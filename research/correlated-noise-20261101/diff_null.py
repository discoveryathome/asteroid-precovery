"""Different-epoch template subtraction for XE036 using the 1953 XE035 plate of the same sky.

XE035 (1953-02-18, 103aE, 40 min) is resampled onto the XE036 pixel grid through both plates'
affine-corrected WCS, PSF-matched and flux-scaled; the difference keeps the 1955 trail and
removes static sky. The pipeline's unchanged fit_trail + screen then runs on (a) the AN10 red
positive, (b) the same null sites as the raw-pixel XE036 sample, (c) the same injections.
Usage (repo root, PYTHONPATH=src): diff_null.py OUT.json
"""
import json
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
import numpy as np
from scipy.ndimage import map_coordinates, gaussian_filter
sys.path.insert(0, 'src')
sys.path.insert(0, str(Path(__file__).resolve().parent))
from precovery import pipeline as P
from precovery.verify import strict_time, screen
import labeled_sample as LS

SCALE = 'stars' if '--scale-stars' in sys.argv else 'all-pixels'


def resample_reference(e36, e35):
    ny, nx = e36['data'].shape
    yy, xx = np.mgrid[:ny, :nx]
    p = np.column_stack([xx.ravel(), yy.ravel()]).astype(float)
    corrected = e36['coeff'][0] + p @ e36['coeff'][1:]
    sky = e36['w'].all_pix2world(corrected, 0)
    raw35 = e35['w'].all_world2pix(sky, 0)
    p35 = (raw35 - e35['coeff'][0]) @ np.linalg.inv(e35['coeff'][1:])
    ref = map_coordinates(e35['data'], [p35[:, 1], p35[:, 0]], order=3, cval=np.nan, mode='constant')
    return ref.reshape(ny, nx)


def build_difference(e36, e35):
    ref = resample_reference(e36, e35)
    covered = np.isfinite(ref)
    img = e36['data'].copy()
    s36, s35 = e36['psf'], e35['psf']
    fill = np.where(covered, ref, np.nanmedian(ref))
    if s35 < s36:
        fill = gaussian_filter(fill, np.sqrt(s36 ** 2 - s35 ** 2))
        blur = 'reference'
    else:
        img = gaussian_filter(img, np.sqrt(s35 ** 2 - s36 ** 2))
        blur = 'science'
    # robust linear flux/background match on covered pixels, iterative 3-sigma clipping
    a = np.column_stack([np.ones(covered.sum()), fill[covered]])
    b = img[covered]
    keep = np.ones(len(b), bool)
    for _ in range(6):
        coef = np.linalg.lstsq(a[keep], b[keep], rcond=None)[0]
        r = b - a @ coef
        s = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        keep = np.abs(r) < 3 * s
    if SCALE == 'stars':
        # slope from pixels clearly above sky on BOTH plates (5 to 40 sigma): sky-dominated fits dilute it
        mi, mr = np.median(img[covered]), np.median(fill[covered])
        si = 1.4826 * np.median(np.abs(img[covered] - mi))
        sr = 1.4826 * np.median(np.abs(fill[covered] - mr))
        both = covered & (img - mi > 5 * si) & (fill - mr > 5 * sr) & (img - mi < 40 * si) & (fill - mr < 40 * sr)
        slope = float(np.median((img[both] - mi) / (fill[both] - mr)))
        coef = np.array([mi - slope * mr, slope])
    diff = img - (coef[0] + coef[1] * fill) + np.median(img)
    diff[~covered] = np.nan
    return diff, covered, dict(flux_scale_method=SCALE, psf_xe036=s36, psf_xe035=s35, blurred=blur, scale=float(coef[1]),
                               offset=float(coef[0]), covered_fraction=float(covered.mean()),
                               diff_mad=float(1.4826 * np.nanmedian(np.abs(diff - np.nanmedian(diff)))),
                               raw_mad=float(1.4826 * np.median(np.abs(e36['data'] - np.median(e36['data'])))))


def window_ok(diff, c, v):
    hx, hy = int(abs(v[0]) / 2 + 20), int(abs(v[1]) / 2 + 20)
    win = diff[max(0, int(c[1]) - hy):int(c[1]) + hy + 1, max(0, int(c[0]) - hx):int(c[0]) + hx + 1]
    return bool(np.isfinite(win).all())


def main(out):
    e36, e35 = LS.setup('XE036'), LS.setup('XE035')
    gaia = LS.gaia_all()
    diff, covered, info = build_difference(e36, e35)
    print('difference image', json.dumps(info), flush=True)
    result = dict(info=info)
    # (a) the documented positive on the difference image
    controls = json.loads(Path('data/controls.json').read_text())
    c = controls['positive_controls'][0]
    eph = P.ephemeris()
    s, e = map(strict_time, [c['start'], c['end']])
    pred = [P.prediction(t, eph) for t in [s, s + (e - s) / 2, e]]
    pix = np.array([P.predicted_pixel([p['ra_deg'], p['dec_deg']], e36['w'], e36['coeff']) for p in pred])
    pos = {}
    for name, data in [('raw', e36['data']), ('difference', np.where(np.isfinite(diff), diff, np.nanmedian(diff)))]:
        m, _ = P.fit_trail(data, pix[1], pix[2] - pix[0], e36['cal']['stars'])
        sky = P.measured_sky([m['x'], m['y']], e36['w'], e36['coeff'])
        res = float(np.linalg.norm(P.offsets(sky, [pred[1]['ra_deg'], pred[1]['dec_deg']])))
        near = min(np.linalg.norm(P.offsets(sky, [q[1], q[2]])) for q in e36['cat'])
        m.update(residual_arcsec=res, stationary_match=bool(near < 4), epoch_matches=True, time_provenance=True,
                 inside_coverage=True, calibration_valid=True, artifact_flag=False, edge_flag=False)
        v = screen(m)
        pos[name] = dict(snr=m['snr'], dbic=m['delta_bic_vs_point'], length_ratio=m['length_ratio'],
                         angle_error_deg=m['angle_error_deg'], residual_arcsec=res, noise=m['noise'],
                         screen_pass=v['screen_pass'], reasons=v['reasons'])
        print('positive', name, json.dumps({k: (round(x, 3) if isinstance(x, float) else x) for k, x in pos[name].items()}), flush=True)
    result['positive_an10_red'] = pos
    filled = np.where(np.isfinite(diff), diff, np.nanmedian(diff))
    # (b) paired null: same sites as the raw XE036 sample, only where the window is fully covered
    raw_null = {(r['rot_deg'], tuple(np.round(r['center'], 2))): r
                for r in json.loads(Path(f'{HERE}/sample-XE036-null.json').read_text())['rows']}
    pairs = []
    for rot, cc, v in LS.trials(e36, 'XE036'):
        if not window_ok(diff, cc, v):
            continue
        r = LS.evaluate(e36, filled, cc, v, gaia)
        raw = raw_null[(rot, tuple(np.round(cc, 2)))]
        pairs.append(dict(rot_deg=rot, center=cc.tolist(), diff=r, raw_screen_pass=raw['screen_pass'],
                          raw_snr=raw['snr'], raw_dbic=raw['dbic'], raw_reasons=raw['reasons']))
    result['null'] = pairs
    n = len(pairs)
    print('null sites', n, 'raw pass', sum(p['raw_screen_pass'] for p in pairs),
          'diff pass', sum(p['diff']['screen_pass'] for p in pairs),
          'both', sum(p['raw_screen_pass'] and p['diff']['screen_pass'] for p in pairs), flush=True)
    # (c) paired injections: same sites, amplitudes and centres as the raw XE036 injection sample
    inj = json.loads(Path(f'{HERE}/sample-XE036-inject.json').read_text())['rows']
    ipairs = []
    for row in inj:
        cc = np.array(row['center'])
        rot = np.deg2rad(row['rot_deg'])
        v = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]]) @ e36['vec']
        if not window_ok(diff, cc, v):
            continue
        T = LS.trail_image(e36['data'].shape, np.array(row['injected_center']), v, e36['psf'])
        # inject into the 1955 pixels before subtraction, then rebuild the difference image locally:
        # the reference and fitted scale are unchanged by a faint injection, so add the PSF-matched trail.
        Tm = T if info['blurred'] == 'reference' else gaussian_filter(T, np.sqrt(info['psf_xe035'] ** 2 - info['psf_xe036'] ** 2))
        r = LS.evaluate(e36, filled + row['injected_amplitude'] * Tm, cc, v, gaia, inject_center=np.array(row['injected_center']))
        ipairs.append(dict(true_snr=row['injected_true_snr'], raw_screen_pass=row['screen_pass'], raw_snr=row['snr'],
                           diff=r))
    result['inject'] = ipairs
    for lvl in sorted({p['true_snr'] for p in ipairs}):
        pp = [p for p in ipairs if p['true_snr'] == lvl]
        print(f'inject true SNR {lvl}: n={len(pp)} raw recall {np.mean([p["raw_screen_pass"] for p in pp]):.2f}'
              f' diff recall {np.mean([p["diff"]["screen_pass"] for p in pp]):.2f}'
              f' raw snr {np.median([p["raw_snr"] for p in pp]):.1f} diff snr {np.median([p["diff"]["snr"] for p in pp]):.1f}', flush=True)
    Path(out).write_text(json.dumps(result, indent=1) + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
