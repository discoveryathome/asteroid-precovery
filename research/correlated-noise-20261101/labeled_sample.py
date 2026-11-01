"""Empirical null and injection-recovery sample for the trail screen, on real plate pixels.

Usage (repo root, PYTHONPATH=src): labeled_sample.py PLATE OUT.json [--inject]
PLATE is XE035 (1953 neighbouring field, no AN10 trail), XE036 or XO036 (AN10 trail masked).
Uses the pipeline's own calibrate/fit_trail/screen unchanged; reads frozen inputs only.
"""
import csv
import json
import sys
import time
from pathlib import Path
HERE = Path(__file__).resolve().parent
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.special import erf
sys.path.insert(0, 'src')
from precovery import pipeline as P
from precovery.verify import strict_time, screen

ROT_DEG = (0, 45, 90, 135)        # template rotations relative to the predicted motion
STEP = {'XE035': 24, 'XE036': 24, 'XO036': 30}
INJECT_TRUE_SNR = (3, 4, 5, 6, 8, 12)
NOISE = json.loads(Path(f'{HERE}/noise-model.json').read_text())
F_KEY = 'L15_tapered'


def gaia_all():
    rows = []
    for r in csv.DictReader(open('data/raw/gaia-dr3.csv')):
        rows.append((float(r['ra']), float(r['dec']), float(r['phot_g_mean_mag'] or 99)))
    return np.array(rows)


def setup(region):
    controls = json.loads(Path('data/controls.json').read_text())
    eph, cat = P.ephemeris(), P.catalog()
    ctrl = {c['region']: c for c in controls['positive_controls']}
    vec = None
    trail_mask = []
    for reg in ['XE036', 'XO036']:
        c = ctrl[reg]
        d, h = fits.getdata(c['file'], header=True)
        w = WCS(h)
        coeff, cal = P.calibrate(d, w, cat, abs(h['XPIXELSZ'] * h['PLTSCALE'] / 1000))
        s, e = map(strict_time, [c['start'], c['end']])
        pix = np.array([P.predicted_pixel([p['ra_deg'], p['dec_deg']], w, coeff)
                        for p in [P.prediction(t, eph) for t in [s, s + (e - s) / 2, e]]])
        if reg == region or (region == 'XE035' and reg == 'XE036'):
            vec = pix[2] - pix[0]
        if reg == region:   # chain of pseudo-stars masking the real trail (radius 8 px) plus 10 px beyond its ends
            u = (pix[2] - pix[0]) / np.linalg.norm(pix[2] - pix[0])
            for t in np.arange(-10, np.linalg.norm(pix[2] - pix[0]) + 10, 3):
                q = pix[0] + u * t
                trail_mask.append(dict(x=float(q[0]), y=float(q[1]), sigma_pixel=8 / 3))
            trail_seg = (pix[0], pix[2])
    f = 'data/raw/irsa-XE035.fits' if region == 'XE035' else ctrl[region]['file']
    d, h = fits.getdata(f, header=True)
    w = WCS(h)
    scale = abs(h['XPIXELSZ'] * h['PLTSCALE'] / 1000)
    coeff, cal = P.calibrate(d, w, cat, scale)
    psf = float(np.median([s['sigma_pixel'] for s in cal['stars'] if s['used']]))
    return dict(file=f, data=d.astype(float), w=w, scale=scale, coeff=coeff, cal=cal, cat=cat,
                vec=vec, trail_mask=trail_mask, trail_seg=trail_seg if region != 'XE035' else None, psf=psf)


def trail_image(shape, center, vec, sig):
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    L = np.linalg.norm(vec)
    ang = np.arctan2(vec[1], vec[0])
    along = (xx - center[0]) * np.cos(ang) + (yy - center[1]) * np.sin(ang)
    cross = -(xx - center[0]) * np.sin(ang) + (yy - center[1]) * np.cos(ang)
    prof = .5 * (erf((along + L / 2) / (np.sqrt(2) * sig)) - erf((along - L / 2) / (np.sqrt(2) * sig)))
    return prof * np.exp(-cross ** 2 / (2 * sig ** 2))


def seg_dist(p, a, b):
    t = np.clip(np.dot(p - a, b - a) / np.dot(b - a, b - a), 0, 1)
    return float(np.linalg.norm(p - (a + t * (b - a))))


def trials(env, region):
    d = env['data']
    out = []
    for rot in ROT_DEG:
        r = np.deg2rad(rot)
        v = np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]]) @ env['vec']
        hx, hy = int(abs(v[0]) / 2 + 20), int(abs(v[1]) / 2 + 20)
        for cy in np.arange(hy + 2, d.shape[0] - hy - 2, STEP[region]):
            for cx in np.arange(hx + 2, d.shape[1] - hx - 2, STEP[region]):
                c = np.array([cx + .37, cy + .61], float)   # off-grid sub-pixel centres
                if env['trail_seg'] is not None:
                    # the predicted trail itself must not lie on the trial's own track
                    a, b = c - v / 2, c + v / 2
                    if min(seg_dist(a, *env['trail_seg']), seg_dist(b, *env['trail_seg']),
                           seg_dist(env['trail_seg'][0], a, b), seg_dist(env['trail_seg'][1], a, b)) < 15:
                        continue
                out.append((rot, c, v))
    return out


def evaluate(env, data, c, v, gaia, inject_center=None):
    m, (patch, fitted, valid) = P.fit_trail(data, c, v, env['cal']['stars'] + env['trail_mask'])
    sky = P.measured_sky([m['x'], m['y']], env['w'], env['coeff'])
    near_cat = min(np.linalg.norm(P.offsets(sky, [s[1], s[2]])) for s in env['cat'])
    dra = ((gaia[:, 0] - sky[0] + 180) % 360 - 180) * np.cos(np.deg2rad(sky[1])) * 3600
    ddec = (gaia[:, 1] - sky[1]) * 3600
    sep = np.hypot(dra, ddec)
    j = int(np.argmin(sep))
    if inject_center is None:
        residual = 0.0      # conditioned on a candidate lying exactly at the prediction
    else:
        residual = float(np.hypot(m['x'] - inject_center[0], m['y'] - inject_center[1]) * env['scale'])
    m.update(residual_arcsec=residual, stationary_match=bool(near_cat < 4), epoch_matches=True,
             time_provenance=True, inside_coverage=True, calibration_valid=True,
             artifact_flag=False, edge_flag=False)
    v_ = screen(m)
    morph = [k for k in v_['reasons'] if k in ('snr', 'delta_bic_vs_point', 'length_ratio', 'angle_error_deg', 'fit_converged')]
    # linear matched filter at the fixed trial centre and orientation (no fitting)
    x0, x1, y0, y1 = m['pixel_bounds']
    T = trail_image(patch.shape, (c[0] - x0, c[1] - y0), v, env['psf'])
    Tv = T[valid] - T[valid].mean()
    a_hat = float(np.sum(Tv * (patch[valid] - np.median(patch[valid]))) / np.sum(Tv * Tv))
    return dict(snr=m['snr'], dbic=m['delta_bic_vs_point'], length_ratio=m['length_ratio'],
                angle_error_deg=m['angle_error_deg'], sigma_pixel=m['sigma_pixel'], noise=m['noise'],
                amplitude=m['amplitude'], converged=m['fit_converged'], residual_arcsec=residual,
                fit_xy=[m['x'], m['y']], near_cat_arcsec=float(near_cat), near_gaia_any_arcsec=float(sep[j]),
                near_gaia_any_g=float(gaia[j, 2]), screen_pass=v_['screen_pass'], reasons=v_['reasons'],
                morphology_pass=not morph, a_hat=a_hat, a_white_sigma=float(m['noise'] / np.sqrt(np.sum(Tv * Tv))),
                unmasked=m['unmasked_pixels'])


def main(region, out, inject):
    env = setup(region)
    gaia = gaia_all()
    tr = trials(env, region)
    noise = NOISE[region]
    F = noise['pricing'][F_KEY]['inflation_white']
    sigma = noise['sigma_std']
    rng = np.random.default_rng(int(sum(map(ord, region))) + (1000 if inject else 0))
    print(region, 'trials', len(tr), 'psf', round(env['psf'], 3), 'F', round(F, 3), 'sigma', round(sigma, 1), flush=True)
    rows = []
    t0 = time.time()
    if not inject:
        for i, (rot, c, v) in enumerate(tr):
            r = evaluate(env, env['data'], c, v, gaia)
            r.update(rot_deg=rot, center=c.tolist())
            rows.append(r)
    else:
        # one injection per trial site, cycling the true-SNR grid; true SNR uses the measured F and sigma
        for i, (rot, c, v) in enumerate(tr):
            snr_true = INJECT_TRUE_SNR[i % len(INJECT_TRUE_SNR)]
            ci = c + rng.uniform(-3, 3, 2)
            T = trail_image(env['data'].shape, ci, v, env['psf'])
            hx, hy = int(abs(v[0]) / 2 + 20), int(abs(v[1]) / 2 + 20)
            win = T[int(c[1]) - hy:int(c[1]) + hy + 1, int(c[0]) - hx:int(c[0]) + hx + 1]
            Tp = win - win.mean()
            amp = snr_true * sigma * F / np.sqrt(np.sum(Tp * Tp))
            r = evaluate(env, env['data'] + amp * T, c, v, gaia, inject_center=ci)
            r.update(rot_deg=rot, center=c.tolist(), injected_center=ci.tolist(), injected_true_snr=snr_true,
                     injected_amplitude=float(amp), injected_amp_over_sigma=float(amp / sigma))
            rows.append(r)
    print(region, 'done', len(rows), 'in', round(time.time() - t0, 1), 's', flush=True)
    Path(out).write_text(json.dumps(dict(region=region, file=env['file'], inject=inject, F=F, F_key=F_KEY,
                                         sigma_std=sigma, psf_sigma_pixel=env['psf'], rotations=ROT_DEG,
                                         step=STEP[region], rows=rows), indent=1) + '\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], '--inject' in sys.argv)
