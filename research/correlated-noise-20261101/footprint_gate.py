"""Price and test a footprint version of the stationary gate, plus an along-track uniformity statistic.

The pipeline's stationary gate asks whether a catalogue star lies within 4 arcsec of the fitted
trail's CENTRE. This re-fits every sample row that passed the unchanged screen and asks the same
question along the whole fitted segment:
  A: pipeline catalogue (P.catalog(), same cuts), within 4 arcsec of the fitted segment
  B: every Gaia row, magnitude-dependent radius (8/14/24/40 arcsec for G >17 / 14-17 / 11-14 / <11)
It also measures U, the largest of six along-track segment amplitudes over the median of the
other five. No threshold is tuned on the two documented positives.
Usage (repo root, PYTHONPATH=src): footprint_gate.py OUT.json
"""
import json
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
import numpy as np
sys.path.insert(0, 'src')
sys.path.insert(0, str(Path(__file__).resolve().parent))
from precovery import pipeline as P
from precovery.verify import screen, strict_time
import labeled_sample as LS

K = 6


def radius_b(g):
    return 8 if g > 17 else 14 if g > 14 else 24 if g > 11 else 40


def seg_dist_arcsec(sky_pts, a, b):
    """Distances (arcsec) from catalogue points to the segment a-b, in the local tangent plane at a."""
    d = np.array([P.offsets(s, a) for s in sky_pts])
    e = P.offsets(b, a)
    t = np.clip(d @ e / max(e @ e, 1e-12), 0, 1)
    return np.linalg.norm(d - t[:, None] * e[None, :], axis=1)


def refit(env, data, c, v, gaia, inject_center=None):
    m, (patch, fitted, valid) = P.fit_trail(data, c, v, env['cal']['stars'] + env['trail_mask'])
    w, coeff, scale = env['w'], env['coeff'], env['scale']
    sky = P.measured_sky([m['x'], m['y']], w, coeff)
    near_cat = min(np.linalg.norm(P.offsets(sky, [s[1], s[2]])) for s in env['cat'])
    residual = 0.0 if inject_center is None else float(np.hypot(m['x'] - inject_center[0], m['y'] - inject_center[1]) * scale)
    m.update(residual_arcsec=residual, stationary_match=bool(near_cat < 4), epoch_matches=True, time_provenance=True,
             inside_coverage=True, calibration_valid=True, artifact_flag=False, edge_flag=False)
    verdict = screen(m)
    # fitted segment ends, in sky
    u = np.array([np.cos(m['angle_rad']), np.sin(m['angle_rad'])])
    ends = [np.array([m['x'], m['y']]) + s * u * m['length_pixel'] / 2 for s in (-1, 1)]
    a, b = [P.measured_sky(e, w, coeff) for e in ends]
    cat_pts = np.array([[s[1], s[2]] for s in env['cat']])
    dA = seg_dist_arcsec(cat_pts, a, b)
    gmask = (np.abs(gaia[:, 1] - sky[1]) < .1)
    gp = gaia[gmask]
    dB = seg_dist_arcsec(gp[:, :2], a, b) if len(gp) else np.array([np.inf])
    rB = np.array([radius_b(g) for g in gp[:, 2]]) if len(gp) else np.array([0])
    gate_a = bool(dA.min() < 4)
    gate_b = bool(np.any(dB < rB))
    # along-track uniformity
    x0, x1, y0, y1 = m['pixel_bounds']
    yy, xx = np.mgrid[y0:y1, x0:x1]
    along = (xx - m['x']) * u[0] + (yy - m['y']) * u[1]
    cross = -(xx - m['x']) * u[1] + (yy - m['y']) * u[0]
    wgt = np.exp(-cross ** 2 / (2 * m['sigma_pixel'] ** 2))
    edges = np.linspace(-m['length_pixel'] / 2, m['length_pixel'] / 2, K + 1)
    amps = []
    for k in range(K):
        sel = valid & (along >= edges[k]) & (along < edges[k + 1]) & (np.abs(cross) < 3 * m['sigma_pixel'])
        ww = wgt[sel]
        amps.append(float(np.sum((patch[sel] - m['background']) * ww) / max(np.sum(ww * ww), 1e-9)) if sel.sum() else np.nan)
    amps = np.array(amps)
    j = int(np.nanargmax(amps))
    rest = np.delete(amps, j)
    med = np.nanmedian(rest)
    U = float(amps[j] / med) if med > 0 else None    # None: the other segments carry no positive flux
    return dict(screen_pass=verdict['screen_pass'], reasons=verdict['reasons'], snr=m['snr'], dbic=m['delta_bic_vs_point'],
                gate_a_reject=gate_a, gate_b_reject=gate_b, min_dist_catalogue_arcsec=float(dA.min()),
                uniformity_U=U, segment_amplitudes=amps.tolist(), fit_amplitude=m['amplitude'])


def positives(gaia):
    controls = json.loads(Path('data/controls.json').read_text())
    eph = P.ephemeris()
    out = {}
    for c in controls['positive_controls']:
        env = LS.setup(c['region'])
        env['trail_mask'] = []          # the positive itself must not be masked
        s, e = map(strict_time, [c['start'], c['end']])
        pred = [P.prediction(t, eph) for t in [s, s + (e - s) / 2, e]]
        pix = np.array([P.predicted_pixel([p['ra_deg'], p['dec_deg']], env['w'], env['coeff']) for p in pred])
        r = refit(env, env['data'], pix[1], pix[2] - pix[0], gaia)
        out[c['id']] = r
        print('positive', c['id'], {k: (round(x, 3) if isinstance(x, float) else x) for k, x in r.items() if k != 'segment_amplitudes'}, flush=True)
    return out


def main(out):
    gaia = LS.gaia_all()
    result = dict(positives=positives(gaia), plates={})
    for region in ['XE035', 'XE036', 'XO036']:
        env = LS.setup(region)
        rows = []
        for kind in ['null', 'inject']:
            sample = json.loads(Path(f'{HERE}/sample-{region}-{kind}.json').read_text())['rows']
            for row in sample:
                if not row['screen_pass']:
                    continue
                rot = np.deg2rad(row['rot_deg'])
                v = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]]) @ env['vec']
                c = np.array(row['center'])
                if kind == 'null':
                    r = refit(env, env['data'], c, v, gaia)
                else:
                    ci = np.array(row['injected_center'])
                    T = LS.trail_image(env['data'].shape, ci, v, env['psf'])
                    r = refit(env, env['data'] + row['injected_amplitude'] * T, c, v, gaia, inject_center=ci)
                    r['true_snr'] = row['injected_true_snr']
                r.update(kind=kind, rot_deg=row['rot_deg'], center=row['center'], reproduced=r['screen_pass'] == row['screen_pass'])
                rows.append(r)
        result['plates'][region] = rows
        nul = [r for r in rows if r['kind'] == 'null']
        inj = [r for r in rows if r['kind'] == 'inject']
        print(region, 'null passes', len(nul), 'reproduced', sum(r['reproduced'] for r in nul),
              'gateA rejects', sum(r['gate_a_reject'] for r in nul), 'gateB rejects', sum(r['gate_b_reject'] for r in nul),
              '| injected passes', len(inj), 'gateA rejects', sum(r['gate_a_reject'] for r in inj),
              'gateB rejects', sum(r['gate_b_reject'] for r in inj), flush=True)
    Path(out).write_text(json.dumps(result, indent=1) + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
