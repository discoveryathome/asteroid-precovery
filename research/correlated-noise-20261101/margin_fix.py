"""Test a wider fit window so the length gate can see the ends of over-length features.

fit_trail_margin is the pipeline's own fit_trail source with the fixed 20 px window margin made a
parameter (two literal edits, asserted). Margin rule tested: 20 + 0.25 * predicted length, i.e.
the pipeline's 20 px clearance kept for the longest feature the length gate accepts (1.5x).
Runs: control (margin 20 must reproduce the pipeline exactly), both documented positives,
the over-length probe, the null grids and the injections, all with the wider margin.
Usage (repo root, PYTHONPATH=src): margin_fix.py PART OUT.json   (PART: control | XE035 | XE036 | XO036)
"""
import inspect
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

SRC = inspect.getsource(P.fit_trail)
assert SRC.count('/2+20)') == 2 and SRC.startswith('def fit_trail(data, center, vector, mask_stars=()):')
SRC = SRC.replace('def fit_trail(data, center, vector, mask_stars=()):',
                  'def fit_trail_margin(data, center, vector, mask_stars=(), margin=20):').replace('/2+20)', '/2+margin)')
NS = dict(vars(P))
exec(SRC, NS)
fit_trail_margin = NS['fit_trail_margin']
ORIGINAL = P.fit_trail


def use_margin(margin):
    P.fit_trail = (lambda d, c, v, m=(): fit_trail_margin(d, c, v, m, margin=margin)) if margin != 20 else ORIGINAL


def margin_for(vec):
    return int(round(20 + 0.25 * np.linalg.norm(vec)))


def positives(margin_rule):
    controls = json.loads(Path('data/controls.json').read_text())
    eph, cat = P.ephemeris(), P.catalog()
    out = {}
    for c in controls['positive_controls']:
        from astropy.io import fits
        from astropy.wcs import WCS
        d, h = fits.getdata(c['file'], header=True)
        w = WCS(h)
        coeff, cal = P.calibrate(d, w, cat, abs(h['XPIXELSZ'] * h['PLTSCALE'] / 1000))
        s, e = map(strict_time, [c['start'], c['end']])
        pred = [P.prediction(t, eph) for t in [s, s + (e - s) / 2, e]]
        pix = np.array([P.predicted_pixel([p['ra_deg'], p['dec_deg']], w, coeff) for p in pred])
        margin = margin_rule(pix[2] - pix[0])
        m, _ = fit_trail_margin(d, pix[1], pix[2] - pix[0], cal['stars'], margin=margin)
        sky = P.measured_sky([m['x'], m['y']], w, coeff)
        res = float(np.linalg.norm(P.offsets(sky, [pred[1]['ra_deg'], pred[1]['dec_deg']])))
        near = min(np.linalg.norm(P.offsets(sky, [q[1], q[2]])) for q in cat)
        m.update(residual_arcsec=res, stationary_match=bool(near < 4), epoch_matches=True, time_provenance=True,
                 inside_coverage=True, calibration_valid=True, artifact_flag=False, edge_flag=False)
        v = screen(m)
        out[c['id']] = dict(margin=margin, snr=m['snr'], dbic=m['delta_bic_vs_point'], length_ratio=m['length_ratio'],
                            angle_error_deg=m['angle_error_deg'], residual_arcsec=res, screen_pass=v['screen_pass'], reasons=v['reasons'])
    return out


def control():
    """margin=20 copy must equal the pipeline on both positives, bit for bit."""
    base = positives(lambda v: 20)
    ref = {m['id']: m for m in json.loads(Path('data/derived/measurements.json').read_text())}
    # exact equality holds on this platform's float path; the committed values came from macOS arm64
    ok = all(abs(base[k]['snr'] - ref[k]['snr']) < 1e-6 and abs(base[k]['dbic'] - ref[k]['delta_bic_vs_point']) < 1e-4 for k in base)
    wide = positives(margin_for)
    return dict(control_reproduces_pipeline=bool(ok), margin20=base, wide=wide)


def plate(region):
    env = LS.setup(region)
    gaia = LS.gaia_all()
    margin = margin_for(env['vec'])
    use_margin(margin)
    null = []
    for rot, c, v in LS.trials(env, region):
        r = LS.evaluate(env, env['data'], c, v, gaia)
        null.append(dict(rot_deg=rot, center=c.tolist(), screen_pass=r['screen_pass'], snr=r['snr'], dbic=r['dbic'], reasons=r['reasons']))
    inj = []
    for row in json.loads(Path(f'{HERE}/sample-{region}-inject.json').read_text())['rows']:
        if row['injected_true_snr'] < 6:
            continue
        rot = np.deg2rad(row['rot_deg'])
        v = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]]) @ env['vec']
        ci = np.array(row['injected_center'])
        T = LS.trail_image(env['data'].shape, ci, v, env['psf'])
        r = LS.evaluate(env, env['data'] + row['injected_amplitude'] * T, np.array(row['center']), v, gaia, inject_center=ci)
        inj.append(dict(true_snr=row['injected_true_snr'], screen_pass=r['screen_pass'], baseline_pass=row['screen_pass'], reasons=r['reasons']))
    probe = []
    tr = [t for t in LS.trials(env, region) if t[0] == 0]
    rng = np.random.default_rng(7 + len(region))      # same sites as morphology_probe.py
    picks = rng.choice(len(tr), size=min(12, len(tr)), replace=False)
    noise = LS.NOISE[region]
    F = noise['pricing'][LS.F_KEY]['inflation_white']
    for i in picks:
        _, c, v = tr[i]
        base = LS.trail_image(env['data'].shape, c, v, env['psf'])
        hx, hy = int(abs(v[0]) / 2 + 20), int(abs(v[1]) / 2 + 20)
        win = base[int(c[1]) - hy:int(c[1]) + hy + 1, int(c[0]) - hx:int(c[0]) + hx + 1]
        amp = 10 * noise['sigma_std'] * F / np.sqrt(np.sum((win - win.mean()) ** 2))
        for mult in (1.0, 1.3, 1.5, 2.0, 4.0, 10.0):
            feat = LS.trail_image(env['data'].shape, c, v * mult, env['psf'])
            r = LS.evaluate(env, env['data'] + amp * feat, c, v, gaia, inject_center=c)
            probe.append(dict(length_mult=mult, screen_pass=r['screen_pass'], length_ratio=r['length_ratio'], reasons=r['reasons']))
    use_margin(20)
    res = dict(region=region, margin=margin, null=null, inject=inj, probe=probe)
    print(region, 'margin', margin, 'null passes', sum(r['screen_pass'] for r in null), '/', len(null), flush=True)
    for lvl in sorted({r['true_snr'] for r in inj}):
        rr = [r for r in inj if r['true_snr'] == lvl]
        print(f'   inject {lvl}: wide recall {np.mean([r["screen_pass"] for r in rr]):.2f} baseline {np.mean([r["baseline_pass"] for r in rr]):.2f}', flush=True)
    print('   probe', ' | '.join(f"x{m}: {sum(r['screen_pass'] for r in probe if r['length_mult'] == m)}/12 L {np.median([r['length_ratio'] for r in probe if r['length_mult'] == m]):.2f}"
                                for m in (1.0, 1.3, 1.5, 2.0, 4.0, 10.0)), flush=True)
    return res


if __name__ == '__main__':
    part, out = sys.argv[1], sys.argv[2]
    res = control() if part == 'control' else plate(part)
    if part == 'control':
        print(json.dumps(res, indent=1))
    Path(out).write_text(json.dumps(res, indent=1) + '\n')
