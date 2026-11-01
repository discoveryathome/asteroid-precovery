"""Can the screen tell an over-length straight feature from the predicted trail?

Injects, into real XE035 / XE036 / XO036 null-site pixels, straight features through the
trial centre at the predicted angle (+/- small offsets) whose length is a multiple of the
predicted trail length, then runs the pipeline's unchanged fit_trail + screen.
Usage (repo root, PYTHONPATH=src): morphology_probe.py OUT.json
"""
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, 'src')
sys.path.insert(0, str(Path(__file__).resolve().parent))
from precovery import pipeline as P
from precovery.verify import screen
import labeled_sample as LS

LENGTH_MULT = (1.0, 1.3, 1.5, 2.0, 4.0, 10.0)
ANGLE_OFF_DEG = (0, 6)
TRUE_SNR = 10          # of the predicted-length trail; longer features keep the same peak amplitude
SITES = 12


def main(out):
    rows = []
    for region in ['XE035', 'XE036', 'XO036']:
        env = LS.setup(region)
        noise = LS.NOISE[region]
        F = noise['pricing'][LS.F_KEY]['inflation_white']
        sigma = noise['sigma_std']
        gaia = LS.gaia_all()
        tr = [t for t in LS.trials(env, region) if t[0] == 0]
        rng = np.random.default_rng(7 + len(region))
        picks = rng.choice(len(tr), size=min(SITES, len(tr)), replace=False)
        for i in picks:
            _, c, v = tr[i]
            L = np.linalg.norm(v)
            base = LS.trail_image(env['data'].shape, c, v, env['psf'])
            hx, hy = int(abs(v[0]) / 2 + 20), int(abs(v[1]) / 2 + 20)
            win = base[int(c[1]) - hy:int(c[1]) + hy + 1, int(c[0]) - hx:int(c[0]) + hx + 1]
            amp = TRUE_SNR * sigma * F / np.sqrt(np.sum((win - win.mean()) ** 2))
            for off in ANGLE_OFF_DEG:
                r = np.deg2rad(off)
                rot = np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
                for mult in LENGTH_MULT:
                    feat = LS.trail_image(env['data'].shape, c, rot @ v * mult, env['psf'])
                    res = LS.evaluate(env, env['data'] + amp * feat, c, v, gaia, inject_center=c)
                    rows.append(dict(region=region, site=c.tolist(), angle_offset_deg=off, length_mult=mult,
                                     fit_length_ratio=res['length_ratio'], snr=res['snr'], dbic=res['dbic'],
                                     angle_error_deg=res['angle_error_deg'], residual_arcsec=res['residual_arcsec'],
                                     screen_pass=res['screen_pass'], reasons=res['reasons']))
        print(region, 'done', flush=True)
    Path(out).write_text(json.dumps(rows, indent=1) + '\n')
    for region in ['XE035', 'XE036', 'XO036']:
        for off in ANGLE_OFF_DEG:
            line = []
            for mult in LENGTH_MULT:
                rr = [r for r in rows if r['region'] == region and r['angle_offset_deg'] == off and r['length_mult'] == mult]
                line.append(f"x{mult}: pass {sum(r['screen_pass'] for r in rr)}/{len(rr)} "
                            f"Lfit {np.median([r['fit_length_ratio'] for r in rr]):.2f}")
            print(region, f'angle+{off}', ' | '.join(line))


if __name__ == '__main__':
    main(sys.argv[1])
