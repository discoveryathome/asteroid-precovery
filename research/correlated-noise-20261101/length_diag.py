"""Why over-length features pass the length gate: fitted length vs window margin, and a fixed-length cost scan.

Three XE036 probe sites (the first three morphology_probe.py sites), a 1x and a 10x straight
feature at true SNR 10. For each: the full fit at margins 20, 54 and 100 px, then the robust cost
with length held at 0.8-1.55x and the other six parameters refitted (margin 54).
Usage (repo root, PYTHONPATH=src): length_diag.py OUT.json
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.special import erf
from scipy.optimize import least_squares
sys.path.insert(0, 'src')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import labeled_sample as LS
import margin_fix as MF

SCAN = (0.8, 1.0, 1.12, 1.3, 1.45, 1.55)


def main(out):
    env = LS.setup('XE036')
    noise = LS.NOISE['XE036']
    F = noise['pricing'][LS.F_KEY]['inflation_white']
    tr = [t for t in LS.trials(env, 'XE036') if t[0] == 0]
    picks = np.random.default_rng(7 + len('XE036')).choice(len(tr), size=12, replace=False)
    mask = env['cal']['stars'] + env['trail_mask']
    rows = []
    for i in picks[:3]:
        _, c, v = tr[i]
        L0 = np.linalg.norm(v)
        base = LS.trail_image(env['data'].shape, c, v, env['psf'])
        hx, hy = int(abs(v[0]) / 2 + 20), int(abs(v[1]) / 2 + 20)
        win = base[int(c[1]) - hy:int(c[1]) + hy + 1, int(c[0]) - hx:int(c[0]) + hx + 1]
        amp = 10 * noise['sigma_std'] * F / np.sqrt(np.sum((win - win.mean()) ** 2))
        for mult in (1.0, 10.0):
            data = env['data'] + amp * LS.trail_image(env['data'].shape, c, v * mult, env['psf'])
            fits_ = {}
            for margin in (20, 54, 100):
                m, _ = MF.fit_trail_margin(data, c, v, mask, margin=margin)
                fits_[str(margin)] = dict(length_ratio=m['length_ratio'], converged=m['fit_converged'])
            m, (patch, _, valid) = MF.fit_trail_margin(data, c, v, mask, margin=54)
            x0, x1, y0, y1 = m['pixel_bounds']
            yy, xx = np.mgrid[y0:y1, x0:x1]

            def model(p, L):
                x, y, ang, sig, a, b = p
                along = (xx - x) * np.cos(ang) + (yy - y) * np.sin(ang)
                cross = -(xx - x) * np.sin(ang) + (yy - y) * np.cos(ang)
                return b + a * .5 * (erf((along + L / 2) / (np.sqrt(2) * sig)) - erf((along - L / 2) / (np.sqrt(2) * sig))) * np.exp(-cross ** 2 / (2 * sig ** 2))
            p0 = [m['x'], m['y'], m['angle_rad'], m['sigma_pixel'], m['amplitude'], m['background']]
            scan = {str(f): float(least_squares(lambda p: ((model(p, f * L0) - patch) / m['noise'])[valid], p0, loss='soft_l1').cost)
                    for f in SCAN}
            rows.append(dict(site=c.tolist(), length_mult=mult, full_fit=fits_, fixed_length_cost_margin54=scan,
                             best_fixed_length=min(scan, key=scan.get)))
            print(rows[-1], flush=True)
    Path(out).write_text(json.dumps(rows, indent=1) + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
