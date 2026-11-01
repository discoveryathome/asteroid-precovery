"""Grain-only null: the false-pass floor if every static sky source were removed perfectly.

Replaces each plate's pixels with stationary Gaussian noise carrying that plate's measured
autocorrelation (noise-model.json, the same kernel used for pricing) and sigma, then runs the
identical null trial grid through the pipeline's unchanged fit_trail + screen.
Usage (repo root, PYTHONPATH=src): grain_null.py PLATE OUT.json
"""
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import labeled_sample as LS


def correlated_noise(shape, R, rng):
    L = R.shape[0] // 2
    n0, n1 = shape[0] + 2 * L + 8, shape[1] + 2 * L + 8
    k = np.zeros((n0, n1))
    k[np.ix_(np.arange(-L, L + 1) % n0, np.arange(-L, L + 1) % n1)] = R
    P = np.clip(np.fft.fft2(k).real, 0, None)
    w = rng.standard_normal((n0, n1))
    x = np.fft.ifft2(np.fft.fft2(w) * np.sqrt(P)).real
    x = x[:shape[0], :shape[1]]
    return x / x.std()


def main(region, out):
    env = LS.setup(region)
    noise = LS.NOISE[region]
    Lmax = 25
    L = 15
    # rebuild the tapered L15 kernel from the stored radial profile is lossy; recompute from the 2D values
    R2 = np.array(noise['R2d']) if 'R2d' in noise else None
    if R2 is None:
        raise SystemExit('noise-model.json lacks R2d; rerun noise_model.py')
    R = R2[Lmax - L:Lmax + L + 1, Lmax - L:Lmax + L + 1]
    yy, xx = np.mgrid[-L:L + 1, -L:L + 1]
    R = R * np.clip(1 - np.hypot(xx, yy) / (L + 1), 0, 1)
    rng = np.random.default_rng(20261101 + len(region))
    synth = np.median(env['data']) + noise['sigma_std'] * correlated_noise(env['data'].shape, R, rng)
    env['data'] = synth
    gaia = LS.gaia_all()
    rows = []
    for rot, c, v in LS.trials(env, region):
        r = LS.evaluate(env, synth, c, v, gaia)
        r.update(rot_deg=rot, center=c.tolist())
        rows.append(r)
    Path(out).write_text(json.dumps(dict(region=region, kind='grain-only synthetic null', rows=rows), indent=1) + '\n')
    snr = np.array([r['snr'] for r in rows])
    print(region, 'grain-only null n', len(rows), 'snr>=5', int((snr >= 5).sum()),
          'morph pass', sum(r['morphology_pass'] for r in rows), 'full pass', sum(r['screen_pass'] for r in rows),
          'synthetic std', round(float(synth.std()), 1), flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
