"""PNG montage of null-sample screen passes (real pixels), for visual inspection.

Usage (repo root, PYTHONPATH=src): montage.py OUT.png PLATE [PLATE ...]
Shows up to 8 full-screen passes per plate, highest nominal SNR first, with the fitted trail
segment overlaid, beside the documented AN10 trail on the same plate.
"""
import json
import os
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
import numpy as np
from astropy.io import fits
os.environ.setdefault('MPLCONFIGDIR', str(Path('.cache/matplotlib').resolve()))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FILES = {'XE035': 'data/raw/irsa-XE035.fits', 'XE036': 'data/raw/XE036.fits', 'XO036': 'data/raw/XO036.fits'}
POS = {'XE036': 'an10-red-1955', 'XO036': 'an10-blue-1955'}


def panel(ax, data, x, y, L, ang, title, half):
    x0, y0 = int(round(x)) - half, int(round(y)) - half
    cut = data[max(0, y0):y0 + 2 * half + 1, max(0, x0):x0 + 2 * half + 1].astype(float)
    lo, hi = np.percentile(cut, [5, 99])
    ax.imshow(cut, origin='lower', cmap='gray_r', vmin=lo, vmax=hi)
    dx, dy = np.cos(ang) * L / 2, np.sin(ang) * L / 2
    cx, cy = x - max(0, x0), y - max(0, y0)
    ax.plot([cx - dx, cx + dx], [cy - dy - 6, cy + dy - 6], 'r-', lw=.8)
    ax.set_title(title, fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])


def main(out, plates):
    rows = []
    meas = {m['id']: m for m in json.loads(Path('data/derived/measurements.json').read_text())}
    for region in plates:
        data = fits.getdata(FILES[region])
        null = json.loads(Path(f'{HERE}/sample-{region}-null.json').read_text())['rows']
        passes = sorted([r for r in null if r['screen_pass']], key=lambda r: -r['snr'])[:8]
        entries = []
        if region in POS:
            m = meas[POS[region]]
            entries.append((m['x'], m['y'], m['length_pixel'], m['angle_rad'], f"{region} AN10 positive\nSNR {m['snr']:.1f} dBIC {m['delta_bic_vs_point']:.0f}"))
        for r in passes:
            # fitted angle is not stored; draw along the trial vector (rotation from the predicted one)
            entries.append((r['fit_xy'][0], r['fit_xy'][1], None, r['rot_deg'],
                            f"{region} null rot {r['rot_deg']}\nSNR {r['snr']:.1f} dBIC {r['dbic']:.0f} L {r['length_ratio']:.2f}"))
        rows.append((region, entries))
    ncol = max(len(e) for _, e in rows)
    fig, axes = plt.subplots(len(rows), ncol, figsize=(1.9 * ncol, 2.2 * len(rows)), squeeze=False)
    vec = {'XE035': (133.62, 31.04), 'XE036': (133.62, 31.04), 'XO036': (60.03, 13.7)}
    for axrow, (region, entries) in zip(axes, rows):
        data = fits.getdata(FILES[region])
        v = np.array(vec[region])
        half = int(np.linalg.norm(v) / 2 + 25)
        for ax, (x, y, L, ang, title) in zip(axrow, entries):
            if L is None:
                rot = np.deg2rad(ang)
                vv = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]]) @ v
                L, ang = np.linalg.norm(vv), np.arctan2(vv[1], vv[0])
            panel(ax, data, x, y, L, ang, title, half)
        for ax in axrow[len(entries):]:
            ax.axis('off')
    fig.suptitle('Null-sample screen passes on real pixels (red line: trial track, offset 6 px below)', fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=100)
    print('wrote', out)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
