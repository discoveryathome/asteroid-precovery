"""Measure each plate's pixel-noise autocorrelation and price a whitened trail filter.

Run from the repo root with PYTHONPATH=src. Reads frozen inputs only; writes one JSON.
"""
import csv
import json
import sys
from datetime import datetime
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.special import erf
from scipy.sparse.linalg import LinearOperator, cg
sys.path.insert(0, 'src')
from precovery import pipeline as P
from precovery.verify import strict_time

LAGS = (5, 10, 15, 25)      # truncation radii tested for the autocorrelation (pixels)
TILE = 64                   # plane-detrending tile (pixels)


def gaia_all(epoch):
    """Every Gaia row, no quality or magnitude cuts, proper motion applied where present."""
    out = []
    for r in csv.DictReader(open('data/raw/gaia-dr3.csv')):
        ra, dec, g = float(r['ra']), float(r['dec']), float(r['phot_g_mean_mag'] or 21)
        if r['pmra'] and r['pmdec']:
            dt = epoch - 2016.0
            ra += float(r['pmra']) * dt / (3.6e6 * np.cos(np.deg2rad(dec)))
            dec += float(r['pmdec']) * dt / 3.6e6
        out.append((r['source_id'], ra, dec, g))
    return out


def source_mask(data, w, scale, gaia, extra=()):
    """True = usable noise pixel. Masks every Gaia source, bright outliers, edges and the trail."""
    ny, nx = data.shape
    yy, xx = np.mgrid[:ny, :nx]
    keep = np.ones(data.shape, bool)
    for _, ra, dec, g in gaia:
        x, y = w.all_world2pix([[ra, dec]], 0)[0]
        rad_arcsec = 8 if g > 17 else 14 if g > 14 else 24 if g > 11 else 40
        rad = rad_arcsec / scale
        if -rad < x < nx + rad and -rad < y < ny + rad:
            keep &= (xx - x) ** 2 + (yy - y) ** 2 > rad ** 2
    med = np.median(data)
    sig = 1.4826 * np.median(np.abs(data - med))
    hot = (data > med + 3 * sig) | (data < med - 5 * sig)
    # dilate the outlier mask by 3 pixels
    grow = hot.copy()
    for dy in range(-3, 4):
        for dx in range(-3, 4):
            grow |= np.roll(np.roll(hot, dy, 0), dx, 1)
    keep &= ~grow
    keep[:4] = keep[-4:] = False
    keep[:, :4] = keep[:, -4:] = False
    for (x0, y0, x1, y1, pad) in extra:   # trail segment(s) to exclude
        d = np.abs((xx - x0) * (y1 - y0) - (yy - y0) * (x1 - x0)) / np.hypot(x1 - x0, y1 - y0)
        t = ((xx - x0) * (x1 - x0) + (yy - y0) * (y1 - y0)) / ((x1 - x0) ** 2 + (y1 - y0) ** 2)
        keep &= ~((d < pad) & (t > -pad / np.hypot(x1 - x0, y1 - y0)) & (t < 1 + pad / np.hypot(x1 - x0, y1 - y0)))
    return keep


def detrend(data, keep):
    """Subtract a plane per tile, fitted to usable pixels only."""
    r = np.zeros(data.shape)
    ny, nx = data.shape
    for y0 in range(0, ny, TILE):
        for x0 in range(0, nx, TILE):
            sl = (slice(y0, min(ny, y0 + TILE)), slice(x0, min(nx, x0 + TILE)))
            k = keep[sl]
            if k.sum() < 200:
                keep[sl] = False
                continue
            yy, xx = np.mgrid[sl]
            A = np.column_stack([np.ones(k.sum()), xx[k], yy[k]])
            c = np.linalg.lstsq(A, data[sl][k].astype(float), rcond=None)[0]
            r[sl] = data[sl] - (c[0] + c[1] * xx + c[2] * yy)
    r[~keep] = 0
    return r, keep


def autocov(r, keep, L):
    """Masked autocovariance for |dx|,|dy| <= L via zero-padded FFTs."""
    ny, nx = r.shape
    shape = (ny + L + 1, nx + L + 1)
    fr = np.fft.rfft2(r * keep, shape)
    fm = np.fft.rfft2(keep.astype(float), shape)
    num = np.fft.irfft2(fr * np.conj(fr), shape)
    den = np.fft.irfft2(fm * np.conj(fm), shape)
    idx = np.arange(-L, L + 1)
    num = num[np.ix_(idx % shape[0], idx % shape[1])]
    den = den[np.ix_(idx % shape[0], idx % shape[1])]
    return num / np.maximum(den, 1), den


def template(shape, center, length, ang, sig):
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    x, y = center
    along = (xx - x) * np.cos(ang) + (yy - y) * np.sin(ang)
    cross = -(xx - x) * np.sin(ang) + (yy - y) * np.cos(ang)
    prof = .5 * (erf((along + length / 2) / (np.sqrt(2) * sig)) - erf((along - length / 2) / (np.sqrt(2) * sig)))
    return prof * np.exp(-cross ** 2 / (2 * sig ** 2))


def conv_same(x, R):
    """Linear convolution of x with the (2L+1)^2 kernel R, cropped to x's shape."""
    L = R.shape[0] // 2
    shape = (x.shape[0] + 2 * L, x.shape[1] + 2 * L)
    out = np.fft.irfft2(np.fft.rfft2(x, shape) * np.fft.rfft2(R, shape), shape)
    return out[L:L + x.shape[0], L:L + x.shape[1]]


def price(T, R):
    """Per unit amplitude and unit pixel sigma: white-filter inflation and best-case whitening gain.

    T is the background-projected template over the fit window; R the normalized
    autocorrelation kernel (zero beyond its radius). Returns the white matched filter's
    true SNR inflation F = sqrt(T'RT/T'T) and the GLS gain sqrt(T'R^-1T * T'RT)/T'T,
    with T'R^-1T solved inside the window (no pixels masked: favorable).
    """
    tt = float(np.sum(T * T))
    trt = float(np.sum(T * conv_same(T, R)))
    n = T.size
    op = LinearOperator((n, n), matvec=lambda v: conv_same(v.reshape(T.shape), R).ravel(), dtype=float)
    z, info = cg(op, T.ravel(), rtol=1e-8, maxiter=5000)
    tri = float(T.ravel() @ z)
    return dict(TT=tt, TRT=trt, TRinvT=tri, cg_info=int(info),
                inflation_white=float(np.sqrt(trt / tt)),
                best_case_whitening_gain=float(np.sqrt(tri * trt) / tt) if tri > 0 else None)


def kernel_psd_min(R):
    """Minimum of the kernel's 2D spectrum on a fine grid (negative = not a valid covariance)."""
    big = np.zeros((256, 256))
    L = R.shape[0] // 2
    big[np.ix_(np.arange(-L, L + 1) % 256, np.arange(-L, L + 1) % 256)] = R
    return float(np.fft.fft2(big).real.min())


def main(out):
    controls = json.loads(Path('data/controls.json').read_text())
    eph, cat = P.ephemeris(), P.catalog()
    plates = {}
    for c in controls['positive_controls']:
        d, h = fits.getdata(c['file'], header=True)
        w = WCS(h)
        scale = abs(h['XPIXELSZ'] * h['PLTSCALE'] / 1000)
        coeff, cal = P.calibrate(d, w, cat, scale)
        s, e = map(strict_time, [c['start'], c['end']])
        mid = s + (e - s) / 2
        pix = np.array([P.predicted_pixel([p['ra_deg'], p['dec_deg']], w, coeff)
                        for p in [P.prediction(t, eph) for t in [s, mid, e]]])
        m, _ = P.fit_trail(d, pix[1], pix[2] - pix[0], cal['stars'])
        plates[c['region']] = dict(file=c['file'], data=d, w=w, scale=scale, epoch=1955.07,
                                   trail=(pix[0, 0], pix[0, 1], pix[2, 0], pix[2, 1], 12),
                                   vector=pix[2] - pix[0], fit=m)
    # The 1953 wrong-epoch plate of the neighbouring field: real pixels, no AN10 trail.
    d, h = fits.getdata('data/raw/irsa-XE035.fits', header=True)
    plates['XE035'] = dict(file='data/raw/irsa-XE035.fits', data=d, w=WCS(h),
                           scale=abs(h['XPIXELSZ'] * h['PLTSCALE'] / 1000), epoch=1953.13,
                           trail=None, vector=plates['XE036']['vector'], fit=plates['XE036']['fit'])
    result = {}
    for region, p in plates.items():
        d = p['data'].astype(float)
        keep = source_mask(d, p['w'], p['scale'], gaia_all(p['epoch']), [p['trail']] if p['trail'] else [])
        r, keep = detrend(d, keep)
        Lmax = max(LAGS)
        C, pairs = autocov(r, keep, Lmax)
        var = C[Lmax, Lmax]
        Rfull = C / var
        resid = r[keep]
        row = dict(file=p['file'], usable_pixels=int(keep.sum()), usable_fraction=float(keep.mean()),
                   sigma_std=float(np.sqrt(var)), sigma_mad=float(1.4826 * np.median(np.abs(resid - np.median(resid)))),
                   pipeline_window_noise=p['fit']['noise'] if region != 'XE035' else None,
                   R_1_0=float(Rfull[Lmax, Lmax + 1]), R_0_1=float(Rfull[Lmax + 1, Lmax]),
                   R_2_0=float(Rfull[Lmax, Lmax + 2]), R_0_2=float(Rfull[Lmax + 2, Lmax]),
                   R_1_1=float(Rfull[Lmax + 1, Lmax + 1]), R_1_m1=float(Rfull[Lmax + 1, Lmax - 1]),
                   R_5_0=float(Rfull[Lmax, Lmax + 5]), R_0_5=float(Rfull[Lmax + 5, Lmax]),
                   R_10_0=float(Rfull[Lmax, Lmax + 10]), R_0_10=float(Rfull[Lmax + 10, Lmax]),
                   R_25_0=float(Rfull[Lmax, Lmax + 25]), R_0_25=float(Rfull[Lmax + 25, Lmax]),
                   min_pairs_at_Lmax=float(pairs[[0, -1, Lmax, Lmax], [Lmax, Lmax, 0, -1]].min()),
                   radial_profile=[float(np.mean([Rfull[Lmax + dy, Lmax + dx] for dy in range(-k, k + 1)
                                                  for dx in range(-k, k + 1) if max(abs(dx), abs(dy)) == k]))
                                   for k in range(Lmax + 1)],
                   R2d=Rfull.tolist(),
                   pricing={})
        f = p['fit']
        length = float(np.linalg.norm(p['vector']))
        ang = float(np.arctan2(p['vector'][1], p['vector'][0]))
        hx, hy = int(abs(p['vector'][0]) / 2 + 20), int(abs(p['vector'][1]) / 2 + 20)
        shape = (2 * hy + 1, 2 * hx + 1)
        sig = f['sigma_pixel']
        T = template(shape, (hx, hy), length, ang, sig)
        Tp = T - T.mean()
        for L in LAGS:
            R = Rfull[Lmax - L:Lmax + L + 1, Lmax - L:Lmax + L + 1].copy()
            yy, xx = np.mgrid[-L:L + 1, -L:L + 1]
            taper = np.clip(1 - np.hypot(xx, yy) / (L + 1), 0, 1)   # radial triangle taper
            for name, kern in [('raw', R * (np.hypot(xx, yy) <= L)), ('tapered', R * taper)]:
                pr = price(Tp, kern)
                pr['kernel_spectrum_min'] = kernel_psd_min(kern)
                pr['inflation_white_unprojected'] = float(np.sqrt(np.sum(T * conv_same(T, kern)) / np.sum(T * T)))
                row['pricing'][f'L{L}_{name}'] = pr
        row['template'] = dict(length_pixel=length, angle_deg=float(np.rad2deg(ang)), sigma_pixel=sig,
                               window=list(shape))
        result[region] = row
        print(region, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()
                                  if k not in ('pricing', 'radial_profile', 'R2d')}))
        for k, v in row['pricing'].items():
            print('   ', k, {x: round(y, 4) if isinstance(y, float) else y for x, y in v.items()
                             if x in ('inflation_white', 'inflation_white_unprojected', 'best_case_whitening_gain', 'cg_info', 'kernel_spectrum_min')})
    Path(out).write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
