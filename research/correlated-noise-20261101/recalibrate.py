"""Post-hoc price of three per-plate SNR calibrations and the footprint gates, on the stored samples.

Only the SNR gate changes in variants 1-2; every other gate keeps the pipeline's verdict.
  0  pipeline as is (fixed factor 2)
  1  grain:     SNR x 2 / F, F from the measured autocorrelation (noise-model.json)
  2  empirical: SNR x 2 / s, s the robust scatter of the linear trail filter over the plate's own
                null sites (pipeline noise units; includes static clutter)
Footprint gates A/B (footprint-gate.json) are applied on top of variant 0.
Thresholds come from the null/injection samples only; the two documented positives are reported, not used.
Usage (repo root): recalibrate.py OUTDIR
"""
import json
import sys
from pathlib import Path
import numpy as np

OUT = Path(sys.argv[1])
noise = json.loads((OUT / 'noise-model.json').read_text())
meas = {m['id']: m for m in json.loads(Path('data/derived/measurements.json').read_text())}
fp = json.loads((OUT / 'footprint-gate.json').read_text()) if (OUT / 'footprint-gate.json').exists() else None
POS = {'XE036': 'an10-red-1955', 'XO036': 'an10-blue-1955'}


def rstd(x):
    x = np.asarray(x, float)
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def passes(row, factor):
    others = [k for k in row['reasons'] if k != 'snr']
    return not others and row['snr'] * factor >= 5


result = {}
for region in ['XE035', 'XE036', 'XO036']:
    null = json.loads((OUT / f'sample-{region}-null.json').read_text())
    inj = json.loads((OUT / f'sample-{region}-inject.json').read_text())['rows']
    F = noise[region]['pricing'][null['F_key']]['inflation_white']
    s = rstd([r['a_hat'] / r['a_white_sigma'] for r in null['rows']])
    variants = {'0-pipeline': 1.0, '1-grain-F': 2 / F, '2-empirical-null': 2 / s}
    res = dict(F_grain=F, empirical_scatter=s, variants={})
    for name, f in variants.items():
        fa = sum(passes(r, f) for r in null['rows'])
        rec = {}
        for lvl in sorted({r['injected_true_snr'] for r in inj}):
            rr = [r for r in inj if r['injected_true_snr'] == lvl]
            rec[str(lvl)] = round(sum(passes(r, f) for r in rr) / len(rr), 3)
        v = dict(snr_factor=f, null_passes=fa, null_trials=len(null['rows']), null_pass_rate=fa / len(null['rows']), recall=rec)
        if region in POS:
            m = meas[POS[region]]
            v['positive_corrected_snr'] = m['snr'] * f
            v['positive_passes'] = bool(m['verdict']['screen_pass'] and m['snr'] * f >= 5)
        res['variants'][name] = v
    if fp:
        rows = fp['plates'][region]
        for gate in ['gate_a_reject', 'gate_b_reject']:
            nul = [r for r in rows if r['kind'] == 'null']
            fa = len(nul) - sum(r[gate] for r in nul)
            rec = {}
            for lvl in sorted({r['injected_true_snr'] for r in inj}):
                n = sum(r['injected_true_snr'] == lvl for r in inj)
                kept = sum(1 for r in rows if r['kind'] == 'inject' and r['true_snr'] == lvl and r['screen_pass'] and not r[gate])
                rec[str(lvl)] = round(kept / n, 3)
            v = dict(null_passes=fa, null_trials=len(null['rows']), null_pass_rate=fa / len(null['rows']), recall=rec)
            if region in POS:
                pr = fp['positives'][POS[region]]
                v['positive_passes'] = bool(pr['screen_pass'] and not pr[gate])
            res['variants']['3-footprint-' + gate[5]] = v
    result[region] = res
    print(region, f"F={F:.2f} empirical s={s:.2f}")
    for name, v in res['variants'].items():
        extra = f" positive {'PASS' if v.get('positive_passes') else 'FAIL'}" if region in POS else ''
        if 'positive_corrected_snr' in v:
            extra += f" (SNR {v['positive_corrected_snr']:.1f})"
        print(f"   {name:22s} null {v['null_passes']:3d}/{v['null_trials']} ({100 * v['null_pass_rate']:.1f}%) recall {v['recall']}{extra}")
(OUT / 'recalibration.json').write_text(json.dumps(result, indent=2) + '\n')
