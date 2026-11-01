"""Condense the per-row experiment outputs (grain, difference, footprint, margin, probe) into one JSON.

Usage (repo root): experiments_summary.py DIR   (DIR holds the per-row outputs; writes DIR/experiments-summary.json)
"""
import json
import sys
from pathlib import Path
import numpy as np

D = Path(sys.argv[1])
out = {}


def load(name):
    p = D / name
    return json.loads(p.read_text()) if p.exists() else None


g = {}
for region in ['XE035', 'XE036', 'XO036']:
    rows = load(f'grain-{region}-null.json')
    if rows:
        rows = rows['rows']
        z = np.array([r['a_hat'] / r['a_white_sigma'] for r in rows])
        g[region] = dict(trials=len(rows), screen_passes=sum(r['screen_pass'] for r in rows),
                         linear_filter_scatter=float(1.4826 * np.median(np.abs(z - np.median(z)))))
out['grain_only_null'] = g

d = {}
for name, f in [('all-pixel-scale', 'diff-XE036.json'), ('star-pixel-scale', 'diff-XE036-scalestars.json')]:
    x = load(f)
    if not x:
        continue
    rec = {}
    for lvl in sorted({p['true_snr'] for p in x['inject']}):
        pp = [p for p in x['inject'] if p['true_snr'] == lvl]
        rec[str(lvl)] = dict(n=len(pp), raw=float(np.mean([p['raw_screen_pass'] for p in pp])),
                             difference=float(np.mean([p['diff']['screen_pass'] for p in pp])))
    d[name] = dict(info=x['info'], positive=x['positive_an10_red'], null_sites=len(x['null']),
                   null_raw_passes=sum(p['raw_screen_pass'] for p in x['null']),
                   null_difference_passes=sum(p['diff']['screen_pass'] for p in x['null']), recall=rec)
out['difference_imaging_XE036'] = d

fp = load('footprint-gate.json')
if fp:
    f = dict(positives={k: {a: v[a] for a in ('gate_a_reject', 'gate_b_reject', 'min_dist_catalogue_arcsec', 'uniformity_U')}
                        for k, v in fp['positives'].items()}, plates={})
    for region, rows in fp['plates'].items():
        nul = [r for r in rows if r['kind'] == 'null']
        inj = [r for r in rows if r['kind'] == 'inject']
        U = lambda rr: np.array([np.inf if r['uniformity_U'] is None else r['uniformity_U'] for r in rr])
        thr = float(np.percentile(U([r for r in inj if r['true_snr'] >= 8]), 95))
        f['plates'][region] = dict(null_passes=len(nul), injected_passes=len(inj),
                                   refits_reproduce_screen=all(r['reproduced'] for r in rows),
                                   gate_a_rejects_null=sum(r['gate_a_reject'] for r in nul), gate_a_rejects_injected=sum(r['gate_a_reject'] for r in inj),
                                   gate_b_rejects_null=sum(r['gate_b_reject'] for r in nul), gate_b_rejects_injected=sum(r['gate_b_reject'] for r in inj),
                                   uniformity_threshold_p95_injected_snr_ge_8=thr,
                                   uniformity_rejects_null=int((U(nul) > thr).sum()), uniformity_rejects_injected=int((U(inj) > thr).sum()),
                                   uniformity_null_median=float(np.median(U(nul))), uniformity_injected_median=float(np.median(U(inj))))
    out['footprint_and_uniformity'] = f

probe = load('morphology-probe.json')
if probe:
    pr = {}
    for region in ['XE035', 'XE036', 'XO036']:
        for off in (0, 6):
            pr[f'{region}_angle+{off}'] = {str(m): dict(passes=sum(r['screen_pass'] for r in probe if r['region'] == region and r['angle_offset_deg'] == off and r['length_mult'] == m),
                                                        sites=sum(1 for r in probe if r['region'] == region and r['angle_offset_deg'] == off and r['length_mult'] == m),
                                                        median_fit_length_ratio=float(np.median([r['fit_length_ratio'] for r in probe if r['region'] == region and r['angle_offset_deg'] == off and r['length_mult'] == m])))
                                           for m in (1.0, 1.3, 1.5, 2.0, 4.0, 10.0)}
    out['over_length_probe_margin20'] = pr

mg = {}
ctrl = load('margin-control.json')
if ctrl:
    mg['control_reproduces_pipeline'] = ctrl['control_reproduces_pipeline']
    mg['positives_wide'] = ctrl['wide']
for region in ['XE035', 'XE036', 'XO036']:
    x = load(f'margin-{region}.json')
    if not x:
        continue
    rec = {}
    for lvl in sorted({r['true_snr'] for r in x['inject']}):
        rr = [r for r in x['inject'] if r['true_snr'] == lvl]
        rec[str(lvl)] = dict(wide=float(np.mean([r['screen_pass'] for r in rr])), baseline=float(np.mean([r['baseline_pass'] for r in rr])))
    mg[region] = dict(margin=x['margin'], null_passes=sum(r['screen_pass'] for r in x['null']), null_trials=len(x['null']), recall=rec,
                      probe={str(m): dict(passes=sum(r['screen_pass'] for r in x['probe'] if r['length_mult'] == m),
                                          median_fit_length_ratio=float(np.median([r['length_ratio'] for r in x['probe'] if r['length_mult'] == m])))
                             for m in (1.0, 1.3, 1.5, 2.0, 4.0, 10.0)})
out['wide_margin'] = mg
ld = load('length-diag.json')
if ld:
    out['optimizer_stall_XE036'] = ld
(D / 'experiments-summary.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=1)[:6000])
