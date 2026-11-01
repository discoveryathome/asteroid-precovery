"""Summarize the null and injection samples into one JSON and a printed table."""
import json
import sys
from pathlib import Path
import numpy as np

OUT = Path(sys.argv[1])
NOISE = json.loads((OUT / 'noise-model.json').read_text())


def robust_std(x):
    x = np.asarray(x, float)
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def pct(x, q):
    return [round(float(v), 2) for v in np.percentile(x, q)]


summary = {}
for region in ['XE035', 'XE036', 'XO036']:
    s = {}
    null = json.loads((OUT / f'sample-{region}-null.json').read_text())
    rows = null['rows']
    snr = np.array([r['snr'] for r in rows])
    dbic = np.array([r['dbic'] for r in rows])
    morph = [r for r in rows if r['morphology_pass']]
    full = [r for r in rows if r['screen_pass']]
    z = np.array([r['a_hat'] / r['a_white_sigma'] for r in rows])
    noise_ratio = np.median([r['noise'] for r in rows]) / NOISE[region]['sigma_std']
    s['null'] = dict(
        trials=len(rows), rotations=null['rotations'], step_pixel=null['step'],
        snr_p50_p90_p99_max=pct(snr, [50, 90, 99, 100]), dbic_p50_p90_p99_max=pct(dbic, [50, 90, 99, 100]),
        snr_ge_5=int(np.sum(snr >= 5)), dbic_ge_10=int(np.sum(dbic >= 10)),
        morphology_pass=len(morph), morphology_pass_rate=len(morph) / len(rows),
        screen_pass_with_stationary_gate=len(full),
        screen_pass_on_any_gaia_within_4arcsec=sum(r['near_gaia_any_arcsec'] < 4 for r in full),
        screen_pass_on_gaia_outside_gate_mag_range=sum(r['near_gaia_any_arcsec'] < 4 and not 12.5 < r['near_gaia_any_g'] < 18.5 for r in full),
        morphology_pass_rejected_by_stationary_gate=sum(r['stationary_match'] if 'stationary_match' in r else ('stationary_match' in r['reasons']) for r in morph),
        passes=[dict(rot=r['rot_deg'], center=[round(c, 1) for c in r['center']], fit_xy=[round(c, 1) for c in r['fit_xy']],
                     snr=round(r['snr'], 2), dbic=round(r['dbic'], 1), length_ratio=round(r['length_ratio'], 3),
                     angle_error_deg=round(r['angle_error_deg'], 2), near_cat_arcsec=round(r['near_cat_arcsec'], 2),
                     near_gaia_any_arcsec=round(r['near_gaia_any_arcsec'], 2), near_gaia_any_g=round(r['near_gaia_any_g'], 2),
                     screen_pass=r['screen_pass'], reasons=r['reasons']) for r in morph],
        linear_filter_z_robust_std=robust_std(z),
        median_window_noise_over_sigma_std=float(noise_ratio),
        predicted_z_std_from_autocorrelation=float(NOISE[region]['pricing'][null['F_key']]['inflation_white'] / noise_ratio),
    )
    pos = {'XE036': 'an10-red-1955', 'XO036': 'an10-blue-1955'}.get(region)
    if pos:
        m = next(x for x in json.loads(Path('data/derived/measurements.json').read_text()) if x['id'] == pos)
        s['null']['positive'] = dict(id=pos, snr=m['snr'], dbic=m['delta_bic_vs_point'],
                                     null_full_passes_with_snr_ge=sum(r['snr'] >= m['snr'] for r in full),
                                     null_full_passes_with_dbic_ge=sum(r['dbic'] >= m['delta_bic_vs_point'] for r in full),
                                     null_full_pass_snr_max=max([r['snr'] for r in full], default=None),
                                     null_full_pass_dbic_max=max([r['dbic'] for r in full], default=None))
    s['null']['full_pass_by_rotation'] = {str(rot): [sum(r['screen_pass'] for r in rows if r['rot_deg'] == rot),
                                                     sum(r['rot_deg'] == rot for r in rows)] for rot in null['rotations']}
    inj = json.loads((OUT / f'sample-{region}-inject.json').read_text())
    rec = {}
    for level in sorted({r['injected_true_snr'] for r in inj['rows']}):
        rr = [r for r in inj['rows'] if r['injected_true_snr'] == level]
        fails = {}
        for r in rr:
            for k in r['reasons']:
                fails[k] = fails.get(k, 0) + 1
        rec[str(level)] = dict(n=len(rr), recall=sum(r['screen_pass'] for r in rr) / len(rr),
                               morphology_recall=sum(r['morphology_pass'] for r in rr) / len(rr),
                               median_pipeline_snr=float(np.median([r['snr'] for r in rr])),
                               median_dbic=float(np.median([r['dbic'] for r in rr])),
                               median_amp_over_sigma=float(np.median([r['injected_amp_over_sigma'] for r in rr])),
                               median_position_error_arcsec=float(np.median([r['residual_arcsec'] for r in rr])),
                               fail_reasons=dict(sorted(fails.items(), key=lambda kv: -kv[1])))
    s['inject'] = rec
    pipe = np.array([r['snr'] for r in inj['rows']])
    true = np.array([r['injected_true_snr'] for r in inj['rows']])
    s['inject_pipeline_snr_over_true_median'] = float(np.median(pipe / true))
    summary[region] = s
    n = s['null']
    print(f"{region}: null n={n['trials']} snr p50/90/99/max={n['snr_p50_p90_p99_max']} dbic={n['dbic_p50_p90_p99_max']}"
          f" snr>=5:{n['snr_ge_5']} dbic>=10:{n['dbic_ge_10']} morph-pass {n['morphology_pass']} full-pass {n['screen_pass_with_stationary_gate']}"
          f" (on any Gaia<4\": {n['screen_pass_on_any_gaia_within_4arcsec']}; outside gate mags: {n['screen_pass_on_gaia_outside_gate_mag_range']})")
    if 'positive' in n:
        print(f"   positive {n['positive']}")
    print(f"   full pass by rotation {n['full_pass_by_rotation']}")
    print(f"   linear z robust std {n['linear_filter_z_robust_std']:.3f} vs predicted {n['predicted_z_std_from_autocorrelation']:.3f}"
          f" (window noise/sigma {n['median_window_noise_over_sigma_std']:.3f}); pipeline SNR/true median {s['inject_pipeline_snr_over_true_median']:.3f}")
    for k, v in rec.items():
        print(f"   inject true SNR {k}: n={v['n']} recall={v['recall']:.2f} morph={v['morphology_recall']:.2f}"
              f" pipeSNR={v['median_pipeline_snr']:.2f} dBIC={v['median_dbic']:.1f} poserr={v['median_position_error_arcsec']:.2f}\" fails={v['fail_reasons']}")
(OUT / 'sample-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
