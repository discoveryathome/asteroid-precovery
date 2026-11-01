"""Stdlib-only check of this package: pinned input hashes, output hashes and headline numbers.

Run from the repository root: python3 research/correlated-noise-20261101/check_package.py
Exits non-zero on any mismatch. Does not rerun the reduction (see README for that).
"""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    m = json.loads((HERE / 'manifest.json').read_text())
    bad = []
    for row in m['inputs']:
        if sha256(ROOT / row['path']) != row['sha256']:
            bad.append(row['path'])
    # inputs must also agree with the project's own frozen manifest where it lists them
    frozen = {r['path']: r['sha256'] for r in json.loads((ROOT / 'data/manifest.json').read_text())['files']}
    for row in m['inputs']:
        if row['path'] in frozen and frozen[row['path']] != row['sha256']:
            bad.append(row['path'] + ' (differs from data/manifest.json)')
    for row in m['outputs']:
        if sha256(HERE / row['path']) != row['sha256']:
            bad.append(row['path'])
    if bad:
        sys.exit('FAIL: hash mismatch: ' + ', '.join(bad))
    noise = json.loads((HERE / 'noise-model.json').read_text())
    key = m['headline']['pricing_key']
    for region, want in m['headline']['inflation_white'].items():
        got = noise[region]['pricing'][key]['inflation_white']
        if abs(got - want) > 5e-4:
            sys.exit(f'FAIL: {region} inflation {got} != {want}')
    for region, want in m['headline']['best_case_whitening_gain'].items():
        got = noise[region]['pricing'][key]['best_case_whitening_gain']
        if abs(got - want) > 5e-4:
            sys.exit(f'FAIL: {region} whitening gain {got} != {want}')
    summ = json.loads((HERE / 'sample-summary.json').read_text())
    for region, want in m['headline']['null_screen_passes'].items():
        got = summ[region]['null']['screen_pass_with_stationary_gate']
        if got != want:
            sys.exit(f'FAIL: {region} null screen passes {got} != {want}')
    print(f"PASS: {len(m['inputs'])} input and {len(m['outputs'])} output hashes; headline numbers match")


if __name__ == '__main__':
    main()
