"""Write manifest.json: SHA-256 of the frozen inputs and code this package used, of every file it ships,
and the headline numbers check_package.py re-reads. Run from the repository root after the outputs exist."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = ['data/raw/XE036.fits', 'data/raw/XO036.fits', 'data/raw/irsa-XE035.fits', 'data/raw/gaia-dr3.csv',
          'data/raw/horizons-an10-1955.txt', 'data/controls.json', 'data/derived/measurements.json',
          'src/precovery/pipeline.py', 'src/precovery/verify.py', 'src/precovery/acquire.py']
SKIP = {'manifest.json', 'check.log'}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


noise = json.loads((HERE / 'noise-model.json').read_text())
summ = json.loads((HERE / 'sample-summary.json').read_text())
key = 'L15_tapered'
manifest = dict(
    schema='research-package-manifest/v1',
    base_commit='9795ae8e87ba525d57b7cb6a135abb8a9eb34c56',
    inputs=[dict(path=p, sha256=sha(p)) for p in INPUTS],
    outputs=[dict(path=f.name, sha256=sha(f)) for f in sorted(HERE.iterdir())
             if f.is_file() and f.name not in SKIP and f.suffix in ('.py', '.json', '.md', '.png')],
    headline=dict(pricing_key=key,
                  inflation_white={r: round(noise[r]['pricing'][key]['inflation_white'], 4) for r in noise},
                  best_case_whitening_gain={r: round(noise[r]['pricing'][key]['best_case_whitening_gain'], 4) for r in noise},
                  null_screen_passes={r: summ[r]['null']['screen_pass_with_stationary_gate'] for r in summ}))
(HERE / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print('manifest:', len(manifest['inputs']), 'inputs,', len(manifest['outputs']), 'outputs')
