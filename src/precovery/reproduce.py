"""One-command offline certificate check and full pixel reduction."""
import json
import re
import sys
from pathlib import Path
from astropy.io import fits
import numpy as np
from .acquire import ROOT
from .pipeline import run
from .verify import check_manifest


def validate_inputs():
    check_manifest(ROOT,json.loads((ROOT/'data/manifest.json').read_text()))
    c=json.loads((ROOT/'data/controls.json').read_text())
    source=(ROOT/c['time_override']['snapshot']).read_text()
    pattern=r'1955 01 (26\.\d+) (\d+) (\d+) ([\d.]+) \+(\d+) (\d+) ([\d.]+)\s+675'
    rows=re.findall(pattern,source)
    if len(rows)!=4: raise ValueError('Published control observations missing')
    for i,row in enumerate(rows):
        day,hh,mm,ss,dd,dm,ds=map(float,row)
        ra=(hh+mm/60+ss/3600)*15;dec=dd+dm/60+ds/3600
        ctrl=c['positive_controls'][i//2]
        if not np.allclose([ra,dec],ctrl['reference_endpoints_deg'][i%2],rtol=0,atol=1e-10):
            raise ValueError('Control coordinate differs from original circular')
        if abs(day-ctrl['reference_fractional_days'][i%2])>1e-8: raise ValueError('Control time differs from circular')
        from datetime import datetime
        t=datetime.fromisoformat(ctrl[['start','end'][i%2]])
        fraction=26+(t.hour*3600+t.minute*60+t.second)/86400
        if abs(fraction-day)*86400>.5: raise ValueError('Time override does not match rounded published endpoints')
    for region in ['XE036','XO036']:
        h=fits.getheader(ROOT/f'data/raw/{region}.fits')
        if h['REGION'].strip()!=region: raise ValueError('Wrong plate region')
        # Actual plate ID is recorded in the version lock; never guess it from region.
        lock=json.loads((ROOT/'data/versions.json').read_text())
        if h['PLATEID'].strip()!=lock['images'][region]['plate_id']:raise ValueError('Wrong plate ID')
    return c


def main():
    validate_inputs()
    metrics=run(Path(sys.argv[1]) if len(sys.argv)>1 else None)
    if metrics['positives_passing']!=2 or metrics['negatives_rejected']!=4 or not metrics['pair_motion_consistent']:
        raise SystemExit('FAIL: control or motion regression')
    if metrics['scientific_acceptance'] or metrics['mpc_submission']:
        raise SystemExit('FAIL: publication gate violated')
    print('PASS: 2 documented exposures reproduced; 4 negative controls rejected; review pending; no submission.')


if __name__=='__main__': main()
