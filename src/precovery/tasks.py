"""Cheap offline orientation checks for the four Discovery@Home tracks."""
import json
import sys
from astropy.table import Table
from .acquire import ROOT
from .verify import check_manifest


def main(track):
    check_manifest(ROOT,json.loads((ROOT/'data/manifest.json').read_text()))
    if track=='archive':
        table=Table.read(ROOT/'research/sources/irsa-images.tbl',format='ascii.ipac')
        files=[str(r['fname']) for r in table if 'XE036' in str(r['fname']) or 'XO036' in str(r['fname'])]
        assert len(files)==2
        print(json.dumps({'coverage_rows':len(table),'documented_plate_matches':files},indent=2))
    elif track=='astrometry':
        c=json.loads((ROOT/'data/derived/calibration.json').read_text())
        print(json.dumps({k:{x:v[x] for x in ['inliers','loo_rms_per_axis_arcsec']} for k,v in c.items()},indent=2))
        print('Inspect data/controls.json time_override and data/derived/reference-comparison.json')
    elif track=='pipeline':
        from .reproduce import validate_inputs
        validate_inputs()
        print('PASS: immutable inputs, published endpoint transcription, image IDs and time overrides verified')
    elif track=='review':
        r=json.loads((ROOT/'review/record.json').read_text())
        assert not r['scientifically_accepted'] and not r['submission_allowed'] and r['reviewer'] is None
        print(json.dumps({'status':r['status'],'open_obligations':len(r['open_obligations']),'reviewer':r['reviewer']},indent=2))
    else: raise ValueError(track)


if __name__=='__main__':main(sys.argv[1])
