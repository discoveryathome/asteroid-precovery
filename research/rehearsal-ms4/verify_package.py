"""Offline integrity, exposure-identity and cutout-coordinate checks; no detection verdict."""
from pathlib import Path
import hashlib,json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
R=Path(__file__).resolve().parent
m=json.loads((R/'manifest.json').read_text())
for row in m['files']:
 p=(R/row['path']).resolve();assert p.is_relative_to(R)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],p
records=json.loads((R/'exposure-inventory.json').read_text())
assert len(records)==7
counts={}
for c in records:
 counts[c['identity_status']]=counts.get(c['identity_status'],0)+1
 h=fits.getheader(R/c['archive_file'])
 assert h['REGION'].strip()==c['archive_region'].strip()
 assert h['PLATEID'].strip()==c['archive_plate_id'].strip()
 assert h['DATE-OBS']==c['archive_date_obs']
 assert not c['scientifically_accepted'] and not c['image_measured']
 if c['identity_status']=='exact-published-plate-id':
  assert h['PLATEID'].strip()==c['published_plate_id']
  assert h['REGION'].strip()==c['published_region']
  assert h['EXPOSURE']==c['exposure_minutes']
 elif c['identity_status']=='wrong-default-exposure':
  assert h['PLATEID'].strip()!=c['published_plate_id']
  assert not h['DATE-OBS'].startswith(c['date'])
assert counts=={'provisional-alias-visually-matching-field':2,'exact-published-plate-id':3,'wrong-default-exposure':2}
for region in ['XE740','XE741','ER729']:
 p=R/'data/raw'/region
 parent=fits.Header.fromstring(p.with_suffix('.original-header.txt').read_text(),sep='\n')
 cut=fits.getheader(p.with_suffix('.fits'))
 meta=json.loads(p.with_suffix('.provenance.json').read_text())
 x0,x1,y0,y1=meta['pixel_bounds_zero_based']
 assert fits.getdata(p.with_suffix('.fits')).shape==(y1-y0,x1-x0)
 xy=np.array([[0,0],[x1-x0-1,y1-y0-1],[(x1-x0-1)/2,(y1-y0-1)/2]])
 sky1=WCS(parent).all_pix2world(xy+np.array([x0,y0]),0)
 sky2=WCS(cut).all_pix2world(xy,0)
 assert np.max(np.abs(sky1-sky2))*3600<1e-5
review=json.loads((R/'review.json').read_text())
assert review['scientifically_accepted'] is False and review['submission_allowed'] is False
print(f"PASS: {len(m['files'])} file hashes; 7 epoch records; 3 exact plate IDs; 2 provisional aliases; 2 wrong-epoch matches rejected; native WCS preserved. Scientific review pending.")
