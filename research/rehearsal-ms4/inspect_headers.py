import io,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
from astropy.io import fits
from precovery.range_cutout import get_range
R=Path(__file__).resolve().parent
rows=json.loads((R/'coverage.json').read_text());out=[]
for name in sorted(set(x['fname'] for x in rows)):
 url='https://irsa.ipac.caltech.edu/data/DSS/'+name
 p=R/'data/raw'/(Path(name).stem+'.header.bin')
 if not p.exists():
  data,headers=get_range(url,0,28799);p.write_bytes(data)
  p.with_suffix('.provenance.json').write_text(json.dumps({'url':url,'retrieved_at':datetime.now(timezone.utc).isoformat(),'sha256':hashlib.sha256(data).hexdigest(),'range':[0,28799],'response_headers':headers},indent=2)+'\n')
 h=fits.Header.fromfile(io.BytesIO(p.read_bytes()))
 a={'archive_file':name,'headers':{k:h.get(k) for k in ['REGION','PLATEID','DATE-OBS','EXPOSURE','PLTLABEL','TELESCOP','NAXIS1','NAXIS2']}}
 out.append(a);print(a,flush=True)
(R/'header-audit.json').write_text(json.dumps(out,indent=2)+'\n')
