from acquire_case import get,ROOT
from astropy.io import fits
import json
out=[]
for c in json.loads((ROOT/'additional-coverage.json').read_text()):
 if c['region_filename_matches']:continue
 try:
  p=get(c['date']+'-stsci-blue.fits','https://archive.stsci.edu/cgi-bin/dss_search',dict(v='poss2ukstu_blue',r=str(c['ra_deg']),d=str(c['dec_deg']),e='J2000',h='4.25',w='4.25',f='fits',c='none'))
  h=fits.getheader(p)
  result={'date':c['date'],'expected_region':c['published_region'],'expected_plate_id':c['published_plate_id'],'region':h.get('REGION'),'plate_id':h.get('PLATEID'),'date_obs':h.get('DATE-OBS'),'exposure':h.get('EXPOSURE'),'file':str(p.relative_to(ROOT))}
  result['plate_id_matches']=str(h.get('PLATEID','')).strip()==c['published_plate_id']
 except Exception as e:result={'date':c['date'],'error':str(e)}
 out.append(result);print(result,flush=True)
 (ROOT/'stsci-probes.json').write_text(json.dumps(out,indent=2)+'\n')
