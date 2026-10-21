"""Generalize the existing bounded native DSS range reader to an explicit archive URL."""
import io,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from precovery.range_cutout import get_range
from acquire_case import ROOT

def retrieve(url,region,ra,dec,size_arcsec=255):
 prefix=ROOT/'data/raw'/region
 if prefix.with_suffix('.fits').exists():raise FileExistsError(prefix)
 raw,hh=get_range(url,0,28799)
 h=fits.Header.fromfile(io.BytesIO(raw))
 assert h['BITPIX']==16 and h['REGION'].strip()==region
 w=WCS(h);x,y=w.all_world2pix([[ra,dec]],0)[0]
 scale=abs(h['XPIXELSZ']*h['PLTSCALE']/1000)
 half=int(np.ceil(size_arcsec/scale/2));x0,x1=int(x)-half,int(x)+half+1;y0,y1=int(y)-half,int(y)+half+1
 nx,ny=h['NAXIS1'],h['NAXIS2'];assert min(x0,y0)>=0 and x1<=nx and y1<=ny
 headlen=len(h.tostring().encode('ascii'));start,end=headlen+y0*nx*2,headlen+y1*nx*2-1
 assert end-start+1<40_000_000
 rawrows,rh=get_range(url,start,end)
 data=np.frombuffer(rawrows,dtype='>i2').reshape(y1-y0,nx)[:,x0:x1].copy()
 out=h.copy();out['CNPIX1']=h.get('CNPIX1',0)+x0;out['CNPIX2']=h.get('CNPIX2',0)+y0
 fits.PrimaryHDU(data,header=out).writeto(prefix.with_suffix('.fits'))
 prefix.with_suffix('.original-header.txt').write_text(h.tostring(sep='\n'))
 meta=dict(url=url,retrieved_at=datetime.now(timezone.utc).isoformat(),parent_bytes=int(rh['Content-Range'].split('/')[-1]),parent_checksum=None,parent_checksum_reason='Full plate intentionally not downloaded',header_range=[0,28799],header_range_sha256=hashlib.sha256(raw).hexdigest(),pixel_range=[start,end],pixel_range_sha256=hashlib.sha256(rawrows).hexdigest(),pixel_bounds_zero_based=[x0,x1,y0,y1],cutout_size_arcmin=size_arcsec/60,native_pixels=True,etag=rh.get('ETag'),last_modified=rh.get('Last-Modified'),sha256=hashlib.sha256(prefix.with_suffix('.fits').read_bytes()).hexdigest())
 prefix.with_suffix('.provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
 print(region,data.shape,h['DATE-OBS'],h['PLATEID'],h['EXPOSURE'],flush=True)
 return h
if __name__=='__main__':
 records=json.loads((ROOT/'additional-coverage.json').read_text())
 for c in records:
  for name in c['region_filename_matches']:
   retrieve('https://irsa.ipac.caltech.edu/data/DSS/'+name,c['published_region'],c['ra_deg'],c['dec_deg'],size_arcsec=150)
