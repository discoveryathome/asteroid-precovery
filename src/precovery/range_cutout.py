"""Retrieve bounded FITS byte ranges; preserve native pixels and DSS plate solution."""
import io
import json
import hashlib
from datetime import datetime, timezone
import requests
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from .acquire import ROOT


def get_range(url, start, end):
    with requests.get(url, headers={'Range': f'bytes={start}-{end}'},
                      timeout=(15, 90), stream=True) as r:
        r.raise_for_status()
        if r.status_code != 206 or not r.headers.get('Content-Range', '').startswith(f'bytes {start}-{end}/'):
            raise ValueError('Server did not honor the byte range; refusing full image')
        body = bytearray()
        for chunk in r.iter_content(65536):
            body.extend(chunk)
            if len(body) > end-start+1: raise ValueError('Range length exceeded')
        if len(body) != end-start+1: raise ValueError('Truncated range')
        return bytes(body), dict(r.headers)


def retrieve(band, region, ra=129.5, dec=72.665, size_arcmin=10):
    url = f'https://irsa.ipac.caltech.edu/data/DSS/images/dss1{band}/dss1{band}_{region}.fits'
    prefix = ROOT / 'data/raw' / region
    if prefix.with_suffix('.fits').exists(): raise FileExistsError(prefix)
    raw_header, headers = get_range(url, 0, 28799)
    h = fits.Header.fromfile(io.BytesIO(raw_header))
    header_len = len(h.tostring().encode('ascii'))
    if h['BITPIX'] != 16: raise ValueError('Only signed 16-bit uncompressed DSS supported')
    w = WCS(h)
    x, y = w.all_world2pix([[ra, dec]], 0)[0]
    scale = abs(float(h['XPIXELSZ'])*float(h['PLTSCALE'])/1000)
    half = int(np.ceil(size_arcmin*60/scale/2))
    x0,x1 = int(x)-half,int(x)+half+1
    y0,y1 = int(y)-half,int(y)+half+1
    nx,ny=h['NAXIS1'],h['NAXIS2']
    if min(x0,y0)<0 or x1>nx or y1>ny: raise ValueError('Cutout crosses plate edge')
    start,end=header_len+y0*nx*2,header_len+y1*nx*2-1
    if end-start+1>40_000_000: raise ValueError('Range exceeds 40 MB budget')
    rows, row_headers = get_range(url,start,end)
    data = np.frombuffer(rows,dtype='>i2').reshape(y1-y0,nx)[:,x0:x1].copy()
    # CNPIX is the plate-coordinate origin used by the DSS polynomial WCS.
    # CD/CRPIX TAN is absent in these native DSS plates.
    out = h.copy(); out['CNPIX1']=h.get('CNPIX1',0)+x0;out['CNPIX2']=h.get('CNPIX2',0)+y0
    fits.PrimaryHDU(data,header=out).writeto(prefix.with_suffix('.fits'))
    prefix.with_suffix('.original-header.txt').write_text(h.tostring(sep='\n'))
    meta = dict(url=url,retrieved_at=datetime.now(timezone.utc).isoformat(),
                parent_bytes=int(row_headers['Content-Range'].split('/')[-1]),
                parent_checksum=None,parent_checksum_reason='Full plate intentionally not downloaded',
                header_range=[0,28799],header_range_sha256=hashlib.sha256(raw_header).hexdigest(),
                pixel_range=[start,end],pixel_range_sha256=hashlib.sha256(rows).hexdigest(),
                pixel_bounds_zero_based=[x0,x1,y0,y1],cutout_size_arcmin=size_arcmin,
                native_pixels=True,etag=row_headers.get('ETag'),last_modified=row_headers.get('Last-Modified'),
                sha256=hashlib.sha256(prefix.with_suffix('.fits').read_bytes()).hexdigest())
    prefix.with_suffix('.provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(region, data.shape, h.get('DATE-OBS'), h.get('EXPOSURE'), 'bytes transferred',len(rows)+len(raw_header),flush=True)


if __name__=='__main__':
    retrieve('red','XE036')
    retrieve('blue','XO036')
