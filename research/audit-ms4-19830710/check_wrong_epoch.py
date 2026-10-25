"""Offline, standard-library audit of the MS4 1983-07-10 wrong-default STScI response.

Reads the frozen rehearsal-ms4 inputs in place and changes nothing. Explains one recorded
wrong-epoch response; it makes no detection, astrometric or scientific-acceptance claim.
"""
from pathlib import Path
import hashlib,json,math,re

PKG=Path(__file__).resolve().parents[1]/'rehearsal-ms4'
RECORD='1983-07-10'
EXACT='1983-06-14'  # exact published plate-ID match for S728/00TR
PROBES=['1982-08-09','1983-06-14','1983-07-01','1983-07-10']  # every STScI default request in the package
# Input bytes at base 9795ae8; also checked against the package's manifest and provenance records.
PINNED={
 'data/raw/1982-08-09-stsci-blue.fits':'4922e0b224b72e4e7276245285e793b893cead4f8515f4688e1ed945d07adee7',
 'data/raw/1983-06-14-stsci-blue.fits':'95443c45ad3c73ed53754b67b1ded7bddbdff0c9b9d00e33ab0f62c164bbc995',
 'data/raw/1983-07-01-stsci-blue.fits':'426fa03209f73c469f99459316bf2469b84f1e95a78bf6ee7cbaa258d7b4de37',
 'data/raw/1983-07-10-stsci-blue.fits':'3cfe5e717f7606f53de77e52fa79e91cc5a9cf09f6ff594a28487548bf31fe15',
 'data/raw/1982-08-09-stsci-blue.fits.provenance.json':'f05de065160710ad2353c12b876df84a7bd80421f942298eb013e411c3933cc8',
 'data/raw/1983-06-14-stsci-blue.fits.provenance.json':'83f307b846af99ef6c98e0454863203b3adc03dba09c709f8970f2b822e755a4',
 'data/raw/1983-07-01-stsci-blue.fits.provenance.json':'4eb2b474c570431084119ef386b7c12b0b7220f6477e061cc8b4170513a28f89',
 'data/raw/1983-07-10-stsci-blue.fits.provenance.json':'69b84fb9e7ae0a173b253453a664e6cae02ae4542f5e9f65f9be9eb3980e3e68',
 'data/raw/1983-06-14-horizons.txt':'7c511e822f9cc5ffcd72b8bb382f6cc31eeeab493e907584ce618ff91313c59e',
 'additional-coverage.json':'ed9a628c0e02e099889fa83b72c4bb7454da114e7c307f78049cc186c5e8c6e7',
 'exposure-inventory.json':'a5e5c368c2468b875eccd16217475e80f0dae1d239883e68f719b39ad4da8c5f',
 'stsci-probes.json':'75148ff2cec7639339c69193c26a38f648c3553fbeb3847a241c876abb7b673e',
}
# Keywords that describe the cutout rather than the plate; every other card must be identical.
CUTOUT_KEYS={'DATE','CNPIX1','CNPIX2','CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2',
             'DATAMIN','DATAMAX','OBJCTRA','OBJCTDEC','OBJCTX','OBJCTY'}
IDENTITY_KEYS={'SURVEY','REGION','PLATEID','PLTLABEL','DATE-OBS','EXPOSURE','PLATERA','PLATEDEC'}

def sha(b):return hashlib.sha256(b).hexdigest()

def load(rel,manifest):
 b=(PKG/rel).read_bytes();h=sha(b)
 assert h==PINNED[rel],f'{rel}: bytes differ from pinned base'
 assert manifest[rel]==h,f'{rel}: manifest disagrees'
 return b

def cards(b):
 """Primary-header cards as {keyword: raw 80-character card}, stopping at END."""
 out={}
 for i in range(0,len(b),80):
  c=b[i:i+80].decode('ascii');k=c[:8].strip()
  if k=='END':return out
  if k:
   assert k not in out,f'duplicate keyword {k}'
   out[k]=c
 raise ValueError('no END card')

def value(card):
 v=card[10:]
 if v.lstrip().startswith("'"):return re.match(r"\s*'((?:[^']|'')*)'",v).group(1).replace("''","'").rstrip()
 v=v.split('/')[0].strip()
 return {'T':True,'F':False}.get(v) if v in 'TF' else float(v)

def sep_deg(ra1,de1,ra2,de2):
 r=math.radians;a=math.sin(r(de2-de1)/2)**2+math.cos(r(de1))*math.cos(r(de2))*math.sin(r(ra2-ra1)/2)**2
 return math.degrees(2*math.asin(math.sqrt(a)))

def tangent(ra,de,ra0,de0):
 """Gnomonic offsets (deg) of (ra,de) about a field centre."""
 r=math.radians;d=r(ra-ra0);c=math.sin(r(de0))*math.sin(r(de))+math.cos(r(de0))*math.cos(r(de))*math.cos(d)
 x=math.cos(r(de))*math.sin(d)/c;y=(math.cos(r(de0))*math.sin(r(de))-math.sin(r(de0))*math.cos(r(de))*math.cos(d))/c
 return math.degrees(x),math.degrees(y)

def classify(req,h):
 """Accept a response only if region, plate ID and UT date all match the published exposure."""
 ok=(value(h['REGION'])==req['published_region'] and value(h['PLATEID'])==req['published_plate_id']
     and value(h['DATE-OBS']).startswith(req['date']))
 return 'exact' if ok else 'wrong-exposure'

def requested_position(prov):
 url=prov['url']
 q=dict(p.split('=',1) for p in url.split('?',1)[1].split('&'))
 assert q['v']=='poss2ukstu_blue' and q['e']=='J2000'
 return float(q['r']),float(q['d']),float(q['h'])

def main():
 m={f['path']:f['sha256'] for f in json.loads((PKG/'manifest.json').read_text())['files']}
 raw={rel:load(rel,m) for rel in PINNED}
 prov={d:json.loads(raw[f'data/raw/{d}-stsci-blue.fits.provenance.json']) for d in PROBES}
 for d in PROBES:  # the provenance record must describe the same bytes
  assert prov[d]['sha256']==PINNED[f'data/raw/{d}-stsci-blue.fits']
 inv={r['date']:r for r in json.loads(raw['exposure-inventory.json'])}
 hdr={d:cards(raw[f'data/raw/{d}-stsci-blue.fits']) for d in PROBES}
 req=inv[RECORD];got=hdr[RECORD];ref=hdr[EXACT]
 out={'record':RECORD,'requested':{k:req[k] for k in ('published_region','published_plate_id','published_calendar_time')},
      'returned':{k:value(got[k]) for k in ('REGION','PLATEID','PLTLABEL','DATE-OBS','SURVEY')}}

 # 1. Identity: the response is rejected, and the exact 1983-06-14 match is accepted.
 assert classify(req,got)=='wrong-exposure' and classify(inv[EXACT],ref)=='exact'

 # 2. Same physical plate: every plate-level card is byte-identical to the exact 06-14 match.
 assert set(got)==set(ref)
 plate=sorted(set(got)-CUTOUT_KEYS)
 diff=[k for k in plate if got[k]!=ref[k]]
 assert IDENTITY_KEYS<=set(plate) and not diff,diff
 assert [k for k in CUTOUT_KEYS-{'DATE'} if got[k]!=ref[k]],'cutouts should differ in position'
 out['same_plate_as_exact_match']={'plate_cards_compared':len(plate),'differing':diff}

 # 3. Selection: the default plate for each request is the UK Schmidt field whose centre is nearest.
 centres={}
 for rec in json.loads(raw['additional-coverage.json']):
  for row in rec['coverage']:
   if row['facility_name']=='UK Schmidt':
    c=(float(row['crval1']),float(row['crval2']))
    assert centres.setdefault(int(row['tile']),c)==c
 calib={}
 for d in PROBES:  # blue field centres seen in headers versus the red-survey centres from IRSA
  t=int(value(hdr[d]['REGION'])[1:]);calib[t]=round(sep_deg(value(hdr[d]['PLATERA']),value(hdr[d]['PLATEDEC']),*centres[t])*60,2)
 assert set(calib)=={727,728} and max(calib.values())<1.0,calib
 half=value(got['PLTSIZEX'])*value(got['PLTSCALE'])/3600/2  # nominal plate half-width, deg
 rows=[]
 for d in PROBES:
  ra,de,_=requested_position(prov[d])
  dist={t:sep_deg(ra,de,*c) for t,c in centres.items()}
  edge={t:half-max(map(abs,tangent(ra,de,*c))) for t,c in centres.items()}
  near=min(dist,key=dist.get);inner=max(edge,key=edge.get)
  ret=int(value(hdr[d]['REGION'])[1:]);pub=int(inv[d]['published_region'][1:])
  rows.append({'probe':d,'returned':ret,'published':pub,'nearest_centre':near,'largest_edge_margin':inner,
               'deg_to_centre':{t:round(v,3) for t,v in sorted(dist.items())},
               'edge_margin_deg':{t:round(v,3) for t,v in sorted(edge.items())}})
  assert near==ret and inner==ret,(d,near,inner,ret)
  assert edge[pub]>0  # the published field does cover the requested position
 out['selection']={'blue_vs_red_centre_arcmin':calib,'probes':rows}

 # 4. Exclusion without header IDs: where the object was on this plate at its exposure epoch.
 scale=value(got['XPIXELSZ'])/1000*value(got['PLTSCALE'])  # arcsec per scan pixel
 x0,y0=value(got['CNPIX1']),value(got['CNPIX2']);nx,ny=value(got['NAXIS1']),value(got['NAXIS2'])
 ox,oy=value(ref['OBJCTX']),value(ref['OBJCTY'])  # STScI's own plate pixel for the published 06-14 position
 def outside(x,y):return math.hypot(max(x0-x,0,x-(x0+nx)),max(y0-y,0,y-(y0+ny)))
 gap=outside(ox,oy)
 assert gap>0
 # Cross-check with the frozen Horizons row (its orbit already includes the 1983 observations).
 hz=[l for l in raw['data/raw/1983-06-14-horizons.txt'].decode().split('$$SOE')[1].split('$$EOE')[0].splitlines() if l.strip()][0].split(',')
 assert hz[0].strip()=='1983-Jun-14 13:03:45'
 hra,hde=float(hz[3]),float(hz[4])
 hsep=sep_deg(hra,hde,value(got['CRVAL1']),value(got['CRVAL2']))*3600
 halfdiag=math.hypot(nx,ny)/2*scale
 assert hsep>halfdiag
 out['exclusion']={'returned_epoch':value(got['DATE-OBS']),'cutout_box_px':[x0,x0+nx,y0,y0+ny],
                   'object_px_at_returned_epoch':[ox,oy],'gap_outside_box_px':round(gap,1),
                   'gap_outside_box_arcmin':round(gap*scale/60,2),
                   'horizons_to_cutout_centre_arcmin':round(hsep/60,2),'cutout_half_diagonal_arcmin':round(halfdiag/60,2)}

 # Failure controls: each check must be able to fail.
 bad=bytearray(raw['data/raw/1983-07-10-stsci-blue.fits']);bad[-1]^=1
 assert sha(bytes(bad))!=PINNED['data/raw/1983-07-10-stsci-blue.fits']
 patched=dict(got)
 for k,v in (('REGION','S655'),('PLATEID','03YK'),('DATE-OBS','1983-07-10T11:46')):
  patched[k]=f"{k:<8}= '{v:<8}'".ljust(80)
 assert classify(req,patched)=='exact'
 assert [k for k in plate if patched[k]!=ref[k]]
 assert outside(x0+nx/2,y0+ny/2)==0
 out['controls']='flipped byte fails hash; S655/03YK/1983-07-10 header accepted; patched header breaks same-plate; object inside box gives zero gap'

 print(json.dumps(out,indent=1))
 print(f"PASS: {RECORD} request for {req['published_region']}/{req['published_plate_id']} returned "
       f"{out['returned']['REGION']}/{out['returned']['PLATEID']} ({out['returned']['PLTLABEL']}, {out['returned']['DATE-OBS']}); "
       f"same plate as the exact {EXACT} match ({len(plate)} non-cutout header cards identical); nearest-centre and edge-margin "
       f"selection both reproduce all {len(PROBES)} STScI defaults; object {out['exclusion']['gap_outside_box_arcmin']}' outside "
       f"the returned cutout at the returned epoch. No detection or scientific-acceptance verdict.")

if __name__=='__main__':main()
