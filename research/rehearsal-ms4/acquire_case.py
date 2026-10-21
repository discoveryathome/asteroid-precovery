"""Bounded acquisition of a documented second control; no scientific acceptance."""
import argparse, json, hashlib
from pathlib import Path
from datetime import datetime, timezone
from xml.etree import ElementTree as E
from precovery import acquire
from astropy.table import Table
ROOT=Path(__file__).resolve().parent
acquire.ROOT=ROOT
CASES=[
 {'id':'1954','date':'1954-04-08','published_fractional_day':8.37535,'ra_deg':(14+15/60+39.35/3600)*15,'dec_deg':-(13+39/60+52.5/3600),'published_plate':'se0741','exposure_minutes':45},
 {'id':'1955','date':'1955-04-16','published_fractional_day':16.32917,'ra_deg':(14+19/60+36.64/3600)*15,'dec_deg':-(13+38/60+22.7/3600),'published_plate':'se0742','exposure_minutes':50}
]
def get(name,url,params=None):
 p=ROOT/'data/raw'/name
 if p.exists():
  m=json.loads(p.with_suffix(p.suffix+'.provenance.json').read_text())
  assert hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256']
  return p
 return acquire.fetch(name,url,params,limit=10_000_000)
def run(stage):
 if stage=='sources':
  get('observer-report.html','https://andrew-lowe.ca/2002ms4.htm')
  get('case-index.html','https://andrew-lowe.ca/tno-precoveries.htm')
  get('discovery-mpec.html','https://www.minorplanetcenter.net/iau/mpec/K02/K02W27.html')
 elif stage=='coverage':
  rows=[]
  for c in CASES:
   x=get(c['id']+'-atlas.xml','https://irsa.ipac.caltech.edu/cgi-bin/Atlas/nph-atlas',dict(mission='DSS',locstr=f"{c['ra_deg']} {c['dec_deg']} eq",regSize='0.08',covers='on',mode='PI'))
   tree=E.parse(x); assert tree.getroot().get('status')=='ok'
   table=get(c['id']+'-images.tbl',tree.find('./images/metadata').text.strip())
   t=Table.read(table,format='ascii.ipac')
   print(c['id'],len(t),'coverage rows',t.colnames,flush=True)
   for r in t:
    row={k:str(r[k]) for k in t.colnames};rows.append({'case':c['id'],**row});print(row,flush=True)
  (ROOT/'coverage.json').write_text(json.dumps(rows,indent=2)+'\n')
 elif stage=='cutouts':
  for c in CASES:
   x=get(c['id']+'-cutouts.xml','https://irsa.ipac.caltech.edu/cgi-bin/Cutouts/nph-cutouts',dict(mission='DSS',locstr=f"{c['ra_deg']} {c['dec_deg']} eq",sizeX='255',units='arcsec',mode='PI',ntable_cutouts='3',cutouttbl1='images_red',cutouttbl2='images_blue',cutouttbl3='images_ir',outtype='single'))
   print(x.read_text(),flush=True)
 if not (ROOT/'cases.json').exists():
  (ROOT/'cases.json').write_text(json.dumps({'object':'307261 (2002 MS4)','source':'https://andrew-lowe.ca/2002ms4.htm','cases':CASES,'timing_warning':'Observer labels dates UTC; historical UT1 convention and exact header timing require independent review.','scientifically_accepted':False},indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['sources','coverage','cutouts']);run(p.parse_args().stage)
