"""Inventory the other five published epochs; coverage is not identification."""
from acquire_case import get,ROOT
from astropy.table import Table
from xml.etree import ElementTree as E
import json
CASES=[
 ('1982-08-09',9.38380,(16,9,39.06),(-12,20,0.8),'S727','044C'),
 ('1983-06-14',14.54428,(16,16,21.79),(-12,16,28.0),'S728','00TR'),
 ('1983-07-01',1.46809,(16,15,17.21),(-12,14,19.9),'S656','02C3'),
 ('1983-07-10',10.49060,(16,14,47.94),(-12,13,44.0),'S655','03YK'),
 ('1992-06-23',23.57704,(16,53,16.75),(-11,13,38.2),'ER729','A2JL')]
out=[]
for day,fraction,ra,dec,region,plate in CASES:
 ra=(ra[0]+ra[1]/60+ra[2]/3600)*15;dec=dec[0]-dec[1]/60-dec[2]/3600
 x=get(day+'-atlas.xml','https://irsa.ipac.caltech.edu/cgi-bin/Atlas/nph-atlas',dict(mission='DSS',locstr=f'{ra} {dec} eq',regSize='0.08',covers='on',mode='PI'))
 tree=E.parse(x);assert tree.getroot().get('status')=='ok'
 table=get(day+'-images.tbl',tree.find('./images/metadata').text.strip())
 t=Table.read(table,format='ascii.ipac')
 rows=[{k:str(r[k]) for k in t.colnames} for r in t]
 record=dict(date=day,published_fractional_day=fraction,ra_deg=ra,dec_deg=dec,published_region=region,published_plate_id=plate,coverage=rows,region_filename_matches=[r['fname'] for r in rows if r['fname'].endswith('_'+region+'.fits')])
 out.append(record);print(day,region,len(rows),record['region_filename_matches'],flush=True)
 (ROOT/'additional-coverage.json').write_text(json.dumps(out,indent=2)+'\n')
