from acquire_case import ROOT,CASES,get
from datetime import datetime,timedelta
from astropy.io import fits
import json
records=[]
for c,region in zip(CASES,['XE740','XE741']):
 records.append({**c,'archive_file':f'data/raw/{region}.fits','identity_status':'provisional-alias-visually-matching-field','published_region':c['published_plate'],'published_plate_id':None,'observer_image':f"observer-{c['date'].replace('-','')}.gif",'observer_station':'675'})
probes={c['date']:c for c in json.loads((ROOT/'stsci-probes.json').read_text())}
for c in json.loads((ROOT/'additional-coverage.json').read_text()):
 record={k:v for k,v in c.items() if k not in ['coverage','region_filename_matches']}
 record.update(id=c['date'],exposure_minutes=60,observer_station='413',observer_image=f"observer-{c['date'].replace('-','')}.gif")
 if c['region_filename_matches']:
  record.update(archive_file=f"data/raw/{c['published_region']}.fits",identity_status='exact-published-plate-id')
 else:
  p=probes[c['date']]
  record.update(archive_file=p.get('file'),identity_status='exact-published-plate-id' if p.get('plate_id_matches') else 'wrong-default-exposure')
 records.append(record)
for r in records:
 h=fits.getheader(ROOT/r['archive_file'])
 r['archive_region']=h.get('REGION');r['archive_plate_id']=h.get('PLATEID');r['archive_date_obs']=h.get('DATE-OBS');r['archive_exposure_minutes']=h.get('EXPOSURE')
 t=datetime.fromisoformat(r['date'])+timedelta(days=r['published_fractional_day']%1)
 r['published_calendar_time']=t.isoformat();r['published_time_scale_label']='UTC';r['working_time_scale']='UT1' if t.year<1962 else 'UTC'
 r['header_minus_published_seconds_uninterpreted']=(datetime.fromisoformat(h['DATE-OBS'])-t).total_seconds()
 r['scientifically_accepted']=False;r['image_measured']=False
 if r['identity_status']=='exact-published-plate-id':
  assert h['PLATEID'].strip()==r['published_plate_id']
  assert h['EXPOSURE']==r['exposure_minutes']
  get(r['id']+'-horizons.txt','https://ssd.jpl.nasa.gov/api/horizons.api',dict(format='text',COMMAND="'307261;'",CENTER="'413'",MAKE_EPHEM='YES',EPHEM_TYPE='OBSERVER',START_TIME="'"+t.strftime('%Y-%m-%d %H:%M:%S')+"'",STOP_TIME="'"+(t+timedelta(minutes=1)).strftime('%Y-%m-%d %H:%M:%S')+"'",STEP_SIZE="'1 m'",QUANTITIES="'1,3,9,36,37'",CSV_FORMAT='YES',ANG_FORMAT='DEG',EXTRA_PREC='YES',TIME_DIGITS='SECONDS'))
  get(r['observer_image'],'https://andrew-lowe.ca/2002ms4-'+r['date'].replace('-','')+'.gif')
 print(r['date'],r['identity_status'],r['archive_region'],r['archive_plate_id'],r['header_minus_published_seconds_uninterpreted'],flush=True)
(ROOT/'exposure-inventory.json').write_text(json.dumps(records,indent=2)+'\n')
