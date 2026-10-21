from acquire_case import get,ROOT,CASES
from datetime import datetime,timedelta
from astropy.time import Time
import json
for c in CASES:
 date=datetime.fromisoformat(c['date'])+timedelta(days=c['published_fractional_day']%1)
 # Horizons calendar dates before 1962 are UT1; label discrepancy is preserved.
 get(c['id']+'-horizons.txt','https://ssd.jpl.nasa.gov/api/horizons.api',dict(format='text',COMMAND="'307261;'",CENTER="'675'",MAKE_EPHEM='YES',EPHEM_TYPE='OBSERVER',START_TIME="'"+date.strftime('%Y-%m-%d %H:%M:%S')+"'",STOP_TIME="'"+(date+timedelta(minutes=1)).strftime('%Y-%m-%d %H:%M:%S')+"'",STEP_SIZE="'1 m'",QUANTITIES="'1,3,9,36,37'",CSV_FORMAT='YES',ANG_FORMAT='DEG',EXTRA_PREC='YES',TIME_DIGITS='SECONDS'))
for date in ['19540408','19550416']:
 get('observer-'+date+'.gif','https://andrew-lowe.ca/2002ms4-'+date+'.gif')
get('sbdb-ms4.json','https://ssd-api.jpl.nasa.gov/sbdb.api',{'sstr':'307261','cov':'mat','full-prec':'true','discovery':'true'})
print('Frozen predictions, covariance and observer comparison images saved.')
