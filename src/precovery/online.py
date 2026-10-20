"""Optional public-data refresh into a NEW directory; frozen baseline stays unchanged."""
import argparse
import json
from pathlib import Path
from xml.etree import ElementTree as E
from astropy.io import fits
from . import acquire, range_cutout


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--most-only',action='store_true');args=p.parse_args()
    out=args.output.resolve()
    if out.exists(): raise SystemExit('Output must be a new directory')
    out.mkdir(parents=True);acquire.ROOT=out;range_cutout.ROOT=out
    if args.most_only:
        acquire.fetch('most-an10-2010.tbl','https://irsa.ipac.caltech.edu/cgi-bin/MOST/nph-most',
            dict(catalog='wise_merge',input_type='name_input',obj_name='137108',obs_begin='2010 01 01',
                 obs_end='2010 02 01',output_mode='Brief',ephem_step='0.05'))
        print('MOST response saved; coverage only, not a detection or precovery.');return
    acquire.horizons('horizons-an10-1955.txt','1955-01-26 09:00','1955-01-26 11:00')
    acquire.fetch('sbdb-an10.json','https://ssd-api.jpl.nasa.gov/sbdb.api',
                  {'sstr':'137108','cov':'mat','full-prec':'true','discovery':'true'})
    # Exact public query from the baseline; each refreshed response has new provenance.
    baseline=Path(__file__).resolve().parents[2]
    gaia=json.loads((baseline/'data/raw/gaia-dr3.csv.provenance.json').read_text())
    acquire.fetch('gaia-dr3.csv',gaia['url'])
    coverage=acquire.fetch('irsa-dss.xml','https://irsa.ipac.caltech.edu/cgi-bin/Atlas/nph-atlas',
        dict(mission='DSS',locstr='129.5 72.665 eq',regSize='0.2',covers='on',mode='PI'))
    xml=E.parse(coverage)
    if xml.getroot().get('status')!='ok':raise ValueError('Archive search failed')
    table=acquire.fetch('irsa-images.tbl',xml.find('./images/metadata').text.strip())
    from astropy.table import Table
    names=[str(r['fname']) for r in Table.read(table,format='ascii.ipac')]
    for band,region in [('red','XE036'),('blue','XO036')]:
        if not any(region in f for f in names):raise ValueError(f'Coverage lost for {region}')
        range_cutout.retrieve(band,region)
        h=fits.getheader(out/f'data/raw/{region}.fits')
        if h['REGION'].strip()!=region:raise ValueError('Archive image identity mismatch')
    print('New snapshot saved. Review orbit/version/timing changes before using it; baseline was not modified.')


if __name__=='__main__': main()
