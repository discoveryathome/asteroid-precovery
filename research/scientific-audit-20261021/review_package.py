#!/usr/bin/env python3
"""Build/verify a hash-bound, bounded scientific-review request. No acceptance writes."""
import argparse,copy,hashlib,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
FILES=['data/controls.json','data/versions.json','data/raw/XE036.fits','data/raw/XE036.original-header.txt','data/raw/XE036.provenance.json','data/raw/irsa-XE036.fits','data/raw/irsa-XE036.fits.provenance.json','research/sources/mpec-1999-n21.html','research/sources/irsa-images.tbl','review/record.json']

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def build():
    controls=json.loads((ROOT/'data/controls.json').read_text());red=controls['positive_controls'][0]
    assert red['id']=='an10-red-1955'
    h=(ROOT/'data/raw/XE036.original-header.txt').read_text()
    archive=re.search(r"DATE-OBS\s*=\s*'([^']+)'",h).group(1)
    source=(ROOT/controls['time_override']['snapshot']).read_text()
    lines=[x.strip() for x in source.splitlines() if 'J99A10N 5 1955' in x]
    assert len(lines)==4 and '26.40208' in lines[0] and '26.43333' in lines[1]
    return {'schema':'bounded-scientific-review/v1','id':'an10-XE036-timing-20261021','status':'pending-qualified-reviewer','scope':'Assess whether the published start/end timing for the single XE036 exposure can be adopted for retrospective replication, and what original-log evidence is still needed. No astrometric, identity, new-observation or submission acceptance.', 'evidence':[{ 'path':p,'sha256':digest(ROOT/p)} for p in FILES], 'observations':{'archive_DATE_OBS':archive,'published_lines':lines[:2],'replication_start':red['start'],'replication_end':red['end'],'working_timescale':red['timescale'],'time_sigma_seconds':red['time_sigma_seconds'],'sigma_status':'policy floor; not measured clock uncertainty','known_conflict':controls['time_override']['conflict']},'reviewer':None,'consent':False,'reviewed_at':None,'verdict':'pending','rationale':None,'independent_checks':[],'additional_evidence':[],'unresolved':['Original plate/observing-log timing and date convention not yet established.','Published fractional days do not independently authenticate an archive-to-plate alias.','30 second uncertainty floor requires scientific justification; not clock calibration.'],'scientifically_accepted':False,'submission_allowed':False,'publication_allowed':False}

def check(record):
    if record.get('schema')!='bounded-scientific-review/v1':raise ValueError('wrong schema')
    expected={p:digest(ROOT/p) for p in FILES}
    actual={x['path']:x['sha256'] for x in record['evidence']}
    if len(record['evidence'])!=len(expected) or actual!=expected:raise ValueError('evidence mismatch')
    original=build()
    for key in ['id','scope','observations']:
        if record[key]!=original[key]:raise ValueError('altered request or observations')
    if any(record[k] is not False for k in ['scientifically_accepted','submission_allowed','publication_allowed']):raise ValueError('bounded record cannot authorize scientific acceptance or external action')
    verdict=record['verdict']
    if verdict not in ['pending','supported-for-retrospective-use','needs-more-evidence','rejected']:raise ValueError('unknown verdict')
    if verdict!='pending':
        reviewer=record.get('reviewer')
        if not isinstance(reviewer,dict) or not all(isinstance(reviewer.get(k),str) and reviewer[k].strip() for k in ['name','qualification','independence_disclosure']):raise ValueError('identified qualified reviewer required')
        if record.get('consent') is not True or not record.get('reviewed_at'):raise ValueError('consent and review date required')
        if not record.get('rationale') or not record.get('independent_checks'):raise ValueError('reasoned, reproducible decision required')
    return {'status':'PASS','verdict':verdict,'evidence_files':len(expected),'scientific_identity_verification':'not performed by software','scientific_acceptance':False}

def main():
    p=argparse.ArgumentParser();p.add_argument('--build',action='store_true');p.add_argument('--record',type=Path,default=HERE/'review-request.json');args=p.parse_args()
    if args.build:
        if args.record.exists():raise ValueError('refusing to replace existing review record')
        args.record.write_text(json.dumps(build(),indent=2)+'\n')
    print(json.dumps(check(json.loads(args.record.read_text())),indent=2))
if __name__=='__main__':main()
