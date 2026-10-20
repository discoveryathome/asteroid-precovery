"""Explicit, bounded public-data acquisition. No submission or publishing methods."""
import json
from datetime import datetime, timezone
from pathlib import Path
import requests
from .verify import sha256

ROOT = Path(__file__).resolve().parents[2]


def fetch(name, url, params=None, limit=20_000_000):
    path = ROOT / 'data/raw' / name
    if path.exists():
        raise FileExistsError(f'{path}: snapshots are immutable; use a new name')
    with requests.get(url, params=params, timeout=(15, 120), stream=True) as r:
        r.raise_for_status()
        chunks = []; size = 0
        for chunk in r.iter_content(65536):
            size += len(chunk)
            if size > limit: raise ValueError('Download exceeded byte budget')
            chunks.append(chunk)
        body = b''.join(chunks)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
        meta = dict(path=str(path.relative_to(ROOT)), url=r.url,
                    retrieved_at=datetime.now(timezone.utc).isoformat(),
                    sha256=sha256(path), bytes=size, content_type=r.headers.get('Content-Type'),
                    etag=r.headers.get('ETag'), last_modified=r.headers.get('Last-Modified'))
        path.with_suffix(path.suffix+'.provenance.json').write_text(json.dumps(meta, indent=2)+'\n')
    return path


def horizons(name, start, stop, step='1 m'):
    return fetch(name, 'https://ssd.jpl.nasa.gov/api/horizons.api', dict(
        format='text', COMMAND="'137108;'", CENTER="'675'", MAKE_EPHEM='YES',
        EPHEM_TYPE='OBSERVER', START_TIME=f"'{start}'", STOP_TIME=f"'{stop}'",
        STEP_SIZE=f"'{step}'", QUANTITIES="'1,3,9,36,37'", CSV_FORMAT='YES',
        ANG_FORMAT='DEG', EXTRA_PREC='YES', TIME_DIGITS='SECONDS'))
