"""Fail-closed screening and provenance checks, separate from scientific acceptance."""
import hashlib
import math
from datetime import datetime
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_manifest(root, manifest):
    for row in manifest['files']:
        path = (Path(root) / row['path']).resolve()
        if not path.is_relative_to(Path(root).resolve()):
            raise ValueError('Manifest path escapes project')
        if sha256(path) != row['sha256']:
            raise ValueError(f"Checksum mismatch: {row['path']}")


def strict_time(value):
    # datetime deliberately rejects e.g. 09:65; FITS libraries may normalize it.
    return datetime.fromisoformat(value)


def screen(m):
    reasons = []
    thresholds = [('snr', 5.0, 'min'), ('residual_arcsec', 12.0, 'max'),
                  ('length_ratio', 0.6, 'min'), ('length_ratio', 1.5, 'max'),
                  ('angle_error_deg', 12.0, 'max'), ('delta_bic_vs_point', 10.0, 'min')]
    for key, limit, direction in thresholds:
        v = m.get(key)
        if v is None or not math.isfinite(v) or (v < limit if direction == 'min' else v > limit):
            reasons.append(key)
    for key in ['epoch_matches', 'time_provenance', 'inside_coverage', 'calibration_valid', 'fit_converged']:
        if m.get(key) is not True:
            reasons.append(key)
    for key in ['stationary_match', 'artifact_flag', 'edge_flag']:
        if m.get(key) is not False:
            reasons.append(key)
    return {'screen_pass': not reasons, 'reasons': sorted(set(reasons)),
            'scientifically_accepted': False, 'submission_allowed': False,
            'status': 'awaiting-independent-review' if not reasons else 'rejected-by-screen'}
