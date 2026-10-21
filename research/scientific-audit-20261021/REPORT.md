# Bounded scientific review and fresh setup — 21 October 2026

A reviewer-ready timing dossier is prepared for **one exposure, XE036 (red)**.
Scientific signoff remains pending: the project owner confirmed that a qualified
reviewer is still pending. No reviewer is named, no review is fabricated, and no
astrometry is submitted. This package makes the next scientific decision concrete.

## The decision to request

Assess whether the published start/end timing can be used for retrospective
replication of this exposure, and identify the original-log evidence still needed.
This decision alone does not accept the astrometry, asteroid identity, new
observational findings or an MPC submission.

| Evidence | Recorded value | Interpretation |
|---|---|---|
| Native IRSA DATE-OBS | 1955-01-25T10:05:00 | Archive metadata; conflicts with publication |
| Published first endpoint | 1955 Jan 26.40208 | Approximately Jan 26 09:39 |
| Published second endpoint | 1955 Jan 26.43333 | Approximately Jan 26 10:24 |
| Existing replication timescale | UT1 | Must not silently be treated as modern UTC |
| Existing uncertainty floor | 30 seconds | Policy assumption, not calibrated plate-clock error |
| Earlier STScI transcription | Jan 25 09:65 | Invalid clock string; not silently normalized |

`review-request.json` binds ten input/source/review files by SHA-256, includes the
actual two published lines from the cached MPEC, and records the unresolved timing
question. It references the native and independently retrieved image and metadata,
image-version lock, published circular, and existing review obligations. The
published source is [MPEC 1999-N21](https://www.minorplanetcenter.net/iau/mpec/J99/J99N21.html).
No original plate log has been acquired in this pass. The archive-versus-publication
conflict is documented, not resolved by selecting the time giving smaller residuals.

## Proposed review process

1. The maintainer obtains consent from a qualified astronomer/astrometrist and
   agrees the single timing obligation and review capacity. The role is currently
   unassigned. No contact has been sent on the project's behalf.
2. The reviewer starts with `review-request.json`, runs the hash check below,
   inspects the original files and independently checks the time transcription.
   They request original observing/plate-log evidence where necessary; new evidence
   gets a source URL, retrieval details, raw file and hash in `additional_evidence`.
3. In a **copy** of the request, the reviewer records their identity, relevant
   qualifications, consent, independence disclosure, review date, independent
   commands/checks, rationale, unresolved issues and one scoped verdict:
   `supported-for-retrospective-use`, `needs-more-evidence`, or `rejected`.
   “Needs more evidence” is a useful completed review decision, not a failure.
4. The maintainer checks provenance and the reviewer's identity/qualifications,
   then decides how that bounded verdict affects the outstanding obligations in
   `review/record.json`. Software only validates record completeness and hashes;
   it cannot authenticate a reviewer or grant scientific acceptance.
5. New observation acceptance or any external reporting remains a separate,
   expressly authorized scientific decision. This dossier cannot enable either.

```sh
python3 research/scientific-audit-20261021/review_package.py
# After a real review, validate its copied record:
python3 research/scientific-audit-20261021/review_package.py --record /path/to/reviewer-decision.json
python3 -m unittest discover -s research/scientific-audit-20261021 -p 'test_*.py' -v
```

The five regression tests reject altered evidence, silent time changes, anonymous
non-pending verdicts and automatic scientific acceptance. The original review
record remains pending and unchanged. A structurally valid decision is not proof
of the reviewer's credentials or correctness.

## Fresh setup results

A new isolated Python 3.14.6 environment was installed from `requirements.lock`
in this existing worktree. All nine project tests passed, both documented positive
exposures reproduced, all four negative controls were rejected, and the MS4
archive package verifier passed. **All 12 generated outputs were byte-identical
to the frozen derived files.** `cold-receipt.json` records hashes and comparison
results; installation, tests, replay and MS4 logs are retained alongside it.

This was a coordinator-run setup rehearsal, with ordinary package caching allowed.
It is not a new independent donor session, a usability study, or a scientific
review. It establishes that a fresh environment can execute the bounded technical
work. The broader scientific pipeline, calibration and frozen time overrides were
not changed.

Reproduce from the repository root:

```sh
uv venv --python 3.14 .replay/new-review/venv
uv pip install --python .replay/new-review/venv/bin/python -r requirements.lock
PYTHONPATH=src .replay/new-review/venv/bin/python -m unittest discover -s tests -v
PYTHONPATH=src .replay/new-review/venv/bin/python -m precovery.reproduce .replay/new-review/output
.replay/new-review/venv/bin/python research/rehearsal-ms4/verify_package.py
```

Keep participation limited to bounded evidence work while reviewer capacity is
unassigned. The next real milestone is a qualified person's documented decision
on this dossier, followed by an independent contributor trial with orientation
time, completeness, review effort and acceptance latency recorded. No population
false-alarm rate, new observation, improved orbit, or scientific acceptance is
claimed from these controls.
