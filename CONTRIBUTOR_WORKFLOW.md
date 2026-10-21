# Contributor workflow — 21 October 2026

Start with `python3 tools/pick_task.py`. An ordinary draw excludes the qualified
expert-review track; a qualified reviewer may opt in with `--expert`. The maintainer
of record is @Mister-Meeseeks. Scientific review is still a distinct requirement.

Install the pinned environment as described in README, then reproduce the original
AN10 control into a separate output directory:

```sh
PYTHONPATH=src .venv/bin/python -m precovery.reproduce .replay/session
.venv/bin/python research/rehearsal-ms4/verify_package.py
```

The second command checks the imported MS4 rehearsal package: 106 hashes, three
exact plate-ID matches, two provisional aliases and two rejected wrong-epoch
responses. All seven epoch records remain evidence about published observations.
No new MS4 image astrometry, arc extension or discovery is claimed. The package's
original manifest is preserved byte for byte.

A first 30-minute session should tackle one exposure, one alias, one documentary
timestamp discrepancy, or one rejected-match audit. Reserve approximately ten
minutes for verification and handoff. Acquisition failure is completable when the
request, returned identity, error and bounded fallback are recorded precisely.

An exposure receipt must contain the primary observation citation; requested and
returned image/plate IDs and epochs; raw headers; URLs/request parameters and access
date; cutout bounds and native coordinate metadata; file hashes; and unresolved
issues. Never accept a FITS response solely because the HTTP request succeeded.
Wrong-epoch images are useful negative controls, not detections. The rehearsal's
`native_cutout.py` and provenance records document the bounded range-download route
used when the IRSA cutout endpoint failed; use it only for a source whose format
and access terms match, and record byte ranges and original headers.

Technical package completion means reproducible acquisition/identity evidence.
Scientific acceptance requires a qualified reviewer to assess timing, calibration,
uncertainty, morphology, alternative explanations and independent linkage where
applicable. Observation reporting is a later, separately authorized decision.
The existing scripts have no submission path. Matching fields do not prove matching
exposures; do not relax the long-trail classifier for MS4's nearly pointlike motion.

## Session receipt and review

Record the base commit, selected task and seed, input and output hashes, commands,
budget and elapsed time, successful and failed attempts, verification output,
and what remains unresolved. Check the live Handoffs Discussion before claiming
work. Post or open a focused PR only within the contributor's authorization.
A task draw is advisory; it does not reserve work.

A useful outcome can be a new witness, an independently checked result, a
source-backed correction, or a reproducible blocked handoff. Repeating a known
positive control is onboarding evidence. Review records must distinguish code
replay, independent implementation, scientific review and maintainer acceptance.

For the next monitored cohort, record cold-start orientation time, first useful
artifact, reviewer minutes, acceptance latency, duplicate work and unsuccessful
sessions. Later sessions inherit every known valid result. A lack of improvements
is evidence to revise the task pool, not a reason to reset it to weaker baselines.

## Bounded timing review and fresh setup

[The 21 October review package](research/scientific-audit-20261021/REPORT.md)
prepares one hash-bound XE036 timing obligation and a review decision template.
A fresh environment passed the controls and reproduced all 12 derived outputs
byte-for-byte. Qualified reviewer assignment and independent donor onboarding
remain pending. The technical setup pass cannot grant scientific acceptance.
