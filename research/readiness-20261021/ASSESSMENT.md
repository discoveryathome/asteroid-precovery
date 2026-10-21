# Readiness assessment — asteroid-precovery — 21 October 2026

Recommendation: **Bounded evidence-collection beta with qualified scientific review**.

The original broad archive assignment remained incomplete despite useful evidence. Ordinary task draws could assign expert scientific review, and the tested second-case evidence lived outside the project.

## Completed changes

- Imported the existing MS4 rehearsal package with its original 106-file manifest, attribution and scientific nonacceptance record.
- Replaced the broad archive task with one-exposure acquisition, one-alias resolution and one rejected-match audit; narrowed the timing task to one plate discrepancy.
- Added a seeded ordinary picker that excludes expert review unless explicitly opted into; specified technical package completion separately from scientific acceptance.
- Documented wrong-epoch controls, bounded archive fallback and a compact evidence receipt.

## Validation

Original AN10 offline replay and nine tests passed; MS4 hashes, seven epoch records, three exact IDs, two provisional aliases, two wrong-epoch rejections and native WCS preservation passed; 120 ordinary draws never selected expert review; schema validation passed. Reused the existing pinned Python environment; no cold installation or fresh donor session was measured.

All 6 active task cards' first-step commands passed; repeated identical commands were cached within the check.
Executed outputs are in [checks/](checks/). This is coordinator review and replay,
not an independent cold contributor study or external scientific signoff. No new
research-success rate is inferred. Current workflow: [CONTRIBUTOR_WORKFLOW.md](../../CONTRIBUTOR_WORKFLOW.md).

## Next scientific handoff

A qualified reviewer should accept or reject one bounded scientific obligation and define how technical packages enter the research workflow. Test cold onboarding on the revised cards.

Acceptance test: A documented scientific review path, a completed cold task with independently inspectable evidence, and explicit unresolved timing/identity/calibration issues. No automated observation reporting.

```sh
PYTHONPATH=src .venv/bin/python -m precovery.reproduce .replay/reviewer
.venv/bin/python research/rehearsal-ms4/verify_package.py
python3 tools/pick_task.py
```

## Integration state

Base: `438be093db77dd3f8ca955775b0f59ee6b742928`. Remote main at fetch: `1f92397b2c227ea6fab475af9d20c3f1b3519415`.
This work builds on the newer local launch/coordination branch. The existing
[launch PR #1](https://github.com/discoveryathome/asteroid-precovery/pull/1) was still open when
checked, and its Handoffs thread had no comments. Local launch decisions are not
proof that remote main or the public launcher has adopted them. No remote write,
merge, announcement, scientific acceptance or observation submission occurred.
