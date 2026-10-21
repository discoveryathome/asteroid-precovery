# Donor rehearsal: visible progress is achievable, with smaller tasks

The seeded archive-search task produced five archival cutouts for a second documented object, 2002 MS4, plus two explicitly rejected wrong-epoch responses. The first two cutouts arrived about **6 minutes 37 seconds** after the session started; all five arrived by **10 minutes 16 seconds**. This supports an accessible initial contribution queue. It does not establish the rate of finding previously unreported observations.

The session had a 30-minute maximum and ended after **17 minutes 52 seconds**, after packaging a concrete result. See `session.json` for the final timestamps and duration. No baseline changes were merged, and nothing was published, submitted or announced.

## What the contributor produced

Open the [evidence package](../rehearsal-ms4/README.md) or the [five-exposure comparison](../rehearsal-ms4/five-exposure-preview.png).

- Seven documented historical epochs inventoried against public archive responses.
- Five small FITS cutouts: three exact published region/plate-ID matches, and two visually matching fields with provisional legacy-name crosswalks.
- Two requests that silently selected the wrong exposure, retained with IDs and dates so a future verifier can reject them.
- Source URLs, request parameters, response metadata, cutout and range checksums, original headers, frozen JPL predictions/covariance and solution versions.
- An offline verifier that passed 106 file hashes, identity/epoch checks and native WCS preservation checks.
- Twenty-two additional literature links to triage, beyond the investigated case. These are leads, not validated available controls.

This is **new project material reproducing access to old, published science**. No new image astrometry was measured for MS4, no scientific acceptance occurred, no observation arc was extended, and no orbital-uncertainty improvement was demonstrated. The original AN10 reproduction remains the pilot's completed scientific-control replay.

## Observed timeline

Times are measured from 2026-10-20 06:21:47 UTC, including reading, commentary and orchestration. Research proceeded while dependencies installed and the baseline checks ran.

| Elapsed | Milestone |
|---|---|
| 4m 29s | Fresh dependency installation completed; installation itself took about 89 seconds |
| 6m 26s | Baseline reproduction passed |
| 6m 30s | All nine baseline tests and the selected archive starter completed successfully |
| 6m 37s | Both early Palomar cutouts and their provenance were saved |
| 7m 48s | First side-by-side image preview saved |
| 10m 16s | All five selected cutouts acquired |
| About 16m | Complete package integrity/identity/WCS checks passed |

Evidence: `events.jsonl`, raw acquisition timestamps, `baseline.log`, `tests.log`, `starter.log`, and `package-verification.log`. The session's bounded work and final packaging took less than 30 minutes; it was not padded with idle time.

## Friction that matters for launch

1. **The task is too broad for a short-session success criterion.** The selected card asks for a second case and all overlapping exposures. Five useful images still leave it formally partial: two exact requested plates and a full timing/identity audit remain unresolved. A donor dashboard should count bounded acquisition, measured rejection and review-ready evidence as distinct useful outputs.
2. **Archive defaults can be wrong.** Two STScI requests returned a real FITS image of the wrong epoch. The absence of an HTTP error is insufficient. The explicit plate checks made these responses productive negative examples.
3. **An archive request failed.** The IRSA cutout request returned HTTP 502. Its cause was not diagnosed; bounded byte-range retrieval worked. A cold contributor should find that fallback in the task instructions instead of having to discover it in source code.
4. **Cold plotting has visible setup cost.** Matplotlib built a fresh font cache, making the first baseline run close to two minutes. The process completed; future sessions should show progress so this does not look like a hang.
5. **Some scientific work remains difficult.** Legacy plate aliases, date discrepancies, uncertainty calibration and independent review were not resolved. They should be separate advanced tasks. An unqualified agent cannot satisfy the expert-review card by running more code.
6. **The new object needs an appropriate model.** Its old exposures contain only about two arcseconds of predicted motion. The long-trail classifier should not be relaxed to force a pass. Point-source/blend/injection studies provide meaningful next work.

## Proposed launch shape

Give an ordinary first contributor one concrete image or one documented source, an executable first experiment, a budget, and a result format. Accept a useful bounded result such as an exact exposure crosswalk, a frozen cutout with provenance, an empirical recovery plot, or a precisely demonstrated wrong match. Keep scientific promotion and expert approval separate.

Draft follow-ups are in [proposed-starter-tasks.json](../rehearsal-ms4/proposed-starter-tasks.json): quantify an image crosswalk, measure seeded point-source injection recovery, and acquire one exposure from one additional documented case. Their durations are suggestions, not tested guarantees. They have not replaced the project's existing catalogue.

A donor-facing progress display could report: **5 cutouts acquired; 3 exact plate IDs matched; 2 provisional aliases awaiting review; 2 wrong-epoch responses rejected; 0 new observations claimed.** This is visible, honest progress.

## Limits of this rehearsal

This was the original author's **warm-context** rehearsal in a fresh local checkout and virtual environment. It used one OS-random seed, the real task picker and the pinned baseline commit; the task was not cherry-picked. The seed, scope key and catalogue/frontier hashes are in `session.json`. Existing project context and familiarity with the range downloader probably reduced friction. It was not an independent blind donor test, a representative success-rate estimate, or an end-to-end test of the future public hub and launcher.

My recommendation is to keep this as an inaugural-project candidate, revise the starter catalogue into bounded tasks, and run a cold contributor test before promising a uniformly easy experience. Review capacity remains a separate launch requirement: the existing project is still paused with no consenting scientific maintainer assigned.

The isolated contribution is commit `48a1f3638bfdefc57453c27d0b136bd4355f5fd6` on `rehearsal/archive-ms4`, with a verified incremental [Git bundle](contribution.bundle) requiring the baseline commit recorded in `session.json`.

This reference copy relocates evidence-package links. The original report and full
session receipts remain in asteroid-precovery/donor-rehearsal-20261020 in the
incubatedPilots workspace; contribution source commit: 48a1f3638bfdefc57453c27d0b136bd4355f5fd6.
