# Research question

Can a bounded, reproducible archive-to-evidence workflow recover documented pre-discovery asteroid trails, reject specified false matches, and retain enough timing, calibration and provenance evidence for an independent scientific reviewer?

The initial setting is one asteroid, two photographic exposures, four published trail-end observations, and four negative controls. It is a small engineering/scientific replication, not a survey completeness or false-discovery-rate measurement. A precovery is an observation made before recognized discovery and identified retrospectively.

## Frontier and evidence

The historical result is the 1955 Palomar precovery of 1999 AN10 in MPEC 1999-N21. That result already exists. The useful contribution here is a runnable local replication and verifier. The local reduction reproduces two image trails and rejects four controls. `discover/frontier.json` is authoritative about this limited claim. Its evidence is author-checked, not independently checked or maintainer reviewed. Sources and claim audits are in research/BRIEF.md.

Progress, in order: (1) frozen-input replication and negative-control rejection; (2) reviewed archive cases with resolved timing and calibrated uncertainties; (3) observation-arc extension from genuinely additional accepted epochs; (4) uncertainty improvement from controlled orbit refits. Throughput and false-positive rates matter only with an independently labeled larger sample. A small residual alone never establishes identity or discovery.

## Verifier and workflow

`PYTHONPATH=src .venv/bin/python -m precovery.reproduce` checks input hashes and published endpoint transcription, predicts historical topocentric positions and uncertainty, verifies coverage, fits stellar calibration and asteroid trails, assesses alternative explanations and writes evidence. Tests independently exercise coordinate wrapping, covariance conversion, invalid timing, maliciously perfect-residual alternatives, corrupted files, injected displaced trails and byte-identical pixels obtained by two archive retrieval routes.

The online acquisition command is optional and writes to a new directory. Offline verification never contacts an archive. Native cutouts preserve the DSS polynomial WCS. The native headers, query results, JPL solution, covariance and Gaia subset are frozen. No new astrometric submission file is produced.

## Bottleneck

The principal bottleneck is scientific review: inconsistent plate timestamps, photographic systematics, trail-end definition and confirmation against the original records. Historical clock errors can dominate a fast-moving object's along-track residual. Gaia propagation is a first-order local correction with quality cuts; perspective acceleration and a full reference-star covariance treatment remain review/improvement tasks. The nominal Gaussian SNR and robust BIC-like statistic are screening summaries, not calibrated discovery probabilities.

## Tracks

The catalogue has archive searches (30), astrometric verification (30), image-pipeline improvements (20) and expert review (20), totaling 100. All begin with a cheap executable check against the frozen case. Search tasks prioritize public single exposures, exact image identity, complete motion/uncertainty coverage and bounded cutouts. Archive responses can be useful even when they expose wrong epochs or missing coverage.

## Settled questions and avoided dead ends

- Position-only archive matches can select the wrong epoch: the real XE035 1953 control demonstrates this.
- Catalogue/coadd coverage is not evidence of a moving source in a single exposure.
- A fit to an orbit already containing the candidate observations is not independent validation.
- An ephemeris covariance is not an image/timing uncertainty; combine the distinct terms without silently treating missing covariance as zero.
- The four published trail ends are two exposures; correlated endpoints do not supply four independent confirmations.
- Reproducing this known case cannot create new arc extension or a new discovery claim.

## Review and contribution policy

No maintainer has agreed; review capacity is undeclared, not fabricated. `review.paused` remains true and the manifest's UNASSIGNED entry is an explicitly vacant role required by the v0 schema. No hubs are authorized. The project can be inspected and developed locally but cannot accept donated scientific contributions yet.

A future reviewer must inspect immutable inputs, calibration and timing provenance, source morphology, motion between exposures, stationary and artifact alternatives, and known-object confusion. Recompute astrometry independently and resolve the plate-log discrepancy before discussing any submission. A current-orbit fit is a consistency check; proposed new observations require appropriate independent linkage/orbit analysis. Record identity, date, evidence hashes, decision, reasons and unresolved obligations. A human approval is a separate decision; software has no path to grant it.

Changes to the verifier or frontier deserve the greatest care. Review capacity, not donated compute, determines the pace.

Credit original observers, archive providers and catalogue teams in every evidence package. Record contributors and reviewers by consent in CONTRIBUTING.md; authorship requires a substantial scientific contribution and agreement. Original code is MIT, original documentation CC BY 4.0; external materials retain their rights.
