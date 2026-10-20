# Independent review guide

Status: local replication complete; scientific acceptance pending. No reviewer or maintainer has consented. There is no astrometric submission.

1. Run `PYTHONPATH=src .venv/bin/python -m precovery.reproduce .replay/review`. Input hashes, original-circular transcription and image identities must pass before reduction.
2. Compare the two raw native FITS cutouts and original headers with their source URLs and byte-range provenance. The tests compare their overlapping pixels to independently retrieved IRSA service cutouts exactly. Inspect `research/sources/irsa-images.tbl` for coverage and plate alternatives.
3. Read `data/controls.json` before interpreting times. The 1955 dates/times are a **published-record override** of conflicting archive headers. Do not silently normalize an invalid `09:65` string or treat UT1 as UTC. Confirm the original observing log independently.
4. Inspect `data/derived/evidence.png` and the individual positive-control FITS cutouts. Red is faint and elongated; blue is stronger. The plotted model is fitted to pixels, initialized from the predicted trail. Check possible blends, plate defects, nearby saturated stars, scan artifacts and whether fitted trail ends are reliable. The fit uses catalogue-star masks and a robust loss, so model-residual panels alone are not a goodness-of-fit proof.
5. Recompute the astrometry using another reduction method. Review Gaia IDs, proper-motion cuts, affine coefficients, rejected stars and leave-one-out residuals in `calibration.json`. Systematic floors and time errors should be reassessed, not optimized away.
6. Compare measured and predicted midpoints, two-exposure motion and all four published endpoints. Current JPL#100 already uses the historical observations; this is an explicitly non-blind check. No before/after orbit fit has been made.
7. Confirm all four negative controls fail. The stationary-source and shifted-search controls use actual archive pixels. The artifact is an explicitly synthetic single-pixel injection with a fixed seed. The wrong-epoch control is rejected from real image metadata before any candidate measurement. None estimates a population false-discovery rate.
8. Record a decision in `review/record.json` with reviewer consent, identity, evidence hashes, independent commands and unresolved issues. A software pass cannot change scientific acceptance or allow a submission.

## Interpretation of screening statistics

The source model is a uniform line convolved with a circular Gaussian plus background. Center, length, angle, width, amplitude and background are fitted within bounds around the prediction. A point-source alternative is fitted separately. `delta_bic_vs_point` uses robust deviance and a parameter penalty; it is BIC-like, not a calibrated Bayes factor. SNR has a factor-two noise inflation for correlated plate pixels; empirical noise/false-alarm calibration remains open.

Residuals use wrapped RA times cos(dec). Measurement covariance combines the fit centroid covariance, a stellar-calibration floor and 30 seconds of timing uncertainty projected along motion. The frozen Horizons error ellipse is then added for a consistency diagnostic. Acceptance also needs morphology, independent corroboration and reviewed timing; no residual threshold is a discovery criterion.

The cached native WCS has GSC1 polynomial coefficients. It is retained, then locally corrected against Gaia DR3. The original published comparison used USNO-A2.0. These are distinct calibrations and must not be conflated.
