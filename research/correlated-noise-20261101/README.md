# Correlated noise, sky clutter and the trail screen (2026-11-01)

The per-exposure trail screen was run unchanged on real plate pixels away from the AN10 trail
(masked on XE036 and XO036; absent from the 1953 XE035 plate of the same sky), and on synthetic
trails injected into those pixels. With the candidate placed at the prediction it passes 2.1%
(red XE036), 6.5% (blue XO036) and 16.4% (XE035) of blank trials. Both documented positives still
rank above every one of those passes on their own plate. Seven screen changes were priced on this
sample. Two help, each at a recall cost: a per-plate noise inflation (XO036 false passes 80 → 16,
recall lost below true SNR 8) and an along-track uniformity statistic (most false passes removed on
all three plates, 8–19% of injected passes lost). The rest fail. On XE036 the length gate passes
straight features up to 10× the predicted length at 9 of 12 probe sites.

Nothing under `src/`, `data/`, `review/` or `discover/` changes. Every number here is
author-checked only. No result is scientifically accepted, and nothing here is a population or
survey false-alarm rate.

## Task card

- **Card:** `correlated-noise-and-trails` (track `image-pipeline-improvements`, priority_revision 2),
  drawn with `python3 tools/pick_task.py --seed 83fdb98fa8b2f62bfc23bebe32115b3f`. Scope not changed.
- **Base:** `9795ae8e87ba525d57b7cb6a135abb8a9eb34c56`.
- **Focus:** explore. The new approach (a whitened trail statistic) was priced against its
  favourable best case before any implementation.
- **Plates:** XE036 (red, 103aE, 45 min, 1.70″/px), XO036 (blue, 103aO, 12 min, 1.01″/px) and
  `irsa-XE035.fits` (1953-02-18, 103aE, 40 min, same sky, no AN10 trail). All frozen inputs.
- **Pipeline code used unchanged:** `calibrate`, `fit_trail`, `screen`, `predicted_pixel`,
  `measured_sky`. One experiment (`margin_fix.py`) re-executes `fit_trail`'s own source with its
  20 px window margin made a parameter; at 20 px it reproduces both positives exactly.

## 1. Noise inflation: measured, not 2

`noise_model.py` masks every Gaia source (no magnitude cut), 3σ outliers, plate edges and the AN10
trail, removes a plane per 64 px tile, and measures the pixel autocorrelation R to 25 px.
Correlation is gone by lag 3 (XE036) or 4 (XO036); XE035 keeps a weak tail (0.03 at lag 4).

For each plate's predicted trail template, F = sqrt(TᵀRT / TᵀT) is how much a white-noise SNR
overstates the true SNR of the same filter. The pipeline hard-codes F = 2.

| Plate | R(1,0) / R(0,1) | F (measured) | F (grain-only synthetic check) |
|---|---|---|---|
| XE036 | 0.33 / 0.41 | 1.61 | 1.62 |
| XO036 | 0.74 / 0.73 | 2.85 | 3.01 |
| XE035 | 0.36 / 0.45 | 2.02 | 2.09 |

F uses R truncated at 15 px with a triangle taper (`L15_tapered`). Other truncations give
1.58–1.67 (XE036) and 2.66–2.90 (XO036), excluding kernels that aren't valid covariances. The last
column is the scatter of a fixed linear trail filter over the same trial grid on synthetic Gaussian noise built from each plate's measured R
(`grain_null.py`). It checks the computation by simulation; it reuses R, so it isn't independent.

Against grain noise, then, the factor 2 is conservative on the red plate and anti-conservative on
the blue one. Against the real trail-scale scatter (section 3) it is too small on all three.
Correcting both positives with their own F gives SNR 19.5 (red) and 31.7 (blue); both still pass.

## 2. The new statistic, priced first: stopped

A whitened (generalized least squares) trail filter is optimal under correlated noise. Its best
case over the white filter, with a perfectly known, stationary covariance and no masked pixels, is
sqrt(TᵀR⁻¹T · TᵀRT) / TᵀT:

| Plate | Best-case gain (valid kernels) |
|---|---|
| XE036 | 1.03–1.05 |
| XO036 | 1.03–1.04 |
| XE035 | 1.10–1.22 |

On the two plates that contain AN10 the best case is under 5%, below the 10% line set before
pricing, so it wasn't implemented. The correlation length (2–3 px) is short against the trail, so
the noise spectrum is nearly flat over the template's band. XE035, with a longer tail, clears the
line, so a plate like it could justify a whitened filter. Untapered kernels that are not valid
covariances (negative spectrum) were excluded; the XO036 ones also failed to converge.

Section 3 suggests a further limit, not tested here: on real pixels the noise that matters at trail
scale is mostly static sky, which a stationary grain covariance doesn't describe. A whitener built
from a clutter-inclusive covariance wasn't priced.

## 3. Labeled sample on real pixels

`labeled_sample.py` places the predicted trail vector, rotated by 0°, 45°, 90° and 135°, on a grid
of sites (24 px step on the red plates, 30 px on XO036), masks calibration stars as the pipeline
does for positives, masks the AN10 trail on XE036 and XO036 (and skips sites whose track comes
within 15 px of it), and runs `fit_trail` + `screen`. The residual gate is set as satisfied
(residual 0), so each number is the chance that a blank site **at the prediction** passes every
other gate. Neighbouring trials overlap, so they're correlated.

| Plate | Trials (site × rotation) | Screen passes | Grain-only synthetic passes | Trail-filter scatter, real / grain-only |
|---|---|---|---|---|
| XE036 | 283 | 6 (2.1%) | 0 | 3.62 / 1.62 |
| XO036 | 1230 | 80 (6.5%) | 41 (3.3%) | 3.89 / 3.01 |
| XE035 | 372 | 61 (16.4%) | 2 (0.5%) | 4.50 / 2.09 |

Scatter is the robust standard deviation of a fixed linear trail filter, in the pipeline's own
noise units. None of the 147 passes has any Gaia source, at any magnitude, within 4″ of its fitted
centre. `null-passes-montage.png` shows up to 8 of the strongest per plate. Most of those are
trial tracks with a bright, often saturated star near one end, or a chain of faint stars along the
track. The trail model absorbs that light while the point-source alternative is held within ±12 px of the centre,
and the stationary gate only checks the trail's centre. On XO036 a grain-only version of the plate
still gives 41 passes against 80 real ones: there the nominal SNR ≥ 5 cut (about 10 in the linear
filter's units) sits at about 3.3σ of the grain-only scatter (10 / 3.01).

**The positives against their own plate's null.** Red AN10: SNR 15.7, ΔBIC 394; the largest
XE036 null pass is SNR 9.6, ΔBIC 125. Blue AN10: SNR 45.2, ΔBIC 1767; the largest XO036 null pass
is SNR 21.0, ΔBIC 427. On XE035, which has half XE036's grain noise, 8 null passes exceed SNR 15.7
and 5 exceed ΔBIC 394.

**Recall.** One synthetic trail per trial (pipeline profile, plate PSF, predicted length and angle,
up to 3 px off the site, true SNR cycling through the six levels), added to the real pixels, with
amplitude set by a true SNR computed from the measured F and grain σ. A pass requires every gate including a residual under 12″.

| True SNR | 3 | 4 | 5 | 6 | 8 | 12 |
|---|---|---|---|---|---|---|
| XE036 (47–48 each) | 0.02 | 0.04 | 0.17 | 0.47 | 0.81 | 0.87 |
| XO036 (205 each) | 0.43 | 0.66 | 0.82 | 0.90 | 0.95 | 0.95 |
| XE035 (62 each) | 0.11 | 0.21 | 0.27 | 0.50 | 0.63 | 0.82 |

The red screen is strict and the blue one loose, which matches their F. At true SNR 8–12 the
losses differ by plate: on XE035 mostly position capture (residual > 12″), on XE036 a few each of
residual and ΔBIC, on XO036 a spread of angle, ΔBIC, length, residual and stationary-match failures.

## 4. Over-length straight features pass the length gate

`morphology_probe.py` adds straight features through 12 sites per plate, at the predicted angle
and at +6°, with the trail's peak amplitude (true SNR 10 for the 1× case) and 1.0–10× the
predicted length. A 10× feature crosses the whole cutout, like a meteor, scratch or emulsion
defect.

| Plate, angle +0° | 1× | 1.3× | 1.5× | 2× | 4× | 10× |
|---|---|---|---|---|---|---|
| XE036 passes /12 | 11 | 9 | 9 | 9 | 9 | 9 |
| XE036 median fitted length ratio | 1.02 | 1.14 | 1.14 | 1.14 | 1.14 | 1.14 |
| XE035 passes /12 | 11 | 7 | 6 | 6 | 6 | 6 |
| XO036 passes /12 | 12 | 10 | 5 | 3 | 3 | 3 |

At +6° XE036 still passes 3 of 12 at every length from 1.3×. Widening the fit window (section 5)
doesn't change this. It points to the optimizer rather than the window. `length_diag.py` takes
three of the XE036 probe sites and a 10× feature: with the length held fixed and the other six
parameters refitted, the cost keeps falling to the edge of the scan (1.55×, the fit's upper bound)
on two sites and is lowest at 1.45× on the third, but the full fit stops at 1.00–1.29× for margins
of 20, 54 and 100 px. The fit starts at the predicted length, and at this per-pixel signal the cost
surface along length is shallow, so the fitted length mostly reports where the fit began. Three
sites only.

## 5. Fixes priced on this sample

Thresholds come from the null and injection samples only. The two positives are reported, never
used to tune.

| Fix | Null passes (XE036 / XO036 / XE035) | Recall cost | Both positives |
|---|---|---|---|
| None (pipeline) | 6 / 80 / 61 | n/a | pass |
| SNR × 2/F, measured per plate | 6 / **16** / 61 | XO036: 0.95 → 0.94 at true SNR 8, but 0.82 → 0.57 at 5 and 0.90 → 0.79 at 6; XE036 recall rises (0.47 → 0.74 at 6) | pass (19.5, 31.7) |
| SNR × 2/s, s = real-pixel filter scatter | 1 / 8 / 18 | large: XE036 recall 0.55 at true SNR 12 | pass (8.7, 23.3) |
| Along-track uniformity U ≤ p95 of injected passes at true SNR ≥ 8, per plate (2.71 / 2.26 / 3.0) | 0 / 34 / 14 | 8% (XE036), 13% (XO036), 19% (XE035) of injected passes at all SNRs lost | pass (U = 1.27, 1.07) |
| Stationary gate on the whole fitted segment | 6 / 73 / 32 (catalogue, 4″) | 15% / 7% / 22% of injected passes | **blue fails**: a catalogue star lies 0.1″ from its fitted segment |
| Same, every Gaia star, magnitude radius | 0 / 52 / 3 | 63% / 29% / 68% | **both fail** |
| Wider window (20 px + 0.25 × length) | 6 / 67 / 47 | at true SNR 6: XE036 0.47 → 0.38, XO036 0.90 → 0.88, XE035 0.50 → 0.42; over-length features still pass | pass (15.8, 45.9) |
| Different-epoch subtraction (XE036 − XE035) | 175 covered trials: 1 → 0 (all-pixel flux scale), 1 → 9 (star-pixel scale) | all-pixel: 0.46 → 0.36 at SNR 6; star-pixel: 0.91 → 0.74 at SNR 8 | red passes (15.5, 14.1) |

U is the largest of six along-track segment amplitudes over the median of the other five. A real
trail is uniform; a star at one end isn't. Null passes have median U of 3–600 by plate, injected
passes 1.5–1.8. The catalogue-proximity vetoes fail because real trails cross stars: the blue AN10
trail passes 0.1″ from a catalogue star that the pipeline masks. Subtraction with one global
linear flux scale fails: fitted on all pixels, the scale is diluted toward the sky (0.33) and stars
stay unsubtracted; fitted on star pixels (1.18), the difference noise rises 16% and new false
passes appear, likely from star residuals on density images. The two PSFs already match
(σ 1.02 and 1.00 px), so PSF matching isn't the cause.

## Replay

From the repository root, in the pinned environment (Python 3.14.9 here, single-threaded BLAS,
offline). Each script reads `data/` and `src/` and writes into this folder.

```sh
P=research/correlated-noise-20261101
export PYTHONPATH=src OPENBLAS_NUM_THREADS=1
.venv/bin/python $P/noise_model.py $P/noise-model.json
for p in XO036 XE035 XE036; do .venv/bin/python $P/labeled_sample.py $p $P/sample-$p-null.json; done
for p in XO036 XE035 XE036; do .venv/bin/python $P/labeled_sample.py $p $P/sample-$p-inject.json --inject; done
.venv/bin/python $P/summarize.py $P
.venv/bin/python $P/recalibrate.py $P
python3 $P/check_package.py
```

That's about 12 minutes on two cores. The other experiments write per-row files that aren't
shipped; `experiments_summary.py` condenses them into `experiments-summary.json`:

```sh
for p in XO036 XE036 XE035; do .venv/bin/python $P/grain_null.py $p $P/grain-$p-null.json; done
.venv/bin/python $P/morphology_probe.py $P/morphology-probe.json
.venv/bin/python $P/diff_null.py $P/diff-XE036.json
.venv/bin/python $P/diff_null.py $P/diff-XE036-scalestars.json --scale-stars
.venv/bin/python $P/footprint_gate.py $P/footprint-gate.json
for x in control XE036 XO036 XE035; do .venv/bin/python $P/margin_fix.py $x $P/margin-$x.json; done
.venv/bin/python $P/length_diag.py $P/length-diag.json
.venv/bin/python $P/experiments_summary.py $P
.venv/bin/python $P/montage.py $P/null-passes-montage.png XE035 XE036 XO036
```

`check_package.py` is stdlib-only. It checks the SHA-256 of the 10 inputs (against
`data/manifest.json` where listed) and of every shipped file, then re-reads the headline numbers.

## Not checked

- Injection adds linear signal to photographic density. That's a fair approximation near sky
  level only; bright-trail recall isn't tested.
- True SNR uses grain-only F and σ. Clutter raises the real scatter 1.3–2.2× above that, so the
  recall tables flatter a blind search.
- Trials overlap and share sky, so pass counts aren't independent and carry no binomial error bars.
- The pair-motion check (12″) wasn't combined with the per-plate null; joint red + blue false-pass
  rates are not measured.
- Uniformity U was thresholded on these injections and checked on the same plates; no held-out
  plate. It's untested on trails that cross stars by design.
- The optimizer-stall diagnosis covers three XE036 sites only.
- Subtraction used one global linear flux scale and a Gaussian PSF match; a spatially varying or
  nonlinear density match wasn't tried.
- Replay across platforms: outputs here are from Linux x86-64 with OpenBLAS. The pipeline's own
  committed outputs (macOS arm64) replay here with differences up to 1e-5″.

## Open next questions

1. **Make the length gate measure length.** Profile the fit over length (or start it at 1.5×
   and 0.6× as well and keep the best), then rerun `morphology_probe.py`. Success: features ≥ 2×
   fail at every site while both positives and injected 1× trails still pass.
2. **Check the two useful changes off-sample.** F and the uniformity threshold were measured and
   judged on the same three plates. Adapt `labeled_sample.py` to another frozen plate used as a
   null (the MS4 rehearsal package has seven frozen FITS cutouts). Success: F measured there, and the U
   threshold from these plates, give similar false-pass and recall numbers.
3. **Joint pair null.** Pair red and blue null sites by the predicted displacement and measure how
   often both screens and the 12″ pair-motion test pass together.

## Credit

Plates: Palomar Observatory Sky Survey (POSS-I), 48-inch Oschin Schmidt. The 1955 AN10 images
were found by A. Gnädig and A. Doppler and measured by G. V. Williams (MPEC 1999-N21).
Digitization: STScI Digitized Sky Survey; access: IRSA. Catalogue: ESA Gaia DR3. Ephemeris:
NASA/JPL Horizons (JPL#100). Pipeline, controls and frozen inputs: this repository's pilot.

An AI agent (Claude Code) did this work in a donated session for VibeCodyH. It has had no
qualified scientific review and makes no detection, astrometric or acceptance claim.
