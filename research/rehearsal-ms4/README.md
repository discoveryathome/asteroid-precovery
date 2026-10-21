# Second-control acquisition: 2002 MS4

This local contributor rehearsal acquired **five archival cutouts for a second documented object**, and retained **two wrong-epoch archive responses**. This is an acquisition and provenance contribution. It does not yet reproduce calibrated astrometry for this object or add observations to an orbit.

[Five-exposure preview](five-exposure-preview.png) compares the investigator's images with the retrieved pixels. Orange circles mark published positions; blue crosses mark a frozen current-orbit prediction. Neither marker is a newly measured source position.

| Published epoch | Retrieved region / plate ID | Current assessment |
|---|---|---|
| 1954-04-08 | XE740 / 0FGQ | Visually matching field; published alias se0741 and one-day header discrepancy require review |
| 1955-04-16 | XE741 / 0FGR | Visually matching field; published alias se0742 and one-day header discrepancy require review |
| 1982-08-09 | S727 / 044C | Exact published region and plate ID |
| 1983-06-14 | S728 / 00TR | Exact published region and plate ID |
| 1983-07-01 | Default request returned S728 / 00TR | Wrong epoch; requested S656 / 02C3 remains outstanding |
| 1983-07-10 | Default request returned S728 / 00TR | Wrong epoch; requested S655 / 03YK remains outstanding |
| 1992-06-23 | ER729 / A2JL | Exact published region and plate ID |

The primary historical source is [Andrew Lowe's precovery report](https://andrew-lowe.ca/2002ms4.htm). The [MPC discovery announcement](https://www.minorplanetcenter.net/iau/mpec/K02/K02W27.html) and a frozen JPL SBDB response document the object. Its discovery date is 2002-06-18. [The investigator's index](https://andrew-lowe.ca/tno-precoveries.htm) supplies 23 linked case reports including this one; the other 22 are **untriaged leads**, not verified ready tasks. All source snapshots have acquisition URLs and checksums.

## Reproduce the local checks

From the repository root, using its installed environment:

```sh
.venv/bin/python research/rehearsal-ms4/verify_package.py
MPLCONFIGDIR=.replay/mpl .venv/bin/python research/rehearsal-ms4/preview_all.py
```

The first command is offline and checks file hashes, exact exposure IDs, wrong-epoch rejection, and native cutout coordinate preservation. The second rebuilds the visual comparison from frozen data. A new Matplotlib configuration directory may spend time building its initial font cache; that happened in the rehearsal. Use an existing cache for quick subsequent plotting. The package manifest hashes the delivered image; rendering in another environment can change PNG bytes without changing the input pixels.

The baseline remains separately reproducible:

```sh
PYTHONPATH=src .venv/bin/python -m precovery.reproduce .replay/output
PYTHONPATH=src .venv/bin/python -m unittest discover -s tests -v
```

The baseline tests exercise the original AN10 pipeline. They are not validation of MS4 measurements, which have not been made.

## Evidence and limits

`exposure-inventory.json` preserves seven published epochs, exact returned header identities, timing differences, and review status. `coverage.json` and `additional-coverage.json` retain the IRSA results returned for all seven sky positions. This is not an exhaustive search of all providers or every overlapping historical exposure. The four STScI requests are in `stsci-probes.json`.

`data/raw/` holds native FITS cutouts, original parent headers where range extraction was used, observer comparison images, archive response tables, JPL predictions, and acquisition provenance. Full parent plates were intentionally not downloaded; range and cutout checksums do not masquerade as full-plate checksums. `manifest.json` pins delivered files, including original code used in this contribution.

`versions.json` pins JPL solution 25, its 2026-03-14 solution date, the Horizons/covariance epoch and the distinct SBDB elements epoch. The present orbit already contains historical observations. Residuals in `prediction-audit.json` compare **published coordinates** with current predictions, not new image astrometry or independent orbit validation. Pre-1962 Horizons output is UT1; the original report labels its dates UTC. Header date and timing discrepancies remain visible and unresolved.

No reference-star catalogue for these fields has been acquired. The original pilot's Gaia subset covers another field and must not be reused. Approximate two-arcsecond motion during the oldest exposures makes the existing long-trail classifier unsuitable without separate validation. No screening threshold has been weakened to accommodate this object.

See `review.json` for explicit open obligations. The full archive task is partial: five cutouts and a documented second case are ready for follow-up; two requested exposures and a complete identity/timing audit remain outstanding. Newly contributed observation-arc extension is zero and orbital-uncertainty improvement is unmeasured.

## Acquisition scripts

`acquire_case.py`, `inspect_headers.py`, `additional_coverage.py`, `stsci_probe.py`, `native_cutout.py`, `fetch_predictions.py`, and `build_inventory.py` preserve the executed retrieval approach. Existing source snapshots are immutable. Range extraction refuses to replace an existing FITS cutout. Network reruns belong in a separate copy with new raw-data storage; live replies and temporary archive URLs may differ. The immutable package is the reproducibility target.

The initial IRSA cutout request returned HTTP 502; its invocation and traceback are preserved in the session transcript. The underlying cause was not diagnosed. Existing range extraction recovered the two old Palomar cutouts; a URL-parameterized variant recovered ER729. STScI's default blue survey returned two exact plate matches and two incorrect defaults, which were retained rather than accepted.

## Credit and rights

Credit Andrew Lowe for the published precovery work and reference illustrations, Caltech/Palomar and the UK Schmidt Telescope for the photographic observations, STScI and IRSA for DSS access, the MPC for historical records, and NASA/JPL for ephemerides and orbital data. External reports, illustrations, scans and metadata retain their respective rights; the project's CC BY license applies only to original prose, and MIT only to original code. No external publication, astrometric submission or new-observation announcement occurred.
