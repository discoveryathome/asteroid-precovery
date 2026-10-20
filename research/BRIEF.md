# Research brief — 2026-10-20

The user approved this local replication scope and local execution on 2026-10-20. The goal is an auditable archival-image workflow, first validated on a documented result. It is not to rediscover a scientifically unknown observation by optimizing residuals.

## Claim audit and starting frontier

| Claim | Status and evidence |
|---|---|
| Historical images can contain observations predating recognized discovery | Confirmed for this case by [MPEC 1999-N21](https://www.minorplanetcenter.net/iau/mpec/J99/J99N21.html) and [JPL](https://cneos.jpl.nasa.gov/news/news018a.html). |
| MOST supports precovery searches | Confirmed by [IRSA instructions](https://irsa.ipac.caltech.edu/applications/MOST/instructions.html); its output identifies coverage, not a verified asteroid detection. |
| Suitable public pixels can be retrieved at manageable size | Confirmed experimentally here: two native DSS cutouts, independently cross-checked against the IRSA cutout service. |
| These controls improve the present orbit | Not claimed: they are already in the modern JPL fit. A future before/after refit is required to quantify incremental improvement. |
| A small residual is sufficient for discovery | Rejected by project design and by the artifact control, which lies close to the prediction but fails morphology/significance checks. |

MPEC 1999-N21 publishes four trail-end positions for January 26, 1955, from 45-minute red and 12-minute blue exposures. Its account flags uncertain trail boundaries and minute-rounded timing. The original reduction used USNO-A2.0; this pilot uses a separately frozen Gaia DR3 subset. JPL supplies contemporary context and an image of the trail. Neither source is presented as endorsement of this new reduction.

The frontier is a **case-specific local replication**, not a claim to outperform established precovery systems. The already achieved historical identification is the positive control. The new engineering question is whether strangers can cheaply reproduce the image evidence and diagnose failure cases.

## MOST and archive choices

[MOST](https://irsa.ipac.caltech.edu/applications/MOST/) combines moving-object ephemerides with coverage for selected IRSA collections, including WISE/NEOWISE, Spitzer, 2MASS, PTF, ZTF, SOFIA and SPHEREx. It accepts object identifiers or orbital elements. Its [API](https://irsa.ipac.caltech.edu/applications/MOST/MOSTProgramInterface.html) and [official Python tutorial](https://caltech-ipac.github.io/irsa-tutorials/most-queries/) provide automation paths. Query bounded time intervals, retain metadata and exact URLs promptly, and avoid bulk image tarballs until the hit list is understood. Some collections include restricted data; use public records only. An intersection is a search lead. The pilot's 1955 DSS plates are obtained through IRSA Atlas/Cutouts, not through MOST, whose supported datasets differ.

| Public access route | Useful products | Main limitation for precovery |
|---|---|---|
| [IRSA DSS Atlas](https://irsa.ipac.caltech.edu/data/DSS/) and [STScI DSS](https://archive.stsci.edu/dss/) | Photographic plates, polynomial WCS, overlapping fields | Plate selection and historical timing require checks; scans are not calibrated photon counts. |
| [CADC SSOIS](https://www3.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/ssois/documentation.html) | Moving-object coverage across many observatories, links to providers | Access and calibration differ by archive; a footprint does not guarantee a detection. |
| [NOIRLab SIA](https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/DataAccessInterfaces/SimpleImageAccessSIA/SimpleImageAccessSIA.html) | Public image metadata and cutouts, including single-epoch survey products | Choose individual exposures; coadds can suppress or distort moving objects. |
| [IRSA ZTF access](https://ztf.ipac.caltech.edu/page/dr8) | Single-exposure products and server-side cutouts | Verify release, exposure ID, public availability, mask, reference image and timing. |

The local experiment found 12 overlapping IRSA DSS products. The nearest default STScI request selected a 1953 plate. IRSA returned both target 1955 fields. The STScI plate-specific service returned a server-side temporary-directory error; that failure is retained in the acquisition log. The full parent plates would be hundreds of MB to over 1 GB. Bounded HTTP ranges transferred about 10 MB and 28 MB and produced sub-MB native cutouts. Parent-file hashes are explicitly unavailable; retrieved-byte and cutout hashes are recorded instead.

## Position prediction, covariance and calibration

[JPL Horizons](https://ssd.jpl.nasa.gov/horizons/manual.html) supplies topocentric astrometric ICRF coordinates and formal uncertainty ellipses. The saved response pins JPL#100, its solution date, osculating epoch, DE441, observer 675 and the full request. Three-sigma semiaxes become one-sigma covariance by division by three before squaring. Formal covariance can be optimistic under systematics or nonlinear extrapolation. Unavailable uncertainty must not become zero. Pre-1962 output uses UT1.

[OrbFit](https://adams.dm.unipi.it/orbfit/) is a candidate for controlled orbit propagation/refitting and covariance experiments. It is not installed or validated here. Any future refit must freeze observations, force model, weighting/debiasing, solver version, parameters and comparison epoch; include/exclude the precovery while keeping everything else comparable. Osculating elements alone are not a substitute for the full perturbation model.

[Gaia DR3](https://www.cosmos.esa.int/web/gaia/dr3) uses reference epoch 2016.0 and ICRS astrometry. The pilot freezes 191 catalogue rows and query metadata, applies proper motions to 1955, filters problematic stars, then fits a local affine correction to the native WCS using 47 inliers per plate. Leave-one-out stellar residuals and a 0.5-arcsecond systematic floor enter the evidence. Linear propagation and missing full star-error correlations are limitations, not precision guarantees.

## MPC reporting and review

[MPC reporting guidance](https://docs.minorplanetcenter.net/mpc-ops-docs/astrometry/reporting-observations/) prefers ADES and distinguishes midpoint timing from trail-end timing. One must not report a trail midpoint and both ends as separate measures of the same exposure. [ADES documentation](https://docs.minorplanetcenter.net/mpc-ops-docs/observations/ades-format/) describes XML/PSV validation; [allowed values](https://docs.minorplanetcenter.net/mpc-ops-docs/observations/valid-ades-values/) define catalogue and instrument codes. A future reviewed export would need the accepted object designation, station, time convention, frame, position, mode, catalogue, uncertainties, reduction provenance and appropriate observers/measurers. Unknown photometry must not be invented. This pilot deliberately has no submission/export path.

For this case, archive DATE-OBS values conflict with the original publication. Published trail-end times supply a visible, documented replication override with a conservative 30-second timing floor. Original plate logs remain an expert-review obligation. The measured inter-exposure displacement, source morphology and catalogue checks support reproduction; they do not resolve every possible confuser or establish calibrated discovery significance.

## Metrics and scope limits

Report positive/negative control outcomes, calibration residuals, trail-versus-point discrimination, pair motion and complete provenance. The historical start-date extension is 16057.60 days; the newly contributed extension is zero. Orbital-uncertainty improvement remains null. A useful future metric is the ratio of propagated sky-ellipse areas or covariance eigenvalues at the same evaluation epoch after consistent with/without fits; the determinant of a mixed-unit covariance alone is hard to interpret.

Two exposures of one target cannot estimate survey recall, completeness or a reliable false-alarm rate. Expand only with documented independent controls. Expert review, not compute volume, is the limiting resource.
