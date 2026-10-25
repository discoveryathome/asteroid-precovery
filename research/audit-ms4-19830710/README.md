# Wrong-epoch audit: 2002 MS4, 1983-07-10

This package explains one recorded wrong-epoch archive response from the
[MS4 rehearsal](../rehearsal-ms4/README.md): the 1983-07-10 exposure. It reads the
frozen inputs in place and changes nothing in `rehearsal-ms4/`. It makes **no detection,
astrometric or scientific-acceptance claim**. The checks are technical and
author-checked; no qualified scientific reviewer has assessed them.

| | Region / plate ID | Plate label | Epoch |
|---|---|---|---|
| Published ([Lowe](https://andrew-lowe.ca/2002ms4.htm)) | S655 / 03YK | not given | 1983-07-10.49060 UTC |
| Returned by the STScI default request | S728 / 00TR | J8626 (SERC-EJ) | 1983-06-14T13:03 |

## Findings

1. **The rejection stands.** The response is a cutout from UK Schmidt plate J8626
   (field S728, plate ID 00TR), taken 26 days before the published exposure.
2. **It's the same physical plate as the exact 1983-06-14 match.** All 146 header
   cards outside the cutout keywords are byte-identical between the 1983-07-10 and
   1983-06-14 files. These include survey, region, plate ID, plate label, DATE-OBS,
   field centre, scan geometry and the full GSC plate-solution coefficients. Only the
   14 cutout keywords differ (corner pixels, local WCS, requested position, data range).
3. **The response follows from position-only selection.** The default request
   carries only a survey, a position and a size. The requested position lies on four
   UK Schmidt fields (655, 656, 727, 728). Field 728 is both the field whose centre is
   nearest (2.874° against 4.131° for the published field 655) and the field with the
   largest margin to the plate edge (1.181° against 0.381°). Both rules also reproduce
   the other three STScI default responses in the package: two exact matches
   (1982-08-09, 1983-06-14) and the other wrong default (1983-07-01). This data can't
   separate the two rules. STScI's [DSS description](https://gsss.stsci.edu/SkySurveys/DSS%20Description.htm)
   says the software narrows candidate plates and picks the best one, without defining
   "best". So this is a consistent explanation, not a documented rule.
4. **The exclusion doesn't rely on the REGION/PLATEID strings.** It uses the plate
   geometry shared with the exact match (finding 2), and it assumes Lowe's 1983-06-14
   position was measured on this plate. At the returned plate's epoch, the object was
   at plate pixel (10217.5, 2405.3): STScI's own conversion (`OBJCTX`/`OBJCTY`) of the
   published 1983-06-14 position. The returned cutout spans X 10952–11102 and
   Y 2423–2573, so the object lies 735 pixels (20.8′) outside it. The frozen Horizons position for that epoch is 23.1′ from the
   cutout centre, against a 3.0′ half-diagonal.

**Consequence.** Under either rule consistent with the data, repeating a default
position-only STScI request at the published 1983-07-10 position returns S728 again.
Acquiring S655/03YK needs a route that selects by plate or field. This audit doesn't
claim that no such route exists.

## Replay

```sh
python3 research/audit-ms4-19830710/check_wrong_epoch.py
```

The check uses only the Python standard library (run with 3.12.3). It is offline and
needs no project environment. It pins the SHA-256 of its twelve inputs and
cross-checks them against `rehearsal-ms4/manifest.json`. The four STScI FITS hashes
are also checked against their `.provenance.json` records, and the requested
positions are read from the pinned provenance URLs. The expected
output is in `check.log`.

Failure controls run on every invocation. A flipped byte fails the hash check. A
header patched to S655 / 03YK / 1983-07-10 is accepted by the identity predicate, so
the predicate isn't vacuous. The same patched header breaks the same-plate comparison.
An object placed inside the cutout gives zero gap. Separately, two scratch copies made
the script exit with `AssertionError: ... bytes differ from pinned base`: one with a
byte flipped in the 1983-07-10 FITS file, and one with the last digit of the requested
RA changed in its provenance URL.

## Assumptions and what was not checked

- The blue (SERC-EJ) centres of fields 655 and 656 were not observed. They are taken
  from IRSA's red (ER) plates on the same UK Schmidt field grid. Blue and red centres
  agree within 0.41′ for both fields where both are in the package (727, 728). Edge
  margins use the nominal 355 mm plate and the returned plate's scale for every field.
- STScI's actual selection algorithm is not documented. It is inferred from four probes.
- `rehearsal-ms4/verify_package.py`, the AN10 baseline and the unit tests were **not
  run** in this session: no project environment was installed and the session wasn't
  sandboxed. The hashes were checked independently instead.
- No image pixels were inspected and no source was measured.
- For the sibling record 1983-07-01, only the selection prediction was tested. Its
  same-plate and exclusion checks were not run.
- The published S655/03YK plate was not acquired. No archive was queried to see
  whether STScI or IRSA hold it.
- The timing of the published 1983-07-10 observation was not reviewed.
- The Horizons solution (JPL #25) already includes the 1983 observations, so the
  Horizons cross-check is not independent. The pixel-geometry exclusion doesn't use it.

## Packaging finding (not changed here)

`rehearsal-ms4/exposure-inventory.json` gives `observer_image` values
`observer-19830701.gif` and `observer-19830710.gif` for the two wrong-default records.
Neither file is in `rehearsal-ms4/data/raw/` or its `manifest.json`. The frozen
`observer-report.html` links both (`http://andrew-lowe.ca/2002ms4-19830701.gif`,
`http://andrew-lowe.ca/2002ms4-19830710.gif`). The rehearsal package is immutable, so
this is reported rather than fixed.

## Open next questions

1. **Acquire S655/03YK.** Find an archive route that selects by plate or field rather
   than by position alone, and record the request, the returned identity and the
   headers. Success: a returned header with `REGION = 'S655'`, `PLATEID = '03YK'` and a
   DATE-OBS on 1983-07-10.
2. The same for S656/02C3 (1983-07-01), after running this package's same-plate and
   exclusion checks on that record.
3. Distinguish the two selection rules with a probe position where they disagree.
   This is only needed if future searches rely on predicting the STScI default.

## Credit

Andrew Lowe identified MS4 on these plates and published the precovery report. The
photographic plates are from the UK Schmidt Telescope (AAO/ROE copyright, per the
headers), digitized and served by STScI/MAST as the Digitized Sky Survey. IRSA supplied
the coverage tables, and NASA/JPL Horizons the frozen ephemeris. This audit builds
directly on the `rehearsal-ms4` package, which acquired and froze every input it reads.

## Session record

- Base commit `9795ae8e87ba525d57b7cb6a135abb8a9eb34c56`.
- Task `audit-one-rejected-match` (track `archive-searches`, priority_revision 2), drawn
  with `python3 tools/pick_task.py --seed 06d8f5aa0a9b9d47d32fad91ba2a9228`.
- The card has two wrong-default records. This one is the second in `stsci-probes.json`
  order, selected by the seed's third 8-hex group (`d32fad91` mod 2 = 1).
- Work done by an AI agent (Claude Code) in a donated session. Author-checked only.
