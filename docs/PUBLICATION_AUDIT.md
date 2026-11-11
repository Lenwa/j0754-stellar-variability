# Public release audit

## Radio/UV analysis follow-up, 10 November 2026

Only source, tests and derived scientific results/plots were added; no new
archive files or account state were required. All 159 tests pass. Gitleaks and
private account-path/email checks found no issues. Radio, UVOT exposure-map
and GALEX count plots were visually reviewed. Source, protocol and input hashes
are retained. Scientific nulls and the UVOT coverage failure are explicit.

## Seven-archive inventory follow-up, 10 November 2026

The public query responses, 460 GALEX aperture photon events and seven small
GALEX/Swift/VLASS FITS products were added with exact retrieval receipts. All
154 tests pass. Gitleaks reports zero findings; private local-account path and
personal-email pattern checks found no matches. Radio cutout WCS positions
place the target inside all three finite image arrays. These checks verify
data access and packaging, not a radio detection or eROSITA association.

## Expanded flare-search follow-up, 10 November 2026

The new source, search products and five public MAST pixel cutouts were checked
with the same Gitleaks configuration (zero findings), plus private local-path
and personal-email pattern checks (no matches). All 151 numerical tests pass.
New plots were visually reviewed. Download receipts contain public archive
URLs, relative paths and hashes. The new data are calibrated 7x7 cutouts totaling
31.75 MB; the original private workspace remains excluded. The new findings are
exploratory source-centered flare candidates, with rejected hits retained.

Prepared 10 November 2026 for publication of the astronomy work only.

## Scope

The release exports the astronomy source component at original source commit
`aefff19` into a fresh standalone Git repository. The unrelated parent project,
its branches, prior commit metadata, conversations, databases, account settings,
credentials and private operational receipts are not copied into this history.
The public commit identity uses a GitHub no-reply email address.

The release includes the research source and tests, specialist documentation,
the updated illustrated findings report, scientific tables and structured
outputs, analysis figures and selected compact archive products. Large image
and calibration files are documented in the native-data manifest. Public source
retrieval URLs are retained separately from private job access routes.

## Checks performed

- Gitleaks 8.30.1 was downloaded from its official release with its published
  checksum verified. The directory scan passes using the checked-in config.
- Initial generic-key findings were individually reviewed. They were SHA-256
  hashes of public API documentation and explicit TESS data-product labels.
  Only those narrow scientific forms are allowlisted; credential rules remain
  enabled. The allowlist is visible in `.gitleaks.toml`.
- Text, binary products and the members of DOCX/NPZ archives were checked for
  local home/account paths and personal email addresses. No remaining matching
  private identifiers were found. Local paths were removed from the export;
  original observations and the private source checkout were preserved.
- DOCX creator, last-modified-by, company and manager metadata were cleared.
  Analysis PNGs were re-encoded without original ancillary metadata. A contact
  sheet of the 51 exported analysis PNGs was visually reviewed; they are
  astronomical plots and cutouts, not private browser or desktop captures.
- No archived conversation, user instruction/profile files, browser HTML state,
  job access receipts, copied third-party paper pages or original Git objects
  were intentionally included.
- The exported source passes all 147 numerical tests. A shorter-window rotation
  analysis also runs against the included data and reproduces the reported
  measurements. Source modifications concern workspace paths and public
  packaging; scientific fitting algorithms were retained.

The checked-in audit records what was checked, not an absolute guarantee against
every possible sensitive string or an independent scientific validation. Earlier
result files keep their original scientific provenance hashes; files edited for
path/privacy cleanup will differ from their original stored byte hashes.
`PUBLIC_FILES_SHA256.json` records the public snapshot's file hashes.

## Rechecking locally

```bash
gitleaks dir . --config .gitleaks.toml --redact
gitleaks git . --config .gitleaks.toml --redact
python -B -m unittest discover -s tools/sky-plates -p 'test_*.py'
```

Report a suspected credential through GitHub's private security reporting route
where available, rather than pasting it into a public issue. Scientific errors
and reproducibility issues can be filed in ordinary repository issues.

## 11 November GALEX/eROSITA follow-up

Public detector calibration, photon/background queries and the target subset of
the official Gaia counterpart catalogue are included with request receipts.
The full 855 MB counterpart catalogue was streamed, not copied into the repo.
The original duplicate-column header and unused shutter API errors are retained.
The corrected analysis, figure, protocol and 166-test suite were checked locally.

## 11 November common flare comparison and evidence ledger

Existing public native pixels and response curves were reused without changes.
The four-event comparison, matched image controls, recovery trials, protocol,
two figures and conditional rate illustration are included. A read-only review
caught unmatched image-control sampling; the corrected controls match exact
retained before/event/after timestamp patterns. A numerical boundary artifact
in the new offset fitter was corrected and pinned by a real-template regression.
Apparent hourly drift and substantial baseline uncertainty remain in the result.

The curated 36-row evidence ledger links 44 local source files and inventories
37 retained result families, distinguishing instruments from reused reductions
and physical hypotheses. A second read-only review found no important scientific
overclaim; catalogue/decade/thousands formatting was corrected. All 176 tests
passed. Scientific figures were visually inspected. Built-in browser security
blocked a local-file HTML preview; the HTML structure and source/image links
were checked offline. No claim of browser-rendered visual QA is made.

## 11 November giant comparison screening

The full 86-row Olah catalogue and its original classification/rotation tables,
five-star physical shortlist, Gaia/MAST responses, 12 reduced Kepler FITS and
seven explicitly unreduced sibling TESS/K2 FITS are retained with public request
receipts. The unavailable Dixon sample and timed-out ZTF queries remain failures,
not screened stars or nondetections. No account pages, cookies or credentials
were retained. The protected journal response was not saved.

A fresh read-only review independently reproduced the physical ranking, all
24 SAP/PDC period maxima and powers, and the historical annual summaries.
Review findings were corrected: cached catalogue reads no longer overwrite
enriched results; published period pairs and their original text are preserved;
plot lines stop at missing annual bins. The source hash manifest was regenerated
after those edits. All 188 component tests passed, including 12 new safeguards.
Scientific plots were inspected. Full-repository Gitleaks and targeted private
path/contact checks were clean. Native survey files remain unchanged.

## 11 November Tang sibling analysis

The seven retained sibling FITS are now analyzed without altering their native
bytes. One additional compact public K2SC FITS is included with its retrieval
receipt and SHA-256. The EVEREST response exceeded the existing 12 MB retention
bound; its failed-size receipt remains visible without retaining the product.
All previous native-manifest and public-request records are preserved.

The initial pre-fit plan snapshot is unchanged. The protocol separates later
centroid/K2SC follow-up and post-review QLP spacecraft-orbit reference offsets.
A fresh read-only review reproduced the native fits, sector holdouts and
corrected zero-excursion result. Native-cadence bridging, local noise, numerical
ties and flat-null scans were corrected with regression tests. The report
distinguishes fitted periods from confirmed orbits, alternative reductions
from new observations, and 2021 modulation from historical fading.

All 205 component tests pass, including 17 new sibling safeguards. Four figures
were visually inspected. New source/native hashes match; the original plan
snapshot hash matches. Full-repository Gitleaks scanned 198.92 MB without
findings, and targeted private path/contact/key checks were clean. Public
request receipts contain archive URLs and catalogue identifiers, with no
account pages or credentials. The publication manifest hashes canonical Git
blobs; native FITS retain their exact archive bytes.

## 11 November J0830 direct-pixel check

Three public 21x21 TESSCut calibrated FITS and a Gaia field response are included
with original byte hashes and request receipts. The geometry-selected control
file is preserved byte-for-byte using a scoped Git attribute. Previous manifest
and public-request records remain unchanged. No private account records or
credentials are present in the new source, figures or public archive products.

The pixel result is explicitly partial: primary/large aperture modulation
survives a bright-control correction; S46 is spatially consistent, while S44–45
remain mixed with image changes. The narrow registered aperture and unstable
faint control are retained. One motion-aware check is labelled as later work,
and no further correction family was added. The report connects this boundary
to the main fading question instead of claiming orbital or causal confirmation.

Fresh read-only review verified native hashes/time headers, both controls'
eligibility, common cadences and independent least-squares results. All 212
component tests pass, including seven new numerical safeguards. Four figures
were visually inspected; eight input/source hashes match. Full-repository
Gitleaks scanned 295.32 MB without findings; targeted private path/contact/key
checks were clean. The canonical Git-blob manifest is updated with this release.
