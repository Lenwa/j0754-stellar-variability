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
