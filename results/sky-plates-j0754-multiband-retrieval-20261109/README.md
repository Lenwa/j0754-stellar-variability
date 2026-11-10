# J0754 multiband retrieval, 9 November 2026

## Result
- ASAS-SN g: 4,284 native measurements, 4,220 passing the existing project's detection screen. Coverage 2017-10-07.4514018 through 2026-11-08.2998711. Re-downloaded from the completed 8 November coordinate-photometry job. The file is byte-identical to the previously retained snapshot; this is not newly acquired post-8-November coverage.
- ZTF r/i: the fresh, unfiltered 3-arcsecond query returned HTTP 200 and zero rows. An independent all-band 10-arcsecond query returned 575 rows, all g, and is byte-identical to the existing target file. Thus no public catalogue r/i lightcurve was retrieved, not a measured nondetection or proof that no images exist.
- ATLAS c/o: not retrieved. The current forced-photometry service requires registration. The historical ATLAS-VAR DR1 detection table requires MAST CasJobs login. No accounts, sign-ins or jobs were created.

## Files
J0754_ASAS-SN_g_native_subset.csv contains every native g row with values and line bytes preserved. asassn_v1_native.csv is the complete original V/g export. J0754_ASAS-SN_g_screening.csv supplies row-level screening decisions without altering the raw data. The empty ZTF r/i response and full 10-arcsecond response are retained alongside timestamped request metadata, hashes and numerical summary.

## Native conventions and limits
ASAS-SN times remain HJD and the supplied UT Date, flux remains mJy, and cameras remain separate. Its export has no explicit image-quality flag. The optional screening uses the existing project rules and is not a full systematic-quality validation. No brightness clipping, zero-point adjustment, nightly averaging, or time-scale conversion was made. ZTF native hjd/mjd, catflags, object IDs, filter codes and all other columns remain unchanged.

The existing ASAS-SN job already extends through 8 November 2026. Obtaining later epochs requires a new recomputation, which has not been submitted. These data do not by themselves identify the fading mechanism.

## Verified public sources
- ASAS-SN result: https://asas-sn.osu.edu/sky-patrol/coordinate/d2c2ad05-c9a6-4362-93ac-523c2584d094
- ASAS-SN export: https://asas-sn.osu.edu/sky-patrol/coordinate/d2c2ad05-c9a6-4362-93ac-523c2584d094/export.csv
- ZTF API: https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html
- ATLAS service and source: https://fallingstar-data.com/forcedphot/ and https://github.com/lukeshingles/atlasserver
- ATLAS-VAR archive: https://archive.stsci.edu/hlsp/atlas-var

## Preservation
The original project snapshots remain unchanged. The new folder is a deliberate retrieval snapshot in the existing project's results area. No maintained source, Git index, commit, branch, or remote publication was changed.
