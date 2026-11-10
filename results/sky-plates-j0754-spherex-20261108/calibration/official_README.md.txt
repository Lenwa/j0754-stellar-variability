# Spectrophotometry Tool instructions
```
photometry-tool perform-photometry --help
Usage: photometry-tool perform-photometry [OPTIONS]

Options:
  --result-dir DIRECTORY          Directory in which all the intended UWS
                                  results of the job will be deposited.
                                  Defaults to '.', the current working
                                  directory. The tool will create a symlink
                                  'result' in that directory pointing to the
                                  'principal' result  [default: .]
  --err-dir DIRECTORY             Directory in which an error document for the
                                  job, if any, will be deposited by the tool.
                                  Defaults to “.”, the current working
                                  directory.  [default: .]
  --runid TEXT                    Allows specifying a UWS 'runid', essentially
                                  a job name, which may be used in log
                                  messages or other metadata. Defaults to a
                                  short string indicating the specific science
                                  tool being run (e.g., 'photo') followed by
                                  an ISO-8601 timestamp, e.g.,
                                  'photo-20250307T172132'. This is not the
                                  same as the UWS 'Job Identifier' passed via
                                  --jobid. [optional]
  --jobid TEXT                    Allows the UWS system to pass in the opaque
                                  UWS 'Job Identifier' by which the job is
                                  referenceable in the REST API. Required; no
                                  default is supported. This is not the same
                                  as the UWS 'Run ID' which is set by the user
                                  and passed to the tool via the --runid
                                  option.  [required]
  --tmp-dir DIRECTORY             Directory available for use as temporary
                                  scratch space for the job. Defaults to “.”,
                                  the current working directory. When
                                  executing in the UWS system, tmp-dir should
                                  never be the same as result-dir, unless it's
                                  deliberately intended for temporary files to
                                  be returned to the user as results.
                                  [default: .]
  --log-dir DIRECTORY             Directory into which any log files from the
                                  job will be deposited.Defaults to “.”, the
                                  current working directory. When executing in
                                  the UWS system, log-dir should never be the
                                  same as result-dir, unless it's deliberately
                                  intended for any logs to be returned to the
                                  user as results.
  --fake [normal|error-short|error-detail|error-signal]
                                  Run a simulated job. Modes: "normal" =
                                  simulate successful run, "error-short" =
                                  create errorSummary.yaml, "error-detail" =
                                  create errorSummary.yaml and symlink "error"
                                  to actual error file, "error-signal" = exit
                                  with SIGINT. [optional flag]
  --point TEXT                    Point-like target in format ra,dec[,name].
                                  [One of  --point OR --elliptical is
                                  required.]
  --elliptical TEXT               Elliptical target
                                  ra,dec,sersic,ratio,pa,effrad[,name].
                                  
                                  [One of  --point OR --elliptical is
                                  required.]
  --time TEXT                     Time range as min,max (MJD). [optional]
  --detector TEXT                 Detector label: D1, D2, D3, D4, D5, or D6
                                  [optional]
  --bkg-region INTEGER            Background region in pixels (default: 15)
                                  [optional]
  --fit-positions / --no-fit-positions
                                  Allow positions to float in Tractor fit.
                                  (Not implemented.) [optional flag]
  --verbose / --no-verbose        Turns on additional verbosity in a job's
                                  logging; should not affect the output
                                  visible to the user. Defaults to --no-
                                  verbose. Translated within the science tools
                                  to logging control options.  [default: no-
                                  verbose]
  --use-simulator-bandpass / --no-use-simulator-bandpass
                                  If True, simulator bandpass model is used to
                                  populate measurement wavelength and
                                  bandwidth. Otherwise, spectral_wcs
                                  calibration product is used.  [default: no-
                                  use-simulator-bandpass]
  --add-image-ids / --no-add-image-ids
                                  If True, image ids are included in the
                                  resulting VOTable. Defaults to False.
                                  [default: no-add-image-ids]
```


# Environment variables
The following environment variables should be set before running `photometry-tool perform-photometry`
- SPXTOOL_BUTLER - Butler datastore root
- SPXTOOL_NCORES: Maximum number of cores that the script is allowed to use for its processing. The
  tool is allowed to take advantage of this either by spawning subprocesses or by multithreading
- SPXPHOTO_INSTRUMENT - Instrument name: `spx` or `simu`. Defaults to `spx`.
- SPXPHOTO_L2_DATASET_TYPE: Dataset type for level2 images. Defaults to `level2`.
- SPXPHOTO_INPUT_COLL - Input collections that include level2 images used by photometry. Defaults to `IRSA-qr1-visible`.
- SPXPHOTO_CALIBRATION_COLL - Input collections that include calibrations used by photometry. Defaults to `IRSA-qr1-tools`.
- SPXPHOTO_QGRAPH_TIMEOUT - Qraph (execution plan) build timeout in seconds. Defaults to 1 hour.
- SPXPHOTO_L3_TIMEOUT - Pipeline execution timeout in seconds. Defaults to 12 hours.
- SPXPHOTO_EXTRA_RESOURCE: Path to IRSA-provided file containing service descriptor XML content.

# Error codes:
- 11 → Butler does not exist
- 12 → Butler exists but there is no level2 data in it.
- 13 → No level2 image matches for the provided sources.
- 14 → Subprocess executing one of the pipelines failed.
- 15 → Unexpected error when running photometry.
- 16 → No combined all sky catalog found.
- 17 → Failed to add service descriptor to VOTable XML content.
- 18 → Unexpected failure during VOTable file creation or modification.
