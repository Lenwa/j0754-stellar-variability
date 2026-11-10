## SPHEREx Pipelines 6.5 (2025-12-11)


### New Features

- Photometry tool v13: Major version of the photometry tool. Changes include CLI, butler-related changes, and bug fixes. [[#SP-1661](https://jira.ipac.caltech.edu/browse/SP-1661)]
- Photometry tool: Add `obs_publisher_did` column to photometry tool VOTable output. [[#SP-1772](https://jira.ipac.caltech.edu/browse/SP-1772)]
- Introduced an afterburner script to update the spectral WCS of Level 2 products, writing the results to a new collection. [[#SP-1823](https://jira.ipac.caltech.edu/browse/SP-1823)]
- Mosaic tool CLI and sample docker delivery [[#SP-1850](https://jira.ipac.caltech.edu/browse/SP-1850)]
- Use spectral WCS calibration product instead of simulator bandpass model to populate measurement wavelength and bandwidth for photometry tool v13.1 release. [[#SP-1862](https://jira.ipac.caltech.edu/browse/SP-1862)]
- Added github actions workflow to build and deploy Docker container for spherex-pipelines. [[#SP-1902](https://jira.ipac.caltech.edu/browse/SP-1902)]


### Changes
- Ingested updated `bad pixel map` from 10/19/25 instrument delivery. [[#SP-1857](https://jira.ipac.caltech.edu/browse/SP-1857)]
- Ingested `SIP` from qr-2 to use as initial SIP guess. [[#SP-1469](https://jira.ipac.caltech.edu/browse/SP-1469)]

- Identification algorithm and flagging for source blooming, cosmic-ray snowball and halo, (BSH Flagging) artifacts: Aded new flag definitions and optional diagnostic output (to be used for validation). [[#SP-1582](https://jira.ipac.caltech.edu/browse/SP-1582)]
- Created a dedicated module for source bloom and cosmic-ray halo flagging and added it to the beginning of L2a pipeline; preserves L1 flags as int16 and L2 flags as int32. [[#SP-1886](https://jira.ipac.caltech.edu/browse/SP-1886)]
- Apply BSH flagging and update persistence current only to active pixel area of detectors. Persistence should ignore bad pixels. Streamlined BSH status information in header. [[#SP-1893](https://jira.ipac.caltech.edu/browse/SP-1893)]


- Update Persistence: update persistence current with bloom, snowball, halo flagging; this includes distinguishing cosmic rays from bright sources.  `SSDC–MS–012 v2.2.1` [[#SP-1581](https://jira.ipac.caltech.edu/browse/SP-1581)]
- Persistence Module Piecewise Function: Added new piecewise function, modified flagging rules, and 
  updated persistence parameters for R6.5. [[#SP-1882](https://jira.ipac.caltech.edu/browse/SP-1882)]
- Update persistence current should be applied to only active area. [[#SP-1890](https://jira.ipac.caltech.edu/browse/SP-1890)]


- Solve for off-center detector position-angles (PA) using boresight and PA from planned observation. [[#SP-1549](https://jira.ipac.caltech.edu/browse/SP-1549)]
- R6.5 astrometry improvements: Added sytematic shift corrections to coarse astro initial guess. [[#SP-1868](https://jira.ipac.caltech.edu/browse/SP-1868)]
- Astrometry diagnostics: diagnostic capability built into the Fine Astrometry module for routine data collection during pipeline execution. [[#SP-1831](https://jira.ipac.caltech.edu/browse/SP-1831)]
- Added `FINAST` quality metric: `Qs_nmatch_var` determined from the number of matched sources in each quadrant of the image. [[#SP-1897](https://jira.ipac.caltech.edu/browse/SP-1897)]

- Dichroic Beam-Splitter (DBS) ghost masking: Added functionality to mask `main` and `edge` DBS ghosts for detectors 3-6 inclusive. `SSDC-MS-019 v3.3`. [[#SP-1858](https://jira.ipac.caltech.edu/browse/SP-1858)]
- Validtion DQA check for corners of `GHOST_EXT` region of interest (ROI), logs a failure if there are no sources in the corners. [[#SP-1878](https://jira.ipac.caltech.edu/browse/SP-1878)]


- Rewrite Query Reference Catalog task to use planned observations as inputs instead of L1 images. This mitigates loss of coverage of srccat when there are missing L1 images. [[#SP-1722](https://jira.ipac.caltech.edu/browse/SP-1722)]
- Keep only necessary columns as float64 in srccat (ra, dec, times), others may use float32. [[#SP-1885](https://jira.ipac.caltech.edu/browse/SP-1885)]
- Added a DQA check for srccat coverage of least 250 pixel padding past image active area region, on all sides. [[#SP-1892](https://jira.ipac.caltech.edu/browse/SP-1892)]
- Query Reference Catalog Photometry Flag: Added `for_l3phot` flag to filter down galactic plane sources
  photometered in L3 without affecting the source mask. [[#SP-1881](https://jira.ipac.caltech.edu/browse/SP-1881)]
- Query Custom Catalog Task: Created a new Query Reference Catalog-like task to use input reference catalogs with custom sky coverage for deep-field. [[#SP-1895](https://jira.ipac.caltech.edu/browse/SP-1895)]

- Photometry Flagging: Added `BLOOM`, `SNOWBALL`, `HALO`, and `SATELLITE_HALO` to the list of flags to mask when
  estimating the diffuse background, the list of flags for Tractor to ignore, and the list of bad photometry
  flags. [[#SP-1876](https://jira.ipac.caltech.edu/browse/SP-1876)]


- L3 outputs: Based on the refcat version in the QRC srccat outputs, include file-level metadata in the Level 3 output Parquet files documenting the refcat version used. [[#SP-1810](https://jira.ipac.caltech.edu/browse/SP-1810)]


- Decoupling Butler for Photometry Tool: When extracting products from the main butler to a temporary one, 
  a new process is launched to avoid storage class mismatches between them. [[#SP-1845](https://jira.ipac.caltech.edu/browse/SP-1845)]




### Bug Fixes

- Fixed `SPHERExImage` layer order on read to match the order of the corresponding FITS extensions. [[#SP-1826](https://jira.ipac.caltech.edu/browse/SP-1826)]
- Forced Photometry Flagging: Fixed source region mask to exclude non-source
  pixels within the bounding box. [[#SP-1833](https://jira.ipac.caltech.edu/browse/SP-1833)]
- Nonfunc flag count set to integer from numpy integer for pyyaml compatibility. [[#SP-1839](https://jira.ipac.caltech.edu/browse/SP-1839)]
- Repeating UUIDs for CAOM Plane Table Fix: Fixes an indentation error in the logic that caused the same UUID 
  to be used as the obsid for multiple plane ids. [[#SP-1846](https://jira.ipac.caltech.edu/browse/SP-1846)]
- Photometry tool v13 bugfix: Subprocess hangs when all sources have no coverage as a consequence of multiprocessing queue blocking. [[#SP-1847](https://jira.ipac.caltech.edu/browse/SP-1847)]
- Photometry tool v13: `obs_pub`lisher_did` should reflect qr2 by default. [[#SP-1848](https://jira.ipac.caltech.edu/browse/SP-1848)]
- Addresses the suppression of the `local_bkg_flg` column in the Firefly display of Spectrophotometry Tool outputs. [[#SP-1856](https://jira.ipac.caltech.edu/browse/SP-1856)]
- Perform Forced Photometry Bugfix: Added persistence flag to bad pixel mask and mask of pixels for Tractor to
  ignore. [[#SP-1874](https://jira.ipac.caltech.edu/browse/SP-1874)]
- Astometry bug fixes r6.5. [[#SP-1891](https://jira.ipac.caltech.edu/browse/SP-1891)]


### Enhancements

- Estimate Diffuse Background Optimization: Reduced memory allocated by the run
  method by adding chunked spline evaluation and removing outdated null
  background objects. [[#SP-1818](https://jira.ipac.caltech.edu/browse/SP-1818)]
- Level 3a Optimization: Reduced memory allocated from reading parquet files by enabling column filtering and
  initializing `astropy.table.Table` objects with shallow copies. [[#SP-1819](https://jira.ipac.caltech.edu/browse/SP-1819)]
- Level 3 Photometry Improvements: Added more flexible pass-through options to save memory while reading `srccat`
  parquet files during L3a production mode. Refactored L3 unit tests for easier testing of different config
  values. [[#SP-1854](https://jira.ipac.caltech.edu/browse/SP-1854)]


### Other

- `click` package pin to less than breaking release 8.3.0. [[#SP-1834](https://jira.ipac.caltech.edu/browse/SP-1834)]

## SPHEREx Pipelines 6.4.1 (2025-09-02)


### Changes

- Nonfunctional pixel map: updated input calibration product handling for nonfunc pixels, including, formatter, DetectorPixelData type, dataset type, fake products and updated scripts to ingest into the butler. [[#SP-1801](https://jira.ipac.caltech.edu/browse/SP-1801)]
- Query Reference Catalog: Include reference catalog version in srccat tiles. [[#SP-1807](https://jira.ipac.caltech.edu/browse/SP-1807)]


### Bug Fixes

- L2a Date Header Fix: Creation time of FITS file based on creation of L2a file rather than
  inheriting date from L1. [[#SP-1806](https://jira.ipac.caltech.edu/browse/SP-1806)]
- Perform Forced Photometry Bugfix: Fixed central source filtering to correctly
  index the inverse variance. [[#SP-1811](https://jira.ipac.caltech.edu/browse/SP-1811)]


### Enhancements

- L3a Memory Optimization: Improved memory consumption by reducing array copying, enhancing window boundary
  object representation, and adding a manual garbage collection option. [[#SP-1764](https://jira.ipac.caltech.edu/browse/SP-1764)]

## SPHEREx Pipelines 6.4 (2025-08-26)


### New Features

- Various CLI and system exit code/error handling changes for photometry tool [[#SP-1769](https://jira.ipac.caltech.edu/browse/SP-1769)]
- README for L2a Configuration Override sample: Include in README instructions
  on how to use the sample config override file for L2a. [[#SP-1798](https://jira.ipac.caltech.edu/browse/SP-1798)]
- Variance map handling for overflow case and mitigation of outliers. [[#SP-1800](https://jira.ipac.caltech.edu/browse/SP-1800)]
- Update variance handling of transients and check electronic noise: Implemented updated module algorithm v4.1 [[ssdc-ms-011](https://spherex-docs.ipac.caltech.edu/ssdc-ms-011/v/main/)] [[#SP-1588](https://jira.ipac.caltech.edu/browse/SP-1588)]
- Set default behavior to add position and velocity to image header, in the case where an input SPK observation snippet is missing, the error is logged and processing continued. [[#SP-1761](https://jira.ipac.caltech.edu/browse/SP-1761)]

### Changes

- Overflow threshold obtained from image header should be scaled by 4. [[#SP-1505](https://jira.ipac.caltech.edu/browse/SP-1505)]
- Added position and velocity reference frames in image header. [[#SP-1766](https://jira.ipac.caltech.edu/browse/SP-1766)]
- Reference frame and FITS standard info in image header. [[#SP-1763](https://jira.ipac.caltech.edu/browse/SP-1763)]
- SPHEREx Reference Catalog v0.73: Updated scripts to support SPHEREx Reference Catalog v0.73. [[#SP-1776](https://jira.ipac.caltech.edu/browse/SP-1776)]
- Persistence module parameters: update to detector-dependent parameters ip\_sat, flu\_sat, r\_flu. [[#SP-1794](https://jira.ipac.caltech.edu/browse/SP-1794)]
- Level 3a: Added support for using Instrument data to obtain center-pixel wavelength and bandwidth,
  in addition to the wavelength model from the Simulator. [[#SP-1796](https://jira.ipac.caltech.edu/browse/SP-1796)]
- Average PSF per-detector list for config parameters.  [[#SP-1592](https://jira.ipac.caltech.edu/browse/SP-1592)]


### Bug Fixes

- opencv-python dependency moved to opencv-python-headless [[#SP-1781](https://jira.ipac.caltech.edu/browse/SP-1781)]
- Patching Timing bugs in Persistence Module: Utilizing 'TELAPSE' in place of 'TINT' for total integration time, reset to measure
  starts to the end of previous exposure to middle of evaluating exposure, and decay time of previous persistence is end of
  previous previous exposure to most previous exposure. [[#SP-1782](https://jira.ipac.caltech.edu/browse/SP-1782)]
- Perform Forced Photometry Bugfix: Subtract global background in place to save
  memory. [[#SP-1783](https://jira.ipac.caltech.edu/browse/SP-1783)]
- Bugfix for source masking, fixed application of proper motion correction. [[#SP-1784](https://jira.ipac.caltech.edu/browse/SP-1784)]
- Fixed warnings related to missing HIERARCH cards when using FITS header keywords longer than 8 characters. [[#SP-1797](https://jira.ipac.caltech.edu/browse/SP-1797)]
- Observer in SPICE routine call using SPK.  [[#SP-1754](https://jira.ipac.caltech.edu/browse/SP-1754)]
- unit of SUR should be adu. [[#SP-1482](https://jira.ipac.caltech.edu/browse/SP-1482)]


### Enhancements

- Perform Forced Photometry Improvements: Enhanced method to retrieve the interpolated PSF. Optimized Lanczos
  resampling and downsampling of the high-resolution PSF used by Tractor. [[#SP-1729](https://jira.ipac.caltech.edu/browse/SP-1729)]
- Photometry performance improvements: Using numpy array indexing and JIT Numba caching to compute window mask [[#SP-1773](https://jira.ipac.caltech.edu/browse/SP-1773)]
- Ensuring SPK coverage. [[#SP-1756](https://jira.ipac.caltech.edu/browse/SP-1756)]

## SPHEREx Pipelines 6.3.1 (2025-08-05)


### New Features

- Evolved Stellar Positions: Modified Query Reference Catalog to update stellar
  positions based on GAIA proper motions and reference epoch. Added "ref_ra" and
  "ref_dec" as pass-through column options in Perform Forced Photometry to
  report reference catalog positions in output. Corrected Tractor utility functions
  to use the cutout origin instead of the cutout center when querying the PSF
  from shifted positions. [[#SP-1537](https://jira.ipac.caltech.edu/browse/SP-1537)]
- Tagged collection for Primary Calibrators: Create a tagged collection for primary calibrators. [[#SP-1733](https://jira.ipac.caltech.edu/browse/SP-1733)]
- New `SPECIAL` SPSFlag for Targeted observations of comet 3I/Atlas. [[#SP-1755](https://jira.ipac.caltech.edu/browse/SP-1755)]
- Photometry tool version 11: Updated documentation/CLI help commands, and software versioning. [[#SP-1758](https://jira.ipac.caltech.edu/browse/SP-1758)]
- Validated position and velocity extracted from SPK with JPL Horizons. [[#SP-1744](https://jira.ipac.caltech.edu/browse/SP-1744)]
- Reinstate ODIR: created a minimal ODIR product that only includes the reobservations based on N Transients. [[#SP-1732](https://jira.ipac.caltech.edu/browse/SP-1732)]
- ODIR connection to L1: Adding ODIR to L1 output and creating aggregate ODIR pipeline yaml. [[#SP-1742](https://jira.ipac.caltech.edu/browse/SP-1742)]


### Changes
- Rewrite Estimate PSF module from scratch using Adreas Faisst's notebook model for one image with one psf zone. Includes RLD number of iterations optimization based on photometry fit quality metric. [[#SP-1738](https://jira.ipac.caltech.edu/browse/SP-1738)]
- Level 3 Product Versions: Added `spherex-pipelines` version metadata to the
  Perform Forced Photometry and Create All Sky Catalog output products. [[#SP-1617](https://jira.ipac.caltech.edu/browse/SP-1617)]
- Spectrophotometry Tool: Updated the workflow to use a Quantum-Backed Butler for executing all processes
  in the three-task L3Pipeline. This change eliminates temporary Butler registry modifications, preventing
  conflicts when multiple processes run concurrently. [[#SP-1695](https://jira.ipac.caltech.edu/browse/SP-1695)]
- Spectrophotometry Tool: Input collections for science images and calibration products are now configured separately,
  each via its own environment variable. Instead of querying calibrations dynamically for matching exposures,
  all relevant calibration datasets are pre-exported to the temporary Butler used during a single run. [[#SP-1713](https://jira.ipac.caltech.edu/browse/SP-1713)]
- Spectrophotometry Tool: Added missing RA/DEC column `FIELDref` elements
  required for cutout service. [[#SP-1760](https://jira.ipac.caltech.edu/browse/SP-1760)]
- Return WCS solution from Estimate Coarse Astrometry if FINAST > 1. [[#SP-1734](https://jira.ipac.caltech.edu/browse/SP-1734)]

### Bug Fixes

- Corrections to FITS headers in Level 2 (L2b) images. [[#SP-1646](https://jira.ipac.caltech.edu/browse/SP-1646)]
- Removed the obsolete pre-TVAC wavelength model that `get_wavelength_mesh` returned when `npix_x` or `npix_y` were specified. [[#SP-1712](https://jira.ipac.caltech.edu/browse/SP-1712)]
- Spectrophotometry Tool: Fixed dimensions for Combine Catalog Tiles task. [[#SP-1714](https://jira.ipac.caltech.edu/browse/SP-1714)]
- Estimate Fine Astrometry Bugfix: Modified inequalities in conditional
  statements to ensure "FINAST" header keyword gets written for all daxis values. [[#SP-1716](https://jira.ipac.caltech.edu/browse/SP-1716)]
- Renamed the WCS coordinate `Bandpass` to `Bandwidth` to reflect correct terminology.
  The previous name was inaccurate: `bandpass` refers to the full spectral response function
  (e.g., the filter transmission curve), whereas `bandwidth` refers to a scalar metric (e.g., FWHM)
  that characterizes its width. [[#SP-1720](https://jira.ipac.caltech.edu/browse/SP-1720)]
- Photometry tool inputs: Time is fixed to be optional. [[#SP-1726](https://jira.ipac.caltech.edu/browse/SP-1726)]
- Spectrophotometry Tool Bugfix: Fixed indexing error to enable proper insertion
  of external service descriptors in XML file. [[#SP-1737](https://jira.ipac.caltech.edu/browse/SP-1737)]
- Create All Sky Catalog: Removed "No primary calibrators..." warning to enable
  cleaner logs and better align with the expected primary calibrator distribution
  (most HEALPix tiles will not have any). [[#SP-1740](https://jira.ipac.caltech.edu/browse/SP-1740)]
- Perform Forced Photometry Bugfix: Corrected input positions when fetching the blended PSF to account for the
  x/y-axis shift in the SPHEREx detector pixel coordinates. [[#SP-1741](https://jira.ipac.caltech.edu/browse/SP-1741)]
- Optional Detector in Photometry CLI: Makes detector command line input optional when running photometry tool. [[#SP-1748](https://jira.ipac.caltech.edu/browse/SP-1748)]
- Create All Sky Catalog: Added string column length keyword to pyarrow metadata
  to prevent `astropy` warnings when reading parquet files. [[#SP-1757](https://jira.ipac.caltech.edu/browse/SP-1757)]


### Enhancements

- Added an LRU cache (`maxsize=1`) to `get_wavelength_mesh` to speed up repeated calls. [[#SP-1715](https://jira.ipac.caltech.edu/browse/SP-1715)]
- Replaced the `region_msk` dictionary with a lazily evaluated iterable to avoid allocating all masks at once. [[#SP-1718](https://jira.ipac.caltech.edu/browse/SP-1718)]

## SPHEREx Pipelines 6.2.2 (2025-06-24)


### New Features

- Operations Support: Added tag_releasable_data.py script to populate a TAGGED collection with data that pass DQA
  and are eligible for release. [[#SP-1692](https://jira.ipac.caltech.edu/browse/SP-1692)]

### Changes

- Exclude non-survey data from SSDC calibration pipelines. [[#SP-1210](https://jira.ipac.caltech.edu/browse/SP-1210)]
- Exclude in-sband observations from SSDC calibration pipelines. [[#SP-1236](https://jira.ipac.caltech.edu/browse/SP-1236)]
- Estimate Diffuse Background: Implemented estimation for spectral background
  component (`SSDC-MS-032 v3.0`). Added check on sigma-clipped median of unmasked
  pixels after subtraction. Streamlined unit tests for level 3 modules with
  shared inputs. [[#SP-1437](https://jira.ipac.caltech.edu/browse/SP-1437)]
- Config switch to do aberration correction in estimate fine astrometry. [[#SP-1703](https://jira.ipac.caltech.edu/browse/SP-1703)]

### Bug Fixes

- Fixed an issue where the software version of the input image was incorrectly saved in the VERSION keyword
  of the IMAGE extension. [[#SP-1688](https://jira.ipac.caltech.edu/browse/SP-1688)]
- Degenerate Polygons in Updated Exposure Regions: By filtering for good fine astrometry results
  in Level 2b images (FINAST == 0), it prevents forming degenerate polygons when updating
  exposure regions. [[#SP-1690](https://jira.ipac.caltech.edu/browse/SP-1690)]
- Fixed a failure in active area determination when some reference rows were missing.
  Active area boundaries are now derived from image readout parameters instead of pixel flags,
  since `REFERENCE` and `PHANTOM` flags are only set on non-missing pixels. [[#SP-1697](https://jira.ipac.caltech.edu/browse/SP-1697)]
- Fixed an issue that prevented the `detector` dimension table from being updated when using `butler register-instrument --update`.
  This command must be used explicitly to update basic dimension tables when needed. [[#SP-1708](https://jira.ipac.caltech.edu/browse/SP-1708)]


## SPHEREx Pipelines 6.2.1 (2025-06-11)


### New Features

- IRSA Tools Interface `SSDC-IF-008`:
  - Added dry-run support to the IRSA photometry tool. [[#SP-1586](https://jira.ipac.caltech.edu/browse/SP-1586)]
  - Unified CLI interface between the photometry tool and IRSA tool. [[#SP-1618](https://jira.ipac.caltech.edu/browse/SP-1618)]
  - Provided a UWS photometry tool template with `photometry-tool perform-photometry` interface for IRSA UWS. [[#SP-1424](https://jira.ipac.caltech.edu/browse/SP-1424)]
- Spectrophotometry Tool `SSDC-MS-051`:
  - Added functionality to run L3a pipeline and primary catalog generation. [[#SP-1120](https://jira.ipac.caltech.edu/browse/SP-1120)]
  - Added support for generating a source catalog from the user’s source list and registering it with the LVF-matching data IDs in a temporary Butler repository for use in the forced photometry pipeline. [[#SP-1519](https://jira.ipac.caltech.edu/browse/SP-1519)]
  - Added support for creating a temporary Butler repository containing all inputs required to run the forced photometry pipeline. [[#SP-1520](https://jira.ipac.caltech.edu/browse/SP-1520)]
  - Added `lvf_query` function to find LVF (Level2) images that contain sources from a user-provided source list. [[#SP-1532](https://jira.ipac.caltech.edu/browse/SP-1532)]
- Data Delivery:
  - Created IRSA Load ID as specified in `SSDC-IF-003, Section 3.2.2` to handle product delivery. [[#SP-1281](https://jira.ipac.caltech.edu/browse/SP-1281)]
  - ObsCore file generation: Added `obs_title` and units to relevant columns, and CAOM SPS plan string support. [[#SP-1383](https://jira.ipac.caltech.edu/browse/SP-1383)]
  - Data Integrity and Checksums: Implemented a script to parallelize MD5 checksum calculations before transferring datasets outside the NMDC, storing the results in the `file_datastore_records` table within the Butler-associated database for integrity verification. [[#SP-1349](https://jira.ipac.caltech.edu/browse/SP-1349)]
- L1 GDS Processing and Association [[#SP-1563](https://jira.ipac.caltech.edu/browse/SP-1563)]:
  - Enabled ingestion of `sgt` and `spk` files in L1.
  - Added support for optional selective file ingestion in GDS association.
  - Included start, middle, and end velocities and positions in SPK header ingest.
  - Fixed minor pipeline bugs related to solar system object computation and pinned `kete==1.0.6`.
- Data Model:
  - Added a sequence flag dimension to track incomplete GDS data fragments. [[#SP-1450](https://jira.ipac.caltech.edu/browse/SP-1450)]

### Changes

- FITS Headers and Metadata:
  - Updated names of spacecraft position and velocity vectors in the SPHEREx image header. [[#SP-1143](https://jira.ipac.caltech.edu/browse/SP-1143)]
  - Removed L1 pixel counting headers from level2 images. [[#SP-1655](https://jira.ipac.caltech.edu/browse/SP-1655)]
  - Added "L2 N_OUTLIER" keyword to flag header and improved L2 DQA unit tests to verify consistency of flag counts. [[#SP-1569](https://jira.ipac.caltech.edu/browse/SP-1569)]
- Spectral WCS:
  - Added support for MWIR flip in spectral WCS data. [[#SP-658](https://jira.ipac.caltech.edu/browse/SP-658)]
  - Added descriptive comments for spectral WCS FITS headers. [[#SP-1648](https://jira.ipac.caltech.edu/browse/SP-1648)]
- Spectral Calibration:
  - Implemented initial version of "spec cal lite" for a tabular bandpass model and central wavelength pixel map option in the mesh grid LVF module. [[#SP-1224](https://jira.ipac.caltech.edu/browse/SP-1224)]
- Fine Astrometry:
  - Revised parameterization to support Galactic-coordinate-dependent processing. Added background subtraction, error handling, and ran refinement twice to meet module requirements. `SSDC-MS-016, v2.5` [[#SP-1294](https://jira.ipac.caltech.edu/browse/SP-1294)]
  - Improved source position comparison between refinement iterations and selected the solution with lower relative difference. [[#SP-1397](https://jira.ipac.caltech.edu/browse/SP-1397)]
  - Added `FINAST` flag classification, indexing support for astrometry.net, and reinstated position angle from PC matrix. [[#SP-1541](https://jira.ipac.caltech.edu/browse/SP-1541)]
  - Enabled PV-to-SIP WCS conversion as a configurable option. [[#SP-1553](https://jira.ipac.caltech.edu/browse/SP-1553)]
  - Cleaned up temporary files generated by the fine astrometry module. [[#SP-1574](https://jira.ipac.caltech.edu/browse/SP-1574)]
- PSF Modeling:
  - Corrected PSF zone interpolation edge calculations to reduce residuals for blended sources. [[#SP-1357](https://jira.ipac.caltech.edu/browse/SP-1357)]
  - Updated average PSF parameter estimation in IOC. [[#SP-1539](https://jira.ipac.caltech.edu/browse/SP-1539)]
  - Parallelized PSF measurement and disabled jitter by default; added FWHM to plots. [[#SP-1578](https://jira.ipac.caltech.edu/browse/SP-1578)]
- Absolute Gain Matrix:
  - Normalized flat field and multiplied by previous flat before computing correction. [[#SP-1542](https://jira.ipac.caltech.edu/browse/SP-1542)]
  - Updated dark current correction to only apply when units are compatible with L2 images. [[#SP-1440](https://jira.ipac.caltech.edu/browse/SP-1440)]
  - Moved dark current subtraction to the L2 module where conversion to physical units occurs. [[#SP-1625](https://jira.ipac.caltech.edu/browse/SP-1625)]
  - Added provenance metadata for input calibrations in the absolute gain matrix product. [[#SP-1585](https://jira.ipac.caltech.edu/browse/SP-1585)]
  - Added optional method using slope fitting to estimate flat field and dark current. `SSDC-MS-026 v5.2` [[#SP-1599](https://jira.ipac.caltech.edu/browse/SP-1599)]
- Nonlinearity Corrections:
  - Corrected the omega pixel correction in the 1D gain function. [[#SP-1367](https://jira.ipac.caltech.edu/browse/SP-1367)]
  - Updated exposure time calculations and added handling for transient pixels. `SSDC-MS-010 v3.1.2` [[#SP-1370](https://jira.ipac.caltech.edu/browse/SP-1370)]
- Persistence Modelling:
  - Updated the persistence model by replacing the lookup table with a piecewise function and tuning parameters. [[#SP-1502](https://jira.ipac.caltech.edu/browse/SP-1502)]
  - Deprecated storage classes and formatters for per-pixel persistence parameters.
    Removed the corresponding per-pixel calibration products from the processing workflow. [[#SP-1445](https://jira.ipac.caltech.edu/browse/SP-1445)]
- Source Masking and Modeling:
  - Implemented radial functional fit for scaling source masking radius. [[#SP-1518](https://jira.ipac.caltech.edu/browse/SP-1518)]
  - Added new rendering methods to `SPHERExTractorPSF` to improve galaxy modeling. [[#SP-1466](https://jira.ipac.caltech.edu/browse/SP-1466)]
- Background Subtraction:
  - Updated the L3.2 module name from "Subtract Diffuse Background" to "Estimate Diffuse Background". [[#SP-1385](https://jira.ipac.caltech.edu/browse/SP-1385)]
  - Added local-background flags to photometry measurement tables and all-sky catalogs. [[#SP-1436](https://jira.ipac.caltech.edu/browse/SP-1436)]
- Forced Photometry:
  - Implemented all columns required for primary catalog construction and reduced memory usage by limiting columns passed to fitters. [[#SP-1346](https://jira.ipac.caltech.edu/browse/SP-1346)]
  - Filtered ignored pixels for Tractor, and converted the `flags` column into a bitmap for source regions. [[#SP-1377](https://jira.ipac.caltech.edu/browse/SP-1377)]
  - Changed the default value of the `combinedOutput` config parameter to `False`. [[#SP-1438](https://jira.ipac.caltech.edu/browse/SP-1438)]
  - Removed outdated utility functions and improved code readability. [[#SP-1443](https://jira.ipac.caltech.edu/browse/SP-1443)]
  - Adjusted galaxy position angle calculations based on cutout orientation and PC matrix determinant. [[#SP-1523](https://jira.ipac.caltech.edu/browse/SP-1523)]
  - Enabled class-based filtering to photometer only those windows containing sources matching the filter. [[#SP-1643](https://jira.ipac.caltech.edu/browse/SP-1643)]
- Photometry Catalogs:
  - Added substrings to external survey IDs to encode survey and data release. [[#SP-1375](https://jira.ipac.caltech.edu/browse/SP-1375)]
  - Added a boolean class filter to the "Create Primary Catalog" task. Removed the filtered-only output. [[#SP-1413](https://jira.ipac.caltech.edu/browse/SP-1413)]
  - Added a dedicated pipeline YAML for the primary calibrator catalog and changed default dataset type. [[#SP-1414](https://jira.ipac.caltech.edu/browse/SP-1414)]
  - Added a pipeline and task to combine all-sky catalogs across HEALPix tiles. [[#SP-1530](https://jira.ipac.caltech.edu/browse/SP-1530)]
- Spectrophotometry Tool Support:
  - Moved `external_source_id` column derivation to the Query Reference Catalog module. [[#SP-1555](https://jira.ipac.caltech.edu/browse/SP-1555)]
  - Added derived `ab_ratio` and `phi` columns to the Query Reference Catalog. [[#SP-1556](https://jira.ipac.caltech.edu/browse/SP-1556)]
  - Added config parameters for "pass-through" columns in photometry and all-sky catalog configurations. [[#SP-1535](https://jira.ipac.caltech.edu/browse/SP-1535)]
  - Added functions to create VOTables and embed IRSA resource links in the XML. [[#SP-1631](https://jira.ipac.caltech.edu/browse/SP-1631)]
  - Added input source catalog validator for L3a pipeline and improved unit tests for column compatibility. [[#SP-1526](https://jira.ipac.caltech.edu/browse/SP-1526)]
- Reference Catalogs:
  - Updated scripts for SPHEREx Reference Catalog v0.71. [[#SP-1388](https://jira.ipac.caltech.edu/browse/SP-1388)]
  - Updated scripts for SPHEREx Reference Catalog v0.72. [[#SP-1628](https://jira.ipac.caltech.edu/browse/SP-1628)]
  - Included missing magnitude in query reference catalog and omitted Pan-STARRS. [[#SP-1548](https://jira.ipac.caltech.edu/browse/SP-1548)]
- IDQR and DQA:
  - Added configuration for maximum transient pixels in IDQR modules. Marked images as bad and flagged observations for re-observation based on count. [[#SP-1570](https://jira.ipac.caltech.edu/browse/SP-1570)]
- CAOM and IRSA Delivery:
  - Added Load ID argument to CAOM generation for IRSA Butler export. [[#SP-1384](https://jira.ipac.caltech.edu/browse/SP-1384)]
  - Included target type in CAOM exports. [[#SP-1483](https://jira.ipac.caltech.edu/browse/SP-1483)]
  - Added `dark` and `solid_angle_pixel_map` products to IRSA deliveries. [[#SP-1554](https://jira.ipac.caltech.edu/browse/SP-1554)]
  - Implemented chained collection support to avoid re-delivering calibrations to IRSA. [[#SP-1456](https://jira.ipac.caltech.edu/browse/SP-1456)]
- Butler and Pipeline Utilities:
  - Created a script to generate tagged collections for primary calibrators with HEALPix6 ID lookup. [[#SP-1412](https://jira.ipac.caltech.edu/browse/SP-1412)]
  - Modified Butler queries to align with updated query methods. [[#SP-1419](https://jira.ipac.caltech.edu/browse/SP-1419)]
  - Added a `where` clause to enable export of Butler dataset subsets. [[#SP-1403](https://jira.ipac.caltech.edu/browse/SP-1403)]
  - Allowed filtering of observations by `sps_plan` in level-2a processing. [[#SP-1684](https://jira.ipac.caltech.edu/browse/SP-1684)]
  - Added `L2a_PERSISTENCE_STATE` tagged collection to track persistence model state. [[#SP-1463](https://jira.ipac.caltech.edu/browse/SP-1463)]
  - Added options to `list_dataset.py` to query by planning period and document code. [[#SP-1568](https://jira.ipac.caltech.edu/browse/SP-1568)]
- Pipeline Infrastructure:
  - Updated the nightly script to process all detectors with ingested calibrations: level-2a sequentially, level-2b in parallel. [[#SP-1378](https://jira.ipac.caltech.edu/browse/SP-1378)]
  - Revised calibration ingest script to use initial and instrument calibrations and removed test calibrations. [[#SP-1417](https://jira.ipac.caltech.edu/browse/SP-1417)]
  - Determined complete L1 history range for use in level-2a pipeline input selection. [[#SP-1462](https://jira.ipac.caltech.edu/browse/SP-1462)]

### Bug Fixes

- Fixed bugs found during L1 integration, including a return value issue in `CorrectPreamplifierTask` that blocked HISTORY comment generation. [[#SP-1076](https://jira.ipac.caltech.edu/browse/SP-1076)]
- Modified numpy array initialization in photometry utilities unit tests to avoid equality comparisons with NaNs. [[#SP-1369](https://jira.ipac.caltech.edu/browse/SP-1369)]
- Fixed cumulative source flag calculation in Perform Forced Photometry to only count bad pixel flags. [[#SP-1377](https://jira.ipac.caltech.edu/browse/SP-1377)]
- Corrected DQA behavior to create variance maps with valid pixels when comparing to the reference pixel distribution. [[#SP-1398](https://jira.ipac.caltech.edu/browse/SP-1398)]
- Updated connection type for solid angle pixel map in the photometry pipeline to enable proper calibration lookup. [[#SP-1405](https://jira.ipac.caltech.edu/browse/SP-1405)]
- Updated Postgres DB compatibility by converting galactic coordinate types from `np.float64` to native `float`. [[#SP-1406](https://jira.ipac.caltech.edu/browse/SP-1406)]
- Updated `_get_hires_psf` utility in Perform Forced Photometry to support L3a image shapes and PSF grids. [[#SP-1408](https://jira.ipac.caltech.edu/browse/SP-1408)]
- Added workaround to restore weak reference to uncertainty’s parent data after `SPHERExImage` is garbage collected. [[#SP-1425](https://jira.ipac.caltech.edu/browse/SP-1425)]
- Ordered calibration provenance HISTORY comments by dataset type to aid image comparison. [[#SP-1426](https://jira.ipac.caltech.edu/browse/SP-1426)]
- Fixed NaN propagation issue in per-exposure PSF measurement due to incorrect PSF flag selection. [[#SP-1431](https://jira.ipac.caltech.edu/browse/SP-1431)]
- Fixed CAOM ingestion when calibration collections are empty. [[#SP-1442](https://jira.ipac.caltech.edu/browse/SP-1442)]
- Added missing `SOURCE` and `CROSSTALK` flag definitions in Level 2 image headers. [[#SP-1507](https://jira.ipac.caltech.edu/browse/SP-1507)]
- Fixed inclusion operator bug in persistence module, correcting discontinuities in current estimation. [[#SP-1517](https://jira.ipac.caltech.edu/browse/SP-1517)]
- Updated CAOM generation with gain factor calibration level, datetime format fixes, and `metaRelease` for observation table. [[#SP-1527](https://jira.ipac.caltech.edu/browse/SP-1527)]
- Added missing `spsplan` dimension for image-to-catalog matching. [[#SP-1528](https://jira.ipac.caltech.edu/browse/SP-1528)]
- Corrected source catalog connection types to `Input` for tasks requiring multiple matched images. [[#SP-1529](https://jira.ipac.caltech.edu/browse/SP-1529)]
- Fixed deep field qualification logic for observations. [[#SP-1543](https://jira.ipac.caltech.edu/browse/SP-1543)]
- Treated `band` argument in IRSA photometry tool as a list of tuples. [[#SP-1559](https://jira.ipac.caltech.edu/browse/SP-1559)]
- Added function to update YAML metadata in Create All Sky Catalog; added unit tests for Astropy-parquet interoperability. [[#SP-1575](https://jira.ipac.caltech.edu/browse/SP-1575)]
- Fixed flag selection logic in Update Persistence to match truth table. [[#SP-1576](https://jira.ipac.caltech.edu/browse/SP-1576)]
- Preserved FITS keywords from L1 that Astropy removes. [[#SP-1577](https://jira.ipac.caltech.edu/browse/SP-1577)]
- Deployed `setuptools-scm` to assign SSDC pipeline version based on Git tag per [PEP 440](https://peps.python.org/pep-0440/). [[#SP-1579](https://jira.ipac.caltech.edu/browse/SP-1579)]
- Updated IDQR re-observe logic to exclude SAA fields and include flagged all-sky observations. [[#SP-1590](https://jira.ipac.caltech.edu/browse/SP-1590)]
- Applied correct PA calculation using CD matrix from fine astrometry; updated L1 detector PA handling. [[#SP-1598](https://jira.ipac.caltech.edu/browse/SP-1598)]
- Filtered out L2 images failing fine astrometry from photometry pipeline input. [[#SP-1601](https://jira.ipac.caltech.edu/browse/SP-1601)]
- Updated flux and wavelength conversions in Perform Forced Photometry to align index ordering; reindexed pass-through columns. [[#SP-1606](https://jira.ipac.caltech.edu/browse/SP-1606)]
- Refactored aggregate-IDQR task and pipeline to operate directly on per-image IDQR datasets. [[#SP-1608](https://jira.ipac.caltech.edu/browse/SP-1608)]
- Fixed missing metadata for `nyquist_comp` field in All Sky Catalog pyarrow schema. [[#SP-1619](https://jira.ipac.caltech.edu/browse/SP-1619)]
- Fixed persistence module slew time estimation for cases with one or two prior exposures. [[#SP-1627](https://jira.ipac.caltech.edu/browse/SP-1627)]
- Removed dark calibration from L3a inputs due to its move to L2b. [[#SP-1645](https://jira.ipac.caltech.edu/browse/SP-1645)]
- Reset WCS to `RADESYS = ICRS` after fine astrometry and removed obsolete `EQUINOX`. [[#SP-1649](https://jira.ipac.caltech.edu/browse/SP-1649)]
- Reassigned zero-valued variance pixels to mean variance in Perform Forced Photometry. [[#SP-1651](https://jira.ipac.caltech.edu/browse/SP-1651)]
- Fixed mapping logic for filtered photometry measurement tables with fewer rows. [[#SP-1656](https://jira.ipac.caltech.edu/browse/SP-1656)]
- Fixed detector ID column type in primary All Sky Catalogs to use `uint8`. [[#SP-1658](https://jira.ipac.caltech.edu/browse/SP-1658)]
- Filtered non-finite pixels missed by flags in Perform Forced Photometry. [[#SP-1659](https://jira.ipac.caltech.edu/browse/SP-1659)]
- Renamed L3a pipeline definition file for spectrophotometry tool; moved tool-specific data to `photometry_tool` subdirectory. [[#SP-1663](https://jira.ipac.caltech.edu/browse/SP-1663)]
- Excluded NaNs and Infs from least squares fits. [[#SP-1668](https://jira.ipac.caltech.edu/browse/SP-1668)]
- Handled masking of bright sources near image edges. [[#SP-1671](https://jira.ipac.caltech.edu/browse/SP-1671)]
- Corrected `CRVAL2` in 0-based active pixel WCS of Level 2 images. [[#SP-1679](https://jira.ipac.caltech.edu/browse/SP-1679)]
- Increased Butler query limit in exposure region update task from 20,000 to 30,000. [[#SP-1689](https://jira.ipac.caltech.edu/browse/SP-1689)]
- Photometry tool: Fail only if all tasks fail. [[#SP-1662](https://jira.ipac.caltech.edu/browse/SP-1662)]


### Enhancements

- PSF Interpolation Speedup: Optimized the blended PSF method (and its helper
  functions) used during Tractor fitting to consume less CPU time and memory. [[#SP-1427](https://jira.ipac.caltech.edu/browse/SP-1427)]

### Other

- Enabled a workaround for `astropy/astropy#17810` to prevent parsing errors with `SPHERExImage` in Astropy 7.0.1. [[#SP-1250](https://jira.ipac.caltech.edu/browse/SP-1250)]
- Cached Simulator data in the GitHub workflow to accelerate unit tests and reduce redundant downloads. [[#SP-1376](https://jira.ipac.caltech.edu/browse/SP-1376)]
- Added a Slack bot that uploads files to a specified channel with token validation. [[#SP-1458](https://jira.ipac.caltech.edu/browse/SP-1458)]
- Implemented retry logic for Slack notifier to avoid webhook overloads and consumer disconnection.
- Updated Slack notification formatting to reduce line wrapping. [[#SP-1454](https://jira.ipac.caltech.edu/browse/SP-1454)]
- Added a shebang line to the L1 processing script for dotenv compatibility. [[#SP-1550](https://jira.ipac.caltech.edu/browse/SP-1550)]


## SPHEREx Pipelines 6.1.0 (2025-02-14)


### New Features

- Initial implementation of photometry flag and fit quality metric: Added columns
  for photometry flag and fit quality metric to output photometry table. Updated
  unit test data to include sample flags. Cleaned up unit tests to remove
  redundant code. [[#SP-1171](https://jira.ipac.caltech.edu/browse/SP-1171)]
- Data Model:
  - Added spatial dimension to exposure table;
    added `butler update-exposure-region` command that populates spatial information for an exposure. [[#SP-1048](https://jira.ipac.caltech.edu/browse/SP-1048)]
  - Added `Non-Survey Data` and `Diagnostic mode` flags to observation table. [[#SP-1209](https://jira.ipac.caltech.edu/browse/SP-1209)]
- IRSA Interface: 
  - Developed a proof-of-concept notebook demonstrating SPHEREx content creation for the Common Archive Observation Model (CAOM) interface in IRSA. [[#SP-1315](https://jira.ipac.caltech.edu/browse/SP-1315)]
  - CAOM and ObsCore table generation: modified scripts to meet IRSA requirements for CAOM, newer approach for ObsCore table independent of `lsst-dax-obscore`. [[#SP-1329](https://jira.ipac.caltech.edu/browse/SP-1329)]


### Changes
- GDS image association:
  - Added GDS image association to pipeline. [[#SP-1318](https://jira.ipac.caltech.edu/browse/SP-1318)]
  - Added SSDC sband telemetry aggregate to GDS image association; the data is included as a binary table extension to the image. `SSDC-MS-002 v2.1.1` https://spherex-docs.ipac.caltech.edu/ssdc-ms-002/ [[#SP-1270](https://jira.ipac.caltech.edu/browse/SP-1270)]
- Initial SIP solution: preserve original PC matrix and CDELT after initial sip solution. [[#SP-1289](https://jira.ipac.caltech.edu/browse/SP-1289)]
- Removed duplicate output dataset type SIP: the SIP is included in the output image's WCS. Added functionality to handle overriding the output WCS with a given WCS. [[#SP-1195](https://jira.ipac.caltech.edu/browse/SP-1195)]
- Detect Outlier Pixels: Median filter implementation with DQA check. `SSDC-MS-020, v0.4` [[#SP-818](https://jira.ipac.caltech.edu/browse/SP-818)]
- Source masking: Implemented updated Source masking version, Crosstalk masking, DBS ghost and Frame edge ghost masking.
  `SSDC-MS-019, v3.0` https://spherex-docs.ipac.caltech.edu/ssdc-ms-019/ [[#SP-981](https://jira.ipac.caltech.edu/browse/SP-981)]
- Level 2a Pipeline: Added an option to save the persistence model dataset only at the end of a l2a run. [[#SP-1255](https://jira.ipac.caltech.edu/browse/SP-1255)]
- Including RA and Dec for source flagged as Astrometric or PSF with no Gaia or PanStarrs data [[#SP-1263](https://jira.ipac.caltech.edu/browse/SP-1263)]
- Position angle: position angle calculation from PC matrix if not available in image header. [[#SP-1265](https://jira.ipac.caltech.edu/browse/SP-1265)]
- Update Photometric Calibration: Parallel processing by obseration for calculation of flux and color correction; needed to improve runtime. [[#SP-1311](https://jira.ipac.caltech.edu/browse/SP-1311)]
- Perform Forced Photometry Pixel Scale: Modified flux conversion function to use pixel scale calculated
  from the solid angle pixel map. Added column for extended source flag. Updated unit tests to include
  sample solid angle pixel map. [[#SP-1321](https://jira.ipac.caltech.edu/browse/SP-1321)]
- Level 1 and Level 2 Data Quality Assessment: 
  - Added flag statistics to image header. [[#SP-1340](https://jira.ipac.caltech.edu/browse/SP-1340)]
  - Added a DQA utility function to throw an error if any nan pixels are found that are not flagged as `MISSING_DATA`. This function is called in both L1 and L2 dqa modules. [[#SP-1254](https://jira.ipac.caltech.edu/browse/SP-1254)]
  - Added data quality checks on level 2 module methods and outputs. [[#SP-1268](https://jira.ipac.caltech.edu/browse/SP-1268)]
  - DQA on Estimate Variance module with reference pixels; mean of distribution around zero and sigma comparable to square root of variance. [[#SP-646](https://jira.ipac.caltech.edu/browse/SP-646)]
- Image Data Quality Report (IDQR):
  - Aggregate IDQR: Now creating single IDQR for all detectors. [[#SP-1266](https://jira.ipac.caltech.edu/browse/SP-1266)]
  - Revised the IDQR to include ecliptic coordinates, added `sps_target` and `sps_flags` fields. [[#SP-768](https://jira.ipac.caltech.edu/browse/SP-768)]
- Photometry Interpolated PSF: Updated Tractor fitting to use a blended PSF
  that smoothly varies across the detector. [[#SP-1341](https://jira.ipac.caltech.edu/browse/SP-1341)]
- PSF zone mapping: improved unit test to verify that the asymmetric PSF zones are mapped to the correct location. Zone numbering is from left to right, bottom to top. [[#SP-1252](https://jira.ipac.caltech.edu/browse/SP-1252)]
- Photometry Unit Tests: Use Tractor installation from IPAC-SW fork to perform
  unit tests. [[#SP-1303](https://jira.ipac.caltech.edu/browse/SP-1303)]
- Perform Forced Photometry: Added `spherex_class` column to output parquet table, which will enable Create All
  Sky Catalog to filter out primary calibrators. [[#SP-1362](https://jira.ipac.caltech.edu/browse/SP-1362)]
- All-sky Primary Catalog: 
  - Added columns `pmra`, `pmdec`, `parallax`, `nyquist_comp`, and function to calculate the Nyquist completeness. [SSDC-MS-035](https://spherex-docs.ipac.caltech.edu/ssdc-ms-035/) [[#SP-1271](https://jira.ipac.caltech.edu/browse/SP-1271)]
  - Select primary calibrators, calspec sources, and secondary calibrators from photometry measurements. Output an additional table of calibrator measurements in addition to primary all-sky catalog. See `SSDC-MS-035, v3.0`. [[#SP-1325](https://jira.ipac.caltech.edu/browse/SP-1325)]
- Solid angle pixel map: added output dataset type for solid angle pixel map calibration product. [[#SP-1322](https://jira.ipac.caltech.edu/browse/SP-1322)]
- Now ingesting solid angle pixel map as calibration product. [[#SP-1350](https://jira.ipac.caltech.edu/browse/SP-1350)]
- Calibration pipeline products: add history to SPHEREx pipeline calibration derived products. [[#SP-1323](https://jira.ipac.caltech.edu/browse/SP-1323)]
- Data Model:
  - File Templates: Updated the directory hierarchy for SSDC-produced files to facilitate volume management and data deliveries organized by survey planning period. [[#SP-1299](https://jira.ipac.caltech.edu/browse/SP-1299)]
  - Replaced instances of combined `observation` and `detector` dimensions with the `exposure` dimension.
    This change simplifies the creation of ObsCore tables using the `lsst-dax-obscore` package. [[#SP-1304](https://jira.ipac.caltech.edu/browse/SP-1304)]
  - Added `sps_target` and `in_sband` fields to the observation table. These fields represent the SPS Target ID and the IN_SBAND indicator from the predicted observation record, respectively. [[#SP-1310](https://jira.ipac.caltech.edu/browse/SP-1310)]
  - Renamed `CHKSUM_ERR` flag to `SUR_ERROR`. [[#SP-1338](https://jira.ipac.caltech.edu/browse/SP-1338)]
  - Renamed `SURERROR` flag to `SUR_MISMATCH`. [[#SP-1358](https://jira.ipac.caltech.edu/browse/SP-1358)]
  - Changed default Butler dataset type for Average PSF to `average_psf`. [[#SP-1348](https://jira.ipac.caltech.edu/browse/SP-1348)]
- Butler Configuration: separated the registry from the butler.yaml file by creating a registry.yaml file. This approach will enable you to retain default configurations while allowing the registry.yaml file to be overwritten for different registries. [[#SP-1278](https://jira.ipac.caltech.edu/browse/SP-1278)]
- SPHEREx Reference Catalog v0.7: Updated scripts to support SPHEREx Reference Catalog v0.7. [[#SP-1277](https://jira.ipac.caltech.edu/browse/SP-1277)]
- Data Export to IRSA: Modified the transferring process of datasets to remote hosts. [[#SP-1281](https://jira.ipac.caltech.edu/browse/SP-1281)]


### Bug Fixes

- Photometry Flagging Bug Fix: Fixed photometry flag to only consider pixels used
  for calculating a source's fit quality (instead of all pixels in flag cutout). [[#SP-1171](https://jira.ipac.caltech.edu/browse/SP-1171)]
- Perform Forced Photometry Bug Fix: Added further check using window cutout WCS to
  filter out sources with minor inconsistencies in predicted/actual window position. [[#SP-1328](https://jira.ipac.caltech.edu/browse/SP-1328)]
- Photometry Central Wavelength Calculation: Modified code to use WCS-converted
  array indices instead of rounded pixel positions when calculating the central
  wavelength. Fixed bug related to out of bounds indices for positions along image
  boundaries. [[#SP-1235](https://jira.ipac.caltech.edu/browse/SP-1235)]
- Photometry Unit Test Bug Fix: Fixed bug in unit tests due to numpy array
  comparison with NaN values from empty array initialization. [[#SP-1251](https://jira.ipac.caltech.edu/browse/SP-1251)]
- SPHEREx Simulator dependencies must be optional. [[#SP-1264](https://jira.ipac.caltech.edu/browse/SP-1264)]
- L1 DQA: Check percentage of NaN pixels is not greater than configurable threshold, default is 80%. [[#SP-1186](https://jira.ipac.caltech.edu/browse/SP-1186)]
- Updating Exposure Regions: Fixed an issue where updating exposure records failed due to missing
  Level 2 products. [[#SP-1296](https://jira.ipac.caltech.edu/browse/SP-1296)]
- Error on Simulated Image Ingest: Due to slight variations in simulated observation start/end time in images from different detectors for the same simulated observation, an update is passed to ensure it is found within the schema. [[#SP-1324](https://jira.ipac.caltech.edu/browse/SP-1324)]
- Astropy Date Warning: Since a standard function within astropy called `_generate_wcs_and_update_header` removes certain header keys such as DATE-BEG and DATE-END. To address this, the keys were saved and added post-function call. [[#SP-1327](https://jira.ipac.caltech.edu/browse/SP-1327)]
- BrokenPipeError Fix: Due to piping out the observation list to head and the receiver closing early. Modified output to handle exception and prevent it from closing early. [[#SP-1337](https://jira.ipac.caltech.edu/browse/SP-1337)]
- Increased the slew group counter range to 9. The previous limit of 4 observations per slew group was insufficient for IOC. [[#SP-1363](https://jira.ipac.caltech.edu/browse/SP-1363)]


### Enhancements

- Perform Forced Photometry Optimization: Optimized measurement table construction to reduce time and memory
  consumption. Modified to use fixed-length numpy arrays instead of lists. Eliminated pipeline-level
  configuration parameter for keeping frame sources. [[#SP-1251](https://jira.ipac.caltech.edu/browse/SP-1251)]


## SPHEREx Pipelines 6.0.4 (2024-12-03)


### Bug Fixes

- Fine Astrometry:
  - Added handling for non-converging subprocesses by introducing timeouts and excluding images that did not pass Level 1 DQA from processing. [[#SP-1185](https://jira.ipac.caltech.edu/browse/SP-1185)]
  - Added a configurable threshold for the minimum number of astrometric calibrators. [[SP-1177](https://jira.ipac.caltech.edu/browse/SP-1177)]
  - Enhanced temporary file handling to ensure cleanup after module completion. [[#SP-1178](https://jira.ipac.caltech.edu/browse/SP-1178)]
- Exposure Averaged PSF: Added handling to ignore NaNs in the resulting summed stack. [[SP-1179](https://jira.ipac.caltech.edu/browse/SP-1179)]
- Per-Exposure PSF: The output PSF is now the measured jitter convolved with the average PSF. [[SP-1246](https://jira.ipac.caltech.edu/browse/SP-1246)]
- SUR Validation: Replace all occurrences of extension name SUR with RAWDATA for L1 full sample data. [[#SP-1242](https://jira.ipac.caltech.edu/browse/SP-1242)]
- Level 2a: Added the `OBSID` keyword if missing in the input image; this keyword is necessary for generating the Image Data Quality Report snippet. [[#SP-1175](https://jira.ipac.caltech.edu/browse/SP-1175)]
- Level 2b: Resolved a processing failure due to an absent upstream source catalog. Updated `PrerequisiteInput` connection type for dataset types that might be unavailable to `Input`. [[#SP-1181](https://jira.ipac.caltech.edu/browse/SP-1181)]
- L1 DQA bugfix: Updated flag depth crtierion, removed COSMIC RAYS flag in lieu of TRANSIENT, removed simulator based NCR header and fixed flag count calculation. [[#SP-1193](https://jira.ipac.caltech.edu/browse/SP-1193)]
- L1 DQA: Task now returns pipeBase.Struct as the only value. [[#SP-1240](https://jira.ipac.caltech.edu/browse/SP-1240)]
- L2 DQA: Fixed an issue where a missing `OMEGA_MEDIAN` keyword in the header caused errors. [[#SP-1180](https://jira.ipac.caltech.edu/browse/SP-1180)]
- Diffuse Background Subtraction: Corrected the calculation for the total number of subregions/boxes. [[#SP-1184](https://jira.ipac.caltech.edu/browse/SP-1184)]
- Forced Photometry: Fixed X and Y image positions input to each window by finding indices that satisfy both dimension boundaries. [[#SP-1248](https://jira.ipac.caltech.edu/browse/SP-1248)]
- Astropy SIP Warning: Updated the process for constructing `SpectralWCS` from a FITS header. The new approach directly extracts spectral WCS-related keys from the original header without initializing a WCS object. This resolves the SIP warning generated by Astropy during `SpectralWCS` initialization. [[#SP-1205](https://jira.ipac.caltech.edu/browse/SP-1205)]


### Enhancements

- Added a function to transform between GWCS and WCS formalism. [[(SP-1126)](https://jira.ipac.caltech.edu/browse/SP-1126)]
- Perform Forced Photometry Memory Utilization [[#SP-1199](https://jira.ipac.caltech.edu/browse/SP-1199)]
  - Optimized cutout initialization for image data, variance layer, and local background estimation to prevent memory allocation errors.
  - Modified the `run_tractor_window` call to only pass the required image attributes instead of the whole image object.
  - Updated variable and parameter names for array numbers to `detector`.
  - Added local background subtraction cases to unit tests.
  - Improved conditional statements and header information related to fit objects.
- Photometry Enhancements [[#SP-1235](https://jira.ipac.caltech.edu/browse/SP-1235)] 
  - Added the central wavelength, X/Y detector pixel position, and detector ID for each source to the output photometry measurement table.
  - Introduced a configuration parameter to enable or disable modeling galaxies with shapes.
  - Filtered out sources outside the image boundary from being passed to each window.
- Integrated Data Quality Report: Improved task management for large datasets. [[SP-1241](https://jira.ipac.caltech.edu/browse/SP-1241)]
- Data Analysis Support: Added functionality to calculate the Median PSF FWHM. Median PSF FWHM and astrometry RMS are now included in the image header. [[SP-1243](https://jira.ipac.caltech.edu/browse/SP-1243)]


## SPHEREx Pipelines 6.0.0 (2024-10-30)


### New Features

- Module implementation:
  - L1.4 `Associate Images with GDS Data`: Initial implementation of S-Band data association to SPHEREx images; data is aggregated and added as a binary table extension in the image. See `SSDC-MS-002, v2.1.1`. [[#SP-1134](https://jira.ipac.caltech.edu/browse/SP-1134)]
  - L1.11 `Level 1 Data Quality Assessment`: Initial implementation of `SSDC-MS-009, v1.3`. Added utility functions and a Pydantic model for DQA records. [[#SP-583](https://jira.ipac.caltech.edu/browse/SP-583)]
  - L2.15 `Level 2 Data Quality Assessment`: Initial implementation of `SSDC-MS-025, v1.2`. Added a Pydantic model for DQA records. [[#SP-1055](https://jira.ipac.caltech.edu/browse/SP-1055)]
  - L2.16 `Estimate Flat Field`: Added functionality to correct for solar elongation. See `SSDC-MS-026, v5.0`. [[#SP-1088](https://jira.ipac.caltech.edu/browse/SP-1088)]
  - L3.3 `Perform Forced Photometry`: Initial implementation of `SSDC-MS-033, v3.0`. [[#SP-963](https://jira.ipac.caltech.edu/browse/SP-963)]
  - L3.7 `Update Photometric Calibration`: Initial implementation of the algorithm using spectral calibrators and simulated SPHEREx spectral calibrator measurements. `SSDC-MS-037, v1.5`. [[#SP-1099](https://jira.ipac.caltech.edu/browse/SP-1099)] 
- GDS Reports:
  - Initial implementation of the Image Data Quality Report (IDQR), which includes pointing solutions and data quality indicators.  [[#SP-1127](https://jira.ipac.caltech.edu/browse/SP-1127)]
- Provenance:
  - References to the calibrations used in L1/L2 image processing are now recorded in the HISTORY comments of the image primary header. [[#SP-585](https://jira.ipac.caltech.edu/browse/SP-585)]
  - Modules applied to L1/L2 images are now recorded in the HISTORY comments of the image primary header. [[#SP-586](https://jira.ipac.caltech.edu/browse/SP-586)]

  
### Changes

- L1.5 `Estimate Coarse Astrometry`:
  - Added the roll angle from the `ckernel` observation-associated table and position and velocity vectors from the `spk` observation-associated table. [[#SP-1125](https://jira.ipac.caltech.edu/browse/SP-1125)]
  - Ingested ray-tracing tabular data as an input calibration. This data is used to derive the initial distortion model and output the initial SIP solution. [[#SP-1133](https://jira.ipac.caltech.edu/browse/SP-1133)]
- L2.7 `Query Reference Catalog`: Updated the module to limit selected sources to those that match flags, magnitude, and color criteria. [[#SP-1150](https://jira.ipac.caltech.edu/browse/SP-1150)]
- L2.16 `Estimate Absolute Gain Matrix`: Added output for the absolute gain matrix and consolidated the `Apply Diffuse Flat Field Correction` with the `Apply Photometric Calibration` module. `SSDC-MS-014, v1.4`. [[#SP-1147](https://jira.ipac.caltech.edu/browse/SP-1147)]
- L3.2 `Estimate Diffuse Background`:
  - Updated task, config, and unit tests based on `SSDC-MS-32, v2.1`. [[#SP-1006](https://jira.ipac.caltech.edu/browse/SP-1006)]
  - Switched to the default background estimator instead of using the median value per box. Replaced the absolute galactic latitude cut with a check on the number of boxes below the masking threshold. Increased exclusion percentile and number of sigma clipping iterations. Added a fake source mask to the unit test for the galactic center case. [[#SP-1142](https://jira.ipac.caltech.edu/browse/SP-1142)]
- L3.3 `Perform Forced Photometry`:
  - Galaxy Photometry Implementation: Updated Query Reference Catalog to include DESI Legacy Survey morphology columns. Modified Tractor utility function to model point sources or Sersic galaxies based on Sersic radius. Added test data with an extended source for unit testing. Adjusted local background configuration parameters to prevent `Background2D` ValueErrors. Added an optional output connection to save Tractor fit results for local development/notebooks. [[#SP-1114](https://jira.ipac.caltech.edu/browse/SP-1114)]
  - Forced Photometry Pipeline: Fixed units returned by the local background estimation function. Added a unit test assertion to check that local background and image cutouts have compatible units. Modified `L3a` and `PerformForcedPhotometry` tasks to run with the `pipetask` command. [[#SP-1128](https://jira.ipac.caltech.edu/browse/SP-1128)]
  - `PSFZoneData` class: Stored the precomputed header as an attribute. Modified the `getPsf` method to use the header attribute for quicker PSF zone retrieval during Tractor window fitting. [[#SP-1130](https://jira.ipac.caltech.edu/browse/SP-1130)]
  - Perform Forced Photometry Speedups: Removed logic related to parallel window fitting. Changed PSF initialization from per-window to a single-time setup. Saved actual PSF zone centers as attributes for quicker PSF zone retrieval. [[#SP-1132](https://jira.ipac.caltech.edu/browse/SP-1132)]
- FITS Image Header:
  - Added ecliptic latitude, longitude, PA, and equatorial RA and Dec to the SPHEREx image header. [[#SP-1105](https://jira.ipac.caltech.edu/browse/SP-1105)]
  - Replaced `EXPSTART` and `EXPSTOP` headers with the standard FITS headers `MJD-BEG` and `MJD-END`. [[#SP-1144](https://jira.ipac.caltech.edu/browse/SP-1144)]


### Bug Fixes

- Coarse Astrometry: Corrected the initial WCS solution to apply to the active pixel area, with adjustments to Level 1 format using CRPIX. [[#SP-1104](https://jira.ipac.caltech.edu/browse/SP-1104)]
- Estimate Fine Astrometry: 
  - SIP distortion in the input image for SCAMP should be retained. [[#SP-1126](https://jira.ipac.caltech.edu/browse/SP-1126)]
- Perform Forced Photometry: Adjusted local background estimation parameters to prevent `Background2D` from raising `ValueErrors`. Assigned high values to NaN elements in the variance map to avoid issues with inverse variance calculations in Tractor. [[#SP-1148](https://jira.ipac.caltech.edu/browse/SP-1148)]
- Photometry Array IDs: Replaced hardcoded `ARRAY` keyword references in the image header with the `getArrayID` method, allowing photometry code to process images from either simulator or SPHEREx instruments. [[#SP-1149](https://jira.ipac.caltech.edu/browse/SP-1149)]
- Simulated Level 1 Image Naming: Updated regex expressions to reflect the use of "detector_" instead of "array-" in one-year simulated Level 1 image names. [[#SP-1165](https://jira.ipac.caltech.edu/browse/SP-1165)]
- General Fixes and Improvements [[#SP-1129](https://jira.ipac.caltech.edu/browse/SP-1129)]
  - Modified FITS headers to ensure `fitsverify` passes without errors for Level 1 and Level 2 images.
  - Corrected zero-based WCS for active pixels in Level 2 images.
  - Fixed column types in the spectral WCS lookup table.


### Other

- Dependency management: Added tractor and astrometry.net dependencies to run the Perform Forced Photometry module. [[#SP-1092](https://jira.ipac.caltech.edu/browse/SP-1092)], [[#SP-1102](https://jira.ipac.caltech.edu/browse/SP-1102)]
- Conflict resolution: Modified solar elongation utilities to lazily initialize MPI communicator objects. This prevents conflicts with TACC's message-passing network, which uses Python 3.9.7. [[#SP-1137](https://jira.ipac.caltech.edu/browse/SP-1137)]
- SPHEREx Reference Catalog v0.6: Updated scripts to support SPHEREx Reference Catalog v0.6. [[#SP-1151](https://jira.ipac.caltech.edu/browse/SP-1151)]
- Support for grouping photometry jobs by size: Added the source count from the `Query Reference Catalog` output to the task metadata. [[#SP-1226](https://jira.ipac.caltech.edu/browse/SP-1226)]


## SPHEREx Pipelines 5.1.0 (2024-08-22)


### New Features

- GDS data: 
  - Added a storage classes for GDS and PEF data; GDS data is serialized as an `ecsv` file and PEF as `json`. [[#SP-998](https://jira.ipac.caltech.edu/browse/SP-998)]
  - Added GDS data inputs to Coarse Astrometry module and Level 1 Pipeline. [[#SP-1070](https://jira.ipac.caltech.edu/browse/SP-1070)]
- Initial WCS with SIP in Coarse Astrometry: Obtains initial WCS from telemetry data and known layout of focal plane. Obtains initial distortion solution from ray tracing model. [[#SP-1005](https://jira.ipac.caltech.edu/browse/SP-1005)]
- Aggregate GDS data: Added functions to get mean, first, and last, values from GDS associated data. [[#SP-1057](https://jira.ipac.caltech.edu/browse/SP-1057)]
- Image metadata: Added support for saving the content of image metadata blocks in Butler. [[#SP-912](https://jira.ipac.caltech.edu/browse/SP-912)]


### Changes

- Overflow threshold handling: Added a utility function, `get_overflow_threshold`, that reads a given image header. The overflow threshold, in units of ADU, is returned, with handling for instrument delivered or simulated images. [[#SP-969](https://jira.ipac.caltech.edu/browse/SP-969)]
- L2.8 `Estimate Fine Astrometry`: Added a function to calculate the residual pixel map and output the result. [[#SP-970](https://jira.ipac.caltech.edu/browse/SP-970)]
- L2.18 `Update Field Dictortion`, see `SSDC-MS-028 v1.3`: 
  - Added functions to initialize polynomial fit for SIP coefficients and use standard library methods. [[#SP-986](https://jira.ipac.caltech.edu/browse/SP-986)]
  - Solid angle pixel map: Added a function to calculate the solid angle pixel map, and output the result. [[#SP-971](https://jira.ipac.caltech.edu/browse/SP-971)]
  - Assessment for pixel solid angle map: added functions to calculate the median and standard deviation of the solid angle pixel map. [[#SP-1087](https://jira.ipac.caltech.edu/browse/SP-1087)]
- Lookup table for persistence: Added initial calibration constant parameters and lookup table input and functionality. [[#SP-982](https://jira.ipac.caltech.edu/browse/SP-982)]
- Electronic gains calibration: Added provenance information. [[#SP-987](https://jira.ipac.caltech.edu/browse/SP-987)]
- Level 3 photometry pipeline: Converted `Task` subclasses to `PipelineTask` subclasses for first 3 module stubs
  (L3.1 `Correct Dark Current`, L3.2 `Subtract Diffuse Background`, L3.3 `Perform Forced Photometry`.) [[#SP-1032](https://jira.ipac.caltech.edu/browse/SP-1032)]
- Optional photometry output connections: Added a configuration parameter which allows two output modes (single
  combined table or per-HEALPix tables) for the L3.3 `Perform Forced Photometry` measurements. [[#SP-1053](https://jira.ipac.caltech.edu/browse/SP-1053)]
- Level 1 DQA support: L1.7 `Validate SUR` now calculates median and sigma without outliers. [[#SP-946](https://jira.ipac.caltech.edu/browse/SP-946)]
- Level 2 DQA Support:
  - Extra outputs for nonlinearity validation [[#SP-944](https://jira.ipac.caltech.edu/browse/SP-944)]
  - Saved the non-linearity, persistence, and cross-talk
    thresholds in output image headers so Level 2 DQA can check flagging. [[#SP-1091](https://jira.ipac.caltech.edu/browse/SP-1091)]
- End-to-End test: Edits for End-to-End test; added `instrument` option to Butler client commands. [[#SP-1077](https://jira.ipac.caltech.edu/browse/SP-1077)]


### Bug Fixes

- Bugfix to correct nonlinearity flagging. [[#SP-950](https://jira.ipac.caltech.edu/browse/SP-950)]
- Query reference catalog should use WCS `world2pix` and `pix2world` with SIP. [[#SP-1080](https://jira.ipac.caltech.edu/browse/SP-1080)]
- Tests: Fixed unit test for parquet table Butler formatter. [[#SP-754](https://jira.ipac.caltech.edu/browse/SP-754)]
- Dependencies: Added `matplotlib` as an explicit dependency to `pyproject.toml`.
- Fixed `ModuleNotFoundError` when running build_and_test GitHub action. [[#SP-984](https://jira.ipac.caltech.edu/browse/SP-984)]


### Other

- Added scripts that illustrate a possible way to transfer data from one Butler repository to another,
  using `butler import` functionality. [[#SP-934](https://jira.ipac.caltech.edu/browse/SP-934)]
- Butler export instructions for IRSA. [[SP-935](https://jira.ipac.caltech.edu/browse/SP-935)]
- Documentation cleanup: Fixed in-code documentation errors related to type references and block quotes. [[#SP-989](https://jira.ipac.caltech.edu/browse/SP-989)]


## SPHEREx Pipelines 5.0.0 (2024-05-14)


### New Features

- L1.6 `Validate Sample-up-the-Ramp`: Initial implementation of the algorithm. See `SSDC-MS-008, v2.5.1`. [[#SP-871](https://jira.ipac.caltech.edu/browse/SP-871)]
- L2.17 `Estimate Flat-Field Correction`: Implemented updated flat field correction algorithm using stacking on spectrally flattened images. Removed planar tilt and image pair differencing technique. See `SSDC-MS-026`, `v4.0.0`. [[#SP-870](https://jira.ipac.caltech.edu/browse/SP-870)]
- Level 2a: Added butler command to detect prior observation with persistence state and run Level 2a on all observations after it. [[#SP-909](https://jira.ipac.caltech.edu/browse/SP-909)]
- Plan Tools: Created `plan-tools` script to calculate IVOA ObsCore components from observation plan files. [[#SP-902](https://jira.ipac.caltech.edu/browse/SP-902)]

### Changes

- Level 2 pipeline automation: Level 2b pipeline is now executed via the `pipetask` command for streamlined configuration, resource consumption tracking, and provenance management. [[#SP-882](https://jira.ipac.caltech.edu/browse/SP-882)]
- Added alternative WCS to `SPHERExImage` extensions: flags, variance and zodi. [[#SP-824](https://jira.ipac.caltech.edu/browse/SP-824)]
- L2.17 `Estimate Flat-Field Correction`: Applied source masking to data processed in the algorithm, and changed calculation of `ff_delta` from mean to sigma clipped (nominally 2) of the stack. See `SSDC-MS-026`, `v4.0.0`. [[#SP-947](https://jira.ipac.caltech.edu/browse/SP-947)]
- Reference point handling in `Estimate Fine Astrometry`: Added a function to compute linear WCS fit while keeping CRPIX constant. Validated that CRPIX1 and CRPIX2 remain unchanged. [[#SP-887](https://jira.ipac.caltech.edu/browse/SP-887)]
- Distortion handling in `Estimate Fine Astrometry`: SIP coefficients for distortion are now applied to input for pattern matching in SCAMP. [[#SP-943](https://jira.ipac.caltech.edu/browse/SP-943)]
- Data Model Improvements:
  - Revised Butler schema: It is now possible to query observations by Survey Plan ID. [[#SP-903](https://jira.ipac.caltech.edu/browse/SP-903)]
  - Observation ID: Fully implemented exposure and observation ID conventions described in the `Binary Form of Observation ID` document. [[#SP-827](https://jira.ipac.caltech.edu/browse/SP-827)]
  - Observation Metadata: Added Survey Planning Software flags to observation metadata. [[#SPS-908](https://jira.ipac.caltech.edu/browse/SPS-908)]
- Detector pixel-level data model (`DetectorPixelData`):
  - L2.3 `Update Persistence Mask`: Implemented the `DetectorPixelData` data model for persistence image  and time-evolving coefficients stored per-pixel. See `SSDC-MS-012`, `v1.1.2`. [[#SP-817](https://jira.ipac.caltech.edu/browse/SP-817)]
  - L2.3 `Update Persistence Mask`: Deployed the `DetectorPixelData` data model to pass per-pixel persistence parameters and state. See `SSDC-MS-012`, `v1.1.2`. [[#SP-814](https://jira.ipac.caltech.edu/browse/SP-814)]
  - L2.1 `Correct Non-Linearity`: Deployed the `DetectorPixelData` data model for non-linearity parameters. See `SSDC-MS-10`, `v3.0` [[#SP-820](https://jira.ipac.caltech.edu/browse/SP-820)]
  - L2.2 `Create Variance Map`: Implemented the `DetectorPixelData` data model for read noise calibration parameters. See `SSDC-MS-011`, `v3.2`. [[#SP-892](https://jira.ipac.caltech.edu/browse/SP-892)]

### Other

- Nightly Integration: Added `Dark Current Estimate`, `Update Field Distortion`, and `Estimate Average PSF` to the nightly integration script. [[#SP-800](https://jira.ipac.caltech.edu/browse/SP-800)]
- Calibration Ingest: Updated code to enable single calibration file ingest. [[#SP-940](https://jira.ipac.caltech.edu/browse/SP-940)]
- Testing Tools: Updated `setup-base-butler-repo` script to facilitate ingest of the simulated data. [[#SP-400](https://jira.ipac.caltech.edu/browse/SP-400)]


## SPHEREx Pipelines 4.1.0 (2023-12-21)


### New Features

- Reference Catalog with Magnitude: Added estimated magnitude calculation for each dichroic pair of detectors in `Query Reference Catalog` module. See `SSDC-MS-018 v2.0`. [[#SP-648](https://jira.ipac.caltech.edu/browse/SP-648)]
- New Persistance Model: Implemented new algorithm for `Update Persistence Mask` module that models persistence with additive parameter coefficients and lookup table. See `SSDC-MS-012 v1.1.2`. [[#SP-653](https://jira.ipac.caltech.edu/browse/SP-653)]
- Estimate Fine Astrometry: Implemented initial version of `Estimate Fine Astrometry` algorithm; module writes solved WCS to output image. See `SSDC-MS-016 v2.4`. [[#SP-776](https://jira.ipac.caltech.edu/browse/SP-776)]
- Measure Per-Exposure PSF: Implemented initial version of `Measure Per-Exposure Point Spread Function` (PSF) algorithm; added jitter difference kernel (JDK) function and convolved JDK with estimated PSF for output. See `SSDC-MS-023 v3.2`. [[#SP-668](https://jira.ipac.caltech.edu/browse/SP-668)]
- Average PSF with Multiple Exposures: Added functionality for multiple exposure processing in the `Estimate Exposure-Averaged Point Spread Function` module; tested in EstimatePSFPipeline. See `SSDC-MS-039 v1.2` [[#SP-802](https://jira.ipac.caltech.edu/browse/SP-802)]
- Estimate PSF Pipeline: Added a new pipeline for the estimate PSF module. [[#SP-786](https://jira.ipac.caltech.edu/browse/SP-786)]
- Butler Integration: Added a storage class for detector pixel-level data, storage class and formatter for non-linearity parameters. [[#SP-671](https://jira.ipac.caltech.edu/browse/SP-671)]

### Changes

- Example Reference Catalog: Updated pipeline code to work with Example Reference Catalog v0.4. Added scripts, located in `scripts/ref_catalog`, to find and copy catalog tiles that overlap with a set of images. Introduced options for catalog preprocessing and added functionality to include a per-image catalog in the iteration over images. [[#SP-650](https://jira.ipac.caltech.edu/browse/SP-650)]
- Source Mask Improvements: Added the lookup table for masking radius as a function of detector and magnitude in `Create Source Mask` module. See `SSDC-MS-019 v2.1`. [[#SP-719](https://jira.ipac.caltech.edu/browse/SP-719)]
- Gain Units: Aligned the use of gain factors such that the units are in electrons per adu. [[#SP-755](https://jira.ipac.caltech.edu/browse/SP-755)]
- Handling of Extra Reference Row Reads: Added funtionality to represent input calibration data, given in physical active and reference pixel area, with extra reference row reads. [[#SP-673](https://jira.ipac.caltech.edu/browse/SP-673)]

### Other

- Technical Debt: Fixed deprecation and resource warnings in unit tests. [[#SP-679](https://jira.ipac.caltech.edu/browse/SP-679)]
- Technical Debt: Removed pydantic model duplication and consolidated custom types. [[#SP-767](https://jira.ipac.caltech.edu/browse/SP-767)]

## SPHEREx Pipelines 4.0.0 (2023-08-27)


### New Features

- Spectral WCS: SPHEREx Level 1 and Level 2 images now contain approximated central wavelength and bandpass for each active pixel as defined in "Representations of spectral coordinates in FITS", [Greisen et al. 2006](https://www.aanda.org/articles/aa/pdf/2006/05/aa3818-05.pdf), section 6, "Coordinates by table lookup". [[#SP-461](https://jira.ipac.caltech.edu/browse/SP-461)]
- L2.2 `Create Variance Map`: Refined algorithm takes into account per-pixel read noise maps and pixel status flags. See `SSDC-MS-011 v3.1`. [[#SP-601](https://jira.ipac.caltech.edu/browse/SP-601)]
- L2.1 `Correct Non-Linearity`: Working implementation of the algorithm. See `SSDC-MS-010 v3.0`. [[#SP-644](https://jira.ipac.caltech.edu/browse/SP-644)]
- L2.22 `Estimate Exposure-Averaged Point Spread Function`: Initial implementation of the algorithm. See `SSDC-MS-039 v1.2`. [[#SP-645](https://jira.ipac.caltech.edu/browse/SP-645)]
- Image Data Quality Report (IDQR): Created a Pydantic model for IDQR. [[#SP-666](https://jira.ipac.caltech.edu/browse/SP-666)]


### Changes

- Format change to support checking variance on reference pixels: The change from L1 to L2 image format, which includes the removal of reference and phantom pixels and was previously the first step in the L2 pipeline, now occurs immediately after L2.3 `Update Persistence Mask`. [[#SP-647](https://jira.ipac.caltech.edu/browse/SP-647)]
- Readout order relative to LVF: Changed the extra reference pixel row reads to be on the lower wavelength side of the image readout to match the diagram from 2023/02/16 instrument delivery. [[#SP-580](https://jira.ipac.caltech.edu/browse/SP-580)]
- Revamped gain units and gain application procedure: Electronic gain per channel is now expressed in `ADU/e-` units. The conversion process from the instrument-reported units of `2E−16ADU/SURstep` to the engineering units of `e-/s` has been restructured.  Specifically, L1.3 `Assemble pixel data and decompress` now outputs an image in the units of `ADU/s`. Subsequently, L1.8 `Convert to Engineering Units` applies the per-channel gain factors. [[#SP-670](https://jira.ipac.caltech.edu/browse/SP-670)]
- Gain factors calibration product format: Transitioned `spx_data.GainFactors` storage class from `JSON` to `YAML` format. [[#SP-594](https://jira.ipac.caltech.edu/browse/SP-594)]
- Numeric precision controls: Added pipeline-level support to set the precision of data and per-pixel flags. [[#SP-630](https://jira.ipac.caltech.edu/browse/SP-630)]
- L2.4 `Apply Bad Pixel Mask`: Removed the assignment of `DICHROIC` and `NONFUNC` flags in the bad pixel mask module, and kept the assignment of `HOT` and `COLD` flags. [[#SP-634](https://jira.ipac.caltech.edu/browse/SP-634)]
- L1.10 `Create Level 1 Image Data Files`: Added placeholders for setting `NONFUNC`/`DICHROIC` flags from instrument ground products. [[#SP-643](https://jira.ipac.caltech.edu/browse/SP-643)]
- L2.1 `Correct Non-Linearity`: Updated the estimated time for an observation when an overlow or transient is encountered. Also updated the module algorithm with up-to-date equation. [[#SP-651](https://jira.ipac.caltech.edu/browse/SP-651)]
- Modified the calculation of PSF zone centers to allow for floating point values, in agreement with the FITS convention. [[#SP-680](https://jira.ipac.caltech.edu/browse/SP-680)]


### Bug Fixes

- Duplicate calibration products in Butler registry: Separated per-instrument from per-detector calibrations to avoid storing multiple copies of the same calibration in Butler registry. [[#SP-582](https://jira.ipac.caltech.edu/browse/SP-582)]
- L2.7 Query Reference Catalog: Fixed bugs and optimized algorithm. The Query Reference Catalog task returns reference catalog sources that overlap with the active area of the images separated by small slews. [[#SP-649](https://jira.ipac.caltech.edu/browse/SP-649)]

### Other

- Test Coverage Report: Integrated `pytest --cov` into GitHub Actions for measuring and tracking code coverage. [[#SP-588](https://jira.ipac.caltech.edu/browse/SP-588)]
- Tooling for release notes: Added infrastructure to produce useful, summarized change logs with `towncrier`. [[#SP-589](https://jira.ipac.caltech.edu/browse/SP-589)]
- Package and Tools Configuration: Migrated from `setup.cfg` to `pyproject.toml`, which provides a standardized way to define project metadata, build requirements, and tool configurations in Python projects. [[#SP-600](https://jira.ipac.caltech.edu/browse/SP-600)]
- Gain Factors: Added validation for Gain Factors ground calibration product. [[#SP-603](https://jira.ipac.caltech.edu/browse/SP-603)]
- Variance propagation: Added unit tests to confirm the variance is propagated as expected when changes to the data are made. [[#SP-639](https://jira.ipac.caltech.edu/browse/SP-639)]
- Pipeline Framework: Support for propagating gain and read noise calibration products through pipelines. [[#SP-641](https://jira.ipac.caltech.edu/browse/SP-641)]
