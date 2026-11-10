# Historical extension audit

This existing project's historical-coverage helper retrieves legacy ASAS data,
audits the unchanged DASCH/APPLAUSE quality cuts, and inventories APPLAUSE plate
footprints without treating coverage as photometry. Maintained source remains
in `tools/sky-plates`; raw responses and results remain
in outer `results/sky-plates-j0754-models-20261108/history`.

The existing astronomy runtime supplies Python, requests, astropy and matplotlib.
No service, model or database needs to be started. From the outer workspace:

```powershell
$skyPython = 'python'
$skySource = 'tools/sky-plates'
$skyHistory = 'results/sky-plates-j0754-models-20261108/history'
& $skyPython -B -m unittest discover -s $skySource -p test_historical_extension.py
& $skyPython -B "$skySource/historical_extension.py" --out $skyHistory --recon 'results/sky-plates-recon-20261107' --pilot 'results/sky-plates-j0754-20261107'
```

The offline audit expects `asas_raw.txt`, the original DASCH reconnaissance and
pilot data. It uses fresh APPLAUSE source/exposure snapshots and bounded footprint
inventory when present. Original URLs, exact ADQL, retrieval times and hashes are
in `query_receipts.jsonl`. To deliberately refresh one public response into a new
snapshot directory:

```powershell
& $skyPython -B "$skySource/fetch_historical_extension.py" --out $skyHistory --name asas_raw.txt --url 'https://www.astrouw.edu.pl/cgi-asas/asas_cgi_get_data?075446+1641.7,asas3'
```

The fetcher saves error responses and does not imply that HTTP success means a
measurement was returned. It permits read-only TAP queries with `--tap-query`.
Do not overwrite a retained snapshot without deliberately recording a new capture.

The saved audit recovers 422 ASAS rows, 390 A/B MAG_0 rows and 381 distinct HJD/frame
combinations spanning 2002–2009. It preserves four field zero points, all five
apertures, sentinels, duplicates and grades. Frame dispersions are not stellar
uncertainties. HJD has not been converted to BJD or assumed to have a more precise
time standard. The initial retrieval's annual shape is descriptive; the separate
comparison-star follow-up below supports a candidate smaller episode.

The APPLAUSE inventory covers 295 physical plates from a bounded 25-degree-center
search. It is not complete for all conceivable wider fields and provides no new
pixel photometry. Source measurements remain in their individual plate natural
systems. X/Y scans are not independent physical plates. Catalina and NSVS access
results are explicitly distinguished from recovered photometric data in
`history_audit.txt`. No flags are relaxed to fill historical gaps.

APPLAUSE provider documentation and acknowledgements are retained in the raw
snapshot, and the main component README supplies the full project attribution.
Legacy ASAS data and documentation are provided by the All Sky Automated Survey
at Warsaw University Observatory. NSVS metadata are from KASI's analysis of the
LANL survey; Catalina coverage is from its public CRTS data service.

## Same-exposure comparison follow-up

`historical_extension_controls_fetch.py` preserves requests for the original
sixteen Gaia controls. `historical_extension_controls.py` reads those snapshots
offline, applies the position/A/B/MAG_0 gates, centers each field/dataset on 2004,
and matches exact HJD/frame pairs. Repeated Gaia datasets and shared actual ASAS
catalogue apertures receive one vote. At least three control groups are needed.

All 338 target main-field exposures have 7–16 eligible control groups. Corrected
2007−2004 fading is +0.147 mag; 2009−2007 recovery is −0.114 mag. Apertures 1/2
and leave-one-control-out checks retain the shape. A root reconstruction restricted
to main-field controls gives +0.143125/−0.114 mag. Two individual controls show
similar-sized changes, so this is comparison-supported target variation without
a uniqueness, period, independent-survey replication or cause claim. Correlated
nearby apertures and target-specific systematics remain relevant.

```powershell
& $skyPython -B "$skySource/historical_extension_controls.py"
& $skyPython -B -m unittest discover -s $skySource -p test_historical_extension_controls.py
```

The fetch/analysis helpers use the fixed current result paths recorded in source.
The preserved control protocol predates downloads; aperture/jackknife checks are
post-protocol sensitivities. Three tests guard dataset baselines, exact frame
identity, repeated-star votes and shared native-aperture deduplication. Full raw
data, position audits, output and limitations are under `history/controls/`.
