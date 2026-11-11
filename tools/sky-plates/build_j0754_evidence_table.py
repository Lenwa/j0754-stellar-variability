"""Render the curated evidence ledger and inventory retained result families."""
import csv,hashlib,html,json,re
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DOCS=ROOT/'docs'
SOURCE=DOCS/'j0754_evidence_records.json'
FLARES=ROOT/'results/sky-plates-j0754-flare-comparison-20261111'

def readable(text):
    protected={}
    def hold(match):
        token=chr(0xe000+len(protected));protected[token]=match.group();return token
    # Identifiers, decades, thousands separators and scientific notation are
    # semantic tokens; spacing for compact prose must not alter them.
    text=re.sub(r'\b\d+(?:\.\d+)?e[+-]?\d+',hold,text)
    text=re.sub(r'\b(?:[123]eRASS|2SXPS|J\d{6}\.\d+[+-]\d{6}|\d{4}s|\d{1,3}(?:,\d{3})+)\b',hold,text)
    text=re.sub(r'(?<=[A-Za-z])(?=\d)|(?<=\d)(?=[A-Za-z])',' ',text)
    text=re.sub(r'(\d) e ([+-]?\d+)',r'\1e\2',text)
    text=re.sub(r'\b(SB|DR|W|E|S|K) (\d+)',r'\1\2',text)
    text=text.replace('cm 2','cm2')
    text=text.replace('2 SXPS','2SXPS')
    text=re.sub(r'(?<=[,;:])(?=[A-Za-z0-9])',' ',text)
    text=re.sub(r'(?<=[a-z0-9])\(', ' (',text)
    text=re.sub(r'\)(?=[A-Za-z])',') ',text)
    for token,value in protected.items():
        # Separate a protected numerical token from units without changing it.
        text=re.sub(re.escape(token)+r'(?=[A-Za-z])',token+' ',text)
        text=text.replace(token,value)
    return text

def main():
    data=json.loads(SOURCE.read_text());records=data['records'];ids=[r['id'] for r in records]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate evidence ID')
    allowed={'measured','candidate','inconclusive','not-recovered','model-compatible','model-rejected'}
    references={}
    for row in records:
        if row['status'] not in allowed:raise ValueError(row['status'])
        for name in ['finding','meaning','limit','independence']:
            if not row.get(name):raise ValueError(f'Missing {name}')
        for source in row['sources']:
            if source.startswith('https://'):continue
            p=(ROOT/source).resolve()
            if not p.is_relative_to(ROOT.resolve()) or not p.is_file():raise ValueError(f'Missing/invalid source {source}')
            references[source]=hashlib.sha256(p.read_bytes()).hexdigest()
    inventory=[]
    for directory in sorted((ROOT/'results').iterdir()):
        if not directory.is_dir() or not ('j0754' in directory.name or 'recon' in directory.name):continue
        files=[p for p in directory.rglob('*') if p.is_file()]
        inventory.append(dict(directory=directory.relative_to(ROOT).as_posix(),files=len(files),bytes=sum(p.stat().st_size for p in files),
            file_types=dict(sorted(Counter(p.suffix.lower() for p in files).items())),
            protocols=[p.relative_to(ROOT).as_posix() for p in files if 'protocol' in p.name or 'specification' in p.name],
            reports=[p.relative_to(ROOT).as_posix() for p in files if p.name.startswith('report')]))
    (DOCS/'J0754_dataset_inventory_2026-11-11.json').write_text(json.dumps(inventory,indent=2)+'\n',encoding='utf8',newline='\n')
    columns=['id','area','dataset','status','finding','meaning','limit','independence','sources']
    with (DOCS/'J0754_evidence_2026-11-11.csv').open('w',newline='',encoding='utf8') as fp:
        writer=csv.DictWriter(fp,fieldnames=columns);writer.writeheader()
        for r in records:writer.writerow({k:' | '.join(r[k]) if k=='sources' else readable(r[k]) for k in columns})
    escape=lambda s:html.escape(readable(str(s)))
    css='body{font:15px/1.5 system-ui,sans-serif;color:#172c3b;background:#f7f9fb;margin:0}main{max-width:1700px;margin:auto;padding:32px}h1{font-size:30px}h2{margin-top:36px}p{max-width:1100px}a{color:#12568e}table{width:100%;border-collapse:collapse;background:white}th,td{text-align:left;vertical-align:top;padding:12px;border:1px solid #d8e0e6}th{background:#16364d;color:white}tbody tr:nth-child(even){background:#f0f5f8}small{display:block;color:#42586a;margin-top:8px}figure{margin:24px 0}img{max-width:1100px;width:100%;height:auto}caption{text-align:left;padding:10px 0;font-weight:bold}.scroll{overflow:auto}.tag{display:inline-block;font-size:12px;padding:2px 6px;background:#e1ebf0;border-radius:3px;margin:5px 0}ul{max-width:1050px}@media print{body{font-size:10px;background:white}main{padding:0}th{color:black;background:#e5edf2}tr{break-inside:avoid}img{max-height:900px;width:auto;max-width:100%}}'
    parts=['<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',f'<title>J0754 evidence | {data["as_of"]}</title><style>{css}</style></head><body><main>',
        '<h1>J0754: the evidence together</h1>',f'<p>{escape(data["target"])}<br>Updated {data["as_of"]} · {len(records)} evidence entries · {len(inventory)} retained result families</p>',
        f'<p><strong>{escape(data["synthesis"])}</strong></p>',
        '<p>This is a source-linked research ledger, not a peer-reviewed classification. It incorporates the four-event TESS comparison before consolidating the wider evidence. Earlier reports remain historical snapshots.</p>',
        '<p><a href="J0754_evidence_2026-11-11.csv">CSV table</a> · <a href="j0754_evidence_records.json">Curated records</a> · <a href="J0754_dataset_inventory_2026-11-11.json">Dataset inventory</a> · <a href="../results/sky-plates-j0754-flare-comparison-20261111/report.txt">Four-event report</a></p>',
        '<h2>The four TESS excursions</h2><p>Common model measurements. ED is equivalent duration: integrated relative excess, expressed as seconds of quiescent stellar output. The model includes tails not fully observed. Face-on patch energies use the shared 9000 K flare, 5206 K photosphere and 54.64 solar-luminosity assumptions; isotropic equivalents are four times larger.</p>',
        '<div class="scroll"><table><thead><tr><th>Date</th><th>Peak excess</th><th>Model FWHM</th><th>Model ED / baseline–aperture range</th><th>Nominal face-on energy</th></tr></thead><tbody>']
    flare=json.loads((FLARES/'comparison.json').read_text())
    for event in flare['events']:
        f=event['photometry']['central'];lo,hi=event['energy']['model_ED_baseline_aperture_range_s']
        parts.append(f'<tr><td>{event["date_utc"][:10]}</td><td>{100*f["amplitude_fractional"]:.2f}%</td><td>{f["fwhm_h"]:.2f} h</td><td>{f["equivalent_duration_s"]:.0f} s / {lo:.0f}–{hi:.0f} s</td><td>{event["energy"]["face_on_patch_erg"]:.2e} erg</td></tr>')
    parts+=['</tbody></table></div><p>The broad November 16 event and the 2023 event have substantial baseline sensitivity. November 22 and December 23 show unusual apparent hourly position changes under the matched-control test. Image-gradient terms do not remove their brightness excess, but stationary-source astrometry remains unresolved. The data do not separate a close binary component.</p>',
        '<h2>Evidence table</h2><p>Status refers to the specific row. “Measured” does not certify its proposed physical explanation; “not recovered” does not imply universal absence. The lineage column prevents counting a new reduction or model fit as another instrument.</p>',
        '<div class="scroll"><table><thead><tr><th>Dataset / test</th><th>Finding</th><th>What it supports</th><th>Practical limit</th><th>Lineage / sources</th></tr></thead><tbody>']
    for row in records:
        links=[]
        for i,source in enumerate(row['sources'],1):
            url=source if source.startswith('https://') else '../'+source
            links.append(f'<a href="{html.escape(url,quote=True)}">Source {i}</a>')
        parts.append(f'<tr id="{row["id"]}"><td><strong>{row["id"]} · {escape(row["dataset"])}</strong><small>{escape(row["area"])}</small><span class="tag">{escape(row["status"])}</span></td><td>{escape(row["finding"])}</td><td>{escape(row["meaning"])}</td><td>{escape(row["limit"])}</td><td>{escape(row["independence"])}<small>{" · ".join(links)}</small></td></tr>')
    parts+=['</tbody></table></div><h2>What would change the picture</h2><ul>',
        '<li>Repeatable high-resolution velocities and line-shape/activity measurements across several 35.92-day cycles could test the companion and synchronization.</li>',
        '<li>Temperature-sensitive spectra plus simultaneous multicolor photometry at different brightness states could separate changing cool surface coverage from attenuation.</li>',
        '<li>An official-PRF/attitude and moving-source check could resolve the remaining TESS position diagnostics; new well-localized flares would help assign the emitting source.</li>',
        '</ul><h2>Plots and measurement record</h2><figure><img src="../results/sky-plates-j0754-flare-comparison-20261111/four_event_measurements.png" alt="Four-event light curves, common profile fits and baseline/aperture sensitivity"><figcaption>Common measurements; uncertainty from baseline choice is shown explicitly.</figcaption></figure>',
        '<figure><img src="../results/sky-plates-j0754-flare-comparison-20261111/four_event_localization.png" alt="Signed event-difference images and matched-injection location scales"><figcaption>Difference images and empirical-template position checks. These scales omit official-PRF/color systematics.</figcaption></figure>',
        '<p>Rebuild with <code>python -B tools/sky-plates/build_j0754_evidence_table.py</code>. Input source hashes and validation are retained beside this table. No telescope inputs were modified and no message was sent to an astronomer.</p></main></body></html>']
    (DOCS/'J0754_evidence_2026-11-11.html').write_text('\n'.join(parts)+'\n',encoding='utf8',newline='\n')
    summary=['J0754 — CONSOLIDATED EVIDENCE','11 November 2026','',data['synthesis'],'',
        'Read the formatted table: J0754_evidence_2026-11-11.html','Machine-readable table: J0754_evidence_2026-11-11.csv','',
        'The strongest measured results are repeat fading/recovery, independent optical replication, wavelength-dependent dimming, a recurring roughly36-day signal, published chromospheric emission and archival X-rays.',
        'The four TESS excursions remain positive under common measurements. Their conditional face-on energies are about2.0,2.5,7.1,2.1e37erg. November16/November28 have broad baseline envelopes. Apparent hourly position drift in November22/December23 remains an open attribution diagnostic; the tested pointing terms do not remove the brightness excess.',
        'Gaia gives the companion candidate. Activity detections are not independent orbit confirmations. Spots fit the color changes, but nonuniform obscuration survives. Radio images add limits; UVOT exposure edges prevent a flare inference; GALEX timing supplies only a weak variability hint.',
        'Next decisive measurements: calibrated high-resolution RV/activity spectra, simultaneous multicolor/temperature measurements, and the remaining TESS astrometry checks. Other-star screening is the next separately scheduled stage.','']
    for row in records:summary += [f"{row['id']} [{row['status']}] {readable(row['dataset'])}",readable(row['finding']), 'Supports: '+readable(row['meaning']),'Limit: '+readable(row['limit']),'Sources: '+' | '.join(row['sources']),'']
    summary=[readable(line) if not line.startswith('Sources:') else line for line in summary]
    (DOCS/'J0754_evidence_summary_2026-11-11.txt').write_text('\n'.join(summary),encoding='utf8',newline='\n')
    verification=dict(evidence_entries=len(records),result_families=len(inventory),local_source_files=len(references),all_local_references_exist=True,
        source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),local_source_sha256=references,
        qualification='Reference/hash verification does not independently validate every scientific interpretation. Curated claims retain their original measurement/model qualifications.')
    (DOCS/'J0754_evidence_verification_2026-11-11.json').write_text(json.dumps(verification,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({k:v for k,v in verification.items() if k!='local_source_sha256'},indent=2))

if __name__=='__main__':main()
