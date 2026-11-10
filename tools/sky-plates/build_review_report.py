"""Build the compact J0754 scientific-review PDF and offline supporting packet.

Run from repository root with the bundled document Python (reportlab and Pillow).
No network access or changes to telescope observations.
"""
import csv
import hashlib
import html
import json
from pathlib import Path
import shutil
import zipfile
from xml.sax.saxutils import escape
from reportlab.pdfgen import canvas
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from reportlab.lib.utils import ImageReader
from PIL import Image

ROOT = Path.cwd()
SOURCE = ROOT/'tools/sky-plates'
OUT = ROOT/'results/sky-plates-j0754-review-report-20261108'
W, H = 612, 792
MARGIN = 48
BODY = ParagraphStyle('body', fontName='Times-Roman', fontSize=11, leading=13.5,
                      spaceAfter=5)
SMALL = ParagraphStyle('small', fontName='Helvetica', fontSize=8.3, leading=11)
HEAD = ParagraphStyle('head', fontName='Helvetica-Bold', fontSize=11, leading=14)


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def paragraph(c, text, y, style=BODY, width=W-2*MARGIN):
    p=Paragraph(escape(text),style)
    _, height=p.wrap(width,H)
    p.drawOn(c,MARGIN,y-height)
    return y-height


def furniture(c, page, label):
    c.setFillColorRGB(0,0,0)
    c.setFont('Helvetica',8)
    c.drawString(MARGIN, H-28, 'J0754  |  Archival research report  |  8 November 2026')
    c.drawRightString(W-MARGIN, 26, f'{label}  |  {page}')


def image_box(c, path, x, top, width, height):
    with Image.open(path) as im: iw,ih=im.size
    scale=min(width/iw,height/ih)
    ww,hh=iw*scale,ih*scale
    c.drawImage(ImageReader(str(path)),x+(width-ww)/2,top-hh,
                width=ww,height=hh,mask='auto')
    return hh


FIGURES = [
 ('1','independent-20261108/asassn_replication.png','Modern fading and recovery',
  'Independent ASAS-SN camera measurements reproduce the ZTF fading and recovery. Camera and filter identities are retained; these streams are not freely offset to manufacture a common event.'),
 ('2','dust-20261108/optical_infrared.png','Optical and infrared behavior',
  'NEOWISE W1 and W2 fade less than optical g and partly recover through March 2024. Control corrections and near-contemporaneous matching support an optical-to-infrared color change. Static catalog colors are not event colors.'),
 ('3','models-20261108/historical_modern_pattern.png','Historical and modern episodes',
  'DASCH photographic measurements, legacy ASAS V and modern optical data retain their own passbands. The 1930s event and smaller 2000s trough are visible without waveform stretching or cross-survey event alignment.'),
 ('4','physical-tests-20261108/tess/sector46_pixel_localization.png','TESS event and image localization',
  'Sector 46 calibrated pixels show a compact excess near the target during the 2021 December 23 excursion. Offset is 0.177 pixel (3.7 arcsec). This moment-based localization has no formal PRF astrometric uncertainty; an unresolved alternative host is not excluded.'),
 ('5','optical-overlay-20261108/optical_flare_overlay.png','Optical flare shape comparison',
  'Peak- and FWHM-normalized profiles compare J0754 with the Davenport empirical optical flare template and an approximate digitization of a published solar white-light kernel (Hao et al. 2017). Shapes, not absolute energies or durations, are compared. No statistical template fit was performed.'),
 ('6a','physical-tests-20261108/fast/night_EW.png','FAST H alpha equivalent widths',
  'Night-averaged pseudo-equivalent widths use absorption-positive convention. November 2009 shows reduced absorption relative to October 2008. Plotted formal errors omit continuum/extraction systematics; a 1.675% differential continuum bias can account for the primary contrast.'),
 ('6b','fast-bluewing-20261108/registered_Halpha_bluewing.png','FAST registered line profiles',
  'Nearby metal-line registration and same-grating comparisons do not reveal a convincing detached blueshifted absorption component. The limited resolution and sparse epochs do not provide a stringent general limit on outflow.'),
 ('7','seasonal-rotation-20261108/seasonal_search.png','Rotation with changing seasonal phase',
  'Gray curves show individual seasons and blue curves their equally weighted mean fit improvement. Both surveys favor approximately 36 days; orange marks the Gaia candidate orbital period. Phase and amplitude are independent between seasons. Peak widths are not formal period confidence intervals.'),
 ('8','amplitude-brightness-20261108/amplitude_vs_brightness.png','Rotational strength and the long dip',
  'The broad fade does not track rotational amplitude monotonically. Nearby-peak amplitudes are subject to period-selection bias. The weak whole-season 2024/25 component partly masks stronger early modulation recovered in shorter windows; see the text and short-window result table.'),
 ('9','models-20261108/model_predictions.png','Competing color predictions',
  'Atmosphere-based model families fit g and W1 and predict additional colors. Curves and envelopes depend on the stated atmosphere/geometry grids and tolerance. W2 was not used to select the illustrative spot example. Overlap precludes a unique cause from these bands alone.'),
 ('10','spherex-20261108/photosphere_check.png','Recovery era infrared spectrum',
  'Local and official-tool SPHEREx reductions use the same observations. A photosphere fitted at shorter wavelengths accounts for most 2.45-5 micron flux. Residuals are sensitive to atmosphere and calibration choices; the plot does not establish a dust-free system.')
]


def build_pdf(content):
    path=OUT/'J0754_report.pdf'
    c=canvas.Canvas(str(path),pagesize=(W,H),pageCompression=1)
    c.setTitle(content['title']); c.setAuthor('Kyle; AI-assisted analysis by Ari')
    metrics=[]
    for index,page in enumerate(content['pages'],1):
        furniture(c,index,'Findings')
        y=H-50
        if index==1:
            title=ParagraphStyle('title',fontName='Helvetica-Bold',fontSize=18,leading=21)
            y=paragraph(c,content['title'],y,title)-6
            y=paragraph(c,content['subtitle'],y,SMALL)-6
            y=paragraph(c,content['credit'],y,SMALL)-12
        y=paragraph(c,page['title'],y,ParagraphStyle('page',parent=HEAD,fontSize=13,leading=16))-9
        for heading,text in page['sections']:
            if heading: y=paragraph(c,heading,y,HEAD)-3
            style=SMALL if heading=='Sources and reproducibility' else BODY
            y=paragraph(c,text,y,style)-6
        if y<44:
            raise ValueError(f'Text page {index} overflows: bottom {y:.1f}')
        metrics.append(dict(page=index,bottom_y=y));c.showPage()
    # Two ordinary plots per page; dense multi-panel figures receive full pages.
    groups=[FIGURES[:2], [FIGURES[2]], [FIGURES[3]], [FIGURES[4]],
            FIGURES[5:7], [FIGURES[7]], [FIGURES[8]], [FIGURES[9]], [FIGURES[10]]]
    for number,group in enumerate(groups,5):
        furniture(c,number,'Figure appendix')
        y=H-52
        for no,rel,title,caption in group:
            y=paragraph(c,f'Figure {no}  {title}',y,HEAD)-8
            cap=Paragraph(escape(caption),SMALL);_,ch=cap.wrap(W-2*MARGIN,H)
            image_height=(270 if len(group)==2 else 590)-ch
            used=image_box(c,ROOT/f'results/sky-plates-j0754-{rel}',MARGIN,y,W-2*MARGIN,image_height)
            y-=used+7;y=paragraph(c,caption,y,SMALL)-18
        if y<44: raise ValueError(f'Figure page {number} overflows: {y}')
        metrics.append(dict(page=number,bottom_y=y));c.showPage()
    c.save()
    (OUT/'layout_metrics.json').write_text(json.dumps(metrics,indent=2))
    return path


def package_support(content):
    # Only this investigation's files; never include private workspace notes,
    # session history, account/job receipts, credentials, or third-party papers.
    gallery=OUT/'figures';gallery.mkdir(exist_ok=True)
    skipped={'page1.png','Tang2011_STScI_poster.png','Tang2012_thesis_page_56.png','Tang2012_thesis_page_64.png'}
    figures=[];seen={}
    for folder in sorted((ROOT/'results').glob('sky-plates-j0754*')):
        if folder==OUT: continue
        for p in sorted(folder.rglob('*.png')):
            if p.name in skipped: continue
            sha=digest(p)
            if sha in seen:
                seen[sha]['duplicate_source_paths'].append(p.relative_to(ROOT).as_posix());continue
            rel=p.relative_to(ROOT/'results')
            target=gallery/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
            record=dict(path='figures/'+rel.as_posix(),source=p.relative_to(ROOT).as_posix(),
                        sha256=sha,duplicate_source_paths=[])
            seen[sha]=record;figures.append(record)
    blocks=['<!doctype html><meta charset="utf-8"><title>J0754 analysis figures</title>',
        '<style>body{font:16px system-ui;max-width:1100px;margin:35px auto;color:#171717}img{max-width:100%;height:auto}section{margin:35px 0}h2{font-size:17px;overflow-wrap:anywhere}p{line-height:1.5}</style>',
        '<h1>J0754 complete analysis figure collection</h1><p>Original resolution figures from this investigation. Read J0754_report.pdf first. Earlier diagnostic figures and failed extraction attempts are preserved for audit; their appearance here does not validate their measurements. The report and current JSON results take precedence over superseded plot labels. Published paper pages are linked in the report rather than reproduced.</p>']
    for i,r in enumerate(figures,1):
        name=html.escape(r['source'])
        status=''
        if any(s in name for s in ('precomputation','/probe/','applause','scan_check')):
            status='<p>Exploratory or incomplete diagnostic; not an accepted independent variability measurement.</p>'
        if 'fixed_period_diagnostic' in name:
            status='<p>Earlier stable-waveform test. Its failure is distinct from the later changing-season rotation detection.</p>'
        blocks.append(f'<section><h2>{i}. {name}</h2>{status}<a href="{r["path"]}"><img loading="lazy" src="{r["path"]}" alt="{name}"></a></section>')
    (OUT/'figure_gallery.html').write_text('\n'.join(blocks),encoding='utf8')
    (OUT/'figure_manifest.json').write_text(json.dumps(figures,indent=2),encoding='utf8')
    selected={
        '20261107':['pilot_results.json','protocol.json','selected_controls.json'],
        'checks-20261108':['modern_differential.json','modern_summary.json','checks_summary.json','scan_comparison_summary.json'],
        'independent-20261108':['asassn_replication.json','replication_analysis_specification.json','protocol.json','references.json','report.txt'],
        'dust-20261108':['dust_analysis.json','event_comparison.json','evidence_dossier.txt','analysis_specification.json','protocol.json'],
        'models-20261108':['model_comparison.json','event_comparison.json','report.txt','protocol.json','model_review.txt'],
        'spherex-20261108':['report.txt','spectrum_comparison.json','method_audit.txt','optical_state_match.json','analysis_specification.json'],
        'physical-tests-20261108':['report.txt','protocol.json','independent_review.txt'],
        'period-20261108':['period_results.json','report.txt','protocol.json'],
        'activity-rhythm-20261108':['rhythm_check.json','protocol.txt'],
        'seasonal-rotation-20261108':['results.json','protocol.txt'],
        'amplitude-brightness-20261108':['results.json','protocol.txt'],
        'short-rotation-20261108':['results.json','protocol.txt'],
        'fast-bluewing-20261108':['bluewing_results.json','protocol.txt'],
        'flare-energy-20261108':['energy_estimate.json'],
        'spots-check-20261108':['spots_only_check.json'],
        'cloud-budget-20261108':['conditional_budget.json'],
        'swift-catalog-check-20261108':['lightcurve_check.json','live_obslc.txt','live_sslc.txt'],
        'optical-overlay-20261108':['comparison_receipt.json','source.json'],
        'gemini-archives-20261108':['archive_audit.txt','ps1_assessment.json','tess_coverage_assessment.json']
    }
    support=OUT/'supporting_results'
    for suffix,names in selected.items():
        folder=ROOT/f'results/sky-plates-j0754-{suffix}'
        for name in names:
            p=folder/name
            if not p.is_file():raise FileNotFoundError(p)
            target=support/folder.name/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
    # All small extracted numerical tables, plus metadata at the deepest key reductions.
    for folder in sorted((ROOT/'results').glob('sky-plates-j0754*')):
        if folder==OUT:continue
        for p in folder.rglob('*.csv'):
            if p.stat().st_size>10_000_000:continue
            target=support/p.relative_to(ROOT/'results');target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
    for rel in ['physical-tests-20261108/tess/sector46_pixel_assessment.json',
                'models-20261108/short_period/fixed_period_results.json',
                'models-20261108/short_period/report.md',
                'dust-20261108/literature/nss_full.json']:
        p=ROOT/f'results/sky-plates-j0754-{rel}'
        if p.is_file():
            target=support/p.relative_to(ROOT/'results');target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
    code=OUT/'source';code.mkdir(exist_ok=True)
    for p in SOURCE.iterdir():
        if p.is_file() and p.suffix in ('.py','.md','.json','.txt'):
            shutil.copy2(p,code/p.name)
    native=[]
    for folder in sorted((ROOT/'results').glob('sky-plates-j0754*')):
        if folder==OUT:continue
        for p in folder.rglob('*'):
            if p.is_file() and any(s in p.name.lower() for s in ('.fits','.fit.gz')):
                native.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=digest(p)))
    (OUT/'native_FITS_inventory.json').write_text(json.dumps(native,indent=2))
    readme=(f'J0754 SCIENTIFIC REVIEW PACKET\n\nStart with J0754_report.pdf: four pages of findings followed by nine figure pages. '
        f'J0754_report.txt contains the same narrative and captions as editable text.\n'
        f'figure_gallery.html provides all {len(figures)} unique original-resolution analysis PNGs; duplicate source paths are in figure_manifest.json. '
        'Some figures record earlier failed or superseded reductions. Published third-party poster/paper pages are excluded; source links remain in the report.\n'
        'supporting_results contains selected analysis JSON/protocol/report files and extracted CSV tables. source contains maintained research scripts and documentation. '
        'Native telescope FITS are listed with SHA256 in native_FITS_inventory.json rather than duplicated in this compact packet. '
        'This is a review packet, not a complete standalone raw-data archive; analysis commands refer to the original outer results/ layout. '
        'Calibration files, native IPAC tables and other downloaded spectra/images remain local and are not all bundled.\n'
        'Report prepared 8 November 2026 local time. No contact, submission or publication has occurred. '
        'Scientific statements are the report authors\' analysis of public observations, not an independent peer review.\n')
    (OUT/'README.txt').write_text(readme,encoding='utf8')
    manifest=[]
    for p in sorted(OUT.rglob('*')):
        if p.is_file() and p.name not in ('packet_manifest.json','J0754_review_packet.zip') and 'qa' not in p.relative_to(OUT).parts:
            manifest.append(dict(path=p.relative_to(OUT).as_posix(),bytes=p.stat().st_size,sha256=digest(p)))
    (OUT/'packet_manifest.json').write_text(json.dumps(manifest,indent=2))
    archive=OUT/'J0754_review_packet.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for r in manifest:z.write(OUT/r['path'],r['path'])
        z.write(OUT/'packet_manifest.json','packet_manifest.json')
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    print(json.dumps(dict(output=str(OUT),figures=len(figures),native_fits_files=len(native),
                         packaged_files=len(manifest)+1,zip_megabytes=round(archive.stat().st_size/1e6,2))))


def main():
    OUT.mkdir(exist_ok=True)
    content=json.loads((SOURCE/'j0754_report_content.json').read_text(encoding='utf8'))
    txt=[content['title'],content['subtitle'],content['credit']]
    for i,p in enumerate(content['pages'],1):
        txt.extend(['',f'PAGE {i}  {p["title"]}'])
        for head,body in p['sections']:txt.extend(['',head,body])
    txt.append('\nFIGURE CAPTIONS')
    for no,rel,title,caption in FIGURES:txt.extend(['',f'Figure {no}  {title}',caption])
    (OUT/'J0754_report.txt').write_text('\n'.join(txt),encoding='utf8')
    build_pdf(content)
    package_support(content)


if __name__=='__main__':main()
