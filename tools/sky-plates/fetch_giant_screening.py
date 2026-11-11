"""Bounded public catalogue/photometry retrieval for the comparison-giant screen."""
import argparse
import hashlib
import io
import json
import re
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from html.parser import HTMLParser
from pathlib import Path
import requests
from screen_giant_analogs import ROOT, OUT, RECON, shortlist


class TableReader(HTMLParser):
    def __init__(self):
        super().__init__();self.tables=[];self.table=None;self.row=None;self.cell=None;self.skip=0
    def handle_starttag(self,tag,attrs):
        if tag=='table':self.table=[]
        elif tag=='tr' and self.table is not None:self.row=[]
        elif tag in ('td','th') and self.row is not None:self.cell=[]
        elif tag in ('annotation','annotation-xml'):self.skip+=1
    def handle_endtag(self,tag):
        if tag in ('annotation','annotation-xml'):self.skip=max(0,self.skip-1)
        elif tag in ('td','th') and self.cell is not None:
            self.row.append(' '.join(''.join(self.cell).split()));self.cell=None
        elif tag=='tr' and self.row is not None:
            self.table.append(self.row);self.row=None
        elif tag=='table' and self.table is not None:
            self.tables.append(self.table);self.table=None
    def handle_data(self,data):
        if self.cell is not None and not self.skip:self.cell.append(data)


def scientific_tables(raw):
    parser=TableReader();parser.feed(raw);return parser.tables


def dump(path,value):path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n',encoding='utf8')


def fetch(name,url,params=None,transform=None,method='GET'):
    path=OUT/name;receipt_path=OUT/(name+'.receipt.json')
    if path.exists() and receipt_path.exists():
        receipt=json.loads(receipt_path.read_text(encoding='utf8'))
        if hashlib.sha256(path.read_bytes()).hexdigest()!=receipt.get('sha256'):raise ValueError('Cached file hash mismatch: '+name)
        return receipt
    record=dict(file=name,url=url,params=params,method=method,retrieved_utc=datetime.now(timezone.utc).isoformat())
    if receipt_path.exists():
        previous=json.loads(receipt_path.read_text(encoding='utf8'))
        history=previous.pop('previous_attempts',[])
        record['previous_attempts']=history+[previous]
    try:
        response=requests.request(method,url,params=params if method=='GET' else None,json=params if method=='POST' else None,timeout=(12,40))
        record.update(http_status=response.status_code,response_sha256=hashlib.sha256(response.content).hexdigest(),response_bytes=len(response.content))
        response.raise_for_status()
        if len(response.content)>12_000_000:raise ValueError('Response exceeds 12 MB bound')
        raw=response.content
        if transform:raw=transform(response.text)
        if isinstance(raw,str):raw=raw.encode('utf8')
        if isinstance(raw,(dict,list)):raw=(json.dumps(raw,indent=2,ensure_ascii=False)+'\n').encode('utf8')
        if raw is not None:
            path.write_bytes(raw);record.update(sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),state='saved')
        else:record.update(state='metadata_only_no_download')
        if b'Error=Table or Catalog not found' in response.content:record['state']='catalogue_not_found'
    except (requests.RequestException,ValueError,OSError) as error:
        record.update(state='failed',error=str(error))
    dump(receipt_path,record);print(name,record['state'],record.get('bytes',''),flush=True)
    return record


def public_tables(raw):
    if 'captcha' in raw.lower() or 'bot manager' in raw.lower():raise ValueError('Public article protected by human verification; no bypass attempted')
    tables=scientific_tables(raw)
    if not tables:raise ValueError('No scientific tables in article response')
    return dict(tables=tables)


def catalogues():
    OUT.mkdir(parents=True,exist_ok=True)
    protocol=OUT/'protocol.json'
    if not protocol.exists():
        plan=ROOT/'docs/superpowers/plans/2026-11-11-giant-screening.md'
        dump(protocol,dict(created_utc=datetime.now(timezone.utc).isoformat(),status='Exploratory catalogue selection; J0754 properties and published benchmark already known.',plan=plan.relative_to(ROOT).as_posix(),plan_sha256=hashlib.sha256(plan.read_bytes()).hexdigest(),selection='Five nearest published giant-labelled rows in Teff/radius/luminosity; archive outcomes excluded from rank.',limitations=['Not a prevalence study or a claim of neglected/unpublished stars.','Dixon 7286 sample is screened only if actual rows are retrieved.','2022 catalogue carries later physical parameters alongside earlier luminosity classifications; retain conflicts.'] ))
    base='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
    jobs=[
        ('olah2022_stars.tsv',base,{'-source':'J/A+A/668/A101/table1','-out.max':'10000','-out.all':'1'},None),
        ('olah2022_flares.tsv',base,{'-source':'J/A+A/668/A101/table2','-out.max':'10000','-out.all':'1'},None),
        ('dixon2025_vizier.tsv',base,{'-source':'J/AJ/169/309','-out.max':'10000'},None),
        ('olah2021_tables.json','https://arxiv.org/html/2010.07623v2',None,public_tables),
        ('benchmark2026_tables.json','https://arxiv.org/html/2610.11799v1',None,public_tables),
        ('dixon2025_journal_tables.json','https://iopscience.iop.org/article/10.3847/1538-3881/adc92a',None,public_tables)]
    with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda j:fetch(*j),jobs))
    rows,selected=shortlist();print('CATALOGUE',len(rows),'SHORTLIST',[(r['KIC'],round(r['physical_score'],3)) for r in selected],flush=True)


def targets():
    _,selected=shortlist()
    result=[dict(name='KIC'+r['KIC'],ra=float(r['_RA']),dec=float(r['_DE']),rank=r['physical_rank']) for r in selected]
    result.extend([dict(name='J0736',ra=114.02711824252,dec=21.23634538304),dict(name='J0830',ra=127.66033651549,dec=14.12030855354)])
    dump(OUT/'targets.json',result);return result


def archives(mast_only=False):
    jobs=[]
    for target in targets():
        name,ra,dec=target['name'],target['ra'],target['dec']
        query=f"SELECT TOP 200 source_id,ra,dec,pmra,pmdec,parallax,parallax_error,phot_g_mean_mag,bp_rp,ruwe FROM gaiadr3.gaia_source WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{60/3600}))=1"
        jobs.append((name+'_gaia.json','https://gea.esac.esa.int/tap-server/tap/sync',dict(REQUEST='doQuery',LANG='ADQL',FORMAT='json',QUERY=query),None))
        request=dict(service='Mast.Caom.Cone',params=dict(ra=ra,dec=dec,radius=5/3600),format='json',pagesize=1000,page=1)
        jobs.append((name+'_mast.json','https://mast.stsci.edu/api/v0/invoke',dict(request=json.dumps(request)),None))
        jobs.append((name+'_ztf.csv','https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves',dict(CIRCLE=f'{ra} {dec} {2/3600}',FORMAT='csv'),None))
    # Existing Tang orbit identifiers are checked again with all diagnostic columns.
    query='SELECT * FROM gaiadr3.nss_two_body_orbit WHERE source_id IN (865283696503303936,651437519650162304,667401913089980544)'
    jobs.append(('tang_gaia_orbits.json','https://gea.esac.esa.int/tap-server/tap/sync',dict(REQUEST='doQuery',LANG='ADQL',FORMAT='json',QUERY=query),None))
    if mast_only:jobs=[job for job in jobs if job[0].endswith('_mast.json') and not (OUT/job[0]).exists()]
    with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda j:fetch(*j),jobs))


def products():
    for target in targets():
        path=OUT/(target['name']+'_mast.json')
        if not path.exists():continue
        metadata=json.loads(path.read_text(encoding='utf8'))
        rows=metadata.get('data',[])
        kepler=[r for r in rows if r.get('obs_collection')=='Kepler' and r.get('dataproduct_type')=='timeseries' and '_lc_' in r.get('obs_id','')]
        if not kepler:continue
        # One LC observation groups all quarter products; short-cadence products excluded later.
        record=sorted(kepler,key=lambda r:str(r.get('obsid')))[0]
        request=dict(service='Mast.Caom.Products',params=dict(obsid=str(record['obsid'])),format='json',pagesize=1000,page=1)
        receipt=fetch(target['name']+'_kepler_products.json','https://mast.stsci.edu/api/v0/invoke',dict(request=json.dumps(request)))
        if receipt['state']!='saved':continue
        data=json.loads((OUT/receipt['file']).read_text(encoding='utf8'))
        lc=sorted({r['dataURI']:r for r in data.get('data',[]) if str(r.get('productFilename','')).endswith('_llc.fits') and int(r.get('size') or 0)<5_000_000}.values(),key=lambda r:r['productFilename'])
        choices=sorted(set([min(4,len(lc)-1),min(16,len(lc)-1)])) if lc else []
        for index in choices:
            r=lc[index];file=r['productFilename']
            fetch(file,'https://mast.stsci.edu/api/v0.1/Download/file',dict(uri=r['dataURI']))


def coverage():
    """Post-screen metadata follow-up when final Q17 cannot cover two known cycles."""
    from astropy.io import fits
    import numpy as np
    from screen_giant_analogs import literature_rows
    periods=literature_rows();selected=[]
    for target in targets():
        if not target['name'].startswith('KIC'):continue
        kic=target['name'][3:];period=periods[kic]['published_rotation_days']
        if period is None:continue
        short=False
        for path in OUT.glob('kplr'+kic.zfill(9)+'-*_llc.fits'):
            with fits.open(path) as h:
                if h[0].header.get('QUARTER')==17:
                    t=h[1].data['TIME'];t=t[np.isfinite(t)]
                    short=len(t)>0 and np.ptp(t)<2*period
        if not short:continue
        data=json.loads((OUT/(target['name']+'_kepler_products.json')).read_text())['data']
        row=next((r for r in data if r.get('description','').endswith(' - Q15') and r['productFilename'].endswith('_llc.fits')),None)
        if row:
            fetch(row['productFilename'],'https://mast.stsci.edu/api/v0.1/Download/file',dict(uri=row['dataURI']));selected.append(target['name'])
    dump(OUT/'coverage_followup.json',dict(reason='After initial inspection, append Q15 when final Q17 cannot cover two cycles at the already published rotation period. Retain Q17. This is an explicit metadata/coverage follow-up, not a blinded selection.',targets=selected,selection='Q15 long-cadence products from previously downloaded product tables'))


def siblings():
    """Retain the immediately useful sibling light curves for the next analysis."""
    selected=[]
    for name in ['J0736','J0830']:
        metadata=json.loads((OUT/(name+'_mast.json')).read_text())
        for row in metadata.get('data',[]):
            uri=str(row.get('dataURL') or '')
            qlp=uri.startswith('mast:HLSP/qlp/') and row.get('sequence_number') in [44,45,46]
            k2=name=='J0830' and str(row.get('obs_id'))=='ktwo211569759-c05_lc'
            if not (qlp or k2):continue
            file=uri.rsplit('/',1)[-1]
            receipt=fetch(file,'https://mast.stsci.edu/api/v0.1/Download/file',dict(uri=uri))
            selected.append(dict(target=name,file=file,uri=uri,collection=row.get('obs_collection'),sequence=row.get('sequence_number'),state=receipt['state'],analysis_status='Downloaded, not reduced in this screening'))
    dump(OUT/'sibling_lightcurve_products.json',selected)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['catalogues','archives','retry-mast','products','coverage','siblings']);args=parser.parse_args()
    if args.stage=='retry-mast':archives(mast_only=True)
    else:globals()[args.stage]()
