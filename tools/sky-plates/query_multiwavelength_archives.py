"""Public metadata reconnaissance at J0754; no accounts and no bulk download.

Every response is saved with exact query, retrieval time and hash. Empty tables,
catalogue matches, pointing coverage and failed requests remain distinct.
"""
import argparse,csv,gzip,hashlib,io,json,time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
import requests
from astropy.io.votable import parse_single_table
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-multiwavelength-20261110'
RA,DEC=118.69117410685,16.694621231895
HEA='https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync'

def json_row_count(data):
    if isinstance(data,list):
        tables=[t for t in data if isinstance(t,dict) and t.get('TableName')!='SqlQuery' and 'Rows' in t]
        return sum(len(t['Rows']) for t in tables) if tables else None
    if isinstance(data,dict):
        rows=data.get('data')
        if isinstance(rows,list):return len(rows)
        if isinstance(rows,dict) and 'Tables' in rows:return sum(len(t.get('Rows',[])) for t in rows['Tables'])
    return None

def tap(query):return dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query)

def cone(table,radius):
    return f"SELECT TOP 300 * FROM {table} WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},{radius}))=1"

def jobs():
    result=[]
    def add(name,url,params):result.append(dict(name=name,url=url,params=params))
    add('heasarc_tables',HEA,tap("SELECT table_name,description FROM TAP_SCHEMA.tables WHERE table_name LIKE '%uvot%' OR table_name LIKE '%xmm%' OR table_name LIKE '%chan%' OR table_name LIKE '%csc%' OR table_name LIKE '%eros%'"))
    for table in ['xmmmaster','chanmaster','swiftmastr']:
        add(table+'_pointings_1deg',HEA,tap(cone(table,1.)))
    for table in ['xmmssc','csc','swuvotssc']:
        add(table+'_sources_60arcsec',HEA,tap(cone(table,60/3600)))
    for table in ['xmmslewcln','xmmslewful']:
        add(table+'_sources_60arcsec',HEA,tap(cone(table,60/3600)))
    add('heasarc_slew_radio_tables',HEA,tap("SELECT table_name,description FROM TAP_SCHEMA.tables WHERE table_name LIKE '%slew%' OR table_name LIKE '%vlass%'"))
    for release in ['DR1_Main','DR1_Supp','DR2_Main','DR2_Hard']:
        add('erosita_'+release,'https://erosita.mpe.mpg.de/erodat/catalogue/SCS',dict(CAT=release,RA=RA,DEC=DEC,SR=120/3600,VERB=3))
    for release in ['DR1_eRASS1','DR2_eRASSc3']:
        add('erosita_limit_'+release,'https://erosita.mpe.mpg.de/erodat/upperlimit/service',dict(ra=RA,dec=DEC,band='024',dr_survey=release))
    for release,table in [('dr17','apogeeStar'),('dr19','apogee_drp_allstar')]:
        query=f'SELECT TOP 20 * FROM {table} WHERE ra BETWEEN {RA-.003} AND {RA+.003} AND dec BETWEEN {DEC-.003} AND {DEC+.003}'
        add('apogee_'+release,'https://skyserver.sdss.org/'+release+'/SkyServerWS/SearchTools/SqlSearch',dict(cmd=query,format='json'))
    for name,service,params in [
        ('mast_observations','Mast.Caom.Cone',dict(ra=RA,dec=DEC,radius=3/3600)),
        ('galex_products','Mast.Caom.Products',dict(obsid='23671')),
        ('galex_catalog','Mast.Galex.Catalog',dict(ra=RA,dec=DEC,radius=30/3600))]:
        add(name,'https://mast.stsci.edu/api/v0/invoke',dict(request=json.dumps(dict(service=service,params=params,format='json',pagesize=1000,page=1))))
    query=f"SELECT TOP 100 o.collection,o.observationID,p.productID,p.planeID,p.time_bounds_lower,p.time_bounds_upper,p.position_resolution FROM caom2.Observation AS o JOIN caom2.Plane AS p ON o.obsID=p.obsID WHERE o.collection LIKE 'VLASS%' AND INTERSECTS(p.position_bounds,CIRCLE('ICRS',{RA},{DEC},0.001))=1"
    add('vlass_cadc_coverage','https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync',tap(query))
    add('gphoton_nuv_ranges','https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest',
        dict(query=f"select startTimeRange,endTimeRange from GPFCore.dbo.fGetTimeRanges(1,10000000000000,{RA},{DEC}) where band='NUV'",format='extjs'))
    add('gphoton_nuv_aspect','https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest',
        dict(query=f"select distinct time from GPFCore.dbo.fGetNearbyAspectEq({RA},{DEC},37.5,1000,1200000000000) where band='NUV' or band='FUV/NUV' order by time",format='extjs'))
    add('swift_uvot_image_listing','https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/2009_12/00031554001/uvot/image/',{})
    # Adaptive follow-ups use IDs returned by the first-stage public queries.
    aspect=OUT/'gphoton_nuv_aspect.json'
    if aspect.exists():
        data=json.loads(aspect.read_text())
        if data.get('status')=='COMPLETE':
            times=[int(r[0]) for r in data['data']['Tables'][0]['Rows']]
            if times:
                add('gphoton_nuv_aperture_photons','https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest',
                    dict(query=f'select time,ra,dec,xi,eta,x,y,flag,q from GPFCore.dbo.fGetNearbyObjEqNUVAllColumns({RA},{DEC},{6/3600},{min(times)-1000},{max(times)+1000},0)',format='extjs'))
    coverage=OUT/'vlass_cadc_coverage.native'
    if coverage.exists():
        rows=list(csv.DictReader(coverage.read_text().splitlines()))
        ids=[r['planeID'] for r in rows if 'planeID' in r]
        if ids:
            quoted=','.join("'"+v+"'" for v in ids)
            add('vlass_artifacts','https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync',tap(f'SELECT planeID,uri,contentType,contentLength FROM caom2.Artifact WHERE planeID IN ({quoted})'))
    return result

def fetch(job):
    job=dict(job,params=dict(job['params']))
    if job['url']==HEA:
        # This TAP service rejects FORMAT=csv; default VOTable is documented.
        job['params'].pop('FORMAT',None)
    name=job['name'];record=dict(job,retrieved_utc=datetime.now(timezone.utc).isoformat())
    path=OUT/(name+'.native')
    try:
        response=requests.get(job['url'],params=job['params'],timeout=65)
        record.update(status_code=response.status_code,url_resolved=response.url,bytes=len(response.content))
        if len(response.content)>8_000_000:raise ValueError('Metadata response exceeds8MB bound')
        path.write_bytes(response.content);record['file']=path.relative_to(ROOT).as_posix();record['sha256']=hashlib.sha256(response.content).hexdigest()
        response.raise_for_status()
        raw=response.content.lstrip()
        if raw.startswith((b'{',b'[')):
            data=response.json();(OUT/(name+'.json')).write_text(json.dumps(data,indent=2),encoding='utf8')
            if isinstance(data,dict) and str(data.get('status','')).upper()=='ERROR':
                raise ValueError('Archive JSON error: '+str(data.get('msg',data)))
            record['state']='json_response';record['json_keys']=list(data) if isinstance(data,dict) else None
            count=json_row_count(data)
            if count is not None:record['rows']=count
        elif b'<VOTABLE' in raw[:1500] or b'<vot:VOTABLE' in raw[:1500]:
            root=ET.fromstring(raw)
            errors=[node.text for node in root.iter() if node.tag.endswith('INFO') and node.attrib.get('name')=='QUERY_STATUS' and node.attrib.get('value')=='ERROR']
            if errors:raise ValueError('Archive query error: '+'; '.join(errors))
            table=parse_single_table(io.BytesIO(raw)).to_table();table.write(OUT/(name+'.ecsv'),format='ascii.ecsv',overwrite=True)
            record.update(state='table',rows=len(table),columns=list(table.colnames))
        elif b'<html' in raw[:1000].lower() or b'<!doctype html' in raw[:1000].lower():
            record.update(state='html_response_not_query_result',preview=response.text[:250])
        else:
            lines=response.text.splitlines();parsed=list(csv.DictReader(lines));record.update(state='csv_or_text',rows=len(parsed),preview=response.text[:250])
        print(name,record['state'],record.get('rows'),flush=True)
    except (requests.RequestException,ValueError,OSError,IndexError,ET.ParseError) as error:
        record.update(state='failed',error=str(error));print(name,'FAILED',str(error)[:120],flush=True)
    (OUT/(name+'_receipt.json')).write_text(json.dumps(record,indent=2),encoding='utf8')
    return record

def fetch_small_assets():
    """Small timing/calibrated files and 1.2-arcmin radio cutouts only."""
    assets=[]
    products=json.loads((OUT/'galex_products.json').read_text())['data']
    timing={r['productFilename']:r for r in products if 'scst' in r['productFilename']}
    (OUT/'galex_small_timing_products.json').write_text(json.dumps(list(timing.values()),indent=2),encoding='utf8')
    for product in timing.values():
        assets.append((product['productFilename'],'https://mast.stsci.edu/api/v0.1/Download/file',dict(uri=product['dataURI'])))
    for suffix in ['uw1_sk.img.gz','uw1_ex.img.gz']:
        name='sw00031554001'+suffix
        assets.append((name,'https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/2009_12/00031554001/uvot/image/'+name,{}))
    for artifact in csv.DictReader((OUT/'vlass_artifacts.native').read_text().splitlines()):
        uri=artifact['uri']
        if uri.endswith('.fits') and '.rms.' not in uri:
            name=uri.split('/')[-1].replace('.subim.fits','.j0754-cutout.fits')
            assets.append((name,'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync',dict(id=uri,CIRCLE=f'{RA} {DEC} 0.01')))
    directory=OUT/'native_files';directory.mkdir(exist_ok=True);records=[]
    for name,url,params in assets:
        record=dict(name=name,url=url,params=params);path=directory/name
        try:
            if not path.exists():
                with requests.get(url,params=params,timeout=65,stream=True) as response:
                    response.raise_for_status();pieces=[];size=0
                    for piece in response.iter_content(131072):
                        size+=len(piece)
                        if size>8_000_000:raise ValueError('Download exceeds8MB bound')
                        pieces.append(piece)
                    payload=b''.join(pieces)
                validation=gzip.decompress(payload) if payload.startswith(b'\x1f\x8b') else payload
                with fits.open(io.BytesIO(validation)) as hdus:hdus.verify('exception')
                path.write_bytes(payload)
            with fits.open(path) as hdus:
                record['hdus']=[dict(name=h.name,shape=list(h.data.shape) if h.data is not None else None,date_obs=h.header.get('DATE-OBS'),exposure=h.header.get('EXPOSURE')) for h in hdus]
            record.update(state='downloaded_or_cached',file=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            print(name,record['state'],record['bytes'],flush=True)
        except (requests.RequestException,ValueError,OSError) as error:
            record.update(state='failed',error=str(error));print(name,'FAILED',str(error)[:150],flush=True)
        records.append(record)
    (OUT/'native_downloads.json').write_text(json.dumps(records,indent=2),encoding='utf8')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--only',default='');parser.add_argument('--assets',action='store_true');parser.add_argument('--summarize-cached',action='store_true');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);todo=jobs()
    if args.only:todo=[j for j in todo if j['name'] in args.only.split(',')]
    if not args.summarize_cached:
        with ThreadPoolExecutor(max_workers=5) as pool:
            records=[r.result() for r in as_completed([pool.submit(fetch,j) for j in todo])]
    else:
        for receipt in OUT.glob('*_receipt.json'):
            r=json.loads(receipt.read_text());p=OUT/(r['name']+'.json')
            if p.exists() and r.get('state')!='failed':
                count=json_row_count(json.loads(p.read_text()))
                if count is not None:r['rows']=count
                receipt.write_text(json.dumps(r,indent=2),encoding='utf8')
    allrecords=[json.loads(p.read_text()) for p in sorted(OUT.glob('*_receipt.json'))]
    (OUT/'query_inventory.json').write_text(json.dumps(allrecords,indent=2),encoding='utf8')
    if args.assets:fetch_small_assets()

if __name__=='__main__':main()
