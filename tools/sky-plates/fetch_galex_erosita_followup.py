"""Public follow-up inputs; stream the large Gaia catalogue, retaining target rows only."""
import argparse, gzip, hashlib, io, json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import requests
from astropy.io import fits
from astropy.table import Table, vstack

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-galex-erosita-20261111'
GP='https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest'
UID=311707200169

def save_request(name,url,params=None):
    path=OUT/name
    if path.exists(): return path
    r=requests.get(url,params=params,timeout=180);r.raise_for_status()
    path.write_bytes(r.content)
    receipt=dict(url=url,params=params,retrieved_utc=datetime.now(timezone.utc).isoformat(),bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest())
    (OUT/(name+'.receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    return path

def galex():
    query='select time,ra,dec,xi,eta,x,y,flag,q from GPFCore.dbo.fGetNearbyObjEqNUVAllColumns(118.69117410685,16.694621231895,0.008333333333333333,822929964995,822935977995,0)'
    save_request('photons_30arcsec.json',GP,dict(query=query,format='extjs'))
    save_request('NUV_flat.fits','https://archive.stsci.edu/prepds/gphoton/cal/cal/NUV_flat.fits')
    for visit,(a,b) in enumerate([(822929964995,822930063995),(822935894995,822935978995)]):
        query=f'select time/50 as tick,count(*) as n from (select time from GPFCore.dbo.NUVPhotonsV where time >= {a} and time < {b} and flag=0 union all select time from GPFCore.dbo.NUVPhotonsNULLV where time >= {a} and time < {b}) as events group by time/50 order by tick'
        save_request(f'global_ticks_{visit}.json',GP,dict(query=query,format='extjs'))
    # One scalar query per 20s bin: gPhoton shutter time and empirical global rate.
    for visit,(start,stop) in enumerate([(822929965.495,822930062.495),(822935895.495,822935977.495)]):
        edges=np.arange(start,stop,20).tolist()+[stop]
        if edges[-1]-edges[-2]<10: edges.pop(-2)
        for i,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
            a,b=int(a*1000),int(b*1000)
            query=f'select shutter*0.05 from GPFCore.dbo.fGetNUVShutter({a},{b})'
            save_request(f'shutter_{visit}_{i}.json',GP,dict(query=query,format='extjs'))
            query=f'select (select count(*) from GPFCore.dbo.NUVPhotonsV where time >= {a} and time < {b} and flag=0)+(select count(*) from GPFCore.dbo.NUVPhotonsNULLV where time >= {a} and time < {b})'
            save_request(f'global_{visit}_{i}.json',GP,dict(query=query,format='extjs'))

def read_exact(stream,n):
    parts=[];left=n
    while left:
        b=stream.read(left)
        if not b: raise EOFError(f'Truncated FITS: missing {left} bytes')
        parts.append(b);left-=len(b)
    return b''.join(parts)

def read_header(stream):
    blocks=[]
    while True:
        block=read_exact(stream,2880);blocks.append(block)
        if any(block[i:i+8]==b'END     ' for i in range(0,2880,80)):break
    return fits.Header.fromstring(b''.join(blocks).decode('ascii'))

def erosita():
    dest=OUT/'erosita_gaia_counterparts.ecsv'
    if dest.exists():return
    url='https://erosita.mpe.mpg.de/dr2/AllSkySurveyData_dr2/Catalogues_dr2/RamosM_DR2/eRASSc3_Main_GDR3_Public_27Jul2026.fits.gz'
    # Avoid storing 855MB of unrelated catalogue entries. Scan every row for this UID.
    with requests.get(url,stream=True,timeout=180) as r:
        r.raise_for_status();r.raw.decode_content=False
        with gzip.GzipFile(fileobj=r.raw) as stream:
            primary=read_header(stream)
            size=abs(primary['BITPIX'])//8
            if primary.get('NAXIS',0)==0:size=0
            else:
                for axis in range(1,primary['NAXIS']+1):size*=primary[f'NAXIS{axis}']
            if size:read_exact(stream,size+(-size)%2880)
            header=read_header(stream)
            if header['XTENSION'].strip()!='BINTABLE':raise ValueError('Expected binary table')
            (OUT/'erosita_gaia_original_header.txt').write_text(header.tostring(sep='\n'),encoding='ascii')
            # Published file repeats GDR3_source_id. Keep both with explicit suffix.
            seen=set()
            for i in range(1,header['TFIELDS']+1):
                key=f'TTYPE{i}';name=header[key]
                if name in seen:header[key]=name+'_duplicate'
                seen.add(name)
            rowlen,nrows=header['NAXIS1'],header['NAXIS2'];found=[];scanned=0
            while scanned<nrows:
                n=min(10000,nrows-scanned);payload=read_exact(stream,n*rowlen)
                h=header.copy();h['NAXIS2']=n
                prefix=fits.PrimaryHDU().header.tostring().encode()+h.tostring().encode()
                with fits.open(io.BytesIO(prefix+payload+b'\0'*((-len(payload))%2880)),memmap=False) as hdus:
                    rows=hdus[1].data;mask=np.asarray(rows['UID'])==UID
                    if np.any(mask):found.append(Table(rows[mask]))
                scanned+=n
                if scanned%500000==0:print('eROSITA rows scanned',scanned,flush=True)
            if not found:raise ValueError('No target UID in complete Gaia counterpart catalogue')
            table=vstack(found);table.write(dest,format='ascii.ecsv',overwrite=True)
    (OUT/'erosita_gaia_receipt.json').write_text(json.dumps(dict(url=url,retrieved_utc=datetime.now(timezone.utc).isoformat(),filter_UID=UID,rows_scanned=scanned,rows_selected=len(table),subset_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),scope='All binary-table rows scanned; retained matching UID rows only.'),indent=2)+'\n',encoding='utf8')
    print('eROSITA counterpart rows',len(table),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('archive',choices=['galex','erosita']);args=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    {'galex':galex,'erosita':erosita}[args.archive]()
