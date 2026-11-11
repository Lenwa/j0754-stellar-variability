"""Bounded public MAST inventory refresh and 7x7 TESSCut downloads.

No accounts or tokens. Reuse cached products; never replace prior FITS.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile
import requests
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-flare-search-20261110'
RA,DEC=118.69117410685,16.694621231895

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--sectors',default='45,72');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);pixelout=OUT/'pixels';pixelout.mkdir(exist_ok=True)
    records=[]
    receipt=OUT/'archive_refresh.json'
    old=json.loads(receipt.read_text()) if receipt.exists() else []
    for name,url,params in [
        ('current_tesscut_coverage','https://mast.stsci.edu/tesscut/api/v0.1/sector',dict(ra=RA,dec=DEC)),
        ('current_mast_cone','https://mast.stsci.edu/api/v0/invoke',dict(request=json.dumps(dict(service='Mast.Caom.Cone',params=dict(ra=RA,dec=DEC,radius=3/3600),format='json',pagesize=1000,page=1))))]:
        try:
            response=requests.get(url,params=params,timeout=60);response.raise_for_status();data=response.json()
            (OUT/(name+'.json')).write_text(json.dumps(data,indent=2),encoding='utf8')
            records.append(dict(kind=name,url=response.url,status='retrieved',rows=len(data.get('data',data.get('results',[])))))
            print(name,records[-1]['rows'],flush=True)
        except (requests.RequestException,ValueError) as error:
            records.append(dict(kind=name,status='failed',error=str(error)));print(name,str(error),flush=True)
    for sector in map(int,args.sectors.split(',')):
        existing=list(pixelout.glob(f'tess-s{sector:04d}-*astrocut.fits'))
        if existing:
            payload=None;paths=existing
        else:
            url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut'
            params=dict(ra=RA,dec=DEC,x=7,y=7,units='px',sector=sector)
            try:
                response=requests.get(url,params=params,timeout=180);response.raise_for_status()
                if len(response.content)>40_000_000:
                    raise ValueError('Unexpected cutout exceeds 40 MB bound')
                with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
                    paths=[]
                    for info in archive.infolist():
                        if not info.filename.endswith('.fits'):
                            continue
                        payload=archive.read(info)
                        with fits.open(io.BytesIO(payload)) as hdus:
                            if int(hdus[0].header['SECTOR'])!=sector or hdus[1].data['FLUX'].shape[1:]!=(7,7):
                                raise ValueError('Unexpected sector or pixel dimensions')
                        path=pixelout/Path(info.filename).name
                        if not path.exists():path.write_bytes(payload)
                        paths.append(path)
                if not paths:raise ValueError('No FITS in archive response')
            except (requests.RequestException,ValueError,zipfile.BadZipFile,OSError) as error:
                records.append(dict(kind='pixels',sector=sector,status='failed',error=str(error)));print(sector,str(error),flush=True);continue
        for path in paths:
            records.append(dict(kind='pixels',sector=sector,status='downloaded_or_cached',file=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),url=f'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={RA}&dec={DEC}&x=7&y=7&units=px&sector={sector}'))
            print(sector,path.name,path.stat().st_size,flush=True)
    # Preserve earlier successful downloads when fetching additional sectors.
    requested=set(map(int,args.sectors.split(',')))
    records.extend(r for r in old if r.get('kind')=='pixels' and r.get('sector') not in requested)
    receipt.write_text(json.dumps(records,indent=2),encoding='utf8')

if __name__=='__main__':main()
