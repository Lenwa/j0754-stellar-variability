"""Bounded public TESSCut retrieval for the approved J0830 source test."""
import hashlib
import io
import json
import zipfile
from datetime import datetime,timezone
from pathlib import Path
import requests
from astropy.io import fits
from analyze_tang_siblings import ROOT,TARGETS
import fetch_giant_screening as retrieval

OUT=ROOT/'results/j0830-pixel-check-20261111'
RA,DEC=TARGETS['J0830']['ra'],TARGETS['J0830']['dec']


def main():
    OUT.mkdir(parents=True,exist_ok=True);pixels=OUT/'pixels';pixels.mkdir(exist_ok=True)
    protocol=OUT/'protocol.json'
    if not protocol.exists():
        protocol.write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),question='Does the 38–39-day signal originate at J0830 rather than nearby stars or shared field effects?',prediction='Signed periodic image concentrates at target; fixed apertures and comparison-corrected target retain the wave; controls do not explain it.',weakening_outcome='Displaced source, a shared field wave, or no direct-pixel reproduction.',sectors=[44,45,46],size_pixels=21,primary_period_days=TARGETS['J0830']['period'],aperture_radii_pixels=[1.,1.5,2.],background='Per-frame plane excluding catalogue sources, robust residual clipping',control_selection='G<=14.5; >=4 pixels from target, >=3 pixels from other G<=16.5 sources, >=2.5 pixels from cutout edge; choose up to five by G proximity before time-series inspection.',quality='QUALITY=0, finite positive aperture flux; common target/control cadences',native_size_bound_bytes=60000000,scope='Three sectors only; no broad new archive search, orbit refit or significance calibration.'),indent=2)+'\n',encoding='utf8')
    retrieval.OUT=OUT
    query=f"SELECT TOP 5000 source_id,ra,dec,pmra,pmdec,phot_g_mean_mag,bp_rp FROM gaiadr3.gaia_source WHERE phot_g_mean_mag<17.5 AND CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},0.1))=1"
    retrieval.fetch('gaia_field.json','https://gea.esac.esa.int/tap-server/tap/sync',dict(REQUEST='doQuery',LANG='ADQL',FORMAT='json',QUERY=query))
    records=[]
    for sector in [44,45,46]:
        cached=list(pixels.glob(f'tess-s{sector:04d}-*astrocut.fits'))
        receipt_path=OUT/f'sector{sector}_pixels.receipt.json'
        if cached and receipt_path.exists():
            record=json.loads(receipt_path.read_text())
            for p in cached:
                if hashlib.sha256(p.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('Cached cutout hash mismatch')
            records.append(record);print(sector,'cached',flush=True);continue
        url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut';params=dict(ra=RA,dec=DEC,x=21,y=21,units='px',sector=sector)
        record=dict(sector=sector,url=url,params=params,method='GET',retrieved_utc=datetime.now(timezone.utc).isoformat())
        try:
            response=requests.get(url,params=params,timeout=(12,120));record.update(http_status=response.status_code,response_bytes=len(response.content),response_sha256=hashlib.sha256(response.content).hexdigest());response.raise_for_status()
            if len(response.content)>60_000_000:raise ValueError('Cutout response exceeds 60 MB bound')
            with zipfile.ZipFile(io.BytesIO(response.content)) as z:
                members=[i for i in z.infolist() if i.filename.endswith('.fits')]
                if len(members)!=1 or members[0].file_size>60_000_000:raise ValueError('Unexpected product count/size')
                payload=z.read(members[0])
            with fits.open(io.BytesIO(payload)) as h:
                if int(h[0].header['SECTOR'])!=sector or h[1].data['FLUX'].shape[1:]!=(21,21):raise ValueError('Wrong sector or shape')
            path=pixels/Path(members[0].filename).name
            if path.exists():raise ValueError('Refusing to overwrite native file')
            path.write_bytes(payload);record.update(state='saved',file=path.relative_to(ROOT).as_posix(),bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest())
        except (requests.RequestException,ValueError,OSError,zipfile.BadZipFile) as e:record.update(state='failed',error=str(e))
        receipt_path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf8');records.append(record);print(sector,record['state'],record.get('bytes',record.get('error')),flush=True)
    (OUT/'pixel_products.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf8')


if __name__=='__main__':main()
