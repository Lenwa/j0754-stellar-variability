"""Fetch bounded public CK04/SVO inputs; native data stay outside source control."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import datetime
import json
from pathlib import Path
import requests

FILTERS={'g':('filter_g.xml','SLOAN/SDSS.g'),'W1':('filter_W1.xml','WISE/WISE.W1'),
         'W2':('filter_W2.xml','WISE/WISE.W2'),'B':('filter_B.xml','Generic/Johnson.B'),
         'V':('filter_V.xml','Generic/Johnson.V'),'r':('filter_r.xml','SLOAN/SDSS.r'),
         'I':('filter_CousinsI.xml','Generic/Cousins.I'),'J':('filter_J.xml','2MASS/2MASS.J'),
         'H':('filter_H.xml','2MASS/2MASS.H'),'Ks':('filter_Ks.xml','2MASS/2MASS.Ks'),
         'i':('filter_sdss_i.xml','SLOAN/SDSS.i'),'z':('filter_z.xml','SLOAN/SDSS.z'),
         'ztfg':('filter_ztfg.xml','Palomar/ZTF.g')}


def main(out):
    out.mkdir(exist_ok=True,parents=True)
    urls={f'{m}_{t}.fits':f'https://ssb.stsci.edu/trds/grid/ck04models/{m}/{m}_{t}.fits'
          for m in ['ckm10','ckm05','ckp00'] for t in range(3500,6501,250)}
    urls.update({name:'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID='+fid for name,fid in FILTERS.values()})
    def one(item):
        name,url=item;path=out/name
        if not path.exists():
            r=requests.get(url,timeout=45)
            (out/(name+'.request.json')).write_text(json.dumps(dict(url=r.url,status=r.status_code,
                bytes=len(r.content),retrieved_utc=datetime.datetime.now(datetime.UTC).isoformat()),indent=2)+'\n')
            r.raise_for_status();path.write_bytes(r.content)
        return dict(file=name,bytes=path.stat().st_size)
    records=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(one,item) for item in urls.items()]):records.append(future.result())
    (out/'downloads.json').write_text(json.dumps(sorted(records,key=lambda r:r['file']),indent=2)+'\n')
    print(len(records),'native atmosphere/filter files ready')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    main(p.parse_args().out)
