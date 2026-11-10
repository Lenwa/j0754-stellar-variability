"""Bounded native SPHEREx cutouts from an already captured SIA2 inventory."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import csv
import datetime
import hashlib
import json
from pathlib import Path

import requests


def run(metadata,out,size=.06):
    out.mkdir(exist_ok=True,parents=True)
    rows=list(csv.DictReader(metadata.open()))
    def one(r):
        det=int(r['energy_bandpassname'].split('D')[-1])
        name=f"{r['obs_id']}D{det}"
        path=out/(name+'.fits');receipt=out/(name+'.request.json')
        record=dict(obs_id=r['obs_id'],detector=det,parent_url=r['access_url'],
                    t_min=r['t_min'],t_max=r['t_max'],file=path.name)
        if path.exists():
            raw=path.read_bytes()
            if not raw.startswith(b'SIMPLE  ='):raise ValueError('Invalid cached FITS: '+str(path))
            return dict(record,status='CACHED',bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
        if receipt.exists():
            prior=out/(name+'.prior_'+datetime.datetime.now(datetime.UTC).strftime('%Y%m%dT%H%M%S')+'.json')
            prior.write_bytes(receipt.read_bytes())
        try:
            response=requests.get(r['access_url'],params=dict(center='118.69117410685,16.694621231895',size=size),stream=True,timeout=(15,90))
            record.update(url=response.url,http_status=response.status_code,
                          retrieved_utc=datetime.datetime.now(datetime.UTC).isoformat())
            response.raise_for_status();parts=[];total=0
            for part in response.iter_content(262144):
                total+=len(part)
                if total>15_000_000:raise ValueError('Cutout exceeds 15MB cap')
                parts.append(part)
            raw=b''.join(parts)
            if not raw.startswith(b'SIMPLE  ='):raise ValueError('Response is not FITS')
            path.write_bytes(raw)
            record.update(status='DOWNLOADED',bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
        except (requests.RequestException,ValueError) as exc:
            record.update(status='FAILED',error=str(exc))
        receipt.write_text(json.dumps(record,indent=2)+'\n')
        return record
    records=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for i,f in enumerate(as_completed([pool.submit(one,r) for r in rows]),1):
            r=f.result();records.append(r)
            if i%12==0 or r['status']=='FAILED':print(i,'/',len(rows),r['file'],r['status'],flush=True)
    (out/'downloads.json').write_text(json.dumps(sorted(records,key=lambda r:r['file']),indent=2)+'\n')
    print('Completed',len(records),'cutouts;',sum(r['status']=='FAILED' for r in records),'failures',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--metadata',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--size',type=float,default=.06)
    a=p.parse_args();run(a.metadata,a.out,a.size)
