"""Seasonal check at the fixed Gaia period; reuse preserved daily photometry."""
import csv
import json
from pathlib import Path
import numpy as np
from astropy.time import Time
from orbit_photometry import fit, PERIOD

OUT=Path('results/sky-plates-j0754-activity-rhythm-20261108')
INPUT=Path('results/sky-plates-j0754-models-20261108/short_period/asassn_g_daily.csv')


def main():
    OUT.mkdir(exist_ok=True)
    protocol='Exploratory follow-up: fixed35.9204566d, separate July-June observing seasons. Fit camera offsets and within-camera linear trends jointly with sine/cosine using existing fit(); errors floored0.02mag. At least30 daily groups and90day span. No period selection or photometric clipping. Annual bq-camera medians show slow trend separately. White-noise-only errors do not establish seasonal detections or physical temperature. Earlier frozen-waveform held-out test controls persistent-clock claim. No claim of long activity-cycle period from one modern dip.'
    (OUT/'protocol.txt').write_text(protocol,encoding='utf8')
    rows=list(csv.DictReader(INPUT.open()))
    t=np.array([float(r['t']) for r in rows]);y=np.array([float(r['y']) for r in rows])
    e=np.maximum(.02,np.array([float(r['err']) for r in rows]));g=np.array([r['group'] for r in rows])
    years=np.array([int(r['day'][:4])-(int(r['day'][5:7])<7) for r in rows])
    results=[]
    for year in sorted(set(years)):
        m=years==year
        if sum(m)<30 or np.ptp(t[m])<90:continue
        a=fit(t[m],y[m],e[m],g[m],mode='group_linear')
        bq=m&(g=='bq');q=a['parameters'][0]
        results.append(dict(season=f'{year}/{str(year+1)[2:]}',n=int(sum(m)),
            start=Time(min(t[m]),format='jd').isot[:10],stop=Time(max(t[m]),format='jd').isot[:10],
            bq_camera_median_mag=float(np.median(y[bq])) if sum(bq)>=10 else None,
            bq_groups=int(sum(bq)),**q,
            retention=a['phase_energy_retained_range'],residual_rms_mag=a['residual_rms_mag']))
    p=json.load(open('results/sky-plates-j0754-models-20261108/short_period/fixed_period_results.json'))
    r=dict(period_days=PERIOD,input=str(INPUT),method=protocol,seasonal=results,
        earlier_heldout=p['asassn_g']['spline365']['heldout'] if 'heldout' in p['asassn_g']['spline365'] else 'See preserved fixed_period_results.json and report.md',
        limits='No independent temperature measurements or companion confirmation. Seasonal amplitude is a fit at an imposed period, not a searched/detected rotation period. Phase is poorly defined when amplitude is small. Long dip and evolving spot model are compatible but not causally identified.')
    (OUT/'rhythm_check.json').write_text(json.dumps(r,indent=2),encoding='utf8')
    print(json.dumps(r,indent=2))


if __name__=='__main__':main()
