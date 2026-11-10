"""Check 2024/25 rotational cancellation in shorter, predeclared windows."""
import csv
import json
from pathlib import Path
import numpy as np
from astropy.time import Time
from orbit_photometry import daily, fit, PERIOD

OUT = Path('results/sky-plates-j0754-short-rotation-20261108')


def check(rows, name):
    t = np.array([float(r['t']) for r in rows]); y = np.array([float(r['y']) for r in rows])
    e = np.maximum(.02, [float(r['err']) for r in rows]); g = np.array([r['group'] for r in rows])
    start = Time('2024-09-01').jd; mid = Time('2025-01-15').jd; end = Time('2025-06-01').jd
    windows = [('whole', start, end), ('first_half', start, mid), ('second_half', mid, end)]
    for width in (120, 90):
        for offset in range(0, int(end-start-width)+1, 30):
            windows.append((f'rolling{width}_{offset}', start+offset, start+offset+width))
    result = []
    rng = np.random.default_rng(20261109)
    for label, lo, hi in windows:
        m = (t >= lo) & (t < hi)
        entry = dict(window=label, start=Time(lo,format='jd').isot[:10],
                     stop=Time(hi,format='jd').isot[:10], n=int(sum(m)))
        if sum(m) < 12 or len(np.unique(np.floor(t[m]-.5))) < 10 or np.ptp(t[m]) < 50:
            entry['status'] = 'insufficient_cadence'; result.append(entry); continue
        a = fit(t[m],y[m],e[m],g[m],mode='group_linear')
        if sum(m)-a['nuisance_rank']-2 < 8 or min(a['phase_energy_retained_range']) < .3:
            entry['status'] = 'poor_geometry'; result.append(entry); continue
        p = a['parameters'][0]
        # Shared signs on seven-day blocks preserve simultaneous camera residuals.
        block = np.floor((t[m]-lo)/7).astype(int)
        signs = rng.choice([-1.,1.],size=(block.max()+1,512))
        residual = a['_null_residual']
        beta_null = a['_op']@(residual[:,None]*signs[block])
        amp_null = 2*np.linalg.norm(beta_null,axis=0)
        entry.update(status='fit', peak_to_peak_mag=2*p['semi_amplitude_mag'],
            peak_to_peak_white_sigma=2*p['amplitude_white_sigma'],
            phase_cycles=p['cosine_phase_cycles'], beta=a['beta'],
            null_peak_to_peak_p95=float(np.quantile(amp_null,.95)),
            null_fraction_ge_observed=float(np.mean(amp_null >= 2*p['semi_amplitude_mag'])),
            phase_energy_retained=a['phase_energy_retained_range'])
        result.append(entry)
    return result


def main():
    OUT.mkdir(exist_ok=True)
    protocol = ('Exploratory follow-up: fixed35.9204566days, no period optimization. '
        'Whole2024Sep1â€“2025Jun1 and two halves split2025Jan15; additionally90/120day '
        'windows startingSep1 every30days. Joint camera offsets/slopes, errors floor.02mag, '
        '>=12groups,>=10nights,>=50days,>=8residualdof,minimum phase energy retention.3. '
        'Compare amplitudes/phases; strong short-window amplitudes with offset phases '
        'would support cancellation of a changing signal. Weak/unconstrained shorter fits '
        'do not establish rapid spot evolution. Shared7day block-sign512draw nulls are '
        'diagnostic only, not calibrated FAP; overlapping windows are not independent '
        'and pointwise tails are not trial-corrected. White errors omit correlated noise. '
        'No model of arbitrary faster-than-window evolution or upper limit on all variability.')
    (OUT/'protocol.txt').write_text(protocol,encoding='utf8')
    with open('results/sky-plates-j0754-models-20261108/short_period/asassn_g_daily.csv') as f:
        asas=list(csv.DictReader(f))
    raw=json.loads(Path('results/sky-plates-j0754-checks-20261108/modern_differential.json').read_text())['measurements']
    ztf=daily([dict(t=r['mjd']+2400000.5,y=r['corrected_relative_mag'],err=.02,group='ZTF') for r in raw])
    result=dict(protocol=protocol,ASASSN=check(asas,'ASASSN'),ZTF=check(ztf,'ZTF'))
    (OUT/'results.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    for name in ('ASASSN','ZTF'):
        print(name)
        for r in result[name]:
            print(r['window'],r['n'],r['status'],
                {k:round(r[k],4) for k in ('peak_to_peak_mag','peak_to_peak_white_sigma','phase_cycles','null_peak_to_peak_p95','null_fraction_ge_observed') if k in r})


if __name__=='__main__': main()
