"""Compact, source-linked candidate catalogue and comparison figure."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
from astropy.time import Time
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-flare-search-20261110'

def main():
    followup=json.loads((OUT/'pixel_followup.json').read_text());catalogue=[]
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    for event,ax in zip(followup['events'],axes.ravel()):
        sector=event['sector'];center=np.mean(event['event_window_btjd'])
        data=np.genfromtxt(OUT/f'sector{sector}_direct_pixel_20min.csv',delimiter=',',names=True)
        near=abs(data['btjd']-center)<.2;index=np.flatnonzero(near)[np.argmax(data['residual'][near])]
        peak=float(data['btjd'][index]);date=Time(peak+2457000,format='jd',scale='tdb').utc.isot
        use=abs(data['btjd']-peak)<.65;x=(data['btjd'][use]-peak)*24;y=data['residual'][use]*100
        gap=np.flatnonzero(np.diff(x)>1.)+1
        ax.plot(np.insert(x,gap,np.nan),np.insert(y,gap,np.nan),'.-',ms=4,lw=.8)
        known=sector==46;amplitude=float(data['residual'][index]);offset=event['difference_centroid_offset_pixels']
        ax.axhline(0,color='gray',lw=.6);ax.set(xlabel='Hours from observed 20-min peak',ylabel='Direct aperture excess (%)',
            title=f'{date[:10]} | S{sector} | {"original event" if known else "additional candidate"}\nPeak {amplitude*100:.2f}% | excess offset {offset:.2f} pixel')
        catalogue.append(dict(date_utc=date,sector=sector,peak_btjd=peak,
            peak_20min_direct_fractional=amplitude,positive_excess_offset_pixels=offset,
            positive_excess_offset_arcsec=offset*21,
            event_mean_stationary_fractional_scale=event['stationary_only_fractional_scale'],
            event_mean_scale_with_motion=event['stationary_plus_motion_model']['fractional_stationary_scale'],
            status='Previously localized event' if known else 'Additional source-centered flare-like candidate',
            note='Peak depends on aperture, baseline and binning; excess moment has no calibrated position uncertainty.'))
    fig.suptitle('J0754: three additional source-centered flare-like brightenings\nSame direct pixel aperture and 2-day median baseline; each event has its own vertical scale')
    fig.savefig(OUT/'four_event_lightcurves.png',dpi=160);plt.close(fig)
    (OUT/'candidate_catalogue.json').write_text(json.dumps(catalogue,indent=2),encoding='utf8')
    with (OUT/'candidate_catalogue.csv').open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fieldnames=list(catalogue[0]));writer.writeheader();writer.writerows(catalogue)
    # Add every retrieved pixel file, including files cached in earlier requests.
    receipt=json.loads((OUT/'archive_refresh.json').read_text())
    files={r.get('file') for r in receipt if r.get('kind')=='pixels' and r.get('status')=='downloaded_or_cached'}
    for path in sorted((OUT/'pixels').glob('*astrocut.fits')):
        relative=path.relative_to(ROOT).as_posix()
        if relative not in files:
            with fits.open(path) as hdus:sector=int(hdus[0].header['SECTOR'])
            receipt.append(dict(kind='pixels',sector=sector,status='downloaded_or_cached',file=relative,bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                url=f'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=118.69117410685&dec=16.694621231895&x=7&y=7&units=px&sector={sector}'))
    (OUT/'archive_refresh.json').write_text(json.dumps(receipt,indent=2),encoding='utf8')
    print(json.dumps(catalogue,indent=2))

if __name__=='__main__':main()
