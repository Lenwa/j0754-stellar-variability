"""Cross-check the complete saved 7x7 TESS pixel series with a direct aperture."""
import hashlib,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import analyze_tess_variability as tv
import search_tess_flares as sf

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-flare-search-20261110'

def main():
    paths=sorted((OUT/'pixels').glob('*astrocut.fits'))+list((ROOT/'results/sky-plates-j0754-physical-tests-20261108/tess/pixels').glob('*astrocut.fits'))
    records=[];fig,axes=plt.subplots(len(paths),1,figsize=(14,3*len(paths)),layout='constrained')
    for path,ax in zip(paths,np.atleast_1d(axes)):
        with fits.open(path,memmap=False) as hdus:
            d=hdus[1].data;sector=int(hdus[0].header['SECTOR']);t=np.asarray(d['TIME'],float);cube=np.asarray(d['FLUX'],float)
            xy=WCS(hdus[2].header).all_world2pix([[118.69117410685,16.694621231895]],0)[0]
            yy,xx=np.indices(cube.shape[1:]);edge=(xx==0)|(yy==0)|(xx==6)|(yy==6);core=(xx-xy[0])**2+(yy-xy[1])**2<=1.5**2
            flux=np.sum(cube[:,core],axis=1)-np.median(cube[:,edge],axis=1)*sum(core.ravel())
            mask=(d['QUALITY']==0)&np.isfinite(t)&np.isfinite(flux)&(flux>0)&np.all(np.isfinite(cube),axis=(1,2));good=np.flatnonzero(mask);good=good[np.argsort(t[good])]
            nt=t[good];ny=flux[good]/np.median(flux[good])-1;b=tv.bin_series(nt,ny)
            detections=[]
            for method,tt,y in [('native',nt,ny),('20min',b['t'],b['y'])]:
                scan=sf.screen(tt,y)
                if method=='20min':bscan=scan
                for event in scan['events']:
                    detections.append(dict(event,method=method,start_btjd=float(tt[event['start_index']]),stop_btjd=float(tt[event['stop_index']-1])))
            events=sf.merged_events(detections)
            for event in events:
                peak=event['peak_btjd'];window=[max(event['start_btjd'],peak-.08),min(event['stop_btjd'],peak+.08)]
                if window[1]-window[0]<.02:window=[peak-.04,peak+.04]
                try:
                    diagnostic=tv.pixel_difference(t,cube,d['QUALITY'],window,xy,guard=.15,side=.3)
                    event['pixel_offset']=diagnostic['difference_centroid_offset_pixels']
                    event['stationary_scale']=diagnostic['stationary_only_fractional_scale']
                except ValueError as error:event['pixel_failure']=str(error)
                gaps=np.diff(nt)>2.5*np.median(np.diff(nt));boundary=np.r_[nt[0],nt[-1],nt[:-1][gaps],nt[1:][gaps]]
                event['gap_or_edge_hours']=float(min(abs(boundary-peak))*24)
                event['known_screen_match']=next((label for center,label in [(2541.3366,'new_S45'),(2572.0677,'original_S46'),(3277.3968,'new_S72')] if abs(center-peak)<.3),None)
            ax.plot(b['t'],bscan['residual']*100,'.',ms=2)
            for e in events:ax.axvline(e['peak_btjd'],color='orange',lw=.8)
            ax.set(title=f'S{sector} direct 1.5-pixel core minus outer-ring background | q=0',ylabel='Residual (%)')
            records.append(dict(sector=sector,file=path.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),accepted_rows=len(nt),noise_20min=bscan['noise'],events=events))
            np.savetxt(OUT/f'sector{sector}_direct_pixel_20min.csv',np.column_stack([b['t'],b['y'],bscan['baseline'],bscan['residual']]),delimiter=',',header='btjd,relative_flux,baseline,residual',comments='')
    axes[-1].set_xlabel('BTJD; orange: screen candidates including systematic effects')
    fig.savefig(OUT/'direct_pixel_search.png',dpi=130);plt.close(fig)
    (OUT/'direct_pixel_search.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),products=records),indent=2,allow_nan=False),encoding='utf8')
    print(json.dumps([dict(sector=r['sector'],events=r['events']) for r in records],indent=2))

if __name__=='__main__':main()
