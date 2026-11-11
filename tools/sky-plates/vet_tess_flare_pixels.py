"""Post-screen pixel follow-up: fixed windows, sideband sensitivity and maps."""
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import analyze_tess_variability as tv

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-flare-search-20261110'
EVENTS=[(45,[2534.52,2534.60]),(45,[2541.30,2541.38]),(46,[2572.0260444438954,2572.1371604654987]),(72,[3277.30,3277.48])]

def main():
    records=[]
    fig,axes=plt.subplots(len(EVENTS),4,figsize=(17,4*len(EVENTS)),layout='constrained')
    for (sector,event),row in zip(EVENTS,axes):
        directory=OUT/'pixels' if sector!=46 else ROOT/'results/sky-plates-j0754-physical-tests-20261108/tess/pixels'
        path=next(directory.glob(f'tess-s{sector:04d}-*astrocut.fits'))
        with fits.open(path,memmap=False) as hdus:
            data=hdus[1].data; t=np.array(data['TIME'],float);cube=np.array(data['FLUX'],float);q=np.array(data['QUALITY'])
            xy=WCS(hdus[2].header).all_world2pix([[118.69117410685,16.694621231895]],0)[0]
            base=tv.pixel_difference(t,cube,q,event,xy,guard=.15,side=.3)
            sensitivity=[]
            for guard,side in [(.1,.1),(.15,.3),(.3,.4),(.5,.4)]:
                try:
                    r=tv.pixel_difference(t,cube,q,event,xy,guard=guard,side=side)
                    sensitivity.append(dict(guard_days=guard,side_days=side,offset_pixels=r['difference_centroid_offset_pixels'],
                        centroid_xy=r['difference_positive_centroid_xy'],stationary_scale=r['stationary_only_fractional_scale'],
                        scale_with_motion=r['stationary_plus_motion_model']['fractional_stationary_scale']))
                except ValueError as error:sensitivity.append(dict(guard_days=guard,side_days=side,error=str(error)))
            maps=base['_maps']; yy,xx=np.indices(cube.shape[1:]);edge=(xx==0)|(yy==0)|(xx==xx.max())|(yy==yy.max());core=(xx-xy[0])**2+(yy-xy[1])**2<=1.5**2
            good=base['_good'];peak=np.mean(event);local=good&(abs(t-peak)<1.)
            before,eventmask,after=base['_masks'];times=base['sideband_mean_btjd']
            coreflux=np.sum(cube[:,core],axis=1)-np.median(cube[:,edge],axis=1)*np.sum(core)
            border=np.median(cube[:,edge],axis=1)
            anchor=[float(np.mean(coreflux[before])),float(np.mean(coreflux[after]))]
            line=anchor[0]+(anchor[1]-anchor[0])*(t-times[0])/(times[2]-times[0])
            rel=coreflux/line-1
            centroidseries=[]
            # One-hour frame averages reduce positive-moment noise bias.
            for center in np.arange(peak-.3,peak+.301,1/24):
                use=good&(abs(t-center)<1/48)
                if sum(use)<3:continue
                actual=float(np.mean(t[use]));alpha=(actual-times[0])/(times[2]-times[0])
                baseline=(1-alpha)*maps['before']+alpha*maps['after'];difference=np.mean(cube[use],axis=0)-baseline
                plane=np.column_stack([np.ones(xx.size),xx.ravel(),yy.ravel()]);co=np.linalg.lstsq(plane[edge.ravel()],difference[edge],rcond=None)[0]
                corrected=difference-(plane@co).reshape(xx.shape);pos=np.maximum(corrected,0);total=float(pos.sum())
                centroid=[float((xx*pos).sum()/total),float((yy*pos).sum()/total)] if total>0 else None
                centroidseries.append(dict(btjd=actual,frames=int(sum(use)),relative_core_flux=float(np.mean(rel[use])),positive_moment_xy=centroid))
            for ax,name in zip(row[:2],['baseline','plane_subtracted_difference']):
                m=maps[name];im=ax.imshow(m,origin='lower',cmap='viridis' if name=='baseline' else 'RdBu_r')
                ax.plot(*xy,'wx',ms=11,mew=2)
                if name!='baseline':ax.plot(*base['difference_positive_centroid_xy'],'k+',ms=12,mew=2)
                ax.set(title=f'S{sector}: {name.replace("_"," ")}',xlabel='Pixel x',ylabel='Pixel y');fig.colorbar(im,ax=ax,label='e-/s')
            row[2].plot((t[local]-peak)*24,rel[local]*100,'.',ms=2,label='Core, border-subtracted')
            row[2].axvspan((event[0]-peak)*24,(event[1]-peak)*24,color='orange',alpha=.25)
            row[2].set(xlabel='Hours from selected window center',ylabel='Local core change (%)',title=f'Excess offset {base["difference_centroid_offset_pixels"]:.2f} pixels');row[2].legend(fontsize=7)
            row[3].plot((t[local]-peak)*24,border[local],'.',ms=2)
            row[3].set(xlabel='Hours from selected window center',ylabel='Outer-ring median (e-/s/pixel)',title='Local pixel background')
            output={k:v for k,v in base.items() if not k.startswith('_')}
            output.update(sector=sector,event_window_btjd=event,file=path.relative_to(ROOT).as_posix(),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),sideband_sensitivity=sensitivity,
                local_core_border_correlation=tv.correlation(rel[local],border[local]),
                hourly_positive_moments=centroidseries,
                note='Post-screen follow-up. White cross is target WCS; black plus is positive excess moment. Large offsets undermine target association. Positive-moment centroids have noise bias and no formal calibrated error.')
            records.append(output)
            table=np.column_stack([t[local],q[local],coreflux[local],border[local],rel[local]])
            np.savetxt(OUT/f'sector{sector}_{int(peak)}_local_pixel_series.csv',table,delimiter=',',header='btjd,quality,core_flux_e_per_s,border_e_per_s_per_pixel,relative_core_flux',comments='')
    fig.suptitle('New screen candidates compared with the known S46 event | same pixel procedure')
    fig.savefig(OUT/'pixel_candidate_comparison.png',dpi=150);plt.close(fig)
    (OUT/'pixel_followup.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),events=records),indent=2,allow_nan=False),encoding='utf8')
    print(json.dumps([dict(sector=r['sector'],offset_pixels=r['difference_centroid_offset_pixels'],xy=r['difference_positive_centroid_xy'],target=r['target_wcs_pixel_xy'],scale=r['stationary_only_fractional_scale'],sensitivity=r['sideband_sensitivity']) for r in records],indent=2))

if __name__=='__main__':main()
