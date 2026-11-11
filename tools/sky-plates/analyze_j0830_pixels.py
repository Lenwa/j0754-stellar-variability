"""Direct calibrated-image extraction and bounded J0830 attribution test."""
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_tang_siblings import ROOT,TARGETS,PERIODS
from screen_giant_analogs import write_csv,gap_series
from tang_sibling_core import absolute_time,bin_series,periodic_fit,scan_periods,clean
from j0830_pixel_core import background_planes,aperture_flux,differential_flux,periodic_image,spatial_components,register_on_control

OUT=ROOT/'results/j0830-pixel-check-20261111'
PERIOD=TARGETS['J0830']['period']


def world_positions(rows):
    return [[r['ra']+(r.get('pmra') or 0)*5.8/3.6e6/np.cos(np.deg2rad(r['dec'])),r['dec']+(r.get('pmdec') or 0)*5.8/3.6e6] for r in rows]


def gaussian_center(image,xy):
    """Descriptive spatial fit; no calibrated positional confidence interval."""
    yy,xx=np.indices(image.shape);mask=(xx-xy[0])**2+(yy-xy[1])**2<=3.5**2
    x,y,z=xx[mask],yy[mask],image[mask]
    def residual(p):
        a,cx,cy,sx,sy,base,gx,gy=p
        return a*np.exp(-.5*((x-cx)**2/sx**2+(y-cy)**2/sy**2))+base+gx*(x-xy[0])+gy*(y-xy[1])-z
    scale=max(float(np.max(z)-np.min(z)),1e-9)
    fit=least_squares(residual,[scale,*xy,.9,.9,float(np.median(z)),0,0],bounds=([0,xy[0]-2,xy[1]-2,.35,.35,-np.inf,-np.inf,-np.inf],[np.inf,xy[0]+2,xy[1]+2,2.5,2.5,np.inf,np.inf,np.inf]),x_scale='jac')
    return dict(xy=fit.x[1:3].tolist(),offset_from_wcs_pixels=float(np.linalg.norm(fit.x[1:3]-xy)),rms=float(np.sqrt(np.mean(fit.fun**2))),success=bool(fit.success),width_xy_pixels=fit.x[3:5].tolist())


def extract(selection):
    records=[];parts=[];products=json.loads((OUT/'pixel_products.json').read_text())
    for product in products:
        if product['state']!='saved':continue
        path=ROOT/product['file']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=product['sha256']:raise ValueError('Native hash mismatch')
        with fits.open(path,memmap=False) as h:
            d=h[1].data;sector=int(h[0].header['SECTOR']);group='S'+str(sector);wcs=WCS(h[2].header)
            target_xy=wcs.all_world2pix(world_positions([selection['target']]),0)[0]
            control_xy=wcs.all_world2pix(world_positions(selection['selected_controls']),0)
            field_xy=wcs.all_world2pix(world_positions(selection['all_field_sources']),0)
            t=absolute_time(d['TIME'],h[1].header);cube=np.asarray(d['FLUX'],float)
            good=(np.asarray(d['QUALITY'])==0)&np.isfinite(t)&np.isfinite(cube).all(axis=(1,2));t=t[good];cube=cube[good]
            plane=background_planes(cube,field_xy);corrected=cube-plane
            control=np.column_stack([aperture_flux(corrected,xy,1.5) for xy in control_xy])
            target={str(radius):aperture_flux(corrected,target_xy,radius) for radius in [1.,1.5,2.]}
            common=np.isfinite(control).all(axis=1)&(control>0).all(axis=1)
            for values in target.values():common&=np.isfinite(values)&(values>0)
            t,cube,corrected,plane,control=[v[common] for v in [t,cube,corrected,plane,control]];target={k:v[common] for k,v in target.items()}
            g=np.full(len(t),group);variants={}
            aligned,offsets=register_on_control(corrected,control_xy[0])
            if np.max(abs(offsets))>1.5:raise ValueError('Motion exceeds declared crop margin; attribution remains unresolved')
            registered_control=aperture_flux(aligned,control_xy[0],1.5)
            if np.any(~np.isfinite(registered_control)) or np.any(registered_control<=0):raise ValueError('Aligned comparison unavailable')
            for radius,flux in target.items():
                differential,ensemble=differential_flux(flux,control)
                variants['aperture_r'+radius]=flux/np.median(flux)-1
                variants['differential_r'+radius]=differential/np.median(differential)-1
                registered=aperture_flux(aligned,target_xy,float(radius));registered_diff,_=differential_flux(registered,registered_control[:,None])
                variants['registered_r'+radius]=registered/np.median(registered)-1
                variants['registered_differential_r'+radius]=registered_diff/np.median(registered_diff)-1
                if radius=='1.5':
                    variants['ensemble']=ensemble/np.median(ensemble)-1
                    for j in range(control.shape[1]):
                        one,_=differential_flux(flux,control[:,[j]])
                        variants['only_control_'+str(j+1)]=one/np.median(one)-1
            for j in range(control.shape[1]):variants['control_'+str(j+1)]=control[:,j]/np.median(control[:,j])-1
            variants['background']=np.median(plane,axis=(1,2))
            binned={name:bin_series(t,values,g) for name,values in variants.items()}
            # The identical native rows and bins are used for target, controls and maps.
            image_bins=bin_series(t,np.zeros(len(t)),g,extra=corrected.reshape(len(t),-1))
            bt=image_bins['t'];images=image_bins['extra'].reshape(-1,*cube.shape[1:])
            aligned_bins=bin_series(t,np.zeros(len(t)),g,extra=aligned[:,2:-2,2:-2].reshape(len(t),-1));aligned_images=aligned_bins['extra'].reshape(-1,17,17)
            if not np.all(np.isfinite(aligned_images)):raise ValueError('Registered map has invalid interior pixels')
            for name,b in binned.items():
                if not np.array_equal(b['t'],bt):raise ValueError('Mismatched control/target sampling')
            record=dict(sector=sector,file=path.relative_to(ROOT).as_posix(),native_sha256=product['sha256'],camera=int(h[0].header['CAMERA']),ccd=int(h[0].header['CCD']),quality_zero_finite=int(sum(good)),common_accepted=len(t),bins=len(bt),first_bjd=float(min(t)),last_bjd=float(max(t)),target_wcs_xy=target_xy.tolist(),control_wcs_xy=control_xy.tolist(),target_median_flux_e_per_s={k:float(np.median(v)) for k,v in target.items()},control_median_flux_e_per_s=np.median(control,axis=0).tolist(),registration_max_abs_offset_pixels=np.max(abs(offsets),axis=0).tolist(),registration_source_id=selection['selected_controls'][0]['source_id'])
            rows=[dict(bjd=float(ti),sector=group,**{name:float(b['y'][i]) for name,b in binned.items()}) for i,ti in enumerate(bt)]
            write_csv(OUT/f'sector{sector}_binned_photometry.csv',rows)
            parts.append(dict(t=bt,g=np.full(len(bt),group),variants={name:b['y'] for name,b in binned.items()},images=images,registered_images=aligned_images,registered_reference=np.median(aligned[:,2:-2,2:-2],axis=0),median_image=np.median(corrected,axis=0),raw_median=np.median(cube,axis=0),xy=target_xy,controls=control_xy,sector=sector))
            records.append(record)
    if len(parts)!=3:raise ValueError('The bounded test requires all three requested sectors')
    return records,parts


def main():
    selection=json.loads((OUT/'field_selection.json').read_text());records,parts=extract(selection)
    t=np.concatenate([p['t'] for p in parts]);g=np.concatenate([p['g'] for p in parts]);variants={name:np.concatenate([p['variants'][name] for p in parts]) for name in parts[0]['variants']}
    fits_out={};curves={}
    for name,y in variants.items():
        if name=='background':continue
        fits_out[name]={}
        for mode in ['global_linear','sector_linear','arc_offset']:
            fit=periodic_fit(t,y,g,PERIOD,mode);peak,power=scan_periods(t,y,g,PERIODS,mode)
            fits_out[name][mode]=dict(at_gaia=clean(fit),free_peak=peak)
            write_csv(OUT/f'{name}_{mode}_periodogram.csv',[dict(period_days=float(pi),power=float(vi)) for pi,vi in zip(PERIODS,power)])
            curves[(name,mode)]=fit
    reference=curves[('aperture_r1.5','sector_linear')]['beta'];spatial=[];fig,axes=plt.subplots(3,3,figsize=(12,11),layout='constrained')
    for p,record,row in zip(parts,records,axes):
        maps=periodic_image(p['t'],p['images'],p['g'],PERIOD,reference,'sector_linear');image=maps['signed'];xy=p['xy'];center=gaussian_center(image,xy);baseline_center=gaussian_center(p['median_image'],xy)
        stationary=spatial_components(p['median_image'],image,xy)
        spatial.append(dict(sector=p['sector'],gaussian_periodic_center=center,gaussian_reference_center=baseline_center,periodic_minus_reference_pixels=float(np.linalg.norm(np.array(center['xy'])-baseline_center['xy'])),spatial_components=stationary,signal_condition=maps['signal_condition']))
        for ax,im,title in zip(row,[p['median_image'],image,maps['quadrature']],['Background-subtracted reference','Signed 38.389-day modulation','Quadrature modulation']):
            limit=np.max(abs(im));artist=ax.imshow(im,origin='lower',cmap='viridis' if title.startswith('Background') else 'RdBu_r',vmin=None if title.startswith('Background') else -limit,vmax=None if title.startswith('Background') else limit)
            ax.plot(*xy,'wx',ms=10,mew=2)
            for j,cxy in enumerate(p['controls']):ax.plot(*cxy,'o',mfc='none',mec='yellow',ms=7);ax.text(cxy[0]+.3,cxy[1],str(j+1),color='yellow',fontsize=8)
            if title.startswith('Signed'):ax.plot(*center['xy'],'k+',ms=10,mew=2)
            ax.set(title=f'S{p["sector"]}: {title}',xlabel='Pixel x',ylabel='Pixel y');fig.colorbar(artist,ax=ax,shrink=.7,label='e-/s/pixel')
        np.savez_compressed(OUT/f'sector{p["sector"]}_maps.npz',reference=p['median_image'],signed=image,quadrature=maps['quadrature'],cosine=maps['cosine'],sine=maps['sine'])
    fig.suptitle('J0830 spatial test: white cross = Gaia position; circles = fixed comparison stars\nSigned maps use the same stellar phase; black plus = descriptive Gaussian center')
    fig.savefig(OUT/'spatial_attribution.png',dpi=150);plt.close(fig)
    reference_registered=curves[('registered_r1.5','sector_linear')]['beta'];registered_spatial=[]
    fig,axes=plt.subplots(3,2,figsize=(9,11),layout='constrained')
    for p,row in zip(parts,axes):
        maps=periodic_image(p['t'],p['registered_images'],p['g'],PERIOD,reference_registered,'sector_linear');xy=p['xy']-2;template=p['registered_reference'];center=gaussian_center(maps['signed'],xy);basecenter=gaussian_center(template,xy)
        registered_spatial.append(dict(sector=p['sector'],gaussian_periodic_center=center,gaussian_reference_center=basecenter,periodic_minus_reference_pixels=float(np.linalg.norm(np.array(center['xy'])-basecenter['xy'])),spatial_components=spatial_components(template,maps['signed'],xy),signal_condition=maps['signal_condition']))
        for ax,im,title in zip(row,[template,maps['signed']],['Aligned reference','Signed 38.389-day modulation']):
            limit=np.max(abs(im));artist=ax.imshow(im,origin='lower',cmap='viridis' if title.startswith('Aligned') else 'RdBu_r',vmin=None if title.startswith('Aligned') else -limit,vmax=None if title.startswith('Aligned') else limit);ax.plot(*xy,'wx',ms=10,mew=2)
            if title.startswith('Signed'):ax.plot(*center['xy'],'k+',ms=10,mew=2)
            ax.set(title=f'S{p["sector"]}: {title}',xlabel='Cropped pixel x',ylabel='Cropped pixel y');fig.colorbar(artist,ax=ax,shrink=.75,label='e-/s/pixel')
        np.savez_compressed(OUT/f'sector{p["sector"]}_registered_maps.npz',reference=template,signed=maps['signed'],quadrature=maps['quadrature'],cosine=maps['cosine'],sine=maps['sine'])
    fig.suptitle('One motion-aware check: align images using bright comparison 1 only\nWhite cross = target Gaia WCS; black plus = descriptive modulation center')
    fig.savefig(OUT/'registered_spatial_attribution.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(11,9),layout='constrained')
    for ax,names,title in [(axes[0],['aperture_r1.5','differential_r1.5'],'Direct and comparison-corrected J0830'),(axes[1],['control_1','control_2','ensemble'],'Comparison stars on identical cadences'),(axes[2],['only_control_1','only_control_2'],'Each comparison used separately')]:
        for name in names:
            fit=curves[(name,'sector_linear')]
            for sector in sorted(set(g)):
                m=g==sector;ax.plot(*gap_series(t[m]-min(t),(variants[name]-fit['_nuisance'])[m]*100,.5),'.',ms=2,label=name if sector==min(g) else None)
            ax.plot(t-min(t),fit['_periodic']*100,lw=1,alpha=.6)
        ax.set(title=title,xlabel='Days within 2021 epoch',ylabel='Sector-linear-corrected flux (%)');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Independent aperture extraction of the same TESS exposures; no QLP reference-image offsets')
    fig.savefig(OUT/'direct_and_comparison_lightcurves.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(11,7),layout='constrained')
    for ax,name in zip(axes,['registered_r1.5','registered_differential_r1.5']):
        fit=curves[(name,'sector_linear')]
        for sector in sorted(set(g)):
            m=g==sector;ax.plot(*gap_series(t[m]-min(t),(variants[name]-fit['_nuisance'])[m]*100,.5),'.',ms=3,label=sector);ax.plot(*gap_series(t[m]-min(t),fit['_periodic'][m]*100,.5),color='black',lw=1)
        ax.set(title=name.replace('_',' '),xlabel='Days within 2021 epoch',ylabel='Sector-linear-corrected flux (%)');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('After comparison-star image alignment: direct and bright-control differential flux')
    fig.savefig(OUT/'registered_lightcurves.png',dpi=150);plt.close(fig)
    result=dict(question='Is the five-week signal spatially associated with J0830?',records=records,control_source_ids=[r['source_id'] for r in selection['selected_controls']],models=fits_out,spatial=spatial,post_initial_motion_check=dict(status='Added after unregistered maps showed shared dipoles; one bounded check, not an optimized family',method='Independent bright-control positive-moment centroid within fixed 2.5px radius; cubic image translation; crop two edge pixels for maps',registered_spatial=registered_spatial),reference_phase_beta=reference,limitations=['Calibrated TESSCut full-frame pixels are another reduction of the same TESS exposures, not an independent instrument.','Only about two proposed cycles; free scans do not supply formal period precision.','Spatial centers are descriptive fits, not calibrated position intervals; unresolved blends remain possible.','Controls were selected from catalogue geometry before light-curve inspection; both retained even if variable.','Signed spatial maps are conditional on the fixed period and primary aperture phase.','Interpolation can change undersampled peak pixels; aperture and spatial results are model-dependent.'])
    (OUT/'results.json').write_text(json.dumps(clean(result),indent=2,allow_nan=False)+'\n',encoding='utf8')
    sources=[Path(__file__),Path(__file__).with_name('j0830_pixel_core.py'),Path(__file__).with_name('fetch_j0830_pixels.py'),Path(__file__).with_name('tang_sibling_core.py'),OUT/'field_selection.json']+[ROOT/r['file'] for r in records]
    (OUT/'input_source_hashes.json').write_text(json.dumps({p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},indent=2)+'\n',encoding='utf8')
    for name in ['aperture_r1.5','differential_r1.5','control_1','control_2','only_control_1','only_control_2']:
        print(name,[(m,round(v['free_peak']['best_period_days'],3),round(v['at_gaia']['harmonic_semiamplitudes'][0]*100,3),round(v['at_gaia']['power'],3)) for m,v in fits_out[name].items()],flush=True)
    print('SPATIAL',[(r['sector'],r['gaussian_periodic_center']['offset_from_wcs_pixels'],r['periodic_minus_reference_pixels'],r['spatial_components']['stationary_fraction']) for r in spatial],flush=True)
    print('REGISTERED_SPATIAL',[(r['sector'],r['gaussian_periodic_center']['offset_from_wcs_pixels'],r['periodic_minus_reference_pixels'],r['spatial_components']['stationary_fraction']) for r in registered_spatial],flush=True)
    for name in ['registered_r1.5','registered_differential_r1.5']:print(name,[(m,round(v['free_peak']['best_period_days'],3),round(v['at_gaia']['harmonic_semiamplitudes'][0]*100,3),round(v['at_gaia']['power'],3)) for m,v in fits_out[name].items()],flush=True)


if __name__=='__main__':main()
