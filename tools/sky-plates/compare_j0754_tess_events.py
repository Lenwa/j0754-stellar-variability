"""Common finite-cadence measurements of the four saved J0754 TESS excursions.

Models describe shapes; they do not certify flare physics or identify an
unresolved component. See protocol.json for the exploratory measurement rules.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.ndimage import shift
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import analyze_tess_variability as tv
import search_tess_flares as search
from stellar_models import read_bandpass

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-flare-comparison-20261111'
PRIOR=ROOT/'results/sky-plates-j0754-flare-search-20261110'
EVENTS=[(45,2534.564726409669,[2534.52,2534.60]),
        (45,2541.3292121620702,[2541.30,2541.38]),
        (46,2572.0817252024544,[2572.0260444438954,2572.1371604654987]),
        (72,3277.4457208103445,[3277.30,3277.48])]

def exposure_profile(t,width,peak,rise,decay):
    """Unit-peak two-sided exponential averaged across finite exposures; hours."""
    t,width=np.broadcast_arrays(np.asarray(t,float),np.asarray(width,float))
    if rise<=0 or decay<=0 or np.any(width<=0):raise ValueError('Positive timescales required')
    a=t-peak-width/2;b=t-peak+width/2
    left=rise*(np.exp(np.minimum(b,0)/rise)-np.exp(np.minimum(a,0)/rise))
    right=decay*(np.exp(-np.maximum(a,0)/decay)-np.exp(-np.maximum(b,0)/decay))
    return (left+right)/width

def baseline(t,params,degree):
    return np.polynomial.polynomial.polyval(np.asarray(t)/24,np.asarray(params)[4:5+degree])

def prediction(t,width,params,degree):
    return baseline(t,params,degree)+params[0]*exposure_profile(t,width,*params[1:4])

def fit_event(t,y,width,degree=1,initial=None):
    t,y,width=map(lambda x:np.asarray(x,float),(t,y,width))
    bounds=([0.,-4.,.05,.05]+[-.2]*(degree+1),[.3,4.,16.,24.]+[.2]*(degree+1))
    starts=[initial] if initial is not None else [np.r_[max(.001,float(np.max(y)-np.median(y))),0.,r,d,np.median(y),np.zeros(degree)] for r,d in [(.5,4.),(3.,3.),(1.,10.)]]
    fits_=[least_squares(lambda p:prediction(t,width,p,degree)-y,p,bounds=bounds,ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=400) for p in starts]
    best=min(fits_,key=lambda f:f.fun@f.fun);p=best.x;quiet=1+baseline(p[1],p,degree)
    return dict(params=p.tolist(),degree=degree,success=bool(best.success),rms=float(np.sqrt(np.mean(best.fun**2))),
        amplitude_fractional=float(p[0]/quiet),peak_offset_h=float(p[1]),rise_tau_h=float(p[2]),decay_tau_h=float(p[3]),
        fwhm_h=float(np.log(2)*(p[2]+p[3])),model_1percent_duration_h=float(np.log(100)*(p[2]+p[3])),
        equivalent_duration_s=float(p[0]*(p[2]+p[3])*3600/quiet),
        parameter_at_bound=bool(np.any(np.isclose(p,bounds[0],rtol=0,atol=1e-4)) or np.any(np.isclose(p,bounds[1],rtol=0,atol=1e-4))))

def observed_area(t_hours,relative,max_gap_hours):
    return tv.equivalent_duration(np.asarray(t_hours)/24,np.asarray(relative),max_gap_hours/24)

def energy_factor(wave,response,star_temperature,flare_temperature):
    def band(temp):return np.trapezoid(response*wave**-5/np.expm1(.014387768775039337/(wave*temp)),wave)
    return float((flare_temperature/star_temperature)**4*band(star_temperature)/band(flare_temperature))

def pixel_products():
    dirs=[PRIOR/'pixels',ROOT/'results/sky-plates-j0754-physical-tests-20261108/tess/pixels']
    products={}
    for path in [p for d in dirs for p in d.glob('*astrocut.fits')]:
        with fits.open(path,memmap=False) as h:
            data=h[1].data;cube=np.array(data['FLUX'],float);t=np.array(data['TIME'],float);q=np.array(data['QUALITY']);wcs=WCS(h[2].header)
            good=(q==0)&np.isfinite(t)&np.all(np.isfinite(cube),axis=(1,2));xy=wcs.all_world2pix([[118.69117410685,16.694621231895]],0)[0]
            yy,xx=np.indices(cube.shape[1:]);edge=(xx==0)|(yy==0)|(xx==xx.max())|(yy==yy.max())
            products[int(h[0].header['SECTOR'])]=dict(t=t,cube=cube,quality=q,good=good,xy=xy,xx=xx,yy=yy,edge=edge,
                cadence_days=float(h[1].header['TIMEDEL']),path=path,wcs=wcs,header=h[2].header.copy(),dead_fraction_header=float(h[1].header.get('DEADC',1)))
    return products

def aperture_series(p,radius):
    core=(p['xx']-p['xy'][0])**2+(p['yy']-p['xy'][1])**2<=radius*radius
    flux=np.sum(p['cube'][:,core],axis=1)-np.median(p['cube'][:,p['edge']],axis=1)*np.sum(core)
    mask=p['good']&np.isfinite(flux)&(flux>0)
    return p['t'][mask],flux[mask]

def measurements(p,peak,draws):
    variants=[];central=None;series=None
    for radius in [1.,1.5,2.]:
        t,flux=aperture_series(p,radius);xh=(t-peak)*24;reference=np.median(flux[abs(xh)<36]);y=flux/reference-1
        for degree in [1,2]:
            for window in [[-24,30],[-36,36]]:
                inside=(xh>=window[0])&(xh<=window[1]);x=xh[inside];value=y[inside];width=np.full(len(x),p['cadence_days']*24)
                fit=fit_event(x,value,width,degree);fit.update(aperture_pixels=radius,context_h=window,native_samples=int(sum(inside)))
                base=baseline(x,fit['params'],degree);relative=(value-base)/(1+base)
                a=(x>=-12)&(x<=18)
                fit['observed_signed_ED_s']=observed_area(x[a],relative[a],2.5*p['cadence_days']*24)
                fit['observed_integration_window_h']=[-12,18]
                fit['max_retained_gap_h']=float(np.max(np.diff(x[a])))
                fit['sampled_model_fraction']=float(observed_area(x[a],fit['amplitude_fractional']*exposure_profile(x[a],width[a],*fit['params'][1:4]),2.5*p['cadence_days']*24)/fit['equivalent_duration_s'])
                variants.append(fit)
                if radius==1.5 and degree==1 and window==[-24,30]:
                    central=fit;series=dict(x=x,y=value,width=width,baseline=base,relative=relative,reference_flux=float(reference))
    rng=np.random.default_rng(int(peak*100));x,y,width=[series[k] for k in ['x','y','width']]
    model=prediction(x,width,central['params'],1);residual=y-model;block=max(1,int(round(1/(p['cadence_days']*24))))
    boot=[]
    for _ in range(draws):
        starts=rng.integers(0,len(x),size=int(np.ceil(len(x)/block)))
        index=np.concatenate([(np.arange(block)+s)%len(x) for s in starts])[:len(x)]
        fit=fit_event(x,model+residual[index],width,1,central['params'])
        if fit['success']:boot.append([fit[k] for k in ['amplitude_fractional','fwhm_h','equivalent_duration_s','rise_tau_h','decay_tau_h']])
    return dict(central=central,variants=variants,
        bootstrap=dict(method='One-hour circular residual blocks; conditional on selected profile/baseline, not astrophysical false-alarm calibration.',draws=draws,successful=len(boot),
            columns=['amplitude_fractional','fwhm_h','equivalent_duration_s','rise_tau_h','decay_tau_h'],p16_p50_p84=np.quantile(boot,[.16,.5,.84],axis=0).tolist()),_series=series)

def fit_image_offset(image,template):
    """Signed empirical-template shift plus a planar background; pixels."""
    image,template=np.asarray(image,float),np.asarray(template,float)
    yy,xx=np.indices(image.shape);plane=np.column_stack([np.ones(image.size),xx.ravel(),yy.ravel()])
    q=np.linalg.qr(plane)[0];data=image.ravel()-q@(q.T@image.ravel())
    def solve(offset):
        # grid-constant interpolates continuously through the boundary. 'constant'
        # drops an edge row abruptly at a tiny nonzero shift and can pin the fit.
        shifted=shift(template,[offset[1],offset[0]],order=3,mode='grid-constant',cval=0).ravel();v=shifted-q@(q.T@shifted)
        amp=max(0.,float(v@data/max(v@v,1e-30)));return amp,v
    opt=least_squares(lambda off:data-solve(off)[0]*solve(off)[1],[0.,0.],bounds=(-1.25,1.25),ftol=1e-9,xtol=1e-9,gtol=1e-9,max_nfev=100)
    amp,v=solve(opt.x);noise=np.linalg.norm(data-amp*v)
    return dict(offset_xy=opt.x.tolist(),offset_pixels=float(np.linalg.norm(opt.x)),amplitude=float(amp),localized=bool(amp*np.linalg.norm(v)>1e-7),
                at_bound=bool(np.any(abs(opt.x)>1.24)),residual_norm=float(noise))

def target_template(p,maps):
    # Keep the observed target profile; remove the outer-ring plane. It is not an official TESS PRF.
    image=maps['baseline'];xx,yy=p['xx'],p['yy'];plane=np.column_stack([np.ones(image.size),xx.ravel(),yy.ravel()])
    return image-(plane@np.linalg.lstsq(plane[p['edge'].ravel()],image[p['edge']],rcond=None)[0]).reshape(image.shape)

def track_velocity(track):
    if len(track)<3:return None
    time=np.array([r['btjd'] for r in track])*24;time-=np.mean(time)
    weights=np.array([r['amplitude'] for r in track])**2
    design=np.column_stack([np.ones(len(time)),time])*np.sqrt(weights[:,None])
    offsets=np.array([r['offset_xy'] for r in track])*np.sqrt(weights[:,None])
    beta=np.linalg.lstsq(design,offsets,rcond=None)[0]
    return dict(velocity_xy_pixels_per_h=beta[1].tolist(),speed_pixels_per_h=float(np.linalg.norm(beta[1])),samples=len(track))

def matched_indices(p,desired_times):
    """Match an observed sampling pattern without filling gaps or flagged frames."""
    t=p['t'];desired=np.asarray(desired_times)
    finite=np.flatnonzero(np.isfinite(t));times=t[finite]
    after=np.clip(np.searchsorted(times,desired),0,len(times)-1);before=np.maximum(0,after-1)
    k=np.where(abs(times[before]-desired)<abs(times[after]-desired),before,after);index=finite[k]
    if np.any(abs(t[index]-desired)>.3*p['cadence_days']) or not np.all(p['good'][index]) or len(np.unique(index))!=len(index):
        raise ValueError('Control does not retain the complete observed time pattern')
    return index

def matched_difference(p,center,peak,observed):
    indices=[matched_indices(p,center+p['t'][m]-peak) for m in observed['_masks']]
    times=[float(np.mean(p['t'][i])) for i in indices];images=[np.mean(p['cube'][i],axis=0) for i in indices]
    alpha=(times[1]-times[0])/(times[2]-times[0]);ref=(1-alpha)*images[0]+alpha*images[2]
    return dict(sideband_mean_btjd=times,_maps=dict(before=images[0],after=images[2],baseline=ref,difference=images[1]-ref),
                counts=[len(i) for i in indices])

def gradient_amplitudes(image,template):
    yy,xx=np.indices(image.shape);plane=np.column_stack([np.ones(image.size),xx.ravel(),yy.ravel()]);gy,gx=np.gradient(template)
    fixed=np.column_stack([template.ravel(),plane]);moving=np.column_stack([template.ravel(),gx.ravel(),gy.ravel(),plane])
    a=np.linalg.lstsq(fixed,image.ravel(),rcond=None)[0];b=np.linalg.lstsq(moving,image.ravel(),rcond=None)[0]
    return dict(fixed_amplitude=float(a[0]),motion_nuisance_amplitude=float(b[0]),
        amplitude_relative_change=float((b[0]-a[0])/a[0]) if abs(a[0])>1e-12 else None,
        derivative_coefficients_pixels=b[1:3].tolist())

def localization(p,peak,event,draws):
    raw=tv.pixel_difference(p['t'],p['cube'],p['quality'],event,p['xy'],guard=.15,side=.3)
    template=target_template(p,raw['_maps']);observed=fit_image_offset(raw['_maps']['difference'],template)
    injected_scale=raw['stationary_only_fractional_scale'];controls=[];rng=np.random.default_rng(int(peak*1000))
    for center in rng.permutation(np.linspace(np.nanmin(p['t'])+1,np.nanmax(p['t'])-1,600)):
        if any(abs(center-e[1])<1 for e in EVENTS):continue
        window=np.asarray(event)-np.mean(event)+center
        try:blank=matched_difference(p,center,peak,raw)
        except ValueError:continue
        base=target_template(p,blank['_maps']);fit=fit_image_offset(blank['_maps']['difference']+injected_scale*base,base)
        controls.append(dict(btjd=float(center),offset_xy=fit['offset_xy'],at_bound=fit['at_bound'],localized=fit['localized'],matched_counts=blank['counts'],_blank=blank,_template=base))
        if len(controls)>=draws:break
    if len(controls)<20:raise ValueError('Insufficient localization injection controls')
    valid=[r for r in controls if r['localized']];offsets=np.array([r['offset_xy'] for r in valid]);bias=np.median(offsets,axis=0)
    radial=np.linalg.norm(offsets-bias,axis=1);radius95=float(np.quantile(radial,.95));corrected=np.array(observed['offset_xy'])-bias
    # Difference positions in short slices: a descriptive displacement check only.
    track=[]
    for a in np.arange(peak-.18,peak+.181,1/24):
        select=p['good']&(p['t']>=a)&(p['t']<a+1/24)
        if sum(select)<3:continue
        time=float(np.mean(p['t'][select]));alpha=(time-raw['sideband_mean_btjd'][0])/(raw['sideband_mean_btjd'][2]-raw['sideband_mean_btjd'][0])
        ref=(1-alpha)*raw['_maps']['before']+alpha*raw['_maps']['after'];diff=np.mean(p['cube'][select],axis=0)-ref
        fit=fit_image_offset(diff,template)
        if fit['amplitude']>=.3*injected_scale and fit['localized'] and not fit['at_bound']:
            track.append(dict(btjd=time,**fit,motion_nuisance=gradient_amplitudes(diff,template),_frame_offsets=(p['t'][select]-peak).tolist()))
    velocity=track_velocity(track);speeds=[]
    if velocity is not None:
        for control in controls:
            blank=control['_blank'];base=control['_template'];fake=[]
            for hour in track:
                time=control['btjd']+hour['btjd']-peak
                try:indices=matched_indices(p,control['btjd']+np.array(hour['_frame_offsets']))
                except ValueError:
                    fake=[];break
                actual=float(np.mean(p['t'][indices]));alpha=(actual-blank['sideband_mean_btjd'][0])/(blank['sideband_mean_btjd'][2]-blank['sideband_mean_btjd'][0])
                ref=(1-alpha)*blank['_maps']['before']+alpha*blank['_maps']['after']
                difference=np.mean(p['cube'][indices],axis=0)-ref+hour['amplitude']*base
                fit=fit_image_offset(difference,base)
                if fit['localized'] and not fit['at_bound']:fake.append(dict(btjd=actual,**fit))
            v=track_velocity(fake)
            if v is not None:speeds.append(v['speed_pixels_per_h'])
        velocity.update(stationary_control_count=len(speeds),stationary_speed95_pixels_per_h=float(np.quantile(speeds,.95)) if speeds else None,
            empirical_stationary_reference_p=float((1+np.count_nonzero(np.array(speeds)>=velocity['speed_pixels_per_h']))/(len(speeds)+1)))
    return dict(empirical_template_fit=observed,event_motion_nuisance=gradient_amplitudes(raw['_maps']['difference'],template),
        observed_before_event_after_counts=[int(sum(m)) for m in raw['_masks']],
        calibration_bias_pixels=bias.tolist(),calibrated_offset_pixels=float(np.linalg.norm(corrected)),
        injection_radius95_pixels=radius95,within_empirical95_radius=bool(np.linalg.norm(corrected)<=radius95),
        controls=controls,control_count=len(controls),control_fit_bound_fraction=float(np.mean([r['at_bound'] for r in controls])),
        prior_positive_moment_offset_pixels=raw['difference_centroid_offset_pixels'],hourly_track=track,hourly_motion=velocity,
        motion_assessment='Amplitude-weighted hourly linear displacement compared with blank-time stationary injections at the observed hourly amplitudes. Conditional noise diagnostic, not a blanket asteroid exclusion or orbital proper motion.',
        interpretation='Offset relative to the observed target template, calibrated with stationary injections at blank times. Conditional on template and empirical noise, not an association posterior or resolved-component identification.',
        _maps=raw['_maps'],_template=template)

def sensitivity(products,events):
    rows=[]
    for sector,p in sorted(products.items()):
        t,flux=aperture_series(p,1.5);y=flux/np.median(flux)-1;scan=search.screen(t,y);b=tv.bin_series(t,y);bs=search.screen(b['t'],b['y'])
        centers=t.min()+(np.arange(48)+.5)/48*np.ptp(t)
        for event in events:
            fit=event['photometry']['central'];amp=fit['amplitude_fractional'];rise=fit['rise_tau_h'];decay=fit['decay_tau_h'];records=[]
            for peak in centers:
                shape=amp*exposure_profile((t-peak)*24,np.full(len(t),p['cadence_days']*24),0,rise,decay)
                visibility=float(np.sum(shape)*p['cadence_days']*86400/fit['equivalent_duration_s'])
                residual=tv.injected_residual(t,y,scan['baseline'],shape,2.)
                added_b=tv.bin_values(shape,b['inverse'],len(b['t']));res_b=tv.injected_residual(b['t'],b['y'],bs['baseline'],added_b,2.)
                tolerance=max(fit['fwhm_h']/24,1/24)
                detected=search.threshold_runs(t,residual,scan['noise'])+search.threshold_runs(b['t'],res_b,bs['noise'])
                preexisting=any(abs(r['peak_btjd']-peak)<tolerance for r in scan['events']+bs['events'])
                recovered=any(abs(r['peak_btjd']-peak)<tolerance for r in detected) and not preexisting
                records.append(dict(peak_btjd=float(peak),sampled_model_fraction=visibility,new_recovery=bool(recovered),preexisting_match=bool(preexisting)))
            available=[r for r in records if r['sampled_model_fraction']>=.8]
            recovery=sum(r['new_recovery'] for r in records)/len(records);span=float(np.ptp(t))
            rows.append(dict(sector=sector,event_date=event['date_utc'][:10],trials=len(records),calendar_span_days=span,
                new_recoveries=sum(r['new_recovery'] for r in records),well_sampled_trials=len(available),
                recovered_well_sampled=sum(r['new_recovery'] for r in available),effective_search_days=span*recovery,
                nominal_accepted_cadence_days=float(len(t)*p['cadence_days']),trials_detail=records))
        print('Direct-aperture sensitivity finished',sector,flush=True)
    return dict(rows=rows,profiles='Measured finite-cadence asymmetric profiles;48 uniform peak times per interval, including gaps. Native OR20min screen; observed noise held fixed, running baseline recomputed.',
        qualification='Effective search time is model-specific, finite-grid exploratory recovery exposure after FFI calibration. It is not a universal completeness correction or a population flare-frequency distribution.')

def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items() if not k.startswith('_')}
    if isinstance(value,list):return [clean(v) for v in value]
    return value

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--bootstrap',type=int,default=160);parser.add_argument('--positions',type=int,default=80);parser.add_argument('--sensitivity',action='store_true');parser.add_argument('--reuse-sensitivity',action='store_true');args=parser.parse_args()
    cached_sensitivity=None
    if args.reuse_sensitivity:
        previous=(OUT/'comparison.json').read_bytes();old=json.loads(previous)
        cached_sensitivity=old['sensitivity']
        if cached_sensitivity is None:raise ValueError('No saved sensitivity experiment to reuse')
        cached_sensitivity['reused_from_result_sha256']=hashlib.sha256(previous).hexdigest()
        cached_sensitivity['original_analysis_source_sha256']=old['source_sha256']
    products=pixel_products();events=[]
    bp=read_bandpass(ROOT/'results/sky-plates-j0754-flare-energy-20261108/TESS_Red.xml','TESS/TESS.Red')
    if bp.detector!=0:raise ValueError('Expected retained energy response')
    wave=bp.wave*1e-10;nominal_factor=energy_factor(wave,bp.transmission,5205.7764,9000);nominal_lum=54.640694*3.828e33
    fig,axes=plt.subplots(4,2,figsize=(13,15),layout='constrained')
    mapsfig,mapaxes=plt.subplots(2,2,figsize=(10,9),layout='constrained')
    for i,(sector,peak,window) in enumerate(EVENTS):
        p=products[sector];photo=measurements(p,peak,args.bootstrap);loc=localization(p,peak,window,args.positions)
        fit=photo['central'];ed=fit['equivalent_duration_s'];variant_ed=[r['equivalent_duration_s'] for r in photo['variants']]
        grid=[duration*lum*3.828e33*energy_factor(wave,bp.transmission,ts,tf) for duration in [min(variant_ed),max(variant_ed)] for lum in [30,60] for ts in [4500,5250] for tf in [8000,12000]]
        record=dict(sector=sector,selected_peak_btjd=peak,date_utc=Time(peak+2457000,format='jd',scale='tdb').utc.isot,photometry=photo,localization=loc,
            energy=dict(isotropic_equivalent_erg=ed*nominal_lum*nominal_factor,face_on_patch_erg=ed*nominal_lum*nominal_factor/4,
                isotropic_parameter_sensitivity_range_erg=[min(grid),max(grid)],face_on_parameter_sensitivity_range_erg=[min(grid)/4,max(grid)/4],
                model_ED_baseline_aperture_range_s=[min(variant_ed),max(variant_ed)],temperature_geometry_scope='Same nominal L=54.640694Lsun,Tstar=5205.7764K,Tflare=9000K for every event; sensitivity envelope is not a confidence interval.'),
            native_input=p['path'].relative_to(ROOT).as_posix())
        events.append(record);s=photo['_series'];x=s['x'];b=tv.bin_series(x/24,s['relative']);axes[i,0].plot(b['t']*24,b['y']*100,'.-',label='Direct aperture,20min')
        fine=np.linspace(-24,30,800);curve=fit['amplitude_fractional']*exposure_profile(fine,np.full_like(fine,.001),*fit['params'][1:4]);axes[i,0].plot(fine,curve*100,label='Asymmetric profile');axes[i,0].axhline(0,color='gray',lw=.7)
        axes[i,0].set(title=f"{record['date_utc'][:10]} | FWHM {fit['fwhm_h']:.2f}h | ED {ed:.0f}s",xlabel='Hours from selected peak',ylabel='Excess (%)');axes[i,0].legend(fontsize=8)
        for variant in photo['variants']:
            axes[i,1].plot(variant['aperture_pixels'],variant['equivalent_duration_s'],'o',color='tab:blue' if variant['degree']==1 else 'tab:orange',alpha=.7)
        axes[i,1].axhline(ed,color='black',lw=.8);axes[i,1].set(title='Blue:linear; orange:quadratic baseline',xlabel='Aperture radius (pixels)',ylabel='Full-model ED (seconds)')
        ax=mapaxes.ravel()[i];im=ax.imshow(loc['_maps']['plane_subtracted_difference'],origin='lower',cmap='RdBu_r');ax.plot(*p['xy'],'wx',ms=10,mew=2)
        ax.set(title=f"{record['date_utc'][:10]} | template shift {loc['calibrated_offset_pixels']:.2f}px\nInjected95% radius {loc['injection_radius95_pixels']:.2f}px",xlabel='Pixel x',ylabel='Pixel y');mapsfig.colorbar(im,ax=ax,label='e-/s')
        print(record['date_utc'][:10],clean(dict(photometry=fit,energy=record['energy'],offset=loc['calibrated_offset_pixels'],radius95=loc['injection_radius95_pixels'])),flush=True)
    fig.suptitle('J0754: common TESS event measurements | local baseline removed; signed data')
    fig.savefig(OUT/'four_event_measurements.png',dpi=150);plt.close(fig)
    mapsfig.suptitle('Signed event differences; target WCS marked | empirical-template astrometry')
    mapsfig.savefig(OUT/'four_event_localization.png',dpi=160);plt.close(mapsfig)
    result=dict(events=events,energy_nominal_conversion_per_ED_second_erg=nominal_lum*nominal_factor,
        assumptions=dict(stellar_luminosity_Lsun=54.640694,stellar_temperature_K=5205.7764,flare_temperature_K=9000,
            extinction='No additional extinction or epoch-dependent luminosity correction; systematic uncertainty retained.',geometry='E_iso=4*E_face_on_patch; neither estimates kinetic/magnetic energy.'),
        sensitivity=sensitivity(products,events) if args.sensitivity else cached_sensitivity,
        input_sha256={p['path'].relative_to(ROOT).as_posix():hashlib.sha256(p['path'].read_bytes()).hexdigest() for p in products.values()},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),protocol_sha256=hashlib.sha256((OUT/'protocol.json').read_bytes()).hexdigest())
    (OUT/'comparison.json').write_text(json.dumps(clean(result),indent=2,allow_nan=False)+'\n',encoding='utf8',newline='\n')
    with (OUT/'event_measurements.csv').open('w',newline='',encoding='utf8') as fp:
        writer=csv.DictWriter(fp,fieldnames=['date','sector','peak_percent','fwhm_h','rise_tau_h','decay_tau_h','model_ED_s','observed_signed_ED_s','isotropic_energy_erg','face_on_patch_energy_erg','template_offset_pixels','injection_radius95_pixels']);writer.writeheader()
        for r in events:
            f=r['photometry']['central'];writer.writerow(dict(date=r['date_utc'][:10],sector=r['sector'],peak_percent=100*f['amplitude_fractional'],fwhm_h=f['fwhm_h'],rise_tau_h=f['rise_tau_h'],decay_tau_h=f['decay_tau_h'],model_ED_s=f['equivalent_duration_s'],observed_signed_ED_s=f['observed_signed_ED_s'],isotropic_energy_erg=r['energy']['isotropic_equivalent_erg'],face_on_patch_energy_erg=r['energy']['face_on_patch_erg'],template_offset_pixels=r['localization']['calibrated_offset_pixels'],injection_radius95_pixels=r['localization']['injection_radius95_pixels']))

if __name__=='__main__':main()
