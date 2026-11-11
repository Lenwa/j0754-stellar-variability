"""Offline fixed-position radio photometry and bounded recovered-UV diagnostics."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy.stats import chi2,binomtest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

ROOT=Path(__file__).resolve().parents[2]
INPUT=ROOT/'results/sky-plates-j0754-multiwavelength-20261110'
OUT=ROOT/'results/sky-plates-j0754-radio-uv-20261110'
TARGET=SkyCoord(118.69117410685,16.694621231895,unit='deg')

def mad(y):
    y=np.asarray(y);y=y[np.isfinite(y)]
    return float(1.4826*np.median(abs(y-np.median(y))))

def tangent_offsets(world):
    icrs=world.transform_to(TARGET.frame)
    # WCS SkyCoords can retain extra equinox/obstime attributes after conversion.
    # Reconstruct ICRS coordinates so strict offset frame checks compare like frames.
    clean=SkyCoord(icrs.ra,icrs.dec,frame='icrs')
    east,north=TARGET.spherical_offsets_to(clean)
    return east.arcsec,north.arcsec

def beam(east,north,major,minor,pa):
    angle=np.deg2rad(pa)
    along=np.asarray(east)*np.sin(angle)+np.asarray(north)*np.cos(angle)
    across=np.asarray(east)*np.cos(angle)-np.asarray(north)*np.sin(angle)
    return np.exp(-4*np.log(2)*((along/major)**2+(across/minor)**2))

def forced_fit(east,north,image,cx,cy,major,minor,pa):
    dx=np.asarray(east)-cx;dy=np.asarray(north)-cy
    use=(dx*dx+dy*dy<=8**2)&np.isfinite(image)
    design=np.column_stack([beam(dx[use],dy[use],major,minor,pa),np.ones(sum(use.ravel())),dx[use],dy[use]])
    return float(np.linalg.lstsq(design,np.asarray(image)[use],rcond=None)[0][0])

def poisson_deviance(counts,exposure):
    counts=np.asarray(counts,float);exposure=np.asarray(exposure,float)
    expected=exposure*sum(counts)/sum(exposure);use=counts>0
    return float(2*np.sum(counts[use]*np.log(counts[use]/expected[use])))

def radio():
    rows=[];fig,axes=plt.subplots(1,3,figsize=(14,5),layout='constrained');files=sorted((INPUT/'native_files').glob('VLASS*.fits'))
    for path,ax in zip(files,axes):
        with fits.open(path) as hdus:
            h=hdus[0].header;data=hdus[0].data.squeeze().astype(float)
            if h['BUNIT'].strip().lower()!='jy/beam':raise ValueError('Expected Jy/beam')
            wcs=WCS(h).celestial;yy,xx=np.indices(data.shape);world=wcs.pixel_to_world(xx,yy)
            east,north=tangent_offsets(world);r=np.hypot(east,north)
            major,minor,pa=h['BMAJ']*3600,h['BMIN']*3600,h['BPA']
            amplitude=forced_fit(east,north,data,0,0,major,minor,pa)
            blank=[]
            for bx in np.arange(-24.,24.01,3.):
                for by in np.arange(-24.,24.01,3.):
                    if 10<=np.hypot(bx,by)<=24:
                        blank.append(forced_fit(east,north,data,bx,by,major,minor,pa))
            blank=np.asarray(blank);offset=float(np.median(blank));sigma=mad(blank);flux=amplitude-offset
            annulus=(r>=10)&(r<=30);pixel_rms=mad(data[annulus]);background=float(np.median(data[annulus]))
            local=(r<5);indices=np.flatnonzero(local);peak=indices[np.argmax(data[local])];iy,ix=np.unravel_index(peak,data.shape)
            near=[]
            for dx in [-1.,0.,1.]:
                for dy in [-1.,0.,1.]:near.append(forced_fit(east,north,data,dx,dy,major,minor,pa)-offset)
            record=dict(file=path.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),date=h['DATE-OBS'],frequency_Hz=h['CRVAL3'],
                beam_major_minor_pa_arcsec_deg=[major,minor,pa],flux_mJy=flux*1000,empirical_sigma_mJy=sigma*1000,
                snr=flux/sigma,pixel_rms_mJy_per_beam=pixel_rms*1000,blank_fits=len(blank),blank_offset_mJy=offset*1000,
                blank_std_mJy=float(np.std(blank,ddof=1))*1000,source_detected_5sigma=bool(flux>5*sigma),
                conservative_3sigma_image_bound_mJy=(max(0,flux)+3*sigma)*1000,
                nominal_5sigma_image_sensitivity_mJy=5*sigma*1000,
                max_pixel_within5arcsec_mJy_per_beam=float((data[iy,ix]-background)*1000),
                max_pixel_offset_arcsec=float(r[iy,ix]),position_shift_1arcsec_flux_range_mJy=[min(near)*1000,max(near)*1000])
            rows.append(record)
            show=ax.imshow(data*1000,origin='lower',cmap='RdBu_r',vmin=-.5,vmax=.5)
            xy=wcs.world_to_pixel(TARGET);ax.plot(*xy,'k+',ms=14,mew=2);ax.set_xlim(xy[0]-16,xy[0]+16);ax.set_ylim(xy[1]-16,xy[1]+16)
            ax.add_patch(Ellipse((xy[0]-11,xy[1]-11),minor,major,angle=pa,facecolor='none',edgecolor='black'))
            ax.set(title=f"{h['DATE-OBS'][:10]} | 3 GHz\n{flux*1000:+.3f} ± {sigma*1000:.3f} mJy",xlabel='Pixel x (1 arcsec/pixel)',ylabel='Pixel y')
    fig.colorbar(show,ax=axes,label='mJy/beam');fig.suptitle('VLASS: fixed J0754 position marked; common image color scale')
    fig.savefig(OUT/'radio_three_epochs.png',dpi=150);plt.close(fig)
    flux=np.array([r['flux_mJy'] for r in rows]);error=np.array([r['empirical_sigma_mJy'] for r in rows]);weight=1/error**2
    mean=float(sum(flux*weight)/sum(weight));stat=float(sum(((flux-mean)/error)**2))
    result=dict(epochs=rows,constant_flux_mJy=mean,constant_flux_sigma_mJy=float(1/np.sqrt(sum(weight))),
        constant_flux_chi2=stat,degrees_of_freedom=2,constant_flux_p=float(chi2.sf(stat,2)),
        zero_flux_chi2=float(sum((flux/error)**2)),zero_flux_p=float(chi2.sf(sum((flux/error)**2),3)),
        caveat='Empirical blank-fit noise includes correlated image pixels. Blanks overlap; scatter estimation is finite-field. QL calibration/flux biases and astrometric errors are not fully modelled. Limits are approximate image-level bounds, not precise physical-flux exclusions.')
    return result

def uvot():
    records=[];fig,axes=plt.subplots(2,2,figsize=(11,10),layout='constrained')
    with fits.open(INPUT/'native_files/sw00031554001uw1_sk.img.gz') as images,fits.open(INPUT/'native_files/sw00031554001uw1_ex.img.gz') as exposures:
        for i,(image,exposure) in enumerate(zip(images[1:],exposures[1:])):
            h=image.header;data=image.data.astype(float);exp=exposure.data.astype(float);w=WCS(h).celestial;xy=np.array(w.world_to_pixel(TARGET))
            cx,cy=np.rint(xy).astype(int);sl=(slice(cy-35,cy+36),slice(cx-35,cx+36));d=data[sl];e=exp[sl];yy,xx=np.indices(d.shape);xx+=cx-35;yy+=cy-35
            world=w.pixel_to_world(xx,yy);rad=world.separation(TARGET).arcsec;ap=rad<=5;bg=(rad>=15)&(rad<=25)&(e>0)
            background_rate=float(np.median(d[bg]/e[bg]));counts=float(np.sum(d[ap]));bgcounts=float(np.sum(e[ap])*background_rate)
            effective=float(np.median(e[ap]));lo,hi=np.quantile(e[ap],[0,1]);uniform=bool(lo>.95*hi and hi>0)
            record=dict(extension=image.name,date=h['DATE-OBS'],header_field_exposure_s=h['EXPOSURE'],target_center_exposure_s=float(exp[cy,cx]),
                aperture_exposure_min_median_max_s=[float(lo),effective,float(hi)],uniform_exposure_aperture=uniform,
                aperture_recorded_counts=counts,estimated_background_recorded_counts=bgcounts,
                uncorrected_rate_per_s=(counts-bgcounts)/effective if effective>0 and uniform else None,
                use_for_epoch_comparison=uniform,
                note='Recorded sky-image counts; resampled pixels and approximate background. No coincidence-loss/aperture/large-scale-sensitivity correction. Nonuniform exposure cannot be summarized by a single aperture exposure without PSF modelling.')
            records.append(record)
            for ax,array,title in [(axes[i,0],d,'Recorded counts'),(axes[i,1],e,'Exposure map (seconds)')]:
                im=ax.imshow(array,origin='lower',cmap='viridis',norm=matplotlib.colors.AsinhNorm(linear_width=max(1,np.nanmax(array)/200)))
                ax.plot(xy[0]-(cx-35),xy[1]-(cy-35),'r+',ms=12,mew=2);ax.set(title=f'{h["DATE-OBS"]}\n{title}',xlabel='Cutout x',ylabel='Cutout y');fig.colorbar(im,ax=ax)
    fig.suptitle('Swift UVW1: target cross and actual per-pixel coverage')
    fig.savefig(OUT/'uvot_exposure_diagnostic.png',dpi=130);plt.close(fig)
    return records

def galex():
    photons=json.loads((INPUT/'gphoton_nuv_aperture_photons.json').read_text())['data']['Tables'][0]
    fields=[f['name'] for f in photons['Fields']];p=np.asarray(photons['Rows'],float);times=p[:,fields.index('time')]/1000
    aspect=json.loads((INPUT/'gphoton_nuv_aspect.json').read_text())['data']['Tables'][0]['Rows'];at=np.array([float(r[0])/1000 for r in aspect]);cuts=np.r_[0,np.flatnonzero(np.diff(at)>2)+1,len(at)]
    visits=[];allcounts=[];allexp=[];fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for j,(a,b) in enumerate(zip(cuts[:-1],cuts[1:])):
        start=at[a]-.5;stop=at[b-1]+.5;inside=(times>=start)&(times<stop);n=int(sum(inside));seconds=float(stop-start)
        edges=np.r_[np.arange(start,stop,20.),stop];edges=np.unique(edges);counts,_=np.histogram(times[inside],edges);exposure=np.diff(edges)
        allcounts.extend(counts.tolist());allexp.extend(exposure.tolist());axes[j].errorbar((edges[:-1]+edges[1:])/2-start,counts/exposure,yerr=np.sqrt(counts)/exposure,fmt='o');axes[j].axhline(n/seconds,color='gray');axes[j].set(title=f'Visit{j+1}: {n} photons / {seconds:.0f}s',xlabel='Seconds into raw aspect run',ylabel='Uncalibrated aperture count rate / s')
        records=[]
        for path in sorted((INPUT/'native_files').glob('*scst*')):
            with fits.open(path) as hdus:
                d=hdus[1].data
                if j==int('0002' in path.name):
                    for col in ['t_dead_nuv','status_flag','hvnom_nuv','owmoving']:
                        vals=np.asarray(d[col],float);records.append(dict(column=col,unit=hdus[1].columns[col].unit,minimum=float(np.min(vals)),median=float(np.median(vals)),maximum=float(np.max(vals))))
        visits.append(dict(start_galex_s=start,stop_galex_s=stop,nominal_aspect_duration_s=seconds,aperture_photons=n,raw_rate_per_s=n/seconds,
            photon_detector_xy_mean=np.mean(p[inside][:,[fields.index('x'),fields.index('y')]],axis=0).tolist(),spacecraft_diagnostics=records,
            bins=[dict(seconds_from_start=float((aa+bb)/2-start),counts=int(c),duration_s=float(bb-aa)) for aa,bb,c in zip(edges[:-1],edges[1:],counts)]))
    fig.suptitle('GALEX NUV, 2 February2006 | 6-arcsec aperture | no response/deadtime correction')
    fig.savefig(OUT/'galex_raw_counts.png',dpi=150);plt.close(fig)
    n=[r['aperture_photons'] for r in visits];ex=[r['nominal_aspect_duration_s'] for r in visits]
    return dict(visits=visits,between_visit_exact_binomial_p=float(binomtest(n[0],sum(n),ex[0]/sum(ex)).pvalue),
        binned_constant_rate_deviance=poisson_deviance(allcounts,allexp),binned_degrees_of_freedom=len(allcounts)-1,
        binned_chi2_approx_p=float(chi2.sf(poisson_deviance(allcounts,allexp),len(allcounts)-1)),
        caveat='Raw source+background photon arrival diagnostics. Assumes constant throughput for statistical reference; deadtime, flat response, aperture losses, spacecraft flags and background not calibrated. This is not a physical UV-variability significance.')

def xray():
    t=Table.read(INPUT/'erosita_DR2_Main.ecsv');star_sep=TARGET.separation(SkyCoord(t['ra'],t['dec'],unit='deg')).arcsec[0]
    uv=Table.read(INPUT/'swuvotssc_sources_60arcsec.ecsv');center=SkyCoord(t['ra'][0],t['dec'][0],unit='deg');distance=center.separation(SkyCoord(uv['ra'],uv['dec'],unit='deg')).arcsec
    return dict(erosita_name=str(t['iauname'][0]),offset_arcsec=float(star_sep),catalogue_1sigma_pos_err_arcsec=float(t['pos_err'][0]),
        offset_in_catalogue_pos_err_units=float(star_sep/t['pos_err'][0]),
        catalogued_uvot_objects=[dict(name=str(r['name']),offset_from_xray_arcsec=float(sep),uvw1_vega_mag=float(r['uvw1_mag']),uvw1_quality_flags=int(r['uvw1_quality_flags'])) for r,sep in zip(uv,distance)],
        flags={c:str(t[c][0]) for c in t.colnames if c.startswith('flag_')},
        assessment='Consistent candidate, not uniquely identified: a faint UVOT catalogue object is also nearby. No posterior association probability or Xray flare timing inferred.')

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    result=dict(radio=radio(),uvot=uvot(),galex=galex(),erosita=xray(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),protocol_sha256=hashlib.sha256((OUT/'protocol.json').read_bytes()).hexdigest())
    inputs=list((INPUT/'native_files').glob('*'))+[INPUT/name for name in ['gphoton_nuv_aperture_photons.json','gphoton_nuv_aspect.json','erosita_DR2_Main.ecsv','swuvotssc_sources_60arcsec.ecsv']]
    result['input_sha256']={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs if p.is_file()}
    (OUT/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf8')
    with (OUT/'radio_measurements.csv').open('w',newline='',encoding='utf8') as f:
        cols=['date','flux_mJy','empirical_sigma_mJy','snr','pixel_rms_mJy_per_beam','conservative_3sigma_image_bound_mJy','source_detected_5sigma'];writer=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');writer.writeheader();writer.writerows(result['radio']['epochs'])
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
