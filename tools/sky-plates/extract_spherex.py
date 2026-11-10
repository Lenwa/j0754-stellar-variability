"""Local ePSF extraction using preserved native cutouts and matched calibrations."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

from spherex_photometry import render_epsf,fit_linear_source,surface_to_ujy,calibration_value,native_xy,pixel_mask,fit_centroid_sensitivity

COORD=SkyCoord(118.69117410685*u.deg,16.694621231895*u.deg)


def load_calibration(path):
    lookup={r['uri'].split('/')[-1]:r['planeid'] for r in csv.DictReader((path/'calibration_uuid_map.csv').open())}
    result={}
    for det in range(1,7):
        e=next(path.glob(f'epsf_D{det}_*.fits'));w=next(path.glob(f'spectral_wcs_D{det}_*.fits'));s=next(path.glob(f'solid_angle_pixel_map_D{det}_*.fits'))
        with fits.open(e) as h:
            table=h['EPSF'].data.copy();over=(int(h['EPSF'].header['OVSMPY']),int(h['EPSF'].header['OVSMPX']))
        with fits.open(w) as h:wave=h['CWAVE'].data.copy();width=h['CBAND'].data.copy()
        with fits.open(s) as h:area=h['IMAGE'].data.copy()
        result[det]=dict(table=table,oversampling=over,wave=wave,width=width,area=area,
                         spectral_uuid=lookup[w.name],files=[str(e),str(w),str(s)])
    return result


def extract(path,calibration,diagnostic=False):
    with fits.open(path,memmap=False) as h:
        header=h['IMAGE'].header;image=np.asarray(h['IMAGE'].data,float)
        variance=np.asarray(h['VARIANCE'].data,float);flags=np.asarray(h['FLAGS'].data,np.int64)
        det=int(header['DETECTOR']);c=calibration[det]
        x,y=WCS(header).world_to_pixel(COORD);x=float(x);y=float(y)
        nx,ny=native_xy(header,x,y)
        history='\n'.join(str(line) for line in header.get('HISTORY',[]))
        ids=set(re.findall(r'\[CALIB\] spectral_wcs ([a-f0-9-]+)',history))
        if ids!={c['spectral_uuid']}:raise ValueError('Unmatched spectral calibration: '+str(ids))
        wave=calibration_value(c['wave'],nx,ny);width=calibration_value(c['width'],nx,ny);area=calibration_value(c['area'],nx,ny)
        table=c['table'];i=int(np.argmin((table['XCENTER']-nx)**2+(table['YCENTER']-ny)**2))
        psf=render_epsf(table['EPSF'][i],image.shape,x,y,c['oversampling'])
        yy,xx=np.indices(image.shape);r2=(xx-x)**2+(yy-y)**2
        bad=pixel_mask(flags)
        # Other catalogued sources must not set the fitted local background.
        bad|=((flags&(1<<21))!=0)&(r2>3.5**2)
        fit=fit_linear_source(image,variance,psf,bad,x,y)
        factor=surface_to_ujy(1,area)
        core=r2<=2.5**2
        combined_flags=int(np.bitwise_or.reduce(flags[core]))
        flux=fit['amplitude']*factor;error=fit['amplitude_error']*factor
        zodi_fit=fit_linear_source(image-np.asarray(h['ZODI'].data,float),variance,psf,bad,x,y)
        constant=fit_linear_source(image,variance,psf,bad,x,y,plane=False)
        native_bad=pixel_mask(flags,False)|(((flags&(1<<21))!=0)&(r2>3.5**2))
        native_fit=fit_linear_source(image,variance,psf,native_bad,x,y)
        record=dict(file=path.name,detector=det,mjd=float(header.get('MJD-OBS',header.get('MJD-BEG'))),
            season='2025' if float(header.get('MJD-OBS',header.get('MJD-BEG')))<61041 else '2026',
            x=x,y=y,native_x=nx,native_y=ny,epsf_row=i,wavelength_um=wave,width_um=width,pixel_area_arcsec2=area,
            flux_ujy=flux,error_ujy=error,background_mjy_sr=fit['background'],
            fit_ql=fit['fit_ql'],reduced_chi2=fit['reduced_chi2'],retained_psf_fraction=fit['retained_psf_fraction'],
            core_flags=combined_flags,core_masked_fraction=float(np.mean(bad[core])),
            l1_dqa=header.get('L1DQAFLG'),l2_dqa=header.get('L2DQAFLG'),fine_astrometry=header.get('FINAST'),
            image_version=h[0].header.get('VERSION'),spectral_calibration_uuid=c['spectral_uuid'],
            flux_constant_background_ujy=constant['amplitude']*factor,
            flux_zodi_subtracted_plus_plane_ujy=zodi_fit['amplitude']*factor,
            flux_native_mask_ujy=native_fit['amplitude']*factor,
            source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        for radius in [2.,3.]:
            good=(r2<=radius**2)&(~bad)&np.isfinite(image)&np.isfinite(variance)&(variance>0)
            background=fit['background']+fit['background_x']*(xx-x)+fit['background_y']*(yy-y)
            fraction=psf[good].sum()
            record[f'aperture_{int(radius)}_ujy']=float(np.sum((image-background)[good])/fraction*factor) if fraction>0 else None
        shifts=[]
        for dx in [-.3,0.,.3]:
            for dy in [-.3,0.,.3]:
                q=render_epsf(table['EPSF'][i],image.shape,x+dx,y+dy,c['oversampling'])
                f=fit_linear_source(image,variance,q,bad,x,y)
                shifts.append((f['chi2'],dx,dy,f['amplitude']*factor))
        best=min(shifts)
        record.update(shift_x=best[1],shift_y=best[2],flux_free_centroid_grid_ujy=best[3])
        fine=fit_centroid_sensitivity(table['EPSF'][i],image,variance,bad,x,y,c['oversampling'])
        record.update(fine_centroid_dx=fine['dx'],fine_centroid_dy=fine['dy'],
                      fine_centroid_flux_ujy=fine['amplitude']*factor,
                      fine_centroid_chi2=fine['chi2'],fine_centroid_success=fine['success'],
                      fine_centroid_boundary=fine['boundary'])
        record['primary_quality']=(record['l1_dqa']=='Pass' and record['l2_dqa']=='Pass' and record['fine_astrometry']==0 and
                                   record['retained_psf_fraction']>=.8 and record['core_masked_fraction']<=.2 and fit['fit_ql'] is not None and fit['fit_ql']<5)
        if diagnostic:
            fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
            for ax,a,title in zip(axes,[image,fit['model'],fit['residual']],['Calibrated cutout','ePSF + local plane','Residual']):
                im=ax.imshow(a,origin='lower');ax.plot(x,y,'+',c='cyan');ax.set_title(title);fig.colorbar(im,ax=ax,shrink=.7)
            fig.suptitle(f"{path.name}; Î»={wave:.4f} Âµm; flux={flux:.0f} Â± {error:.0f} ÂµJy; fit quality={fit['fit_ql']:.2f}")
            return record,fig
        return record,None


def run(science,calpath,out,one=None):
    out.mkdir(exist_ok=True,parents=True);cal=load_calibration(calpath)
    paths=[one] if one else sorted(science.glob('*.fits'))
    records=[];failures=[]
    for j,path in enumerate(paths):
        try:
            r,fig=extract(path,cal,diagnostic=j==0);records.append(r)
            if fig:fig.savefig(out/'extraction_diagnostic.png',dpi=150);plt.close(fig)
        except (ValueError,KeyError,np.linalg.LinAlgError) as exc:failures.append(dict(file=path.name,error=str(exc)))
    if records:
        with (out/'local_spectrophotometry.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
        fig,axes=plt.subplots(2,1,figsize=(12,8),sharex=True,layout='constrained')
        for season in ['2025','2026']:
            data=[r for r in records if r['season']==season and r['primary_quality']]
            axes[0].errorbar([r['wavelength_um'] for r in data],[r['flux_ujy']/1000 for r in data],
                            yerr=[r['error_ujy']/1000 for r in data],fmt='.',label=season,alpha=.7)
            axes[1].plot([r['wavelength_um'] for r in data],[r['aperture_3_ujy']/r['flux_ujy'] if r['flux_ujy']!=0 else np.nan for r in data],'.',label=season)
        axes[0].set_ylabel('Extracted flux density (mJy)');axes[0].legend();axes[0].set_title('Local ePSF extraction; conditional errors exclude PSF/astrometry/calibration systematics')
        axes[1].set_ylabel('Aperture / ePSF flux');axes[1].set_xlabel('Wavelength (Âµm)');axes[1].axhline(1,c='gray',lw=.6)
        for ax in axes:ax.grid(alpha=.2)
        fig.savefig(out/'local_spectrum.png',dpi=170);plt.close(fig)
    report=dict(rows=len(records),quality_rows=sum(r['primary_quality'] for r in records),failed=failures,
                records=records,calibration_files=sorted({f for c in cal.values() for f in c['files']}),
                method='Current ePSF sampled on native pixel centers, normalized full stamp; fixed Gaia position; inverse-variance source plus plane fit; conservative masks; CWAVE/CBAND and solid-angle calibrations matched by UUID.')
    (out/'local_extraction.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('Extracted',report['rows'],'quality',report['quality_rows'],'failures',len(failures))
    for r in records[:3]:print({k:r[k] for k in ['file','wavelength_um','flux_ujy','error_ujy','fit_ql','primary_quality']})
    if failures:print(failures[:5])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--science',type=Path,required=True);p.add_argument('--calibration',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--one',type=Path)
    a=p.parse_args();run(a.science,a.calibration,a.out,a.one)
