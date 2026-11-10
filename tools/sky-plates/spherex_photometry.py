"""SPHEREx ePSF sampling and weighted point-source extraction contracts."""
import numpy as np
from astropy.io.fits import Header
from astropy.wcs import WCS
from photutils.psf import ImagePSF
from scipy.optimize import minimize

BAD_BITS=[0,1,2,6,7,9,10,11,14,15,17,19,20,22,23,24,25,26,27,28,29]
NATIVE_BAD_BITS=[2,6,9,10,11,15,17,26,27,28,29]


def render_epsf(epsf,shape,x,y,oversampling=(5,5)):
    """Sample the already pixel-convolved ePSF, then normalize native stamp.

    Matches the current IRSA tutorial's ImagePSF.evaluate implementation. Do not
    integrate the ePSF a second time or confuse its oversampled sum with the
    native source amplitude. Normalize on an enclosing grid, so detector-edge
    cutouts preserve missing wings instead of redistributing them inwards.
    """
    model=ImagePSF(np.asarray(epsf,float),oversampling=oversampling)
    yy,xx=np.indices(shape)
    p=model.evaluate(xx,yy,flux=1.,x_0=float(x),y_0=float(y))
    half_y=int(np.ceil(np.asarray(epsf).shape[0]/oversampling[0]))+2
    half_x=int(np.ceil(np.asarray(epsf).shape[1]/oversampling[1]))+2
    gy,gx=np.mgrid[int(np.floor(y))-half_y:int(np.floor(y))+half_y+1,
                   int(np.floor(x))-half_x:int(np.floor(x))+half_x+1]
    total=model.evaluate(gx,gy,flux=1.,x_0=float(x),y_0=float(y)).sum()
    if not np.isfinite(total) or total<=0:raise ValueError('Invalid sampled ePSF normalization')
    return p/total


def pixel_mask(flags,conservative=True):
    bits=BAD_BITS if conservative else NATIVE_BAD_BITS
    return (np.asarray(flags,np.int64)&sum(1<<b for b in bits))!=0


def fit_centroid_sensitivity(epsf,image,variance,bad,x,y,oversampling=(5,5)):
    """Diagnostic bounded centroid fit; the fixed-position primary is unchanged.

    Keep the primary pixel selection/background coordinates fixed while moving
    only the source template. This diagnoses small astrometric/model errors,
    not total accuracy or an alternate astrometric catalogue measurement.
    """
    def evaluate(offset):
        p=render_epsf(epsf,image.shape,x+offset[0],y+offset[1],oversampling)
        return fit_linear_source(image,variance,p,bad,x,y)
    result=minimize(lambda offset:evaluate(offset)['chi2'],[0.,0.],method='L-BFGS-B',
                    bounds=[(-.15,.15),(-.15,.15)],options={'ftol':1e-12,'gtol':1e-7,'maxiter':80})
    fit=evaluate(result.x)
    return dict(dx=float(result.x[0]),dy=float(result.x[1]),amplitude=fit['amplitude'],
                chi2=fit['chi2'],success=bool(result.success),boundary=bool(np.any(np.abs(result.x)>.149)))


def native_xy(header,x,y):
    # Astropy otherwise attaches the primary celestial SIP polynomial even to
    # the alternate LINEAR WCS. Isolate its cards before constructing the map.
    prefixes=('WCSAXES','CTYPE','CUNIT','CRPIX','CRVAL','CDELT','PC','CD','CROTA')
    linear=Header({k:v for k,v in header.items() if k.endswith('A') and k.startswith(prefixes)})
    a,b=WCS(linear,key='A').all_pix2world(float(x),float(y),0)
    return float(a),float(b)


def calibration_value(image,x,y):
    a=np.asarray(image)
    if not 0<=x<=a.shape[1]-1 or not 0<=y<=a.shape[0]-1:raise ValueError('Outside calibration image')
    x0=int(np.floor(x));y0=int(np.floor(y));x1=min(x0+1,a.shape[1]-1);y1=min(y0+1,a.shape[0]-1)
    dx=x-x0;dy=y-y0
    return float((1-dy)*((1-dx)*a[y0,x0]+dx*a[y0,x1])+dy*((1-dx)*a[y1,x0]+dx*a[y1,x1]))


def surface_to_ujy(amplitude,area_arcsec2):
    return float(amplitude*area_arcsec2*1e12/206264.80624709636**2)


def fit_linear_source(image,variance,psf,bad,x,y,radius=5.,plane=True):
    image=np.asarray(image,float);variance=np.asarray(variance,float)
    yy,xx=np.indices(image.shape);dx=xx-x;dy=yy-y
    good=np.isfinite(image)&np.isfinite(variance)&(variance>0)&(~bad)&(dx**2+dy**2<=radius**2)
    columns=[psf,np.ones_like(psf)]
    if plane:columns.extend([dx,dy])
    if good.sum()<max(12,len(columns)+3):raise ValueError('Insufficient valid source/background pixels')
    a=np.column_stack([c[good] for c in columns]);w=1/np.sqrt(variance[good])
    aw=a*w[:,None];yw=image[good]*w
    beta,_,rank,_=np.linalg.lstsq(aw,yw,rcond=None)
    if rank<len(columns):raise ValueError('Rank-deficient point-source fit')
    covariance=np.linalg.inv(aw.T@aw)
    model=sum(b*c for b,c in zip(beta,columns));resid=image-model
    core=good&(dx**2+dy**2<=2.5**2)
    chi=float(np.sum((resid[good]*w)**2));dof=max(1,int(good.sum())-len(columns))
    return dict(amplitude=float(beta[0]),amplitude_error=float(np.sqrt(covariance[0,0])),
                background=float(beta[1]),background_x=float(beta[2]) if plane else 0.,
                background_y=float(beta[3]) if plane else 0.,n_pixels=int(good.sum()),
                chi2=chi,reduced_chi2=chi/dof,
                fit_ql=float(np.mean(np.abs(resid[core])/np.sqrt(variance[core]))) if np.any(core) else None,
                retained_psf_fraction=float(psf[good].sum()),model=model,residual=resid,good=good)
