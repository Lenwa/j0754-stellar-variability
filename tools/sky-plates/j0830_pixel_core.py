"""Small numerical tools for direct-pixel variability attribution."""
import numpy as np
from scipy.ndimage import shift
from tang_sibling_core import nuisance_matrix,orthogonal_columns,signal_matrix


def background_planes(cube,positions,mask_radius=2.8):
    cube=np.asarray(cube,float);yy,xx=np.indices(cube.shape[1:]);allowed=np.ones(xx.shape,bool)
    for x,y in positions:allowed&=(xx-x)**2+(yy-y)**2>mask_radius**2
    design=np.column_stack([np.ones(xx.size),xx.ravel(),yy.ravel()]);planes=[]
    for frame in cube:
        mask=allowed.ravel()&np.isfinite(frame.ravel());values=frame.ravel()
        if sum(mask)<20:planes.append(np.full(frame.shape,np.nan));continue
        for _ in range(3):
            beta=np.linalg.lstsq(design[mask],values[mask],rcond=None)[0];residual=values-design@beta
            center=np.median(residual[mask]);scale=1.4826*np.median(abs(residual[mask]-center))
            if scale<1e-10:break
            next_mask=mask&(abs(residual-center)<=4*scale)
            if sum(next_mask)<20:break
            mask=next_mask
        planes.append((design@beta).reshape(frame.shape))
    return np.asarray(planes)


def aperture_flux(cube,xy,radius):
    yy,xx=np.indices(np.asarray(cube).shape[1:]);mask=(xx-xy[0])**2+(yy-xy[1])**2<=radius**2
    if not np.any(mask):raise ValueError('Empty aperture')
    return np.sum(np.asarray(cube)[:,mask],axis=1)


def register_on_control(cube,xy,radius=2.5):
    """Align translations using only an independent bright comparison star.

    Offsets are measured x,y relative to its median measured position. NaN
    edges from interpolation must be cropped before pixel-by-pixel fitting.
    """
    cube=np.asarray(cube,float);yy,xx=np.indices(cube.shape[1:]);mask=(xx-xy[0])**2+(yy-xy[1])**2<=radius**2
    weights=np.maximum(cube[:,mask],0);total=np.sum(weights,axis=1)
    if np.any(total<=0) or not np.all(np.isfinite(total)):raise ValueError('Comparison centroid unavailable')
    centroids=np.column_stack([weights@xx[mask]/total,weights@yy[mask]/total]);offsets=centroids-np.median(centroids,axis=0)
    aligned=np.stack([shift(frame,(-offset[1],-offset[0]),order=3,mode='constant',cval=np.nan) for frame,offset in zip(cube,offsets)])
    return aligned,offsets


def differential_flux(target,controls):
    target,controls=np.asarray(target,float),np.asarray(controls,float)
    if controls.ndim!=2 or controls.shape[1]<1:raise ValueError('At least one comparison star required')
    if not np.all(np.isfinite(controls)) or np.any(controls<=0):raise ValueError('Controls must be finite and positive on common cadences')
    ensemble=np.median(controls/np.median(controls,axis=0),axis=1)
    return target/ensemble,ensemble


def periodic_image(t,cube,groups,period,reference_beta,mode='sector_linear'):
    cube=np.asarray(cube,float);q=orthogonal_columns(nuisance_matrix(t,groups,mode));s=signal_matrix(t,period);s-=q@(q.T@s)
    values=cube.reshape(len(t),-1);values=values-q@(q.T@values)
    beta=np.linalg.lstsq(s,values,rcond=None)[0];direction=np.asarray(reference_beta,float);direction=direction/np.linalg.norm(direction)
    signed=(direction@beta).reshape(cube.shape[1:]);quadrature=(np.array([-direction[1],direction[0]])@beta).reshape(cube.shape[1:])
    return dict(signed=signed,quadrature=quadrature,cosine=beta[0].reshape(signed.shape),sine=beta[1].reshape(signed.shape),signal_condition=float(np.linalg.cond(s)))


def spatial_components(template,signed,xy,radius=3.):
    """Fit stationary light, its two motion derivatives and a background plane."""
    template,signed=np.asarray(template,float),np.asarray(signed,float);yy,xx=np.indices(template.shape);dy,dx=np.gradient(template)
    mask=((xx-xy[0])**2+(yy-xy[1])**2<=radius**2)&np.isfinite(template)&np.isfinite(signed)
    matrix=np.column_stack([template[mask],dx[mask],dy[mask],np.ones(sum(mask.ravel())),xx[mask],yy[mask]])
    beta=np.linalg.lstsq(matrix,signed[mask],rcond=None)[0];residual=signed[mask]-matrix@beta
    return dict(stationary_fraction=float(beta[0]),motion_coefficients=beta[1:3].tolist(),plane_coefficients=beta[3:].tolist(),residual_rms=float(np.sqrt(np.mean(residual**2))),condition=float(np.linalg.cond(matrix)),pixels=int(sum(mask.ravel())))
