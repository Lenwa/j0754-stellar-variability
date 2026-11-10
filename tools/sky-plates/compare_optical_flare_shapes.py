"""Peak/FWHM-normalized visual comparison; run from the repository root."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT = Path('results/sky-plates-j0754-optical-overlay-20261108')
INPUT = Path('results/sky-plates-j0754-physical-tests-20261108/tess/sector46_pixel_event_timeseries.csv')


def normalize(t, y):
    """Use nearest bracketing half-height crossings around the observed peak."""
    k = int(np.argmax(y))
    y = y / y[k]
    left = np.flatnonzero(y[:k] < .5)[-1]
    right = k + np.flatnonzero(y[k:] < .5)[0]
    def cross(a, b):
        return float(t[a] + (.5-y[a])*(t[b]-t[a])/(y[b]-y[a]))
    lo, hi = cross(left, left+1), cross(right-1, right)
    width = hi-lo
    assert width > 0 and lo < t[k] < hi
    return (t-t[k])/width, y, width, float(t[k])


def main():
    OUT.mkdir(exist_ok=True)
    a = np.genfromtxt(INPUT, delimiter=',', names=True, dtype=None, encoding='utf8')
    t = (a['btjd']-2572.0782528500076)*24
    y = a['relative_core_flux']
    tx, ty, tw, _ = normalize(t, y)
    # Approximate manual digitization of RED optical-continuum kernel in
    # Hao et al. 2017 Fig.3c, original 660x856 JPEG. Not raw telescope samples.
    # Chosen as an impulsive solar example, not optimized against J0754.
    xp = np.array([97,118,140,162,184,207,229,250,273,297,319,341,363,385,407,429,451,473,495,513,535,557,579])
    yp = np.array([585,592,602,609,605,602,600,598,594,512,543,539,550,562,590,600,596,591,586,584,623,619,616])
    st = (xp-162)*16/(513-162)+3  # minutes after 03:13, from 03:16 and 03:32 ticks
    sy = (594-yp)/330  # panel c: y=594 is 0; y=528 is 0.2 contrast
    sx, sy, sw, _ = normalize(st, sy)
    np.savetxt(OUT/'solar_red_kernel_digitized.csv', np.c_[xp,yp,st,sy], delimiter=',',
               header='figure_x_pixel,figure_y_pixel,minutes_after_0313,peak_normalized_contrast', comments='')
    # Published empirical white-light template, Davenport 2014 Eqs.1,4.
    u = np.linspace(-3,8,11001)
    f = np.zeros_like(u)
    rise = (u > -1) & (u <= 0)
    f[rise] = 1+1.941*u[rise]-.175*u[rise]**2-2.246*u[rise]**3-1.125*u[rise]**4
    decay = u > 0
    f[decay] = .689*np.exp(-1.600*u[decay])+.303*np.exp(-.2783*u[decay])
    dx, dy, _, _ = normalize(u, f)
    fig, ax = plt.subplots(figsize=(10,5.7))
    ax.plot(tx,ty,'.-',color='#275a99',lw=1,ms=3,label=f'J0754 candidate â€” TESS optical (FWHM {tw:.1f} h)')
    ax.plot(sx,sy,'o-',color='#d55232',lw=1.5,ms=4,label=f'Solar flare, 10 Mar 2015 â€” 425 nm kernel (~{sw:.1f} min)')
    ax.plot(dx,dy,'--',color='#19806b',lw=2,label='Davenport 2014 empirical optical stellar-flare template')
    ax.axhline(0,color='.6',lw=.7); ax.axvline(0,color='.8',lw=.7)
    ax.set(xlim=(-1.6,3.5),ylim=(-.42,1.13),xlabel='Time from peak / each curveâ€™s full width at half maximum',
           ylabel='Excess brightness / peak excess', title='Does J0754 resemble an optical flare?')
    ax.legend(fontsize=9,loc='upper right');ax.grid(alpha=.15)
    fig.text(.09,.025,'Peaks and widths aligned; no shape fitting. Solar points approximated from Hao et al. (2017), Fig. 3c.\nDifferent optical passbands and spatial scales. This compares shapes, not energies or source identity.',fontsize=9)
    fig.tight_layout(rect=(0,.09,1,1));fig.savefig(OUT/'optical_flare_overlay.png',dpi=180)
    receipt = {'input':str(INPUT),'normalization':'Original TESS relative-flux baseline retained; solar plotted zero contrast retained; divide excess by peak, align peaks, divide time by nearest linearly interpolated FWHM. No independent rise/decay stretches or fit.',
        'tess_fwhm_hours':tw,'solar_fwhm_minutes_approximate':sw,
        'solar_source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5736556/',
        'solar_method':'Approximate manual digitization, red kernel Fig3c 4250 Angstrom; coordinates preserved; about 2 pixels reading precision, no formal errors.',
        'stellar_template_source':'https://arxiv.org/html/1411.3723','stellar_template':'Equations 1 and 4; empirical GJ1243 classical flare template, not a separate named FRED-superflare catalogue.',
        'limits':'Visual comparison only. TESS candidate is on/near J0754; blend/systematic alternatives remain. No superflare energy classification.'}
    (OUT/'comparison_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf8')
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
