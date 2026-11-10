"""Descriptive full-record comparison: retain gaps, bands and unequal amplitudes.

This is post-hoc inspection of previously analyzed data, with no optimized epoch
shift, stretch, replacement period or claim that the episodes share a cause.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time

from asassn_replication import parse_row


def main(historical,asassn,out,legacy_asas=None,legacy_controls=None):
    old=list(csv.DictReader(historical.open()))
    grouped=defaultdict(list)
    for r in old:grouped[int(float(r['year']))].append(r)
    annual=[dict(year=y,n=len(g),mag=float(np.median([float(r['corrected_mag']) for r in g])),
                 series=sorted({r['series'] for r in g})) for y,g in sorted(grouped.items())]
    optical=[]
    for raw in csv.DictReader(asassn.open()):
        r=parse_row(raw)
        if r['accepted'] and r['band']=='g':optical.append(r)
    baseline=np.median([r['mag'] for r in annual if r['n']>=3 and (1910<=r['year']<1930 or 1950<=r['year']<1970)])
    note=dict(status='DESCRIPTIVE_POST_HOC',historical_annual=annual,
        input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [historical,asassn]},
        historical_baseline=float(baseline),sparse_1892_1893_points=sum(r['n'] for r in annual if 1892<=r['year']<=1893),
        points_1973_1982=sum(r['n'] for r in annual if 1973<=r['year']<=1982),
        interpretation='Historical 1930s fade/recovery and modern multiband fade/recovery are a recorded pattern. The 1892 feature is sparse. An ~88-89-year spacing is a candidate recurrence interval, not a fitted period. A half-interval hypothesis must confront the actually observed 1970s record. No common waveform or mechanism is established.',
        caveats=['Blue-sensitive historical natural magnitudes and modern g have different responses; depths cannot be equated directly.',
                 'The historical onset is sparsely sampled, so ingress/egress symmetry is not securely measurable.',
                 'No shift, amplitude normalization or time stretching was optimized to make events match.',
                 'APPLAUSE independent historical photometric corroboration remains inconclusive.'])
    legacy=list(csv.DictReader(legacy_asas.open())) if legacy_asas else []
    if legacy_asas:
        note['input_sha256'][str(legacy_asas)]=hashlib.sha256(legacy_asas.read_bytes()).hexdigest()
        note['legacy_asas_annual']=legacy
        note['legacy_asas_qualification']='Separate V system/dataset zero points; new main-field annual shape is not yet comparison-star controlled.'
    controls=list(csv.DictReader(legacy_controls.open())) if legacy_controls else []
    if legacy_controls:
        note['input_sha256'][str(legacy_controls)]=hashlib.sha256(legacy_controls.read_bytes()).hexdigest()
        note['legacy_asas_same_frame_annual']=controls
        note['legacy_asas_qualification']='Separate V system/dataset zero points; exact HJD/frame common-control correction supports the smaller episode, without uniqueness or cause established.'
    (out/'event_comparison.json').write_text(json.dumps(note,indent=2)+'\n')
    fig,axes=plt.subplots(3 if legacy else 2,1,figsize=(12,10 if legacy else 7),layout='constrained')
    modern_ax=axes[-1]
    axes[0].scatter([float(r['year']) for r in old],[float(r['corrected_mag'])-baseline for r in old],s=8,alpha=.2,c='steelblue',label='Same-exposure differential detections')
    solid=[r for r in annual if r['n']>=3];sparse=[r for r in annual if r['n']<3]
    axes[0].plot([r['year']+.5 for r in solid],[r['mag']-baseline for r in solid],'o',color='navy',ms=4,label='Annual median: â‰¥3 detections')
    axes[0].plot([r['year']+.5 for r in sparse],[r['mag']-baseline for r in sparse],'o',mfc='none',mec='navy',ms=5,label='Sparse year: 1â€“2 detections')
    axes[0].axvspan(1930,1940,alpha=.1,color='red');axes[0].axvspan(1973,1983,alpha=.08,color='orange')
    axes[0].set_title('Historical record: 1930s event; sparse early feature; retain the 1970s coverage')
    axes[0].set_ylabel('DASCH relative magnitude\n(blue-sensitive plates)');axes[0].set_xlabel('Calendar year')
    axes[0].legend(fontsize=8,loc='upper right');axes[0].set_xlim(1890,1991)
    if legacy:
        main_field=[r for r in legacy if r['dataset']=='2']
        axes[1].plot([int(r['year'])+.5 for r in main_field],[float(r['median_mag0']) for r in main_field],'o-',color='purple')
        axes[1].set_ylabel('Legacy ASAS V magnitude\n(main field, A/B MAG_0)')
        axes[1].set_title('New 2002â€“2009 coverage: smaller annual-median trough; comparison-star check pending')
        axes[1].set_xlabel('Calendar year');axes[1].set_xlim(2001.5,2010)
        if controls:
            axes[1].plot([int(r['year'])+.5 for r in controls],[float(r['target_corrected_median']) for r in controls],'s--',color='teal',label='Exact-frame comparison corrected')
            axes[1].lines[0].set_label('Raw main-field annual median')
            axes[1].set_title('New ASAS coverage: smaller candidate fading/recovery survives exact-frame comparison')
            axes[1].legend(fontsize=8)
    for camera in sorted({r['camera'] for r in optical}):
        g=[r for r in optical if r['camera']==camera];base=[r['mag'] for r in g if r['year']==2018]
        if len({r['day'] for r in g if r['year']==2018})<10:continue
        buckets=defaultdict(list)
        for r in g:buckets[int((r['jd']-2458119.5)//30)].append(r)
        bins=[v for v in buckets.values() if len({r['day'] for r in v})>=3]
        modern_ax.plot([float(Time(np.median([r['jd'] for r in v]),format='jd').decimalyear) for v in bins],
                     [np.median([r['mag'] for r in v])-np.median(base) for v in bins],'.-',ms=3,lw=.8,label=camera)
    modern_ax.axvspan(2018,2026.86,alpha=.06,color='red');modern_ax.set_title('Modern record: camera-separated optical fading and partial recovery')
    modern_ax.set_ylabel('ASAS-SN g change\n(relative to each cameraâ€™s 2018 median)');modern_ax.set_xlabel('Calendar year')
    modern_ax.legend(ncol=5,fontsize=8);modern_ax.set_xlim(2017.5,2027)
    for ax in axes:ax.invert_yaxis();ax.grid(alpha=.2)
    for ax in [axes[0],modern_ax]:ax.axhline(0,c='gray',lw=.6)
    fig.suptitle('J0754 has repeated fading-and-recovery episodes\nDifferent bands and sampling; no forced alignment or measured long recurrence period')
    fig.savefig(out/'historical_modern_pattern.png',dpi=170);plt.close(fig)
    print('Historical points',len(old),'1973â€“1982',note['points_1973_1982'],'1892â€“1893',note['sparse_1892_1893_points'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--historical',type=Path,required=True)
    p.add_argument('--asassn',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--legacy-asas',type=Path,help='Optional per-dataset annual CSV, kept in its own V panel')
    p.add_argument('--legacy-controls',type=Path,help='Optional exact-frame-corrected legacy annual CSV')
    a=p.parse_args();main(a.historical,a.asassn,a.out,a.legacy_asas,a.legacy_controls)
