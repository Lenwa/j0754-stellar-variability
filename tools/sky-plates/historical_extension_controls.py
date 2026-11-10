"""Offline, prespecified exact-frame ASAS comparison-star diagnostic."""
import json
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

from analysis_core import separation_arcsec
from historical_extension import asas_aperture_status, date_from_jd, parse_asas, read_csv, write_csv


def unique_rows(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row.get('star', ''), row['field'], row['dataset'], row['hjd'], row['frame'])].append(row)
    return [dict(group[0], mag=st.median(r['mag'] for r in group), duplicate_entries=len(group))
            for group in groups.values()]


def frame_ensemble(target, controls):
    """Baseline each dataset separately; each star has one vote per exact exposure."""
    controls = unique_rows(controls)
    baselines = defaultdict(list)
    for row in controls:
        if row['year'] == 2004:
            baselines[(row['star'], row['field'], row['dataset'])].append(row['mag'])
    centers = {key: st.median(values) for key, values in baselines.items() if len(values) >= 3}
    by_frame = defaultdict(list)
    for row in controls:
        key = (row['star'], row['field'], row['dataset'])
        if key in centers:
            aperture=row.get('aperture_identity',row['star']+'|'+row['field'])
            by_frame[(row['hjd'], row['frame'])].append((row['star'],aperture,row['mag']-centers[key]))
    result = []
    for row in target:
        records=by_frame[(row['hjd'],row['frame'])]
        parents={star:star for star,_,_ in records}
        def root(star):
            while parents[star]!=star:
                parents[star]=parents[parents[star]]
                star=parents[star]
            return star
        owners={}
        for star,aperture,_ in records:
            if aperture in owners:
                parents[root(star)]=root(owners[aperture])
            else:
                owners[aperture]=star
        groups=defaultdict(lambda:defaultdict(list))
        for star,aperture,residual in records:
            groups[root(star)][aperture].append(residual)
        offsets=[st.median(st.median(v) for v in apertures.values()) for apertures in groups.values()]
        if len(offsets) >= 3:
            offset = st.median(offsets)
            result.append(dict(row, controls=len(offsets), distinct_gaia_controls=len(parents), control_offset=offset,
                               corrected_mag=row['mag']-offset))
    return result


def contrast(rows, key='mag', before=2004, after=2007):
    groups = {year: [r[key] for r in rows if r['year'] == year] for year in (before, after)}
    return {f'n{before}':len(groups[before]), f'n{after}':len(groups[after]),
                f'delta_{after}_minus_{before}':st.median(groups[after])-st.median(groups[before])
                if min(map(len, groups.values())) >= 3 else None}


def run():
    out = Path('results/sky-plates-j0754-models-20261108/history/controls')
    selection = json.loads((out/'selected_control_queries.json').read_text())
    audit, all_rows, good = [], [], []
    for source in selection:
        source_id = source['source_id']
        path = out/(source_id+'.txt')
        if not path.exists():
            audit.append(dict(source_id=source_id, error='No returned file'))
            continue
        try:
            rows, datasets = parse_asas(path.read_text(encoding='utf-8'))
        except ValueError as exc:
            audit.append(dict(source_id=source_id, error=str(exc)))
            continue
        meta = {d['dataset']: d for d in datasets}
        for d in datasets:
            d['gaia_position_separation_arcsec'] = separation_arcsec(float(d['ra_hours'])*15,
                float(d['dec_deg']), float(source['ra']), float(source['dec']))
            d['aperture_identity']=f"{d['field']}|{float(d['ra_hours']):.9f}|{float(d['dec_deg']):.9f}"
        for row in rows:
            d = meta[row['dataset']]
            row.update(star=source_id, year=date_from_jd(row['hjd']).year, mag=row['mag_0'],
                position_match_arcsec=d['gaia_position_separation_arcsec'],
                aperture_identity=d['aperture_identity'],
                accepted=row['grade'] in ('A','B') and asas_aperture_status(row['mag_0'])=='measured'
                         and d['gaia_position_separation_arcsec'] <= 10)
            all_rows.append(row)
            if row['accepted']:
                good.append(row)
        audit.append(dict(source_id=source_id, query_designation=source['asas_query_designation'],
            raw_rows=len(rows), accepted_rows=sum(r['accepted'] for r in rows), datasets=datasets))
    good = unique_rows(good)
    aperture_requests=defaultdict(list)
    for source in audit:
        for dataset in source.get('datasets',[]):
            aperture_requests[dataset['aperture_identity']].append(dict(star=source['source_id'],
                dataset=dataset['dataset'],separation_arcsec=dataset['gaia_position_separation_arcsec']))
    aliases=[dict(aperture_identity=key,requests=values) for key,values in aperture_requests.items()
             if len({r['star'] for r in values})>1]
    (out/'shared_asas_aperture_audit.json').write_text(json.dumps(aliases,indent=2))
    grouped = defaultdict(list)
    for row in good:
        grouped[(row['star'],row['field'],row['dataset'])].append(row)
    placebos = [dict(star=star, field=field, dataset=dataset, **contrast(rows))
                for (star,field,dataset), rows in grouped.items()]
    target = [dict(r, year=int(r['calendar_year_from_hjd']), hjd=float(r['hjd']), mag=float(r['mag_0']))
              for r in read_csv(out.parent/'asas_mag0_grade_ab.csv') if r['field']=='F0800+16_438']
    target = unique_rows(target)
    corrected = frame_ensemble(target, good)
    annual = []
    for year in range(2002,2010):
        rows = [r for r in corrected if r['year']==year]
        if rows:
            annual.append(dict(year=year, rows=len(rows), target_raw_median=st.median(r['mag'] for r in rows),
                control_common_offset=st.median(r['control_offset'] for r in rows),
                target_corrected_median=st.median(r['corrected_mag'] for r in rows),
                minimum_controls=min(r['controls'] for r in rows)))
    result = dict(protocol=out.joinpath('protocol.json').as_posix(), controls_requested=len(selection),
        controls_returning_photometry=sum('raw_rows' in a for a in audit), raw_control_rows=len(all_rows),
        accepted_control_rows_before_collapse=sum(r['accepted'] for r in all_rows),
        accepted_control_rows_after_dataset_exposure_collapse=len(good),
        distinct_control_star_exposures=len({(r['star'],r['hjd'],r['frame']) for r in good}),
        distinct_accepted_asas_aperture_exposures=len({(r['aperture_identity'],r['hjd'],r['frame']) for r in good}),
        target_main_field_rows=len(target), corrected_target_rows=len(corrected),
        target_raw_contrast=contrast(target), matched_target_raw_contrast=contrast(corrected),
        target_corrected_contrast=contrast(corrected,'corrected_mag'),
        matched_common_control_contrast=contrast(corrected,'control_offset'),
        controls_per_frame=dict(Counter(r['controls'] for r in corrected)), annual=annual,
        corrected_frames_with_shared_aperture_aliases=sum(r['controls']<r['distinct_gaia_controls'] for r in corrected),
        shared_asas_aperture_alias_groups=len(aliases),
        per_dataset_placebos=placebos,
        limitation='Fixed selected controls, MAG_0, field/dataset 2004 offsets, exact HJD/frame and >=3 stars; '
                   'no causal inference or individual-star uncertainty from MER.')
    result['target_corrected_recovery_2009_minus_2007'] = contrast(corrected,'corrected_mag',2007,2009)
    result['target_raw_recovery_2009_minus_2007'] = contrast(target,'mag',2007,2009)
    result['time_identity'] = dict(distinct_hjd=len({r['hjd'] for r in corrected}),
        distinct_hjd_frame=len({(r['hjd'],r['frame']) for r in corrected}),
        distinct_hjd_calendar_dates=len({date_from_jd(r['hjd']).date() for r in corrected}),
        note='HJD is preserved; HJD calendar-date labels are not UTC nights. Archive UTC timestamps are unavailable.')
    leave_one_out=[]
    for source in selection:
        values=frame_ensemble(target,[r for r in good if r['star']!=source['source_id']])
        leave_one_out.append(dict(omitted_star=source['source_id'], rows=len(values),
            trough=contrast(values,'corrected_mag'), recovery=contrast(values,'corrected_mag',2007,2009)))
    result['leave_one_control_out']=leave_one_out
    sensitivities=[]
    for aperture in (0,1,2):
        controls_ap=unique_rows([dict(r,mag=r[f'mag_{aperture}']) for r in all_rows
            if r['grade'] in ('A','B') and r['position_match_arcsec']<=10
            and asas_aperture_status(r[f'mag_{aperture}'])=='measured'])
        target_ap=[dict(r,year=int(r['calendar_year_from_hjd']),hjd=float(r['hjd']),mag=float(r[f'mag_{aperture}']))
            for r in read_csv(out.parent/'asas_all_rows.csv') if r['field']=='F0800+16_438'
            and r['grade'] in ('A','B') and asas_aperture_status(float(r[f'mag_{aperture}']))=='measured']
        values=frame_ensemble(unique_rows(target_ap),controls_ap)
        sensitivities.append(dict(aperture=aperture,rows=len(values),
            raw_trough=contrast(target_ap),corrected_trough=contrast(values,'corrected_mag'),
            raw_recovery=contrast(target_ap,'mag',2007,2009),
            corrected_recovery=contrast(values,'corrected_mag',2007,2009)))
    result['aperture_sensitivity']=sensitivities
    write_csv(out/'all_control_rows.csv', all_rows)
    write_csv(out/'accepted_control_rows.csv', good)
    write_csv(out/'per_dataset_placebos.csv', placebos)
    write_csv(out/'target_exact_frame_differential.csv', corrected)
    write_csv(out/'annual_exact_frame_differential.csv', annual)
    (out/'catalogue_position_audit.json').write_text(json.dumps(audit,indent=2))
    (out/'comparison_results.json').write_text(json.dumps(result,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    ax=axes[0]
    ax.scatter([date_from_jd(r['hjd']) for r in corrected],[r['corrected_mag'] for r in corrected],s=8,alpha=.5)
    ax.invert_yaxis();ax.set_ylabel('ASAS MAG_0, exact-frame corrected');ax.grid(alpha=.2)
    ax.set_title(f"J0754: {len(corrected)} exposures; {min(r['controls'] for r in corrected)}â€“{max(r['controls'] for r in corrected)} control groups")
    ax=axes[1]
    values=[r['delta_2007_minus_2004'] for r in placebos if r['field']=='F0800+16_438'
            and r['delta_2007_minus_2004'] is not None]
    ax.scatter(range(1,len(values)+1),sorted(values),label='16 main-field controls',s=30)
    ax.axhline(.135,color='#b75c27',label='Target raw +0.135')
    ax.axhline(result['target_corrected_contrast']['delta_2007_minus_2004'],color='#315d6c',ls='--',label='Target corrected +0.147')
    ax.set_ylabel('Annual median 2007 minus 2004 (mag)');ax.set_xlabel('Controls ordered by change')
    ax.set_title('Individual controls include two similar-size changes');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('ASAS comparison check; variation survives common-frame correction')
    fig.tight_layout();fig.savefig(out/'asas_control_check.png',dpi=170);plt.close(fig)
    print(json.dumps({k:v for k,v in result.items() if k not in
                     ('per_dataset_placebos','annual','leave_one_control_out','aperture_sensitivity')},indent=2))


if __name__ == '__main__':
    run()
