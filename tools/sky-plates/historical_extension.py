"""Offline historical-coverage audit; preserves survey systems and rejection reasons."""
import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from analysis_core import number, quality_dasch, separation_arcsec

TARGET_RA = 118.69117410685
TARGET_DEC = 16.694621231895
TARGET_GAIA = '667401913089980544'


def parse_asas(text):
    rows, datasets = [], []
    current, pending_ndata = None, None
    for lineno, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            if line.startswith('#ndata='):
                pending_ndata = line.split('=', 1)[1].strip()
            elif line.startswith('#dataset='):
                label, rest = line.split('=', 1)[1].split(';', 1)
                fields = rest.split()
                current = dict(dataset=label.strip(), entry=fields[0], field=fields[-1],
                               ndata=pending_ndata)
                datasets.append(current)
            elif current is not None and '=' in line:
                key, value = line[1:].split('=', 1)
                if key.strip() in ('ra', 'cra'):
                    current[key.strip() + '_hours'] = value.split()[0]
                elif key.strip() in ('dec', 'cdec'):
                    current[key.strip() + '_deg'] = value.split()[0]
                else:
                    current[key.strip()] = value.strip()
            continue
        values = line.split()
        if current is None or len(values) != 13 or values[11] not in ('A', 'B', 'C', 'D'):
            raise ValueError(f'Unexpected ASAS row at line {lineno}: {line}')
        row = dict(dataset=current['dataset'], field=current['field'], entry=current['entry'],
                   source_line=lineno, hjd_minus_2450000=float(values[0]),
                   hjd=2450000 + float(values[0]), grade=values[11], frame=values[12])
        for aperture in range(5):
            row[f'mag_{aperture}'] = float(values[1 + aperture])
            row[f'frame_dispersion_{aperture}'] = float(values[6 + aperture])
        rows.append(row)
    return rows, datasets


def asas_aperture_status(value):
    if not math.isfinite(value):
        return 'invalid'
    if value == 29.999:
        return 'below_threshold_no_numeric_limit'
    if value == 99.999:
        return 'negative_flux'
    return 'measured' if 0 < value < 29.999 else 'invalid'


def dasch_reasons(row):
    reasons = []
    mag, error, limit = (number(row.get(k)) for k in
                         ('magcal_magdep', 'magcal_local_rms', 'limiting_mag_local'))
    if not math.isfinite(mag):
        reasons.append('no_photometric_detection')
    if not math.isfinite(error):
        reasons.append('missing_local_rms')
    elif error <= 0:
        reasons.append('nonpositive_local_rms')
    elif error > .3:
        reasons.append('error_above_0.3')
    if not math.isfinite(limit):
        reasons.append('missing_local_limit')
    elif math.isfinite(mag) and limit - mag < .5:
        reasons.append('within_0.5mag_of_limit')
    positions = [number(row.get(k)) for k in
                 ('ra_deg', 'dec_deg', 'ra_cat_corrected', 'dec_cat_corrected')]
    if not all(math.isfinite(v) for v in positions):
        reasons.append('missing_astrometry')
    elif separation_arcsec(*positions) > 10:
        reasons.append('position_offset_above_10arcsec')
    flags = int(row.get('aflags') or 0)
    for bit in (6, 11, 12, 13, 15):
        if flags & (1 << bit):
            reasons.append(f'aflag_bit_{bit}')
    return reasons


def applause_reasons(row):
    """The unchanged pilot APPLAUSE screen; calibration outliers are not vetoed."""
    reasons = []
    for key in ('natmag', 'cat_natmag'):
        if not math.isfinite(number(row.get(key))):
            reasons.append('missing_' + key)
    if number(row.get('phot_range_flags')) != 0:
        reasons.append('outside_calibrated_depth_range')
    if not number(row.get('model_prediction')) >= .9:
        reasons.append('artifact_probability_gate')
    if not 0 < number(row.get('natmag_error')) <= .3:
        reasons.append('magnitude_error_gate')
    if not number(row.get('annular_bin')) <= 8:
        reasons.append('plate_rim')
    return reasons


def votable_rows(path):
    from astropy.io.votable import parse_single_table
    table = parse_single_table(path).to_table()
    return [{key: '' if getattr(row[key], 'mask', False) else str(row[key])
             for key in table.colnames} for row in table]


def spherical_polygon_contains(stc_polygon, ra_deg, dec_deg):
    """Inside a convex, less-than-hemisphere spherical plate polygon, including edges."""
    values = stc_polygon.split()
    if len(values) < 8 or values[:2] != ['Polygon', 'ICRS'] or len(values[2:]) % 2:
        raise ValueError('Expected a convex STC Polygon ICRS footprint')
    coords = list(map(float, values[2:]))

    def vector(ra, dec):
        ra, dec = math.radians(ra), math.radians(dec)
        return (math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec))

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    vertices = [vector(ra, dec) for ra, dec in zip(coords[::2], coords[1::2])]
    center = tuple(sum(v[axis] for v in vertices) for axis in range(3))
    target = vector(ra_deg, dec_deg)
    for a, b in zip(vertices, vertices[1:] + vertices[:1]):
        normal = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                  a[0] * b[1] - a[1] * b[0])
        if dot(normal, target) * dot(normal, center) < -1e-14:
            return False
    return True


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def date_from_jd(jd):
    return dt.datetime(2000, 1, 1, 12) + dt.timedelta(days=jd - 2451545)


def span(rows, key):
    if not rows:
        return None
    low, high = min(float(r[key]) for r in rows), max(float(r[key]) for r in rows)
    return dict(first=low, last=high, first_calendar=date_from_jd(low).isoformat(),
                last_calendar=date_from_jd(high).isoformat())


def run(out, recon, pilot):
    raw_path = out / 'asas_raw.txt'
    asas, datasets = parse_asas(raw_path.read_text(encoding='utf-8'))
    grouped, annual, identities = defaultdict(list), defaultdict(list), Counter()
    for row in asas:
        for aperture in range(5):
            row[f'mag_{aperture}_status'] = asas_aperture_status(row[f'mag_{aperture}'])
        row['accepted_mag_0'] = row['grade'] in ('A', 'B') and row['mag_0_status'] == 'measured'
        row['calendar_year_from_hjd'] = date_from_jd(row['hjd']).year
        identity = (row['hjd'], row['frame'])
        identities[identity] += 1
        grouped[row['dataset']].append(row)
        if row['accepted_mag_0']:
            annual[(row['dataset'], row['calendar_year_from_hjd'])].append(row['mag_0'])
    for row in asas:
        row['same_hjd_frame_count_all_datasets'] = identities[(row['hjd'], row['frame'])]
    summaries = []
    for dataset in datasets:
        rows = grouped[dataset['dataset']]
        if len(rows) != int(dataset['ndata']):
            raise ValueError(f'ASAS declared count differs for dataset {dataset["dataset"]}')
        ra = float(dataset['ra_hours']) * 15
        dec = float(dataset['dec_deg'])
        accepted = [r for r in rows if r['accepted_mag_0']]
        summaries.append(dict(dataset=dataset, rows=len(rows),
                              separation_arcsec=separation_arcsec(ra, dec, TARGET_RA, TARGET_DEC),
                              grades=dict(Counter(r['grade'] for r in rows)),
                              accepted_mag_0=len(accepted), all_span=span(rows, 'hjd'),
                              accepted_span=span(accepted, 'hjd'),
                              accepted_median_mag0=statistics.median(r['mag_0'] for r in accepted)
                              if accepted else None))
    write_csv(out / 'asas_all_rows.csv', asas)
    write_csv(out / 'asas_mag0_grade_ab.csv', [r for r in asas if r['accepted_mag_0']])
    write_csv(out / 'asas_annual_by_dataset.csv',
              [dict(dataset=d, year=y, count=len(v), median_mag0=statistics.median(v))
               for (d, y), v in sorted(annual.items())])
    dasch = read_csv(recon / 'J0754_apass_lightcurve.csv')
    differential_keys = {r['key'] for r in read_csv(pilot / 'dasch_differential.csv')}
    decades, rejected = defaultdict(list), []
    for row in dasch:
        reasons = dasch_reasons(row)
        if (not reasons) != quality_dasch(row):
            raise AssertionError('Reason audit disagrees with the unchanged primary quality gate')
        year = date_from_jd(float(row['date_jd'])).year
        key = ':'.join(str(row.get(k) or '0') for k in
                       ('series', 'plate_number', 'mosaic_number', 'solution_number'))
        audited = dict(row, calendar_year=year, primary_quality_pass=not reasons,
                       rejection_reasons=';'.join(reasons), key=key,
                       differential_retained=key in differential_keys,
                       differential_reason='primary_quality_fail' if reasons else
                       ('retained' if key in differential_keys else 'fewer_than_5_eligible_same_exposure_controls'))
        rejected.append(audited)
        decades[year // 10 * 10].append(audited)
    write_csv(out / 'dasch_all_rows_quality_audit.csv', rejected)
    write_csv(out / 'dasch_before1900_quality_audit.csv', [r for r in rejected if r['calendar_year'] < 1900])
    decade_summary = []
    for decade, rows in sorted(decades.items()):
        decade_summary.append(dict(decade=decade, exposures=len(rows),
                                   detections=sum(bool(r['magcal_magdep']) for r in rows),
                                   primary_quality=sum(r['primary_quality_pass'] for r in rows),
                                   differential=sum(r['differential_retained'] for r in rows),
                                   reasons=dict(Counter(reason for r in rows for reason in
                                                        r['rejection_reasons'].split(';') if reason))))
    fresh_applause = out / 'applause_fresh_source_rows.csv'
    fresh_dates = out / 'applause_fresh_exposures.xml'
    applause = read_csv(fresh_applause if fresh_applause.exists() else recon / 'J0754_applause_cone.csv')
    dates = votable_rows(fresh_dates) if fresh_dates.exists() else read_csv(recon / 'applause_exposure_dates.csv')
    write_csv(out / 'applause_fresh_exposures.csv', dates)
    date_map = defaultdict(list)
    for row in dates:
        date_map[row['plate_id']].append(row)
    plate_map = defaultdict(list)
    for row in applause:
        row['target_separation_arcsec'] = separation_arcsec(number(row['ra_icrs']), number(row['dec_icrs']),
                                                           TARGET_RA, TARGET_DEC)
        row['exact_target_gaia_match'] = row.get('gaiaedr3_id') == TARGET_GAIA
        row['primary_quality_reasons'] = ';'.join(applause_reasons(row))
        row['primary_quality_pass'] = not row['primary_quality_reasons']
        row['ut_mid'] = ';'.join(r['ut_mid'] for r in date_map[row['plate_id']])
        plate_map[row['plate_id']].append(row)
    plate_inventory = []
    for plate_id, rows in plate_map.items():
        target = [r for r in rows if r['exact_target_gaia_match'] or r['target_separation_arcsec'] <= 3]
        natmag = [r for r in target if math.isfinite(number(r.get('natmag')))]
        plate_inventory.append(dict(plate_id=plate_id, cone_rows=len(rows),
                                    target_candidate_rows=len(target), natmag_rows=len(natmag),
                                    calibrated_non_depth_flag_rows=sum(int(r.get('phot_range_flags') or 0) == 0
                                                                      for r in natmag),
                                    minimum_separation_arcsec=min(r['target_separation_arcsec'] for r in rows),
                                    target_gaia_ids=';'.join(sorted(set(r.get('gaiaedr3_id', '') for r in target))),
                                    ut_mid=';'.join(r['ut_mid'] for r in date_map[plate_id]),
                                    depth_flags=';'.join(r.get('phot_range_flags', '') for r in natmag),
                                    calibration_flags=';'.join(r.get('phot_calib_flags', '') for r in natmag)))
    write_csv(out / 'applause_all_candidate_plates.csv', plate_inventory)
    exact_target = [r for r in applause if r['exact_target_gaia_match']]
    good_target = [r for r in exact_target if r['primary_quality_pass']]
    write_csv(out / 'applause_target_rows_quality.csv', exact_target)
    write_csv(out / 'applause_good_natural_photometry.csv', good_target)
    summary = dict(asas=dict(raw_rows=len(asas), datasets=summaries,
                             accepted_mag0_grade_ab=sum(r['accepted_mag_0'] for r in asas),
                             accepted_distinct_hjd_frame=len({(r['hjd'], r['frame']) for r in asas
                                                             if r['accepted_mag_0']}),
                             same_hjd_frame_duplicate_groups=sum(v > 1 for v in identities.values()),
                             notes=['HJD is preserved; calendar dates are labels, with no HJD/BJD conversion.',
                                    'MER columns are frame dispersion, not individual stellar uncertainties.',
                                    'MAG_0 follows the provider rule of thumb for V>=12; all five apertures retained.',
                                    'Datasets remain separate; no cross-survey or cross-field zero-point alignment.',
                                    'Grade D retained in all-rows file; only A/B measured MAG_0 used in summary.']),
                   dasch=dict(decades=decade_summary, unchanged_primary_gate_verified=True,
                              additional_rows_promoted=0),
                   applause=dict(cone_rows=len(applause), physical_plates=len(plate_map),
                                 natmag_physical_plates=sum(r['natmag_rows'] > 0 for r in plate_inventory),
                                 exact_target_natural_mag_rows=len(exact_target),
                                 primary_quality_rows=len(good_target),
                                 primary_quality_plates=len({r['plate_id'] for r in good_target}),
                                 inventory=plate_inventory,
                                 note='Refreshed full-epoch cone inventory; X/Y scans of a physical plate are correlated. '
                                      'No new pixel photometry inferred from catalogue or image availability.'))
    (out / 'historical_extension_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    candidates_path = out / 'applause_solution_candidates.xml'
    if candidates_path.exists():
        candidates = votable_rows(candidates_path)
        covered = [r for r in candidates if spherical_polygon_contains(r['stc_polygon'], TARGET_RA, TARGET_DEC)]
        write_csv(out / 'applause_covering_solution_inventory.csv', covered)
        summary['applause']['footprint_search'] = dict(center_search_radius_deg=25,
            candidate_solution_rows=len(candidates), covering_solution_rows=len(covered),
            covering_physical_plates=len({r['plate_id'] for r in covered}),
            maximum_candidate_fov1=max(number(r['fov1']) for r in candidates),
            maximum_candidate_fov2=max(number(r['fov2']) for r in candidates),
            no_catalogue_target_detection_plates=sorted({r['plate_id'] for r in covered} -
                                                        {r['plate_id'] for r in exact_target}, key=int),
            limitation='Bounded to solution centers within 25 degrees; geometric footprint is coverage, '
                       'not a detection or validated photometric limit.')
        meta_paths = [out / ('applause_covering_' + kind + '.xml')
                      for kind in ('plates', 'exposures', 'calibration')]
        if all(path.exists() for path in meta_paths):
            plates_meta, exposure_meta, calibration_meta = [votable_rows(path) for path in meta_paths]
            for kind, values in zip(('plates', 'exposures', 'calibration'),
                                    (plates_meta, exposure_meta, calibration_meta)):
                write_csv(out / ('applause_covering_' + kind + '.csv'), values)
            inventory = []
            for plate in plates_meta:
                pid = plate['plate_id']
                exps = [r for r in exposure_meta if r['plate_id'] == pid]
                calibrations = [r for r in calibration_meta if r['plate_id'] == pid]
                faint_limits = [number(r['faint_limit']) for r in calibrations
                                if math.isfinite(number(r['faint_limit']))]
                colors = [number(r['color_term']) for r in calibrations
                          if math.isfinite(number(r['color_term']))]
                inventory.append(dict(plate_id=pid, archive_id=plate['archive_id'],
                    ut_mid=';'.join(r['ut_mid'] for r in exps),
                    exposure_seconds=';'.join(r['exptime'] for r in exps),
                    telescope=plate['telescope'], emulsion=plate['emulsion'], filter=plate['filter'],
                    catalogue_target_rows=sum(r['plate_id'] == pid for r in exact_target),
                    primary_quality_target_rows=sum(r['plate_id'] == pid for r in good_target),
                    calibration_rows=len(calibrations),
                    calibration_faint_limit_min=min(faint_limits) if faint_limits else '',
                    calibration_faint_limit_max=max(faint_limits) if faint_limits else '',
                    color_term_min=min(colors) if colors else '', color_term_max=max(colors) if colors else '',
                    evidence_type='target_catalogue_measurement' if any(r['plate_id'] == pid for r in exact_target)
                                  else 'geometric_coverage_without_exact_target_catalogue_measurement'))
            write_csv(out / 'applause_covering_plate_metadata.csv', inventory)
            all_dates = [r['ut_mid'] for r in exposure_meta if r['ut_mid'].strip()]
            summary['applause']['footprint_search'].update(
                exposure_metadata_rows=len(exposure_meta), calibration_metadata_rows=len(calibration_meta),
                first_exposure=min(all_dates), last_exposure=max(all_dates),
                exposure_decades=dict(Counter(int(number(r['year_mid'])) // 10 * 10 for r in exposure_meta
                                              if math.isfinite(number(r['year_mid'])))),
                archives=dict(Counter(r['archive_id'] for r in inventory)))
            inputs_metadata = meta_paths
        else:
            inputs_metadata = []
        (out / 'historical_extension_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    inputs = [raw_path, recon / 'J0754_apass_lightcurve.csv', recon / 'J0754_applause_cone.csv',
              recon / 'applause_exposure_dates.csv', pilot / 'dasch_differential.csv']
    inputs.extend(p for p in (fresh_applause, fresh_dates, candidates_path) if p.exists())
    if candidates_path.exists():
        inputs.extend(inputs_metadata)
    (out / 'analysis_input_sha256.json').write_text(json.dumps(
        {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}, indent=2), encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax, dataset in zip(axes.flat, datasets):
        rows = grouped[dataset['dataset']]
        good = [r for r in rows if r['accepted_mag_0']]
        bad = [r for r in rows if not r['accepted_mag_0'] and r['mag_0_status'] == 'measured']
        ax.scatter([date_from_jd(r['hjd']) for r in good], [r['mag_0'] for r in good],
                   s=10, color='#1c7086', label='A/B')
        ax.scatter([date_from_jd(r['hjd']) for r in bad], [r['mag_0'] for r in bad],
                   s=18, marker='x', color='#b96431', label='Other grades')
        ax.set_title(f'Dataset {dataset["dataset"]}: {dataset["field"]}; {len(good)} A/B rows')
        ax.set_ylim(12.65, 11.8)
        ax.set_ylabel('ASAS V, MAG_0')
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle('J0754: legacy ASAS observations; independent field zero points preserved')
    fig.text(.5, .015, 'Calendar labels from archive HJD. MER is frame scatter, not individual-star uncertainty.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0,.035,1,.96))
    fig.savefig(out / 'asas_legacy_fields.png', dpi=170)
    plt.close(fig)
    print(json.dumps(dict(asas_rows=len(asas), asas_grade_ab=summary['asas']['accepted_mag0_grade_ab'],
                         asas_unique=summary['asas']['accepted_distinct_hjd_frame'],
                         applause_plates=len(plate_map), applause_good_rows=len(good_target),
                         applause_good_plates=summary['applause']['primary_quality_plates']), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--recon', required=True, type=Path)
    parser.add_argument('--pilot', required=True, type=Path)
    args = parser.parse_args()
    run(args.out, args.recon, args.pilot)
