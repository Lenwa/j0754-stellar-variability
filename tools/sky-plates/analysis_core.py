"""Numerical contracts for the experimental historical-plate investigation.

DASCH flag meanings follow the MIT-licensed pkgw/daschlab photometry module;
the code here independently implements the documented numeric selection.
"""
import math
import random
import statistics as st

BASELINES = ((1910, 1930), (1950, 1970))
EVENT = (1930, 1940)
STANDARD_BAD_AFLAGS = sum(1 << n for n in (6, 11, 12, 13, 15))


def number(value, default=float("nan")):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def separation_arcsec(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(math.radians, (ra1, dec1, ra2, dec2))
    a = math.sin((dec2-dec1)/2)**2 + math.cos(dec1)*math.cos(dec2)*math.sin((ra2-ra1)/2)**2
    return math.degrees(2*math.asin(min(1, math.sqrt(a)))) * 3600


def in_baseline(year):
    return any(lo <= year < hi for lo, hi in BASELINES)


def eligible_catalog(row, target):
    values = [number(row.get(k)) for k in ("stdmag", "color", "num_matches")]
    return (row["ref_number"] != target["ref_number"]
            and all(math.isfinite(v) for v in values)
            and values[2] >= 300
            and abs(values[0]-number(target["stdmag"])) <= 1
            and abs(values[1]-number(target["color"])) <= .3)


def eligible_applause(row, target):
    """Match static catalogue colour and synthetic magnitude in the same image."""
    values = [number(obj.get(k)) for obj in (row,target) for k in ("cat_natmag","gaiaedr3_bp_rp")]
    return (all(math.isfinite(v) for v in values)
            and abs(values[0]-values[2]) <= 1
            and abs(values[1]-values[3]) <= .3)


def quality_dasch(row, strict=False):
    mag = number(row.get("magcal_magdep"))
    err = number(row.get("magcal_local_rms"))
    limit = number(row.get("limiting_mag_local"))
    values = [number(row.get(k)) for k in ("ra_deg", "dec_deg", "ra_cat_corrected", "dec_cat_corrected")]
    if not all(math.isfinite(v) for v in [mag, err, limit, *values]):
        return False
    flags = int(row.get("aflags") or 0)
    mask = STANDARD_BAD_AFLAGS
    if strict:
        mask |= (1 << 7) | (1 << 8) | (1 << 10) | (1 << 25)
        if int(row.get("reject_flag") or 0):
            return False
    return flags & mask == 0 and 0 < err <= .3 and limit-mag >= .5 and separation_arcsec(*values) <= 10


def differential_curve(target, controls, minimum_controls=5):
    """Use an ensemble of same-exposure control residuals, without clipping target."""
    by_key = {}
    for label, measurements in controls.items():
        baseline = [r["mag"] for r in measurements if in_baseline(r["year"])]
        if len(baseline) < 3:
            continue
        center = st.median(baseline)
        grouped = {}
        for row in measurements:
            grouped.setdefault(row["key"], []).append(row["mag"] - center)
        for key, offsets in grouped.items():
            by_key.setdefault(key, {})[label] = st.median(offsets)
    result = []
    for row in target:
        offsets = list(by_key.get(row["key"], {}).values())
        if len(offsets) < minimum_controls:
            continue
        offset = st.median(offsets)
        scatter = 1.4826 * st.median(abs(v-offset) for v in offsets)
        result.append(dict(row, corrected_mag=row["mag"]-offset,
                           control_offset=offset, control_scatter=scatter,
                           controls=len(offsets)))
    return result


def annual_contrast(rows, replicates=2000, minimum_per_year=3, value_key="corrected_mag"):
    """Median event year minus median baseline year; bootstrap calendar years."""
    grouped = {}
    for row in rows:
        grouped.setdefault(math.floor(row["year"]), []).append(row[value_key])
    annual = {year: st.median(values) for year, values in grouped.items() if len(values) >= minimum_per_year}
    event = [v for y, v in annual.items() if EVENT[0] <= y < EVENT[1]]
    base = [v for y, v in annual.items() if in_baseline(y)]
    result = {"event_years": len(event), "baseline_years": len(base),
              "event_points": sum(len(v) for y,v in grouped.items() if EVENT[0] <= y < EVENT[1]),
              "baseline_points": sum(len(v) for y,v in grouped.items() if in_baseline(y)),
              "delta_mag": None, "bootstrap_95_percent": None}
    if len(event) < 2 or len(base) < 2:
        return result
    rng = random.Random(20261107)
    boots = sorted(st.median(rng.choices(event, k=len(event))) - st.median(rng.choices(base, k=len(base))) for _ in range(replicates))
    result.update(delta_mag=st.median(event)-st.median(base),
                  bootstrap_95_percent=[boots[int(.025*(replicates-1))], boots[int(.975*(replicates-1))]])
    return result
