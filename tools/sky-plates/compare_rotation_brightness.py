"""Compare seasonal rotational amplitude with brightness using retained data."""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from orbit_photometry import nuisance, daily, PERIOD

OUT = Path('results/sky-plates-j0754-amplitude-brightness-20261108')
BASE = Path('results/sky-plates-j0754-seasonal-rotation-20261108')


def amplitude(t, y, e, groups, period):
    phase = 2*np.pi*(t-2457000)/period
    n = nuisance(t, groups, 'group_linear')
    design = np.column_stack([n, np.cos(phase), np.sin(phase)])
    weighted = design/e[:, None]
    beta = np.linalg.lstsq(weighted, y/e, rcond=None)[0]
    return float(2*np.linalg.norm(beta[-2:]))


def main():
    OUT.mkdir(exist_ok=True)
    protocol = ('Exploratory comparison specified before calculation: use same July-June '
        'seasons retained by previous rotation search. Compare peak-to-peak sinusoidal '
        'amplitude with median brightness, using ASAS-SN bq only for brightness to avoid '
        'camera changes, and independent differential ZTF g. Amplitude joint-fits camera '
        'offsets/slopes at fixed35.9204566days and separately at previously selected '
        'seasonal30-42day peak. Selected-period amplitudes have selection bias. Spearman '
        'rank correlation is descriptive across few non-independent seasons, no significance '
        'claim. Simple prediction: increasing rotational amplitude as star fades, decreasing '
        'during recovery. Failure weakens this simple relation, not all spot geometries. '
        'Amplitude measures longitudinal asymmetry, not total spot coverage.')
    (OUT/'protocol.txt').write_text(protocol, encoding='utf8')
    prior = json.loads((BASE/'results.json').read_text())
    with open('results/sky-plates-j0754-models-20261108/short_period/asassn_g_daily.csv') as f:
        asas = list(csv.DictReader(f))
    raw = json.loads(Path('results/sky-plates-j0754-checks-20261108/modern_differential.json').read_text())['measurements']
    ztf = daily([dict(t=r['mjd']+2400000.5, y=r['corrected_relative_mag'],
                      err=.02, group='ZTF') for r in raw])
    results = dict(protocol=protocol)
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex='col')
    for col, (name, rows) in enumerate((('ASASSN', asas), ('ZTF', ztf))):
        t = np.array([float(r['t']) for r in rows])
        y = np.array([float(r['y']) for r in rows])
        e = np.maximum(.02, [float(r['err']) for r in rows])
        g = np.array([r['group'] for r in rows])
        dates = Time(t, format='jd').isot
        years = np.array([int(d[:4])-(int(d[5:7]) < 7) for d in dates])
        records = []
        for item in prior[name]['seasons']:
            year = int(item['season'][:4]); m = years == year
            bright = m & (g == 'bq') if name == 'ASASSN' else m
            if sum(bright) < 10:
                continue
            records.append(dict(season=item['season'], year=year+1, n=int(sum(m)),
                brightness_n=int(sum(bright)), median_mag=float(np.median(y[bright])),
                fixed_period_peak_to_peak_mag=amplitude(t[m], y[m], e[m], g[m], PERIOD),
                selected_period_days=item['near_period_days'],
                selected_period_peak_to_peak_mag=amplitude(t[m], y[m], e[m], g[m], item['near_period_days'])))
        b = [r['median_mag'] for r in records]
        fixed = [r['fixed_period_peak_to_peak_mag'] for r in records]
        chosen = [r['selected_period_peak_to_peak_mag'] for r in records]
        results[name] = dict(seasons=records,
            descriptive_spearman_magnitude_vs_fixed_amplitude=float(spearmanr(b, fixed).statistic),
            descriptive_spearman_magnitude_vs_selected_amplitude=float(spearmanr(b, chosen).statistic))
        xx = [r['year'] for r in records]
        axes[0,col].plot(xx, b, 'o-', color='black')
        axes[0,col].invert_yaxis()
        axes[0,col].set_title(name+' (bq camera brightness)' if name == 'ASASSN' else 'ZTF differential g')
        axes[0,col].set_ylabel('Median g magnitude\nbrighter toward top')
        axes[1,col].plot(xx, fixed, 'o-', label='35.92 days')
        axes[1,col].plot(xx, chosen, 's--', label='Seasonal nearby peak')
        axes[1,col].set(ylabel='Rotational peak-to-peak (mag)', xlabel='Observing season ending year', ylim=(0,None))
        axes[1,col].legend(fontsize=8)
    fig.suptitle('Does the rotational variation grow as J0754 fades?')
    fig.tight_layout(); fig.savefig(OUT/'amplitude_vs_brightness.png', dpi=160); plt.close(fig)
    (OUT/'results.json').write_text(json.dumps(results, indent=2), encoding='utf8')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
