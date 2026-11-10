"""Exploratory rotation search with independent amplitudes/phases per season."""
import csv
import json
from pathlib import Path
import numpy as np
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from orbit_photometry import nuisance, daily

OUT = Path('results/sky-plates-j0754-seasonal-rotation-20261108')
FREQ = np.linspace(1/100, 1/10, 1801)
PERIODS = 1/FREQ


def search(rows, name):
    t = np.array([float(r['t']) for r in rows])
    y = np.array([float(r['y']) for r in rows])
    e = np.maximum(.02, [float(r['err']) for r in rows])
    g = np.array([r['group'] for r in rows])
    dates = Time(t, format='jd').isot
    season = np.array([int(d[:4])-(int(d[5:7]) < 7) for d in dates])
    records, powers, simulations = [], [], {d: [] for d in (3, 7, 14)}
    rng = np.random.default_rng(20261109)
    for year in sorted(set(season)):
        m = season == year
        if sum(m) < 30 or np.ptp(t[m]) < 120:
            continue
        tt, yy, ee, gg = t[m], y[m], e[m], g[m]
        n = nuisance(tt, gg, 'group_linear')/ee[:, None]
        u, s, _ = np.linalg.svd(n, full_matrices=False)
        q = u[:, :sum(s > max(n.shape)*np.finfo(float).eps*s[0])]
        residual = yy/ee-q@(q.T@(yy/ee))
        # At each frequency, project both sinusoids through the same nuisance
        # model and orthonormalize. Power is fractional null-SSE reduction.
        basis = []
        for f in FREQ:
            angle = 2*np.pi*(tt-2457000)*f
            a = np.column_stack((np.cos(angle), np.sin(angle)))/ee[:, None]
            a -= q@(q.T@a)
            b, _ = np.linalg.qr(a, mode='reduced')
            basis.append(b.T)
        basis = np.array(basis)
        power = np.sum((basis@residual)**2, axis=1)/(residual@residual)
        powers.append(power)
        best = int(np.argmax(power))
        near = (PERIODS >= 30) & (PERIODS <= 42)
        ni = np.flatnonzero(near)[np.argmax(power[near])]
        records.append(dict(season=f'{year}/{str(year+1)[2:]}', n=int(sum(m)),
                            best_period_days=float(PERIODS[best]),
                            best_power=float(power[best]),
                            near_period_days=float(PERIODS[ni]),
                            near_power=float(power[ni]),
                            power_at_gaia=float(np.interp(1/35.9204566, FREQ, power))))
        for days in simulations:
            block = np.floor((tt-tt.min())/days).astype(int)
            signs = rng.choice([-1., 1.], size=(block.max()+1, 256))
            z = residual[:, None]*signs[block]
            z -= q@(q.T@z)
            denom = np.sum(z*z, axis=0)
            sim = np.sum((basis@z)**2, axis=1)/denom[None, :]
            simulations[days].append(sim)
    combined = np.mean(powers, axis=0)
    ni = np.flatnonzero((PERIODS >= 30) & (PERIODS <= 42))
    best = int(np.argmax(combined)); nb = ni[np.argmax(combined[ni])]
    diagnostics = {}
    for days, sims in simulations.items():
        sim = np.mean(sims, axis=0)
        diagnostics[days] = dict(draws=256,
            fraction_null_max_anywhere_ge_observed=float(np.mean(sim.max(axis=0) >= combined[best])),
            fraction_null_max_near_ge_observed=float(np.mean(sim[ni].max(axis=0) >= combined[nb])),
            null_near_max_p95=float(np.quantile(sim[ni].max(axis=0), .95)))
    np.savetxt(OUT/f'{name}_powers.csv', np.column_stack([PERIODS, combined, *powers]),
               delimiter=',', header='period_days,mean_power,'+','.join(r['season'] for r in records))
    return dict(seasons=records, best_period_days=float(PERIODS[best]),
                best_mean_power=float(combined[best]), near_period_days=float(PERIODS[nb]),
                near_mean_power=float(combined[nb]), block_sign_diagnostics=diagnostics), powers, combined


def main():
    OUT.mkdir(exist_ok=True)
    protocol = ('Exploratory follow-up specified before this search: 10â€“100 day frequency grid, '
        '30â€“42 day neighborhood of Gaia assessed explicitly. Independent sine/cosine in each '
        'Julyâ€“June season; joint per-camera offsets and linear trends, >=30 daily camera groups '
        'and >=120 day span. Equal-season mean fractional SSE reduction. No clipping or downloads. '
        'Block-sign null diagnostics use shared signs across cameras in 3/7/14-day blocks and '
        '256 draws, scanning the same range; these are sensitivity diagnostics, not calibrated '
        'astrophysical false-alarm probabilities. A recurring peak near 36 days across seasons '
        'and independent ZTF would support rotation; broad/discordant peaks weaken identification. '
        'Season-long stationarity still assumed; this is not a fully evolving spot model.')
    (OUT/'protocol.txt').write_text(protocol, encoding='utf8')
    with open('results/sky-plates-j0754-models-20261108/short_period/asassn_g_daily.csv') as f:
        rows = list(csv.DictReader(f))
    raw = json.loads(Path('results/sky-plates-j0754-checks-20261108/modern_differential.json').read_text())['measurements']
    ztf = daily([dict(t=r['mjd']+2400000.5, y=r['corrected_relative_mag'], err=.02, group='ZTF') for r in raw])
    results = dict(protocol=protocol)
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for ax, name, data in zip(axes, ('ASASSN', 'ZTF'), (rows, ztf)):
        result, powers, mean = search(data, name)
        results[name] = result
        for p in powers:
            ax.plot(PERIODS, p, color='gray', alpha=.3, lw=.7)
        ax.plot(PERIODS, mean, color='navy', label='Mean across seasons')
        ax.axvline(35.9204566, color='darkorange', ls='--', label='Gaia 35.92 d')
        ax.set(title=name, ylabel='Fractional improvement', xlim=(10, 100))
        ax.legend(fontsize=8)
    axes[-1].set_xlabel('Trial period (days); phase and amplitude free each season')
    fig.tight_layout(); fig.savefig(OUT/'seasonal_search.png', dpi=160); plt.close(fig)
    (OUT/'results.json').write_text(json.dumps(results, indent=2), encoding='utf8')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
