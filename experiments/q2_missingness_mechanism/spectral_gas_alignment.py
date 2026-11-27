#!/usr/bin/env python3
"""Test whether spectral–flue-gas lag structure distinguishes shared from independent row loss."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.signal import correlate


def max_lag_corr(x, y, max_lag=60):
    x = np.diff(np.asarray(x, float))
    y = np.diff(np.asarray(y, float))
    n = min(len(x), len(y))
    x, y = x[:n], y[:n]
    cross = correlate(y, x, mode="full", method="fft")
    px, py = np.r_[0.0, np.cumsum(x)], np.r_[0.0, np.cumsum(y)]
    px2, py2 = np.r_[0.0, np.cumsum(x*x)], np.r_[0.0, np.cumsum(y*y)]
    out = []
    for lag in range(-max_lag, max_lag + 1):
        if lag < 0:
            x0, x1, y0, y1 = -lag, n, 0, n + lag
        else:
            x0, x1, y0, y1 = 0, n - lag, lag, n
        count = x1 - x0
        if count <= 20:
            continue
        sx, sy = px[x1]-px[x0], py[y1]-py[y0]
        sxx, syy = px2[x1]-px2[x0], py2[y1]-py2[y0]
        sxy = float(cross[n - 1 + lag])
        cov = sxy - sx*sy/count
        vx, vy = max(sxx-sx*sx/count, 0.0), max(syy-sy*sy/count, 0.0)
        c = cov / np.sqrt(vx*vy) if vx > 0 and vy > 0 else np.nan
        if np.isfinite(c):
            out.append((float(c), lag))
    if not out:
        return np.nan, np.nan, np.nan
    corr0 = next((c for c, lag in out if lag == 0), np.nan)
    best, lag = max(out, key=lambda z: abs(z[0]))
    return corr0, lag, abs(best)


def series_from_combined_cache(z):
    X = z['X'].astype(float)
    wl = z['wl'].astype(float)
    row_sum = X.sum(axis=1)
    return {
        'spec_integral': np.trapezoid(X, wl, axis=1),
        'spec_centroid': (X @ wl) / np.maximum(row_sum, 1e-12),
        'spec_p95': np.quantile(X, .95, axis=1),
        'spec_max': X.max(axis=1),
        'CO2': z['CO2'].astype(float),
        'CO': z['CO'].astype(float),
    }


def load_repo_inputs(spectra_dir: Path, process_csv: Path):
    proc = pd.read_csv(process_csv)
    data = {}
    for p in sorted(spectra_dir.glob('raw_6C*.npz')):
        heat = p.stem.replace('raw_', '')
        with np.load(p) as z:
            X = z['x'].astype(float)
            wl = z['wl'].astype(float)
            seq = z['sequence'].astype(int)
        g = proc[proc.heat == heat].sort_values('row')
        if len(g) != len(seq) or not np.array_equal(g.row.to_numpy(int), seq):
            raise ValueError(f'Row alignment check failed for {heat}: spectra={len(seq)}, process={len(g)}')
        row_sum = X.sum(axis=1)
        data[heat] = {
            'spec_integral': np.trapezoid(X, wl, axis=1),
            'spec_centroid': (X @ wl) / np.maximum(row_sum, 1e-12),
            'spec_p95': np.quantile(X, .95, axis=1),
            'spec_max': X.max(axis=1),
            'CO2': g.CO2.to_numpy(float),
            'CO': g.CO.to_numpy(float),
        }
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache-dir', type=Path, help='Combined per-heat NPZ caches for local runs')
    ap.add_argument('--spectra-dir', type=Path, help='Repository raw spectra NPZ directory')
    ap.add_argument('--process-csv', type=Path, help='Repository process.csv with row, CO2, CO, heat')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--rep', type=int, default=25)
    ap.add_argument('--max-lag', type=int, default=60)
    ap.add_argument('--seed', type=int, default=20261128)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    if args.cache_dir:
        data = {p.stem: series_from_combined_cache(np.load(p)) for p in args.cache_dir.glob('6C*.npz')}
    elif args.spectra_dir and args.process_csv:
        data = load_repo_inputs(args.spectra_dir, args.process_csv)
    else:
        ap.error('Provide --cache-dir or both --spectra-dir and --process-csv')

    target = '6C02057'
    if len(data) != 6 or target not in data:
        raise ValueError(f'Expected six 6C heats, found {sorted(data)}')
    controls = {h: d for h, d in data.items() if h != target}
    rng = np.random.default_rng(args.seed)
    pairs = [(s, g) for s in ['spec_integral','spec_centroid','spec_p95','spec_max'] for g in ['CO2','CO']]
    rows = []

    for sf, gas in pairs:
        c0, lag, cm = max_lag_corr(data[target][sf], data[target][gas], args.max_lag)
        rows.append({'heat':target,'scenario':'target','rep':0,'spectral_feature':sf,'gas':gas,
                     'corr0':c0,'best_lag_rows':lag,'max_abs_corr':cm})

    for h, d in controls.items():
        n = len(d['CO2']); m = min(775, n-1)
        for scenario in ['global','independent']:
            for rep in range(args.rep):
                if scenario == 'global':
                    ix = np.sort(rng.choice(n, size=m, replace=False))
                    ix_s = ix_g = ix
                else:
                    ix_s = np.sort(rng.choice(n, size=m, replace=False))
                    ix_g = np.sort(rng.choice(n, size=m, replace=False))
                for sf, gas in pairs:
                    c0, lag, cm = max_lag_corr(d[sf][ix_s], d[gas][ix_g], args.max_lag)
                    rows.append({'heat':h,'scenario':scenario,'rep':rep,'spectral_feature':sf,'gas':gas,
                                 'corr0':c0,'best_lag_rows':lag,'max_abs_corr':cm})

    R = pd.DataFrame(rows)
    R.to_csv(args.out/'spectral_gas_alignment_simulations.csv', index=False)
    target_rows = R[R.heat == target]
    summary = []
    for (sf, gas), t in target_rows.groupby(['spectral_feature','gas']):
        tv, tl = float(t.max_abs_corr.iloc[0]), float(t.best_lag_rows.iloc[0])
        for scenario in ['global','independent']:
            q = R[(R.heat != target) & (R.scenario == scenario) & (R.spectral_feature == sf) & (R.gas == gas)]
            per_heat = q.groupby('heat').max_abs_corr.median().dropna()
            pooled = q.max_abs_corr.dropna().to_numpy()
            summary.append({'spectral_feature':sf,'gas':gas,'scenario':scenario,
                            'target_max_abs_corr':tv,'target_best_lag_rows':tl,
                            'control_heat_median_of_medians':float(per_heat.median()),
                            'control_heat_min':float(per_heat.min()),'control_heat_max':float(per_heat.max()),
                            'pooled_sim_p10':float(np.quantile(pooled,.1)),
                            'pooled_sim_p90':float(np.quantile(pooled,.9)),
                            'target_inside_p10_p90':bool(np.quantile(pooled,.1)<=tv<=np.quantile(pooled,.9)),
                            'n_control_heats':int(len(per_heat))})
    S = pd.DataFrame(summary)
    S.to_csv(args.out/'spectral_gas_alignment_summary.csv', index=False)
    print(S.to_string(index=False))

if __name__ == '__main__':
    main()
