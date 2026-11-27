"""Registered heat curves and cross-fitted residual predictors for causal gas forecasts."""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.ndimage import gaussian_filter1d
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler

GASES = ("CO2", "CO")
GRID = np.arange(0., 6001., 10.)


def bands(wl, x):
    edges = np.linspace(345, 1015, 17)
    variables = []
    for lo, hi in zip(edges, edges[1:]):
        v = x[:, (wl >= lo) & (wl < hi)]
        if not v.shape[1]:
            raise ValueError("Missing physical wavelength band")
        variables.append(np.log1p(v.mean(axis=1)))
        high, median = np.percentile(v, [95, 50], axis=1)
        variables.append(np.log1p(high)-np.log1p(median))
    return np.column_stack(variables)


def trailing(a, n):
    prefix = np.vstack([np.zeros(a.shape[1]), np.cumsum(a, axis=0)])
    stop = np.arange(1, len(a)+1)
    start = np.maximum(0, stop-n)
    return (prefix[stop]-prefix[start])/(stop-start)[:, None]


@dataclass
class Heat:
    t: np.ndarray
    dose: np.ndarray
    flow: np.ndarray
    spectrum: np.ndarray
    prefix: np.ndarray
    target: np.ndarray
    dt: float


def prepare(oxygen, y, spec, dt):
    t = np.arange(len(y))*dt
    i = np.arange(len(t))
    dose = oxygen-oxygen[0]
    flow = []
    for sec in [10, 30, 60]:
        j = np.maximum(0, i-int(sec/dt))
        flow.append((dose-dose[j])/np.maximum(dt, t-t[j]))
    q = np.column_stack(flow)
    short = trailing(spec, int(5/dt))
    long = trailing(spec, int(30/dt))
    initial = t <= 60
    pref_band = short[initial, ::2].mean(0).reshape(4, 4).mean(1)
    prefix = np.r_[q[initial,0].mean(), q[initial,0].std(), dose[initial][-1]/60,
                   pref_band, spec[initial,1::2].mean(), short[initial,::2].mean(1).std()]
    spectrum = np.column_stack([short, short-long, short-short[initial].mean(0)])
    return Heat(t, dose, q, spectrum, prefix, y, dt)


def get_rows(a, lag):
    end = np.flatnonzero(a.t >= 120.)
    idx = end-int(lag/a.dt)
    if a.t[idx[0]] < 60:
        raise ValueError("Initial spectrum unavailable at forecast issue")
    t = a.t[idx]
    dose = a.dose[idx]
    # Continue the net dose with past mean supply only; preserve phase order.
    phase = np.maximum.accumulate(dose + 1.25*lag*dose/t)
    p = np.column_stack([phase/3000, dose/3000, t/600, a.flow[idx]/5, dose/t/5])
    return end, idx, phase, p


class Experiment:
    def __init__(self, heats):
        self.heats = heats

    @lru_cache(None)
    def lag_stats(self, h, g, d):
        a = self.heats[h]
        end, idx, _, _ = self.rows(h, d)
        end, idx = end[::4], idx[::4]
        u = a.dose[end]/3000
        q = a.flow[idx, 0]/5
        X = np.column_stack([np.ones(len(u)), a.t[end]/600, u,
                             np.maximum(u[:, None]-[.15,.4,.7,.9], 0), q, q*q])
        y = a.target[end, g]
        return X, y, X.T@X/len(X), X.T@y/len(X)

    @lru_cache(None)
    def rows(self, h, lag):
        return get_rows(self.heats[h], lag)

    def lead(self, train, g, scheme):
        if scheme == "30s":
            return 30
        losses = []
        for d in range(6, 61, 6):
            stats = [self.lag_stats(h,g,d) for h in train]
            penalty = np.eye(stats[0][2].shape[0])*.001
            penalty[0,0] = 0.
            A = sum(s[2] for s in stats) + (len(train)-1)*penalty
            b = sum(s[3] for s in stats)
            losses.append(np.mean([np.abs(y-X@np.linalg.solve(A-aa,b-bb)).mean()
                                   for X,y,aa,bb in stats]))
        return 6*(1+int(np.argmin(losses)))

    @lru_cache(None)
    def curve(self, h, g, lag):
        end, _, phase, _ = self.rows(h, lag)
        keep = np.r_[np.diff(phase)>0, True]
        if phase[-1] > GRID[-1]:
            raise ValueError("Increase registered dose grid")
        v = np.interp(GRID, phase[keep], self.heats[h].target[end[keep],g])
        return gaussian_filter1d(v, 6, mode="nearest"), (GRID>=phase[0]) & (GRID<=phase[-1])

    @lru_cache(None)
    def shared(self, train, g, lag):
        curves, valid = map(np.stack, zip(*(self.curve(h,g,lag) for h in train)))
        count = valid.sum(0)
        mean = curves.mean(0)
        np.divide((curves*valid).sum(0), count, where=count>0, out=mean)
        return mean, (curves-mean)*valid

    def shape_prediction(self, train, held, g, lag, use_spectrum):
        mean, dev = self.shared(train,g,lag)
        _, _, vt = np.linalg.svd(dev, full_matrices=False)
        modes = vt[:2]
        z = np.stack([self.heats[h].prefix for h in train])
        new = self.heats[held].prefix
        if not use_spectrum:
            z, new = z[:, :3], new[:3]
        center, scale = z.mean(0), np.maximum(z.std(0),1e-6)
        model = Ridge(alpha=10.).fit((z-center)/scale, dev @ modes.T)
        scores = model.predict(np.clip((new-center)/scale,-3,3)[None])[0]
        return np.interp(self.rows(held,lag)[2],GRID, mean + scores@modes)

    def predict(self, train, held, g, lag):
        end, idx, phase, proc = self.rows(held,lag)
        base = np.interp(phase, GRID, self.shared(train,g,lag)[0])
        pred = {"shared": base, "fpca_process":self.shape_prediction(train,held,g,lag,False),
                "fpca_fusion":self.shape_prediction(train,held,g,lag,True)}
        xs, spec, residuals, weights = [], [], [], []
        for h in train:
            a = self.heats[h]
            e, i, ph, p = self.rows(h,lag)
            # Each training heat's residual uses a template from the other training heats.
            others = tuple(k for k in train if k != h)
            target = a.target[e,g]-np.interp(ph, GRID, self.shared(others,g,lag)[0])
            step = int(2/a.dt)
            xs.append(p[::step])
            spec.append(a.spectrum[i[::step]])
            residuals.append(target[::step])
            weights.append(np.full(len(xs[-1]),1/len(xs[-1])))
        x = np.vstack(xs)
        y = np.concatenate(residuals)
        w = np.concatenate(weights)
        w = w / w.mean()
        for inputs,xtest,xtrain in [
                ("process",proc,x),
                ("fusion",np.column_stack([proc,self.heats[held].spectrum[idx]]),
                 np.column_stack([x,np.vstack(spec)]))]:
            model = ExtraTreesRegressor(n_estimators=96,min_samples_leaf=20,max_features=.65,
                                        random_state=20261127,n_jobs=2)
            model.fit(xtrain,y,sample_weight=np.concatenate(weights))
            pred["tree_"+inputs] = base + .5*model.predict(xtest)
            scaler = StandardScaler().fit(xtrain,sample_weight=w)
            X = np.clip(scaler.transform(xtrain),-5,5)
            y_mean = np.average(y,weights=w)
            y_sd = max(np.sqrt(np.average((y-y_mean)**2,weights=w)), 1e-8)
            net = MLPRegressor(hidden_layer_sizes=(16,8), activation="tanh", alpha=.1,
                               batch_size=len(X),learning_rate_init=.008,random_state=20261127)
            for _ in range(60):
                net.partial_fit(X,(y-y_mean)/y_sd,sample_weight=w)
            correction = y_mean + y_sd*net.predict(np.clip(scaler.transform(xtest),-5,5))
            pred["mlp_"+inputs] = base + .5*correction
        return end, idx, {k:np.maximum(0,v) for k,v in pred.items()}

    def clear(self):
        for cached in [self.rows,self.curve,self.shared,self.lag_stats]:
            cached.cache_clear()
