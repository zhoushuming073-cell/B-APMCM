"""Histogram gradient boosting of cross-fitted curve residuals."""
import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from .curve_transfer import GRID


def predict(experiment, train, held, gas, lead):
    if held in train:
        raise ValueError("Test heat must be excluded")
    end, idx, phase, process = experiment.rows(held, lead)
    heat = experiment.heats[held]
    baseline = np.interp(phase, GRID, experiment.shared(train, gas, lead)[0])
    inputs, targets, weights = [], [], []
    for h in train:
        a = experiment.heats[h]
        e, i, ph, p = experiment.rows(h, lead)
        others = tuple(k for k in train if k != h)
        residual = a.target[e, gas] - np.interp(ph, GRID, experiment.shared(others, gas, lead)[0])
        step = int(2/a.dt)
        inputs.append(np.column_stack([p, a.spectrum[i]])[::step])
        targets.append(residual[::step])
        weights.append(np.full(len(inputs[-1]), 1/len(inputs[-1])))
    X = np.vstack(inputs)
    y = np.concatenate(targets)
    w = np.concatenate(weights)
    w /= w.mean()
    test = np.column_stack([process, heat.spectrum[idx]])
    forecasts = {}
    for name, columns in [("boost_process", slice(0, process.shape[1])), ("boost_fusion", slice(None))]:
        model = HistGradientBoostingRegressor(
            loss="squared_error", learning_rate=.04, max_iter=160,
            max_leaf_nodes=7, min_samples_leaf=30, l2_regularization=10.,
            early_stopping=False, random_state=20261127)
        model.fit(X[:, columns], y, sample_weight=w)
        forecasts[name] = np.maximum(0, baseline + .5*model.predict(test[:, columns]))
    return end, idx, forecasts
