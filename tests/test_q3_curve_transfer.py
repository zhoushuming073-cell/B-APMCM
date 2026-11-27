import numpy as np
from apmcm.common.paths import PROC
from apmcm.q3.curve_transfer import Experiment, get_rows, prepare


def load_compact(path):
    with np.load(path) as a:
        # Sample physical channels for a lightweight leakage regression test.
        s = np.log1p(a["x"][:, :32])
        return a["oxygen"], a["gas"], s


def test_future_process_and_spectrum_cannot_change_earlier_issued_features():
    o, y, s = load_compact(sorted(PROC.glob("raw_*.npz"))[0])
    initial = prepare(o, y, s, .5)
    cutoff = 180
    o[cutoff+1:] += 20.
    s[cutoff+1:] += 2.
    changed = prepare(o, y, s, .5)
    np.testing.assert_array_equal(initial.prefix,changed.prefix)
    for a,b in zip(get_rows(initial,60),get_rows(changed,60)):
        idx = get_rows(initial,60)[1] <= cutoff
        np.testing.assert_array_equal(a[idx],b[idx])
    np.testing.assert_array_equal(initial.spectrum[:cutoff+1],changed.spectrum[:cutoff+1])


def test_test_heat_gas_has_no_effect_on_lag_or_forecast():
    from threadpoolctl import threadpool_limits
    files = sorted(PROC.glob("raw_*.npz"))[:4]
    data = {p.stem:prepare(*load_compact(p),.5) for p in files}
    held, *train = data
    train = tuple(train)
    with threadpool_limits(1):
        experiment = Experiment(data)
        lag = experiment.lead(train,0,"q2")
        before = experiment.predict(train,held,0,lag)[2]
        data[held].target[:] = -12345
        modified = Experiment(data)
        assert modified.lead(train,0,"q2") == lag
        after = modified.predict(train,held,0,lag)[2]
    for method in before:
        np.testing.assert_allclose(before[method],after[method],atol=1e-10)
    experiment.clear()
