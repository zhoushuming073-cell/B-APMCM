import numpy as np

from apmcm.common.paths import DATA
from apmcm.q1.baseline import surface


def test_bilinear_surface_hits_observed_nodes():
    t = np.array([0, .5, 1.])
    w = np.array([350, 370, 400])
    z = 4 * t[:, None] + 2 * w[None, :]
    fit = surface(t, w, z, 1, 0)
    np.testing.assert_allclose(fit(t, w), z)
    np.testing.assert_allclose(fit.ev(.25, 360), 721)


def test_published_repairs_are_traced_to_masks():
    with np.load(DATA / "processed/scenarios/6C02057_scenario.npz") as z:
        approved = z["approved"]
        upper = z["saturation_repaired"]
        assert not (approved & ~z["screen"]).any()
        assert not (upper & ~z["saturated"]).any()
        np.testing.assert_array_equal(z["scenario"][~(approved | upper)], z["raw"][~(approved | upper)])
