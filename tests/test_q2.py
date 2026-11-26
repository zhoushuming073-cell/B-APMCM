import numpy as np
from scipy.integrate import trapezoid

from apmcm.q2.delay import design
from apmcm.q2.reconstruct import AveragedCounter


def test_oxygen_counter_endpoints_and_mass():
    t = np.arange(0, 20.5, .5)
    y = 7 + 2*t + .5*np.maximum(0, t - 4)
    counter = AveragedCounter(t, y, 4)
    fine = np.linspace(0, 20, 2001)
    flow = counter.derivative()(fine)
    np.testing.assert_allclose(counter(t[[0,-1]]), y[[0,-1]])
    assert np.min(flow) >= 0
    np.testing.assert_allclose(trapezoid(flow, fine), y[-1]-y[0])


def test_clock_reduces_to_constant_delay_at_constant_speed():
    t = np.arange(0, 160, 2)
    q = np.full_like(t, 10.)
    stage = np.ones((len(t), 1))
    _, tau = design((t, q, stage, None), ("shift", 24., np.log(2)))
    np.testing.assert_allclose(tau, 12.)
