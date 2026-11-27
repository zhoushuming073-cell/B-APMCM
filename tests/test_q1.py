import numpy as np
import pytest

from apmcm.q1.figures import display_indices
from apmcm.q1.reconstruct import (
    fit_scales, options, repair_sat, score, validate_scenario_variants,
)


def test_multiscale_score_has_expected_shape_and_finite_values():
    rng = np.random.default_rng(17)
    base = np.linspace(800, 1600, 96, dtype=np.float32)
    spectra = np.maximum(
        0, base[None, :] + rng.normal(0, 4, size=(80, 96))
    ).astype(np.float32)
    scales = fit_scales(spectra[:55])
    result = score(spectra, scales)
    assert result.shape == spectra.shape
    assert np.isfinite(result).all()
    assert (result >= 0).all()


def test_options_require_two_sided_template_support():
    frames, channels = 40, 80
    spectra = np.full((frames, channels), 1200.0, dtype=np.float32)
    spectra += np.linspace(-20, 20, channels, dtype=np.float32)[None, :]
    spectra[20, 40] = 2400.0
    mask = np.zeros_like(spectra, dtype=bool)
    mask[20, 40] = True
    value, difference, agreement, gap = options(spectra, mask)
    assert np.isfinite(value[20, 40])
    assert agreement[20, 40]
    assert difference[20, 40] > 0
    assert gap[20, 40] < 1e-5


def test_saturation_template_is_not_misrepresented_as_hard_lower_bound():
    spectra = np.full((5, 20), 1000.0, dtype=np.float32)
    spectra[1] *= 0.95
    spectra[3] *= 1.05
    saturated = np.zeros_like(spectra, dtype=bool)
    saturated[2, 8:10] = True
    spectra[saturated] = 65535.0

    template_scenario, repaired = repair_sat(spectra, saturated)
    assert repaired.sum() == 2
    assert (template_scenario[repaired] < 65535.0).all()
    lower_bound = template_scenario.copy()
    lower_bound[saturated] = np.maximum(65535.0, lower_bound[saturated])
    assert (lower_bound[saturated] >= 65535.0).all()
    validate_scenario_variants(
        spectra, spectra, template_scenario, lower_bound,
        np.zeros_like(saturated), np.zeros_like(saturated),
        saturated, repaired,
    )


def test_scenario_provenance_masks_are_enforced():
    raw = np.full((4, 8), 100.0, dtype=np.float32)
    screen = np.zeros_like(raw, dtype=bool)
    screen[1, 3] = True
    approved = screen.copy()
    saturated = np.zeros_like(raw, dtype=bool)
    saturated[2, 4] = True
    repaired = saturated.copy()
    anomaly_only = raw.copy()
    anomaly_only[1, 3] = 120.0
    template = anomaly_only.copy()
    template[2, 4] = 200.0
    lower_bound = template.copy()
    lower_bound[2, 4] = 65535.0
    validate_scenario_variants(
        raw, anomaly_only, template, lower_bound,
        screen, approved, saturated, repaired,
    )
    with pytest.raises(ValueError):
        validate_scenario_variants(
            raw, anomaly_only, template, lower_bound,
            screen, np.zeros_like(approved), saturated, repaired,
        )


def test_display_indices_preserve_all_wavelength_channels_by_default():
    np.testing.assert_array_equal(display_indices(2047, 1), np.arange(2047))
    assert display_indices(21, 5)[-1] == 20
    with pytest.raises(ValueError):
        display_indices(21, 0)
