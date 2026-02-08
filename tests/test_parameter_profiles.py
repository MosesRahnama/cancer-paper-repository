import pytest

from simulations.parameter_profiles import (
    CALIBRATED_MONTHLY_PROFILE,
    annual_to_step_probability,
)


def test_annual_to_step_probability_monthly_reference():
    monthly = annual_to_step_probability(0.105, steps_per_year=12)
    assert monthly == pytest.approx(0.0092017, rel=1e-5)


def test_annual_to_step_probability_validation():
    with pytest.raises(ValueError):
        annual_to_step_probability(-0.01)
    with pytest.raises(ValueError):
        annual_to_step_probability(1.01)
    with pytest.raises(ValueError):
        annual_to_step_probability(0.1, steps_per_year=0)


def test_calibrated_profile_has_expected_keys():
    kwargs = CALIBRATED_MONTHLY_PROFILE.run_kwargs()
    assert set(kwargs.keys()) == {"p_mutation", "p_progression", "p_healthy_division"}
