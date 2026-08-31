import pytest

import lsdtt3


def test_crn_column_steady_state_round_trip():
    column = lsdtt3.CRNColumn()
    sample = column.steady_state(effective_erosion_rate_g_cm2_yr=0.01)
    inferred = lsdtt3.crn_apparent_erosion_rate("Be10", sample.concentration_be10)
    assert inferred["status"] == "Converged"
    assert inferred["erosion_rate_g_cm2_yr"] == pytest.approx(0.01, rel=1e-6)


def test_crn_column_transient_profile():
    column = lsdtt3.CRNColumn()
    column.set_erosion_history([0.0, 1000.0, -9999.0], [0.0001, 0.0001, 0.0001], [0.0, 0.0, 0.0])
    profile = column.transient([0.0, 1.0])
    assert len(profile) == 2
    assert profile[0].concentration_be10 > profile[1].concentration_be10
