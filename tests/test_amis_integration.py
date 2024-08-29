import copy
import pytest
from trachoma_amis import amis_integration, trachoma_params
import numpy as np
import numpy.testing as npt


def test_build_transmission_model():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[0, 1], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model([1], [(0.2, 0.4)], 0)
    npt.assert_array_equal(result, [[0, 0.40625]])
    assert isinstance(result, np.ndarray)


def test_running_twice_with_different_seed_produces_different_result_at_first_week():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[1], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 2], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)

    assert result[0][0] == 0.40625
    assert result[1][0] == 0.4722222222222222


def test_running_twice_with_different_seed_produces_different_result_after_10_weeks():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[10], initial_infect_frac=0.5, num_cores=2
    )
    result = transmission_model(seeds=[1, 2], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    assert result[0][0] == 0.7575757575757576
    assert result[1][0] == 0.8285714285714286


def test_running_twice_with_different_seed_produces_different_result_after_10_weeks_single_core():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[10], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 2], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    assert result[0][0] == 0.7575757575757576
    assert result[1][0] == 0.8


def test_running_twice_with_same_seed_produces_same_result_after_10_weeks():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[10], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 1], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    assert result[0][0] == result[1][0]


def test_running_with_modified_coverage_eliminates_disease_if_high():
    # This needs to be long enough after the first treatment that
    # all the diseased individuals have recovered
    first_week_after_treatment = 150
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[first_week_after_treatment], initial_infect_frac=1, num_cores=1
    )

    old_params = copy.deepcopy(trachoma_params.params)

    trachoma_params.params["MDA_Eff"] = 1.0
    result = transmission_model(
        seeds=[1, 1], params=[(0.4, 0.9999), (0.4, 0)], n_tims=0
    )
    assert result[0][0] == 0.0
    assert result[1][0] == 1.0

    # Currently the params are a mutable global so modifying them could effect other tests
    # so make sure to restore them
    trachoma_params.params = old_params
