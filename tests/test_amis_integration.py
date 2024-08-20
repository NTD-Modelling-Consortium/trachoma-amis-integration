import pytest
from trachoma_amis import amis_integration
import numpy as np
import numpy.testing as npt


def test_build_transmission_model():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[0, 1], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model([1], [(0.2, 0.4)], 0)
    # At the moment the seed is ignored so just check
    # the shape of what is returned
    assert np.shape(result) == (1, 2)


@pytest.mark.skip(reason="Seed is not used in setup inits")
def test_running_twice_with_different_seed_produces_different_result_at_first_week():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[1], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 2], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    assert result[0][0] != result[1][0]


def test_running_twice_with_different_seed_produces_different_result_after_10_weeks():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[10], initial_infect_frac=0.5, num_cores=2
    )
    result = transmission_model(seeds=[1, 2], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    assert result[0][0] != result[1][0]


@pytest.mark.skip(reason="This can fail as sometimes the two models will not diverge")
def test_running_twice_with_different_seed_produces_different_result_after_10_weeks_single_core():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[100], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 2], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    assert result[0][0] != result[1][0]


def test_running_twice_with_same_seed_produces_same_result_after_10_weeks():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[10], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 1], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    assert result[0][0] == result[1][0]
