import pytest
from trachoma_amis import amis_integration
import numpy as np
import numpy.testing as npt


def test_build_transmission_model():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[0, 1], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model([1], [(0.2, 0.4)], 0)
    # At the moment the seed is ignored so this is different each time
    # hence the wide margin on matching
    npt.assert_almost_equal(result, np.array([[0, 0.5]]), decimal=1)


@pytest.mark.skip(reason="Seed is not used in setup inits")
def test_running_twice_with_different_seed_produces_different_result_at_first_week():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[1], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 2], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    assert result[0][0] != result[1][0]


@pytest.mark.skip(
    reason="When running with multiple cores, the same state is used for each process"
)
def test_running_twice_with_different_seed_produces_different_result_after_10_weeks():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[10], initial_infect_frac=0.5, num_cores=2
    )
    result = transmission_model(seeds=[1, 2], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    assert result[0][0] != result[1][0]


def test_running_twice_with_different_seed_produces_different_result_after_10_weeks_single_core():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[10], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 2], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    # This works only becuase the single core means the state is shared from one run to the next
    assert result[0][0] != result[1][0]


@pytest.mark.skip(reason="Seed is not used to run simulations")
def test_running_twice_with_same_seed_produces_same_result_after_10_weeks():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[10], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 1], params=[(0.2, 0.4), (0.2, 0.4)], n_tims=0)
    assert result[0][0] == result[1][0]
