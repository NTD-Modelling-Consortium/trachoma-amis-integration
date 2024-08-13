from trachoma_amis import amis_integration
import numpy as np
import numpy.testing as npt


def test_build_transmission_model():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[0, 1], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model([1], [0.2], 0)
    # At the moment the seed is ignored so this is different each time
    # hence the wide margin on matching
    npt.assert_almost_equal(result, np.array([[0, 0.5]]), decimal=1)


def test_running_twice_with_different_seed_produces_different_result_at_first_week():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[1], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 2], params=[0.2, 0.2], n_tims=0)
    # TODO: this currently fails because the setup is shared between
    # each of the runs so at 1 week they are all in sync
    assert result[0][0] != result[1][0]


def test_running_twice_with_different_seed_produces_different_result_after_10_weeks():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[10], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 2], params=[0.2, 0.2], n_tims=0)
    # TODO: I don't understand why this works, as far as I can tell each run
    # uses the same copied state for each run
    assert result[0][0] != result[1][0]


def test_running_twice_with_same_seed_produces_same_result_after_10_weeks():
    transmission_model = amis_integration.build_transmission_model(
        fitting_points=[10], initial_infect_frac=0.5, num_cores=1
    )
    result = transmission_model(seeds=[1, 1], params=[0.2, 0.2], n_tims=0)
    assert result[0][0] == result[1][0]
