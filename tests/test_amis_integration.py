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
    npt.assert_almost_equal(result, np.array([[0], [0.5]]), decimal=1)
