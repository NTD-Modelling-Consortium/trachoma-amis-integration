import os
import sys
from pathlib import Path
import numpy as np

# Add trachoma_amis package to path
trachoma_amis_dir = os.getenv("TRACHOMA_AMIS_DIR", "")
if trachoma_amis_dir:
    sys.path.insert(0, trachoma_amis_dir)

from ntdmc_trachoma.trachoma_functions import (
    get_Intervention_times,
    getOutputTimes,
)
from trachoma_amis.trachoma_params import params as parameters
from trachoma_amis.trachoma_params import sim_params, demog
from trachoma_amis.simulation_setup import (
    get_start_date,
    run_trachoma_simulations,
    PATH_TO_FITTING_PREP_ARTEFACTS,
)

__all__ = ["build_transmission_model"]

START_DATE = get_start_date()

id = os.getenv("SLURM_ARRAY_TASK_ID")
PATH_TO_ENDGAME_INPUTS = (
    PATH_TO_FITTING_PREP_ARTEFACTS / "trachoma" / "data" / "coverage" / "endgame_inputs"
)
mda_filepath = PATH_TO_ENDGAME_INPUTS / f"InputMDA_MTP_{id}.csv"
distToUse = "Expo"


def build_transmission_model(
    fitting_points: list[int], initial_infect_frac=0.1, num_cores=-2
):
    """Create a closure for the AMIS to run the trachoma model.

    Parameters
    ----------
    fitting_points
        A list of indexes of weeks to fit to.
    initial_infect_frac
        Fraction of the population initially infected.
    num_cores
        Number of computer cores to spread transmission model across.

    Returns
    -------
    function
        A function of three arguments (seeds, betavals, n_sims).
    """
    def run_trachoma(seeds, params, n_tims):
        # AMIS-specific parameter unpacking
        betas = [amisPars[0:(int(sim_params["timesim"] / 52))] for amisPars in params]
        coverages = [amisPars[-2] for amisPars in params]
        k_parameters = [amisPars[-1] for amisPars in params]
        
        # Use common simulation engine
        results = run_trachoma_simulations(
            seeds=seeds,
            betas=betas,
            coverages=coverages,
            k_parameters=k_parameters,
            mda_filepath=mda_filepath,
            parameters=parameters,
            sim_params=sim_params,
            demog=demog,
            initial_infect_frac=initial_infect_frac,
            num_cores=num_cores,
            timeToStopImportation=-1,  # Disabled for AMIS
            distToUse=distToUse
        )
        # Get the prevalence from the returned values dictionary
        # This is an exemple of post-processing extracting prevalence among
        # children aged 1 to 9 yo.  To be modified by the modelling team.

        # 'results' is a list of 2-tuples (v, r) where
        # - v is a dict containing various recorded data ('vals' in the trachoma model code)
        # - r is a list of outputResults instances, see trachoma.trachoma_functions.py
        # Must return 2d numpy array with:
        # rows = simulations
        # columns = data points
        return [
            np.transpose(
                np.transpose(
                    [
                        # 'v' is a list, so convert it to NumPy array to index it from
                        # a list 'fitting_points'
                        np.asarray(v["True_Prev_Disease_children_1_9"])[fitting_points]
                        for (v, r) in results
                    ]
                )
            ),
            np.transpose(
                np.transpose(
                    [
                        # 'v' is a list, so convert it to NumPy array to index it from
                        # a list 'fitting_points'
                        np.asarray(v["True_Infections_Disease_children_1_9"])[
                            fitting_points
                        ]
                        for (v, r) in results
                    ]
                )
            ),
        ]

    return run_trachoma
