from datetime import date
import copy
import random

from joblib import Parallel, delayed
import numpy as np

from trachoma.trachoma_functions import (
    run_single_simulation,
    Set_inits,
    readPlatformData,
    getInterventionDates,
    get_Intervention_times,
)
from .trachoma_params import params, sim_params, demog

__all__ = ["build_transmission_model"]

START_DATE = date(2019, 1, 1)


def set_start_date(datestr: str):
    global START_DATE
    try:
        START_DATE = date.fromisoformat(datestr)
    except ValueError as e:
        msg = (
            f"{e}.\n"
            "Valid dates formats are YYYY-MM-DD, YYYYMMDD or YYYY-WXX-D."
        )
        raise ValueError(msg)


def setup_mda(cov_filepath, burnin):
    MDAData = readPlatformData(cov_filepath, "MDA")
    MDA_dates = getInterventionDates(MDAData)
    MDA_times = get_Intervention_times(MDA_dates, START_DATE, burnin)
    return MDA_times, MDAData


def setup_vaccine(cov_filepath, burnin):
    VaccData = readPlatformData(cov_filepath, "Vaccine")
    Vaccine_dates = getInterventionDates(VaccData)
    vacc_times = get_Intervention_times(Vaccine_dates, START_DATE, burnin)
    return vacc_times, VaccData


def setup(initial_prevalence: float):
    vals = Set_inits(params, demog, sim_params, np.random.get_state())
    ids = random.sample(range(params["N"]), k=int(initial_prevalence * params["N"]))
    vals["IndI"][ids] = 1
    vals["T_latent"][ids] = vals["Ind_latent"][ids]
    vals["No_Inf"][ids] = 1

    MDA_times, MDAData = setup_mda("scen2c.csv", sim_params["burnin"])
    vacc_times, VaccData = setup_vaccine("scen2c.csv", sim_params["burnin"])
    sim_params["N_MDA"] = len(MDA_times)
    sim_params["N_Vaccines"] = len(vacc_times)

    return (
        vals,
        params,
        sim_params,
        demog,
        MDA_times,
        MDAData,
        vacc_times,
        VaccData,
    )


def build_transmission_model(
        fitting_points: list[int],
        initial_infect_frac=0.01,
        num_cores=-2
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
    (
        init_vals,
        params,
        sim_params,
        demog,
        MDA_times,
        MDAData,
        vacc_times,
        VaccData,
    ) = setup(initial_prevalence)

    def run_trachoma(seeds, betavals, n_sims):
        results: list[tuple[dict, list]]
        results = Parallel(n_jobs=num_cores)(
            delayed(run_single_simulation)(
                pickleData=copy.deepcopy(init_vals),
                params=params,
                timesim=sim_params["timesim"],
                burnin=sim_params["burnin"],
                demog=demog,
                beta=beta[0],
                MDA_times=MDA_times,
                MDAData=MDAData,
                vacc_times=vacc_times,
                VaccData=VaccData,
                outputTimes=[],
                index=i,
                numpy_state=np.random.get_state(),
            )
            for i, beta in enumerate(betavals)
        )
        # Get the prevalence from the returned values dictionary
        # Could probably also compute it again here.
        return np.transpose(
            [
                t[0]['True_Prev_Disease_children_1_9'][fitting_points]
                for t in results
            ]
        )
    return run_trachoma
