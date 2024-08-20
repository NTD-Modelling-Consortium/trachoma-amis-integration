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
    getOutputTimes,
)
from .trachoma_params import params as parameters
from .trachoma_params import sim_params, demog

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
    MDA_times, MDAData = setup_mda("scen2c.csv", sim_params["burnin"])
    vacc_times, VaccData = setup_vaccine("scen2c.csv", sim_params["burnin"])
    sim_params["N_MDA"] = len(MDA_times)
    sim_params["N_Vaccines"] = len(vacc_times)

    return (
        parameters,
        sim_params,
        demog,
        MDA_times,
        MDAData,
        vacc_times,
        VaccData,
    )


def create_initial_population(initial_prevalence: float, MDAData):
    vals = Set_inits(parameters, demog, sim_params, MDAData, np.random.get_state())
    ids = np.random.choice(
        range(parameters["N"]), int(initial_prevalence * parameters["N"]), replace=False
    )
    vals["IndI"][ids] = 1
    vals["T_latent"][ids] = vals["Ind_latent"][ids]
    vals["No_Inf"][ids] = 1
    return vals


def alterMDACoverage(MDAData, coverage):
    """ update the coverage of each MDA for a run to be a given value.
    Parameters
    ----------
    MDAData
        A list of MDA's to be done with date and coverage of the MDA included. 
        The coverage is given by the 4th value within each MDA of MDAData so we update the [3] position.
    coverage    
        The new coverage level for each MDA
    Returns
    -------
    function
        MDAData with updated coverage value
    """
    for MDA in MDAData:
        MDA[3] = coverage
    return MDAData

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
        parameters,
        sim_params,
        demog,
        MDA_times,
        MDAData,
        vacc_times,
        VaccData,
    ) = setup(initial_infect_frac)

    outputTimes = get_Intervention_times(
        getOutputTimes(range(2019, 2041)),
        START_DATE,
        sim_params['burnin'],
    )

    def do_single_run(seed, beta, coverage, i):
        np.random.seed(seed)
        init_vals = create_initial_population(initial_infect_frac, MDAData)
        random_state = np.random.get_state()
        return run_single_simulation(
            pickleData=copy.deepcopy(init_vals),
            params=parameters,
            timesim=sim_params["timesim"],
            burnin=sim_params["burnin"],
            demog=demog,
            beta=beta,
            MDA_times=MDA_times,
            MDAData=alterMDACoverage(MDAData, coverage),
            vacc_times=vacc_times,
            VaccData=VaccData,
            outputTimes=outputTimes,
            index=i,
            numpy_state=random_state,
        )

    def run_trachoma(seeds, params, n_tims):
        results: list[tuple[dict, list]]
        results = Parallel(n_jobs=num_cores)(
            delayed(do_single_run)(seed, beta=amisPars[0], coverage=amisPars[1], i=i)
            # params now will have 2 columns, one for beta and one for coverage
            # iterate over these, naming them amisPars, and amisPars[0] being beta
            # amisPars[1] being the coverage value
            for i, (seed, amisPars) in enumerate(zip(seeds, params))
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
        return np.transpose(
            np.transpose(
                [
                    # 'v' is a list, so convert it to NumPy array to index it from
                    # a list 'fitting_points'
                    np.asarray(v["True_Prev_Disease_children_1_9"])[fitting_points]
                    for (v, r) in results
                ]
            )
        )

    return run_trachoma
