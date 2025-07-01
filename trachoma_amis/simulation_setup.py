"""
Shared functions for trachoma simulations used by both AMIS integration and projections.
"""
from datetime import date
import copy
import numpy as np
from pathlib import Path
import os

from ntdmc_trachoma.trachoma_functions import (
    run_single_simulation,
    Set_inits,
    readPlatformData,
    getInterventionDates,
    get_Intervention_times,
    getOutputTimes,
)
from joblib import Parallel, delayed

# Global start date
START_DATE = date(1996, 1, 1)

# Common paths
PATH_TO_WORKING_DIR = Path(os.getenv("TRACHOMA_AMIS_DIR", ""))
PATH_TO_FITTING_PREP_ARTEFACTS = Path(
    os.getenv(
        "PATH_TO_FITTING_PREP_ARTEFACTS",
        PATH_TO_WORKING_DIR / "fitting-prep/artefacts",
    )
)
PATH_TO_PROJECTIONS_PREP_ARTEFACTS = Path(
    os.getenv(
        "PATH_TO_PROJECTIONS_PREP_ARTEFACTS",
        PATH_TO_WORKING_DIR / "projections-prep/artefacts",
    )
)
PATH_TO_PROJECTIONS_ARTEFACTS = Path(
    os.getenv(
        "PATH_TO_PROJECTIONS_ARTEFACTS", PATH_TO_WORKING_DIR / "projections/artefacts"
    )
)


def set_start_date(datestr: str):
    """Set the global START_DATE from a date string."""
    global START_DATE
    try:
        START_DATE = date.fromisoformat(datestr)
    except ValueError as e:
        msg = f"{e}.\n" "Valid dates formats are YYYY-MM-DD, YYYYMMDD or YYYY-WXX-D."
        raise ValueError(msg)


def get_start_date():
    """Get the current START_DATE value."""
    return START_DATE


def setup_mda(cov_filepath, burnin):
    """Set up MDA data and times from coverage file."""
    MDAData = readPlatformData(cov_filepath, "MDA")
    MDA_dates = getInterventionDates(MDAData)
    MDA_times = get_Intervention_times(MDA_dates, START_DATE, burnin)
    return MDA_times, MDAData


def setup_vaccine(cov_filepath, burnin):
    """Set up vaccine data and times from coverage file."""
    VaccData = readPlatformData(cov_filepath, "Vaccine")
    Vaccine_dates = getInterventionDates(VaccData)
    vacc_times = get_Intervention_times(Vaccine_dates, START_DATE, burnin)
    return vacc_times, VaccData


def setup(mda_filepath, sim_params, parameters, demog):
    """
    Set up simulation parameters, MDA and vaccine data.
    
    Parameters
    ----------
    mda_filepath : Path or str
        Path to MDA coverage file
    sim_params : dict
        Simulation parameters dictionary (will be modified in place)
    parameters : dict
        Model parameters
    demog : dict
        Demographic parameters
        
    Returns
    -------
    tuple
        (parameters, sim_params, demog, MDA_times, MDAData, vacc_times, VaccData)
    """
    MDA_times, MDAData = setup_mda(mda_filepath, sim_params["burnin"])
    vacc_times, VaccData = setup_vaccine(mda_filepath, sim_params["burnin"])
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


def create_initial_population(initial_prevalence: float, MDAData, parameters, demog, sim_params, distToUse="Expo"):
    """
    Create initial population with specified prevalence.
    
    Parameters
    ----------
    initial_prevalence : float
        Initial infection prevalence
    MDAData : list
        MDA data
    parameters : dict
        Model parameters
    demog : dict
        Demographic parameters
    sim_params : dict
        Simulation parameters
    distToUse : str
        Distribution to use (default "Expo")
        
    Returns
    -------
    dict
        Initial population values
    """
    vals = Set_inits(
        parameters,
        demog,
        sim_params,
        MDAData,
        np.random.get_state(),
        distToUse=distToUse,
    )
    ids = np.random.choice(
        range(parameters["N"]), int(initial_prevalence * parameters["N"]), replace=False
    )
    vals["IndI"][ids] = 1
    vals["T_latent"][ids] = vals["Ind_latent"][ids]
    vals["No_Inf"][ids] = 1
    return vals


def alterMDACoverage(MDAData, coverage):
    """
    Update the coverage of each MDA for a run to be a given value.
    
    Parameters
    ----------
    MDAData
        A list of MDA's to be done with date and coverage of the MDA included.
        The coverage is given by the 4th value within each MDA of MDAData so we update the [3] position.
    coverage
        The new coverage level for each MDA
        
    Returns
    -------
    list
        MDAData with updated coverage value
    """
    for MDA in MDAData:
        MDA[3] = coverage
    return MDAData


def create_single_run_function(
    MDAData, 
    initial_infect_frac, 
    parameters, 
    sim_params, 
    demog, 
    MDA_times, 
    vacc_times, 
    VaccData, 
    outputTimes, 
    distToUse="Expo",
    timeToStopImportation=-1
):
    """
    Create a function to run a single simulation with the given parameters.
    
    This returns a closure that captures the common parameters and returns
    a function suitable for parallel execution.
    
    Parameters
    ----------
    MDAData : list
        MDA data
    initial_infect_frac : float
        Initial infection fraction
    parameters : dict
        Model parameters
    sim_params : dict
        Simulation parameters
    demog : dict
        Demographic parameters
    MDA_times : list
        MDA intervention times
    vacc_times : list
        Vaccine intervention times
    VaccData : list
        Vaccine data
    outputTimes : list
        Times to output results
    distToUse : str
        Distribution to use (default "Expo")
    timeToStopImportation : int
        Time to stop importation (default -1, disabled)
        
    Returns
    -------
    function
        Function to run a single simulation
    """
    def do_single_run(seed, beta, coverage, k_parameter, i):
        np.random.seed(seed)
        altered_mda_coverage = alterMDACoverage(MDAData, coverage)
        init_vals = create_initial_population(
            initial_infect_frac, altered_mda_coverage, parameters, demog, sim_params, distToUse
        )
        random_state = np.random.get_state()
        parameters_copy = copy.deepcopy(parameters)
        parameters_copy["infection_risk_shape"] = k_parameter
        return run_single_simulation(
            pickleData=copy.deepcopy(init_vals),
            params=parameters_copy,
            timesim=sim_params["timesim"],
            burnin=sim_params["burnin"],
            demog=demog,
            beta=beta,
            MDA_times=MDA_times,
            MDAData=altered_mda_coverage,
            vacc_times=vacc_times,
            VaccData=VaccData,
            outputTimes=outputTimes,
            index=i,
            numpy_state=random_state,
            doIHMEOutput=False,
            doSurvey=False,
            distToUse=distToUse,
            postMDAImportationReduction=True,
            timeToStopImportation=timeToStopImportation,
        )
    
    return do_single_run


def run_trachoma_simulations(
    seeds, 
    betas, 
    coverages, 
    k_parameters,
    mda_filepath,
    parameters,
    sim_params,
    demog,
    initial_infect_frac=0.1,
    num_cores=-2,
    timeToStopImportation=-1,
    distToUse="Expo"
):
    """
    Core function for running trachoma simulations with given parameters.
    
    Used by both AMIS fitting and projections scripts to eliminate code duplication.
    
    Parameters
    ----------
    seeds : list
        Random seeds for each simulation
    betas : list or array
        Beta transmission parameters for each simulation
    coverages : list
        MDA coverage values for each simulation  
    k_parameters : list
        Infection risk shape parameters for each simulation
    mda_filepath : Path or str
        Path to MDA coverage file
    parameters : dict
        Model parameters
    sim_params : dict
        Simulation parameters
    demog : dict
        Demographic parameters
    initial_infect_frac : float
        Initial infection fraction (default 0.1)
    num_cores : int
        Number of cores for parallel processing (default -2)
    timeToStopImportation : int
        Time to stop importation, -1 to disable (default -1)
    distToUse : str
        Distribution to use (default "Expo")
        
    Returns
    -------
    list
        Raw simulation results from run_single_simulation
    """
    # Common setup logic
    (
        parameters_updated,
        sim_params_updated,
        demog_updated,
        MDA_times,
        MDAData,
        vacc_times,
        VaccData,
    ) = setup(mda_filepath, sim_params, parameters, demog)
    
    # Common output times setup
    outputTimes = get_Intervention_times(
        getOutputTimes(range(1996, 2022)),
        START_DATE,
        sim_params_updated["burnin"],
    )
    
    # Create the single run function
    do_single_run = create_single_run_function(
        MDAData, 
        initial_infect_frac, 
        parameters_updated, 
        sim_params_updated, 
        demog_updated, 
        MDA_times, 
        vacc_times, 
        VaccData, 
        outputTimes, 
        distToUse,
        timeToStopImportation
    )
    
    # Common parallel execution
    results = Parallel(n_jobs=num_cores)(
        delayed(do_single_run)(
            seeds[i],
            betas[i],
            coverages[i],
            k_parameters[i],
            i,
        )
        for i in range(len(seeds))
    )
    
    return results