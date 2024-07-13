from joblib import Parallel, delayed
import numpy as np

from trachoma.trachoma_functions import(
    run_single_simulation,
    Set_inits,
    readPlatformData,
    getInterventionDates,
    get_Intervention_times,
    getOutputTimes,
)

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
    

def get_amis_wrapper(initial_prevalence = 0.01, num_cores=-2):
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

    def run_trachoma(seeds, params, n_tims):
        results, _ = Parallel(n_jobs=num_cores)(
            delayed(run_single_simulation)(
                pickleData = copy.deepcopy(init_vals), 
                params = params, 
                timesim = sim_params['timesim'],
                burnin = sim_params['burnin'],
                demog=demog, 
                beta = beta[0],
                MDA_times = MDA_times, 
                MDAData=MDAData, 
                vacc_times = vacc_times, 
                VaccData = VaccData,
                outputTimes= outputTimes, 
                index = i,
            )
            for i, beta in params
        )
        return 
