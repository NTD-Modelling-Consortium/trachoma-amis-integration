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
from trachoma_params import params, sim_params, demog

START_DATE = date(2019, 1, 1)


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


def get_amis_wrapper(initial_prevalence=0.01, num_cores=-2):
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
            )
            for i, beta in params
        )
        return
