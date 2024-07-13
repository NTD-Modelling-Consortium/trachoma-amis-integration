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


def setup(initial_prevalence: float)
    params = trachoma_params.params
    sim_params = trachoma_params.sim_params
    demog = trachoma_params.demog
    
    vals = Set_inits(params, demog, sim_params)
    ids = random.sample(
        range(popsize),
        k=int(initial_prevalence * sim_params['N'])
    )
    vals['IndI'][ids] = 1
    vals['T_latent'][ids] = vals['Ind_latent'][ids]
    vals['No_Inf'][ids] = 1

    MDA_times, MDAData = setup_mda("scen2c.csv")
    vacc_times, VaccData = setup_vaccine("scen2c.csv")
    sim_params['N_MDA'] = len(MDA_times)
    sim_params['N_Vaccines'] = len(vacc_times)

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
