from trachoma.trachoma_functions import *
import multiprocessing
from joblib import Parallel, delayed
import pickle
import numpy as np
import time
import os
import pandas as pd


folder_id = "source-data-20250130" # needs to be manually changed if multiple batches to be run due to not enough storage. if so, must run some batches, send to cloud, delete, repeat
species = "trachoma"
species_prefix = "Trachoma_"
IU_SLURM = os.getenv("SLURM_ARRAY_TASK_ID")
num_cores = 10


start = time.time()

initial_infect_frac = 0.1
parameters = {'N': 2500,
          'av_I_duration' : 2,
          'av_ID_duration':300/7,
          'inf_red':0.45,
          'min_ID':11, #Parameters relating to duration of infection period, including ID period
          'av_D_duration':200/7,
          'min_D':10.1/7, #Parameters relating to duration of disease period
          'dis_red':0.3,
          'v_1':1,
          'v_2':2.6,
          'phi':1.4,
          'epsilon':0.5,#Parameters relating to lambda function- calculating force of infection
          #Parameters relating to MDA
          'MDA_Cov':0.8,
          'MDA_Eff': 0.85, # Efficacy of treatment
          'rho':0.3,
          'nweeks_year':52,
          'babiesMaxAge':0.5, #Note this is years, need to check it converts to weeks later
          'youngChildMaxAge':9,#Note this is years, need to check it converts to weeks later
          'olderChildMaxAge':15, #Note this is years, need to check it converts to weeks later
          'b1':1,#this relates to bacterial load function
          'ep2':0.114,
          'n_inf_sev':38,
          'TestSensitivity': 0.96,
          'TestSpecificity': 0.965,
          'SecularTrendIndicator': 0,
          'SecularTrendYearlyBetaDecrease': 0.01,
          'vacc_prob_block_transmission':  0.5, 
          'vacc_reduce_bacterial_load': 0, 
          'vacc_reduce_duration': 0,
          'vacc_coverage': 0,  
          'vacc_waning_length': 52 * 5,
          'importation_rate': 1/(52*2500),
          #'importation_rate': 0.00008,
          'importation_reduction_rate': (0.9)**(1/10)}

sim_params = {'timesim':(52*100)-1,
              'burnin': (52*70)-1,
              'N_MDA': 0} # Set by main script later


demog = {'tau': 0.0004807692, 
         'max_age': 3120,
         'mean_age': 1040}



START_DATE = date(1996, 1, 1)

#id = os.getenv("SLURM_ARRAY_TASK_ID")
# for testing

#mda_filepath = 'endgame_inputs/InputMDA_MTP_' + str(id) + '.csv'
#mda_filepath = 'endgame_inputs/InputMDA_MTP_' + str(id) + '.csv'


'''
    functions from amis file to set up the simulations
'''

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


def setup():
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


'''
    function to generate the NTDMC results from the simulations
'''
def getResultsNTDMC(results, Start_date, burnin):
    '''
    Function to collate results for NTDMC
    '''

    
    for i in range(len(results)):
        d = copy.deepcopy(results[i][0])
        prevs = np.array(d['True_Prev_Disease_children_1_9']) 
        start = burnin # get prevalence from the end of the burnin onwards
        step = 52 # step forward 52 weeks
        chosenPrevs = prevs[start::step] 
        if i == 0:
           df = pd.DataFrame(0, range(len(chosenPrevs)), columns= range(len(results)+4))
           df = df.rename(columns={0: "Time", 1: "age_start", 2: "age_end", 3: "measure"}) 
           df.iloc[:, 0] = range(Start_date.year, Start_date.year + len(chosenPrevs))
           df.iloc[:, 1] = np.repeat(1, len(chosenPrevs))
           df.iloc[:, 2] = np.repeat(9, len(chosenPrevs))
           df.iloc[:, 3] = np.repeat("prevalence", len(chosenPrevs))
        df.iloc[:,i+4] = chosenPrevs
    for i in range(len(results)):
        df = df.rename(columns={i+4: "draw_"+ str(i)}) 
    return df

'''
    Define names of files and other information for the runs. 
    This is the section which needs the most editing for each IU
'''

pathCountry = '../Maps/table_iu_idx_trachoma.csv' # some path here
df_IU_country = pd.read_csv(pathCountry)
iu = df_IU_country['IU_ID'].values[int(IU_SLURM)]
country = df_IU_country['country'].values[int(IU_SLURM)]
    
print(country)
print(str(iu).zfill(5))


''' 
    change to be whatever the scenario we are running is
'''

print("Making projections for trachoma for IU " + str(iu).zfill(5) + "...")
mda_filepath = 'endgame_inputs/InputMDA_MTP_projections_' + str(iu) + '.csv'


'''
    file name for IU specific parameters
'''
#ParamFilePath = '~/Documents/trachoma/post_AMIS_analysis/InputPars_MTP_trachoma/InputPars_MTP_' + str(iu) + '.csv'
ParamFilePath = f'projections/{species}/' + str(folder_id) + f'/{country}/{country}' + str(iu).zfill(5) + f'/InputBeta_{species_prefix}{country}' + str(iu).zfill(5) + '.csv'
print(ParamFilePath)
amisparams = pd.read_csv(ParamFilePath)
amisparams.columns = [s.replace(' ', '') for s in amisparams.columns]

# define the lists of random seeds, R0 and k
seeds = amisparams.iloc[:, 0].tolist()
seeds=list(map(int, seeds))
betas = amisparams.iloc[:, 1].tolist()
coverages = amisparams.iloc[:, 2].tolist()

'''
    numSims should be set to 200
'''
numSims = 200

'''
    additional functions to run the simulations
'''

(
        parameters,
        sim_params,
        demog,
        MDA_times,
        MDAData,
        vacc_times,
        VaccData,
    ) = setup()


paramsToAlter = copy.deepcopy(parameters)


outputTimes = get_Intervention_times(
    getOutputTimes(range(1996,2022)),
    START_DATE,
    sim_params['burnin'],
)

def do_single_run(seed, beta, coverage, i):
    np.random.seed(seed)
    altered_mda_coverage = alterMDACoverage(MDAData, coverage)
    init_vals = create_initial_population(initial_infect_frac, altered_mda_coverage)
    random_state = np.random.get_state()
    return run_single_simulation(
        pickleData=copy.deepcopy(init_vals),
        params=parameters,
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
        distToUse = "Exponential",
    )

results = Parallel(n_jobs=num_cores)(
             delayed(do_single_run)(seeds[i], betas[i], coverages[i], i) for i in range(numSims)
         )

simData=[item[0] for item in results]

'''
    Make pickle files and NTDMC data outputs
'''

# want outputs like <ascaris-folder>/AGO/AGO02049/Asc_AGO02049.p
newOutputSimDataFilePath = f'projections/trachoma/' + str(folder_id) + f'/{country}/{country}' + str(iu).zfill(5) + f'/{species_prefix}{country}' + str(iu).zfill(5) + '.p'
#newOutputSimDataFilePath = "pickleETC.p"

# only save subset of data
subset_keys = ["IndI",
        "IndD",
        "No_Inf",
        "T_latent",
        "T_ID",
        "T_D",
        "Ind_latent",
        "Ind_ID_period_base",
        "Ind_D_period_base",
        "bact_load",
        "Age",
        "vaccinated",
        "time_since_vaccinated",
        "treatProbability",
        "MDA_coverage",
        "systematic_non_compliance",
        "ids"]

subset_sim_data = [{key: d[key] for key in subset_keys if key in d} for d in simData]

print("Pickle file name:")
print(newOutputSimDataFilePath)
pickle.dump( subset_sim_data, open( newOutputSimDataFilePath, 'wb' ) )


NTDMC =  getResultsNTDMC(results, START_DATE, sim_params["burnin"])

PrevDatasetFilePath = f'projections/trachoma/' + str(folder_id) + f'/{country}/{country}' + str(iu).zfill(5) + f'/PrevDataset_{species_prefix}{country}' + str(iu).zfill(5) + '.csv'
#PrevDatasetFilePath = "NTDMCETC.csv"
print("PrevDataset_species_iu.csv file name:")
print(PrevDatasetFilePath)
NTDMC.to_csv(PrevDatasetFilePath, index=False)

print('Finished projections for ' + str(species) + ' in IU ' + str(iu) + ".")

#################################################################################################################

end = time.time()
elapsed_time = end-start
print("elapsed time in seconds: " + str(elapsed_time) )

