params = {'N': 2500,
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
          'importation_rate': 0,
          'importation_reduction_rate': (0.9)**(1/10)}

sim_params = {'timesim':(52*96)+1,
              'burnin': (52*70)-1,
              'N_MDA': 0} # Set by main script later


demog = {'tau': 0.0004807692, 
         'max_age': 3120,
         'mean_age': 1040}
