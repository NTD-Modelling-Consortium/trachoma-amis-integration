library(AMISforInfectiousDiseases)
library(dplyr)
library(tidyr)

setwd("../")
load("Maps/trachoma_maps.rds")
load("Maps/trachoma_map_years.rds")

# from Igor's runs
failed_ids = c(123,150,170,196,203,296,315,369,371,373,376,377,379,380,495,502,515,528,530,531,532,533,535,539,541)
ctd_ids = c()

# 30 Jan:  from previous runs
#failed_ids = c(29, # fitting didn't work properly
#	       133,134, # randomly stopped?
#               123,153,196,296,315,322,374,375,376,377,378,379,380,381,495,502,515,527,530,532,533,539,541)
#ctd_ids = c(55,63,86,113,114,121,145,167,168,170,171,173,174,175,179,180,181,191,192,197,198,199,200,
 #203,205,206,208,262,276,298,326,327,362,363,370,371,372,412,435,441,446,447,471,503,505,507,508,510,528,531,535,549) # ids that were run for extra iterations
# failed ctd runs: 150,369

# loading 'iu_task_lookup' (batches-IUs look up table for the fitting)
load("Maps/iu_task_lookup.Rdata")
id_vec = setdiff(1:max(iu_task_lookup$TaskID),c(failed_ids))
M_l <- 200 # just to simulate from normal to run vioplots
species="trachoma"

# load(paste0("../../trajectories/trajectories_",1,"_",species,".Rdata")) # load 'trajectories' of batch 1
# total_num_years <- ncol(trajectories)
# num_sampled_traj <- 500
# 
# year_indices <-  c(18L,29L,38L) # 2002, 2013, 2022 
# 
# map_years <- c(2002, 2013, 2022)  # years in the map samples amis fitted to
# all_years <- 2002:2022

InputPars_MTP_path_above <- paste0("post_AMIS_analysis/InputPars_MTP_",species,"/")
if (!dir.exists(InputPars_MTP_path_above)) {dir.create(InputPars_MTP_path_above)}

#sampled_params_all = c()
for(id in id_vec){
  
  #### Load AMIS output
  if(!id %in% ctd_ids){
    load(paste0("AMIS_output/amis_output",id,".Rdata")) # loads amis_output
    
  } else {
    load(paste0("AMIS_output_ctd/amis_output",id,".Rdata")) # loads amis_output
  }
  
  #iu_names <- as.character(unique(iu_task_lookup_updated$IU_ID[iu_task_lookup_updated$TaskID==id]))
  iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
  ess = amis_output$ess
  #### Sample draws from the posterior
  num_sub_samples_posterior <- 200
  # InputPars_MTP_path <- paste0("../InputPars_MTP/InputPars_MTP_batch_",id,"/"
  # if (!dir.exists(InputPars_MTP_path)) {dir.create(InputPars_MTP_path)}
  
  for (iu in iu_names) {
    if(ess[which(iu_names==iu)] >= 200){
      sampled_params <- sample_parameters(x = amis_output, n_samples = num_sub_samples_posterior, locations = which(iu_names==iu))
      
      # this puts all ius in the same folder
      file_name <- paste0(InputPars_MTP_path_above, paste0("InputPars_MTP_",iu,".csv"))
      write.csv(sampled_params, file=file_name, row.names = F)
      
      # sampled_params_iu = cbind(IU_ID=iu,sampled_params)
      # sampled_params_all = rbind(sampled_params_all,sampled_params_iu)
      
      
    }
  }
  
  if(id%%100==0){cat(paste0("id=",id, "; "))}
  
}

#save(sampled_params_all,file= paste0(InputPars_MTP_path_above,"InputPars_MTP_allIUs.rds"))
#load(paste0(InputPars_MTP_path_above,"InputPars_MTP_allIUs.rds"))


