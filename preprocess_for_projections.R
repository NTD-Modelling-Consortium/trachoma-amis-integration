library(AMISforInfectiousDiseases)
library(dplyr)
library(tidyr)

setwd("../")
load("Maps/trachoma_maps.rds")
load("Maps/trachoma_map_years.rds")


failed_ids = c() # this needs to be updated once we know which batches failed
ctd_ids = c()

# Early Feb: second runs
# from Igor's runs
#failed_ids = c(123,150,170,196,203,296,315,369,371,373,376,377,379,380,495,502,515,528,530,531,532,533,535,539,541)
# from Raiha's runs
#failed_ids = c(123,150,196,203,205,296,315,369,370,371,373,376,379,380,495,502,515,528,530,531,532,533,535,539,541)
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

InputPars_MTP_path_above <- paste0("post_AMIS_analysis/InputPars_MTP_",species,"/")
if (!dir.exists(InputPars_MTP_path_above)) {dir.create(InputPars_MTP_path_above)}

sampled_params_all = c()
for(id in id_vec){
  
  #### Load AMIS output
  load(paste0("AMIS_output/amis_output",id,".Rdata")) # loads amis_output
  
  iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
  ess = amis_output$ess
  #### Sample draws from the posterior
  num_sub_samples_posterior <- 200
  
  for (iu in iu_names) {
    if(ess[which(iu_names==iu)] >= 200){
      sampled_params <- sample_parameters(x = amis_output, n_samples = num_sub_samples_posterior, locations = which(iu_names==iu))
      
      # this puts all ius in the same folder
      file_name <- paste0(InputPars_MTP_path_above, paste0("InputPars_MTP_",iu,".csv"))
      write.csv(sampled_params, file=file_name, row.names = F)
      
      sampled_params_iu = cbind(IU_ID=iu,sampled_params)
      sampled_params_all = rbind(sampled_params_all,sampled_params_iu)
      
      
    }
  }
  
  if(id%%100==0){cat(paste0("id=",id, "; "))}
  
}

save(sampled_params_all,file= paste0(InputPars_MTP_path_above,"InputPars_MTP_allIUs.rds"))


cat(paste0("Produced posterior parameter samples \n"))
