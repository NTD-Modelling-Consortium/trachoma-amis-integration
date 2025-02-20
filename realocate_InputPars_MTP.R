#### Sample draws from the AMIS output

library("AMISforInfectiousDiseases")
library("dplyr")
setwd("../post_AMIS_analysis")
# species <- "haematobium"
# species <- "mansoni_high"
# species <- "mansoni_low"

set.seed(123)
folder_id = "source-data-20250220" # needs to be manually changed if multiple batches to be run due to not enough storage. if so, must run some batches, send to cloud, delete, repeat
yearschange_index = c(2000,2010,2020)-1926+1 # 1926 is data start date (1996) minus 70 year burn in. 
n_timschange = length(yearschange_index)


failed_ids = c() # this needs to be filled in once we know which batches failed
ctd_ids = c()

# Early Feb: second runs
# from Igor's runs
#failed_ids = c(123,150,170,196,203,296,315,369,371,373,376,377,379,380,495,502,515,528,530,531,532,533,535,539,541)
# from Raiha's runs
#failed_ids = c(123,150,196,203,205,296,315,369,370,371,373,376,379,380,495,502,515,528,530,531,532,533,535,539,541)

# 30 Jan:  from first runs
#failed_ids = c(29, # fitting didn't work properly
#              133,134, # randomly stopped?
#               123,153,196,296,315,322,374,375,376,377,378,379,380,381,495,502,515,527,530,532,533,539,541)
#ctd_ids = c(55,63,86,113,114,121,145,167,168,170,171,173,174,175,179,180,181,191,192,197,198,199,200,
 #203,205,206,208,262,276,298,326,327,362,363,370,371,372,412,435,441,446,447,471,503,505,507,508,510,528,531,535,549) # ids that were run for extra iterations
# failed ctd runs: 150,369

species_all <- c("trachoma")

for(species in species_all){
  
  # loading 'iu_task_lookup' just to get the country codes
  df_IU_country <- read.csv("../Maps/table_iu_idx_trachoma.csv")
  countries <- sort(unique(df_IU_country$country))
  
  path_species <- paste0("../ntd-model-trachoma/projections/",species,"/")
  if (!dir.exists(path_species)) {dir.create(path_species)}

  path_folder_id <- paste0("../ntd-model-trachoma/projections/",species,"/",folder_id,"/")
  if (!dir.exists(path_folder_id)) {dir.create(path_folder_id)}

  for(country in countries){
    path_country <- paste0("../ntd-model-trachoma/projections/",species,"/",folder_id,"/",country,"/")
    if (!dir.exists(path_country)) {dir.create(path_country)}
  }
  
}

for(species in species_all){
  species_prefix <- "Trachoma_"
 
  # loading 'iu_task_lookup' (batches-IUs look up table for the fitting)
  load("../Maps/iu_task_lookup.Rdata")
  num_batches <- max(iu_task_lookup$TaskID)
  
  cat(paste0("species: ",species, " \n"))
  cat(paste0("num_batches: ",num_batches, " \n"))
  
  for(id in setdiff(1:num_batches,failed_ids)){
   
    #### Load AMIS output
    if(!id %in% ctd_ids){
      load(paste0("../AMIS_output/amis_output",id,".Rdata")) # loads amis_output
      
    } else {
      load(paste0("../AMIS_output_ctd/amis_output",id,".Rdata")) # loads amis_output
    }
    
    iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
    ess = amis_output$ess
    num_samples <- 200
  
    for (iu in iu_names) {
      if(ess[which(iu_names==iu)] >= 200){
        wh <- which(df_IU_country$IU_ID==iu)
        if(length(wh)!=1){stop("iu must be found exactly once in df_IU_country")}
        country <- df_IU_country[wh, "country"]
        
        iu0 <- sprintf("%05d", as.integer(iu))
        path_iu <- paste0("../ntd-model-trachoma/projections/",species,"/",folder_id,"/",country,"/",country, iu0,"/")
        if (!dir.exists(path_iu)) {dir.create(path_iu)}
        
        file_name_old <- paste0("InputPars_MTP_trachoma/InputPars_MTP_",iu,".csv")
        sampled_params <- read.csv(file_name_old)
        
        sampled_params <- sampled_params[, c("seed", "beta_init", paste0("beta",yearschange_index),"eff_cov", "k_parameter")]
      

	sampled_params_full = matrix(NA,ncol=102,nrow=nrow(sampled_params))
	colnames(sampled_params_full) = c("seed",paste0("beta",1926:2024),"eff_cov", "k_parameter")
	sampled_params_full[,1] = sampled_params[,1]
        sampled_params_full[,101] = sampled_params[,6]       
       	sampled_params_full[,102] = sampled_params[,7]      

        n_reps_1 = yearschange_index[1] - 1
        sampled_params_full[,2:(n_reps_1+1)] = matrix(rep(sampled_params[,2],n_reps_1),ncol=n_reps_1)
        
        n_reps_2 = yearschange_index[2] - yearschange_index[1] 
        sampled_params_full[,(yearschange_index[1]+1):(yearschange_index[2])] = matrix(rep(sampled_params[,2]*sampled_params[,3],n_reps_2),ncol=n_reps_2)
        
        n_reps_3 = yearschange_index[3] - yearschange_index[2] 
        sampled_params_full[,(yearschange_index[2]+1):(yearschange_index[3])] = matrix(rep(sampled_params[,2]*sampled_params[,3]*sampled_params[,4],n_reps_3),ncol=n_reps_3)
        
        n_reps_4 = (2026-1926+1) - yearschange_index[3] 
        sampled_params_full[,(yearschange_index[3]+1):(2026-1926+1)] = matrix(rep(sampled_params[,2]*sampled_params[,3]*sampled_params[,4]*sampled_params[,5],n_reps_4),ncol=n_reps_4)
	
        file_name_new <- paste0(path_iu, paste0("InputBeta_",species_prefix, country, iu0,".csv"))
        write.csv(sampled_params_full, file=file_name_new, row.names = F)
      }
    }
    
    if(id%%100==0){cat(paste0("id=",id, "; "))}
  
  }
  
  cat(paste0("Samples realocated for all IUs in InputPars_MTP_",species, "/ \n"))

}
