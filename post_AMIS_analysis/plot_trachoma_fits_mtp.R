library(truncnorm)
library(invgamma)
library(ggplot2)
library(gridExtra)
library(magrittr)
library(dplyr)
library(tidyr)
library(reshape2)
library(patchwork)
library(AMISforInfectiousDiseases)
library(sf)
library(viridis)

setwd("~/Documents/trachoma-endgame")
load("Maps/trachoma_maps.rds")
load("Maps/trachoma_map_years.rds")
country_iu_lookup = read.csv("Maps/table_iu_idx_trachoma.csv")

failed_ids = c(55,69,113,114,123,150,153,170,180,189,276,296,369,375,376,379,380,381,495,503,510,528,530,532,533,541)
ctd_ids = c()
id_vec = setdiff(1:564,c(failed_ids))

n_samples = 1000 
saveFigures = T 
savePars = T

# define prior hyperparameters and boundaries
beta_rate = 5
beta_lb = 0
beta_ub = Inf

eff_cov_lb = 0
eff_cov_ub = 1
a_eff_cov = 5
b_eff_cov = 2 

k_mean = 6.4 # based on data from Gambia
k_sd = 3
k_lb = 0.01
k_ub = 20 


if(savePars==T){
  
  # save ESS and append all samples
  ess_gt200 = c()
  prev_colnames_all =  paste0("prev_",years_vector)
  colnames_all = c("IU_ID","seed","beta_init","beta75","beta85","beta95","eff_cov","k_parameter",prev_colnames_all)
  sampled_params_all = matrix(ncol=33,nrow=0,dimnames =list(NULL,colnames_all))
  
  for (id in c(id_vec)){
    # load amis output
    load(paste0("AMIS_output/amis_output",id,".Rdata")) # loads amis_output
    
    #load map
    prevalence_map = trachoma_maps
    prevalence_map = lapply(1:length(prevalence_map), function(t){
      output=list(data = as.matrix(prevalence_map[[t]]$data %>%
                                     filter(TaskID==id) %>%
                                     select(max_tf_prev_upr),ncol=1),
                  likelihood = prevalence_map[[t]]$likelihood)
      rownames(output$data) = prevalence_map[[t]]$data$IU_ID[prevalence_map[[t]]$data$TaskID==id]
      return(output)
    })
    # just get maps with data for this batch
    map_bool_id = which(sapply(1:length(prevalence_map), function(t) !all(is.na(prevalence_map[[t]]$data))))
    years_vector_id = years_vector[map_bool_id]
    prevalence_map = lapply(map_bool_id, function(t) prevalence_map[[t]])
    
    #save ess
    iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
    ess = amis_output$ess
    inds_to_include = which(ess >= 200)
    if(length(inds_to_include)>0){
      ess_gt200 = rbind(ess_gt200, cbind(IU_ID = iu_names[inds_to_include], TaskID=id, ess=ess[inds_to_include]))
    }
    
    #if (!id %in% failed_ids){
    for (iu in iu_names[inds_to_include]){
      
      country_code = country_iu_lookup$country[which(country_iu_lookup$IU_ID == iu)]
      iu0 <- sprintf("%05d", as.integer(iu))
      params_to_use <- read.csv(paste0("post_AMIS_analysis/trachoma-source-data-20250225-espen-input-bet-csv/InputBet_",country_code,iu0,".csv"))
      sampled_params = data.frame(IU_ID=iu,
                                  seed=params_to_use$seed,
                                  amis_output$param[params_to_use$seed,],
                                  amis_output$simulated_prevalences[params_to_use$seed,])
      
      #file_name <- paste0("post_AMIS_analysis/InputPars_MTP_trachoma/InputPars_MTP_",iu,".csv")
      #sampled_params = cbind(IU_ID=iu,read.csv(file_name))
      
      #rename columns
      prev_colnames = paste0("prev_",years_vector_id)
      colnames(sampled_params) = c(colnames(sampled_params[1:8]),prev_colnames)
      #add in missing years
      sampled_params[,setdiff(colnames(sampled_params_all), colnames(sampled_params))] <- NA
      sampled_params = sampled_params[,colnames_all]
      #append
      sampled_params_all = rbind(sampled_params_all,sampled_params)
    }
    #}
  }
  save(sampled_params_all,file="sampled_params_all.Rdata")
  
}


# plot trajectories and parameter distributions for good batches and IUs with ess>200 
if (saveFigures==T){
  for (id in c(id_vec)){
  
    load(paste0("AMIS_output/amis_output",id,".Rdata")) # loads amis_output
    load(paste0("trajectories/trajectories_",id,".Rdata"))
    load(paste0("infections/infections",id,".Rdata"))
      
    mda_history = read.csv(paste0("ntd-model-trachoma/trachoma/data/coverage/endgame_inputs/InputMDA_MTP_",id,".csv")) %>%
      select(starts_with("X")) %>% 
      pivot_longer(cols=starts_with("X"), names_to="Year", values_to="ind") %>%
      mutate(Year = as.numeric(substr(Year,2,5))) %>%
      filter(ind==1)
    
    prevalence_map = trachoma_maps
    prevalence_map = lapply(1:length(prevalence_map), function(t){
      output=list(data = as.matrix(prevalence_map[[t]]$data %>% 
                                     filter(TaskID==id) %>% 
                                     select(max_tf_prev_upr),ncol=1), 
                  likelihood = prevalence_map[[t]]$likelihood)
      rownames(output$data) = prevalence_map[[t]]$data$IU_ID[prevalence_map[[t]]$data$TaskID==id]
      return(output)
    })
    
    # just get maps with data for this batch
    map_bool_id = which(sapply(1:length(prevalence_map), function(t) !all(is.na(prevalence_map[[t]]$data))))
    years_vector_id = years_vector[map_bool_id]
    prevalence_map = lapply(map_bool_id, function(t) {
      data = as.data.frame(prevalence_map[[t]]$data) %>%
        mutate(tf_prev_lwr = case_when(
          max_tf_prev_upr == 0.05 ~ 0,
          max_tf_prev_upr == 0.099 ~ 0.05,
          max_tf_prev_upr == 0.299 ~ 0.1,
          max_tf_prev_upr == 0.499 ~ 0.3,
          max_tf_prev_upr == 1 ~ 0.5))
      return(list(data=data,likelihood=prevalence_map[[t]]$likelihood))
    })
    all_years_vector_id = seq(1996,2022,by=0.25)
    
    # prepare data for plot
    cov_values = as.data.frame(amis_output$param) %>%
      mutate(Var1 = row_number()) %>%
      select(Var1,eff_cov)
    
    trajectories_melt = melt(trajectories) %>%
      left_join(cov_values,by="Var1")
    colnames(trajectories_melt) =c("id","time","prevalence","eff_cov")
    
    times_lookup = data.frame(time=unique(trajectories_melt$time),year=all_years_vector_id)
    trajectories_melt = trajectories_melt %>%
      left_join(times_lookup,by=c("time"))
    
    # same for infections
    infections_melt = melt(infections) %>%
      left_join(cov_values,by="Var1")
    colnames(infections_melt) =c("id","time","prevalence","eff_cov")

    infections_melt = infections_melt %>%
      left_join(times_lookup,by="time")
    
    # draws from prior  
    L <- nrow(prevalence_map[[1]]$data)
    panel_ncols <- 4
    panel_nrows <- ceiling(L/panel_ncols)
    
    iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
    
    trajectories_plots = list()
    infections_plots = list()
    for (j in 1:length(rownames(prevalence_map[[1]]$data))){
      iu = iu_names[j]
      samples = sample(1:n_samples,100,replace=F)
      # plot trajectories
      p = ggplot(trajectories_melt[which(trajectories_melt$id %in% samples),], 
                 aes(x = year, y = prevalence, group=id)) + 
        geom_line(aes(alpha = eff_cov), lwd=0.3) + 
        theme_bw() + 
        ggtitle(paste0(iu_names[j])) + 
        scale_y_continuous(limits = c(0,1)) + 
        theme(legend.position="none")
      for (year_ind in 1:length(prevalence_map)){
        p = p + geom_segment(y=prevalence_map[[year_ind]]$data[j,2],yend=prevalence_map[[year_ind]]$data[j,1],x=(years_vector_id[year_ind]+0.75),col="red") 
      }
      for (mda_ind in 1:nrow(mda_history)){
        p = p + geom_vline(xintercept=as.numeric(mda_history$Year[mda_ind],alpha=0.5),col="green")
      }
      trajectories_plots[[j]] = p
      
      # plot infections
      p = ggplot(infections_melt[which(infections_melt$id %in% samples),],
                 aes(x = year, y = prevalence, group=id)) +
        geom_line(aes(alpha = eff_cov), lwd=0.3) +
        theme_bw() +
        ggtitle(paste0(iu_names[j])) +
        scale_y_continuous(limits = c(0,1)) +
        theme(legend.position="none")
      for (mda_ind in 1:nrow(mda_history)){
        p = p + geom_vline(xintercept=as.numeric(mda_history$Year[mda_ind],alpha=0.5),col="green")
      }
      infections_plots[[j]] = p
    }
    
    png(paste0("post_AMIS_analysis/plots/trajectories_prior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
    print(patchwork::wrap_plots(trajectories_plots,nrow = panel_nrows, ncol = panel_ncols))
    dev.off()
    
  
    png(paste0("post_AMIS_analysis/plots/infections_prior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
    print(patchwork::wrap_plots(infections_plots,nrow = panel_nrows, ncol = panel_ncols))
    dev.off()
 
    #### Sample draws from the posterior
    num_sub_samples_posterior <- 200
    ess_all = amis_output$ess
    inds_to_exclude = which(ess_all < 200)
    inds_to_include = setdiff(1:length(rownames(prevalence_map[[1]]$data)),inds_to_exclude)
    
    # draws from posterior
    trajectories_plots = list()
    infections_plots = list()
    #parameters_plots = list()
    parameters_betainit_plots = list()
    parameters_beta2000_plots = list()
    parameters_beta2010_plots = list()
    parameters_beta2020_plots = list()
    parameters_k_plots = list()
    parameters_effcov_plots = list()
    
    if(length(inds_to_include)>0){ # make sure number of IUs with ESS >= 200 is nonzero
      for (j in inds_to_include){
        list_ind = which(inds_to_include==j)
        iu = iu_names[j]
        ess = ess_all[j]
        
        
        country_code = country_iu_lookup$country[which(country_iu_lookup$IU_ID == iu)]
        iu0 <- sprintf("%05d", as.integer(iu))
        
        sampled_params_iu <- read.csv(paste0("post_AMIS_analysis/trachoma-source-data-20250225-espen-input-bet-csv/InputBet_",country_code,iu0,".csv"))

        #sampled_params_iu = read.csv(paste0("post_AMIS_analysis/InputPars_MTP_trachoma/InputPars_MTP_",iu,".csv"))
        
        p = ggplot() +
          geom_density(data=sampled_params_iu,aes(x=beta1996)) +
          geom_density(aes(x=rexp(10000,rate=beta_rate)),linetype="dotted") +
          theme_bw() + 
          ggtitle(paste0(iu_names[j], " (ESS:",round(ess,digits=0),")")) +
          xlim(c(0,0.4))
        
        parameters_betainit_plots[[list_ind]] = p
        
        
        p = ggplot(data=sampled_params_iu,aes(x=(beta2000))) +
          geom_density() +
          theme_bw() + 
          ggtitle(paste0(iu_names[j], " (ESS:",round(ess,digits=0),")")) +
          xlim(c(0,0.4)) +
          xlab("beta_firstchange")
        
        
        parameters_beta2000_plots[[list_ind]] = p
        
        
        p = ggplot(data=sampled_params_iu,aes(x=(beta2010))) +
          geom_density() +
          theme_bw() + 
          ggtitle(paste0(iu_names[j], " (ESS:",round(ess,digits=0),")")) +
          xlim(c(0,0.4)) +
          xlab("beta_secondchange")
        
        parameters_beta2010_plots[[list_ind]] = p
        
        p = ggplot(data=sampled_params_iu,aes(x=(beta2020))) +
          geom_density() +
          theme_bw() + 
          ggtitle(paste0(iu_names[j], " (ESS:",round(ess,digits=0),")")) +
          xlim(c(0,0.4)) +
          xlab("beta_thirdchange")
        
        parameters_beta2020_plots[[list_ind]] = p
        
        p = ggplot() +
          geom_density(data=sampled_params_iu,aes(x=eff_cov)) +
          geom_density(aes(x=rbeta(10000,a_eff_cov,b_eff_cov)),linetype="dotted") +
          theme_bw() + 
          ggtitle(paste0(iu_names[j], " (ESS:",round(ess,digits=0),")")) +
          xlim(c(0,1))
        
        parameters_effcov_plots[[list_ind]] = p
        
        p = ggplot() +
          geom_density(data=sampled_params_iu,aes(x=k_parameter)) +
          geom_density(aes(x=rtruncnorm(10000,a=k_lb, b=k_ub, mean=k_mean, sd=k_sd)),linetype="dotted") +
          theme_bw() + 
          ggtitle(paste0(iu_names[j], " (ESS:",round(ess,digits=0),")")) +
          xlim(c(0,20)) 
        
        parameters_k_plots[[list_ind]] = p
        
        # plot trajectories
        p = ggplot(trajectories_melt[which(trajectories_melt$id %in% sampled_params_iu$seed),], 
                   aes(x = year, y = prevalence, group=id)) + 
          geom_line(aes(alpha = eff_cov), lwd=0.3) + 
          theme_bw() + 
          ggtitle(paste0(iu_names[j], " (ESS:",round(ess,digits=0),")")) + 
          scale_y_continuous(limits = c(0,1)) + 
          theme(legend.position="none")
        for (year_ind in 1:length(prevalence_map)){
          p = p + geom_segment(y=prevalence_map[[year_ind]]$data[j,2],yend=prevalence_map[[year_ind]]$data[j,1],x=(years_vector_id[year_ind]+0.75),col="red") 
        }
        for (mda_ind in 1:nrow(mda_history)){
          p = p + geom_vline(xintercept=as.numeric(mda_history$Year[mda_ind],alpha=0.5),col="green")
        }
        trajectories_plots[[list_ind]] = p
        
        # plot infections
        p = ggplot(infections_melt[which(infections_melt$id %in% sampled_params_iu$seed),],
                   aes(x = year, y = prevalence, group=id)) +
          geom_line(aes(alpha = eff_cov), lwd=0.3) +
          theme_bw() +
          ggtitle(paste0(iu_names[j], " (ESS:",round(ess,digits=0),")")) +
          scale_y_continuous(limits = c(0,1)) +
          theme(legend.position="none")

        for (mda_ind in 1:nrow(mda_history)){
          p = p + geom_vline(xintercept=as.numeric(mda_history$Year[mda_ind],alpha=0.5),col="green")
        }
        infections_plots[[list_ind]] = p
      }
      
      png(paste0("post_AMIS_analysis/plots/trajectories_posterior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
      print(patchwork::wrap_plots(trajectories_plots,nrow = panel_nrows, ncol = panel_ncols))
      dev.off()
      
      png(paste0("post_AMIS_analysis/plots/infections_posterior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
      print(patchwork::wrap_plots(infections_plots,nrow = panel_nrows, ncol = panel_ncols))
      dev.off()
      
      png(paste0("post_AMIS_analysis/plots/betainit_posterior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
      print(patchwork::wrap_plots(parameters_betainit_plots,nrow = panel_nrows, ncol = panel_ncols))
      dev.off()
      
      png(paste0("post_AMIS_analysis/plots/beta2000_posterior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
      print(patchwork::wrap_plots(parameters_beta2000_plots,nrow = panel_nrows, ncol = panel_ncols))
      dev.off()
      
      png(paste0("post_AMIS_analysis/plots/beta2010_posterior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
      print(patchwork::wrap_plots(parameters_beta2010_plots,nrow = panel_nrows, ncol = panel_ncols))
      dev.off()
      
      png(paste0("post_AMIS_analysis/plots/beta2020_posterior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
      print(patchwork::wrap_plots(parameters_beta2020_plots,nrow = panel_nrows, ncol = panel_ncols))
      dev.off()
      
      png(paste0("post_AMIS_analysis/plots/effcov_posterior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
      print(patchwork::wrap_plots(parameters_effcov_plots,nrow = panel_nrows, ncol = panel_ncols))
      dev.off()
      
      png(paste0("post_AMIS_analysis/plots/k_posterior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
      print(patchwork::wrap_plots(parameters_k_plots,nrow = panel_nrows, ncol = panel_ncols))
      dev.off()
      
    }
  } 
}

# plot initial trajectories for failed batches
if (saveFigures==T){
  for (id in c(failed_ids)){
    
    load(paste0("trajectories/trajectories_",id,".Rdata"))
    load(paste0("infections/infections",id,".Rdata"))
    
    mda_history = read.csv(paste0("ntd-model-trachoma/trachoma/data/coverage/endgame_inputs/InputMDA_MTP_",id,".csv")) %>%
      select(starts_with("X")) %>% 
      pivot_longer(cols=starts_with("X"), names_to="Year", values_to="ind") %>%
      mutate(Year = as.numeric(substr(Year,2,5))) %>%
      filter(ind==1)
    
    prevalence_map = trachoma_maps
    prevalence_map = lapply(1:length(prevalence_map), function(t){
      output=list(data = as.matrix(prevalence_map[[t]]$data %>% 
                                     filter(TaskID==id) %>% 
                                     select(max_tf_prev_upr),ncol=1), 
                  likelihood = prevalence_map[[t]]$likelihood)
      rownames(output$data) = prevalence_map[[t]]$data$IU_ID[prevalence_map[[t]]$data$TaskID==id]
      return(output)
    })
    
    # just get maps with data for this batch
    map_bool_id = which(sapply(1:length(prevalence_map), function(t) !all(is.na(prevalence_map[[t]]$data))))
    years_vector_id = years_vector[map_bool_id]
    prevalence_map = lapply(map_bool_id, function(t) {
      data = as.data.frame(prevalence_map[[t]]$data) %>%
        mutate(tf_prev_lwr = case_when(
          max_tf_prev_upr == 0.05 ~ 0,
          max_tf_prev_upr == 0.099 ~ 0.05,
          max_tf_prev_upr == 0.299 ~ 0.1,
          max_tf_prev_upr == 0.499 ~ 0.3,
          max_tf_prev_upr == 1 ~ 0.5))
      return(list(data=data,likelihood=prevalence_map[[t]]$likelihood))
    })
    all_years_vector_id = seq(1996,2022,by=0.25)
    
    
    trajectories_melt = melt(trajectories) 
    colnames(trajectories_melt) =c("id","time","prevalence")
    
    times_lookup = data.frame(time=unique(trajectories_melt$time),year=all_years_vector_id)
    trajectories_melt = trajectories_melt %>%
      left_join(times_lookup,by=c("time"))
    
    # same for infections
    infections_melt = melt(infections) 
    colnames(infections_melt) =c("id","time","prevalence")
    
    infections_melt = infections_melt %>%
      left_join(times_lookup,by="time")
    
    # draws from prior  
    L <- nrow(prevalence_map[[1]]$data)
    panel_ncols <- 4
    panel_nrows <- ceiling(L/panel_ncols)
    
    iu_names <- rownames(prevalence_map[[1]]$data)
    
    trajectories_plots = list()
    infections_plots = list()
    for (j in 1:length(rownames(prevalence_map[[1]]$data))){
      iu = iu_names[j]
      samples = sample(1:n_samples,100,replace=F)
      # plot trajectories
      p = ggplot(trajectories_melt[which(trajectories_melt$id %in% samples),], 
                 aes(x = year, y = prevalence, group=id)) + 
        geom_line(lwd=0.3) + 
        theme_bw() + 
        ggtitle(paste0(iu_names[j])) + 
        scale_y_continuous(limits = c(0,1)) + 
        theme(legend.position="none")
      for (year_ind in 1:length(prevalence_map)){
        p = p + geom_segment(y=prevalence_map[[year_ind]]$data[j,2],yend=prevalence_map[[year_ind]]$data[j,1],x=(years_vector_id[year_ind]+0.75),col="red") 
      }
      for (mda_ind in 1:nrow(mda_history)){
        p = p + geom_vline(xintercept=as.numeric(mda_history$Year[mda_ind],alpha=0.5),col="green")
      }
      trajectories_plots[[j]] = p
      
      # plot infections
      p = ggplot(infections_melt[which(infections_melt$id %in% samples),],
                 aes(x = year, y = prevalence, group=id)) +
        geom_line(lwd=0.3) +
        theme_bw() +
        ggtitle(paste0(iu_names[j])) +
        scale_y_continuous(limits = c(0,1)) +
        theme(legend.position="none")
      for (mda_ind in 1:nrow(mda_history)){
        p = p + geom_vline(xintercept=as.numeric(mda_history$Year[mda_ind],alpha=0.5),col="green")
      }
      infections_plots[[j]] = p
    }
    
    png(paste0("post_AMIS_analysis/plots/failed_trajectories_prior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
    print(patchwork::wrap_plots(trajectories_plots,nrow = panel_nrows, ncol = panel_ncols))
    dev.off()
    
    
    png(paste0("post_AMIS_analysis/plots/failed_infections_prior_",id,".png"),width=10,height=panel_nrows*2.5,units="in",res=900)
    print(patchwork::wrap_plots(infections_plots,nrow = panel_nrows, ncol = panel_ncols))
    dev.off()
  } 
}

# trajectories and infections and betas with all IUs in one PDF
if (saveFigures==T){

  
  #### Sample draws from the posterior
  L <- 2789 # total number of IUs
  panel_ncols <- 4
  panel_nrows <- ceiling(L/panel_ncols)
  num_sub_samples_posterior <- 200
  
  # draws from posterior
  list_counter = 0
  trajectories_plots = list()
  infections_plots = list()
  beta_baseline_first = list()
  beta_baseline_second = list()
  beta_baseline_third = list()
  betadrop_coverage = list()
  iu_order = data.frame(list_ind=rep(NA,L),iu_code=rep(NA,L))
  
  for (id in c(id_vec)){

    load(paste0("AMIS_output/amis_output",id,".Rdata")) # loads amis_output
    load(paste0("trajectories/trajectories_",id,".Rdata"))
    load(paste0("infections/infections",id,".Rdata"))
    
    mda_history = read.csv(paste0("ntd-model-trachoma/trachoma/data/coverage/endgame_inputs/InputMDA_MTP_",id,".csv")) %>%
      select(starts_with("X")) %>% 
      pivot_longer(cols=starts_with("X"), names_to="Year", values_to="ind") %>%
      mutate(Year = as.numeric(substr(Year,2,5))) %>%
      filter(ind==1)
    
    
    prevalence_map = trachoma_maps
    prevalence_map = lapply(1:length(prevalence_map), function(t){
      output=list(data = as.matrix(prevalence_map[[t]]$data %>% 
                                     filter(TaskID==id) %>% 
                                     select(max_tf_prev_upr),ncol=1), 
                  likelihood = prevalence_map[[t]]$likelihood)
      rownames(output$data) = prevalence_map[[t]]$data$IU_ID[prevalence_map[[t]]$data$TaskID==id]
      return(output)
    })
    
    # just get maps with data for this batch
    map_bool_id = which(sapply(1:length(prevalence_map), function(t) !all(is.na(prevalence_map[[t]]$data))))
    years_vector_id = years_vector[map_bool_id]
    prevalence_map = lapply(map_bool_id, function(t) {
      data = as.data.frame(prevalence_map[[t]]$data) %>%
        mutate(tf_prev_lwr = case_when(
          max_tf_prev_upr == 0.05 ~ 0,
          max_tf_prev_upr == 0.099 ~ 0.05,
          max_tf_prev_upr == 0.299 ~ 0.1,
          max_tf_prev_upr == 0.499 ~ 0.3,
          max_tf_prev_upr == 1 ~ 0.5))
      return(list(data=data,likelihood=prevalence_map[[t]]$likelihood))
    })
    #all_years_vector_id = seq(min(years_vector_id),2021,by=0.25)
    all_years_vector_id = seq(1996,2022,by=0.25)
    
    # prepare data for plot
    cov_values = as.data.frame(amis_output$param) %>%
      mutate(Var1 = row_number()) %>%
      select(Var1,eff_cov)
    
    trajectories_melt = melt(trajectories) %>%
      left_join(cov_values,by="Var1")
    colnames(trajectories_melt) =c("id","time","prevalence","eff_cov")
    
    times_lookup = data.frame(time=unique(trajectories_melt$time),year=all_years_vector_id)
    trajectories_melt = trajectories_melt %>%
      left_join(times_lookup,by=c("time"))
    
    # same for infections
    infections_melt = melt(infections) %>%
      left_join(cov_values,by="Var1")
    colnames(infections_melt) =c("id","time","prevalence","eff_cov")
    
    infections_melt = infections_melt %>%
      left_join(times_lookup,by="time")
    
    iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
    ess_all = amis_output$ess
    inds_to_exclude = which(ess_all < 200)
    inds_to_include = setdiff(1:length(rownames(prevalence_map[[1]]$data)),inds_to_exclude)
    
    if(length(inds_to_include)>0){ # make sure number of IUs with ESS >= 200 is nonzero
      
      # get list indices for this batch
      plot_indices = list_counter + (1:length(inds_to_include))
      
      for (j in inds_to_include){
        list_ind = plot_indices[which(inds_to_include==j)]
        iu = iu_names[j]
        ess = ess_all[j]
        
        country_code = country_iu_lookup$country[which(country_iu_lookup$IU_ID == iu)]
        iu0 <- sprintf("%05d", as.integer(iu))
        iu_order[list_ind,] = c(list_ind,paste0(country_code,iu0))
        
        sampled_params_iu <- read.csv(paste0("post_AMIS_analysis/trachoma-source-data-20250225-espen-input-bet-csv/InputBet_",country_code,iu0,".csv"))
        
        # prep beta data for plot
        beta_vals = sampled_params_iu %>%
          select(paste0("beta",1996:2024))
        beta_summary = apply(beta_vals, 2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))
        colnames(beta_summary) = paste0(1996:2024)
        beta_summary %<>% melt() %>% pivot_wider(names_from="Var1",values_from="value")
        colnames(beta_summary) = c("Year","quant_0.025","quant_0.5","quant_0.975")
        beta_summary %<>% filter(Year <=2022)
        
        # beta reductions
        beta_reductions = sampled_params_iu %>%
          select(beta1996,beta2000,beta2010,beta2020,eff_cov) %>%
          mutate(reduction2000 = (beta2000-beta1996)/beta1996,
                 reduction2010 = (beta2010-beta1996)/beta1996,
                 reduction2020 = (beta2020-beta1996)/beta1996) %>%
          select(beta1996,reduction2000,reduction2010,reduction2020,eff_cov)
        
        # plot betas
        p = ggplot(data=beta_reductions,
                   aes(x = beta1996, y = reduction2000)) +
          geom_point(aes(alpha=eff_cov)) +
          theme_bw() +
          ggtitle(paste0(country_code,iu0," (ESS:",round(ess,digits=0),")")) +
          scale_y_continuous(limits = c(-1,0)) +
          scale_x_continuous(limits = c(0.01,2.5),trans = "log10") +
          theme(legend.position="none")
        beta_baseline_first[[list_ind]] = p
        
        p = ggplot(data=beta_reductions,
                   aes(x = beta1996, y = reduction2010)) +
          geom_point(aes(alpha=eff_cov)) +
          theme_bw() +
          ggtitle(paste0(country_code,iu0," (ESS:",round(ess,digits=0),")")) +
          scale_y_continuous(limits = c(-1,0)) +
          scale_x_continuous(limits = c(0.01,2.5),trans = "log10") +
          theme(legend.position="none")
        beta_baseline_second[[list_ind]] = p
        
        p = ggplot(data=beta_reductions,
                   aes(x = beta1996, y = reduction2020)) +
          geom_point(aes(alpha=eff_cov)) +
          theme_bw() +
          ggtitle(paste0(country_code,iu0," (ESS:",round(ess,digits=0),")")) +
          scale_y_continuous(limits = c(-1,0)) +
          scale_x_continuous(limits = c(0.01,2.5),trans = "log10") +
          theme(legend.position="none")
        beta_baseline_third[[list_ind]] = p
        
        p = ggplot(data=beta_reductions,
                   aes(x = eff_cov, y = reduction2020)) +
          geom_point() +
          theme_bw() +
          ggtitle(paste0(country_code,iu0," (ESS:",round(ess,digits=0),")")) +
          scale_y_continuous(limits = c(-1,0)) +
          scale_x_continuous(limits = c(0,1)) +
          theme(legend.position="none")
        betadrop_coverage[[list_ind]] = p
        
        
        # plot infections
        p = ggplot(infections_melt[which(infections_melt$id %in% sampled_params_iu$seed),],
                   aes(x = year, y = prevalence, group=id)) +
          geom_line(aes(alpha = eff_cov), lwd=0.3) +
          theme_bw() +
          ggtitle(paste0(country_code,iu0," (ESS:",round(ess,digits=0),")")) +
          scale_y_continuous(limits = c(0,1)) +
          theme(legend.position="none")

        for (mda_ind in 1:nrow(mda_history)){
          p = p + geom_vline(xintercept=as.numeric(mda_history$Year[mda_ind],alpha=0.5),col="green")
        }
        infections_plots[[list_ind]] = p

        # plot trajectories
        p = ggplot(data=trajectories_melt[which(trajectories_melt$id %in% sampled_params_iu$seed),]) +
          geom_line(aes(x = year, y = prevalence, group=id, alpha = eff_cov), lwd=0.3) +
          theme_bw() +
          ggtitle(paste0(country_code,iu0," (ESS:",round(ess,digits=0),")")) +
          #scale_y_continuous(limits = c(0,1)) +
          theme(legend.position="none")
        for (mda_ind in 1:nrow(mda_history)){
          p = p + geom_vline(xintercept=as.numeric(mda_history$Year[mda_ind],alpha=0.5),col="green")
        }
        beta_multiplier = 1 #max(beta_summary$quant_0.975) # tried doing some scaling of 2nd beta axis to make the plots more legible but it didn't really work...
        p = p +
          geom_ribbon(data=beta_summary, aes(x = Year, ymin = (quant_0.025/beta_multiplier), ymax = (quant_0.975/beta_multiplier),alpha=0.5), col="blue",fill="blue") +
          geom_line(data=beta_summary, aes(x = Year, y = (quant_0.5)/beta_multiplier), col="blue") +
          scale_y_continuous("prevalence", sec.axis = sec_axis(~ .*beta_multiplier , name = "beta (blue ribbon)"), limits = c(0,1.25))
        for (year_ind in 1:length(prevalence_map)){
          p = p + geom_segment(y=prevalence_map[[year_ind]]$data[j,2],yend=prevalence_map[[year_ind]]$data[j,1],x=(years_vector_id[year_ind]+0.75),col="red")
        }

        trajectories_plots[[list_ind]] = p

        # # code to plot beta trajectories separately to the trajectories
        # p_beta = ggplot(beta_summary) +
        #   geom_line(aes(x = Year, y = quant_0.5)) +
        #   geom_ribbon(aes(x = Year, ymin = quant_0.025, ymax = quant_0.975,alpha=0.5)) +
        #   theme_bw() +
        #   ylab("beta")+
        #   theme(legend.position="none") +
        #   scale_y_continuous(limits=c(0,1.25))
        # 
        # p_joined = patchwork::wrap_plots(list(p,p_beta),nrow = 2, ncol = 1)
        # trajectories_plots[[list_ind]] = p_joined
        
      }
      list_counter = list_counter + length(inds_to_include)
    }
  }

  iu_order = iu_order %>% 
    filter(!is.na(iu_code)) %>%
    mutate(ISO3 = substr(iu_code,1,3),
           IU_ID = as.numeric(substr(iu_code,4,nchar(iu_code))))
  
  
  country_codes = unique(iu_order$ISO3)
  for (country in country_codes){
    print(country)

    country_inds = which(iu_order$ISO3==country)
    trajectories_plots_countries = trajectories_plots[country_inds]
    infections_plots_countries = infections_plots[country_inds]

    pdf(paste0("post_AMIS_analysis/plots/trajectories_posterior_",country,".pdf"),width=12,height=(ceiling(length(country_inds)/panel_ncols)*2))
    print(patchwork::wrap_plots(trajectories_plots_countries,nrow = ceiling(length(country_inds)/panel_ncols), ncol = panel_ncols))
    dev.off()

    pdf(paste0("post_AMIS_analysis/plots/infections_posterior_",country,".pdf"),width=10,height=(ceiling(length(country_inds)/panel_ncols)*2.5))
    print(patchwork::wrap_plots(infections_plots_countries,nrow = ceiling(length(country_inds)/panel_ncols), ncol = panel_ncols))
    dev.off()
  }
  
  pdf(paste0("post_AMIS_analysis/plots/beta_baseline_first_allIUs.pdf"),width=12,height=(ceiling(L/panel_ncols)*2))
  print(patchwork::wrap_plots(beta_baseline_first,nrow = ceiling(L/panel_ncols), ncol = panel_ncols))
  dev.off()

  pdf(paste0("post_AMIS_analysis/plots/beta_baseline_second_allIUs.pdf"),width=12,height=(ceiling(L/panel_ncols)*2))
  print(patchwork::wrap_plots(beta_baseline_second,nrow = ceiling(L/panel_ncols), ncol = panel_ncols))
  dev.off()

  pdf(paste0("post_AMIS_analysis/plots/beta_baseline_third_allIUs.pdf"),width=12,height=(ceiling(L/panel_ncols)*2))
  print(patchwork::wrap_plots(beta_baseline_third,nrow = ceiling(L/panel_ncols), ncol = panel_ncols))
  dev.off()

  pdf(paste0("post_AMIS_analysis/plots/betadrop_coverage_allIUs.pdf"),width=12,height=(ceiling(L/panel_ncols)*2))
  print(patchwork::wrap_plots(betadrop_coverage,nrow = ceiling(L/panel_ncols), ncol = panel_ncols))
  dev.off()
  
  # # just plot those that failed follow up surveys
  # load("Maps/trachoma_espen_data.Rdata")
  # 
  # # impact survey failures
  # ius_that_failed_impact_surveys = trachoma_espen_data %>%
  #   filter(type_prevalence %in% c("impact") & !tf_prev_upr == 0.05)
  # iu_impact_failures_inds = which(iu_order$IU_ID %in% ius_that_failed_impact_surveys$IU_ID)
  # 
  # pdf(paste0("post_AMIS_analysis/plots/failed_surveys_impact_trajectories.pdf"),width=12,height=(ceiling(length(iu_impact_failures_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(trajectories_plots[iu_impact_failures_inds],nrow = ceiling(length(iu_impact_failures_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/failed_surveys_impact_beta_baseline_first.pdf"),width=12,height=(ceiling(length(iu_impact_failures_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(beta_baseline_first[iu_impact_failures_inds],nrow = ceiling(length(iu_impact_failures_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/failed_surveys_impact_beta_baseline_second.pdf"),width=12,height=(ceiling(length(iu_impact_failures_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(beta_baseline_second[iu_impact_failures_inds],nrow = ceiling(length(iu_impact_failures_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/failed_surveys_impact_beta_baseline_third.pdf"),width=12,height=(ceiling(length(iu_impact_failures_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(beta_baseline_third[iu_impact_failures_inds],nrow = ceiling(length(iu_impact_failures_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/failed_surveys_impact_betadrop_coverage.pdf"),width=12,height=(ceiling(length(iu_impact_failures_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(betadrop_coverage[iu_impact_failures_inds],nrow = ceiling(length(iu_impact_failures_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # # surveillance survey failures
  # ius_that_failed_surveillance_surveys = trachoma_espen_data %>%
  #   filter(type_prevalence %in% c("surveillance") & !tf_prev_upr == 0.05)
  # iu_surveillance_failures_inds = which(iu_order$IU_ID %in% ius_that_failed_surveillance_surveys$IU_ID)
  # 
  # pdf(paste0("post_AMIS_analysis/plots/failed_surveys_surveillance_trajectories.pdf"),width=12,height=(ceiling(length(iu_surveillance_failures_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(trajectories_plots[iu_surveillance_failures_inds],nrow = ceiling(length(iu_surveillance_failures_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/failed_surveys_surveillance_beta_baseline_first.pdf"),width=12,height=(ceiling(length(iu_surveillance_failures_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(beta_baseline_first[iu_surveillance_failures_inds],nrow = ceiling(length(iu_surveillance_failures_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/failed_surveys_surveillance_beta_baseline_second.pdf"),width=12,height=(ceiling(length(iu_surveillance_failures_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(beta_baseline_second[iu_surveillance_failures_inds],nrow = ceiling(length(iu_surveillance_failures_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/failed_surveys_surveillance_beta_baseline_third.pdf"),width=12,height=(ceiling(length(iu_surveillance_failures_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(beta_baseline_third[iu_surveillance_failures_inds],nrow = ceiling(length(iu_surveillance_failures_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/failed_surveys_surveillance_betadrop_coverage.pdf"),width=12,height=(ceiling(length(iu_surveillance_failures_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(betadrop_coverage[iu_surveillance_failures_inds],nrow = ceiling(length(iu_surveillance_failures_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # 
  # # rebounds
  # ius_that_rebounds = trachoma_espen_data %>%
  #   filter(IU_ID %in% ius_that_failed_surveillance_surveys$IU_ID & (!IU_ID %in% ius_that_failed_impact_surveys$IU_ID))
  # iu_rebounds_inds = which(iu_order$IU_ID %in% ius_that_rebounds$IU_ID)
  # 
  # pdf(paste0("post_AMIS_analysis/plots/rebounds_trajectories.pdf"),width=12,height=(ceiling(length(iu_rebounds_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(trajectories_plots[iu_rebounds_inds],nrow = ceiling(length(iu_rebounds_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/rebounds_beta_baseline_first.pdf"),width=12,height=(ceiling(length(iu_rebounds_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(beta_baseline_first[iu_rebounds_inds],nrow = ceiling(length(iu_rebounds_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/rebounds_beta_baseline_second.pdf"),width=12,height=(ceiling(length(iu_rebounds_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(beta_baseline_second[iu_rebounds_inds],nrow = ceiling(length(iu_rebounds_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/rebounds_beta_baseline_third.pdf"),width=12,height=(ceiling(length(iu_rebounds_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(beta_baseline_third[iu_rebounds_inds],nrow = ceiling(length(iu_rebounds_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  # 
  # pdf(paste0("post_AMIS_analysis/plots/rebounds_betadrop_coverage.pdf"),width=12,height=(ceiling(length(iu_rebounds_inds)/panel_ncols)*2))
  # print(patchwork::wrap_plots(betadrop_coverage[iu_rebounds_inds],nrow = ceiling(length(iu_rebounds_inds)/panel_ncols), ncol = panel_ncols))
  # dev.off()
  
}

# Map median posterior prevalence, cov and beta
shape <- read_sf(dsn = "../ESPEN_IU_2021/", layer = "ESPEN_IU_2021")

# # IUs to analyse
# full_list_IUs = as.data.frame(ess_gt200) %>%
#   mutate_if(is.character, as.numeric)
# write.csv(full_list_IUs, file="trachoma_ius_fitted.csv",row.names=F)
# 
# # load list of IUs to include
# list_of_endemic_ius = read.csv("Maps/trachoma_IU_selection_with_IUcode.csv")
# 
# ess_filtered_ius = as.data.frame(ess_gt200) %>%
#   mutate_if(is.character, as.numeric) %>%
#   filter(IU_ID %in% unique(list_of_endemic_ius$IU_ID)) %>%# Claudio N-A's  list of ius
#   left_join(shape,by="IU_ID") %>%
#   mutate(IU_ID_full = paste0(ADMIN0ISO3,sprintf("%05d", as.integer(IU_ID))))
# 
# write.csv(ess_filtered_ius %>% select(IU_ID_full),file="trachoma_ius_to_run.csv",row.names=F)


# ####### previously treated but not treated in future projections
# load("Maps/trachoma_data.Rdata")
# # fit only to years with new data (based on when max_survey_n_clean changes)
# new_surveys_completed=trachoma_data %>% 
#   arrange(IU_ID,Year) %>%
#   group_by(IU_ID) %>%
#   mutate(max_survey_n_clean = ifelse(is.na(max_survey_n),-1,max_survey_n)) %>%
#   mutate(new_data = (!(max_survey_n_clean == lag(max_survey_n_clean))) | (Year==1996 & max_survey_n==0)) %>%
#   filter(new_data == TRUE)
# 
# 
# # IUs under post MDA surveillance
# post_mda_surv_label = "Under post-MDA surveillance" # set to <5% tf_prevalence
# trachoma_Geoconnect_to_IUID = read.csv("Maps/trachoma_IU_Match_PB.csv") # PB's lookup table fo mapping Geoconnect_ID to IU_ID
# trachoma_data_raw = read.csv("Maps/trachomaComb_IU.csv") %>%
#   left_join(trachoma_Geoconnect_to_IUID %>% select(Geoconnect_ID,IU_ID), by="Geoconnect_ID") %>%
#   filter(!is.na(IU_ID))
# 
# post_mda_surv = trachoma_data_raw %>%
#   filter(Year==2021 & tf_prevalence == post_mda_surv_label)
# 
# # get IUs to project
# treatment_years = trachoma_data %>%
#   filter(PC_in_group == 1) %>%
#   group_by(IU_ID) %>%
#   summarise(last_treated_year = floor(max(Year)),
#             first_treated_year = floor(min(Year))) %>%
#   filter(last_treated_year < 2020) %>%
#   filter(IU_ID %in% list_of_endemic_ius$IU_ID & !IU_ID %in% post_mda_surv$IU_ID) %>%
#   select(IU_ID)
# 
# write.csv(treatment_years,file="no_future_treatment.csv",row.names=F)

####### plot betas

load("sampled_params_all.Rdata")
load("Maps/trachoma_espen_data.Rdata")
load("Maps/trachoma_data.Rdata")

#From Deirdre: Before we make this decision, would it be possible to see the beta trends 
#in the fits by IU - either as a plot of beta over time, 
#or relative beta over time, or as a scatter plot beta1 against 
#bet2 etc so that we can see whether fitting a trend line would be plausible?

### summary plots
# prep data for plot
all_beta_vals = sampled_params_all %>%
  mutate(beta1996 = beta_init,
         beta2000 = beta_init*beta75,
         beta2010 = beta_init*beta75*beta85,
         beta2020 = beta_init*beta75*beta85*beta95) %>%
  select(IU_ID,paste0("beta",c(1996,2000,2010,2020))) 
all_beta_vals = data.frame(cbind(all_beta_vals[,"IU_ID"],
                           matrix(rep(all_beta_vals[,"beta1996"],4),ncol=4),
                           matrix(rep(all_beta_vals[,"beta2000"],10),ncol=10),
                           matrix(rep(all_beta_vals[,"beta2010"],10),ncol=10),
                           matrix(rep(all_beta_vals[,"beta2020"],5),ncol=5)))
colnames(all_beta_vals) = c("IU_ID",paste0(1996:2024))
#save(all_beta_vals,file="all_beta_vals.Rdata")

all_betas_medians = all_beta_vals %>%
  mutate_if(is.character, as.numeric) %>%
  group_by(IU_ID) %>% 
  summarise(across(starts_with(c("19","20")), function(x) quantile(x, prob = 0.5))) 
#save(all_betas_medians,file="all_betas_medians.Rdata")

all_betas_medians_melt = all_betas_medians %>%
  melt(id.vars = "IU_ID") %>%
  mutate(Year = as.numeric(substr(variable,1,4))) %>%
  select(IU_ID,Year,value) %>%
  ungroup()
colnames(all_betas_medians_melt) = c("IU_ID","Year","median")


all_betas_upper = all_beta_vals %>%
  mutate_if(is.character, as.numeric) %>%
  group_by(IU_ID) %>% 
  summarise(across(starts_with(c("19","20")), function(x) quantile(x, prob = 0.95)))
#save(all_betas_upper,file="all_betas_upper.Rdata")

all_betas_lower = all_beta_vals %>%
  mutate_if(is.character, as.numeric) %>%
  group_by(IU_ID) %>% 
  summarise(across(starts_with(c("19","20")), function(x) quantile(x, prob = 0.05))) 
#save(all_betas_lower,file="all_betas_lower.Rdata")

all_betas_lower_min = apply(all_betas_lower[,2:30],2,min)
all_betas_upper_max = apply(all_betas_upper[,2:30],2,max)
all_betas_bounds = data.frame(Year=1996:2024, lower_overall=all_betas_lower_min, upper_overall=all_betas_upper_max)

# trend lines
pdf(paste0("post_AMIS_analysis/plots/betatrend_posterior_summary.pdf"),width=5,height=5)
ggplot() +
  geom_ribbon(data=all_betas_bounds, aes(x = Year, ymin = lower_overall, ymax = upper_overall, alpha=0.5), colour="lightblue",fill="lightblue") +
  geom_line(data=all_betas_medians_melt, aes(x = Year, y = median, group=IU_ID, alpha=0.5)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1.5))  +
  ylab("beta")+
  theme(legend.position="none") 
dev.off()


#### Plot trend in beta for all IUs

L <- 2789 # total number of IUs
panel_ncols <- 4
panel_nrows <- ceiling(L/panel_ncols)
num_sub_samples_posterior <- 200

# draws from posterior
list_counter = 0
beta_plots = list()


for (id in c(id_vec)){

  load(paste0("AMIS_output/amis_output",id,".Rdata")) # loads amis_output
  prevalence_map = trachoma_maps
  prevalence_map = lapply(1:length(prevalence_map), function(t){
    output=list(data = as.matrix(prevalence_map[[t]]$data %>%
                                   filter(TaskID==id) %>%
                                   select(max_tf_prev_upr),ncol=1),
                likelihood = prevalence_map[[t]]$likelihood)
    rownames(output$data) = prevalence_map[[t]]$data$IU_ID[prevalence_map[[t]]$data$TaskID==id]
    return(output)
  })

  # just get maps with data for this batch
  map_bool_id = which(sapply(1:length(prevalence_map), function(t) !all(is.na(prevalence_map[[t]]$data))))
  years_vector_id = years_vector[map_bool_id]
  prevalence_map = lapply(map_bool_id, function(t) {
    data = as.data.frame(prevalence_map[[t]]$data) %>%
      mutate(tf_prev_lwr = case_when(
        max_tf_prev_upr == 0.05 ~ 0,
        max_tf_prev_upr == 0.099 ~ 0.05,
        max_tf_prev_upr == 0.299 ~ 0.1,
        max_tf_prev_upr == 0.499 ~ 0.3,
        max_tf_prev_upr == 1 ~ 0.5))
    return(list(data=data,likelihood=prevalence_map[[t]]$likelihood))
  })

  iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
  ess_all = amis_output$ess
  inds_to_exclude = which(ess_all < 200)
  inds_to_include = setdiff(1:length(rownames(prevalence_map[[1]]$data)),inds_to_exclude)

  if(length(inds_to_include)>0){ # make sure number of IUs with ESS >= 200 is nonzero

    # get list indices for this batch
    plot_indices = list_counter + (1:length(inds_to_include))

    for (j in inds_to_include){

      list_ind = plot_indices[which(inds_to_include==j)]
      iu = iu_names[j]
      ess = ess_all[j]

      country_code = country_iu_lookup$country[which(country_iu_lookup$IU_ID == iu)]
      iu0 <- sprintf("%05d", as.integer(iu))

      sampled_params_iu <- read.csv(paste0("post_AMIS_analysis/trachoma-source-data-20250225-espen-input-bet-csv/InputBet_",country_code,iu0,".csv"))

      # prep data for plot
      beta_vals = sampled_params_iu %>%
        select(paste0("beta",1996:2024))
      beta_summary = apply(beta_vals, 2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))
      colnames(beta_summary) = paste0(1996:2024)
      beta_summary %<>% melt() %>% pivot_wider(names_from="Var1",values_from="value")
      colnames(beta_summary) = c("Year","quant_0.025","quant_0.5","quant_0.975")

      # plot trajectories
      p = ggplot(beta_summary) +
        geom_line(aes(x = Year, y = quant_0.5)) +
        geom_ribbon(aes(x = Year, ymin = quant_0.025, ymax = quant_0.975,alpha=0.5)) +
        theme_bw() +
        ggtitle(paste0(iu_names[j])) +
        ylab("beta")+
        theme(legend.position="none")

      beta_plots[[list_ind]] = p
    }
    list_counter = list_counter + length(inds_to_include)
  }
}

pdf(paste0("post_AMIS_analysis/plots/betatrend_posterior_allIUs.pdf"),width=10,height=panel_nrows*2.5)
print(patchwork::wrap_plots(beta_plots,nrow = panel_nrows, ncol = panel_ncols))
dev.off()


####### make summary plots

# setup code to vary beta throughout simulation
randomWalk = TRUE
# parameters of beta distribution
randomWalk_a = 5
randomWalk_b = 1
# index of the column of params_full defined in transmission_model
yearschange_index = c(2000,2010,2020)-1926+1 # 1926 is data start date (1996) minus 70 year burn in.
n_timschange = length(yearschange_index)

rprior = function(n){

  params = matrix(NA,ncol=6,nrow=n)
  colnames(params) = c("beta_init",paste0("beta",yearschange_index),"eff_cov", "k_parameter")

  if (randomWalk==TRUE){

    params[,1] = rexp(n,rate=beta_rate)
    params[,2:4] = rbeta(n*3, randomWalk_a, randomWalk_b)
    params[,5] = rbeta(n,a_eff_cov,b_eff_cov)
    params[,6] = rtruncnorm(n,a=k_lb, b=k_ub, mean=k_mean, sd=k_sd)
  } else {

    params = matrix(NA,ncol=3,nrow=n)
    params[,1] = matrix(rexp(n,rate=beta_rate),ncol=1)
    params[,2] = rbeta(n,a_eff_cov,b_eff_cov)
    params[,3] = rtruncnorm(n,a=k_lb, b=k_ub, mean=k_mean, sd=k_sd)
  }

  return(params)
}

# prior density
dprior = function(x, log){
  if(randomWalk==TRUE){


    if(!is.matrix(x)){
      if (log){
        d = sum(dexp(x[1],rate=beta_rate,log=T),
                dbeta(x[2:4], randomWalk_a, randomWalk_b,log=T),
                dbeta(x[5],a_eff_cov,b_eff_cov,log=T),
                log(dtruncnorm(x[6],a=k_lb, b=k_ub, mean=k_mean, sd=k_sd)))
      } else {
        d = prod(dexp(x[1],rate=beta_rate,log=F),
                 dbeta(x[2:4], randomWalk_a, randomWalk_b,log=F),
                 dbeta(x[5],a_eff_cov,b_eff_cov,log=F),
                 dtruncnorm(x[6],a=k_lb, b=k_ub, mean=k_mean, sd=k_sd))
      }

    } else {
      if (log){
        d = sum(dexp(x[,1],rate=beta_rate,log=T),
                dbeta(x[,2:4], randomWalk_a, randomWalk_b,log=T),
                dbeta(x[,5],a_eff_cov,b_eff_cov,log=T),
                log(dtruncnorm(x[,6],a=k_lb, b=k_ub, mean=k_mean, sd=k_sd)))
      } else {
        d = prod(dexp(x[,1],rate=beta_rate,log=F),
                 dbeta(x[,2:4], randomWalk_a, randomWalk_b,log=F),
                 dbeta(x[,5],a_eff_cov,b_eff_cov,log=F),
                 dtruncnorm(x[,6],a=k_lb, b=k_ub, mean=k_mean, sd=k_sd))
      }
    }
  } else {

    if(!is.matrix(x)){
      if (log){
        d = sum(dexp(x[1],rate=beta_rate,log=T),
                dbeta(x[2],a_eff_cov,b_eff_cov,log=T),
                log(dtruncnorm(x[3],a=k_lb, b=k_ub, mean=k_mean, sd=k_sd)))
      } else {
        d = prod(dexp(x[1],rate=beta_rate,log=F),
                 dbeta(x[2],a_eff_cov,b_eff_cov,log=F),
                 dtruncnorm(x[3],a=k_lb, b=k_ub, mean=k_mean, sd=k_sd))
      }

    } else {
      if (log){
        d = sum(dexp(x[,1],rate=beta_rate,log=T),
                dbeta(x[,2],a_eff_cov,b_eff_cov,log=T),
                log(dtruncnorm(x[,3],a=k_lb, b=k_ub, mean=k_mean, sd=k_sd)))
      } else {
        d = prod(dexp(x[,1],rate=beta_rate,log=F),
                 dbeta(x[,2],a_eff_cov,b_eff_cov,log=F),
                 dtruncnorm(x[,3],a=k_lb, b=k_ub, mean=k_mean, sd=k_sd))
      }
    }
  }
  return(d)
}

prior = list(rprior=rprior,dprior=dprior)

priorbetamed <- qexp(0.5,rate=beta_rate)
priorbetalwr <- qexp(0.025,rate=beta_rate)
priorbetaupr <- qexp(0.975,rate=beta_rate)

priorbeta_multmed <- qbeta(0.5,randomWalk_a, randomWalk_b)
priorbeta_multlwr <- qbeta(0.025,randomWalk_a, randomWalk_b)
priorbeta_multupr <- qbeta(0.975,randomWalk_a, randomWalk_b)

priorcovmed <- qbeta(0.5,a_eff_cov,b_eff_cov)
priorcovlwr <- qbeta(0.025,a_eff_cov,b_eff_cov)
priorcovupr <- qbeta(0.975,a_eff_cov,b_eff_cov)

priork_parametermed <- qtruncnorm(0.5,a=k_lb, b=k_ub, mean=k_mean, sd=k_sd)
priork_parameterlwr <- qtruncnorm(0.025,a=k_lb, b=k_ub, mean=k_mean, sd=k_sd)
priork_parameterupr <- qtruncnorm(0.975,a=k_lb, b=k_ub, mean=k_mean, sd=k_sd)


#find first 3 maps for each IU
trachoma_histories = trachoma_data %>%
  filter(PC_in_group==1) %>%
  group_by(IU_ID) %>%
  summarise(first_mda_year = min(Year))

first_maps_all = c()
for (iu in unique(sampled_params_all$IU_ID)){
  #load map
  prevalence_map = trachoma_maps
  prevalence_map_iu = sapply(1:length(prevalence_map), function(t){
    output=as.matrix(prevalence_map[[t]]$data %>%
                       filter(IU_ID==iu) %>%
                       select(max_tf_prev_upr),ncol=1)
    return(output)
  })
  first_maps_ind = which(!is.na(prevalence_map_iu))
  first_maps_ind = first_maps_ind[1:min(length(first_maps_ind),3)]

  #first MDA year
  first_mda_year = as.numeric(trachoma_histories %>%
                                filter(IU_ID == iu) %>%
                                select(first_mda_year))

  #join
  first_maps_postMDA = as.numeric(c(years_vector[first_maps_ind], rep(NA,3-length(years_vector[first_maps_ind]))) >= first_mda_year)
  first_maps_prev = c(prevalence_map_iu[first_maps_ind], rep(NA,3-length(years_vector[first_maps_ind])))
  first_maps_all = rbind(first_maps_all,cbind(IU_ID=iu,index=1:3,first_maps_postMDA,first_maps_prev))
}

first_maps_all = as.data.frame(first_maps_all)
first_maps_all$first_maps_postMDA[is.na(first_maps_all$first_maps_postMDA) & !is.na(first_maps_all$first_maps_prev)] = 0

first_maps_prev_wide = first_maps_all %>%
  select(IU_ID,index,first_maps_prev) %>%
  pivot_wider(values_from=first_maps_prev, names_prefix="map_prev_", names_from=index) %>%
  mutate_if(is.character, as.numeric)
first_maps_prev_wide$IU_ID = as.numeric(first_maps_prev_wide$IU_ID)

first_maps_postMDA_wide = first_maps_all %>%
  select(IU_ID,index,first_maps_postMDA) %>%
  pivot_wider(values_from=first_maps_postMDA, names_prefix="postMDA_flag_", names_from=index) %>%
  mutate_if(is.character, as.numeric)
first_maps_postMDA_wide$IU_ID = as.numeric(first_maps_postMDA_wide$IU_ID)

#join all
dat <- sampled_params_all %>%
  mutate_if(is.character, as.numeric) %>%
  #left_join(ess_filtered_ius, by="IU_ID") %>%
  filter(ess >= 200) %>%
  left_join(first_maps_prev_wide,by="IU_ID") %>%
  left_join(first_maps_postMDA_wide, by="IU_ID")


mncov <- aggregate(eff_cov~IU_ID, data=dat, FUN="median")
lwrcov <- aggregate(eff_cov~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprcov <- aggregate(eff_cov~IU_ID, data=dat, FUN=quantile, probs=0.95)

cov <- (cbind(mncov, lwrcov[,2], uprcov[,2]))
colnames(cov) <- c("IU_ID", "cov", "covlwr", "covupr")

mnk_parameter <- aggregate(k_parameter~IU_ID, data=dat, FUN="median")
lwrk_parameter <- aggregate(k_parameter~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprk_parameter <- aggregate(k_parameter~IU_ID, data=dat, FUN=quantile, probs=0.95)

k_parameter <- (cbind(mnk_parameter, lwrk_parameter[,2], uprk_parameter[,2]))
colnames(k_parameter) <- c("IU_ID", "k_parameter", "k_parameterlwr", "k_parameterupr")

mnbeta_init <- aggregate(beta_init~IU_ID, data=dat, FUN="median")
lwrbeta_init <- aggregate(beta_init~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprbeta_init <- aggregate(beta_init~IU_ID, data=dat, FUN=quantile, probs=0.95)

beta_init <- (cbind(mnbeta_init, lwrbeta_init[,2], uprbeta_init[,2]))
colnames(beta_init) <- c("IU_ID", "beta_init", "beta_initlwr", "beta_initupr")

mnbeta75 <- aggregate(beta75~IU_ID, data=dat, FUN="median")
lwrbeta75 <- aggregate(beta75~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprbeta75 <- aggregate(beta75~IU_ID, data=dat, FUN=quantile, probs=0.95)

beta75 <- (cbind(mnbeta75, lwrbeta75[,2], uprbeta75[,2]))
colnames(beta75) <- c("IU_ID", "beta75", "beta75lwr", "beta75upr")

mnbeta85 <- aggregate(beta85~IU_ID, data=dat, FUN="median")
lwrbeta85 <- aggregate(beta85~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprbeta85 <- aggregate(beta85~IU_ID, data=dat, FUN=quantile, probs=0.95)

beta85 <- (cbind(mnbeta85, lwrbeta85[,2], uprbeta85[,2]))
colnames(beta85) <- c("IU_ID", "beta85", "beta85lwr", "beta85upr")

mnbeta95 <- aggregate(beta95~IU_ID, data=dat, FUN="median")
lwrbeta95 <- aggregate(beta95~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprbeta95 <- aggregate(beta95~IU_ID, data=dat, FUN=quantile, probs=0.95)

beta95 <- (cbind(mnbeta95, lwrbeta95[,2], uprbeta95[,2]))
colnames(beta95) <- c("IU_ID", "beta95", "beta95lwr", "beta95upr")


mnmap_prev_1 <- aggregate(map_prev_1~IU_ID, data=dat, FUN="median") %>%
  mutate(lwrmap_prev_1 = case_when(
    map_prev_1 == 0.05 ~ 0,
    map_prev_1 == 0.099 ~ 0.05,
    map_prev_1 == 0.299 ~ 0.1,
    map_prev_1 == 0.499 ~ 0.3,
    map_prev_1 == 1 ~ 0.5,
    is.na(map_prev_1) ~ NA))
map_prev_1 = as.data.frame(beta_init[,1,drop=F]) %>%
  left_join(mnmap_prev_1)
colnames(map_prev_1) <- c("IU_ID", "uprmap_prev_1","lwrmap_prev_1")

mnmap_prev_2 <- aggregate(map_prev_2~IU_ID, data=dat, FUN="median") %>%
  mutate(lwrmap_prev_2 = case_when(
    map_prev_2 == 0.05 ~ 0,
    map_prev_2 == 0.099 ~ 0.05,
    map_prev_2 == 0.299 ~ 0.1,
    map_prev_2 == 0.499 ~ 0.3,
    map_prev_2 == 1 ~ 0.5,
    is.na(map_prev_2) ~ NA))
map_prev_2 = as.data.frame(beta_init[,1,drop=F]) %>%
  left_join(mnmap_prev_2)
colnames(map_prev_2) <- c("IU_ID", "uprmap_prev_2","lwrmap_prev_2")

mnmap_prev_3 <- aggregate(map_prev_3~IU_ID, data=dat, FUN="median") %>%
  mutate(lwrmap_prev_3 = case_when(
    map_prev_3 == 0.05 ~ 0,
    map_prev_3 == 0.099 ~ 0.05,
    map_prev_3 == 0.299 ~ 0.1,
    map_prev_3 == 0.499 ~ 0.3,
    map_prev_3 == 1 ~ 0.5,
    is.na(map_prev_3) ~ NA))
map_prev_3 = as.data.frame(beta_init[,1,drop=F]) %>%
  left_join(mnmap_prev_3)
colnames(map_prev_3) <- c("IU_ID", "uprmap_prev_3","lwrmap_prev_3")


mnpostMDA_flag_1 <- aggregate(postMDA_flag_1~IU_ID, data=dat, FUN="median")
postMDA_flag_1 = as.data.frame(beta_init[,1,drop=F]) %>%
  left_join(mnpostMDA_flag_1)
colnames(postMDA_flag_1) <- c("IU_ID", "postMDA_flag_1")
postMDA_flag_1$postMDA_flag_1 = as.factor(postMDA_flag_1$postMDA_flag_1)


mnpostMDA_flag_2 <- aggregate(postMDA_flag_2~IU_ID, data=dat, FUN="median")
postMDA_flag_2 = as.data.frame(beta_init[,1,drop=F]) %>%
  left_join(mnpostMDA_flag_2)
colnames(postMDA_flag_2) <- c("IU_ID", "postMDA_flag_2")
postMDA_flag_2$postMDA_flag_2 = as.factor(postMDA_flag_2$postMDA_flag_2)


mnpostMDA_flag_3 <- aggregate(postMDA_flag_3~IU_ID, data=dat, FUN="median")
postMDA_flag_3 = as.data.frame(beta_init[,1,drop=F]) %>%
  left_join(mnpostMDA_flag_3)
colnames(postMDA_flag_3) <- c("IU_ID", "postMDA_flag_3")
postMDA_flag_3$postMDA_flag_3 = as.factor(postMDA_flag_3$postMDA_flag_3)

# # Create data frames to plot surveys by year
# # mnprev_1996 <- aggregate(prev_1996~IU_ID, data=dat, FUN="median")
# # lwrprev_1996 <- aggregate(prev_1996~IU_ID, data=dat, FUN=quantile, probs=0.05)
# # uprprev_1996 <- aggregate(prev_1996~IU_ID, data=dat, FUN=quantile, probs=0.95)
# #
# # prev_1996_subset <- cbind(mnprev_1996, lwrprev_1996[,2], uprprev_1996[,2])
# # prev_1996 = as.data.frame(beta[,1,drop=F]) %>%
# #   left_join(prev_1996_subset)
# # colnames(prev_1996) <- c("IU_ID", "prev_1996", "prev_1996lwr", "prev_1996upr")
#
# # mnprev_1997 <- aggregate(prev_1997~IU_ID, data=dat, FUN="median")
# # lwrprev_1997 <- aggregate(prev_1997~IU_ID, data=dat, FUN=quantile, probs=0.05)
# # uprprev_1997 <- aggregate(prev_1997~IU_ID, data=dat, FUN=quantile, probs=0.95)
# #
# # prev_1997_subset <- cbind(mnprev_1997, lwrprev_1997[,2], uprprev_1997[,2])
# # prev_1997 = as.data.frame(beta[,1,drop=F]) %>%
# #   left_join(prev_1997_subset)
# # colnames(prev_1997) <- c("IU_ID", "prev_1997", "prev_1997lwr", "prev_1997upr")
#
# mnprev_1999 <- aggregate(prev_1999~IU_ID, data=dat, FUN="median")
# lwrprev_1999 <- aggregate(prev_1999~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_1999 <- aggregate(prev_1999~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_1999_subset <- cbind(mnprev_1999, lwrprev_1999[,2], uprprev_1999[,2])
# prev_1999 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_1999_subset)
# colnames(prev_1999) <- c("IU_ID", "prev_1999", "prev_1999lwr", "prev_1999upr")
#
# # mnprev_2000 <- aggregate(prev_2000~IU_ID, data=dat, FUN="median")
# # lwrprev_2000 <- aggregate(prev_2000~IU_ID, data=dat, FUN=quantile, probs=0.05)
# # uprprev_2000 <- aggregate(prev_2000~IU_ID, data=dat, FUN=quantile, probs=0.95)
# #
# # prev_2000_subset <- cbind(mnprev_2000, lwrprev_2000[,2], uprprev_2000[,2])
# # prev_2000 = as.data.frame(beta_init[,1,drop=F]) %>%
# #   left_join(prev_2000_subset)
# # colnames(prev_2000) <- c("IU_ID", "prev_2000", "prev_2000lwr", "prev_2000upr")
#
# mnprev_2001 <- aggregate(prev_2001~IU_ID, data=dat, FUN="median")
# lwrprev_2001 <- aggregate(prev_2001~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2001 <- aggregate(prev_2001~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2001_subset <- cbind(mnprev_2001, lwrprev_2001[,2], uprprev_2001[,2])
# prev_2001 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2001_subset)
# colnames(prev_2001) <- c("IU_ID", "prev_2001", "prev_2001lwr", "prev_2001upr")
#
# mnprev_2002 <- aggregate(prev_2002~IU_ID, data=dat, FUN="median")
# lwrprev_2002 <- aggregate(prev_2002~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2002 <- aggregate(prev_2002~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2002_subset <- cbind(mnprev_2002, lwrprev_2002[,2], uprprev_2002[,2])
# prev_2002 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2002_subset)
# colnames(prev_2002) <- c("IU_ID", "prev_2002", "prev_2002lwr", "prev_2002upr")
#
# mnprev_2003 <- aggregate(prev_2003~IU_ID, data=dat, FUN="median")
# lwrprev_2003 <- aggregate(prev_2003~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2003 <- aggregate(prev_2003~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2003_subset <- cbind(mnprev_2003, lwrprev_2003[,2], uprprev_2003[,2])
# prev_2003 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2003_subset)
# colnames(prev_2003) <- c("IU_ID", "prev_2003", "prev_2003lwr", "prev_2003upr")
#
# mnprev_2004 <- aggregate(prev_2004~IU_ID, data=dat, FUN="median")
# lwrprev_2004 <- aggregate(prev_2004~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2004 <- aggregate(prev_2004~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2004_subset <- cbind(mnprev_2004, lwrprev_2004[,2], uprprev_2004[,2])
# prev_2004 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2004_subset)
# colnames(prev_2004) <- c("IU_ID", "prev_2004", "prev_2004lwr", "prev_2004upr")
#
# mnprev_2005 <- aggregate(prev_2005~IU_ID, data=dat, FUN="median")
# lwrprev_2005 <- aggregate(prev_2005~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2005 <- aggregate(prev_2005~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2005_subset <- cbind(mnprev_2005, lwrprev_2005[,2], uprprev_2005[,2])
# prev_2005 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2005_subset)
# colnames(prev_2005) <- c("IU_ID", "prev_2005", "prev_2005lwr", "prev_2005upr")
#
# mnprev_2006 <- aggregate(prev_2006~IU_ID, data=dat, FUN="median")
# lwrprev_2006 <- aggregate(prev_2006~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2006 <- aggregate(prev_2006~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2006_subset <- cbind(mnprev_2006, lwrprev_2006[,2], uprprev_2006[,2])
# prev_2006 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2006_subset)
# colnames(prev_2006) <- c("IU_ID", "prev_2006", "prev_2006lwr", "prev_2006upr")
#
# mnprev_2007 <- aggregate(prev_2007~IU_ID, data=dat, FUN="median")
# lwrprev_2007 <- aggregate(prev_2007~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2007 <- aggregate(prev_2007~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2007_subset <- cbind(mnprev_2007, lwrprev_2007[,2], uprprev_2007[,2])
# prev_2007 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2007_subset)
# colnames(prev_2007) <- c("IU_ID", "prev_2007", "prev_2007lwr", "prev_2007upr")
#
# mnprev_2008 <- aggregate(prev_2008~IU_ID, data=dat, FUN="median")
# lwrprev_2008 <- aggregate(prev_2008~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2008 <- aggregate(prev_2008~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2008_subset <- cbind(mnprev_2008, lwrprev_2008[,2], uprprev_2008[,2])
# prev_2008 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2008_subset)
# colnames(prev_2008) <- c("IU_ID", "prev_2008", "prev_2008lwr", "prev_2008upr")
#
# mnprev_2009 <- aggregate(prev_2009~IU_ID, data=dat, FUN="median")
# lwrprev_2009 <- aggregate(prev_2009~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2009 <- aggregate(prev_2009~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2009_subset <- cbind(mnprev_2009, lwrprev_2009[,2], uprprev_2009[,2])
# prev_2009 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2009_subset)
# colnames(prev_2009) <- c("IU_ID", "prev_2009", "prev_2009lwr", "prev_2009upr")
#
# mnprev_2010 <- aggregate(prev_2010~IU_ID, data=dat, FUN="median")
# lwrprev_2010 <- aggregate(prev_2010~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2010 <- aggregate(prev_2010~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2010_subset <- cbind(mnprev_2010, lwrprev_2010[,2], uprprev_2010[,2])
# prev_2010 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2010_subset)
# colnames(prev_2010) <- c("IU_ID", "prev_2010", "prev_2010lwr", "prev_2010upr")
#
# mnprev_2011 <- aggregate(prev_2011~IU_ID, data=dat, FUN="median")
# lwrprev_2011 <- aggregate(prev_2011~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2011 <- aggregate(prev_2011~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2011_subset <- cbind(mnprev_2011, lwrprev_2011[,2], uprprev_2011[,2])
# prev_2011 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2011_subset)
# colnames(prev_2011) <- c("IU_ID", "prev_2011", "prev_2011lwr", "prev_2011upr")
#
# mnprev_2012 <- aggregate(prev_2012~IU_ID, data=dat, FUN="median")
# lwrprev_2012 <- aggregate(prev_2012~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2012 <- aggregate(prev_2012~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2012_subset <- cbind(mnprev_2012, lwrprev_2012[,2], uprprev_2012[,2])
# prev_2012 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2012_subset)
# colnames(prev_2012) <- c("IU_ID", "prev_2012", "prev_2012lwr", "prev_2012upr")
#
# mnprev_2013 <- aggregate(prev_2013~IU_ID, data=dat, FUN="median")
# lwrprev_2013 <- aggregate(prev_2013~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2013 <- aggregate(prev_2013~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2013_subset <- cbind(mnprev_2013, lwrprev_2013[,2], uprprev_2013[,2])
# prev_2013 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2013_subset)
# colnames(prev_2013) <- c("IU_ID", "prev_2013", "prev_2013lwr", "prev_2013upr")
#
# mnprev_2014 <- aggregate(prev_2014~IU_ID, data=dat, FUN="median")
# lwrprev_2014 <- aggregate(prev_2014~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2014 <- aggregate(prev_2014~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2014_subset <- cbind(mnprev_2014, lwrprev_2014[,2], uprprev_2014[,2])
# prev_2014 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2014_subset)
# colnames(prev_2014) <- c("IU_ID", "prev_2014", "prev_2014lwr", "prev_2014upr")
#
# mnprev_2015 <- aggregate(prev_2015~IU_ID, data=dat, FUN="median")
# lwrprev_2015 <- aggregate(prev_2015~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2015 <- aggregate(prev_2015~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2015_subset <- cbind(mnprev_2015, lwrprev_2015[,2], uprprev_2015[,2])
# prev_2015 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2015_subset)
# colnames(prev_2015) <- c("IU_ID", "prev_2015", "prev_2015lwr", "prev_2015upr")
#
# mnprev_2016 <- aggregate(prev_2016~IU_ID, data=dat, FUN="median")
# lwrprev_2016 <- aggregate(prev_2016~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2016 <- aggregate(prev_2016~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2016_subset <- cbind(mnprev_2016, lwrprev_2016[,2], uprprev_2016[,2])
# prev_2016 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2016_subset)
# colnames(prev_2016) <- c("IU_ID", "prev_2016", "prev_2016lwr", "prev_2016upr")
#
# mnprev_2017 <- aggregate(prev_2017~IU_ID, data=dat, FUN="median")
# lwrprev_2017 <- aggregate(prev_2017~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2017 <- aggregate(prev_2017~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2017_subset <- cbind(mnprev_2017, lwrprev_2017[,2], uprprev_2017[,2])
# prev_2017 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2017_subset)
# colnames(prev_2017) <- c("IU_ID", "prev_2017", "prev_2017lwr", "prev_2017upr")
#
# mnprev_2018 <- aggregate(prev_2018~IU_ID, data=dat, FUN="median")
# lwrprev_2018 <- aggregate(prev_2018~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2018 <- aggregate(prev_2018~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2018_subset <- cbind(mnprev_2018, lwrprev_2018[,2], uprprev_2018[,2])
# prev_2018 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2018_subset)
# colnames(prev_2018) <- c("IU_ID", "prev_2018", "prev_2018lwr", "prev_2018upr")
#
# mnprev_2019 <- aggregate(prev_2019~IU_ID, data=dat, FUN="median")
# lwrprev_2019 <- aggregate(prev_2019~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2019 <- aggregate(prev_2019~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2019_subset <- cbind(mnprev_2019, lwrprev_2019[,2], uprprev_2019[,2])
# prev_2019 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2019_subset)
# colnames(prev_2019) <- c("IU_ID", "prev_2019", "prev_2019lwr", "prev_2019upr")
#
# mnprev_2020 <- aggregate(prev_2020~IU_ID, data=dat, FUN="median")
# lwrprev_2020 <- aggregate(prev_2020~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2020 <- aggregate(prev_2020~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2020_subset <- cbind(mnprev_2020, lwrprev_2020[,2], uprprev_2020[,2])
# prev_2020 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2020_subset)
# colnames(prev_2020) <- c("IU_ID", "prev_2020", "prev_2020lwr", "prev_2020upr")
#
# mnprev_2021 <- aggregate(prev_2021~IU_ID, data=dat, FUN="median")
# lwrprev_2021 <- aggregate(prev_2021~IU_ID, data=dat, FUN=quantile, probs=0.05)
# uprprev_2021 <- aggregate(prev_2021~IU_ID, data=dat, FUN=quantile, probs=0.95)
#
# prev_2021_subset <- cbind(mnprev_2021, lwrprev_2021[,2], uprprev_2021[,2])
# prev_2021 = as.data.frame(beta_init[,1,drop=F]) %>%
#   left_join(prev_2021_subset)
# colnames(prev_2021) <- c("IU_ID", "prev_2021", "prev_2021lwr", "prev_2021upr")


df <- cbind(cov,k_parameter[,2:4],beta_init[,2:4], beta75[,2:4],beta85[,2:4],beta95[,2:4],
            map_prev_1[,2:3],map_prev_2[,2:3],map_prev_3[,2:3],postMDA_flag_1[,2,drop=F],postMDA_flag_2[,2,drop=F],postMDA_flag_3[,2,drop=F])#,
            #prev_1999[,2:4],prev_2001[,2:4],prev_2002[,2:4],
            #prev_2003[,2:4],prev_2004[,2:4],prev_2005[,2:4],prev_2006[,2:4],prev_2007[,2:4],prev_2008[,2:4],
            #prev_2009[,2:4],prev_2010[,2:4],prev_2011[,2:4],prev_2012[,2:4],prev_2013[,2:4],prev_2014[,2:4],
            #prev_2015[,2:4],prev_2016[,2:4],prev_2017[,2:4],prev_2018[,2:4],prev_2019[,2:4],prev_2020[,2:4],prev_2021[,2:4]) #%>%
df <- df[order(df$beta_init),]
df$IUN <- seq(1,nrow(df))


pcov <- ggplot(data=df, aes(y=IUN, x=cov))+
  geom_segment(aes(x=covlwr, xend=covupr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("effective coverage") +
  theme_bw() +
  scale_y_continuous(name="IU")


pk_parameter <- ggplot(data=df, aes(y=IUN, x=k_parameter))+
  geom_segment(aes(x=k_parameterlwr, xend=k_parameterupr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("k parameter") +
  theme_bw() +
  scale_y_continuous(name="IU")

#labels = scales::trans_format("log10", scales::math_format(10^.x)))

pbeta_init <- ggplot(data=df, aes(y=IUN, x=beta_init))+
  geom_segment(aes(x=beta_initlwr, xend=beta_initupr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("baseline beta") +
  theme_bw() +
  scale_y_continuous(name="IU") +
  scale_x_continuous(trans="log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::number_format())


pbeta75 <- ggplot(data=df, aes(y=IUN, x=beta75))+
  geom_segment(aes(x=beta75lwr, xend=beta75upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("beta multiplier 2000") +
  theme_bw() +
  scale_y_continuous(name="IU") +
  scale_x_continuous(trans="log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::number_format(),
                     limits = c(0.15,1))


pbeta85 <- ggplot(data=df, aes(y=IUN, x=beta85))+
  geom_segment(aes(x=beta85lwr, xend=beta85upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("beta multiplier 2010") +
  theme_bw() +
  scale_y_continuous(name="IU") +
  scale_x_continuous(trans="log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::number_format(),
                     limits = c(0.15,1))


pbeta95 <- ggplot(data=df, aes(y=IUN, x=beta95))+
  geom_segment(aes(x=beta95lwr, xend=beta95upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("beta multiplier 2020") +
  theme_bw() +
  scale_y_continuous(name="IU") +
  scale_x_continuous(trans="log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::number_format(),
                     limits = c(0.15,1))




pmap_prev_1 <- ggplot(data=df, aes(y=IUN, x=map_prev_1))+
  geom_segment(aes(x=lwrmap_prev_1, xend=uprmap_prev_1, y=IUN, yend=IUN, col=postMDA_flag_1))+
  #geom_point(aes(col=postMDA_flag_1)) +
  xlab("1st map prevalence") +
  xlim(c(0,1)) +
  theme_bw() +
  theme(legend.position="none") +
  scale_y_continuous(name="IU")


pmap_prev_2 <- ggplot(data=df, aes(y=IUN, x=map_prev_2))+
  geom_segment(aes(x=lwrmap_prev_2, xend=uprmap_prev_2, y=IUN, yend=IUN, col=postMDA_flag_1))+
  #geom_point(aes(col=postMDA_flag_1)) +
  xlab("2nd map prevalence") +
  xlim(c(0,1)) +
  theme_bw() +
  theme(legend.position="none") +
  scale_y_continuous(name="IU")


pmap_prev_3 <- ggplot(data=df, aes(y=IUN, x=map_prev_3))+
  geom_segment(aes(x=lwrmap_prev_3, xend=uprmap_prev_3, y=IUN, yend=IUN, col=postMDA_flag_1))+
  #geom_point(aes(col=postMDA_flag_1)) +
  xlab("3rd map prevalence") +
  xlim(c(0,1)) +
  theme_bw() +
  theme(legend.position="none") +
  scale_y_continuous(name="IU")



png(paste0("post_AMIS_analysis/plots/posteriors_multipletimepts.png") , width=10, height=10,units ="in",res=1500)
grid.arrange(
  pcov+
    annotate(geom="segment", x=priorcovlwr, xend=priorcovupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorcovmed, y=mean(df$IUN), col="red"),
  pk_parameter+
    annotate(geom="segment", x=priork_parameterlwr, xend=priork_parameterupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priork_parametermed, y=mean(df$IUN), col="red"),
  pbeta_init +   annotate(geom="segment", x=priorbetalwr, xend=priorbetaupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorbetamed, y=mean(df$IUN), col="red"),
  pbeta75+   annotate(geom="segment", x=priorbeta_multlwr, xend=priorbeta_multupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorbeta_multmed, y=mean(df$IUN), col="red"),
  pbeta85+   annotate(geom="segment", x=priorbeta_multlwr, xend=priorbeta_multupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorbeta_multmed, y=mean(df$IUN), col="red"),
  pbeta95+   annotate(geom="segment", x=priorbeta_multlwr, xend=priorbeta_multupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorbeta_multmed, y=mean(df$IUN), col="red"),
  pmap_prev_1, pmap_prev_2, pmap_prev_3,
  nrow=3)
dev.off()


# Maps in Evandro's style
th <- theme(plot.title = element_text(size=12, hjust=0.5))
lwd_borders <- 0.1
colour <- "D"

#cov
cov_medians = mncov %>%
  mutate_if(is.character, as.numeric)
colnames(cov_medians) = c("IU_ID","cov")
espen_cov = shape %>%
  dplyr::select(IU_ID,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  left_join(cov_medians, by="IU_ID")
map_cov = st_as_sf(espen_cov)

#k_parameter
k_parameter_medians = mnk_parameter %>%
  mutate_if(is.character, as.numeric)
colnames(k_parameter_medians) = c("IU_ID","k_parameter")
espen_k_parameter = shape %>%
  dplyr::select(IU_ID,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  left_join(k_parameter_medians, by="IU_ID")
map_k_parameter = st_as_sf(espen_k_parameter)


#beta_init
beta_init_medians = mnbeta_init %>%
  mutate_if(is.character, as.numeric)
colnames(beta_init_medians) = c("IU_ID","beta_init")
espen_beta_init = shape %>%
  dplyr::select(IU_ID,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  left_join(beta_init_medians, by="IU_ID")
map_beta_init = st_as_sf(espen_beta_init)

#beta75
beta75_medians = mnbeta75 %>%
  mutate_if(is.character, as.numeric)
colnames(beta75_medians) = c("IU_ID","beta75")
espen_beta75 = shape %>%
  dplyr::select(IU_ID,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  left_join(beta75_medians, by="IU_ID")
map_beta75 = st_as_sf(espen_beta75)

#beta85
beta85_medians = mnbeta85 %>%
  mutate_if(is.character, as.numeric)
colnames(beta85_medians) = c("IU_ID","beta85")
espen_beta85 = shape %>%
  dplyr::select(IU_ID,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  left_join(beta85_medians, by="IU_ID")
map_beta85 = st_as_sf(espen_beta85)

#beta95
beta95_medians = mnbeta95 %>%
  mutate_if(is.character, as.numeric)
colnames(beta95_medians) = c("IU_ID","beta95")
espen_beta95 = shape %>%
  dplyr::select(IU_ID,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  left_join(beta95_medians, by="IU_ID")
map_beta95 = st_as_sf(espen_beta95)

#1st map
map_prev_1_medians = map_prev_1 %>%
  mutate_if(is.character, as.numeric)
colnames(map_prev_1_medians) = c("IU_ID","uprmap_prev_1","lwrmap_prev_1")
espen_map_prev_1 = shape %>%
  dplyr::select(IU_ID,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  left_join(map_prev_1_medians, by="IU_ID")
map_map_prev_1 = st_as_sf(espen_map_prev_1)

#2nd map
map_prev_2_medians = map_prev_2 %>%
  mutate_if(is.character, as.numeric)
colnames(map_prev_2_medians) = c("IU_ID","uprmap_prev_2","lwrmap_prev_2")
espen_map_prev_2 = shape %>%
  dplyr::select(IU_ID,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  left_join(map_prev_2_medians, by="IU_ID")
map_map_prev_2 = st_as_sf(espen_map_prev_2)

#3rd map
map_prev_3_medians = map_prev_3 %>%
  mutate_if(is.character, as.numeric)
colnames(map_prev_3_medians) = c("IU_ID","uprmap_prev_3","lwrmap_prev_3")
espen_map_prev_3 = shape %>%
  dplyr::select(IU_ID,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  left_join(map_prev_3_medians, by="IU_ID")
map_map_prev_3 = st_as_sf(espen_map_prev_3)



# Parameters estimates in Evandro's style
p_cov <- ggplot(map_cov) +
  geom_sf(aes(fill=cov), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median cov")) +
  ggtitle("effective coverage") + th
#p_cov

p_k_parameter <- ggplot(map_k_parameter) +
  geom_sf(aes(fill=k_parameter), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median k_parameter")) +
  ggtitle("k parameter") + th
#p_k_parameter

p_beta_init <- ggplot(map_beta_init) +
  geom_sf(aes(fill=beta_init), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median baseline beta")) +
  ggtitle("baseline beta") + th
#p_beta_init

p_beta75 <- ggplot(map_beta75) +
  geom_sf(aes(fill=beta75), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median beta multiplier 2000")) +
  ggtitle("beta multiplier 2000") + th

p_beta85 <- ggplot(map_beta85) +
  geom_sf(aes(fill=beta85), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median beta multiplier 2010")) +
  ggtitle("beta multiplier 2010") + th

p_beta95 <- ggplot(map_beta95) +
  geom_sf(aes(fill=beta95), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median beta multiplier 2020")) +
  ggtitle("beta multiplier 2020") + th

p_map_prev_1 <- ggplot(map_map_prev_1) +
  geom_sf(aes(fill=uprmap_prev_1), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="upper bound")) +
  ggtitle("1st map prevalence") + th

p_map_prev_2 <- ggplot(map_map_prev_2) +
  geom_sf(aes(fill=uprmap_prev_2), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="upper bound")) +
  ggtitle("2nd map prevalence") + th

p_map_prev_3 <- ggplot(map_map_prev_3) +
  geom_sf(aes(fill=uprmap_prev_3), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="upper bound")) +
  ggtitle("3rd map prevalence") + th


plot_pars <- p_cov + p_k_parameter + p_beta_init + p_beta75 + p_beta85 + p_beta95 + p_map_prev_1 + p_map_prev_2 + p_map_prev_3 & theme(legend.position = "right")
# plot_pars <- plot_pars + plot_layout(guides = "collect")
fileName <- paste0("post_AMIS_analysis/plots/map_parameters_all_countries.png")
ggsave(fileName, plot_pars, width = 400, height = 250, units = "mm")

# only ethiopia
p_cov_eth <- ggplot(map_cov%>% filter(ADMIN0ISO3=="ETH")) +
  geom_sf(aes(fill=cov), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median cov")) +
  ggtitle("effective coverage") + th
#p_cov

p_k_parameter_eth <- ggplot(map_k_parameter%>% filter(ADMIN0ISO3=="ETH")) +
  geom_sf(aes(fill=k_parameter), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median k_parameter")) +
  ggtitle("k parameter") + th
#p_k_parameter

p_beta_init_eth <- ggplot(map_beta_init%>% filter(ADMIN0ISO3=="ETH")) +
  geom_sf(aes(fill=beta_init), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median baseline beta")) +
  ggtitle("baseline beta") + th
#p_beta_init

p_beta75_eth <- ggplot(map_beta75%>% filter(ADMIN0ISO3=="ETH")) +
  geom_sf(aes(fill=beta75), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median beta multiplier 2000")) +
  ggtitle("beta multiplier 2000") + th

p_beta85_eth <- ggplot(map_beta85%>% filter(ADMIN0ISO3=="ETH")) +
  geom_sf(aes(fill=beta85), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median beta multiplier 2010")) +
  ggtitle("beta multiplier 2010") + th

p_beta95_eth <- ggplot(map_beta95%>% filter(ADMIN0ISO3=="ETH")) +
  geom_sf(aes(fill=beta95), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Median beta multiplier 2020")) +
  ggtitle("beta multiplier 2020") + th

p_map_prev_1_eth <- ggplot(map_map_prev_1%>% filter(ADMIN0ISO3=="ETH")) +
  geom_sf(aes(fill=uprmap_prev_1), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="upper bound")) +
  ggtitle("1st map prevalence") + th

p_map_prev_2_eth <- ggplot(map_map_prev_2%>% filter(ADMIN0ISO3=="ETH")) +
  geom_sf(aes(fill=uprmap_prev_2), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="upper bound")) +
  ggtitle("2nd map prevalence") + th

p_map_prev_3_eth <- ggplot(map_map_prev_3%>% filter(ADMIN0ISO3=="ETH")) +
  geom_sf(aes(fill=uprmap_prev_3), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="upper bound")) +
  ggtitle("3rd map prevalence") + th


plot_pars_eth <- p_cov_eth + p_k_parameter_eth + p_beta_init_eth + p_beta75_eth + p_beta85_eth + p_beta95_eth + p_map_prev_1_eth + p_map_prev_2_eth + p_map_prev_3_eth & theme(legend.position = "right")
fileName <- paste0("post_AMIS_analysis/plots/map_parameters_ETH.png")
ggsave(fileName, plot_pars_eth, width = 400, height = 250, units = "mm")


# # plot surveys by year
# pprev1996 <- ggplot(data=df, aes(y=IUN, x=prev_1996))+
#   geom_segment(aes(x=prev_1996lwr, xend=prev_1996upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 1996") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev1997 <- ggplot(data=df, aes(y=IUN, x=prev_1997))+
#   geom_segment(aes(x=prev_1997lwr, xend=prev_1997upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 1997") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev1999 <- ggplot(data=df, aes(y=IUN, x=prev_1999))+
#   geom_segment(aes(x=prev_1999lwr, xend=prev_1999upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 1999") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2000 <- ggplot(data=df, aes(y=IUN, x=prev_2000))+
#   geom_segment(aes(x=prev_2000lwr, xend=prev_2000upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2000") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2001 <- ggplot(data=df, aes(y=IUN, x=prev_2001))+
#   geom_segment(aes(x=prev_2001lwr, xend=prev_2001upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2001") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2002 <- ggplot(data=df, aes(y=IUN, x=prev_2002))+
#   geom_segment(aes(x=prev_2002lwr, xend=prev_2002upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2002") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2003 <- ggplot(data=df, aes(y=IUN, x=prev_2003))+
#   geom_segment(aes(x=prev_2003lwr, xend=prev_2003upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2003") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2004 <- ggplot(data=df, aes(y=IUN, x=prev_2004))+
#   geom_segment(aes(x=prev_2004lwr, xend=prev_2004upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2004") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2005 <- ggplot(data=df, aes(y=IUN, x=prev_2005))+
#   geom_segment(aes(x=prev_2005lwr, xend=prev_2005upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2005") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2006 <- ggplot(data=df, aes(y=IUN, x=prev_2006))+
#   geom_segment(aes(x=prev_2006lwr, xend=prev_2006upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2006") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2007 <- ggplot(data=df, aes(y=IUN, x=prev_2007))+
#   geom_segment(aes(x=prev_2007lwr, xend=prev_2007upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2007") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2008 <- ggplot(data=df, aes(y=IUN, x=prev_2008))+
#   geom_segment(aes(x=prev_2008lwr, xend=prev_2008upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2008") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2009 <- ggplot(data=df, aes(y=IUN, x=prev_2009))+
#   geom_segment(aes(x=prev_2009lwr, xend=prev_2009upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2009") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2010 <- ggplot(data=df, aes(y=IUN, x=prev_2010))+
#   geom_segment(aes(x=prev_2010lwr, xend=prev_2010upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2010") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2011 <- ggplot(data=df, aes(y=IUN, x=prev_2011))+
#   geom_segment(aes(x=prev_2011lwr, xend=prev_2011upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2011") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2012 <- ggplot(data=df, aes(y=IUN, x=prev_2012))+
#   geom_segment(aes(x=prev_2012lwr, xend=prev_2012upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2012") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2013 <- ggplot(data=df, aes(y=IUN, x=prev_2013))+
#   geom_segment(aes(x=prev_2013lwr, xend=prev_2013upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2013") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2014 <- ggplot(data=df, aes(y=IUN, x=prev_2014))+
#   geom_segment(aes(x=prev_2014lwr, xend=prev_2014upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2014") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2015 <- ggplot(data=df, aes(y=IUN, x=prev_2015))+
#   geom_segment(aes(x=prev_2015lwr, xend=prev_2015upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2015") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2016 <- ggplot(data=df, aes(y=IUN, x=prev_2016))+
#   geom_segment(aes(x=prev_2016lwr, xend=prev_2016upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2016") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2017 <- ggplot(data=df, aes(y=IUN, x=prev_2017))+
#   geom_segment(aes(x=prev_2017lwr, xend=prev_2017upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2017") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2018 <- ggplot(data=df, aes(y=IUN, x=prev_2018))+
#   geom_segment(aes(x=prev_2018lwr, xend=prev_2018upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2018") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2019 <- ggplot(data=df, aes(y=IUN, x=prev_2019))+
#   geom_segment(aes(x=prev_2019lwr, xend=prev_2019upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2019") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2020 <- ggplot(data=df, aes(y=IUN, x=prev_2020))+
#   geom_segment(aes(x=prev_2020lwr, xend=prev_2020upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2020") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
#
# pprev2021 <- ggplot(data=df, aes(y=IUN, x=prev_2021))+
#   geom_segment(aes(x=prev_2021lwr, xend=prev_2021upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
#   geom_point() +
#   xlab("Prevalence 2021") +
#   theme(axis.text.y = element_blank()) +
#   theme_bw() +
#   scale_y_continuous(name="IU") +
#   xlim(0,1)
#
# png(paste0("post_AMIS_analysis/plots/posterior_prevs_multipletimepts.png") , width=10, height=10,units ="in",res=1500)
# grid.arrange(
#   pprev1999,pprev2001,pprev2002,pprev2003,pprev2004,pprev2005,pprev2006,pprev2007,pprev2008,pprev2009,
#   pprev2010,pprev2011,pprev2012,pprev2013,pprev2014,pprev2015,pprev2016,pprev2017,pprev2018,pprev2019,pprev2020,pprev2021,
#   nrow=5)
# dev.off()
#
# #Prevalence
# # prev_1996_medians = mnprev_1996 %>%
# #   mutate_if(is.character, as.numeric)
# # colnames(prev_1996_medians) = c("IU_ID","prev_1996")
# # espen_prev_1996 = shape %>%
# #   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
# #   left_join(prev_1996_medians, by="IU_ID")
# # map_prev_1996 = st_as_sf(espen_prev_1996)
# #
# #
# # prev_1997_medians = mnprev_1997 %>%
# #   mutate_if(is.character, as.numeric)
# # colnames(prev_1997_medians) = c("IU_ID","prev_1997")
# # espen_prev_1997 = shape %>%
# #   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
# #   left_join(prev_1997_medians, by="IU_ID")
# # map_prev_1997 = st_as_sf(espen_prev_1997)
#
#
# prev_1999_medians = mnprev_1999 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_1999_medians) = c("IU_ID","prev_1999")
# espen_prev_1999 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_1999_medians, by="IU_ID")
# map_prev_1999 = st_as_sf(espen_prev_1999)
#
#
# # prev_2000_medians = mnprev_2000 %>%
# #   mutate_if(is.character, as.numeric)
# # colnames(prev_2000_medians) = c("IU_ID","prev_2000")
# # espen_prev_2000 = shape %>%
# #   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
# #   left_join(prev_2000_medians, by="IU_ID")
# # map_prev_2000 = st_as_sf(espen_prev_2000)
#
#
# prev_2001_medians = mnprev_2001 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2001_medians) = c("IU_ID","prev_2001")
# espen_prev_2001 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2001_medians, by="IU_ID")
# map_prev_2001 = st_as_sf(espen_prev_2001)
#
#
# prev_2002_medians = mnprev_2002 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2002_medians) = c("IU_ID","prev_2002")
# espen_prev_2002 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2002_medians, by="IU_ID")
# map_prev_2002 = st_as_sf(espen_prev_2002)
#
#
# prev_2003_medians = mnprev_2003 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2003_medians) = c("IU_ID","prev_2003")
# espen_prev_2003 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2003_medians, by="IU_ID")
# map_prev_2003 = st_as_sf(espen_prev_2003)
#
#
# prev_2004_medians = mnprev_2004 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2004_medians) = c("IU_ID","prev_2004")
# espen_prev_2004 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2004_medians, by="IU_ID")
# map_prev_2004 = st_as_sf(espen_prev_2004)
#
#
# prev_2005_medians = mnprev_2005 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2005_medians) = c("IU_ID","prev_2005")
# espen_prev_2005 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2005_medians, by="IU_ID")
# map_prev_2005 = st_as_sf(espen_prev_2005)
#
#
# prev_2006_medians = mnprev_2006 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2006_medians) = c("IU_ID","prev_2006")
# espen_prev_2006 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2006_medians, by="IU_ID")
# map_prev_2006 = st_as_sf(espen_prev_2006)
#
#
# prev_2007_medians = mnprev_2007 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2007_medians) = c("IU_ID","prev_2007")
# espen_prev_2007 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2007_medians, by="IU_ID")
# map_prev_2007 = st_as_sf(espen_prev_2007)
#
#
# prev_2008_medians = mnprev_2008 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2008_medians) = c("IU_ID","prev_2008")
# espen_prev_2008 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2008_medians, by="IU_ID")
# map_prev_2008 = st_as_sf(espen_prev_2008)
#
#
# prev_2009_medians = mnprev_2009 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2009_medians) = c("IU_ID","prev_2009")
# espen_prev_2009 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2009_medians, by="IU_ID")
# map_prev_2009 = st_as_sf(espen_prev_2009)
#
#
# prev_2010_medians = mnprev_2010 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2010_medians) = c("IU_ID","prev_2010")
# espen_prev_2010 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2010_medians, by="IU_ID")
# map_prev_2010 = st_as_sf(espen_prev_2010)
#
#
# prev_2011_medians = mnprev_2011 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2011_medians) = c("IU_ID","prev_2011")
# espen_prev_2011 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2011_medians, by="IU_ID")
# map_prev_2011 = st_as_sf(espen_prev_2011)
#
#
# prev_2012_medians = mnprev_2012 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2012_medians) = c("IU_ID","prev_2012")
# espen_prev_2012 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2012_medians, by="IU_ID")
# map_prev_2012 = st_as_sf(espen_prev_2012)
#
#
# prev_2013_medians = mnprev_2013 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2013_medians) = c("IU_ID","prev_2013")
# espen_prev_2013 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2013_medians, by="IU_ID")
# map_prev_2013 = st_as_sf(espen_prev_2013)
#
#
# prev_2014_medians = mnprev_2014 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2014_medians) = c("IU_ID","prev_2014")
# espen_prev_2014 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2014_medians, by="IU_ID")
# map_prev_2014 = st_as_sf(espen_prev_2014)
#
#
# prev_2015_medians = mnprev_2015 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2015_medians) = c("IU_ID","prev_2015")
# espen_prev_2015 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2015_medians, by="IU_ID")
# map_prev_2015 = st_as_sf(espen_prev_2015)
#
#
# prev_2016_medians = mnprev_2016 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2016_medians) = c("IU_ID","prev_2016")
# espen_prev_2016 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2016_medians, by="IU_ID")
# map_prev_2016 = st_as_sf(espen_prev_2016)
#
#
# prev_2017_medians = mnprev_2017 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2017_medians) = c("IU_ID","prev_2017")
# espen_prev_2017 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2017_medians, by="IU_ID")
# map_prev_2017 = st_as_sf(espen_prev_2017)
#
#
# prev_2018_medians = mnprev_2018 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2018_medians) = c("IU_ID","prev_2018")
# espen_prev_2018 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2018_medians, by="IU_ID")
# map_prev_2018 = st_as_sf(espen_prev_2018)
#
#
# prev_2019_medians = mnprev_2019 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2019_medians) = c("IU_ID","prev_2019")
# espen_prev_2019 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2019_medians, by="IU_ID")
# map_prev_2019 = st_as_sf(espen_prev_2019)
#
#
# prev_2020_medians = mnprev_2020 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2020_medians) = c("IU_ID","prev_2020")
# espen_prev_2020 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2020_medians, by="IU_ID")
# map_prev_2020 = st_as_sf(espen_prev_2020)
#
#
# prev_2021_medians = mnprev_2021 %>%
#   mutate_if(is.character, as.numeric)
# colnames(prev_2021_medians) = c("IU_ID","prev_2021")
# espen_prev_2021 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(prev_2021_medians, by="IU_ID")
# map_prev_2021 = st_as_sf(espen_prev_2021)
#
#
# # p_prev1996 <- ggplot(map_prev_1996) +
# #   geom_sf(aes(fill=prev_1996), lwd=lwd_borders) +
# #   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
# #   theme_bw() + theme(legend.position="right") +
# #   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
# #   ggtitle(as.character(years_vector[which(years_vector==1996)])) + th
# #
# # p_prev1997 <- ggplot(map_prev_1997) +
# #   geom_sf(aes(fill=prev_1997), lwd=lwd_borders) +
# #   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
# #   theme_bw() + theme(legend.position="right") +
# #   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
# #   ggtitle(as.character(years_vector[which(years_vector==1997)])) + th
#
# p_prev1999 <- ggplot(map_prev_1999) +
#   geom_sf(aes(fill=prev_1999), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==1999)])) + th
#
# # p_prev2000 <- ggplot(map_prev_2000) +
# #   geom_sf(aes(fill=prev_2000), lwd=lwd_borders) +
# #   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
# #   theme_bw() + theme(legend.position="right") +
# #   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
# #   ggtitle(as.character(years_vector[which(years_vector==2000)])) + th
#
# p_prev2001 <- ggplot(map_prev_2001) +
#   geom_sf(aes(fill=prev_2001), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2001)])) + th
#
# p_prev2002 <- ggplot(map_prev_2002) +
#   geom_sf(aes(fill=prev_2002), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2002)])) + th
#
# p_prev2003 <- ggplot(map_prev_2003) +
#   geom_sf(aes(fill=prev_2003), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2003)])) + th
#
# p_prev2004 <- ggplot(map_prev_2004) +
#   geom_sf(aes(fill=prev_2004), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2004)])) + th
#
# p_prev2005 <- ggplot(map_prev_2005) +
#   geom_sf(aes(fill=prev_2005), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2005)])) + th
#
# p_prev2006 <- ggplot(map_prev_2006) +
#   geom_sf(aes(fill=prev_2006), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2006)])) + th
#
# p_prev2007 <- ggplot(map_prev_2007) +
#   geom_sf(aes(fill=prev_2007), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2007)])) + th
#
# p_prev2008 <- ggplot(map_prev_2008) +
#   geom_sf(aes(fill=prev_2008), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2008)])) + th
#
# p_prev2009 <- ggplot(map_prev_2009) +
#   geom_sf(aes(fill=prev_2009), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2009)])) + th
#
# p_prev2010 <- ggplot(map_prev_2010) +
#   geom_sf(aes(fill=prev_2010), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2010)])) + th
#
# p_prev2011 <- ggplot(map_prev_2011) +
#   geom_sf(aes(fill=prev_2011), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2011)])) + th
#
# p_prev2012 <- ggplot(map_prev_2012) +
#   geom_sf(aes(fill=prev_2012), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2012)])) + th
#
# p_prev2013 <- ggplot(map_prev_2013) +
#   geom_sf(aes(fill=prev_2013), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2013)])) + th
#
# p_prev2014 <- ggplot(map_prev_2014) +
#   geom_sf(aes(fill=prev_2014), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2014)])) + th
#
# p_prev2015 <- ggplot(map_prev_2015) +
#   geom_sf(aes(fill=prev_2015), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2015)])) + th
#
# p_prev2016 <- ggplot(map_prev_2016) +
#   geom_sf(aes(fill=prev_2016), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2016)])) + th
#
# p_prev2017 <- ggplot(map_prev_2017) +
#   geom_sf(aes(fill=prev_2017), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2017)])) + th
#
# p_prev2018 <- ggplot(map_prev_2018) +
#   geom_sf(aes(fill=prev_2018), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2018)])) + th
#
# p_prev2019 <- ggplot(map_prev_2019) +
#   geom_sf(aes(fill=prev_2019), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2019)])) + th
#
# p_prev2020 <- ggplot(map_prev_2020) +
#   geom_sf(aes(fill=prev_2020), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2020)])) + th
#
# p_prev2021 <- ggplot(map_prev_2021) +
#   geom_sf(aes(fill=prev_2021), lwd=lwd_borders) +
#   scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
#   theme_bw() + theme(legend.position="right") +
#   guides(fill=guide_colorbar(title="Median estimated prevalence")) +
#   ggtitle(as.character(years_vector[which(years_vector==2021)])) + th
#
#
# plot_prevs <- p_prev1999 + p_prev2001 + p_prev2002 +
#   p_prev2003 + p_prev2004 + p_prev2005 + p_prev2006 + p_prev2007 + p_prev2008 +
#   p_prev2009 + p_prev2010 + p_prev2011 + p_prev2012 + p_prev2013 + p_prev2014 + p_prev2015 +
#   p_prev2016 + p_prev2017 + p_prev2018 + p_prev2019 + p_prev2020 + p_prev2021 & theme(legend.position = "right")
# plot_prevs <- plot_prevs + plot_layout(guides = "collect")
#
# fileName <- paste0("post_AMIS_analysis/plots/map_median_prevs_all_countries.png")
# ggsave(fileName, plot_prevs, width = 550, height = 500, units = "mm")

