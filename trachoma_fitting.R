
# on cluster:
id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# for testing:
#id = 11

library(dplyr)
library(reticulate)
library(AMISforInfectiousDiseases)
library(truncnorm)

# setup code to vary beta throughout simulation
randomWalk = TRUE
# parameters of beta distribution
randomWalk_a = 5
randomWalk_b = 1
# index of the column of params_full defined in transmission_model
yearschange_index = c(2000,2010,2020)-1926+1 # 1926 is data start date (1996) minus 70 year burn in. 
n_timschange = length(yearschange_index)
 
# load trachoma model 
# for testing
#setwd("~/Documents/trachoma-endgame/trachoma-amis-integration/")
#reticulate::use_virtualenv("~/Documents/trachoma-endgame/trachoma-venv", required=TRUE)
reticulate::use_virtualenv("../trachoma-venv", required=TRUE)
reticulate::py_config()
amis_int_mod <- import("trachoma_amis")

# Load prevalence map and filter rows for TaskID == id
load(paste0("../Maps/trachoma_maps.rds"))
load("../Maps/trachoma_map_years.rds")
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

# Algorithm parameters
amis_params<-default_amis_params()
amis_params$max_iters=50
amis_params$n_samples=1000
amis_params$target_ess=500
amis_params$sigma=0.0025
amis_params$boundaries=c(-Inf,Inf)

# shell to save trajectories
trajectories = c() # save simulated trajectories as code is running
if (!dir.exists("../trajectories")) {dir.create("../trajectories")}
save(trajectories,file=paste0("../trajectories/trajectories_",id,".Rdata"))

# save infection
infections = c() # save simulated infections as code is running
if (!dir.exists("../infections")) {dir.create("../infections")}
save(infections,file=paste0("../infections/infections",id,".Rdata"))

# Fit to years_vector years
# Index starts at 0
# (70 years burn in, currently I think maps are at the end of the year)
#all_years_vector_id = seq(min(years_vector_id),(max(years_vector_id)+1),by=0.25) 
all_years_vector_id = seq(1996,2022,by=0.25) 
weeks_indices <- as.integer(70*52-1 + (all_years_vector_id - min(years_vector))*52) 
initial_infect = 0.1
num_cores = 8L

model_func <- amis_int_mod$build_transmission_model(
    weeks_indices, initial_infect, num_cores
)

error_function <- function(e) {
    err <- reticulate::py_last_error()
    msg <- sprintf(
        "Exception type %s was raised during the execution of the model:",
        err$type
    )
    message(msg)
    stop(print(err))
}
# Wrap model_func so as to print the Python traceback
# in case of a raised exception during execution of the
# Python model. See
# https://rstudio.github.io/reticulate/reference/py_last_error.html
transmission_model = function(seeds, params, n_tims) {
  # tryCatch(
  #   expr={
       
      # fill in beta years so there is 1 beta column every year
      params_full = matrix(NA,ncol=98,nrow=nrow(params))
      colnames(params_full) = c(paste0("beta",1926:2021),"eff_cov", "k_parameter")
  
      if(randomWalk==T){
        
        n_reps_1 = yearschange_index[1] - 1
        params_full[,1:n_reps_1] = matrix(rep(params[,1],n_reps_1),ncol=n_reps_1)
        
        n_reps_2 = yearschange_index[2] - yearschange_index[1] 
        params_full[,yearschange_index[1]:(yearschange_index[2]-1)] = matrix(rep(params[,1]*params[,2],n_reps_2),ncol=n_reps_2)
        
        n_reps_3 = yearschange_index[3] - yearschange_index[2] 
        params_full[,yearschange_index[2]:(yearschange_index[3]-1)] = matrix(rep(params[,1]*params[,2]*params[,3],n_reps_3),ncol=n_reps_3)
        
        n_reps_4 = (2022-1926+1) - yearschange_index[3] 
        params_full[,yearschange_index[3]:(2022-1926)] = matrix(rep(params[,1]*params[,2]*params[,3]*params[,4],n_reps_4),ncol=n_reps_4)
        
        params_full[,97:98] = params[,5:6]
        
      } else {
        
        n_reps = (2022-1926+1) - 1 
        params_full[,1:(2022-1926)] = matrix(rep(params[,1],n_reps),ncol=n_reps)
        
        params_full[,97:98] = params[,2:3]
        
      }
  
      # run model
      output=model_func(seeds, params_full, n_tims)
      
      #save simulated trajectories 
      load(paste0("../trajectories/trajectories_",id,".Rdata"))
      trajectories =  rbind(trajectories,output[[1]])
      #trajectories =  rbind(trajectories,output)
      save(trajectories, file=paste0("../trajectories/trajectories_",id,".Rdata"))


      #save simulated infections
      load(paste0("../infections/infections",id,".Rdata"))
      infections =  rbind(infections,output[[2]])
      save(infections, file=paste0("../infections/infections",id,".Rdata"))
      
      output_map_years = output[[1]][,which(all_years_vector_id %in% (years_vector_id+0.75)),drop=F]
      #output_map_years = output[,which(all_years_vector_id %in% (years_vector_id+0.75))]
    
      return(output_map_years)
      
  #   },
  #   error=error_function
  # )
}

# # test model
# # shell to save trajectories
# trajectories = c() # save simulated trajectories as code is running
# if (!dir.exists("../test-trajectories")) {dir.create("../test-trajectories")}
# save(trajectories,file=paste0("../test-trajectories/trajectories_",id,".Rdata"))
# # save infection
# infections = c() # save simulated infections as code is running
# if (!dir.exists("../test-infections")) {dir.create("../test-infections")}
# save(infections,file=paste0("../test-infections/infections",id,".Rdata"))
# transmission_model = function(seeds, params, n_tims) {
#   # tryCatch(
#   #   expr={
#   output=model_func(seeds, params, n_tims)
# 
#   #save simulated trajectories
#   load(paste0("../test-trajectories/trajectories_",id,".Rdata"))
#   trajectories =  rbind(trajectories,output[[1]])
#   #trajectories =  rbind(trajectories,output)
#   save(trajectories, file=paste0("../test-trajectories/trajectories_",id,".Rdata"))
# 
# 
#   #save simulated infections
#   load(paste0("../test-infections/infections",id,".Rdata"))
#   infections =  rbind(infections,output[[2]])
#   save(infections, file=paste0("../test-infections/infections",id,".Rdata"))
# 
#   output_map_years = output[[1]][,which(all_years_vector_id %in% (years_vector_id+0.75))]
#   #output_map_years = output[,which(all_years_vector_id %in% (years_vector_id+0.75))]
# 
#   return(output_map_years)
# 
#   #   },
#   #   error=error_function
#   # )
# }
# params = cbind(seed=1:10,rprior(10))
# #sampled_params_test = read.csv("../sampled_params_50886.csv")
# #params = as.matrix(sampled_params_test %>%
# #  select(seed,beta,eff_cov))
# #params[,2] = params[,2]
# #params[,3] = 0.85
# test=transmission_model(as.integer(params[,1]),params[,2:3],45)
# par(mfrow=c(2,1))
# load(paste0("../test-infections/infections",id,".Rdata"))
# plot(x=all_years_vector_id,y=infections[1,],type="l",ylim=c(0,0.6),col="grey")
# for (i in 2:100){lines(x=all_years_vector_id,y=infections[i,],col="grey")}
# lines(x=all_years_vector_id,y=apply(infections,2,mean),lwd=2,col="blue")
# load(paste0("../test-trajectories/trajectories_",id,".Rdata"))
# plot(x=all_years_vector_id,y=trajectories[1,],type="l",ylim=c(0,0.7),col="grey")
# for (i in 2:100){lines(x=all_years_vector_id,y=trajectories[i,],col="grey")}
# lines(x=all_years_vector_id,y=apply(trajectories,2,mean),lwd=2,col="blue")

st<-Sys.time()
amis_output <- AMISforInfectiousDiseases::amis(
    prevalence_map,
    seed=id,
    amis_params=amis_params,
    transmission_model=transmission_model,
    prior
)
en<-Sys.time()
dur_amis<-as.numeric(difftime(en,st,units="mins"))
if (!dir.exists("../AMIS_output")) {dir.create("../AMIS_output")}
save(amis_output,file=paste0("../AMIS_output/amis_output",id,".Rdata"))


# Currently errors - I think because I
# don't know where the weights need to be set
# "No weight on any particles for locations in the active set."
print(amis_output)

# save summary
ess<-amis_output$ess
n_success<-length(which(ess>=amis_params[["target_ess"]]))
failures<-which(ess<amis_params[["target_ess"]])
n_failure<-length(failures)
if (n_failure>0) {cat(paste(failures,id,ess[failures]),file = paste0("../ESS_NOT_REACHED.txt"),sep = "\n", append = TRUE)}
if (!file.exists(paste0("../summary.csv"))) {cat("ID,n_failure,n_success,n_sim,min_ess,duration_amis,durarion_subsampling\n",file=paste0("../summary.csv"))}
cat(id,n_failure,n_success,length(amis_output$seeds),min(ess),dur_amis,NA,"\n",sep=",",file=paste0("../summary.csv"),append=TRUE)

