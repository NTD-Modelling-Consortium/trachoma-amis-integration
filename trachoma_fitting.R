
# on cluster:
#id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# for testing:
id = 556
library(dplyr)
library(reticulate)

# load trachoma model 
# first 3 lines only necessary on local machine - on cluster will install R package, set working dir. and source venv before running script
library(AMISforInfectiousDiseases)
setwd("~/Documents/trachoma-endgame/trachoma-amis-integration/")
reticulate::use_virtualenv("~/Documents/trachoma-endgame/trachoma-venv", required=TRUE)
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
beta_rate = 2
beta_lb = 0
beta_ub = Inf

eff_cov_lb = 0
eff_cov_ub = 1
a_eff_cov = 2
b_eff_cov = 2

# random number generator
rprior = function(n){
  params = matrix(NA,ncol=2,nrow=n)
  params[,1] = rexp(n,rate=beta_rate)
  params[,2] = rbeta(n,a_eff_cov,b_eff_cov)
  colnames(params) = c("beta","eff_cov")
  return(params)
}

# prior density
dprior = function(x, log){
  if(!is.matrix(x)){
    if (log){
      d = sum(dexp(x[1],rate=beta_rate,log=T),
              dbeta(x[2],a_eff_cov,b_eff_cov,log=T))
    } else {
      d = prod(dunif(x[1],rate=beta_rate,log=F),
               dbeta(x[2],a_eff_cov,b_eff_cov,log=F))
    }
    
  } else {
    if (log){
      d = sum(dunif(x[,1],rate=beta_rate,log=T),
              dbeta(x[,2],a_eff_cov,b_eff_cov,log=T))
    } else {
      d = prod(dunif(x[,1],rate=beta_rate,log=F),
               dbeta(x[,2],a_eff_cov,b_eff_cov,log=F))
    }
  }
  return(d)
}

prior = list(rprior=rprior,dprior=dprior)


# Algorithm parameters
amis_params<-default_amis_params()
amis_params$max_iters=10
amis_params$n_samples=500
amis_params$target_ess =500
amis_params$sigma=0.0025
amis_params$boundaries=c(-Inf,Inf)
amis_params$boundaries_param = matrix(c(beta_lb,eff_cov_lb,beta_ub,eff_cov_ub),ncol=2)

# shell to save trajectories
trajectories = c() # save simulated trajectories as code is running
if (!dir.exists("../trajectories")) {dir.create("../trajectories")}
save(trajectories,file=paste0("../trajectories/trajectories_",id,".Rdata"))

# Fit to years_vector years (70 years burn in, currently I think maps are at the start of the year/end of previous year since index starts at 0)
all_years_vector_id = seq(min(years_vector_id),(max(years_vector_id)+1),by=0.25)
weeks_indices <- as.integer(70*52 + (all_years_vector_id - min(years_vector))*52)
initial_infect = 0.5
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
      output=model_func(seeds, params, n_tims)
      
      #save simulated trajectories 
      load(paste0("../trajectories/trajectories_",id,".Rdata"))
      trajectories =  rbind(trajectories,output)
      save(trajectories, file=paste0("../trajectories/trajectories_",id,".Rdata"))
      
      output_map_years = output[,which(all_years_vector_id %in% (years_vector_id+0.75))]
      return(output_map_years)
      
  #   },
  #   error=error_function
  # )
}
#test=transmission_model(1:10,rprior(10),21)

st<-Sys.time()
amis_output <- AMISforInfectiousDiseases::amis(
    prevalence_map,
    seed=id,
    amis_params=amis_params,
    transmission_model=transmission_model,
    prior,
    output_dir = paste0("../amis_intermittant_outputs/batch",id,"_")
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

