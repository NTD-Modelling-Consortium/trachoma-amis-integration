library(reticulate)
amis_int_mod <- import("trachoma_amis")

# Example prevalence map, with two locations, fitting to two time times
# Both locations start at 0.031, and the second time point is 0.021
prevalence_map <- vector("list", 2)
prevalence_map[[1]]$data <- matrix(c(0.031, 0.031))
prevalence_map[[2]]$data <- matrix(c(0.021, 0.021))

#' the "dprior" function
#' Note the second parameter _must_ be called log
#' Unclear what this function does
density_function <- function(parameters, log) {
    return(0.5)
}

#' The "rprior" function that returns a matrix whose columns are the parameters
#' and each row is a sample
rnd_function <- function(num_samples) {
    # set beta to just be 0.04 for each sample
    beta <- rep(0.04, num_samples)
  
    # set coverage to be randomly generated numbers between 0 and 1
    coverage <- runif(num_samples, min = 0, max = 1)
    
    # Combine the two columns to form the params matrix
    result <- cbind(beta, coverage)
    
    return(result)

}

prior <- list("dprior" = density_function, "rprior" = rnd_function)

#' Fit to initial and last week (in this case 100 years)
weeks_indices <- c(0L, 5199L)
initial_infect = 0.01
num_cores = 1L

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least two arguments must be supplied", call.=FALSE)
} else if (length(args)==2) {
  # No data location provided
  args[3] <- NULL
}

mda_coverage_data_filename <- args[1]
vaccine_coverage_data_filename <- args[2]
coverage_data_path <- args[3]
model_func <- amis_int_mod$build_transmission_model(
    weeks_indices, initial_infect, num_cores,
    mda_coverage_data_filename,
    vaccine_coverage_data_filename,
    coverage_data_path
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
amis_output <- AMISforInfectiousDiseases::amis(
    prevalence_map,
    # Wrap model_func so as to print the Python traceback
    # in case of a raised exception during execution of the
    # Python model. See
    # https://rstudio.github.io/reticulate/reference/py_last_error.html
    function(seeds, params, n_tims) {
        tryCatch(
            expr={model_func(seeds, params, n_tims)},
            error=error_function
        )
    },
    prior
)
