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
    return(matrix(c(3, 0.04), ncol = 2, nrow = num_samples, byrow = TRUE))
}

prior <- list("dprior" = density_function, "rprior" = rnd_function)

#' Fit to initial and last week (in this case 100 years)
weeks_indices <- c(0L, 5200L)
initial_infect = 0.01
num_cores = 1L

model_func <- amis_int_mod$build_transmission_model(
    fitting_points, initial_infect, num_cores
)

amis_output <- AMISforInfectiousDiseases::amis(
    prevalence_map, model_func, prior
)