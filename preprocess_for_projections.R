library(AMISforInfectiousDiseases)
library(dplyr)
library(tidyr)
library(optparse)

kPathToWorkingDir <- Sys.getenv("TRACHOMA_AMIS_DIR")
kPathToModel <- Sys.getenv("TRACHOMA_MODEL_DIR")
kPathToMaps <- file.path(kPathToWorkingDir, "Maps")
kPathToPostAmisAnalysis <- file.path(kPathToWorkingDir, "post_AMIS_analysis")
kPathToAmisOutput <- file.path(kPathToWorkingDir, "AMIS_output")

create_directory_structure <- function(species) {
  kPathToInputParsMTP <- file.path(kPathToPostAmisAnalysis, paste0("InputPars_MTP_", species))
  if (!dir.exists(kPathToInputParsMTP)) {
    dir.create(kPathToInputParsMTP, recursive = T)
  }
}

process_all_batches <- function(id_vec) {
  sampled_params_all <- c()

  # Process each batch ID
  for (id in id_vec) {
    #### Load AMIS output
    amis_file_path <- file.path(kPathToAmisOutput, paste0("amis_output", id, ".Rdata"))
    if (!file.exists(amis_file_path)) {
      cat(sprintf("Error: AMIS output file not found for ID: %d\n", id))
      return(NULL)
    }

    load(amis_file_path) # loads amis_output
    sampled_params_all <- rbind(
      sampled_params_all,
      process_batch(id, amis_output)
    )

    if (id %% 100 == 0) {
      cat(paste0("id=", id, "; "))
    }
  }
  sampled_params_all
}

process_batch <- function(id, amis_output_data, species = "trachoma") {
  cat(sprintf("Processing batch ID: %d\n", id))

  iu_names <- rownames(amis_output_data$prevalence_map[[1]]$data)
  ess <- amis_output_data$ess

  #### Sample draws from the posterior
  num_sub_samples_posterior <- 200

  sampled_params_all <- c()
  for (iu in iu_names) {
    if (ess[which(iu_names == iu)] >= 200) {
      sampled_params <- sample_parameters(
        x = amis_output,
        n_samples = num_sub_samples_posterior,
        locations = which(iu_names == iu)
      )

      # this puts all ius in the same folder
      kPathToInputParsMTP <- file.path(kPathToPostAmisAnalysis, paste0("InputPars_MTP_", species))
      file_name <- file.path(kPathToInputParsMTP, paste0("InputPars_MTP_", iu, ".csv"))
      write.csv(sampled_params, file = file_name, row.names = F)

      sampled_params_iu <- cbind(IU_ID = iu, sampled_params)
      sampled_params_all <- rbind(sampled_params_all, sampled_params_iu)
    }
  }

  sampled_params_all
}

# Define command line options
option_list <- list(
  make_option(c("-i", "--id"),
    type = "integer",
    help = "Single ID to process. If not provided, will process all IDs"
  ),
  make_option(c("-s", "--species"),
    type = "character",
    default = "trachoma",
    help = "Species to process [default=%default]"
  ),
  make_option(c("-f", "--failed_ids"),
    type = "character",
    default = "",
    help = paste(
      "Comma-separated list of failed IDs to skip.",
      "Needs to be filled in once we know which batch-ids failed.",
      "Will be ignored when a single ID is specified."
    )
  )
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Process failed IDs if provided
failed_ids <- if (!is.null(opts$failed_ids)) {
  as.numeric(strsplit(opts$failed_ids, ",")[[1]])
} else {
  c()
}


# Early Feb: second runs
# from Igor's runs
# failed_ids = c(123,150,170,196,203,296,315,369,371,373,376,377,379,380,495,502,515,528,530,531,532,533,535,539,541)
# from Raiha's runs
# failed_ids = c(123,150,196,203,205,296,315,369,370,371,373,376,379,380,495,502,515,528,530,531,532,533,535,539,541)
# 30 Jan:  from previous runs
# failed_ids = c(29, # fitting didn't work properly
# 	       133,134, # randomly stopped?
#               123,153,196,296,315,322,374,375,376,377,378,379,380,381,495,502,515,527,530,532,533,539,541)
# ctd_ids = c(55,63,86,113,114,121,145,167,168,170,171,173,174,175,179,180,181,191,192,197,198,199,200,
# 203,205,206,208,262,276,298,326,327,362,363,370,371,372,412,435,441,446,447,471,503,505,507,508,510,528,531,535,549) # ids that were run for extra iterations
# failed ctd runs: 150,369


load(file.path(kPathToMaps, "trachoma_maps.rds"))
load(file.path(kPathToMaps, "trachoma_map_years.rds"))
# loading 'iu_task_lookup' (batches-IUs look up table for the fitting)
load(file.path(kPathToMaps, "iu_task_lookup.Rdata"))

create_directory_structure(species = opts$species)
if (!is.null(opts$id)) {
  # Process single ID
  load(file.path(kPathToAmisOutput, paste0("amis_output", opts$id, ".Rdata")))
  sampled_params_all <- process_batch(opts$id, amis_output)
  if (!is.null(sampled_params_all)) {
    output_file <- paste0("InputPars_MTP_", opts$id, ".rds")
    save(sampled_params_all, file.path(kPathToInputParsMTP, output_file))
    cat(sprintf("Saved posterior parameter samples for single batch [%d] to: %s\n", opts$id, output_file))
  }
} else {
  id_vec <- setdiff(1:max(iu_task_lookup$TaskID), c(failed_ids))
  cat(sprintf("Processing %d batch IDs\n", length(id_vec)))

  sampled_params_all <- process_all_batches(id_vec)
  output_file <- file.path(kPathToInputParsMTP, "InputPars_MTP_allIUs.rds")
  save(sampled_params_all, file = output_file)
  cat(sprintf("Saved posterior parameter samples for all batches to: %s\n", output_file))
}
