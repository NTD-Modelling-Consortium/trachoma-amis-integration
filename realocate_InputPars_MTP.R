#### Sample draws from the AMIS output
library("AMISforInfectiousDiseases")
library("dplyr")
library("optparse")

kPathToWorkingDir <- Sys.getenv("TRACHOMA_AMIS_DIR")
kPathToModel <- Sys.getenv("TRACHOMA_MODEL_DIR")
kPathToMaps <- file.path(kPathToWorkingDir, "Maps")
kPathToPostAmisAnalysis <- file.path(kPathToWorkingDir, "post_AMIS_analysis")
kPathToAmisOutput <- file.path(kPathToWorkingDir, "AMIS_output")

set.seed(123)

randomWalk <- TRUE
yearschange_index <- c(2000, 2010, 2020) - 1926 + 1 # 1926 is data start date (1996) minus 70 year burn in.
n_timschange <- length(yearschange_index)

# Early Feb: second runs
# from Igor's runs
# failed_ids = c(123,150,170,196,203,296,315,369,371,373,376,377,379,380,495,502,515,528,530,531,532,533,535,539,541)
# from Raiha's runs
# failed_ids = c(123,150,196,203,205,296,315,369,370,371,373,376,379,380,495,502,515,528,530,531,532,533,535,539,541)

# 30 Jan:  from first runs
# failed_ids = c(29, # fitting didn't work properly
#              133,134, # randomly stopped?
#               123,153,196,296,315,322,374,375,376,377,378,379,380,381,495,502,515,527,530,532,533,539,541)
# ctd_ids = c(55,63,86,113,114,121,145,167,168,170,171,173,174,175,179,180,181,191,192,197,198,199,200,
# 203,205,206,208,262,276,298,326,327,362,363,370,371,372,412,435,441,446,447,471,503,505,507,508,510,528,531,535,549) # ids that were run for extra iterations
# failed ctd runs: 150,369

kSpeciesAll <- c("trachoma")

create_directory_structure <- function(countries, species_list) {
  kPathToModelProjections <- file.path(kPathToModel, "projections")

  for (s in species_list) {
    for (country in countries) {
      path_country <- file.path(kPathToModelProjections, s, folder_id, country)
      if (!dir.exists(path_country)) {
        dir.create(path_country, recursive = T)
      }
    }
  }
}

process_all_batches <- function(df_iu_country, failed_ids, folder_id) {
  # loading 'iu_task_lookup' (batches-IUs look up table for the fitting)
  load(file.path(kPathToMaps, "iu_task_lookup.Rdata"))
  num_batches <- max(iu_task_lookup$TaskID)
  cat(paste0("num_batches: ", num_batches, " \n"))

  for (species in kSpeciesAll) {
    cat(paste0("species: ", species, " \n"))

    for (id in setdiff(1:num_batches, failed_ids)) {
      # loads amis_output
      load(file.path(kPathToAmisOutput, paste0("amis_output", id, ".Rdata")))
      process_batch(id, amis_output, df_iu_country, folder_id)

      if (id %% 100 == 0) {
        cat(paste0("id=", id, "; "))
      }
    }

    cat(paste0("Samples realocated for all IUs in InputPars_MTP_", species, "/ \n"))
  }
}

process_batch <- function(id, amis_output_data, df_iu_country, folder_id, species = "trachoma") {
  kPathToModelProjections <- file.path(kPathToModel, "projections")
  kPathToFolderId <- file.path(kPathToModelProjections, species, folder_id)

  iu_names <- rownames(amis_output$prevalence_map[[1]]$data)
  ess <- amis_output$ess
  num_samples <- 200

  for (iu in iu_names) {
    if (ess[which(iu_names == iu)] >= 200) {
      wh <- which(df_iu_country$IU_ID == iu)
      if (length(wh) != 1) {
        stop("iu must be found exactly once in df_iu_country")
      }
      country <- df_iu_country[wh, "country"]

      iu0 <- sprintf("%05d", as.integer(iu))
      path_iu <- file.path(kPathToFolderId, country, paste0(country, iu0))
      if (!dir.exists(path_iu)) {
        dir.create(path_iu, recursive = T)
      }

      file_name_old <- file.path(kPathToPostAmisAnalysis, paste0("InputPars_MTP_", iu, ".csv"))
      sampled_params <- read.csv(file_name_old)

      sampled_params <- sampled_params[, c("seed", "beta_init", paste0("beta", yearschange_index), "eff_cov", "k_parameter")]


      sampled_params_full <- matrix(NA, ncol = 102, nrow = nrow(sampled_params))
      colnames(sampled_params_full) <- c("seed", paste0("beta", 1926:2024), "eff_cov", "k_parameter")


      if (randomWalk == T) {
        sampled_params_full[, 1] <- sampled_params[, 1]
        sampled_params_full[, 101] <- sampled_params[, 6]
        sampled_params_full[, 102] <- sampled_params[, 7]

        n_reps_1 <- yearschange_index[1] - 1
        sampled_params_full[, 2:(n_reps_1 + 1)] <- matrix(rep(sampled_params[, 2], n_reps_1), ncol = n_reps_1)

        n_reps_2 <- yearschange_index[2] - yearschange_index[1]
        sampled_params_full[, (yearschange_index[1] + 1):(yearschange_index[2])] <- matrix(rep(sampled_params[, 2] * sampled_params[, 3], n_reps_2), ncol = n_reps_2)

        n_reps_3 <- yearschange_index[3] - yearschange_index[2]
        sampled_params_full[, (yearschange_index[2] + 1):(yearschange_index[3])] <- matrix(rep(sampled_params[, 2] * sampled_params[, 3] * sampled_params[, 4], n_reps_3), ncol = n_reps_3)

        n_reps_4 <- (2026 - 1926) - yearschange_index[3]
        sampled_params_full[, (yearschange_index[3] + 1):(2026 - 1926)] <- matrix(rep(sampled_params[, 2] * sampled_params[, 3] * sampled_params[, 4] * sampled_params[, 5], n_reps_4), ncol = n_reps_4)
      } else {
        sampled_params_full[, 1] <- sampled_params[, 1]
        sampled_params_full[, 101] <- sampled_params[, 3]
        sampled_params_full[, 102] <- sampled_params[, 4]

        n_reps <- (2026 - 1926) - 1
        sampled_params_full[, 2:(2026 - 1926)] <- matrix(rep(sampled_params[, 2], n_reps), ncol = n_reps)
      }

      file_name_new <- file.path(path_iu, paste0("InputBet_", country, iu0, ".csv"))
      write.csv(sampled_params_full, file = file_name_new, row.names = F)
    }
  }
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
  make_option(c("-f", "--failed-ids"),
    type = "character",
    default = "",
    help = paste(
      "Comma-separated list of failed IDs to skip.",
      "Needs to be filled in once we know which batch-ids failed.",
      "Will be ignored when a single ID is specified."
    )
  ),
  make_option(c("-d", "--folder-id"),
    type = "character",
    help = "Example: 'source-data-20250220'."
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

df_iu_country <- read.csv(file.path(kPathToMaps, "table_iu_idx_trachoma.csv"))
countries <- sort(unique(df_iu_country$country))

create_directory_structure(countries, kSpeciesAll)

if (!is.null(opts$id)) {
  # Process single ID
  cat(sprintf("Processing single batch ID: %d\n", opts$id))
  load(file.path(kPathToAmisOutput, paste("amis_output", id, ".Rdata")))
  process_batch(opts$id, amis_output, df_iu_country, opts$folder_id)
} else {
  cat("Processing all batch IDs\n")
  process_all_batches(df_iu_country, failed_ids, opts$folder_id)
}
