library(dplyr)
library(tidyr)
library(optparse)

kPathToWorkingDir <- Sys.getenv("TRACHOMA_AMIS_DIR")
kPathToInputs <- Sys.getenv("PATH_TO_FITTING_PREP_INPUTS")
kPathToArtefacts <- Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS")
kPathToModelOutput <- file.path(kPathToArtefacts, "model_output")
kPathToMaps <- file.path(kPathToArtefacts, "Maps")

other_prev_categories <- c("Endemicity unknown", "Not Suspected", "Suspected Endemic") # ignore these
post_mda_surv <- "Under post-MDA surveillance" # set to <5% tf_prevalence

# Read in data
trachoma_Geoconnect_to_IUID <- read.csv(file.path(kPathToInputs, "trachoma_IU_Match_PB.csv")) # PB's lookup table fo mapping Geoconnect_ID to IU_ID
trachoma_data_raw <- read.csv(file.path(kPathToInputs, "trachomaComb_IU.csv")) %>%
  left_join(trachoma_Geoconnect_to_IUID %>% select(Geoconnect_ID, IU_ID), by = "Geoconnect_ID")

# there are 10 Geoconnect IDs that could not be mapped to IU_ID
# View(unique(trachoma_data_raw$Geoconnect_ID[which(is.na(trachoma_data_raw$IU_ID))]))

# convert tf_prevalence to numerical value
# for cases where many Geoconnnect IDs map to a single IU_ID:
# - TODO: want to take population weighted mode tf_prevalence (i.e. prevalence interval that has the most people)
# - if MDA occurred in any of the Geoconnect_IDs, assume it occurred in the whole IU_ID
trachoma_espen_data <- trachoma_data_raw %>%
  filter(!is.na(trachoma_data_raw$IU_ID)) %>%
  mutate(tf_prevalence = ifelse(tf_prevalence %in% other_prev_categories, NA, tf_prevalence)) %>%
  mutate(tf_prevalence = ifelse(tf_prevalence %in% post_mda_surv, "<5%", tf_prevalence)) %>%
  mutate(survey_n = ifelse(is.na(tf_prevalence), NA, survey_n)) %>%
  mutate(tf_prev_upr = case_when(
    tf_prevalence == "<5%" ~ 0.05,
    tf_prevalence == "5-9.9%" ~ 0.099,
    tf_prevalence == "10-29.9%" ~ 0.299,
    tf_prevalence == "30-49.9%" ~ 0.499,
    tf_prevalence == ">=50%" ~ 1
  ))

save(trachoma_espen_data, file = file.path(kPathToMaps, "trachoma_espen_data.Rdata"))

# select relevant columns
# if type_prevalence changes OR survey_n changes
trachoma_data <- trachoma_espen_data %>%
  group_by(IU_ID) %>%
  mutate(survey_n_clean = ifelse(is.na(survey_n), -1, survey_n)) %>%
  mutate(survey_change = ifelse((survey_n_clean == lag(survey_n_clean) & !(type_prevalence == lag(type_prevalence))), 1, 0)) %>%
  mutate(survey_change = ifelse(is.na(survey_change), 0, survey_change)) %>%
  ungroup() %>%
  group_by(IU_ID, Year, ADMIN0ISO2) %>%
  summarise(
    max_tf_prev_upr = ifelse(length(which(!is.na(tf_prev_upr))) > 0, max(tf_prev_upr, na.rm = T), NA),
    max_survey_n = ifelse(length(which(!is.na(tf_prev_upr))) > 0, max(survey_n, na.rm = T), NA),
    max_survey_change = ifelse(length(which(!is.na(tf_prev_upr))) > 0, max(survey_change, na.rm = T), NA),
    PC_in_group = max(PC)
  ) %>%
  ungroup()

save(trachoma_data, file = file.path(kPathToMaps, "trachoma_data.Rdata"))


# # look at dups (many Geoconnect_ID to ESPEN IU)
# dups = trachoma_espen_data %>%
#   select(IU_ID,Year,tf_prevalence,Geoconnect_ID) %>%
#   group_by(IU_ID,Year,tf_prevalence) %>%
#   tally() %>%
#   filter(n>1) %>%
#   ungroup() %>%
#   group_by(IU_ID,Year) %>%
#   tally() %>%
#   filter(n>1)

# likelihood for trachoma model
# ignore if 'Endemicity unknown, Not Suspected, Suspected Endemic
trachoma_likelihood <- function(data, sim_prev, log) {
  # assuming 250 individuals
  n <- 250

  # tf bands (range of positive #): <5% (0-12), 5-9.9% (13-24) , 10-29.9% (25-74), 30-49.9% (75-124), >=50% (125+)
  # probability of being in band
  prob_0_0.05 <- pbinom(12, n, sim_prev)
  prob_0.05_0.1 <- pbinom(24, n, sim_prev) - pbinom(12, n, sim_prev)
  prob_0.1_0.3 <- pbinom(74, n, sim_prev) - pbinom(24, n, sim_prev)
  prob_0.3_0.5 <- pbinom(124, n, sim_prev) - pbinom(74, n, sim_prev)
  prob_0.5_1 <- pbinom(250, n, sim_prev) - pbinom(124, n, sim_prev)

  if (log) {
    lh <- case_when(
      data == 0.05 ~ log(prob_0_0.05),
      data == 0.099 ~ log(prob_0.05_0.1),
      data == 0.299 ~ log(prob_0.1_0.3),
      data == 0.499 ~ log(prob_0.3_0.5),
      data == 1 ~ log(prob_0.5_1)
    )
  } else {
    lh <- case_when(
      data == 0.05 ~ prob_0_0.05,
      data == 0.099 ~ prob_0.05_0.1,
      data == 0.299 ~ prob_0.1_0.3,
      data == 0.499 ~ prob_0.3_0.5,
      data == 1 ~ prob_0.5_1
    )
  }
  return(lh)
}

# # look at data
# summary = trachoma_data %>%
#   filter(!is.na(tf_prevalence)) %>%
#   group_by(Year,tf_prevalence) %>%
#   tally()
# summary$tf_prevalence = factor(summary$tf_prevalence, levels=c("<5%","5-9.9%","10-29.9%","30-49.9%",">=50%"))
# library(ggplot2)
# ggplot(data=summary) +
#   geom_bar(aes(x=tf_prevalence,y=n),stat="identity") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   facet_wrap(~Year)

# fit only to years with new data (based on when max_survey_n_clean changes)
new_surveys_completed <- trachoma_data %>%
  arrange(IU_ID, Year) %>%
  group_by(IU_ID) %>%
  mutate(max_survey_n_clean = ifelse(is.na(max_survey_n), -1, max_survey_n)) %>%
  mutate(new_data = (!(max_survey_n_clean == lag(max_survey_n_clean))) | (Year == 1996 & max_survey_n == 0) | max_survey_change == 1) %>%
  filter(new_data == TRUE)

# MDA histories
trachoma_histories <- trachoma_data %>%
  select(IU_ID, Year, PC_in_group, ADMIN0ISO2) %>%
  filter(IU_ID %in% new_surveys_completed$IU_ID) %>%
  pivot_wider(names_from = Year, values_from = PC_in_group)

# #I don't think this needs to be done for trachoma
# Overwrite MDA files for IUs with no treatments to add negligible coverage (otherwise python code will break)
# iu_no_mda = unique((trachoma_histories %>%
#                       filter(IU_ID %in% (iu_task_lookup %>%
#                                                    filter(TaskID == 1))$IU_ID))$IU_ID)
# inds = which(trachoma_histories$IU_ID %in% iu_no_mda)
# trachoma_histories[inds, colnames(trachoma_histories) %in% c("1996","1997")] = 1e-11
#

# find unique histories
unique_history <- trachoma_histories %>%
  arrange(IU_ID) %>%
  group_by(across(c(-IU_ID))) %>%
  summarise(ius = list(IU_ID)) %>%
  ungroup() %>%
  mutate(index = row_number())

# initialise iu_task_lookup
iu_task_lookup <- lapply(1:nrow(unique_history), function(i) data.frame(unique_history[["ius"]][[i]], unique_history[["index"]][i]))
iu_task_lookup <- do.call(rbind, iu_task_lookup)
colnames(iu_task_lookup) <- c("IU_ID", "TaskID")
save(iu_task_lookup, file = file.path(kPathToMaps, "iu_task_lookup.Rdata"))

# get task IDs and unique survey years in the data
ius_vector <- iu_task_lookup$IU_ID
years_vector <- sort(unique(new_surveys_completed$Year))

# template to append data to
trachoma_all <- data.frame(
  IU_ID = rep(ius_vector, each = length(years_vector)),
  Year = rep(years_vector, length(ius_vector))
) %>%
  left_join(iu_task_lookup) %>%
  left_join(new_surveys_completed %>% select(IU_ID, Year, max_tf_prev_upr))

# save maps
trachoma_maps <- lapply(years_vector, function(year) {
  list(
    data = trachoma_all %>%
      filter(Year == year) %>%
      select(IU_ID, TaskID, max_tf_prev_upr) %>%
      arrange(IU_ID),
    likelihood = trachoma_likelihood
  )
})
save(years_vector, file = file.path(kPathToMaps, "trachoma_map_years.rds"))
save(trachoma_maps, file = file.path(kPathToMaps, "trachoma_maps.rds"))

# save MDA inputs for python model

# Get into format for python model
coverage_wide <- trachoma_histories %>%
  mutate(
    "Country/Region" = "All",
    "Intervention Type" = "Treatment",
    "Platform Type" = "Campaign",
    Platform = "MDA",
    # Drug = MDA_scheme,
    "Cohort (if not total pop in country/region)" = NA,
    "min age" = 1,
    "max age" = 100
  ) %>%
  select(
    IU_ID, "Country/Region", "Intervention Type", "Platform Type", Platform, # Drug,
    "Cohort (if not total pop in country/region)", "min age", "max age",
    starts_with(c("19", "20"))
  )

# export histories files
kPathToEndgameInputs <- file.path(kPathToArtefacts, "trachoma", "data", "coverage", "endgame_inputs")
if (!dir.exists(kPathToEndgameInputs)) {
  dir.create(kPathToEndgameInputs, recursive = T)
}

process_all_batches <- function(id_vec) {
  num_batches <- max(iu_task_lookup$TaskID)
  cat(sprintf("num_batches: %d\n", num_batches))

  for (id in id_vec) {
    process_batch(id)
  }
}

process_batch <- function(id) {
  iu_file <- file.path(kPathToEndgameInputs, paste0("IUs_MTP_", id, ".csv"))
  ius_per_batch <- iu_task_lookup %>%
    filter(TaskID == id)
  ius <- matrix(ius_per_batch$IU_ID, ncol = 1)
  write.table(ius, file = iu_file, row.names = F, col.names = F, quote = F, sep = ",") # write input parameter file

  for (iu in ius[1]) {
    mda_cov_iu <- coverage_wide %>%
      filter(IU_ID == iu) %>%
      select(-IU_ID)

    mda_cov_iu_xl <- rbind(colnames(mda_cov_iu), mda_cov_iu) %>%
      mutate_at(vars(starts_with(c("19", "20"))), as.numeric)

    # The first year columns must be a numeric integer for it to be read correctly into python!!
    mda_path <- file.path(kPathToEndgameInputs, paste0("InputMDA_MTP_", id, ".csv"))
    write.table(mda_cov_iu_xl, file = mda_path, col.names = F, row.names = F, sep = ",") # write input MDA file
  }
}

option_list <- list(
  make_option(c("-i", "--id"),
    type = "integer",
    help = "Single ID to process. If not provided, will process all IDs"
  )
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)


if (!is.null(opts$id)) {
  cat(sprintf("Processing single batch ID: %d\n", opts$id))
  process_batch(opts$id)
} else {
  id_vec <- 1:max(iu_task_lookup$TaskID)
  cat("Processing all batch IDs\n")
  process_all_batches(id_vec)
}
