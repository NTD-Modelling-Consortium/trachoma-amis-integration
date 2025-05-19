
setwd("./Maps")
library(dplyr)
library(tidyr)

load("trachoma_maps.rds")
load("trachoma_data.Rdata")
load("iu_task_lookup.Rdata")

# fit only to years with new data (based on when max_survey_n_clean changes)
new_surveys_completed=trachoma_data %>% 
  arrange(IU_ID,Year) %>%
  group_by(IU_ID) %>%
  mutate(max_survey_n_clean = ifelse(is.na(max_survey_n),-1,max_survey_n)) %>%
  mutate(new_data = (!(max_survey_n_clean == lag(max_survey_n_clean))) | (Year==1996 & max_survey_n==0) | max_survey_change == 1) %>%
  filter(new_data == TRUE)

# Repeat most recent history if any treatment was given in the previous 3 years
recent_treatments_since_2019 = trachoma_data %>%
  filter(PC_in_group == 1) %>%
  group_by(IU_ID) %>%
  summarise(last_treated_year = floor(max(Year))) %>%
  filter(last_treated_year >= 2019) %>%
  left_join(trachoma_data, by=c("IU_ID","last_treated_year"="Year")) %>%
  select(IU_ID, last_treated_year, PC_in_group, ADMIN0ISO2)

recent_treatments_future_years = cbind(Year = c(rep(2022:2025, each=nrow(recent_treatments_since_2019))),
                                       rbind(recent_treatments_since_2019,recent_treatments_since_2019,recent_treatments_since_2019,recent_treatments_since_2019)) %>%
  select(IU_ID,Year,ADMIN0ISO2,PC_in_group)

# add years with no MDA
ius_no_recent_treatments = unique(trachoma_data$IU_ID)[which(!unique(trachoma_data$IU_ID) %in% recent_treatments_since_2019$IU_ID)]
ius_no_recent_treatments_data = unique(trachoma_data %>% 
  filter(IU_ID %in% ius_no_recent_treatments) %>%
  select(IU_ID,ADMIN0ISO2))
no_recent_treatments_future_years = cbind(IU_ID = rep(ius_no_recent_treatments_data$IU_ID, 4),
                                          Year = c(rep(2022:2025, each=length(ius_no_recent_treatments))),
                                          ADMIN0ISO2 = rep(ius_no_recent_treatments_data$ADMIN0ISO2, 4),
                                          PC_in_group = 0)

# MDA histories
trachoma_data_joined = rbind(trachoma_data %>% select(IU_ID,Year,ADMIN0ISO2,PC_in_group),recent_treatments_future_years,no_recent_treatments_future_years)
projections_histories = trachoma_data_joined %>%
  filter(IU_ID %in% new_surveys_completed$IU_ID) %>%
  pivot_wider(names_from=Year, values_from=PC_in_group)

write.csv(trachoma_data_joined, file=paste0("mda_history_trachoma.csv"))


# save MDA inputs for python model

# Get into format for python model
coverage_wide = projections_histories %>%
  mutate("Country/Region" = "All",
         "Intervention Type"= "Treatment",
         "Platform Type" = "Campaign",
         Platform = "MDA",
         #Drug = MDA_scheme,
         "Cohort (if not total pop in country/region)" = NA,
         "min age" = 1,
         "max age" = 100) %>%
  select(IU_ID,"Country/Region","Intervention Type","Platform Type", Platform, #Drug,
         "Cohort (if not total pop in country/region)","min age","max age",
         starts_with(c("19","20")))

# export histories files
inputs_path = "../ntd-model-trachoma/trachoma/data/coverage/endgame_inputs"
if (!dir.exists(inputs_path)) {dir.create(inputs_path)}

for (id in 1:max(iu_task_lookup$TaskID)){
  
  iu_file<-file.path(inputs_path,paste0("IUs_MTP_",id,".csv"))
  ius_per_batch = iu_task_lookup %>%
    filter(TaskID == id)
  ius = matrix(ius_per_batch$IU_ID,ncol=1)
  write.table(ius,file=iu_file, row.names=F, col.names = F, quote=F, sep=",")# write input parameter file
  
  for (iu in ius) {
    
    mda_cov_iu = coverage_wide %>%
      filter(IU_ID == iu) %>%
      select(-IU_ID) 
    
    mda_cov_iu_xl = rbind(colnames(mda_cov_iu),mda_cov_iu) %>% 
      mutate_at(vars(starts_with(c("19","20"))), as.numeric)
    
    # The first year columns must be a numeric integer for it to be read correctly into python!!
    mda_path<-file.path(inputs_path,paste0("InputMDA_MTP_projections_",iu,".csv"))
    write.table(mda_cov_iu_xl, file=mda_path, col.names=FALSE, row.names=FALSE, sep=",") # write input MDA file
  }
}

# make look up table for country codes
table_iu_idx <- trachoma_maps[[1]]$data[,c("IU_ID","TaskID")]
colnames(table_iu_idx) <- c("IU_ID","TaskID")
rownames(table_iu_idx) <- NULL
library(sf)
df <- read_sf(dsn = "../ESPEN_IU_2021/", layer = "ESPEN_IU_2021") %>%
  filter(IU_ID %in% table_iu_idx$IU_ID)
df <- st_drop_geometry(df[,c("IU_ID","ADMIN0ISO3")])
head(df)
table_iu_idx$country <- NA
for(i in 1:nrow(table_iu_idx)){
  IU <- table_iu_idx[i,"IU_ID"]
  wh <- which(df$IU_ID==IU)[1]
  country <- df[wh,"ADMIN0ISO3"]
  table_iu_idx[i,"country"] <- country
}
table_iu_idx$IU_ID <- as.character(table_iu_idx$IU_ID)
table_iu_idx$TaskID <- as.integer(table_iu_idx$TaskID)
table_iu_idx$country <- as.character(table_iu_idx$country)

# join column that gives last survey year 
last_survey_year = new_surveys_completed %>%
  mutate(IU_ID = as.character(IU_ID)) %>%
  group_by(IU_ID) %>%
  summarise(stop_importation_year = max(Year))

table_iu_idx = table_iu_idx %>%
  left_join(last_survey_year)

write.csv(table_iu_idx,file=paste0("table_iu_idx_trachoma.csv"),row.names=F)




