# Read libraries and fuctions  ==============================================================

source("scripts/000_libraries_and_functions.R")

# load =========================================================================
cat("load\n")

dt_cohort <- qread("input/dt_cohort.qs")

# calculate study end date =====================================================

## most recent pos date ==  pcr_Date in Wales
pcr_date <- dt_cohort %>%
  filter(!is.na(most_recent_pos_date)) %>%
  summarise(max = max(most_recent_pos_date))

vacc_date <- dt_cohort %>%
  filter(Booster_date <= today()) %>%
  summarise(max = max(Booster_date))

if (study_end_date != min(pcr_date$max, vacc_date$max)) {
  cat(
    "You need to update study end date to be\n",
    as.character(min(pcr_date$max, vacc_date$max)), "\n",
    sep = ""
  )
}

rename_vacc <- function(x) {
  case_when(
    x == "AZ"     ~ "AZ",
    x == "Moderna" ~ "MD",
    x == "BNT162b2" ~ "PB",
      )
}

## based on this updated study end date to 17 Feb

###################################################
## Clean the environment
###################################################

# Sanity check for time between doses
#z = dt_cohort_clean %>%
 # mutate(days = Dose2_date - Dose1_date) %>%
  #mutate(vacc_dose1_dose2_diff_week = floor(interval(Dose1_date, Dose2_date) / dweeks(1))) %>%
   #        select(study_id, Dose1_date, Dose2_date, days, vacc_dose1_dose2_diff_week) %>%
  #group_by(vacc_dose1_dose2_diff_week) %>%
  #summarise(
   # min = min(days, na.rm=T),
    #max = max(days, na.rm =T)
  #) %>%
  #arrange(vacc_dose1_dose2_diff_week)

# checking the COVID death indicator
#z = dt_cohort %>%
 # mutate(death_covid_date = if_else(covid_death == 1, date_of_death_rounded, NA_Date_)) %>%
#  mutate(death_noncovid_date = if_else(covid_death == 0, date_of_death_rounded, NA_Date_)) %>%
#  select(covid_death, date_of_death_rounded,death_covid_date, death_noncovid_date) %>%
 # mutate(
#  has_death_noncovid  = as.numeric(!is.na(death_noncovid_date)),
 # has_death_covid     = as.numeric(!is.na(death_covid_date))
#)


dt_cohort_clean = dt_cohort %>%
  ## select only variables we need
  dplyr::select(study_id, age, sex, settlement_band_2015, mdm_decile, Dose1_name, Dose1_date, Dose2_name,
         Dose2_date, Dose3_name, Dose3_date, Booster_name, Booster_date, date_of_death_rounded, covid_death,
         covid_admission_date,most_recent_pos_date, n_tests_gp, prior_infect_wksgrp, BNF_group, lgd) %>% collect() %>%
  mutate(
    # vaccination
    has_vacc_dose1  = as.numeric(!is.na(Dose1_date)), # all have dose 1
    has_vacc_dose2  = as.numeric(!is.na(Dose2_date)), # all have dose 2
    has_vacc_doseb  = as.numeric(!is.na(Dose3_date)), #18056 people at this stage have booster
    Dose1_name = rename_vacc(Dose1_name), # this just renames to AZ/PB/MD
    Dose2_name = rename_vacc(Dose2_name),
    Dose3_name = rename_vacc(Dose3_name),
    Booster_name = rename_vacc(Booster_name)
  ) %>%
  mutate(
    # time between 1st and 2nd dose
    vacc_dose1_dose2_diff_week = floor(interval(Dose1_date, Dose2_date) / dweeks(1)),
    vacc_dose1_dose2_diff_cat = case_when(
      has_vacc_dose1 == 0 | has_vacc_dose2 == 0   ~ "No dose 1 or dose 2",
      between(vacc_dose1_dose2_diff_week,  3,  6) ~ "03-06wk",
      between(vacc_dose1_dose2_diff_week,  7,  8) ~ "07-08wk",
      between(vacc_dose1_dose2_diff_week,  9, 10) ~ "09-10wk",
      between(vacc_dose1_dose2_diff_week, 11, 12) ~ "11-12wk",
      vacc_dose1_dose2_diff_week >= 13            ~ "13+wk"
    ) %>% factor()
    ) %>%
  mutate(
    has_covid_infection = as.numeric(!is.na(most_recent_pos_date)),
    has_covid_hosp      = as.numeric(!is.na(covid_admission_date)),
    # time since most recent infection prior to dose 2
    dose2_prior_infection_week = floor(interval(most_recent_pos_date, Dose2_date) / dweeks()),
    dose2_prior_infection_cat = case_when(
      is.na(most_recent_pos_date) | dose2_prior_infection_week <0   ~ "No prior infection", # included those with positve specimen after Dose 2 as no prior infection
      between(dose2_prior_infection_week,  0, 12) ~ "00-12wk", # 0-2 months
      between(dose2_prior_infection_week, 13, 25) ~ "13-25wk", # 3-5 months
      between(dose2_prior_infection_week, 26, 38) ~ "26-38wk", # 6-8 months
      dose2_prior_infection_week >= 39            ~ "39+wk",   # 9+ months
    ) %>% factor()
  ) %>%
  mutate(
    # death
    death_noncovid_date = if_else(covid_death == 0, date_of_death_rounded, NA_Date_),
    death_covid_date    = if_else(covid_death == 1, date_of_death_rounded, NA_Date_),
    has_death_noncovid  = as.numeric(!is.na(death_noncovid_date)),
    has_death_covid     = as.numeric(!is.na(death_covid_date))
  ) %>%
  ## ##Make sure time between dose 1 and dose 2 is OK
  mutate(days_btn_1_2 = Dose2_date - Dose1_date) %>%  # Yes - starts at 21 
  ## Now create a clean booster which is min of Dose 3 and booster doses
  mutate(B_date_copy = Booster_date) %>% #create back up for booster date / name
  mutate(B_name_copy = Booster_name) %>%
  mutate(Booster_date = pmin(Dose3_date, Booster_date, na.rm = TRUE)) %>%
  mutate(Booster_name = if_else(
    condition = !is.na(Dose3_date) & Dose3_date == Booster_date,
    true = Dose3_name,
    false = Booster_name
  )) %>%
  mutate(Dose4_date = pmax(Dose3_date, B_date_copy)) %>%
  mutate(days_btn_2_b = Booster_date - Dose2_date) %>%
  filter(days_btn_2_b >=21 | is.na(days_btn_2_b)) %>% #apply the 21 day rule
  mutate(days_btn_4_b = Dose4_date - Booster_date) %>%
  filter(days_btn_4_b >=21 | is.na(days_btn_4_b)) %>% #apply the 21 day rule
     mutate(study_end_date = study_end_date) 

dt_cohort_orig = dt_cohort # this is the dataset without the variables / further cleaning of booster

dt_cohort = dt_cohort_clean

##############################################################################
## 17 April -  realised that issues between booster and dose 2 were causing 
## problems with survival analysis so added a cleaning step to remove those with
## booster within 21 days
## Do all cleaning for vaccine doses above
## Want to have dose 1, dose 2, booster (min of dose 3 / booster),
## dose 4 (censoring; where dose 3 and dose 4 are complete)
##############################################################################

# Tidy up factors / labels =====================================================

###################################
## Align to age groups
###################################

#z = dt_cohort %>%
# group_by(age_gp) %>%
#tally()
# 
# dt_cohort <-mutate(dt_cohort,age_gp =
#                      cut(age, breaks = c(-Inf,49,54,59,64,69,74,79,Inf),
#                          labels=c("18-49","50-54",
#                                   "55-59","60-64","65-69",
#                                   "70-74","75-79","80+")))
# 
# # dt_cohort <-mutate(dt_cohort,age_gp =
# #                      cut(age, breaks = c(-Inf,29,34,39,44,49,54,59,64,69,74,79,Inf),
# #                          labels=c("18-29","30-34", "35-39","40-44","45-49","50-54",
# #                                   "55-59","60-64","65-69",
# #                                   "70-74","75-79","80+")))
# 
# dt_cohort$age_gp <- factor(dt_cohort$age_gp, ordered = FALSE,
#                                   levels = c("18-49", "50-54", "55-59","60-64","65-69",
#                                             "70-74","75-79","80+"))
dt_cohort = dt_cohort %>%
  mutate(
    # age group
    age_gp = case_when(
      age <= 17                ~ "17_and_under",
      age %>% between(18,  49) ~ "18-49",
      age %>% between(50,  54) ~ "50-54",
      age %>% between(55,  59) ~ "55-59",
      age %>% between(60,  64) ~ "60-64",
      age %>% between(65,  69) ~ "65-69",
      age %>% between(70,  74) ~ "70-74",
      age %>% between(75,  79) ~ "75-79",
      age >= 80               ~ "80+"
    ) %>%
      factor() %>%
      fct_explicit_na()
  )

label(dt_cohort$age_gp) = "Age groups"

x = dt_cohort %>% count(age_gp)

# # BNF ================================
# # if BNF is NA saying they have no risk groups
# 
# # names(dt_cohort)
# z = dt_cohort %>%
#   count(BNF_group)
# 
# dt_cohort = dt_cohort %>%
#    mutate(BNF_group = ifelse(is.na(BNF_group), 0, BNF_group))
# 
# dt_cohort$BNF_group = factor(dt_cohort$BNF_group, levels=c(0,1,2,3,4,5,6),
#                              labels=c("0","1","2","3","4","5","6+"))
# 
# label(dt_cohort$BNF_group) = "BNF Chapters"


###################################
## Create urban / rural indicator
###################################
## Taken from the NISRA defintion A-e urban,
dt_cohort<-mutate(dt_cohort,ur_combined=
                    if_else(settlement_band_2015=="A"|
                              settlement_band_2015=="B"| 
                              settlement_band_2015=="C"|
                              settlement_band_2015=="D"|
                              settlement_band_2015=="E", 1,
                            if_else(settlement_band_2015=="0",0,
                                    if_else(settlement_band_2015=="NA",9,2))))

## 1 is urban, 2 is rural, 0 wasn't linked, then NA

dt_cohort = dt_cohort %>%
  mutate(ur_combined = ifelse(is.na(ur_combined) | ur_combined == 0, 'unknown', ur_combined))


dt_cohort$ur_combined = factor(dt_cohort$ur_combined, levels=c(1,2, "unknown"),
                               labels=c("Urban", "Rural", "Unknown"))


label(dt_cohort$ur_combined) = "Urban/Rural"

###################################
## Create Derpivation indicator
## 1 is most deprived
###################################


dt_cohort<-mutate(dt_cohort, nimdm =
                    if_else(mdm_decile=="1 " | mdm_decile=="2 ",1,
                            if_else(mdm_decile=="3 " | mdm_decile=="4 ",2,
                                    if_else(mdm_decile=="5 " | mdm_decile=="6 ",3,
                                            if_else(mdm_decile=="7 " | mdm_decile=="8 ",4,
                                                    if_else(mdm_decile=="9 " | mdm_decile=="10",5,0))))))


# z = dt_cohort %>%
#   count(nimdm)

# replace NA with unknown
dt_cohort = dt_cohort %>%
  mutate(nimdm = ifelse(is.na(nimdm) | nimdm==0, 'unknown', nimdm))

dt_cohort$nimdm = factor(dt_cohort$nimdm, levels=c(1,2,3,4,5, "unknown"),
                         labels=c("1 - most deprived","2","3","4","5 - least deprived", "Unknown"))

label(dt_cohort$nimdm) = "Deprivation"

label(dt_cohort$n_tests_gp) = "Number of tests prior to dose 2"

# Sex =========================================================
dt_cohort$sex = as.factor(dt_cohort$sex)

dt_cohort$sex = factor(dt_cohort$sex, levels=c("F","M"),
                               labels=c("Female", "Male"))

label(dt_cohort$sex) = "Gender"

# Save =========================================================================
cat("Save\n")

qsave(
  dt_cohort,
  file = ("input/dt_cohort_clean.qs")
)

message("***Further cleaning and variables generated - dataset to work with for breakthrough***")
# 


