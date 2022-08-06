#######################################################
## working out count of PCR tests before second dose ##
#######################################################

## Keep only PCR tests or where test kit is NA
dt_pcr <- dt_tests[grepl("pcr", TESTKIT, ignore.case = TRUE) == TRUE | is.na(TESTKIT)]

## Pull out tests for people in the study cohort
dt_pcr <- merge(dt_cohort[, .(study_id, Dose2_date)], dt_pcr, by = "study_id", all.x = TRUE)


## drop if specimen date is missing
dt_pcr <- dt_pcr[!is.na(specimen_date)]


## Subset the dt_pcr where the spec date is less than dose 2 date 
dt_pcr_tests_before_dose_2 <- dt_pcr[specimen_date < Dose2_date, 
                                     .(pcr_tests_before_dose_2 = .N), study_id]

dt_cohort <- merge(dt_cohort, dt_pcr_tests_before_dose_2, by = "study_id", all.x = TRUE)

## replace NA w zeros
na_to_zero(dt_cohort, "pcr_tests_before_dose_2")

##################################################################
## Number of PCR tests prior to second dose - Create categorical 
##################################################################

dt_cohort = dt_cohort %>%
  mutate(n_tests_gp = cut(pcr_tests_before_dose_2, breaks=c(-1,0,1,2,4,9,max(pcr_tests_before_dose_2)), 
                          labels=c("0","1","2","3-4","5-9","10+")))


#########################################################
## identifying emergency admissions after a covid test ##
#########################################################

## Subset to emergency admissions only 
dt_emergency_admissions <- dt_admissions[admission_method == 
                                           "emergency" & normal_inpatient == TRUE]

dt_emergency_admissions[, admission_index := .I]

# limit to people currently part of the cohort at this point in the workflow

dt_emergency_admissions <- dt_emergency_admissions[study_id %in% dt_cohort$study_id]

## create a dataset for PCR positive that are in our cohort
dt_pcr_positive <- dt_pcr[result == "positive" & 
                            study_id %in% dt_cohort$study_id, .(study_id, specimen_date)]

# merge admissions and positive PCR tests, and keep only the overlap  -
# but need to add the coded ones back in later.

dt_pos_pcr_and_emergency_admissions <- merge(dt_emergency_admissions, 
                                             dt_pcr_positive, by = "study_id") 

dt_pos_pcr_and_emergency_admissions[, time_between_specimen_and_admission := as.numeric(admission_date - specimen_date)]

dt_pos_pcr_and_emergency_admissions <-
  dt_pos_pcr_and_emergency_admissions[time_between_specimen_and_admission <= 14 &
                                        time_between_specimen_and_admission >= -1 &
                                        (specimen_date <= discharge_date |
                                           is.na(discharge_date))] 

#get first admission determined thorugh this method

dt_first_admission_by_pcr <- dt_pos_pcr_and_emergency_admissions[, .SD[which.min(admission_date)], study_id]

dt_first_admission_by_icd_10 <- dt_emergency_admissions[icd10_covid == TRUE, .SD[which.min(admission_date)], study_id]

# row bind them, which introduces duplicates

dt_first_covid_admission <- rbind(dt_first_admission_by_pcr[, .(study_id, admission_date)], dt_first_admission_by_icd_10[, .(study_id, admission_date)])

# and take the first date to get one per person 

dt_first_covid_admission <- dt_first_covid_admission[, .(covid_admission_date = min(admission_date)), study_id]

dt_cohort <- merge(dt_cohort, dt_first_covid_admission, by = "study_id", all.x = TRUE)


dt_cohort <- dt_cohort[is.na(covid_admission_date) | (covid_admission_date >= booster_start_date & covid_admission_date >= Dose2_date)]


#########################################################
## identifying time since most recent positive 
#########################################################

## Subset to most recent specimen per study ID

pos_pcr_recent = dt_pcr_positive %>%
  group_by(study_id) %>%
  slice(which.max(specimen_date)) %>%
  rename(most_recent_pos_date = specimen_date)

## now link back to dt_covid

dt_cohort = dt_cohort %>%
  left_join(pos_pcr_recent, by="study_id") 

###############################################
## Days between recent pos and second dose
############################################

# Calculate days between positive test and second dose vaccination, 
# categorising any with positives AFTER dose 2 as "No prior infection":

dt_cohort <- dt_cohort %>% 
  mutate(days_covpos_dose2 = as.Date(Dose2_date)+13 - as.Date(most_recent_pos_date))

dt_cohort$days_covpos_dose2 <- as.numeric(dt_cohort$days_covpos_dose2) 

dt_cohort = dt_cohort %>%
  mutate(weeks_covpos = floor(interval(most_recent_pos_date, Dose2_date) / dweeks())) %>%
  mutate(wks_covpos_grp = case_when(
  is.na(weeks_covpos) | weeks_covpos <0    ~ "No prior infection", # those with pos after dose 2 date are recorded as no prior infection
  between(weeks_covpos,  0, 12) ~ "00-12wk", # 0-2 months
  between(weeks_covpos, 13, 25) ~ "13-25wk", # 3-5 months
  between(weeks_covpos, 26, 38) ~ "26-38wk", # 6-8 months
  weeks_covpos >= 39            ~ "39+wk",   # 9+ months
  ) %>% factor())

names(dt_cohort)
names(dt_cohort)[names(dt_cohort) == 'weeks_covpos'] <- 'prior_infect_weeks'
names(dt_cohort)[names(dt_cohort) == 'wks_covpos_grp'] <- 'prior_infect_wksgrp'

############################################
## Days between recent Dose 1 and 2
############################################

dt_cohort = dt_cohort %>%
  mutate(days_btn_1_2 = as.numeric(Dose2_date - Dose1_date)) 


###################################################
## Bring in EPD ##
###################################################

## drop sex from epd to avoid conflict

dt_epd = dt_epd %>%
  dplyr::select(-sex)

dt_cohort <- merge(dt_cohort, dt_epd, by = "study_id", all.x = TRUE) 

## Create BNF categories

dt_cohort = dt_cohort %>%
  mutate(BNF_group = ifelse(BNF_chapters>=6,6,BNF_chapters))

#Check this worked - yes
#z = dt_cohort %>%
 #group_by(BNF_group) %>%
#  summarise(
 #  min = min(BNF_chapters, na.rm=T), 
#  max = max(BNF_chapters, na.rm=T)
# ) %>%
#  arrange(BNF_group) 
 
dt_cohort$BNF_group = factor(dt_cohort$BNF_group, levels=c(0,1,2,3,4,5,6),
                             labels=c("0","1","2","3","4","5","6+"))

label(dt_cohort$BNF_group) = "BNF Chapters"

#z = dt_cohort%>%
 # group_by(BNF_group) %>%

# Save =========================================================================
cat("Save\n")

qsave(
  dt_pcr_positive,
  file = ("input/d_cohort_raw.qs")
)

message("***At this point have vaccine data linked to NHAIS (deprivation,U/R), deaths, emergency admissions,
        tests, EPD. Can remove other datasets. Rule still to review = days between dose 1 and 2***")


