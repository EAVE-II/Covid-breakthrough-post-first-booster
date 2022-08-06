## Set parameters for age and booster 
cohort_min_age = 18
min_days_between_first_and_second_dose = 21 
permitted_booster_doses = c(NA, "BNT162b2", "Moderna")

dt_cohort <- merge(dt_demographics, dt_covid_vaccines_wide, by = "study_id", all.x = TRUE)

## Subset the population to those >= 18
dt_cohort <- dt_cohort %>%
  filter(age >= cohort_min_age)

# include only those with dose 2

dt_cohort <- dt_cohort %>%
  filter(!is.na(Dose2_date))

## Restrict the data to those where primary doses are the same - to account for 
## any potential transcription errors

dt_cohort <- dt_cohort %>%
  filter(Dose1_name == Dose2_name)

## Keep only those where the difference in dose 1 and 2 is greater than 21 days

#dt_cohort <- dt_cohort[Dose2_date - Dose1_date >= min_days_between_first_and_second_dose]

## And keep only those with an mRNA booster
dt_cohort <- dt_cohort %>%
  filter(Booster_name %in% permitted_booster_doses)

# merge deaths and exclude those who died before study start date - need to check how the rounding of date of death affects those near the study start date 

dt_cohort <- dt_cohort %>%
  left_join(dt_deaths, by = "study_id", all.x = TRUE)

## Keep those were date of death rounded is missing or is >= study date

dt_cohort <- dt_cohort %>%
  filter(is.na(date_of_death_rounded) | date_of_death_rounded >= booster_start_date) %>%
  collect()

message("****Cohort rules applied, except 21 day exclusions - link to admissions and tests 102*****")
