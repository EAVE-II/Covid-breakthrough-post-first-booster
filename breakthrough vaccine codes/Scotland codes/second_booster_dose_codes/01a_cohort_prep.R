##########################################################
# Name of file: 01a_cohort_prep.R
# Data release (if applicable):
# Original author(s): Utkarsh Agrawal
# Original date: 7 Sep 2021
# Type of script: Descriptive stats
# Version of R that the script was most recently run on: R 3.5.1
# Description of content: reads in the cohort and merges in vaccination data 
# Approximate run time: Unknown
##########################################################

# 01 Setup ####
#Libraries
library(plyr)
library(tidyverse)
library(survival)
library(lubridate)
#Load data

Location <- "/conf/"  # Server
#Location <- "//isdsf00d03/"  # Desktop

project_path <- paste0(Location,"EAVE/GPanalysis/progs/UA/second_booster_dose_failures")
project_path_vaccine <- paste0(Location,"EAVE/GPanalysis/progs/CR/Vaccine")

EAVE_cohort <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Times2021-07-28.rds"))
EAVE_cohort <- filter(EAVE_cohort, !duplicated(EAVE_LINKNO))

table(EAVE_cohort$death_covid, is.na(EAVE_cohort$NRS.Date.Death), exclude=NULL)
table(EAVE_cohort$icu_death, is.na(EAVE_cohort$date_icu_death), exclude=NULL)
table(EAVE_cohort$hosp_covid, is.na(EAVE_cohort$date_hosp_covid), exclude=NULL)


#a_begin <- as.Date("2021-04-05")  #second dose
a_begin <- as.Date("2020-09-13")  #beginning of vaccination
#remove all who have died before the beginning
EAVE_cohort <- filter(EAVE_cohort, is.na(NRS.Date.Death) | (!is.na(NRS.Date.Death) & NRS.Date.Death > a_begin))
EAVE_cohort <- mutate(EAVE_cohort, ageYear = ageYear +1) # add 1 year to get age at march 2021
#remove under 16s
EAVE_cohort <- filter(EAVE_cohort, ageYear >= 16)
EAVE_cohort <- EAVE_cohort %>% mutate(age_gp = cut(ageYear, breaks=c(-1,seq(19,89,by=5),120), 
                                   labels=c(paste(c(16,seq(20,85, by=5)),seq(19,89, by=5),sep="-"),"90+")))
  

EAVE_Weights <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))
EAVE_cohort  <- EAVE_cohort %>% left_join(EAVE_Weights, by="EAVE_LINKNO")
EAVE_cohort$eave_weight[is.na(EAVE_cohort$eave_weight)] <- mean(EAVE_cohort$eave_weight, na.rm=T)

#adjust inconsistencies in the endpoints and times - all hosp have an admission date
z_max_date_death <- max(EAVE_cohort$NRS.Date.Death, na.rm=T)
z_max_date_icu <- max(EAVE_cohort$date_icu_death, na.rm=T)
#EAVE_cohort <- EAVE_cohort %>% mutate(NRS.Date.Death = case_when(death_covid==1 & is.na(NRS.Date.Death) ~ SpecimenDate + 21,
#                                                       TRUE ~ NRS.Date.Death),
#                            date_icu_death = case_when(icu_death==1 & is.na(date_icu_death) ~ SpecimenDate + 14,
#                                                       TRUE ~ date_icu_death ) ) %>% 
#  mutate(NRS.Date.Death = case_when(NRS.Date.Death > z_max_date_death  ~ z_max_date_death,
#                                    TRUE ~ NRS.Date.Death),
#         date_icu_death = case_when(date_icu_death > z_max_date_icu ~ z_max_date_icu,
#                                    TRUE ~ date_icu_death ) )
EAVE_cohort <- EAVE_cohort %>% mutate(death_covid = case_when(death_covid==1 & is.na(NRS.Date.Death) ~ 0,
                                                       TRUE ~ death_covid),
                            icu_death = case_when(icu_death==1 & is.na(date_icu_death) ~ 0,
                                                       TRUE ~ icu_death ) )


rg <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE.rds"))
rg <- filter(rg, !duplicated(EAVE_LINKNO))
rg <- rg %>% dplyr::select(EAVE_LINKNO:EAVE_CHRONIC_LIVER_DIS, EAVE_CHRONIC_LIVER_DIS:EAVE_DIABETES, 
                           EAVE_HYPERTENSION, n_risk_gps)
z <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds"))
z <- filter(z, !duplicated(EAVE_LINKNO))
z <- z %>% dplyr::select(EAVE_LINKNO, EAVE_Smoking_Status_Worst, EAVE_BP) %>% 
  dplyr::rename(EAVE_Smoke = EAVE_Smoking_Status_Worst)

rg <- rg %>% left_join(z, by="EAVE_LINKNO")

#read in the GP vaccination data
#source("/conf/EAVE/GPanalysis/progs/CR/00_Read_GP_Vaccinations.R")
source("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/00_Read_GP_Vaccinations.R")
#source("/conf/EAVE/GPanalysis/progs/UA/second_dose_failure/00_Read_GP_Vaccinations.R")

#get all positive tests before a_begin
cdw_full  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/CDW_full.rds"))
cdw_full <- cdw_full %>% mutate(date_ecoss_specimen = as_date(date_ecoss_specimen))
z <- cdw_full %>% filter(test_result=="POSITIVE") %>% 
  dplyr::select(EAVE_LINKNO, date_ecoss_specimen) %>% 
  arrange(EAVE_LINKNO, date_ecoss_specimen) %>%
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen)))  #get one positive test per person per day
Positive_Tests <- z #can have duplicate values here
# uncomment when working with infections
#Positive_Tests <- arrange(Positive_Tests,desc(date_ecoss_specimen))
Positive_Tests <- filter(Positive_Tests, !duplicated(EAVE_LINKNO))

#get covid hospitalisations in Study period  - 14 days of a positive test
#first positive test in study period
all_hospitalisations  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/automated_any_hospitalisation_post_01022020.rds"))
summary(all_hospitalisations)
z <- Positive_Tests %>% filter(date_ecoss_specimen >= a_begin) %>% 
  arrange(EAVE_LINKNO, date_ecoss_specimen) %>% 
  filter(!duplicated(EAVE_LINKNO))
#all hospitalisations in study period
z_h <- all_hospitalisations %>% filter(is.na(discharge_date) | discharge_date > a_begin) %>% 
  dplyr::select(-validchi)
z <- z %>% left_join(z_h, by="EAVE_LINKNO") %>% 
  filter(!is.na(admission_date)) #drop those with no admissions
z <- z %>%  mutate(days = as.numeric(admission_date - date_ecoss_specimen)) %>% 
  mutate(event = case_when(days >= -2 & days <= 14 & !is.na(discharge_date) & discharge_date >= date_ecoss_specimen ~ 1,
                           days >= -2 & days <= 14 & is.na(discharge_date)  ~ 1, 
                           TRUE ~ 0)) %>% 
  filter(event==1) # 
z <- z %>% arrange(EAVE_LINKNO, days) %>% filter(!duplicated(EAVE_LINKNO))
z<-filter(z,admission_date<="2022-03-28")
covid_hospitalisations <- z

z <- covid_hospitalisations %>% group_by(admission_date) %>% dplyr::summarise(N=n())
z %>% ggplot(aes(x=admission_date, y=N)) + geom_point()+labs(title="covid hospitalisations")
summary(covid_hospitalisations)

z <- covid_hospitalisations %>% 
  left_join(dplyr::select(cdw_full, EAVE_LINKNO, date_ecoss_specimen, flag_covid_symptomatic, test_result_record_source),
            by=c("EAVE_LINKNO","date_ecoss_specimen"))
z <- z %>% mutate(lab = if_else(test_result_record_source == "ECOSS", "nhs","lh") ) %>% 
  dplyr::select(EAVE_LINKNO, date_ecoss_specimen, admission_date, discharge_date, emergency, flag_covid_symptomatic, lab)
covid_hospitalisations <- z

colnames(covid_hospitalisations)[2]<-"SpecimenDate"
covid_hospitalisations <- select(covid_hospitalisations,EAVE_LINKNO,
                                 SpecimenDate, admission_date)

z <- covid_hospitalisations %>% group_by(admission_date) %>% dplyr::summarise(N=n())
z %>% ggplot(aes(x=admission_date, y=N)) + geom_point()+labs(title="covid hospitalisations")
summary(covid_hospitalisations)
###

#use the EAVE severe cases
z <- EAVE_cohort %>% dplyr::select(EAVE_LINKNO, SpecimenDate, icu_death, date_icu_death, NRS.Date.Death) %>% 
  filter(icu_death==1) %>% 
  filter(date_icu_death > a_begin) %>% 
  dplyr::rename(admission_date = date_icu_death) %>% 
  mutate(admission_date = if_else(is.na(NRS.Date.Death) | !is.na(NRS.Date.Death)&(admission_date <= NRS.Date.Death), admission_date, NRS.Date.Death)) %>% 
  dplyr::select(-icu_death, -NRS.Date.Death)
covid_icu_death <- z
z <- covid_icu_death %>% group_by(admission_date) %>% dplyr::summarise(N=n())
z %>% ggplot(aes(x=admission_date, y=N)) + geom_point()+labs(title="covid icu deaths")
summary(covid_icu_death)


#use the EAVE death cases
z <- EAVE_cohort %>% dplyr::select(EAVE_LINKNO, SpecimenDate, death_covid, NRS.Date.Death) %>% 
  filter(death_covid==1) %>% 
  filter(NRS.Date.Death > a_begin) %>% 
  dplyr::rename(admission_date = NRS.Date.Death) %>% 
  dplyr::select(-death_covid)
covid_death <- z
z <- covid_death %>% group_by(admission_date) %>% dplyr::summarise(N=n())
z %>% ggplot(aes(x=admission_date, y=N)) + geom_point() +labs(title="covid deaths")
summary(covid_death)

all_deaths  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/all_deaths.rds"))
summary(all_deaths)
#get covid death certificate deaths
z <- all_deaths %>%  
  mutate(across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9, ~if_else(. %in% c("U071","U072"), 1,0)))
z <- z %>% rowwise() %>% mutate(rowsum = sum(c_across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9))) %>% 
  mutate(covid_death_cert = if_else(rowsum>=1,1,0)) %>% 
  dplyr::select(EAVE_LINKNO, NRS.Date.Death, covid_death_cert)
all_deaths_covid_dth_cert <- z
z <- all_deaths_covid_dth_cert %>% filter(covid_death_cert==1) %>% 
  group_by(NRS.Date.Death) %>% dplyr::summarise(N=n()) %>% filter(NRS.Date.Death >= a_begin)
z %>% ggplot(aes(x=NRS.Date.Death, y=N)) + geom_point() +labs(title="covid deaths - death certificate")
all_deaths_covid_dth_cert <- all_deaths_covid_dth_cert %>% mutate(NRS.Date.Death = as.Date(NRS.Date.Death))
summary(all_deaths_covid_dth_cert)

#update covid_death for deaths in all_deaths
#covid deaths in all_deaths_covid_dth_cert not in covid deaths
z <- all_deaths_covid_dth_cert %>% filter(covid_death_cert==1) %>% 
  filter(NRS.Date.Death > a_begin) %>% 
  filter(!(EAVE_LINKNO %in% covid_death$EAVE_LINKNO)) %>% 
  dplyr::select(-covid_death_cert)
#most of these are in the recent period now find the specimen date of a positive test if there is one
#the may be multiple tests - get the first one within 90 days of the date of death
z <- z %>% left_join(Positive_Tests, by="EAVE_LINKNO") %>% 
  mutate(date_ecoss_specimen = if_else(!is.na(date_ecoss_specimen) & (date_ecoss_specimen > NRS.Date.Death) & (date_ecoss_specimen <= NRS.Date.Death + 6), NRS.Date.Death, date_ecoss_specimen)) %>% 
  mutate(date_ecoss_specimen = if_else(!is.na(date_ecoss_specimen) & (date_ecoss_specimen > NRS.Date.Death+6) ,as.Date(NA_character_),  date_ecoss_specimen)) %>% 
  mutate(diff = as.numeric(NRS.Date.Death - date_ecoss_specimen)) %>% 
  mutate(date_ecoss_specimen = if_else(!is.na(diff) & diff > 90, as.Date(NA_character_), date_ecoss_specimen)) %>% 
  arrange(EAVE_LINKNO, date_ecoss_specimen)
#z <- z %>%  filter(!duplicated(EAVE_LINKNO))
z <- z %>% group_by(EAVE_LINKNO) %>% dplyr::summarise(NRS.Date.Death = first(NRS.Date.Death),
                                                       date_ecoss_specimen = first(date_ecoss_specimen)) %>% ungroup()
z <- z %>% dplyr::rename(SpecimenDate = date_ecoss_specimen, admission_date = NRS.Date.Death)
z <- z %>% dplyr::select(EAVE_LINKNO, SpecimenDate, admission_date)
z <- bind_rows(covid_death, z)
z<-filter(z,admission_date<="2022-03-28")
covid_death <- z

z <- covid_death %>% group_by(admission_date) %>% dplyr::summarise(N=n())
z %>% ggplot(aes(x=admission_date, y=N)) + geom_point() +labs(title="covid deaths")
summary(covid_death)

#combined hospital death endpoint
z1 <- covid_hospitalisations %>% mutate(hosp=1)
z2 <- covid_death %>% mutate(hosp=2)
z <- bind_rows(z1,z2) %>% 
  arrange(EAVE_LINKNO, hosp) %>%  #take hositalisation first
  filter(!duplicated(EAVE_LINKNO))
covid_hosp_death <- z %>% dplyr::select(-hosp)

#cohort + risk groups
df_cohort <- EAVE_cohort %>% dplyr::select(EAVE_LINKNO:ur6_2016_name, age_gp, eave_weight) %>% 
  left_join(rg, by="EAVE_LINKNO")

#below not done as it is not clear what to measure
#z <- Positive_Tests %>%  mutate(days = as.numeric(specimen_date - a_begin)) %>% 
#  arrange(EAVE_LINKNO, desc)
##  mutate(test_before_dec8 = cut(days, breaks = c((min(days)-1), -28, -21, -14, -7, 0, max(days)),
#                       labels=c("1+m", "4w","3w","2w","0-6d","post-vacc")))
#df_cohort <- df_cohort %>% left_join(dplyr::select(z, EAVE_LINKNO, test_before_dec8), by="EAVE_LINKNO")
#df_cohort <- df_cohort %>% mutate(test_before_dec8 = as.character(test_before_dec8)) %>% 
#  mutate(test_before_dec8 = if_else(is.na(test_before_dec8), "no pos test",test_before_dec8) )

z <- readRDS(paste0(project_path_vaccine,"/output/temp/Qcovid.rds"))
z <- z %>% dplyr::select(-(Sex:age_gp), -Q_BMI)
z <- filter(z, !duplicated(EAVE_LINKNO))
z1 <- df_cohort %>% dplyr::select(-(EAVE_ASTHMA:EAVE_HYPERTENSION), -EAVE_Smoke, -EAVE_BP, -n_risk_gps, -eave_weight) %>% 
  left_join(z, by="EAVE_LINKNO")

df_cohort <- filter(z1, !is.na(eave_weight)) #omit any who - need to fix -  do not match

#z <- df_cohort %>% left_join(dplyr::select(Vaccinations, EAVE_LINKNO, date_vacc_1)) %>% 
#mutate(vacc=if_else(is.na(date_vacc_1),0,1))
#z_sub <- slice_sample(z, n=round(nrow(z)*0.25))
#library(mgcv)
#z1 <- glm(vacc ~ age_gp*Sex + simd2020_sc_quintile + n_tests_gp + n_risk_gps,
#                data=z_sub, family=binomial)
#z_ps <- predict(z1, newdata=z, type="response")
# z <- z %>% mutate(prop_score = z_ps) %>% 
#   mutate(inv_psw = if_else(vacc==1, 1/prop_score, (1-vacc)/(1-prop_score))) %>% 
#   mutate(inv_psw = if_else(inv_psw > quantile(inv_psw, 0.95), quantile(inv_psw, 0.95), inv_psw)) %>% 
#   mutate(inv_psw = inv_psw/mean(inv_psw))

#df_cohort <- z %>% dplyr::select(-date_vacc_1, -vacc)

#update weights
#all known to exist - give a weight of 1 and downweight the rest
z_ids <- c(Vaccinations$EAVE_LINKNO, all_deaths$EAVE_LINKNO, covid_hospitalisations$EAVE_LINKNO, 
           cdw_full$EAVE_LINKNO, all_hospitalisations$EAVE_LINKNO) %>% unique()
#summary(filter(EAVE_cohort, !(EAVE_LINKNO %in% z_ids))$eave_weight)
z_N <- round(sum(df_cohort$eave_weight) )
z_k <- sum(df_cohort$EAVE_LINKNO %in% z_ids)
z_m <- round(sum(filter(df_cohort, (EAVE_LINKNO %in% z_ids))$eave_weight))
z <- df_cohort %>% mutate(ew = if_else(EAVE_LINKNO %in% z_ids, 1, eave_weight*(z_N - z_k)/(z_N - z_m)) )
df_cohort <- z %>% dplyr::select(-eave_weight) %>% dplyr::rename(eave_weight=ew)

z <- read_csv(paste0(Location,"/EAVE/GPanalysis/data/restored/map_files/Datazone2011Lookup.csv")) %>% 
  dplyr::select(DataZone, InterZone, Council, HB)
df_cohort <- df_cohort %>% left_join(z, by="DataZone") %>% 
  mutate(HB = if_else(is.na(HB),"Unknown", HB),
         InterZone = if_else(is.na(InterZone),"Unknown", InterZone),
         Council = if_else(is.na(Council),"Unknown", Council))

#saveRDS(df_cohort, paste0(project_path,"/output/temp/df_cohort.rds"))

covid_hosp_death_df <- covid_hosp_death %>%
  mutate(covid_hosp_death_status=1) %>%
  dplyr::rename(covid_hosp_death_date = admission_date)
z <- covid_death %>%
  mutate(covid_death_status=1) %>%
  dplyr::rename(covid_death_date = admission_date)
covid_hosp_death_df <- left_join(covid_hosp_death_df,z)
z <- covid_hospitalisations %>%
  mutate(covid_hosp_status=1) %>%
  dplyr::rename(covid_hosp_date = admission_date)
covid_hosp_death_df <- left_join(covid_hosp_death_df,z)

z<-left_join(df_cohort,Positive_Tests)
z<-left_join(z,covid_hosp_death_df)
z$SpecimenDate<-coalesce(z$SpecimenDate,z$date_ecoss_specimen)
z$date_ecoss_specimen<-NULL
z<-z %>%
  mutate(SpecimenDate_pos = SpecimenDate)
z$SpecimenDate <- NULL
df_cohort <- z
## now colace the two speciment dates for covid-19
saveRDS(df_cohort, paste0(project_path,"/data/df_cohort_15-03-2022.rds"))
saveRDS(covid_hospitalisations, paste0(project_path,"/data/covid_hospitalisations_15-03-2022.rds"))
saveRDS(covid_death, paste0(project_path,"/data/covid_death_15-03-2022.rds"))
saveRDS(covid_hosp_death, paste0(project_path,"/data/covid_hosp_death_15-03-2022.rds"))
saveRDS(Vaccinations, paste0(project_path,"/data/Vaccinations_15-03-2022.rds"))
#df_cohort1 <- readRDS(paste0(project_path,"/output/temp/df_cohort.rds"))

 #adjustment variables
 variables_hosp <- c("age_gp" , "Sex", "simd2020_sc_quintile", "n_tests_gp", "EAVE_BP","Q_DIAG_DIABETES_2"   , "Q_DIAG_COPD"  ,       
                     "Q_DIAG_CKD_LEVEL","Q_DIAG_DEMENTIA","Q_DIAG_STROKE","Q_LEARN_CAT"   ,      
                     "Q_DIAG_FRACTURE","Q_DIAG_NEURO","Q_DIAG_CCF","Q_DIAG_ASTHMA"    ,   
                     "Q_DIAG_EPILEPSY","Q_DIAG_BLOOD_CANCER","Q_DIAG_VTE","Q_DIAG_CIRRHOSIS" ,   
                     "Q_DIAG_RA_SLE","Q_DIAG_PVD","Q_DIAG_AF","Q_DIAG_PULM_RARE" ,   
                     "Q_DIAG_PARKINSONS","Q_DIAG_DIABETES_1","Q_DIAG_PULM_HYPER")

rm(z_sub, z1)
