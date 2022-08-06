##########################################################
# Name of file: 01b_Get_2Dose_Cohort.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@phs.scot
# Original date: 30 Aug 2021
# Latest update author (if not using version control) - Chris Robertson chris.robertson@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: reads in the cohort and merges in vaccination data 
# Approximate run time: Unknown
#
#   run 01a_Vaccinations_Input.R to read in data
#
##########################################################

library(plyr)
library(tidyverse)
library(survival)
library(lubridate)
library(zoo)

#note df_all will be different for the endpoints as the end dates are different for deaths and hospitalisations
Location <- "/conf/"
project_path <- paste0(Location,"EAVE/GPanalysis/progs/UA/second_booster_dose_failures")
df_cohort <- readRDS(paste0(project_path,"/data/df_cohort_13-12-2021.rds"))
covid_hosp_death <- readRDS(paste0(project_path,"/data/covid_hosp_death_13-12-2021.rds"))
covid_death <- readRDS(paste0(project_path,"/data/covid_death_13-12-2021.rds"))
covid_hospitalisations <- readRDS(paste0(project_path,"/data/covid_hospitalisations_13-12-2021.rds"))
Vaccinations <- readRDS(paste0(project_path,"/data/Vaccinations_13-12-2021.rds"))

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
########

z_df <- Vaccinations %>% filter(!is.na(date_vacc_2)) %>% 
  filter(flag_incon==0) %>% 
  dplyr::select(-vacc_type_2, -flag_incon)  # all inconsistencies omitted so vacc_type_2 is the same as vacc_type

z_event <- covid_hosp_death
z_event <- z_event %>% dplyr::rename(event_date=admission_date) #%>% dplyr::select(-SpecimenDate)
print(nrow(z_event))
#print(table(z_event$EAVE_LINKNO %in% EAVE_cohort$EAVE_LINKNO))
a_end <- max(z_event$event_date)

z <- z_df %>% left_join(z_event, by="EAVE_LINKNO")
z <- z %>% filter(is.na(event_date) | event_date > date_vacc_2) #omit those with event before vaccination
z <- z %>% mutate(event = if_else(is.na(event_date), 0L,1L)) %>% 
  mutate(event_date = if_else(is.na(event_date), a_end ,event_date))
z <- z %>% filter(date_vacc_2 < event_date) #drop any vaccinated after end of study
print(sum(z$event))

a_begin <- min(z$date_vacc_2)


df_all <- z #keep df_all - this has all subjects
#add in the covariates to df_all for descriptives of the cohort
df_all <- df_all %>% left_join(df_cohort, by="EAVE_LINKNO")
df_all <- df_all %>%  filter(!is.na(ageYear))  # drop those who do not link into EAVE
df_all <- df_all %>%  mutate(days_between_vacc = as.numeric(date_vacc_2-date_vacc_1)) %>% 
  mutate(vacc_gap = cut(days_between_vacc, breaks=c(1,48, 62, 76, 90, max(days_between_vacc)),
                        labels=c("<7wk","7-8 wk","9-10 wk","11-12 wk", "13+ wk")))
df_all <- df_all %>%  mutate(days_between_vacc = as.numeric(date_vacc_3-date_vacc_2))
df_all$days_between_vacc[is.na(df_all$days_between_vacc)]<--999
df_all<- mutate(df_all,vacc_gap_b = cut(days_between_vacc, breaks=c(1,120, 150, 180, 210, max(days_between_vacc)),
                          labels=c("<4 months","4-5 months","5-6 months","6-7 months", "8+ months")))
#get positive test before vaccination variable
z <- df_all %>% 
  dplyr::select(EAVE_LINKNO, date_vacc_2)
z <- z %>% 
  left_join(Positive_Tests, by="EAVE_LINKNO") %>% 
  filter(!is.na(date_ecoss_specimen)) %>% 
  filter(date_ecoss_specimen < date_vacc_2) #keep thise sample before vaccination
z <- z %>% 
  group_by(EAVE_LINKNO) %>% 
  dplyr::summarise(N=n()) %>% 
  mutate(pos_before_vacc = 1) %>% dplyr::select(-N)
#now link back
df_all <- df_all %>% 
  left_join(z, by="EAVE_LINKNO") %>% 
  mutate(pos_before_vacc=if_else(is.na(pos_before_vacc), 0,pos_before_vacc))

#get number of tests before vaccination variable
z <- df_all %>% 
  dplyr::select(EAVE_LINKNO, date_vacc_2)
z <- z %>% 
  left_join(dplyr::select(cdw_full, EAVE_LINKNO, date_ecoss_specimen), by="EAVE_LINKNO") %>%
  filter(!is.na(date_ecoss_specimen)) %>% #omit those never tested
  filter(date_ecoss_specimen < date_vacc_2) #keep those sample before vaccination
z <- z %>% group_by(EAVE_LINKNO) %>% dplyr::summarise(N=n()) 
#now link back
df_all <- df_all %>% 
  left_join(z, by="EAVE_LINKNO") %>% 
  mutate(N=if_else(is.na(N), 0L,N))
df_all <- df_all %>% 
  mutate(n_tests_gp = cut(N, breaks=c(-1,0,1,2,4,9,max(N)), 
                          labels=c("0","1","2","3-4","5-9","10+")))
df_all$N <- NULL 

#bmi
df_all <- df_all %>% 
  mutate(bmi_gp = cut(bmi_impute, breaks = c(-1, 19, 24, 29, 34,39, 51),
                      labels=c("<20", "20-24","25-29","30-34","34-39","40+")))

df_all <- df_all %>% 
  mutate(Q_DIAG_DIABETES = if_else(Q_DIAG_DIABETES_1==1,Q_DIAG_DIABETES_1,Q_DIAG_DIABETES_2 ),
         Q_DIAG_CKD = if_else(Q_DIAG_CKD_LEVEL >=1, 1, Q_DIAG_CKD_LEVEL))


print(sum(df_all$event))

df_all<-mutate(df_all,age_gp = 
                 cut(ageYear, breaks = c(-Inf, 29,34,39,44,49,54,59,64,69,74,79,Inf),
                     labels=c("18-29","30-34","35-39","40-44","45-49","50-54",
                              "55-59","60-64","65-69","70-74","75-79","80+")))

df_all<-mutate(df_all,ur_combined=if_else(ur6_2016_name=="1 Large Urban Areas"|
                                            ur6_2016_name=="2 Other Urban Areas"|
                                            ur6_2016_name=="3 Accessible Small Towns",1,2))

df <- filter (df_all,vacc_type != "Mo")
df <- filter (df,vacc_type_3 == "Mo"|vacc_type_3 == "PB"|
                is.na(vacc_type_3))
df <- mutate(df, days_dose_23 = date_vacc_3-date_vacc_2)
df <- filter (df, days_dose_23>28 | is.na(days_dose_23))
df$days_dose_23<-NULL

df <- mutate(df, event_full=event)
df <- mutate(df, event=if_else(event_full==1 & covid_hosp_status==1, 1, 0))
df <- filter(df,!is.na(event))
# Elliotts code comes here

# Calculate days between positive test and second dose vaccination, categorising any with positives AFTER dose 2 as "No prior infection":
df <- df %>% mutate(days_covpos_dose2 = as.Date(date_vacc_2)+13 - as.Date(SpecimenDate_pos))
df$days_covpos_dose2 <- as.numeric(df$days_covpos_dose2)
df <- df %>% mutate(months_covpos = (lubridate::interval(SpecimenDate_pos, (date_vacc_2+13))) %/% months(1))
df <- df %>% mutate(months_covpos_grp = if_else(days_covpos_dose2 < 0,10000,days_covpos_dose2))

# Create variable for number of months between prior infection and dose 2, then convert to factor:
df <- df %>% mutate(months_covpos_grp = case_when(days_covpos_dose2 > 0 & months_covpos %in% 0:2.9 ~ "0-3 month",
                                                  months_covpos %in% 3:5.9 ~ "3-6 month",
                                                  months_covpos %in% 6:8.9 ~ "6-9 month",
                                                  months_covpos >= 9 ~ ">9 month",
                                                  is.na(months_covpos_grp) | months_covpos_grp == 10000 | months_covpos_grp == NA ~ "No prior infection",
                                                  days_covpos_dose2 == 0 ~ "0-3 month"))

# Create variable for number of months between prior infection and dose 2, then convert to factor:
# df <- df %>% mutate(months_covpos_grp = case_when(days_covpos_dose2 > 0 & months_covpos %in% 0:1.9 ~ "0-2 month",
#                                                   months_covpos %in% 2:3.9 ~ "2-4 month",
#                                                   months_covpos %in% 4:5.9 ~ "4-6 month",
#                                                   months_covpos %in% 6:7.9 ~ "6-8 month",
#                                                   months_covpos %in% 8:9.9 ~ "8-10 month",
#                                                   months_covpos >= 10 ~ ">10 month",
#                                                   is.na(months_covpos_grp) | months_covpos_grp == 10000 | months_covpos_grp == NA ~ "No prior infection",
#                                                   days_covpos_dose2 == 0 ~ "0-2 month"))

df$months_covpos_grp <- factor(df$months_covpos_grp, ordered = TRUE,
                               levels = c("No prior infection", "0-3 month", "3-6 month", "6-9 month", ">9 month"))

# Checks - ensure correct groups and that there is no "N/A" group:
addmargins(table(df$months_covpos_grp,df$event))
table(is.na(df$months_covpos_grp))

# Create variable for time of prior infection relative to vaccine doses: 
df <- df %>% mutate(prior_pos_timing = case_when(SpecimenDate_pos <= (date_vacc_1+13) ~ "Pre d1",
                                                 SpecimenDate_pos > (date_vacc_1+13) & SpecimenDate_pos < (date_vacc_2+13) ~ "Between d1 & d2",
                                                 SpecimenDate_pos == (date_vacc_2+13) ~ "Between d1 & d2",
                                                 months_covpos_grp == "No prior infection" ~ "No prior infection"))

# Change names of new variables for clarity:
names(df)
names(df)[names(df) == 'days_covpos_dose2'] <- 'prior_infect_days'
names(df)[names(df) == 'months_covpos'] <- 'prior_infect_months'
names(df)[names(df) == 'months_covpos_grp'] <- 'prior_infect_monthgrp'
str(df)
names(df)

df_all<-df
###########################
seed <- 21021954
controls_per_event <- 100
z_event <- df_all %>% filter(!is.na(covid_hosp_death_date)) %>% 
  dplyr::select(EAVE_LINKNO, SpecimenDate, covid_hosp_death_date) %>% 
  dplyr::rename(event_date = covid_hosp_death_date) 
event_ids <- z_event %>% pull(EAVE_LINKNO) %>% unique()
set.seed(seed)
no_event_ids <- filter(df_all, !(EAVE_LINKNO %in% event_ids)) %>% 
  slice_sample(n=controls_per_event*length(event_ids)) %>%
  pull(EAVE_LINKNO) %>% unique()

z_df <- data.frame(EAVE_LINKNO=c(event_ids, no_event_ids)) %>% 
  left_join(dplyr::select(df_all, EAVE_LINKNO, date_vacc_2), by="EAVE_LINKNO") %>% 
  left_join(z_event, by="EAVE_LINKNO") %>% 
  mutate(event = if_else(is.na(event_date), 0L,1L)) %>% 
  mutate(event_date = if_else(is.na(event_date), a_end,event_date))
z_df<-distinct(z_df,EAVE_LINKNO,date_vacc_2,SpecimenDate,event_date,event)
print(sum(z_df$event))

df<- left_join(z_df,df_all,by=
                c("EAVE_LINKNO","date_vacc_2","SpecimenDate","event_date","event"))
df <- distinct(df,EAVE_LINKNO,date_vacc_2,SpecimenDate,event_date,event,.keep_all = TRUE)

df <- df %>% mutate(weight = 
                     if_else(event==1,1,nrow(df_all)/
                               (controls_per_event*length(event_ids)) ) )  


######
df_cohort_vacc <- filter(df, !is.na(date_vacc_2)) %>%
  mutate(analysis_end_date = covid_hosp_death_date)
a_end<-max(df_cohort_vacc$analysis_end_date,na.rm = TRUE)
df_cohort_vacc<-mutate(df_cohort_vacc,analysis_end_date=if_else(analysis_end_date>a_end,a_end,analysis_end_date))
#df_cohort_vacc$analysis_end_date[is.na(df_cohort_vacc$analysis_end_date)]<-max(df_cohort_vacc$date_vacc_1)
df_cohort_vacc$analysis_end_date[is.na(df_cohort_vacc$analysis_end_date)]<-a_end


#change analysis end date for coxph model
#a_end<-max(df_cohort_vacc$date_vacc_1)
df_cohort_vacc<-mutate(df_cohort_vacc,analysis_end_date_c=if_else(analysis_end_date<date_vacc_2,a_end, analysis_end_date))
df_cohort_vacc<-filter(df_cohort_vacc,date_vacc_1<=a_end)

# variable to flag event of interest for all combined
#df_cohort_vacc<-mutate(df_cohort_vacc,event_all1=if_else(!is.na(covid_hosp_death_date),1L, 0L))
#df_cohort_vacc<-mutate(df_cohort_vacc,event_all1=if_else(analysis_end_date<=date_vacc_2,0L, event_all1))
df_cohort_vacc<-mutate(df_cohort_vacc,event_all1=event)

#merging additional variables
# z<-filter(z_chrt_desc,vacc1==1)
# z<-select(z,EAVE_LINKNO_uv,event,event_hosp,event_death,in_hosp_status,prev_positive_status,
#           care_home_elderly,admission_date)%>%
#   mutate(presence=1)
# df_cohort_vacc<-left_join(df_cohort_vacc,z,
#                           by=c("EAVE_LINKNO"="EAVE_LINKNO_uv"))
# df_cohort_vacc<-filter(df_cohort_vacc,!is.na(presence))
length(which(df_cohort_vacc$event==1))
colnames(df_cohort_vacc)[colnames(df_cohort_vacc)=="event"]<-"event_old"

df_cohort_vacc <- mutate(df_cohort_vacc, date_vacc_3_13 = date_vacc_3 +13)
## variables for time varying cox model
df_cohort_vacc<-mutate(df_cohort_vacc,
                       date_preiod_0014 = if_else(analysis_end_date_c>(date_vacc_2+13),(date_vacc_2+13),analysis_end_date_c))
# df_cohort_vacc<-mutate(df_cohort_vacc,
#                        date_preiod_1428 = if_else(analysis_end_date_c>(date_vacc_2+27),(date_vacc_2+27),analysis_end_date_c))

df_cohort_vacc<-mutate(df_cohort_vacc,
                       date_preiod_1449 = 
                         if_else(analysis_end_date_c>(date_vacc_2+49) &
                                   (date_vacc_3_13>(date_vacc_2+49)|is.na(date_vacc_3_13)),
                                 (date_vacc_2+49),
                                 if_else(analysis_end_date_c<(date_vacc_2+49) &
                                           (date_vacc_3_13>(date_vacc_2+49)|is.na(date_vacc_3_13)),
                                         analysis_end_date_c,
                                         if_else(analysis_end_date_c<(date_vacc_2+49) &
                                                   (analysis_end_date_c<date_vacc_3_13|is.na(date_vacc_3)),
                                                 analysis_end_date_c, 
                                                 if_else(analysis_end_date_c==(date_vacc_2+49), 
                                                         analysis_end_date_c,date_vacc_3_13)))))

df_cohort_vacc<-mutate(df_cohort_vacc,
                       date_preiod_4984 = 
                         if_else(analysis_end_date_c>(date_vacc_2+83) &
                                   (date_vacc_3_13>(date_vacc_2+83)|is.na(date_vacc_3_13)),
                                 (date_vacc_2+83),
                                 if_else(analysis_end_date_c<(date_vacc_2+83) &
                                           (date_vacc_3_13>(date_vacc_2+83)|is.na(date_vacc_3_13)),
                                         analysis_end_date_c,
                                         if_else(analysis_end_date_c<(date_vacc_2+83) &
                                                   (analysis_end_date_c<date_vacc_3_13|is.na(date_vacc_3)),
                                                 analysis_end_date_c, 
                                                 if_else(analysis_end_date_c==(date_vacc_2+83), 
                                                         analysis_end_date_c,date_vacc_3_13)))))

df_cohort_vacc<-mutate(df_cohort_vacc,
                       date_preiod_84168 = 
                         if_else(analysis_end_date_c>(date_vacc_2+167) &
                                   (date_vacc_3_13>(date_vacc_2+167)|is.na(date_vacc_3_13)),
                                 (date_vacc_2+167),
                                 if_else(analysis_end_date_c<(date_vacc_2+167) &
                                           (date_vacc_3_13>(date_vacc_2+167)|is.na(date_vacc_3_13)),
                                         analysis_end_date_c,
                                         if_else(analysis_end_date_c<(date_vacc_2+167) &
                                                   (analysis_end_date_c<date_vacc_3_13|is.na(date_vacc_3)),
                                                 analysis_end_date_c, 
                                                 if_else(analysis_end_date_c==(date_vacc_2+167), 
                                                         analysis_end_date_c,date_vacc_3_13)))))

df_cohort_vacc<-mutate(df_cohort_vacc,
                       date_preiod_168o = 
                         if_else(analysis_end_date_c>(date_vacc_2+167) &
                                   (analysis_end_date_c<date_vacc_3_13|is.na(date_vacc_3_13)),
                                 analysis_end_date_c,
                                 if_else(analysis_end_date_c>(date_vacc_2+167) &
                                           (date_vacc_3_13<analysis_end_date_c|is.na(date_vacc_3_13)),
                                         date_vacc_3_13,
                                         if_else(analysis_end_date_c==date_preiod_84168, 
                                                 analysis_end_date_c,date_vacc_3_13))))

########
# current model has 20 days
# df_cohort_vacc<-mutate(df_cohort_vacc,
#                        date_preiod_b_1428 = if_else(analysis_end_date_c>(date_vacc_3_13+14),(date_vacc_3_13+14),analysis_end_date_c))
# df_cohort_vacc<-mutate(df_cohort_vacc,
#                        date_preiod_b_2842 = if_else(analysis_end_date_c>(date_vacc_3_13+28),(date_vacc_3_13+28),analysis_end_date_c))
# df_cohort_vacc<-mutate(df_cohort_vacc,
#                        date_preiod_b_4256 = if_else(analysis_end_date_c>(date_vacc_3_13+42),(date_vacc_3_13+42),analysis_end_date_c))
# df_cohort_vacc<-mutate(df_cohort_vacc,
#                        date_preiod_b_56o = if_else(is.na(date_vacc_3_13),date_vacc_3_13,analysis_end_date_c))
df_cohort_vacc<-mutate(df_cohort_vacc,
                       date_preiod_b_1428 = if_else(analysis_end_date_c>(date_vacc_3_13+14),(date_vacc_3_13+14),analysis_end_date_c))
df_cohort_vacc<-mutate(df_cohort_vacc,
                       date_preiod_b_28o = if_else(is.na(date_vacc_3_13),date_vacc_3_13,analysis_end_date_c))
########

df_cohort_vacc_g<-gather(df_cohort_vacc,"date_period","end_date",date_preiod_0014,
                         date_preiod_1449,date_preiod_4984,date_preiod_84168,
                         date_preiod_168o,date_preiod_b_1428,date_preiod_b_28o)

df_cohort_vacc_g<-arrange(df_cohort_vacc_g,EAVE_LINKNO,end_date)
df_cohort_vacc_g<-df_cohort_vacc_g %>% 
  group_by(EAVE_LINKNO,end_date) %>%
  filter(!duplicated(EAVE_LINKNO,end_date))
df_cohort_vacc_g<-df_cohort_vacc_g %>% 
  group_by(EAVE_LINKNO) %>%
  mutate(start_date=lag(end_date,1))
df_cohort_vacc_g$start_date<-coalesce(df_cohort_vacc_g$start_date,df_cohort_vacc_g$date_vacc_2)
df_cohort_vacc_g<-mutate(df_cohort_vacc_g,end_date1=if_else(end_date==start_date,end_date+1,end_date))
df_cohort_vacc_g<-df_cohort_vacc_g %>% 
  group_by(EAVE_LINKNO) %>%
  mutate(start_time_date=first(start_date))
df_cohort_vacc_g$start_time<-df_cohort_vacc_g$start_date-df_cohort_vacc_g$start_time_date
df_cohort_vacc_g$end_time<-df_cohort_vacc_g$end_date1-df_cohort_vacc_g$start_time_date
df_cohort_vacc_g$event_old<-as.integer(df_cohort_vacc_g$event_old)
df_cohort_vacc_g <- drop_na(df_cohort_vacc_g, end_time)
df_cohort_vacc_g<-df_cohort_vacc_g %>%
  group_by(EAVE_LINKNO) %>%
  mutate(event=if_else(row_number()==n(),event_old,0L))


#df_cohort_vacc_g1<-df_cohort_vacc_g
# df_cohort_vacc_g <- mutate(df_cohort_vacc_g, date_period = 
#                              if_else(date_period=="date_preiod_b_14o" & vacc_type_3=="Mo", 
#                                      "date_preiod_b_Mo",
#                                      if_else(date_period=="date_preiod_b_14o" & vacc_type_3=="PB",
#                                              "date_preiod_b_PB",date_period)))
df_cohort_vacc_g <- mutate(df_cohort_vacc_g, date_period =
                             if_else(date_period=="date_preiod_b_1428" & vacc_type_3=="Mo",
                                     "date_preiod_b_1428_Mo",
                                     if_else(date_period=="date_preiod_b_28o" & vacc_type_3=="Mo",
                                             "date_preiod_b_28o_Mo",date_period)))

df_cohort_vacc_g <- mutate(df_cohort_vacc_g, date_period =
                             if_else(date_period=="date_preiod_b_1428" & vacc_type_3=="PB",
                                     "date_preiod_b_1428_PB",
                                     if_else(date_period=="date_preiod_b_28o" & vacc_type_3=="PB",
                                             "date_preiod_b_28o_PB",date_period)))

#### adding background level of infection
z_event<-filter(Positive_Tests,date_ecoss_specimen>=a_begin-13)
z_event<-z_event %>%
  group_by(date_ecoss_specimen) %>%
  mutate(num_pos = n())
z_event$EAVE_LINKNO<-NULL
z_event<-filter(z_event,!duplicated(date_ecoss_specimen))
z_event<-arrange(z_event,date_ecoss_specimen)
z_event$num_pos_avg<-rollsumr(z_event$num_pos,k=14,fill=NA)
z_event$num_pos_avg<-z_event$num_pos_avg/14
colnames(z_event)[1]<-"start_date"
z_event$num_pos_avg<-round(z_event$num_pos_avg)
#z <- z_event %>% group_by(date_ecoss_specimen) %>% dplyr::summarise(N=num_pos_avg)
#z %>% ggplot(aes(x=date_ecoss_specimen, y=N)) + geom_point()+labs(title="covid infections")
#z_event<-mutate(z_event,backg_inf=if_else(num_pos_avg<1000,'low',if_else(num_pos_avg>=1000&num_pos_avg<2000,'medium','high')))
z<-data.frame(seq(as.Date("2020-12-30"),as.Date("2021-09-19"),"day"))
colnames(z)[1]<-"start_date"
z<-left_join(z,z_event,by="start_date") %>%
  select(-num_pos)
df_cohort_vacc_g<-left_join(df_cohort_vacc_g,z,by="start_date")

# adding week of event as confounder
z<-data.frame(seq(as.Date("2020-12-08"),as.Date("2022-01-24"),"day"))
colnames(z)[1]<-"start_date"
z<-mutate(z,variant=ifelse(start_date<"2021-05-19","Alpha",
                           ifelse(start_date<"2021-12-15","Delta","Omicron")))
df_cohort_vacc_g<-left_join(df_cohort_vacc_g,z,by="start_date")

#saveRDS(df_cohort_vacc, paste0(project_path,"/data/df_cohort_vacc_13-12-2021_gam.rds"))
#saveRDS(df_cohort_vacc_g, paste0(project_path,"/data/df_cohort_vacc_g_13-12-2021_gam_hosp.rds"))
