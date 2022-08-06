##########################################################
# Name of file: 01a_Vaccinations_Input.R
# Original author(s): Utkarsh Agrawal
# Description of content: Analysis of vaccination failure cohort 
# Approximate run time: Unknown
##########################################################

# 01 Setup ####
#Libraries
library(plyr)
library(tidyverse)
library(survival)
#library(dplyr)
#library(mgcv)
library(tidyr)
library(ggplot2)

Location <- "/conf/"
project_path <- paste0(Location,"EAVE/GPanalysis/progs/UA/second_booster_dose_failures")
#df_cohort <- readRDS(paste0(project_path,"/data/df_cohort_22-09-2021.rds"))
df_cohort_vacc1 <- readRDS(paste0(project_path,"/data/df_cohort_vacc_13-12-2021_gam_full.rds"))
df_cohort_vacc_g1 <- readRDS(paste0(project_path,"/data/df_cohort_vacc_g_13-12-2021_gam.rds"))
covid_hosp_death <- readRDS(paste0(project_path,"/data/covid_hosp_death_15-03-2022.rds"))
#covid_death <- readRDS(paste0(project_path,"/data/covid_death_22-09-2021.rds"))
#covid_hospitalisations <- readRDS(paste0(project_path,"/data/covid_hospitalisations_22-09-2021.rds"))
#Vaccinations <- readRDS(paste0(project_path,"/data/Vaccinations_22-09-2021.rds"))
##################
# data prep for the analysis
df_cohort_vacc1 <- filter(df_cohort_vacc1, !duplicated(EAVE_LINKNO))
tt <-select(df_cohort_vacc_g1, EAVE_LINKNO, date_period, event)
tt<-filter(tt,date_period=="date_preiod_0014"&event==1)
tt<-mutate(tt,presence=1)
tt$date_period<-NULL
tt$event<-NULL
df_cohort_vacc1 <- left_join(df_cohort_vacc1,tt, by="EAVE_LINKNO")
df_cohort_vacc1<-filter(df_cohort_vacc1,is.na(presence))
df_cohort_vacc1$presence<-NULL
rm(tt)
#############

df_cohort_vacc_g_both<-df_cohort_vacc_g1
df_cohort_vacc_both <- df_cohort_vacc1
df_cohort_vacc_g1 <- filter(df_cohort_vacc_g_both,vacc_type=="AZ")
df_cohort_vacc1 <- filter(df_cohort_vacc_both,vacc_type=="AZ")
df_cohort_vacc_g1<-df_cohort_vacc_g_both
df_cohort_vacc1<-df_cohort_vacc_both
df_cohort_vacc_g1 <- filter(df_cohort_vacc_g_both,vacc_type=="PB")
df_cohort_vacc1 <- filter(df_cohort_vacc_both,vacc_type=="PB")

##########
all_deaths  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/all_deaths.rds"))
summary(all_deaths)

# underlying cause deaths
tt<-all_deaths
tt<-filter(tt,UNDERLYING_CAUSE_OF_DEATH=="U071"|UNDERLYING_CAUSE_OF_DEATH=="U072"|
             UNDERLYING_CAUSE_OF_DEATH=="B342"|UNDERLYING_CAUSE_OF_DEATH=="B972")
tt<-select(tt, EAVE_LINKNO)
tt<-mutate(tt,presence=1)
z<-filter(df_cohort_vacc1,event==1&covid_death_status==1)
z<-left_join(z,tt)
z$presence[is.na(z$presence)]<-0
table(z$presence)
#death recorded at any position
tt<-all_deaths
tt<-gather(tt,"cause_posn","death_reason",8:18)
tt<-drop_na(tt,death_reason)
tt<-filter(tt,death_reason=="U071"|death_reason=="U072"|
             death_reason=="B342"|death_reason=="B972")
tt<-select(tt, EAVE_LINKNO)
tt<-mutate(tt,presence=1)
tt<-distinct(tt,EAVE_LINKNO,.keep_all = TRUE)
z<-filter(df_cohort_vacc1,event==1&covid_death_status==1)
z<-left_join(z,tt)
z$presence[is.na(z$presence)]<-0
table(z$presence)
###############

smr01_2022_02_08 <- readRDS("/conf/EAVE/GPanalysis/data/smr01_2022_02_08.rds")
tt<-smr01_2022_02_08
tt<-filter(tt,main_diagnosis_on_admission=="U071"|main_diagnosis_on_admission=="U072"|
             main_diagnosis_on_admission=="B342"|main_diagnosis_on_admission=="B972")
tt<-select(tt, EAVE_LINKNO)
tt<-mutate(tt,presence=1)
tt<-distinct(tt,EAVE_LINKNO,.keep_all = TRUE)
z<-filter(df_cohort_vacc1,event==1&covid_hosp_status==1)
z<-left_join(z,tt)
z$presence[is.na(z$presence)]<-0
table(z$presence)
# hosp recorded any position
tt<-smr01_2022_02_08
tt<-gather(tt,"cause_posn","hosp_reason",4:9)
tt<-drop_na(tt,hosp_reason)
tt<-filter(tt,hosp_reason=="U071"|hosp_reason=="U072"|
             hosp_reason=="B342"|hosp_reason=="B972")
tt<-select(tt, EAVE_LINKNO)
tt<-mutate(tt,presence=1)
tt<-distinct(tt,EAVE_LINKNO,.keep_all = TRUE)
z<-filter(df_cohort_vacc1,event==1&covid_hosp_status==1)
z<-left_join(z,tt)
z$presence[is.na(z$presence)]<-0
table(z$presence)
# covid prior to hosp
z<-filter(df_cohort_vacc1,event==1&covid_hosp_status==1)
z<-mutate(z,days_pos_before_covid=covid_hosp_date-SpecimenDate_pos)
z<-filter(z,days_pos_before_covid>=0&days_pos_before_covid<=14)

#### ARR_BMI_both_vacc
ARR_data <- read.csv("/conf/EAVE/GPanalysis/progs/UA/booster_dose_failures/output/ARR_14-28_days_as_ref.csv")

#3rd_vacc_status
ARR_data1 <- ARR_data[c(2:6),]
names<-c("2_booster_dose_5-6weeks","3_booster_dose_7-8weeks",
         "4_booster_dose_9-10weeks","5_booster_dose_>=11weeks","1_booster_dose_3-4weeks")
ARR_data1[5,2]<-1
ARR_data1[5,3]<-0.001
ARR_data1[5,4]<-1.001
ARR_data1[,1]<-names

ggplot(ARR_data1, aes(x = Adjusted.Rate.Ratio, y = X)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = LCI, xmin = UCI), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") +
  scale_x_continuous(breaks = seq(0.5, 3.5, 0.5), labels = seq(0.5, 3.5, 0.5),
                     limits = c(0.5,3.5)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.text = element_text(face="bold")) +
  ylab("") +
  xlab("Adjusted rate ratio") +
  ggtitle("both vaccines combined")

#2nd_3rd_vacc_status
ARR_data1 <- ARR_data[c(10:20),]
names<-c("Age_30-34_yrs","Age_35-39_yrs","Age_40-44_yrs","Age_45-49_yrs",
         "Age_50-54_yrs","Age_55-59_yrs","Age_60-64_yrs","Age_65-69_yrs",
         "Age_70-74_yrs","Age_75-79_yrs","Age_80+_yrs")
ARR_data1[,1]<-names

ggplot(ARR_data1, aes(x = Adjusted.Rate.Ratio, y = X)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = LCI, xmin = UCI), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") +
  scale_x_continuous(breaks = seq(0.75, 5, 0.5), labels = seq(0.75, 5, 0.5),
                     limits = c(0.9,5)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Adjusted rate ratio") +
  ggtitle("both vaccines combined")
###############