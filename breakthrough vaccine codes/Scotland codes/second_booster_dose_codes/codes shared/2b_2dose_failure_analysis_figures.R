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
covid_hosp_death <- readRDS(paste0(project_path,"/data/covid_hosp_death_22-09-2021.rds"))
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

###############################
z_df <- df_cohort_vacc1 %>%
  dplyr::select(age_gp,Sex, simd2020_sc_quintile, ur6_2016_name, n_risk_gps,
                n_tests_gp,vacc_gap, bmi_gp,prior_infect_monthgrp,event) %>% 
  pivot_longer(cols=age_gp:prior_infect_monthgrp) 

z_df <- z_df %>% 
  group_by(name, value) %>% 
  dplyr::summarise( both_vacc_event = sum(event==1), 
                    total_vacc_adm = sum(!is.na(event))) %>%  
  ungroup()

z_df <- z_df %>% group_by(name) %>% 
  dplyr::mutate( total_event = sum(total_vacc_adm)) %>%  
  ungroup() %>%
  mutate(Percent_of_vacc_adm = round(total_vacc_adm/total_event*100,1))

#############
# df_cohort_vacc_g_pyrs$ur6_2016_name[is.na(df_cohort_vacc_g_pyrs$ur6_2016_name)]<-"unknown"
df_cohort_vacc_g1 <- filter(df_cohort_vacc_g1, date_period!="date_preiod_0014")
chars<-c('age_gp','Sex','simd2020_sc_quintile','ur6_2016_name',
         'n_risk_gps','n_tests_gp','vacc_gap', 'bmi_gp','prior_infect_monthgrp')
event_list<-rbind()
for (i in 1:9){
  p_years<-pyears(Surv(start_time,end_time,event)~get(chars[i]), scale = 365.25 ,
                  data = df_cohort_vacc_g1, data.frame = TRUE,weights = weight)
  z<-data.frame(p_years$data$`get(chars[i])`)
  colnames(z)[1]<-"value"
  z<-mutate(z,name=chars[i])
  z1<-data.frame(round((p_years$data$event)*1000/(p_years$data$pyears),1))
  colnames(z1)[1]<-"rate_per_1000_yrs"
  z<-bind_cols(z,z1)
  event_list<-rbind(event_list,z)
}
z_df <- left_join(z_df,event_list)
z_df <- cbind(z_df, paste(z_df$total_vacc_adm," (",z_df$Percent_of_vacc_adm,")", sep = ""))
z_df <- cbind(z_df, paste(z_df$both_vacc_event," (",z_df$rate_per_1000_yrs,")", sep = ""))

##########
colnames(z_df)[8] <- "total vacc (n, %)"
colnames(z_df)[9] <- "severe outcome (n, rate per 1000 person years)"
results <- select(z_df,name, value,`total vacc (n, %)`,
                  `severe outcome (n, rate per 1000 person years)`)
results_f <- results

colnames(z_df)[8] <- "total vacc AZ (n, %)"
colnames(z_df)[9] <- "severe outcome AZ (n, rate per 1000 person years)"
results <- select(z_df,name, value,`total vacc AZ (n, %)`,
                  `severe outcome AZ (n, rate per 1000 person years)`)
results_f<-left_join(results_f,results)

colnames(z_df)[8] <- "total vacc PB (n, %)"
colnames(z_df)[9] <- "severe outcome PB (n, rate per 1000 person years)"
results <- select(z_df,name, value,`total vacc PB (n, %)`,
                  `severe outcome PB (n, rate per 1000 person years)`)
results_f<-left_join(results_f,results)
write.csv(results_f,paste0(project_path,"/output/cohort_description.csv"))

p_years<-pyears(Surv(start_time,end_time,event)~Sex, scale = 365.25 ,
                data = df_cohort_vacc_g1, data.frame = TRUE, weights = weight)
sum(p_years$data$event)*1000/sum(p_years$data$pyears)

tt<-filter(df_cohort_vacc1,event==1)
table(tt$covid_hosp_status)
table(tt$covid_death_status)
tt<-filter(df_cohort_vacc1,covid_death_status==1)
table(tt$covid_hosp_status)

df_cohort_vacc_g <- df_cohort_vacc_g1

a_end <- as.Date("2021-05-19")
df_cohort_vacc_g1 <- df_cohort_vacc_g1 %>% filter(start_date < a_end)

df_cohort_vacc_g1 <- df_cohort_vacc_g
df_cohort_vacc_g1 <- df_cohort_vacc_g1 %>% filter(start_date >= a_end &
                                                    start_date<"2021-12-15") 

df_cohort_vacc_g1 <- df_cohort_vacc_g
df_cohort_vacc_g1 <- df_cohort_vacc_g1 %>% filter(start_date >= "2021-12-15") 
###############