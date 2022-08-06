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
library(survminer)
#library(dplyr)
#library(mgcv)
library(tidyr)
library(ggplot2)

Location <- "/conf/"
project_path <- paste0(Location,"EAVE/GPanalysis/progs/UA/second_booster_dose_failures")
project_path_vaccine <- paste0(Location,"EAVE/GPanalysis/progs/CR/Vaccine")
#df_cohort <- readRDS(paste0(project_path,"/data/df_cohort_22-09-2021.rds"))
df_cohort_vacc <- readRDS(paste0(project_path,"/data/df_cohort_vacc_15-03-2022_gam.rds"))
df_cohort_vacc_g <- readRDS(paste0(project_path,"/data/df_cohort_vacc_g_15-03-2022_gam.rds"))
#covid_hosp_death <- readRDS(paste0(project_path,"/data/covid_hosp_death_22-09-2021.rds"))
#covid_death <- readRDS(paste0(project_path,"/data/covid_death_22-09-2021.rds"))
#covid_hospitalisations <- readRDS(paste0(project_path,"/data/covid_hospitalisations_22-09-2021.rds"))
#Vaccinations <- readRDS(paste0(project_path,"/data/Vaccinations_22-09-2021.rds"))
##################
# data prep for the analysis

# collapse urban rural for ARR analysis
df_cohort_vacc_g$ur_combined <- as.factor(df_cohort_vacc_g$ur_combined)
df_cohort_vacc_g$ur_combined <- relevel(df_cohort_vacc_g$ur_combined, ref = 1)

#df_cohort_vacc_g$vacc_type_3[is.na(df_cohort_vacc_g$vacc_type_3)]<-"NB"
#df_cohort_vacc_g$vacc_type_3 <- as.factor(df_cohort_vacc_g$vacc_type_3)
#df_cohort_vacc_g$vacc_type_3 <- relevel(df_cohort_vacc_g$vacc_type_3, ref = "NB")
#df_cohort_vacc_g$booster_status <- as.factor(df_cohort_vacc_g$booster_status)
#df_cohort_vacc_g$booster_status <- relevel(df_cohort_vacc_g$booster_status, ref = "nb")
df_cohort_vacc_g$date_period <- as.factor(df_cohort_vacc_g$date_period)
df_cohort_vacc_g$date_period <- relevel(df_cohort_vacc_g$date_period, ref = "date_preiod_1449")
df_cohort_vacc_g$simd2020_sc_quintile <- as.factor(df_cohort_vacc_g$simd2020_sc_quintile)
df_cohort_vacc_g$simd2020_sc_quintile <- relevel(df_cohort_vacc_g$simd2020_sc_quintile, ref = "5-Low")
df_cohort_vacc_g$prior_pos_timing <- as.factor(df_cohort_vacc_g$prior_pos_timing)
df_cohort_vacc_g$prior_pos_timing <- relevel(df_cohort_vacc_g$prior_pos_timing, ref = "No prior infection")

df_cohort_vacc_g <- mutate(df_cohort_vacc_g, p_years=as.numeric(end_time-start_time))
df_cohort_vacc_g$num_pos_avg <- df_cohort_vacc_g$num_pos_avg/100

df_cohort_vacc_g$prior_infect_monthgrp1[df_cohort_vacc_g$prior_infect_monthgrp=="No prior infection"]<-0
df_cohort_vacc_g$prior_infect_monthgrp1[df_cohort_vacc_g$prior_infect_monthgrp=="0-3 month"]<-1
df_cohort_vacc_g$prior_infect_monthgrp1[df_cohort_vacc_g$prior_infect_monthgrp=="3-6 month"]<-2
df_cohort_vacc_g$prior_infect_monthgrp1[df_cohort_vacc_g$prior_infect_monthgrp=="6-9 month"]<-3
df_cohort_vacc_g$prior_infect_monthgrp1[df_cohort_vacc_g$prior_infect_monthgrp==">9 month"]<-4
df_cohort_vacc_g$prior_infect_monthgrp1<-as.factor(df_cohort_vacc_g$prior_infect_monthgrp1)

df_cohort_vacc_g$variant <- as.factor(df_cohort_vacc_g$variant)
df_cohort_vacc_g$variant <- relevel(df_cohort_vacc_g$variant, ref = "Delta")

#df_cohort_vacc_g <- df_cohort_vacc_g %>% 
#  mutate(bmi_gp1 = cut(bmi_impute, breaks = c(-1, 18.5, 24.9, 29.9, 39.9, 51),
#                       labels=c("<18.5", "18.5-24.9","25-29.9","30-39.9","40+")))
df_cohort_vacc_g$bmi_gp <- as.factor(df_cohort_vacc_g$bmi_gp)
df_cohort_vacc_g$bmi_gp <- relevel(df_cohort_vacc_g$bmi_gp, ref = "18.5-24.9")
#########
z <- readRDS(paste0(project_path_vaccine,"/output/temp/Qcovid.rds"))
z <- z %>% dplyr::select(EAVE_LINKNO, Q_BMI)
z <- filter(z, !duplicated(EAVE_LINKNO))
z <- z %>% mutate(bmi_gp1 = cut(Q_BMI, breaks = c(-1, 18.5, 24.9, 29.9, 39.9, 51),
                       labels=c("<18.5", "18.5-24.9","25-29.9","30-39.9","40+")))
z$bmi_gp1<-as.character(z$bmi_gp1)
z$bmi_gp1[is.na(z$bmi_gp1)]<-"NA"
df_cohort_vacc_g <- left_join(df_cohort_vacc_g,z,by="EAVE_LINKNO")
df_cohort_vacc_g$bmi_gp1 <- as.factor(df_cohort_vacc_g$bmi_gp1)
df_cohort_vacc_g$bmi_gp1 <- relevel(df_cohort_vacc_g$bmi_gp1, ref = "18.5-24.9")
######### Modelling
# crude analysis
# sex
# pyears - rate of event per year
# glm_poisson - log rate ratios
################
# df_cohort_vacc_g1 <- df_cohort_vacc_g
# df_cohort_vacc_g <- filter(df_cohort_vacc_g, date_period!="date_preiod_0014")
# 
# glm_pos1 <- glm(event~offset(log(p_years))  +date_period+ Sex + age_gp + vacc_gap+ variant+
#                   prior_infect_monthgrp1+n_risk_gps + n_tests_gp+simd2020_sc_quintile+num_pos_avg+
#                   ur_combined+bmi_gp1 + HB , family=poisson,weights = weight,data=df_cohort_vacc_g)
# z<-as.data.frame(round(exp(glm_pos1$coefficients),3))
# z<-cbind(z,round(exp(glm_pos1$coefficients-coef(summary(glm_pos1))[,2]*1.96),3))
# z<-cbind(z,round(exp(glm_pos1$coefficients+coef(summary(glm_pos1))[,2]*1.96),3))
# z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))
# colnames(z)[1]<-"Adjusted Rate Ratio"
# colnames(z)[2]<-"LCI"
# colnames(z)[3]<-"UCI"
# colnames(z)[4]<-"ARR (LCI-UCI)"
# View(z)
# 
# write.csv(z,paste0(project_path,"/output/ARR_14-28_days_as_ref_b_qbmi.csv"))
# write.csv(z,paste0(project_path,"/output/ARR_AZ_14-28_days_as_ref_b_qbmi.csv"))
# write.csv(z,paste0(project_path,"/output/ARR_PB_14-28_days_as_ref_b_qbmi.csv"))
# 
# df_cohort_vacc_g_both<-df_cohort_vacc_g
# df_cohort_vacc_g <- filter(df_cohort_vacc_g_both,vacc_type=="AZ")
# df_cohort_vacc_g<-df_cohort_vacc_g_both
# df_cohort_vacc_g <- filter(df_cohort_vacc_g_both,vacc_type=="PB")
# ##########################
# 
# conditions <- c('Q_DIAG_ASTHMA','Q_DIAG_CKD','Q_DIAG_CIRRHOSIS','Q_DIAG_NEURO','Q_DIAG_CCF',
#                 'Q_DIAG_DIABETES_1','Q_DIAG_DIABETES_2','Q_DIAG_DEMENTIA','Q_DIAG_CHD')
# ARR <- rbind()
# for (i in 1:1){#length(conditions)
#   temp <- mutate(df_cohort_vacc_g, risk_gp = n_risk_gps)
#   temp$risk_gp <- as.numeric(temp$risk_gp)
#   temp <- mutate(temp, risk_gp=risk_gp-1)
#   temp <- mutate(temp, risk_gp=risk_gp-get(conditions[i]))
#   temp$risk_gp<-as.factor(temp$risk_gp)
#   
#   glm_pos <- glm(event~offset(log(p_years))  +get(conditions[i])+ Sex + age_gp +
#                    simd2020_sc_quintile + risk_gp, 
#                  family=poisson,weights = weight,data=temp)
#   #summary(glm_pos)
#   z<-as.data.frame(round(exp(glm_pos$coefficients),3))
#   z<-cbind(z,round(exp(glm_pos$coefficients-coef(summary(glm_pos))[,2]*1.96),2))
#   z<-cbind(z,round(exp(glm_pos$coefficients+coef(summary(glm_pos))[,2]*1.96),2))
#   colnames(z)[1]<-"Adjusted Rate Ratio"
#   colnames(z)[2]<-"LCI"
#   colnames(z)[3]<-"UCI"
#   rownames(z)[2]<-conditions[i]
#   ARR <- rbind(ARR,z[2,])
#   #print(z)
# }
# ARR<-mutate(ARR, ARR_LCI_UCI = paste(ARR[,1]," (",ARR[,2],"-",ARR[,3],")"))
# colnames(ARR)[4]<-"ARR (LCI-UCI)"
# write.csv(ARR,paste0(project_path,"/output/ARR_conds_14-28_days_as_ref.csv"))
# write.csv(ARR,paste0(project_path,"/output/ARR_conds_14-28_days_as_ref_AZ.csv"))
# write.csv(ARR,paste0(project_path,"/output/ARR_conds_14-28_days_as_ref_PB.csv"))
# 
# 
# p_years1<-pyears(Surv(start_time,end_time,event)~date_period + prior_infect_monthgrp1,
#                  weights = weight , data = df_cohort_vacc_g1, data.frame = TRUE)
# summary(p_years1)
# p_years<-pyears(Surv(start_time,end_time,event)~date_period + prior_infect_monthgrp1,
#                 weights = weight , data = df_cohort_vacc_g, data.frame = TRUE)
# summary(p_years)
# temp<-as.data.frame(p_years1$data)
# write.csv(temp,paste0(project_path,"/output/ARR_pyears_14-28_days_as_ref.csv"))
# 
# # sensitivity check on hospital admission with covid-19
# temp<-filter(df_cohort_vacc_g1, date_period=="date_preiod_0014")
# temp<-filter(temp, event==1) 
# temp<-select(temp, EAVE_LINKNO, event, date_vacc_2, event_date)
# any_hospitalisation_post_01022020 <- readRDS("/conf/EAVE/GPanalysis/data/any_hospitalisation_post_01022020.rds")
# temp<-left_join(temp,any_hospitalisation_post_01022020)
# temp<-mutate(temp, n_days = admission_date- date_vacc_2)
# temp<-filter(temp, n_days <15 & n_days>-15)

# interaction

glm_pos <- glm(event~offset(log(p_years))  +date_period+ bmi_gp+ Sex + age_gp +
                 simd2020_sc_quintile + n_risk_gps + bmi_gp*date_period, 
               family=poisson,weights = weight,data=df_cohort_vacc_g)

z<-as.data.frame(round(exp(glm_pos$coefficients),3))
z<-cbind(z,round(exp(glm_pos$coefficients-coef(summary(glm_pos))[,2]*1.96),2))
z<-cbind(z,round(exp(glm_pos$coefficients+coef(summary(glm_pos))[,2]*1.96),2))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))
colnames(z)[1]<-"Adjusted Rate Ratio"
colnames(z)[2]<-"LCI"
colnames(z)[3]<-"UCI"
colnames(z)[4]<-"ARR (LCI-UCI)"

write.csv(z,paste0(project_path,"/output/interaction_bmi.csv"))
###############################