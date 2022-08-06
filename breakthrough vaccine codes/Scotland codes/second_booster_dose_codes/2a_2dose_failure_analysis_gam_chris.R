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
xdf_full_covid_hosp_death <- readRDS(paste0(project_path,"/data/xdf_full_covid_hosp_death.RDS"))
#df_cohort_vacc <- readRDS(paste0(project_path,"/data/df_cohort_vacc_15-03-2022_gam.rds"))
#df_cohort_vacc_g <- readRDS(paste0(project_path,"/data/df_cohort_vacc_g_15-03-2022_gam.rds"))
##################
# data prep for the analysis
# method 1

df_cohort_vacc_g <- filter(xdf_full_covid_hosp_death, pv_period_f!="uv")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v1_0:3")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v1_4+")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v2_0:1")

# tt1<-filter(df_cohort_vacc_g,event_date>"2021-12-14"&event==1)
# tt<-filter(df_cohort_vacc_g,event==0)
# tt<-rbind(tt,tt1)
# tt<-arrange(tt,EAVE_LINKNO,period_f)
# tt$pv_period_f<-as.character(tt$pv_period_f)
df_cohort_vacc_g$pv_period_f<-as.character(df_cohort_vacc_g$pv_period_f)
df_cohort_vacc_g <- mutate(df_cohort_vacc_g, pv_period_f =
                             if_else(pv_period_f=="v3_2:4" & vacc_type_3=="Mo",
                                     "v3_2:4_Mo",
                                     if_else(pv_period_f=="v3_5:7" & vacc_type_3=="Mo",
                                             "v3_5:7_Mo",
                                             if_else(pv_period_f=="v3_8+" & vacc_type_3=="Mo",
                                                     "v3_8+_Mo",pv_period_f))))

df_cohort_vacc_g <- mutate(df_cohort_vacc_g, pv_period_f =
                             if_else(pv_period_f=="v3_2:4" & vacc_type_3=="PB",
                                     "v3_2:4_PB",
                                     if_else(pv_period_f=="v3_5:7" & vacc_type_3=="PB",
                                             "v3_5:7_PB",
                                             if_else(pv_period_f=="v3_8+" & vacc_type_3=="PB",
                                                     "v3_8+_PB",pv_period_f))))

df_cohort_vacc_g$pv_period_f <- as.factor(df_cohort_vacc_g$pv_period_f)
df_cohort_vacc_g$pv_period_f <- relevel(df_cohort_vacc_g$pv_period_f, ref = "v2_2:9")

df_cohort_vacc_g$simd2020_sc_quintile <- as.factor(df_cohort_vacc_g$simd2020_sc_quintile)
df_cohort_vacc_g$simd2020_sc_quintile <- relevel(df_cohort_vacc_g$simd2020_sc_quintile, ref = "5-Low")

df_cohort_vacc_g$bmi_gp <- as.factor(df_cohort_vacc_g$bmi_gp)
df_cohort_vacc_g$bmi_gp <- relevel(df_cohort_vacc_g$bmi_gp, ref = "18.5-24.9")

df_cohort_vacc_g$vacc_gap <- as.factor(df_cohort_vacc_g$vacc_gap)
df_cohort_vacc_g$vacc_gap <- relevel(df_cohort_vacc_g$vacc_gap, ref = "<7wk")

df_cohort_vacc_g$ur_combined <- as.factor(df_cohort_vacc_g$ur_combined)
df_cohort_vacc_g$ur_combined <- relevel(df_cohort_vacc_g$ur_combined, ref = "2")

df_cohort_vacc_g$age_gp <- as.factor(df_cohort_vacc_g$age_gp)
df_cohort_vacc_g$age_gp <- relevel(df_cohort_vacc_g$age_gp, ref = "18-49")

df_cohort_vacc_g$n_tests_gp <- as.character(df_cohort_vacc_g$n_tests_gp)
df_cohort_vacc_g$n_tests_gp[df_cohort_vacc_g$n_tests_gp=="3"] <- "3-4"
df_cohort_vacc_g$n_tests_gp[df_cohort_vacc_g$n_tests_gp=="4"] <- "3-4"
df_cohort_vacc_g$n_tests_gp[df_cohort_vacc_g$n_tests_gp=="10-19"] <- "10+"
df_cohort_vacc_g$n_tests_gp[df_cohort_vacc_g$n_tests_gp=="20+"] <- "10+"
df_cohort_vacc_g$n_tests_gp <- as.factor(df_cohort_vacc_g$n_tests_gp)
df_cohort_vacc_g$n_tests_gp <- relevel(df_cohort_vacc_g$n_tests_gp, ref = "0")

z.agg <- pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt) ~ 
                  pv_period_f+period_f + Sex + age_gp+n_risk_gps+
                  bmi_gp+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined, 
                data=df_cohort_vacc_g, weight=ew , scale=365.25, data.frame=TRUE)
glm_pos1 <- glm(event ~ offset(log(pyears)) +pv_period_f+period_f + Sex + age_gp+
                  n_risk_gps+bmi_gp+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined, 
                family=poisson, data=z.agg$data)

z<-as.data.frame(round(exp(glm_pos1$coefficients),3))
z<-cbind(z,round(exp(glm_pos1$coefficients-coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z,round(exp(glm_pos1$coefficients+coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))
colnames(z)[1]<-"Adjusted Rate Ratio"
colnames(z)[2]<-"LCI"
colnames(z)[3]<-"UCI"
colnames(z)[4]<-"ARR (LCI-UCI)"
View(z)

z<-as.data.frame(coef(summary(glm_pos1)))
write.csv(z,paste0(project_path,"/output/ARR_14-28_days_as_ref_b.csv"))
z<-as.data.frame(coef(summary(glm_pos1)))
write.csv(z,paste0(project_path,"/output/ARR_AZ_14-28_days_as_ref.csv"))
z<-as.data.frame(coef(summary(glm_pos1)))
write.csv(z,paste0(project_path,"/output/ARR_PB_14-28_days_as_ref.csv"))

df_cohort_vacc_g_both<-df_cohort_vacc_g
df_cohort_vacc_g <- filter(df_cohort_vacc_g_both,vacc_type=="AZ")
df_cohort_vacc_g<-df_cohort_vacc_g_both
df_cohort_vacc_g <- filter(df_cohort_vacc_g_both,vacc_type=="PB")

z<-as.data.frame(coef(summary(glm_pos1)))
write.csv(z,paste0(project_path,"/output/coef_meta_scotland.csv"))
#################

conditions <- c("Q_DIAG_AF", "Q_DIAG_ASTHMA", "Q_DIAG_BLOOD_CANCER", "Q_DIAG_CCF", 
                "Q_DIAG_CEREBRALPALSY","Q_DIAG_CHD","Q_DIAG_CIRRHOSIS", "Q_DIAG_CONGEN_HD",
                "Q_DIAG_COPD" , "Q_DIAG_DEMENTIA","Q_DIAG_DIABETES_1","Q_DIAG_DIABETES_2", 
                "Q_DIAG_EPILEPSY", "Q_DIAG_FRACTURE","Q_DIAG_NEURO","Q_DIAG_PARKINSONS",
                "Q_DIAG_PULM_HYPER","Q_DIAG_PULM_RARE","Q_DIAG_PVD", "Q_DIAG_RA_SLE",
                "Q_DIAG_RESP_CANCER","Q_DIAG_SEV_MENT_ILL","Q_DIAG_SICKLE_CELL","Q_DIAG_STROKE",
                "Q_DIAG_VTE", "Q_HOME_CAT","Q_LEARN_CAT","Q_DIAG_CKD","immuno")
ARR <- rbind()
for (i in 1:length(conditions)){#length(conditions)
  temp <- mutate(df_cohort_vacc_g, risk_gp = n_risk_gps)
  temp$risk_gp <- as.numeric(temp$risk_gp)
  temp <- mutate(temp, risk_gp=risk_gp-1)
  temp <- mutate(temp, risk_gp=risk_gp-get(conditions[i]))
  temp$risk_gp<-as.factor(temp$risk_gp)
  tt<-temp[conditions[i]]
  colnames(tt)<-"cond"
  temp<-cbind(temp,tt)
  temp$cond <- as.factor(temp$cond)
  temp$cond <- relevel(temp$cond, ref = "0")
  z.agg <- pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt) ~ +
                    cond+ period_f+Sex + age_gp+n_risk_gps+
                    bmi_gp+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                    n_tests_gp+ur_combined, 
                  data=temp, weight=ew , scale=365.25, data.frame=TRUE)
  glm_pos <- glm(event ~ offset(log(pyears)) +
                   cond+ period_f+Sex + age_gp+n_risk_gps+
                   bmi_gp+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                   n_tests_gp+ur_combined, 
                  family=poisson, data=z.agg$data)
  #summary(glm_pos)
  z<-as.data.frame(round(exp(glm_pos$coefficients),3))
  z<-cbind(z,round(exp(glm_pos$coefficients-coef(summary(glm_pos))[,2]*1.96),2))
  z<-cbind(z,round(exp(glm_pos$coefficients+coef(summary(glm_pos))[,2]*1.96),2))
  colnames(z)[1]<-"Adjusted Rate Ratio"
  colnames(z)[2]<-"LCI"
  colnames(z)[3]<-"UCI"
  rownames(z)[2]<-conditions[i]
  ARR <- rbind(ARR,z[2,])
  #print(z)
}
ARR<-mutate(ARR, ARR_LCI_UCI = paste(ARR[,1]," (",ARR[,2],"-",ARR[,3],")"))
colnames(ARR)[4]<-"ARR (LCI-UCI)"
write.csv(ARR,paste0(project_path,"/output/ARR_conds_14-28_days_as_ref.csv"))

#################

p_years<-pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt)~Sex, scale = 365.25 ,
                data = df_cohort_vacc_g, data.frame = TRUE, weights = ew)
sum(p_years$data$event)*1000/sum(p_years$data$pyears)

#################

z <- readRDS("/conf/EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid.rds")
z <- filter(z, !duplicated(EAVE_LINKNO))
z <- select(z, EAVE_LINKNO,Q_BMI)
df_cohort_vacc_g <- left_join(df_cohort_vacc_g,z, by="EAVE_LINKNO")

df_cohort_vacc_g <- df_cohort_vacc_g %>% 
  mutate(bmi_gp1 = cut(Q_BMI, breaks = c(-1, 20, 25, 30, 40, 51),
                       labels=c("<18.5", "18.5-24.9","25-29.9","30-39.9","40+")))
df_cohort_vacc_g$bmi_gp1 <- relevel(df_cohort_vacc_g$bmi_gp1, ref = "18.5-24.9")

z.agg <- pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt) ~ 
                  pv_period_f+period_f + Sex + age_gp+n_risk_gps+
                  bmi_gp1+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined, 
                data=df_cohort_vacc_g, weight=ew , scale=365.25, data.frame=TRUE)
glm_pos1 <- glm(event ~ offset(log(pyears)) +pv_period_f+period_f + Sex + age_gp+
                  n_risk_gps+bmi_gp1+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined, 
                family=poisson, data=z.agg$data)

z<-as.data.frame(round(exp(glm_pos1$coefficients),3))
z<-cbind(z,round(exp(glm_pos1$coefficients-coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z,round(exp(glm_pos1$coefficients+coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))
colnames(z)[1]<-"Adjusted Rate Ratio"
colnames(z)[2]<-"LCI"
colnames(z)[3]<-"UCI"
colnames(z)[4]<-"ARR (LCI-UCI)"
View(z)
#################
# # method 2
# z<-data.frame(seq(as.Date("2020-12-08"),as.Date("2022-03-28"),"day"))
# colnames(z)[1]<-"start_date"
# z <- mutate(z,year_week=yearweek(start_date))
# 
# df_cohort_vacc_g<-left_join(df_cohort_vacc_g,z)
# 
# # collapse urban rural for ARR analysis
# df_cohort_vacc_g$ur_combined <- as.factor(df_cohort_vacc_g$ur_combined)
# df_cohort_vacc_g$ur_combined <- relevel(df_cohort_vacc_g$ur_combined, ref = 1)
# 
# #df_cohort_vacc_g$vacc_type_3[is.na(df_cohort_vacc_g$vacc_type_3)]<-"NB"
# #df_cohort_vacc_g$vacc_type_3 <- as.factor(df_cohort_vacc_g$vacc_type_3)
# #df_cohort_vacc_g$vacc_type_3 <- relevel(df_cohort_vacc_g$vacc_type_3, ref = "NB")
# #df_cohort_vacc_g$booster_status <- as.factor(df_cohort_vacc_g$booster_status)
# #df_cohort_vacc_g$booster_status <- relevel(df_cohort_vacc_g$booster_status, ref = "nb")
# df_cohort_vacc_g$date_period <- as.factor(df_cohort_vacc_g$date_period)
# df_cohort_vacc_g$date_period <- relevel(df_cohort_vacc_g$date_period, ref = "date_preiod_1463")
# df_cohort_vacc_g$simd2020_sc_quintile <- as.factor(df_cohort_vacc_g$simd2020_sc_quintile)
# df_cohort_vacc_g$simd2020_sc_quintile <- relevel(df_cohort_vacc_g$simd2020_sc_quintile, ref = "5-Low")
# df_cohort_vacc_g$prior_pos_timing <- as.factor(df_cohort_vacc_g$prior_pos_timing)
# df_cohort_vacc_g$prior_pos_timing <- relevel(df_cohort_vacc_g$prior_pos_timing, ref = "No prior infection")
# 
# df_cohort_vacc_g <- mutate(df_cohort_vacc_g, p_years=as.numeric(end_time-start_time))
# #df_cohort_vacc_g$num_pos_avg <- df_cohort_vacc_g$num_pos_avg/100
# 
# df_cohort_vacc_g$prior_infect_monthgrp1[df_cohort_vacc_g$prior_infect_monthgrp=="No prior infection"]<-0
# df_cohort_vacc_g$prior_infect_monthgrp1[df_cohort_vacc_g$prior_infect_monthgrp=="0-3 month"]<-1
# df_cohort_vacc_g$prior_infect_monthgrp1[df_cohort_vacc_g$prior_infect_monthgrp=="3-6 month"]<-2
# df_cohort_vacc_g$prior_infect_monthgrp1[df_cohort_vacc_g$prior_infect_monthgrp=="6-9 month"]<-3
# df_cohort_vacc_g$prior_infect_monthgrp1[df_cohort_vacc_g$prior_infect_monthgrp==">9 month"]<-4
# df_cohort_vacc_g$prior_infect_monthgrp1<-as.factor(df_cohort_vacc_g$prior_infect_monthgrp1)
# 
# df_cohort_vacc_g$variant <- as.factor(df_cohort_vacc_g$variant)
# df_cohort_vacc_g$variant <- relevel(df_cohort_vacc_g$variant, ref = "Delta")
# 
# #df_cohort_vacc_g <- df_cohort_vacc_g %>% 
# #  mutate(bmi_gp1 = cut(bmi_impute, breaks = c(-1, 18.5, 24.9, 29.9, 39.9, 51),
# #                       labels=c("<18.5", "18.5-24.9","25-29.9","30-39.9","40+")))
# df_cohort_vacc_g$bmi_gp <- as.factor(df_cohort_vacc_g$bmi_gp)
# df_cohort_vacc_g$bmi_gp <- relevel(df_cohort_vacc_g$bmi_gp, ref = "18.5-24.9")
# ######### Modelling
# # crude analysis
# # sex
# # pyears - rate of event per year
# # glm_poisson - log rate ratios
# ################
# df_cohort_vacc_g1 <- df_cohort_vacc_g
# df_cohort_vacc_g <- filter(df_cohort_vacc_g, date_period!="date_preiod_0014")
# df_cohort_vacc_g<-filter(df_cohort_vacc_g,end_time>0)
# 
# z.agg <- pyears(Surv(start_time,end_time,event) ~ 
#                   date_period+year_week + Sex + age_gp+simd2020_sc_quintile, 
#                 data=df_cohort_vacc_g, weight=weight , scale=365.25, data.frame=TRUE)
# glm_pos1 <- glm(event ~ offset(log(pyears)) +date_period+year_week + Sex + age_gp+simd2020_sc_quintile, 
#                 family=poisson, data=z.agg$data)
# 
# glm_pos1 <- glm(event~offset(log(p_years))  +date_period+ Sex + age_gp + vacc_gap+
#                   prior_infect_monthgrp1+n_risk_gps +simd2020_sc_quintile+ start_date+
#                   ur_combined+bmi_gp  , family=poisson,weights = weight,data=df_cohort_vacc_g)
# glm_pos1 <- glm(event~offset(log(p_years))  +date_period+ Sex + age_gp +
#                   year_week +simd2020_sc_quintile, 
#                 family=poisson,weights = weight,data=df_cohort_vacc_g)
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
# write.csv(z,paste0(project_path,"/output/ARR_14-28_days_as_ref_b.csv"))
# write.csv(z,paste0(project_path,"/output/ARR_AZ_14-28_days_as_ref_b.csv"))
# write.csv(z,paste0(project_path,"/output/ARR_PB_14-28_days_as_ref_b.csv"))
# 
# df_cohort_vacc_g_both<-df_cohort_vacc_g
# df_cohort_vacc_g <- filter(df_cohort_vacc_g_both,vacc_type=="AZ")
# df_cohort_vacc_g<-df_cohort_vacc_g_both
# df_cohort_vacc_g <- filter(df_cohort_vacc_g_both,vacc_type=="PB")
###############################