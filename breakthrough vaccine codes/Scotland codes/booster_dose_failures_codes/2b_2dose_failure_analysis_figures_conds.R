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
project_path <- paste0(Location,"EAVE/GPanalysis/progs/UA/booster_dose_failures")
df_cohort_vacc1 <- readRDS(paste0(project_path,"/data/xdf_booster_pop.RDS"))
xdf_full_covid_hosp_death <- readRDS(paste0(project_path,"/data/xdf_full_covid_hosp_death.RDS"))
##################
# data prep for the analysis

df_cohort_vacc_g <- filter(xdf_full_covid_hosp_death, pv_period_f!="uv")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v1_0:3")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v1_4+")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v2_0:3")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v2_4+")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v3_0:1")

# tt1<-filter(df_cohort_vacc_g,event_date>"2021-12-14"&event==1)
# tt<-filter(df_cohort_vacc_g,event==0)
# tt<-rbind(tt,tt1)
# tt<-arrange(tt,EAVE_LINKNO,period_f)
# tt$pv_period_f<-as.character(tt$pv_period_f)

df_cohort_vacc_g$pv_period_f <- as.factor(df_cohort_vacc_g$pv_period_f)
df_cohort_vacc_g$pv_period_f <- relevel(df_cohort_vacc_g$pv_period_f, ref = "v3_2:4")
df_cohort_vacc_g$simd2020_sc_quintile <- as.factor(df_cohort_vacc_g$simd2020_sc_quintile)
df_cohort_vacc_g$simd2020_sc_quintile <- relevel(df_cohort_vacc_g$simd2020_sc_quintile, ref = "5-Low")

df_cohort_vacc_g$bmi_gp <- as.factor(df_cohort_vacc_g$bmi_gp)
df_cohort_vacc_g$bmi_gp <- relevel(df_cohort_vacc_g$bmi_gp, ref = "18.5-24.9")

df_cohort_vacc_g$vacc_gap <- as.factor(df_cohort_vacc_g$vacc_gap)
df_cohort_vacc_g$vacc_gap <- relevel(df_cohort_vacc_g$vacc_gap, ref = "<7wk")

df_cohort_vacc_g$ur_combined <- as.factor(df_cohort_vacc_g$ur_combined)
df_cohort_vacc_g$ur_combined <- relevel(df_cohort_vacc_g$ur_combined, ref = "2")

df_cohort_vacc_g$n_tests_gp <- as.character(df_cohort_vacc_g$n_tests_gp)
df_cohort_vacc_g$n_tests_gp[df_cohort_vacc_g$n_tests_gp=="3"] <- "3-4"
df_cohort_vacc_g$n_tests_gp[df_cohort_vacc_g$n_tests_gp=="4"] <- "3-4"
df_cohort_vacc_g$n_tests_gp[df_cohort_vacc_g$n_tests_gp=="10-19"] <- "10+"
df_cohort_vacc_g$n_tests_gp[df_cohort_vacc_g$n_tests_gp=="20+"] <- "10+"
df_cohort_vacc_g$n_tests_gp <- as.factor(df_cohort_vacc_g$n_tests_gp)
df_cohort_vacc_g$n_tests_gp <- relevel(df_cohort_vacc_g$n_tests_gp, ref = "0")
#####
df_cohort_vacc1$n_tests_gp <- as.character(df_cohort_vacc1$n_tests_gp)
df_cohort_vacc1$n_tests_gp[df_cohort_vacc1$n_tests_gp=="3"] <- "3-4"
df_cohort_vacc1$n_tests_gp[df_cohort_vacc1$n_tests_gp=="4"] <- "3-4"
df_cohort_vacc1$n_tests_gp[df_cohort_vacc1$n_tests_gp=="10-19"] <- "10+"
df_cohort_vacc1$n_tests_gp[df_cohort_vacc1$n_tests_gp=="20+"] <- "10+"
df_cohort_vacc1$n_tests_gp <- as.factor(df_cohort_vacc1$n_tests_gp)

df_cohort_vacc1 <- filter(df_cohort_vacc1, vacc_type!="Mo")
df_cohort_vacc1 <- filter(df_cohort_vacc1, vacc_type_3!="AZ")

#######
df_cohort_vacc_g$age_gp <- as.factor(df_cohort_vacc_g$age_gp)
df_cohort_vacc_g$age_gp <- relevel(df_cohort_vacc_g$age_gp, ref = "18-49")

# df_cohort_vacc_g$age_gp <- NULL
# df_cohort_vacc_g<-filter(df_cohort_vacc_g,ageYear<=74)
# df_cohort_vacc_g<-mutate(df_cohort_vacc_g,age_gp =
#                       cut(ageYear, breaks = c(-Inf, 49,54,59,64,69,74),
#                           labels=c("18-49","50-54","55-59","60-64","65-69",
#                                    "70-74")))
# df_cohort_vacc_g$age_gp <- as.factor(df_cohort_vacc_g$age_gp)
# df_cohort_vacc_g$age_gp <- relevel(df_cohort_vacc_g$age_gp, ref = "18-49")
# 
# df_cohort_vacc1$age_gp <- NULL
# df_cohort_vacc1<-filter(df_cohort_vacc1,ageYear<=74)
# df_cohort_vacc1<-mutate(df_cohort_vacc1,age_gp =
#                            cut(ageYear, breaks = c(-Inf, 49,54,59,64,69,74),
#                                labels=c("18-49","50-54","55-59","60-64","65-69",
#                                         "70-74")))

# df_cohort_vacc_g$age_gp <- NULL
# df_cohort_vacc_g<-filter(df_cohort_vacc_g,ageYear>=75)
# 
# df_cohort_vacc1$age_gp <- NULL
# df_cohort_vacc1<-filter(df_cohort_vacc1,ageYear>=75)
# df_cohort_vacc1<-mutate(df_cohort_vacc1,age_gp =
#                           cut(ageYear, breaks = c(-Inf, Inf),
#                               labels=c("75+")))
###############################
z_df <- df_cohort_vacc1 %>%
  dplyr::select(Q_DIAG_AF, Q_DIAG_ASTHMA, Q_DIAG_BLOOD_CANCER, Q_DIAG_CCF, 
                Q_DIAG_CEREBRALPALSY,Q_DIAG_CHD,Q_DIAG_CIRRHOSIS,Q_DIAG_CONGEN_HD,
                Q_DIAG_COPD , Q_DIAG_DEMENTIA, Q_DIAG_DIABETES_1,Q_DIAG_DIABETES_2, 
                Q_DIAG_EPILEPSY, Q_DIAG_FRACTURE,Q_DIAG_NEURO,Q_DIAG_PARKINSONS,
                Q_DIAG_PULM_HYPER,Q_DIAG_PULM_RARE,Q_DIAG_PVD, Q_DIAG_RA_SLE,
                Q_DIAG_RESP_CANCER,Q_DIAG_SEV_MENT_ILL,Q_DIAG_SICKLE_CELL,
                Q_DIAG_STROKE,Q_DIAG_VTE, Q_HOME_CAT,Q_LEARN_CAT,Q_DIAG_CKD,
                immuno,Q_DIAG_CKD_LEVEL,event) %>% 
  pivot_longer(cols=Q_DIAG_AF:Q_DIAG_CKD_LEVEL) 

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
chars<-c("Q_DIAG_AF", "Q_DIAG_ASTHMA", "Q_DIAG_BLOOD_CANCER", "Q_DIAG_CCF", 
         "Q_DIAG_CEREBRALPALSY","Q_DIAG_CHD","Q_DIAG_CIRRHOSIS", "Q_DIAG_CONGEN_HD",
         "Q_DIAG_COPD" , "Q_DIAG_DEMENTIA","Q_DIAG_DIABETES_1","Q_DIAG_DIABETES_2", 
         "Q_DIAG_EPILEPSY", "Q_DIAG_FRACTURE","Q_DIAG_NEURO","Q_DIAG_PARKINSONS",
         "Q_DIAG_PULM_HYPER","Q_DIAG_PULM_RARE","Q_DIAG_PVD", "Q_DIAG_RA_SLE",
         "Q_DIAG_RESP_CANCER","Q_DIAG_SEV_MENT_ILL","Q_DIAG_SICKLE_CELL","Q_DIAG_STROKE",
         "Q_DIAG_VTE", "Q_HOME_CAT","Q_LEARN_CAT","Q_DIAG_CKD","immuno","Q_DIAG_CKD_LEVEL")
event_list<-rbind()
for (i in 1:length(chars)){
  p_years<-pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt)~get(chars[i]), scale = 365.25 ,
                  data = df_cohort_vacc_g, data.frame = TRUE,weights = ew)
  z<-data.frame(p_years$data$`get(chars[i])`)
  colnames(z)[1]<-"value"
  z<-mutate(z,name=chars[i])
  z1<-data.frame(p_years$data$pyears)
  colnames(z1)[1]<-"person_years"
  z<-bind_cols(z,z1)
  z1<-data.frame(round((p_years$data$event)*1000/(p_years$data$pyears),1))
  colnames(z1)[1]<-"rate_per_1000_yrs"
  z<-bind_cols(z,z1)
  event_list<-rbind(event_list,z)
}
event_list$value <- as.character(event_list$value)
z_df$value <- as.character(z_df$value)
z_df <- left_join(z_df,event_list)

results_f <- select(z_df, name, value, total_vacc_adm, Percent_of_vacc_adm, 
                    both_vacc_event, person_years, rate_per_1000_yrs)

# z_df <- cbind(z_df, paste(z_df$total_vacc_adm," (",z_df$Percent_of_vacc_adm,")", sep = ""))
# z_df <- cbind(z_df, paste(z_df$both_vacc_event," (",z_df$rate_per_1000_yrs,")", sep = ""))
# 
# colnames(z_df)[8] <- "total vacc (n, %)"
# colnames(z_df)[9] <- "severe outcome (n, rate per 1000 person years)"
# results <- select(z_df,name, value,`total vacc (n, %)`,
#                   `severe outcome (n, rate per 1000 person years)`)
# results_f <- results

write.csv(results_f,paste0(project_path,"/output/cohort_description_b_QCovid_all.csv"))
write.csv(results_f,paste0(project_path,"/output/cohort_description_b_QCovid_under_75.csv"))
write.csv(results_f,paste0(project_path,"/output/cohort_description_b_QCovid_over_75.csv"))
##########
#############

p_years<-pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt)~Sex, scale = 365.25 ,
                data = df_cohort_vacc_g, data.frame = TRUE, weights = ew)
sum(p_years$data$event)*1000/sum(p_years$data$pyears)

tt<-filter(df_cohort_vacc1,event==1)
table(tt$covid_hosp)
table(tt$covid_death)
tt<-filter(df_cohort_vacc1,covid_death==1)
table(tt$covid_hosp)
