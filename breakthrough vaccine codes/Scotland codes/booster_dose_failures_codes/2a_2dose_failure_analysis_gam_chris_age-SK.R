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
project_path <- paste0(Location,"EAVE/GPanalysis/progs/UA/booster_dose_failures")
xdf_full_covid_hosp_death <- readRDS(paste0(project_path,"/data/xdf_full_covid_hosp_death.RDS"))
#df_cohort_vacc <- readRDS(paste0(project_path,"/data/df_cohort_vacc_15-03-2022_gam.rds"))
#df_cohort_vacc_g <- readRDS(paste0(project_path,"/data/df_cohort_vacc_g_15-03-2022_gam.rds"))
##################
# data prep for the analysis
# method 1

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

df_cohort_vacc_g = mutate(df_cohort_vacc_g, Q_DIAG_CKD_LEVEL = as.factor(Q_DIAG_CKD_LEVEL))

# df_cohort_vacc_g$age_gp <- NULL
# df_cohort_vacc_g<-filter(df_cohort_vacc_g,ageYear<=74)
# df_cohort_vacc_g<-mutate(df_cohort_vacc_g,age_gp = 
#                       cut(ageYear, breaks = c(-Inf, 49,54,59,64,69,74),
#                           labels=c("18-49","50-54","55-59","60-64","65-69",
#                                    "70-74")))
# df_cohort_vacc_g$age_gp <- as.factor(df_cohort_vacc_g$age_gp)
# df_cohort_vacc_g$age_gp <- relevel(df_cohort_vacc_g$age_gp, ref = "18-49")

# df_cohort_vacc_g$age_gp <- NULL
# df_cohort_vacc_g<-filter(df_cohort_vacc_g,ageYear<=64)
# df_cohort_vacc_g<-mutate(df_cohort_vacc_g,age_gp = 
#                            cut(ageYear, breaks = c(-Inf, 49,54,59,64),
#                                labels=c("18-49","50-54","55-59","60-64")))
# df_cohort_vacc_g$age_gp <- as.factor(df_cohort_vacc_g$age_gp)
# df_cohort_vacc_g$age_gp <- relevel(df_cohort_vacc_g$age_gp, ref = "18-49")

#df_cohort_vacc_g$age_gp <- NULL
#df_cohort_vacc_g<-filter(df_cohort_vacc_g,ageYear<50)

#df_cohort_vacc_g$age_gp <- NULL
#df_cohort_vacc_g<-filter(df_cohort_vacc_g,ageYear<40)

df_cohort_vacc_g$age_gp <- NULL
df_cohort_vacc_g<-filter(df_cohort_vacc_g,ageYear>=75)

df_cohort_vacc_g<-mutate(df_cohort_vacc_g,age_gp =
                      cut(ageYear, breaks = c(75,80, 110),
                          labels=c("75-79","80+"), include.lowest = TRUE))


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

# z<-as.data.frame(coef(summary(glm_pos1)))
# write.csv(z,paste0(project_path,"/output/ARR_14-28_days_as_ref_b_under_75.csv"))
# z<-as.data.frame(coef(summary(glm_pos1)))
# write.csv(z,paste0(project_path,"/output/ARR_14-28_days_as_ref_b_under_65.csv"))
# z<-as.data.frame(coef(summary(glm_pos1)))
# write.csv(z,paste0(project_path,"/output/ARR_14-28_days_as_ref_b_under_50.csv"))
# z<-as.data.frame(coef(summary(glm_pos1)))
# write.csv(z,paste0(project_path,"/output/ARR_14-28_days_as_ref_b_under_40.csv"))
z<-as.data.frame(coef(summary(glm_pos1)))
write.csv(z,paste0(project_path,"/output/ARR_14-28_days_as_ref_b_over_75.csv"))
# z<-as.data.frame(coef(summary(glm_pos1)))
# write.csv(z,paste0(project_path,"/output/ARR_Mo_14-28_days_as_ref_b.csv"))
# z<-as.data.frame(coef(summary(glm_pos1)))
# write.csv(z,paste0(project_path,"/output/ARR_PB_14-28_days_as_ref_b.csv"))
# 
# df_cohort_vacc_g_both<-df_cohort_vacc_g
# df_cohort_vacc_g <- filter(df_cohort_vacc_g_both,vacc_type_3=="Mo")
# df_cohort_vacc_g<-df_cohort_vacc_g_both
# df_cohort_vacc_g <- filter(df_cohort_vacc_g_both,vacc_type_3=="PB")
#################

conditions <- c("Q_DIAG_AF", "Q_DIAG_ASTHMA", "Q_DIAG_BLOOD_CANCER", "Q_DIAG_CCF", 
                "Q_DIAG_CEREBRALPALSY","Q_DIAG_CHD","Q_DIAG_CIRRHOSIS", "Q_DIAG_CONGEN_HD",
                "Q_DIAG_COPD" , "Q_DIAG_DEMENTIA","Q_DIAG_DIABETES_1","Q_DIAG_DIABETES_2", 
                "Q_DIAG_EPILEPSY", "Q_DIAG_FRACTURE","Q_DIAG_NEURO","Q_DIAG_PARKINSONS",
                "Q_DIAG_PULM_HYPER","Q_DIAG_PULM_RARE","Q_DIAG_PVD", "Q_DIAG_RA_SLE",
                "Q_DIAG_RESP_CANCER","Q_DIAG_SEV_MENT_ILL","Q_DIAG_SICKLE_CELL","Q_DIAG_STROKE",
                "Q_DIAG_VTE", "Q_HOME_CAT","Q_LEARN_CAT","Q_DIAG_CKD","immuno","Q_DIAG_CKD_LEVEL")
ARR <- rbind()
for (i in 1:length(conditions)){#length(conditions)
  
  print(i)
  
  tryCatch({
  
  temp <- mutate(df_cohort_vacc_g, risk_gp = n_risk_gps)
  # This is not needed - n_risk_gps is a factor
  # temp$risk_gp <- as.numeric(temp$risk_gp)
  # temp <- mutate(temp, risk_gp=risk_gp-1)
  # temp <- mutate(temp, risk_gp=risk_gp-get(conditions[i]))
  # temp$risk_gp<-as.factor(temp$risk_gp)
  tt<-temp[conditions[i]]
  colnames(tt)<-"cond"
  temp<-cbind(temp,tt)
  temp$cond <- as.factor(temp$cond)
  temp$cond <- relevel(temp$cond, ref = "0")
  z.agg <- pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt) ~ +
                    cond+ period_f+Sex + age_gp+bmi_gp+simd2020_sc_quintile+
                    pos_before_start_u+vacc_gap+n_tests_gp+ur_combined, 
                  data=temp, weight=ew , scale=365.25, data.frame=TRUE)
  glm_pos <- glm(event ~ offset(log(pyears)) +
                   cond+ period_f+Sex + age_gp+bmi_gp+simd2020_sc_quintile+
                   pos_before_start_u+vacc_gap+n_tests_gp+ur_combined, 
                  family=poisson, data=z.agg$data)
  #summary(glm_pos)
  z<-as.data.frame(round(exp(glm_pos$coefficients),3))
  z<-cbind(z,round(exp(glm_pos$coefficients-coef(summary(glm_pos))[,2]*1.96),2))
  z<-cbind(z,round(exp(glm_pos$coefficients+coef(summary(glm_pos))[,2]*1.96),2))
  colnames(z)[1]<-"Adjusted Rate Ratio"
  colnames(z)[2]<-"LCI"
  colnames(z)[3]<-"UCI"
  
  if( conditions[i] == "Q_DIAG_CKD_LEVEL"){
    rownames(z)[2:4]<-c('Q_DIAG_CKD_LEVEL_3', 'Q_DIAG_CKD_LEVEL_4', 'Q_DIAG_CKD_LEVEL_5')
    ARR <- rbind(ARR,z[2:4,])
  } else{
  rownames(z)[2]<-conditions[i]
  ARR <- rbind(ARR,z[2,])
  }
  
  }, error = function(e){})
  
  #print(z)
}
ARR<-mutate(ARR, ARR_LCI_UCI = paste(ARR[,1]," (",ARR[,2],"-",ARR[,3],")"))
colnames(ARR)[4]<-"ARR (LCI-UCI)"
#write.csv(ARR,paste0(project_path,"/output/ARR_conds_14-28_days_as_ref_b_under_75.csv"))
write.csv(ARR,paste0(project_path,"/output/ARR_conds_14-28_days_as_ref_b_over_75.csv"))

#################

p_years<-pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt)~Sex, scale = 365.25 ,
                data = df_cohort_vacc_g, data.frame = TRUE, weights = ew)
sum(p_years$data$event)*1000/sum(p_years$data$pyears)

###############################