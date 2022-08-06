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
df_cohort_vacc1 <- readRDS(paste0(project_path,"/data/xdf_full_cohort.RDS"))
xdf_full_covid_hosp_death <- readRDS(paste0(project_path,"/data/xdf_full_covid_hosp_death.RDS"))
##################
# data prep for the analysis
df_cohort_vacc1$ur6_2016_name[is.na(df_cohort_vacc1$ur6_2016_name)] <- "missing"
xdf_full_covid_hosp_death$ur6_2016_name[is.na(xdf_full_covid_hosp_death$ur6_2016_name)] <- "missing"
df_cohort_vacc1$ur_combined[is.na(df_cohort_vacc1$ur_combined)] <- "missing"
xdf_full_covid_hosp_death$ur_combined[is.na(xdf_full_covid_hosp_death$ur_combined)] <- "missing"

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
#####
df_cohort_vacc1$n_tests_gp <- as.character(df_cohort_vacc1$n_tests_gp)
df_cohort_vacc1$n_tests_gp[df_cohort_vacc1$n_tests_gp=="3"] <- "3-4"
df_cohort_vacc1$n_tests_gp[df_cohort_vacc1$n_tests_gp=="4"] <- "3-4"
df_cohort_vacc1$n_tests_gp[df_cohort_vacc1$n_tests_gp=="10-19"] <- "10+"
df_cohort_vacc1$n_tests_gp[df_cohort_vacc1$n_tests_gp=="20+"] <- "10+"
df_cohort_vacc1$n_tests_gp <- as.factor(df_cohort_vacc1$n_tests_gp)

df_cohort_vacc1 <- filter(df_cohort_vacc1, vacc_type!="Mo")
df_cohort_vacc1 <- filter(df_cohort_vacc1, vacc_type_3=="Mo" |
                             vacc_type_3=="PB" | is.na(vacc_type_3))
###############################
z_df <- df_cohort_vacc1 %>%
  dplyr::select(age_gp,Sex, simd2020_sc_quintile, ur6_2016_name, n_risk_gps,
                n_tests_gp,vacc_gap, bmi_gp,pos_before_start_u,event) %>% 
  pivot_longer(cols=age_gp:pos_before_start_u) 

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
chars<-c('age_gp','Sex','simd2020_sc_quintile','ur6_2016_name',
         'n_risk_gps','n_tests_gp','vacc_gap', 'bmi_gp','pos_before_start_u')
event_list<-rbind()
for (i in 1:9){
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
z_df <- left_join(z_df,event_list)

results_f <- select(z_df, name, value, total_vacc_adm, Percent_of_vacc_adm, 
               both_vacc_event, person_years, rate_per_1000_yrs)

# z_df <- cbind(z_df, paste(z_df$total_vacc_adm," (",z_df$Percent_of_vacc_adm,")", sep = ""))
# z_df <- cbind(z_df, paste(z_df$both_vacc_event," (",z_df$rate_per_1000_yrs,")", sep = ""))
# 
# 
# colnames(z_df)[8] <- "total vacc (n, %)"
# colnames(z_df)[9] <- "severe outcome (n, rate per 1000 person years)"
# results <- select(z_df,name, value,`total vacc (n, %)`,
#                   `severe outcome (n, rate per 1000 person years)`)
# results_f <- results

write.csv(results_f,paste0(project_path,"/output/cohort_description.csv"))
write.csv(results_f,paste0(project_path,"/output/cohort_description_AZ.csv"))
write.csv(results_f,paste0(project_path,"/output/cohort_description_PB.csv"))
##########

df_cohort_vacc_g_both<-df_cohort_vacc_g
df_cohort_vacc_both <- df_cohort_vacc1

df_cohort_vacc_g <- filter(df_cohort_vacc_g_both,vacc_type=="AZ")
df_cohort_vacc1 <- filter(df_cohort_vacc_both,vacc_type=="AZ")
df_cohort_vacc_g<-df_cohort_vacc_g_both
df_cohort_vacc1<-df_cohort_vacc_both
df_cohort_vacc_g <- filter(df_cohort_vacc_g_both,vacc_type=="PB")
df_cohort_vacc1 <- filter(df_cohort_vacc_both,vacc_type=="PB")
##########
##########
p_years<-pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt)~Sex, scale = 365.25 ,
                data = df_cohort_vacc_g, data.frame = TRUE, weights = ew)
sum(p_years$data$event)*1000/sum(p_years$data$pyears)

tt<-filter(df_cohort_vacc1,event==1)
table(tt$covid_hosp)
table(tt$covid_death)
tt<-filter(df_cohort_vacc1,covid_death==1)
table(tt$covid_hosp)

df_cohort_vacc_g <- df_cohort_vacc_g1

a_end <- as.Date("2021-05-19")
df_cohort_vacc_g1 <- df_cohort_vacc_g1 %>% filter(start_date < a_end)

df_cohort_vacc_g1 <- df_cohort_vacc_g
df_cohort_vacc_g1 <- df_cohort_vacc_g1 %>% filter(start_date >= a_end &
                                                    start_date<"2021-12-15") 

df_cohort_vacc_g1 <- df_cohort_vacc_g
df_cohort_vacc_g1 <- df_cohort_vacc_g1 %>% filter(start_date >= "2021-12-15") 
###############