library(tidyverse)
library(plyr)
library(meta)

coef_data_s <- read.csv("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/output/coef_meta_scotland.csv")
coef_data_w <- read.csv("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/output/t_pois_coef_doseb_full.csv")

coef_data_s <- coef_data_s[,c(1,2,3)]
coef_data_w <- coef_data_w[,c(1,2,3,6)]

coef_data_s <- mutate(coef_data_s,country="scotland")
colnames(coef_data_s)[1]<-"term"
colnames(coef_data_s)[2]<-"estimate"
colnames(coef_data_s)[3]<-"std_error"

colnames(coef_data_w)[1]<-"term"
colnames(coef_data_w)[2]<-"estimate"
colnames(coef_data_w)[3]<-"std_error"

coef_data_s <- coef_data_s[c(37:49),]
coef_data_w <- coef_data_w[c(4,9:20),]

names<-c(  "sexM","age_gp50-54","age_gp55-59",
           "age_gp60-64","age_gp65-69","age_gp70-74",
           "age_gp75-79","age_gp80+","n_risk_gps1",
           "n_risk_gps2","n_risk_gps3","n_risk_gps4",
           "n_risk_gps5+")
coef_data_w[,1]<-names

coef_data_ma <- rbind(coef_data_s,coef_data_w)

###############



###############