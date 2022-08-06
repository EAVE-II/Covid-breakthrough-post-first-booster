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

################
temp <- as.data.frame(c(1.23,2.29,3.42,0.62,0.68,0.90,0.84))
categories <- c("49-84","84-168","168o","b_Mo_14","b_Mo_28","b_PB_14","b_PB_28")

z<-as.data.frame(round(exp(glm_pos1$coefficients),3))
z<-cbind(z,round(exp(glm_pos1$coefficients-coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z,round(exp(glm_pos1$coefficients+coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))
colnames(z)[1]<-"Adjusted Rate Ratio"
colnames(z)[2]<-"LCI"
colnames(z)[3]<-"UCI"
colnames(z)[4]<-"ARR (LCI-UCI)"
View(z)
###############################