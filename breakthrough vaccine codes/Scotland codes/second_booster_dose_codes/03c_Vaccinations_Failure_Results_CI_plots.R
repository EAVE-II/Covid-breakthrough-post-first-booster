df_cc <- readRDS("/conf/EAVE/GPanalysis/progs/RM/Vaccine/Vaccine_waning/output/df_cc_ps_matches_m2_death_hosp.rds")
### Cumulative incidence plots ####

df_cohort_vacc_g <- df_cohort_vacc_g1
#library(survminer)
z_title <- "COVID-19 serious events"

surv_fit<-survfit(Surv(start_time,end_time,event)~1, data = df_cohort_vacc_g,weights = weight)
plot(surv_fit,fun="event")
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk")

df_cohort_vacc_g$vacc_type<-as.factor(df_cohort_vacc_g$vacc_type)
cols<-c("#8E063B","#023FA5")
cols<-c("#CC0033","#000099","#CC6600","#003333")
#### by vaccine
z <- survfit(Surv(start_time,end_time, event) ~ vacc_type, data=df_cohort_vacc_g, weights = weight)
plot(z, fun="event", col=cols[1:2], conf.int=T, ylim = c(0,0.005))
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",sub = "by vaccine type")
#legend(x=70,y=0.0002,legend=levels(df_cohort_vacc_g$vacc_type), title = "Vaccine type", lty=1, col=cols[1:2], cex=0.8, bg="gray90")
legend("topleft",legend=c("ChAdOx1","BNT162b2"), title = "Second Dose", lty=1, col=cols[1:2], cex=0.6, bg="gray90")
#summary(survfit(Surv(start_time,end_time, event) ~ vacc_type, data=df_cohort_vacc_g, weights = weight),times=365.25)

tt<-as.data.frame(table(cut(x=z$n.risk,breaks = seq(from=0,to=400,by=50))))
ggsurvplot(z,risk.table = TRUE,data = df_cohort_vacc_g)

z <- survfit(Surv(start_time,end_time, event) ~ vacc_type_3, data=df_cohort_vacc_g, weights = weight)
plot(z, fun="event", col=cols[1:2], conf.int=F, ylim = c(0,0.0012))
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",sub = "by vaccine type")
#legend(x=70,y=0.0002,legend=levels(df_cohort_vacc_g$vacc_type), title = "Vaccine type", lty=1, col=cols[1:2], cex=0.8, bg="gray90")
legend("topleft",legend=c("mRNA-1273","BNT162b2"), title = "Booster dose", lty=1, col=cols[1:2], cex=0.6, bg="gray90")


# plot(surv_fit,fun="event")
# par(new=TRUE)
# plot(z, fun="event", col=cols[1:2], conf.int=T)
# title(main=z_title, xlab="days from vaccination",ylab="cumulative risk")

# AZ
# par(mfrow=c(1,2))
# z_az <- survfit(Surv(start_time,end_time, event) ~ vacc_type, data=df_cohort_vacc_g,subset = vacc_type=="AZ")
# plot(z_az, fun="event", col=cols[1], conf.int=T)
# title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",sub = "AZ")
# legend(x=70,y=0.0002,legend=levels(df_cohort_vacc_g$vacc_type), title = "Vaccine type", lty=1, col=cols[1],lwd=2, cex=0.8, bg="gray90")
# 
#PB
# z_pb <- survfit(Surv(start_time,end_time, event) ~ vacc_type, data=df_cohort_vacc_g,subset = vacc_type=="PB")
# plot(z_pb, fun="event", col=cols[2], conf.int=T)
# title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",sub = "PB")
# legend(x=70,y=0.0002,legend=levels(df_cohort_vacc_g$vacc_type), title = "Vaccine type", lty=1, col=cols[2],lwd=2, cex=0.8, bg="gray90")

#### by age
z_age <- survfit(Surv(start_time,end_time, event) ~ age_gp, data=df_cohort_vacc_g)
plot(z_age, fun="event", col=cols, conf.int=F)
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",sub = "by age groups")
legend("topleft",legend=levels(df_cohort_vacc_g$age_gp), title = "Age Group", lty=1, col=cols, cex=0.6, bg="gray90")


par(mfrow=c(1,2))
plot(surv_fit,fun="event")
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk")

plot(z, fun="event", col=cols[1:2], conf.int=T)
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",sub = "by vaccine type")
#legend("topleft",legend=levels(df_cohort_vacc_g$vacc_type), title = "Vaccine type", lty=1, col=cols, cex=0.8, bg="gray90")

#### by num of risk groups

cols<-c("#CC0033","#000099","#CC6600","#003333","coral4","gray0")
#### by vaccine
z_nrg <- survfit(Surv(start_time,end_time, event) ~ n_risk_gps, data=df_cohort_vacc_g)
plot(z_nrg, fun="event", col=cols, conf.int=F)
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",sub = "by number of risk groups")
#legend(x=70,y=0.0002,legend=levels(df_cohort_vacc_g$vacc_type), title = "Vaccine type", lty=1, col=cols[1:2], cex=0.8, bg="gray90")
legend("topleft",legend=levels(df_cohort_vacc_g$n_risk_gps), title = "Num of Risk Groups", lty=1, col=cols, cex=0.6, bg="gray90")

#### by num of risk groups

cols<-c("#CC0033","#000099","#CC6600","#003333","coral4","gray0")
#### by simd
df_cohort_vacc_g1<-filter(df_cohort_vacc_g,simd2020_sc_quintile!="NA")
df_cohort_vacc_g1$simd2020_sc_quintile<-as.factor(df_cohort_vacc_g1$simd2020_sc_quintile)
z_simd <- survfit(Surv(start_time,end_time, event) ~ simd2020_sc_quintile, data=df_cohort_vacc_g1)
plot(z_simd, fun="event", col=cols, conf.int=F)
axis(1,font=2)
axis(2,font=2)
title(main=z_title, xlab="days from vaccination",ylab="cumulative risk",sub = "by SIMD")
font(side = 1, font=2)
legend("topleft",legend=levels(df_cohort_vacc_g1$simd2020_sc_quintile), title = "SIMD", lty=1, col=cols, cex=0.6, bg="gray90")


###
# histogram of age
EAVE_cohort<-select(EAVE_cohort,EAVE_LINKNO,age_gp)
colnames(EAVE_cohort)[2]<-"age_gp1"
df_cohort_vacc<-left_join(df_cohort_vacc,EAVE_cohort)
df_cohort_vacc_g1<-filter(df_cohort_vacc,event_all==1)
vis_data<-data.frame(table(df_cohort_vacc_g1$age_gp1,df_cohort_vacc_g1$vacc_type))
colnames(vis_data)[1]<-"Age_group"
colnames(vis_data)[2]<-"Vaccine_type"
colnames(vis_data)[3]<-"Number_vaccinated"
vis_data<-filter(vis_data,Age_group!="0-4"&Age_group!="5-9"&Age_group!="10-14")
vis_data$Age_group<-as.character(vis_data$Age_group)
vis_data<-mutate(vis_data,Age_group=if_else(Age_group=="15-19","<=19",Age_group))
black.bold<-element_text(face = "bold", color = "black",size = 13)
ggplot(vis_data)+
  geom_col(aes(x=Age_group,y=Number_vaccinated,fill=Vaccine_type))+
  ylab("COVID-19 hospitalisations and deaths (n)")+
  xlab("Age groups")+
  scale_y_continuous(labels = comma)+
  scale_x_discrete(labels=c("<=19","20","25","30","35","40","45","50","55","60","65","70","75",
                            "80","85","90+"))+
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),axis.text.x = black.bold,
        axis.text.y = black.bold, axis.title = black.bold,
        legend.text = black.bold, legend.title = black.bold)

###
# histogram of age for hospitalisation
df_cohort_vacc_g1<-filter(df_cohort_vacc,event_all==1&!is.na(date_hosp_covid))
vis_data<-data.frame(table(df_cohort_vacc_g1$age_gp1,df_cohort_vacc_g1$vacc_type))
colnames(vis_data)[1]<-"Age_group"
colnames(vis_data)[2]<-"Vaccine_type"
colnames(vis_data)[3]<-"Number_vaccinated"
vis_data<-filter(vis_data,Age_group!="0-4"&Age_group!="5-9"&Age_group!="10-14")
vis_data$Age_group<-as.character(vis_data$Age_group)
vis_data<-mutate(vis_data,Age_group=if_else(Age_group=="15-19","<=19",Age_group))
ggplot(vis_data)+
  geom_col(aes(x=Age_group,y=Number_vaccinated,fill=Vaccine_type))+
  ylab("COVID-19 hospitalisations and deaths (n)")+
  xlab("Age groups")+
  scale_y_continuous(labels = comma)+
  scale_x_discrete(labels=c("<=19","20","25","30","35","40","45","50","55","60","65","70","75",
                            "80","85","90+"))+
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),axis.text.x = black.bold,
        axis.text.y = black.bold, axis.title = black.bold,
        legend.text = black.bold, legend.title = black.bold)

###
# histogram of age for death
df_cohort_vacc_g1<-filter(df_cohort_vacc,event_all==1&!is.na(covid_death_date))
vis_data<-data.frame(table(df_cohort_vacc_g1$age_gp1,df_cohort_vacc_g1$vacc_type))
colnames(vis_data)[1]<-"Age_group"
colnames(vis_data)[2]<-"Vaccine_type"
colnames(vis_data)[3]<-"Number_vaccinated"
vis_data<-filter(vis_data,Age_group!="0-4"&Age_group!="5-9"&Age_group!="10-14")
vis_data$Age_group<-as.character(vis_data$Age_group)
vis_data<-mutate(vis_data,Age_group=if_else(Age_group=="15-19","<=19",Age_group))
ggplot(vis_data)+
  geom_col(aes(x=Age_group,y=Number_vaccinated,fill=Vaccine_type))+
  ylab("COVID-19 hospitalisations and deaths (n)")+
  xlab("Age groups")+
  scale_y_continuous(labels = comma)+
  scale_x_discrete(labels=c("<=19","20","25","30","35","40","45","50","55","60","65","70","75",
                            "80","85","90+"))+
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),axis.text.x = black.bold,
        axis.text.y = black.bold, axis.title = black.bold,
        legend.text = black.bold, legend.title = black.bold)


df_cohort_vacc_g1<-filter(df_cohort_vacc,event_all==1)
vis_data<-data.frame(table(df_cohort_vacc_g1$Sex,df_cohort_vacc_g1$vacc_type))
colnames(vis_data)[1]<-"Sex"
colnames(vis_data)[2]<-"Vaccine_type"
colnames(vis_data)[3]<-"Number_of_events"

ggplot(vis_data, aes(x = Sex, y = Number_of_events, fill = Vaccine_type)) +
  geom_bar(stat="identity",position = position_dodge(),width = 0.5) +
  theme_minimal()+
  labs(x = "Sex", y="Number of events")+
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

##################
names(cohort)
library(tsibble)
today <- Sys.Date()
file_date <- today - 2
file_date <- str_remove_all(file_date, "[-]") %>% as.numeric()
previous_file <- file_date-7
pandemic_start <- as.Date("2020-12-30")


cohort_test <- df_cohort_vacc %>% filter(!is.na(covid_hosp_death_date))
#Add iso_week to allow aggregation
cohort_test$iso_week <- yearweek(cohort_test$covid_hosp_death_date)

# Aggregate by iso_week with basic descriptives
tests_by_week <- cohort_test %>% group_by(iso_week) %>% 
  summarise(wk_comm = first(date(iso_week)),
            tests = n(),
            individuals = n_distinct(EAVE_LINKNO),
            positive = sum(covid_hosp_death_status == 1, na.rm = T),
            min_age = min(ageYear),
            max_age = max(ageYear),
            med_age = median(ageYear))

end_date <- today -3

ggplot(data = tests_by_week, mapping = aes(x = wk_comm, y = tests), color = "blue")+
  geom_bar(stat = "identity", fill = "grey", color = "grey")+
  geom_bar(data = tests_by_week, aes(x = wk_comm, y = positive),
           stat = "identity", fill = "blue", color = "grey")+
  labs(x = "Month", y = "serious covid event",title = "serious covid events over time")+
  scale_x_date(expand = c(0,0), limits = as.Date(c(pandemic_start, end_date)), date_breaks = "months", date_labels = "%b-%y")+
  #geom_vline(xintercept = as.Date("2021-09-14"),linetype="longdash", color = "red", size = 1)+
  geom_vline(xintercept = as.Date("2021-05-19"), linetype="longdash", color = "red", size = 1)+
  geom_vline(xintercept = as.Date("2021-12-15"), linetype="longdash",color = "red", size = 1)+
  geom_label(label = "Delta becomes\ndominant", x = as.Date("2021-04-10"), y = 550, color = "black", fill = "gray", hjust = 0,size=3)+
  #geom_label(label = "Booster dose\nrollout begins", x = as.Date("2021-07-10"), y = 550, color = "black", fill = "gray", hjust = 0,size=3)+
  geom_label(label = "Omicron becomes\ndominant", x = as.Date("2021-10-30"), y = 550, color = "black", fill = "gray", hjust = 0,size=3)
  #theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
  #theme_fivethirtyeight()
Covid_PLOT
# SICSAG_PLOT
ggsave(plot = Covid_PLOT, paste0(filename="/conf/EAVE/GPanalysis/progs/UA/COVID_test_", today, ".png"), width=17, height=7)

