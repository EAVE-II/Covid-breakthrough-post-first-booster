##### Code for Utkarsh #####


##### 0 - Set up #####
#Libraries
library("tidyverse")
library("survival")

## Load data
z_chrt_desc <- readRDS("/conf/EAVE/GPanalysis/analyses/COVID_vaccines_1st_dose/Waning/Covid-vaccine-waning/output/z_chrt_desc.rds")
z_chrt_desc <- filter(df_cohort_vacc1,event==1)
# Subset to vaccination data
z_vacc <- z_chrt_desc %>%
  filter(vacc1 == 1)

# Only those with an event
z_vacc_event <- z_vacc %>%
  mutate(time_to_event = as.numeric(admission_date-date_vacc_1)) %>%
  filter(event == 1)

# Only those with an event >= 14 days
z_vacc_event14 <- z_vacc_event %>%
  filter(time_to_event >= 14)

sum(z_vacc_event14$event==1)
sum(z_vacc_event14$event_hosp==1)
sum(z_vacc_event14$event_death==1)

table(z_vacc_event14$event_death, z_vacc_event14$event_hosp)

# Test#
z_vacc_event_death <- z_vacc_event %>%
  mutate(time_to_death = as.numeric(NRS.Date.Death-date_vacc_1)) %>%
  filter(event_death==1)

which(z_vacc_event_death$time_to_death<0)




#### Utkarsh's code ####
#df_cohort_vacc_g <- z_vacc
# # slecting events 14 days post vacc
df_cohort_vacc_g <- readRDS(paste0("/conf/EAVE/GPanalysis/progs/UA/first_dose_failure/data/df_cohort_coxph_09-4-2021.RDS"))

z1<-filter(df_cohort_vacc_g,event==1)
z2<-filter(z1,date_period=="date_preiod_0007"|date_period=="date_preiod_0714")
df_cohort_vacc_g_orig<-df_cohort_vacc_g
df_cohort_vacc_g <- df_cohort_vacc_g %>%
  filter(!EAVE_LINKNO %in% z2$EAVE_LINKNO)
df_cohort_vacc_g<-filter(df_cohort_vacc_g,date_period!="date_preiod_0007"&
                           date_period!="date_preiod_0714")
# df_cohort_vacc_g<-filter(df_cohort_vacc_g,event==1 & (date_period!="date_preiod_0007"&
#                                                         date_period!="date_preiod_0714"))
df_cohort_vacc <- df_cohort_vacc %>%
  filter(!EAVE_LINKNO %in% z2$EAVE_LINKNO)

length(which(df_cohort_vacc_g$event==1))
length(which(df_cohort_vacc_g$event==1&(df_cohort_vacc_g$date_period=="date_preiod_0007"|
                                          df_cohort_vacc_g$date_period=="date_preiod_0714")))

length(which(df_cohort_vacc_g_orig$event==1))
length(which(df_cohort_vacc_g_orig$event==1&(df_cohort_vacc_g_orig$date_period=="date_preiod_0007"|
                                               df_cohort_vacc_g_orig$date_period=="date_preiod_0714")))

df_cohort_vacc_g$EAVE_Smoke_n[df_cohort_vacc_g$EAVE_Smoke=="Non Smoker"]<-0
df_cohort_vacc_g$EAVE_Smoke_n[df_cohort_vacc_g$EAVE_Smoke=="Smoker"]<-1
df_cohort_vacc_g$EAVE_Smoke_n[df_cohort_vacc_g$EAVE_Smoke=="Ex Smoker"]<-2
df_cohort_vacc_g$EAVE_Smoke_n[df_cohort_vacc_g$EAVE_Smoke=="Unknown"]<-3

# adding week of event as confounder
z<-data.frame(seq(as.Date("2020-12-08"),as.Date("2021-04-18"),"day"))
colnames(z)[1]<-"date_vacc_1"
z<-mutate(z,week_date=cut(z$date_vacc_1,"weeks"))
z$week_date<-as.Date(z$week_date)
z<-mutate(z,week_date=week_date+1)
z$week_date<-lag(z$week_date,1)
z$week_date[1]<-z$week_date[2]
z<-mutate(z,num_weeks_post_vacc=cumsum(!duplicated(week_date)))
z$week_date<-NULL
df_cohort_vacc_g<-left_join(df_cohort_vacc_g,z,by="date_vacc_1")

# BMI cat
df_cohort_vacc_g$bmi_impute<-round(df_cohort_vacc_g$bmi_impute,1)
df_cohort_vacc_g<-mutate(df_cohort_vacc_g,bmi_cat =
                           cut(ageYear, breaks = c(-Inf, 18.5,24.9,29.9,Inf),
                               labels=c("1","0","2","3")))
#underweight - 1
#healthy - 0
#overweight - 2
#obese - 3

# collapse urban rural for ARR analysis
df_cohort_vacc_g<-mutate(df_cohort_vacc_g,
                         ur_combined=if_else(ur6_2016_name=="1 Large Urban Areas"|
                                               ur6_2016_name=="2 Other Urban Areas"|
                                               ur6_2016_name=="3 Accessible Small Towns",1,2))
df_cohort_vacc_g$ur_combined <- as.factor(df_cohort_vacc_g$ur_combined)
df_cohort_vacc_g$ur_combined <- relevel(df_cohort_vacc_g$ur_combined, ref = 1)

# SIMD with numerical cat 
df_cohort_vacc_g$SIMD[df_cohort_vacc_g$simd2020_sc_quintile=="5-Low"]<-5
df_cohort_vacc_g$SIMD[df_cohort_vacc_g$simd2020_sc_quintile=="4"]<-4
df_cohort_vacc_g$SIMD[df_cohort_vacc_g$simd2020_sc_quintile=="3"]<-3
df_cohort_vacc_g$SIMD[df_cohort_vacc_g$simd2020_sc_quintile=="2"]<-2
df_cohort_vacc_g$SIMD[df_cohort_vacc_g$simd2020_sc_quintile=="1 - High"]<-1
df_cohort_vacc_g$SIMD <- as.factor(df_cohort_vacc_g$SIMD)
df_cohort_vacc_g$SIMD <- relevel(df_cohort_vacc_g$SIMD, ref = 5)


# Save
saveRDS(df_cohort_vacc_g, "./df_cohort_vacc_g.rds")


##### 1 - ICD-10 codes associated with hospital admission #####

# Take hosps
z_vacc_event14_hosp <- z_vacc_event14 %>%
  filter(event_hosp ==1)

# Link in SMR data
z_vacc_event14_hosp <- z_vacc_event14_hosp %>%
  left_join(for_UA_joined_smr_cis_flagged, by=c("EAVE_LINKNO_uv"="EAVE_LINKNO")) 

# How many got linked
tbl <- table(z_vacc_event14_hosp$exact_match) # 723 has SMR record

z_vacc_event14_hosp <- z_vacc_event14_hosp%>%
  filter(exact_match==1)

# How many had U071 has the main condition
n <- sum(z_vacc_event14_hosp$covid_acute_main, na.rm=T)
n/nrow(z_vacc_event14_hosp)*100

# Summarise all codes:
tbl2 <- z_vacc_event14_hosp %>%
  summarise_at(vars(covid_acute_main, covid_susp_main, covid_acute, covid_acute_susp),
               sum)
tbl2/nrow(z_vacc_event14_hosp)*100


##### 2 - Plots #####
vacc_type_label <- c("BNT162b2", "ChAdOx1")

### 2.1 Events across time for entire cohort ####
## COVID hospitalisations
p1 <- ggplot(z_chrt_desc) +
  # geom_density(aes(x=hosp_admission_date, fill=age_grp), alpha = 0.5) +
  geom_histogram(aes(x=event_date), binwidth=1, fill="blue") +
  theme_light() +
  #scale_fill_manual(values = c(eave_blue, eave_green, eave_orange)) +
  #scale_color_manual(values = c(eave_blue, eave_green, eave_orange)) +
  labs(x="Date of event", y="Number of events",
       fill="Age group", subtitle = "COVID-19 serious event") +
  geom_vline(xintercept = as.Date("2021-05-19"), linetype="longdash", color = "red", size = 1) +
  annotate("text", x=as.Date("2021-05-21"), y=100, label = "Delta variant \nbecame dominant", hjust=0, size=3.5) +
  geom_vline(xintercept = as.Date("2021-09-14"), linetype="longdash", color = "red", size = 1) +
  annotate("text", x=as.Date("2021-09-14"), y=100, label = "Booster dose \nrollout begins", hjust=-0.05, size=3.5) +
  geom_vline(xintercept = as.Date("2021-12-15"), linetype="longdash", color = "red", size = 1) +
  annotate("text", x=as.Date("2021-11-20"), y=100, label = "Omicron variant \nbecame dominant", hjust=0, size=3.5)


p1

## COVID mortality
p2 <- ggplot(z_chrt_desc) +
  # geom_density(aes(x=hosp_admission_date, fill=age_grp), alpha = 0.5) +
  #geom_histogram(aes(x=NRS.Date.Death, fill=age_grp), position="dodge", binwidth=5) +
  geom_histogram(aes(x=NRS.Date.Death), binwidth=1, fill=eave_blue) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green, eave_orange)) +
  scale_color_manual(values = c(eave_blue, eave_green, eave_orange)) +
  labs(x="Date of COVID-19 death", fill="Age group",
       subtitle = "COVID-19 Death") +
  geom_vline(xintercept = as.Date("2020-12-08"), linetype = 2) +
  #annotate("text", x=as.Date("2020-12-10"), y=220, label = "BNT162b2 introduced", hjust=0, size=3.5) +
  geom_vline(xintercept = as.Date("2021-01-04"), linetype = 2) #+
#annotate("text", x=as.Date("2021-01-06"), y=250, label = "ChAdOx1 introduced", hjust=0, size=3.5)



## COVID hospitalisations and deaths (Composite outcome)
p3 <- ggplot(z_chrt_desc) +
  #geom_histogram(aes(x=admission_date, fill=age_grp), position="dodge", binwidth=5) +
  geom_histogram(aes(x=admission_date), binwidth=1, fill=eave_blue) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green, eave_orange)) +
  scale_color_manual(values = c(eave_blue, eave_green, eave_orange)) +
  labs(x="Date of COVID-19 hospitalisation or death", fill="Age group",
       subtitle = "Composite outcome: COVID-19 hospitalisation or death")+
  geom_vline(xintercept = as.Date("2020-12-08"), linetype = 2) +
  #annotate("text", x=as.Date("2020-12-10"), y=220, label = "BNT162b2 introduced", hjust=0, size=3.5) +
  geom_vline(xintercept = as.Date("2021-01-04"), linetype = 2) #+
#annotate("text", x=as.Date("2021-01-06"), y=250, label = "ChAdOx1 introduced", hjust=0, size=3.5)




#png(file=paste0("./output/final/descriptives/events_age.png"),
#   width = 1000, height=300)
gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
#dev.off()




##### 2.2 - Uptake of vaccinations over time ####
# Labels for vacc type
vacc_type_label <- c("BNT162b2", "ChAdOx1")
names(vacc_type_label) <- c("PB", "AZ")

## Count over time by vaccine type
z_chrt_desc %>%
  ggplot() +
  geom_histogram(aes(x=date_vacc_1, fill=vacc_type), position="dodge", binwidth=2) +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  scale_color_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  labs(x="Date of 1st vaccine (2020/21)", fill="Vaccine type",
       subtitle = paste0("All 1st dose vaccine dates")) +
  scale_y_continuous(labels = function(x) format(x, scientific = F))


## Cumulative plot
z_chrt_desc %>%
  ggplot(aes(x=date_vacc_1, col=vacc_type)) +
  stat_bin(aes(y=cumsum(..count..)), geom="step", binwidth = 1) + 
  #geom_step(stat="count") +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  scale_color_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  labs(x="Date of 1st vaccine (2020/21)", colour="Vaccine type") +
  scale_y_continuous(labels = function(x) format(x, scientific = F))


# Calculate the length?


# Plot
z_chrt_desc %>%
  ggplot(aes(x=date_vacc_1, col=vacc_type)) +
  #stat_bin(aes(y=cumsum(..count..)), geom="step", binwidth = 1) + 
  geom_step(stat="ecdf") +
  theme_light() +
  scale_fill_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  scale_color_manual(values = c(eave_blue, eave_green), labels = vacc_type_label) +
  labs(x="Date of 1st vaccine (2020/21)", colour="Vaccine type") +
  scale_y_continuous(labels = function(x) format(x, scientific = F))




##### 3 - Carry out model #####

#### 3.1 - Check data ####
## Rachel's cohort
z_vacc_cohort <- z_vacc %>%
  mutate(time_to_event = as.numeric(admission_date-date_vacc_1)) %>%
  filter(time_to_event >= 14 | is.na(time_to_event))

# Slightly different numbers - df_cohort_vacc_g

## Check with UAs
length(unique(df_cohort_vacc_g$EAVE_LINKNO))



### Censor events that had hosp with U701 as cause
for_UA_joined_smr_cis_flagged <- readRDS("/conf/EAVE/GPanalysis/data/temp/for_UA_joined_smr_cis_flagged.rds")
# Link in
df_cohort_vacc_g <- df_cohort_vacc_g %>%
  left_join(select(for_UA_joined_smr_cis_flagged,
                   EAVE_LINKNO, covid_acute), by=c("EAVE_LINKNO"="EAVE_LINKNO"))

# Check numbers
sum(df_cohort_vacc_g$event[which(df_cohort_vacc_g$hosp_covid==1)])
length(which(df_cohort_vacc_g$event==1 & df_cohort_vacc_g$hosp_covid==1))
# Not sure what's going on here

table(df_cohort_vacc_g$event, df_cohort_vacc_g$covid_acute)


## Censor
df_cohort_vacc_g2 <- df_cohort_vacc_g %>%
  mutate(event = ifelse(event==1 & covid_acute!=1 & event_hosp == 1, 0, event))

table(df_cohort_vacc_g2$event, df_cohort_vacc_g2$covid_acute)


library("finalfit")

explanatory <- c("date_period", "Sex", "age_gp" , 
                 "prev_positive_status", "in_hosp_status", "care_home_elderly",
                 "SIMD" , "ur_combined", "EAVE_Smoke", "n_risk_gps",  "n_tests_gp", "vacc_type",
                 "HB", "bmi_cat", "num_weeks_post_vacc")

tbl <- df_cohort_vacc_g2 %>%
  mutate(event = as.character(event)) %>%
  summary_factorlist("event", explanatory, p = F)

write.csv(tbl, "./sensitivity_summary_tbl.csv")

## Function
glm_pyears_fn <- function(variable){
  z.fmla <- as.formula(paste("Surv(start_time,end_time,event)",
                             " ~ ",
                             paste(variable, collapse= "+")))
  
  p_years <- pyears(z.fmla, data = df_cohort_vacc_g,
                    data.frame = TRUE)
  
  glm_pos <- glm(event~offset(log(pyears)) + get(variable), family=poisson, data=p_years$data)
  
  round(exp(confint(glm_pos)),2)
}



##### 3.2 - Multivariate model(s) ####


### Person years
p_years<-pyears(Surv(start_time,end_time,event)~date_period + Sex + age_gp  + 
                  prev_positive_status + in_hosp_status + care_home_elderly +
                  SIMD  + ur_combined + EAVE_Smoke + n_risk_gps +  n_tests_gp +vacc_type +
                  HB + bmi_cat + num_weeks_post_vacc, 
                data = df_cohort_vacc_g2, data.frame = TRUE)


### Overall
# Model for all
glm_pos <- glm(event~offset(log(pyears)) +date_period + Sex + age_gp  + 
                 prev_positive_status + in_hosp_status + as.character(care_home_elderly) +
                 SIMD  + ur_combined + EAVE_Smoke + n_risk_gps +  n_tests_gp +
                 HB + bmi_cat + num_weeks_post_vacc, 
               family=poisson,data=p_years$data)


z_glm_est <- as.data.frame(round(exp(cbind(glm_pos$coefficients, confint.default(glm_pos) ) ), 3)) %>%
  rownames_to_column(var ="Var") %>%
  rename(est=2, lwr = 3, upr = 4) %>%
  mutate(est = sprintf("%.2f",round(est,2)),
         upr = sprintf("%.2f",round(upr,2)),
         lwr = sprintf("%.2f",round(lwr,2))) %>%
  mutate(RR_est = paste0(est, " (", lwr, ", ", upr, ")")) %>%  
  select(-c(est, lwr, upr))

z_glm_est

write.csv(z_glm_est, "./sensitivity_both.csv")


### By vacc

glm_multi_vacc <- function(z_vacc_type){
  glm_pos <- glm(event~offset(log(pyears)) +date_period + Sex + age_gp  + 
                   prev_positive_status + in_hosp_status + as.character(care_home_elderly) +
                   SIMD  + ur_combined + EAVE_Smoke + n_risk_gps +  n_tests_gp +
                   HB + bmi_cat + num_weeks_post_vacc, 
                 family=poisson,data=p_years$data,
                 subset = vacc_type == z_vacc_type)
  
  
  z_glm_est <- as.data.frame(round(exp(cbind(glm_pos$coefficients, confint.default(glm_pos) ) ), 3)) %>%
    rownames_to_column(var ="Var") %>%
    rename(est=2, lwr = 3, upr = 4) %>%
    mutate(est = sprintf("%.2f",round(est,2)),
           upr = sprintf("%.2f",round(upr,2)),
           lwr = sprintf("%.2f",round(lwr,2))) %>%
    mutate(RR_est = paste0(est, " (", lwr, ", ", upr, ")")) %>%  
    select(-c(est, lwr, upr))
  
  z_glm_est
}

z_pb <- glm_multi_vacc("PB")
z_az <- glm_multi_vacc("AZ")


# Join
z_glm_est <- z_glm_est %>%
  rename(RR_est_both = RR_est) %>%
left_join(left_join(z_pb,z_az))


write.csv(z_pb, "./sensitivity_pb.csv")
write.csv(z_az, "./sensitivity_az.csv")



##### Utkarsh's adjusted model #####


p_years<-pyears(Surv(start_time,end_time,event)~date_period + Sex + age_gp  + 
                  SIMD  + ur_combined + EAVE_Smoke + bmi_cat + HB + prev_positive_status + 
                  n_risk_gps + in_hosp_status + n_tests_gp + care_home_elderly+num_weeks_post_vacc, 
                data = df_cohort_vacc_g, data.frame = TRUE)
glm_pos <- glm(event~offset(log(pyears)) +date_period+ Sex + age_gp + 
                 SIMD  + ur_combined + EAVE_Smoke + HB + prev_positive_status +
                 n_risk_gps  + in_hosp_status + n_tests_gp + as.character(care_home_elderly)+num_weeks_post_vacc, 
               family=poisson,data=p_years$data)

summary(glm_pos)
#exp(glm_pos$coefficients)
round(exp(glm_pos$coefficients),2)
#round(exp(confint(glm_pos)),2)
round(exp(glm_pos$coefficients-coef(summary(glm_pos))[,2]*1.96),2)
round(exp(glm_pos$coefficients+coef(summary(glm_pos))[,2]*1.96),2)




###### Checking deaths #####

# Checking no. deaths
table(df_cohort_vacc_g$event, df_cohort_vacc_g$event_death)

# All deaths
all_deaths <- readRDS("/conf/EAVE/GPanalysis/data/all_deaths.rds")

death_causes <- unique(sort(all_deaths$CAUSE_OF_DEATH_CODE_0))
death_causes[which(startsWith(death_causes, "U0"))]


all_deaths <- readRDS("/conf/EAVE/GPanalysis/data/all_deaths.rds") %>%
  mutate(main_cod = ifelse(CAUSE_OF_DEATH_CODE_0 == "U071", "U071",
                           ifelse(CAUSE_OF_DEATH_CODE_0=="U072", "U072", "Other")))# %>%
  #mutate(u071_cod = if_any(all_vars() %in% "U071", "U071", "Non-U071")) %>%
  #select(EAVE_LINKNO, main_cod,u071_cod)


all_deaths2 <- all_deaths %>%
  #mutate(u071_cod = if_else(rowSums(all_deaths == "U071") > 0, 1,0))
 # mutate(u071_cod = ifelse(apply(all_deaths == "U071", 1, any), 1, 0))
  #mutate(u071_cod = case_when(any_vars() == "U071"~ 1,TRUE ~ 0))
  mutate(u071_cod = if_else(CAUSE_OF_DEATH_CODE_0 == "U071" | CAUSE_OF_DEATH_CODE_1 == "U071" |
                              CAUSE_OF_DEATH_CODE_2 == "U071" | CAUSE_OF_DEATH_CODE_3 == "U071" |
                              CAUSE_OF_DEATH_CODE_4 == "U071" | CAUSE_OF_DEATH_CODE_5 == "U071" |
                              CAUSE_OF_DEATH_CODE_6 == "U071" | CAUSE_OF_DEATH_CODE_7 == "U071" |
                              CAUSE_OF_DEATH_CODE_8 == "U071" | CAUSE_OF_DEATH_CODE_9 == "U071" , 1,0)) %>%
  mutate(u072_cod = if_else(CAUSE_OF_DEATH_CODE_0 == "U072" | CAUSE_OF_DEATH_CODE_1 == "U072" |
                              CAUSE_OF_DEATH_CODE_2 == "U072" | CAUSE_OF_DEATH_CODE_3 == "U072" |
                              CAUSE_OF_DEATH_CODE_4 == "U072" | CAUSE_OF_DEATH_CODE_5 == "U072" |
                              CAUSE_OF_DEATH_CODE_6 == "U072" | CAUSE_OF_DEATH_CODE_7 == "U072" |
                              CAUSE_OF_DEATH_CODE_8 == "U072" | CAUSE_OF_DEATH_CODE_9 == "U072" , 1,0)) %>%
  select(EAVE_LINKNO, main_cod, u071_cod, u072_cod)



table(all_deaths2$main_cod, all_deaths2$u071_cod)


df_cohort_vacc_g <- df_cohort_vacc_g %>%
  select(-CAUSE_OF_DEATH_CODE_0) %>%
  left_join(all_deaths2) %>%
  mutate(main_cod = replace_na(main_cod, "None")) %>%
  mutate(u071_cod = replace_na(u071_cod, 0)) %>%
  mutate(u072_cod = replace_na(u072_cod, 0))
  

table(df_cohort_vacc_g$event, df_cohort_vacc_g$main_cod)
table(df_cohort_vacc_g$event, df_cohort_vacc_g$u071_cod)
table(df_cohort_vacc_g$event, df_cohort_vacc_g$u072_cod)
