##########################################################
# Name of file: 02a_Get_Full_2Dose_Cohort_Long.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@phs.scot
# Original date: 14 Feb 2022
# Latest update author (if not using version control) - Chris Robertson chris.robertson@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: gets 2 dose cohort in long format
#                         run after 02a_Get_Full_Cohort.R
# Approximate run time: Unknown
#
#   run 01a_Data_Omicron###############################

df_all <- readRDS("/conf/EAVE/GPanalysis/progs/CR/Vaccine_Dose_2_Booster/output/temp/xdf_all_full_u.RDS")
Location <- "/conf/"  # Server
all_deaths  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/all_deaths.rds"))
summary(all_deaths)
z <- all_deaths %>%  
  mutate(across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9, ~if_else(. %in% c("U071","U072"), 1,0)))
z <- z %>% rowwise() %>% mutate(rowsum = sum(c_across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9))) %>% 
  mutate(covid_death_cert = if_else(rowsum>=1,1,0)) 
all_deaths <- z %>% dplyr::select(-rowsum) %>% ungroup()
summary(dplyr::select(all_deaths, NRS.Date.Death, NRS.Reg.Date, covid_death_cert))
#df_all is the full cohort from 02a_Get_Full_Cohort.R
#df_all <- readRDS(paste0("./output/temp/df_all_full.RDS"))
remove(list=ls(pa="^z"))
output_list <- list()

a_begin <- as.Date("2021-09-13")  #beginning of booster period
output_list$a_begin <- a_begin
output_list$endpoint <- "covid_hosp_death"  # "covid_death", "covid_hosp", "covid_hosp_death"

#filter to get the 2 dose cohort only
z_df_all_2d <- df_all %>% filter(!is.na(date_vacc_2))  # only keep those who have had 2 doses
z_df_all_2d  <- z_df_all_2d  %>% filter(ageYear >= 18)  # only keep adults
z_df_all_2d <- z_df_all_2d %>% filter(flag_incon==0)  # only keep consistent vacccine records

z_df_all_2d<-mutate(z_df_all_2d,age_gp = 
                 cut(ageYear, breaks = c(-Inf, 49,54,59,64,69,74,79,Inf),
                     labels=c("18-49","50-54","55-59","60-64","65-69",
                              "70-74","75-79","80+")))

z_df_all_2d<-mutate(z_df_all_2d,ur_combined=
                      if_else(ur6_2016_name=="1 Large Urban Areas"|
                                ur6_2016_name=="2 Other Urban Areas"|
                                ur6_2016_name=="3 Accessible Small Towns",1,2))

z1 <- z_df_all_2d %>% 
  mutate(days_from_previous_pos = as.numeric(a_begin - date_ecoss_specimen_prior))  # all previous tests before a_begin
z1 <- z1 %>% mutate(pos_before_start_u = case_when(is.na(days_from_previous_pos) ~ "not_prev_pos",
                                                 days_from_previous_pos >= 1 & days_from_previous_pos <= 89 ~ "pos_0:2_months",
                                                 days_from_previous_pos >= 90 & days_from_previous_pos <= 179 ~ "pos_3:5_months",
                                                 days_from_previous_pos >= 180 & days_from_previous_pos <= 269 ~ "pos_6:8_months",
                                                 TRUE ~ "pos_9+_months")) 
z_df_all_2d<-z1

z1 <- z_df_all_2d %>%  
  mutate(days_between_vacc = as.numeric(date_vacc_2-date_vacc_1))%>%
  mutate(vacc_gap = cut(days_between_vacc, breaks=c(1,48, 62, 76, 90, max(days_between_vacc, na.rm=T)),
                        labels=c("<7wk","7-8 wk","9-10 wk","11-12 wk", "13+ wk")))
z_df_all_2d<-z1
rm(z1)
#z_df_all_2d1 <-filter(z_df_all_2d,)

if (output_list$endpoint == "covid_death") z_df_all_2d <- z_df_all_2d %>%
  mutate(event=covid_death, event_date=covid_death_date)
if (output_list$endpoint == "covid_hosp") z_df_all_2d <- z_df_all_2d %>%
  mutate(event=covid_hosp, event_date=covid_hosp_date)
if (output_list$endpoint == "covid_hosp_death") z_df_all_2d <- z_df_all_2d %>%
  mutate(event=covid_hosp_death, event_date=covid_hosp_death_date)

output_list$a_end <- max(z_df_all_2d$event_date)

#sample cohort - z_n_controls_per_event non events per event
event_ids <- filter(z_df_all_2d, event==1 & event_date > a_begin)%>% pull(EAVE_LINKNO) %>% unique()
output_list$seed <- 21021954
set.seed(output_list$seed)
output_list$controls_per_event <- 10 
if (output_list$endpoint == "covid_pos") output_list$controls_per_event <- 5 
if (output_list$endpoint == "covid_death") output_list$controls_per_event <- 50 

no_event_ids <- filter(z_df_all_2d, event==0 & event_date > a_begin) %>% 
  filter(if (output_list$endpoint == "covid_pos") pos_before_start != "pos_1:28" else TRUE) %>% 
  slice_sample(n=output_list$controls_per_event*length(event_ids)) %>%
  pull(EAVE_LINKNO) %>% unique()
z_df <- z_df_all_2d %>% dplyr::select(EAVE_LINKNO, date_vacc_1, date_vacc_2, date_vacc_3, event_date, event) %>% #keep minimal columns
  filter(EAVE_LINKNO %in% c(event_ids, no_event_ids))
print(sum(z_df$event))

#z_df <- slice_sample(z_df, n=10000)

#add in deaths and change event_date
#have checked no event_dates after date death
#this bit just checks and does no changes to the z_df data frame
z <- z_df %>% left_join(all_deaths, by="EAVE_LINKNO") %>% 
  mutate(event_date = if_else(!is.na(NRS.Date.Death) & (NRS.Date.Death < event_date), NRS.Date.Death, event_date))
print(sum(z$event))
print(summary(z$event_date))
print(summary(z_df$event_date))
print(summary(z))
table(z$event_date > z$NRS.Date.Death, exclude=NULL)


#set up the survival times
z <- tmerge(z_df,z_df,id=EAVE_LINKNO, endpt = event(event_date, event),tstart=a_begin, tstop=event_date)
z <- z %>% mutate(across(
  .cols = c(tstart, tstop),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
#get the post vacc 1 periods
z_df <- z_df %>% mutate(period_end_date = date_vacc_1)
#z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
z_df <- z_df %>% mutate(across(
  .cols = c(period_end_date),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
z <- tmerge(z,z_df, id=EAVE_LINKNO,
            per1=tdc(period_end_date)
)
names(z)[names(z)=="per1"] <- "pv_uv"

#get the post vacc 1 periods
z_df <- z_df %>% mutate(period_end_date = date_vacc_1 + 27)
#z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
z_df <- z_df %>% mutate(across(
  .cols = c(period_end_date),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
z <- tmerge(z,z_df, id=EAVE_LINKNO,
            per1=tdc(period_end_date)
)
names(z)[names(z)=="per1"] <- "pv_v1_1"

z_df <- z_df %>% mutate(period_end_date = date_vacc_2)
#z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
z_df <- z_df %>% mutate(across(
  .cols = c(period_end_date),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
z <- tmerge(z,z_df, id=EAVE_LINKNO,
            per1=tdc(period_end_date)
)
names(z)[names(z)=="per1"] <- "pv_v1_2"

#get the post vaccination periods - they do not have to be the same length
z_df <- z_df %>% mutate(period_end_date = date_vacc_2 + 13)
#z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
z_df <- z_df %>% mutate(across(
  .cols = c(period_end_date),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
z <- tmerge(z,z_df, id=EAVE_LINKNO,
            per1=tdc(period_end_date)
)
names(z)[names(z)=="per1"] <- "pv_v2_1"

z_df <- z_df %>% mutate(period_end_date = date_vacc_2 + 69)
#z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
z_df <- z_df %>% mutate(across(
  .cols = c(period_end_date),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
z <- tmerge(z,z_df, id=EAVE_LINKNO,
            per1=tdc(period_end_date)
)
names(z)[names(z)=="per1"] <- "pv_v2_2"

z_df <- z_df %>% mutate(period_end_date = date_vacc_2 + 139)
#z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
z_df <- z_df %>% mutate(across(
  .cols = c(period_end_date),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
z <- tmerge(z,z_df, id=EAVE_LINKNO,
            per1=tdc(period_end_date)
)
names(z)[names(z)=="per1"] <- "pv_v2_3"

z_df <- z_df %>% mutate(period_end_date = date_vacc_3)
#z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
z_df <- z_df %>% mutate(across(
  .cols = c(period_end_date),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
z <- tmerge(z,z_df, id=EAVE_LINKNO,
            per1=tdc(period_end_date)
)
names(z)[names(z)=="per1"] <- "pv_v2_4"

#get the post vacc 3 periods
z_df <- z_df %>% mutate(period_end_date = date_vacc_3 + 13)
#z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
z_df <- z_df %>% mutate(across(
  .cols = c(period_end_date),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
z <- tmerge(z,z_df, id=EAVE_LINKNO,
            per1=tdc(period_end_date)
)
names(z)[names(z)=="per1"] <- "pv_v3_1"

z_df <- z_df %>% mutate(period_end_date = date_vacc_3 + 34)
#z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
z_df <- z_df %>% mutate(across(
  .cols = c(period_end_date),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
z <- tmerge(z,z_df, id=EAVE_LINKNO,
            per1=tdc(period_end_date)
)
names(z)[names(z)=="per1"] <- "pv_v3_2"

z_df <- z_df %>% mutate(period_end_date = date_vacc_3 + 55)
#z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
z_df <- z_df %>% mutate(across(
  .cols = c(period_end_date),
  .fns = ~ interval(a_begin, .x) / ddays(1)
))
z <- tmerge(z,z_df, id=EAVE_LINKNO,
            per1=tdc(period_end_date)
)
names(z)[names(z)=="per1"] <- "pv_v3_3"

# z_df <- z_df %>% mutate(period_end_date = date_vacc_3 + 90)
# #z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
# z_df <- z_df %>% mutate(across(
#   .cols = c(period_end_date),
#   .fns = ~ interval(a_begin, .x) / ddays(1)
# ))
# z <- tmerge(z,z_df, id=EAVE_LINKNO,
#             per1=tdc(period_end_date)
# )
# names(z)[names(z)=="per1"] <- "pv_v3_4"

#get the temporal periods
z_period_length <- 7
z_n_periods <- trunc((as.numeric(output_list$a_end-a_begin)/z_period_length))
for (i in 1:z_n_periods){
  #i <- 1
  print(i)
  z_df <- z_df %>% mutate(period_end_date = a_begin + z_period_length*i)
  #z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
  z_df <- z_df %>% mutate(across(
    .cols = c(period_end_date),
    .fns = ~ interval(a_begin, .x) / ddays(1)
  ))
  z <- tmerge(z,z_df, id=EAVE_LINKNO,
              per1=tdc(period_end_date)
  )
  names(z)[names(z)=="per1"] <- paste0("p_",i)
}

z_names <- names(z)[grepl("p_", names(z))]
z1 <- z %>% mutate(period = apply(z[,z_names], 1, sum))
z1 <- z1 %>% dplyr::select(-all_of(z_names))
z_names <- names(z)[grepl("pv_", names(z))]
#z_names <- z_names[!grepl("pv_v3",z_names)] # omit dose 3
z1 <- z1 %>% mutate(pv_period = apply(z[,z_names], 1, sum)) 
z1 <- z1 %>% dplyr::select(-all_of(z_names))
#z_names <- names(z)[grepl("pv_v3", names(z))]
#z1 <- z1 %>% mutate(pv_v3_period = apply(z[,z_names], 1, sum)) 
#z1 <- z1 %>% dplyr::select(-all_of(z_names))
z <- z1

#having got the time dependent periods move the data back to z_df
#z_df in long format
#endpt is the response variable; the origianl one event is duplicated when the periods and pv_periods are added
#don't use event.
#tstart and tstop are the start and stop dates of the intervals
z_df <- z
print(sum(z_df$endpt))


#z_df <- z_df %>% mutate(tstop=as.Date(tstop, origin=as.Date("1970-01-01")))
#calculate the interval lengths 
z_df <- z_df %>% mutate(pyears = as.numeric(tstop-tstart))

#get the dates for the beginning of the periods
z <- z_df %>% group_by(period) %>% 
  dplyr::summarise(date_period_begin = min(tstart), n_events =sum(endpt), pyears=sum(pyears))

#merge in the covariates and calculate other covariates
#if (output_list$endpoint == "covid_hosp") { 
#  df <- z_df %>% 
#    left_join(dplyr::select(z_df_all_2d, EAVE_LINKNO, vacc_type, vacc_type_3, flag_incon, NRS.Date.Death, emergency, covid_mcoa_hosp,
#                            Sex:age_gp, vs:immuno, n_oth_risk_gps), by="EAVE_LINKNO") %>%
#    filter(!is.na(ageYear)) # %>%  #drop those who do not link into EAVE
#} 
if (output_list$endpoint %in% c("covid_hosp", "covid_death", "covid_hosp_death")) { 
  df <- z_df %>% 
    left_join(dplyr::select(z_df_all_2d, EAVE_LINKNO, vacc_type, vacc_type_3, 
                            flag_incon, NRS.Date.Death, Sex:age_gp, vs:immuno, 
                            n_oth_risk_gps,ur_combined,pos_before_start_u), by="EAVE_LINKNO") %>%
    filter(!is.na(ageYear)) # %>%  #drop those who do not link into EAVE
}

z_levs <- sort(unique(df$pv_period))
z_labs <- c("uv", "v1_0:3","v1_4+", "v2_0:1" ,  "v2_2:9" ,  "v2_10:19", "v2_20+", 
            "v3_0:1"  , "v3_2:4" ,  "v3_5:7", "v3_8+")
# z_labs <- c("uv", "v1_0:3","v1_4+", "v2_0:1" ,  "v2_2:9" ,  "v2_10:19", "v2_20+", 
#             "v3_0:1"  , "v3_2:4" ,  "v3_5:8", "v3_9:12", "v3_13+")
#make period and pv_period factors
z_df <- df %>%
  mutate(period_f = factor(period, levels=0:max(period), labels=z$date_period_begin, ordered=FALSE),
         pv_period_f = factor(pv_period, levels=z_levs,labels=z_labs))

df <- z_df  #df is the working data set for the analysis

#calculate the weights for the sampling - event is correct here
df <- df %>% mutate(weight = if_else(event==1,1,nrow(z_df_all_2d)/(output_list$controls_per_event*length(event_ids)) ) )  
#merge the sampling weights and eave_weights
df <- df %>% mutate(ew = eave_weight*weight)



#get dose 3 vaccine type added to pv_period_f
df <- df %>% mutate(pv_period_f_d3b = case_when(is.na(vacc_type_3) ~ as.character(pv_period_f),
                                                !is.na(vacc_type_3) & !grepl("v3", pv_period_f) ~ as.character(pv_period_f),
                                                !is.na(vacc_type_3) & grepl("v3", pv_period_f) ~ paste(as.character(pv_period_f),vacc_type_3,sep="_"),
                                                TRUE ~ "unknown"))
# z_levs <- c(levels(df$pv_period_f)[!grepl("v3", levels(df$pv_period_f))],
#             "v3_0:1_PB", "v3_2:4_PB", "v3_5:8_PB", "v3_9:12_PB", "v3_13+_PB",
#             "v3_0:1_Mo", "v3_2:4_Mo", "v3_5:8_Mo", "v3_9:12_Mo", "v3_13+_Mo",
#             "v3_0:1_AZ", "v3_2:4_AZ", "v3_5:8_AZ", "v3_9:12_AZ", "v3_13+_AZ")
z_levs <- c(levels(df$pv_period_f)[!grepl("v3", levels(df$pv_period_f))],
            "v3_0:1_PB", "v3_2:4_PB", "v3_5:7_PB", "v3_8+_PB",
            "v3_0:1_Mo", "v3_2:4_Mo", "v3_5:7_Mo", "v3_8+_Mo",
            "v3_0:1_AZ", "v3_2:4_AZ", "v3_5:7_AZ", "v3_8+_AZ")
df <- df %>% mutate(pv_period_f_d3b = factor(pv_period_f_d3b, levels=z_levs))

#endpt is the response, pyears is the offset
print(output_list$endpoint)
saveRDS(df, paste0("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/data/xdf_full_",
                   output_list$endpoint,".RDS"))
#saveRDS(z_df_all_2d, paste0("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/data/xdf_full_cohort.RDS"))
remove(list=ls(pa="^z"))
