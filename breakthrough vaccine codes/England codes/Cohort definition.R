# Define the study and sample design 
rm(list = ls())
library(data.table)
library(lubridate)
setwd("/home/cig-RS-jasono/PDrive/DaCVaP/")
source("src/ONS infection survey estimates.R")

chrt <- readRDS(file = "data/FullCohrt.rds")
setDT(chrt)
nrow(chrt)

# first day as per protocol
(first.day <- as.Date("2020-12-07"))

## last day (study close)
(last.day <- as.Date("2022-03-31"))

DT <- chrt[!{Age < 18 | Age > 110} & 
              !{DeRegistrationDateNoNulls < as.Date("2020-12-07")} & 
              !is.na(SecondVaccinationDoseDate) & 
              SecondVaccinationBrand %in% c("AstraZeneca","Pfizer-BioNTech","Moderna") & 
              !{ThirdVaccinationBrand %in% c("Janssen","Novavax","Valneva")} 
           ]
nrow(DT)

## dropping the HESHospitalisation2021 people as they may not have hospitalisation data
DT <- DT[CheckedForHESHospitalisation2021==1]

nrow(DT)

# Easier to define here as opposed to the last filter. 
idx <- DT[DateOfDeath <= (SecondVaccinationDoseDate + 14),StudyID]
length(idx)
DT <- DT[!{StudyID %in% idx}]
nrow(DT)

oldn <- c("SecondVaccinationDoseDate","ThirdVaccinationDoseDate")
newn <- c("ScndVaccDate","ThrdVaccDate")
setnames(DT,old = oldn, new = newn)

################################################################################
#################### hospital records ##########################################
#########- confirmed/suspected hospital admission ##############################

AH <- setDT(readRDS(file = "data/AllHospA.rds"))
names(AH)
## clean the discharge data  - for now assume same as admissionDate
AH[DischargeDate == "1900-01-01", DischargeDate := AdmissionDate + 1]
AH[DischargeDate == "1800-01-01", DischargeDate := AdmissionDate + 1]
AH[DischargeDate < AdmissionDate, DischargeDate := AdmissionDate + 1]

AH[,LOShosp := as.numeric(difftime(DischargeDate,AdmissionDate, unit = "days"))]

# keep if they meet earlier criteria and Admission date is in study time frame
AH1 <- AH[StudyID %in% DT$StudyID &  AdmissionDate >= first.day]

# take just admission date and diagnosis column
AH1 <- AH1[,.(StudyID,AdmissionDate,DischargeDate,LOShosp,COVID19Diagnosis)]

# merge in second and third dose dates - use to place them in second dose or booster phase
tmp <- DT[,.(StudyID,ScndVaccDate,ThrdVaccDate)]
AH2 <- merge(AH1,tmp,by = "StudyID", all.x = T)

# If admission is after Scnd dose run-in phase and is COVID  - 
AH2b <- AH2[AdmissionDate > (ScndVaccDate + 14) & 
              COVID19Diagnosis %in% c("Confirmed","Suspected")]

AH2b[,AdmissionCount:=1:.N,by = StudyID]
AH2b[AdmissionCount==8,]
AH2b[StudyID==1853613]

# order by date, and use the first of these only 
setorder(AH2b,AdmissionDate,StudyID)
AH2b <- AH2b[,.SD[1],by = StudyID]

# if they have a third vaccination and admission is after run-in period call this booster phase
AH2b[!is.na(ThrdVaccDate), phase := ifelse(AdmissionDate > (ThrdVaccDate + 14),"booster","scnd dose")]

# if third dose is missing then it must be second dose phase
AH2b[is.na(ThrdVaccDate), phase := "scnd dose"]

# create label for these particular events
AH2b[,event:=1]

# tabulate these by date to see how the events occur in calendar time.
AH2b[,month:=month(AdmissionDate)]
res <- AH2b[,length(AdmissionDate),by = .(month,phase)]
setorder(res,month,phase)
sum(res$V1)

# 21,447 in the second dose phase, 7269 in the booster period.
DT <- merge(DT,AH2b,by = c("StudyID","ScndVaccDate","ThrdVaccDate"),all.x=T)

DT[,sum(event), by = phase]

#DT[is.na(phase) & event==1,StudyID][2]

################################################################################
##### Events defined by admission within 14,-1 days of a positive test #########
##### merge with hospital data #################################################

## Merge symptom and test data
Test <- setDT(readRDS(file = "data/PosT.rds"))
Test <- Test[StudyID %in% DT$StudyID & TestDate >= first.day]
Symp <- setDT(readRDS(file = "data/Symp.rds"))
Symp <- Symp[StudyID %in% DT$StudyID & TestDate >= first.day]
TS <- merge(Test,Symp,by = c("StudyID","TestDate"),all.x=T)
length(unique(TS$StudyID))

# N = 2,360,620 

# find IDS with test data and hospital admission data
uids <- TS$StudyID[TS$StudyID %in% AH2$StudyID]

# 804,176 have test data and appear in the hospital data set 
length(uids)

# Do the rolling merge
tmp1 <- TS[StudyID %in% uids,.(StudyID,TestResult,TestDate)]

# this allows test done 1 day after admission 
tmp1[,mtime := TestDate-1]
tmp2 <- AH2[StudyID %in% uids,.(StudyID,AdmissionDate,LOShosp,DischargeDate,COVID19Diagnosis,ScndVaccDate,ThrdVaccDate)]
tmp2[,mtime := AdmissionDate]
Mg <- tmp2[tmp1, on = c("StudyID","mtime"), roll = -15]

## If they are not matched then the admission date will be NA 
Mga <- Mg[!is.na(AdmissionDate) & 
            TestResult=="Positive" & 
            COVID19Diagnosis=="No" & 
            AdmissionDate > (ScndVaccDate + 14)]

# take the first occurrence 
setorder(Mga,AdmissionDate,StudyID)
Mga <- Mga[,.SD[1],by = StudyID] 
nrow(Mga)

## only 2,466 of these. 
s <- Mga$StudyID
Mga[StudyID==sample(s,1)]

## Calculate difference between test and admission date 
Mga[, diff := as.numeric(difftime(TestDate,AdmissionDate, units = "days"))]
summary(Mga$diff)

# labels events as in booster or non-booster phases
Mga[!is.na(ThrdVaccDate), phase := ifelse(AdmissionDate > (ThrdVaccDate + 14),"booster","scnd dose")]
Mga[is.na(ThrdVaccDate), phase := "scnd dose"]
Mga[,event:=1]
Mgb <- Mga[,.(StudyID,AdmissionDate,COVID19Diagnosis,event)]
setnames(Mgb, old = c("AdmissionDate","event","COVID19Diagnosis"), new = c("AdmissionDatePosTest","eventPosTest","COVID19DiagnosisPosTest"))

# merge into existing data set 
DT <- merge(DT,Mgb,by = c("StudyID"), all.x=T)

res <- DT[eventPosTest==1,length(AdmissionDatePosTest),by = .(month,phase)]
setorder(res,month,phase)
res

################################################################################
####### Hospital admission and positive test before discharge ##################
################## but not defined as UKHSA COVID cause   ######################
Ng <- tmp2[tmp1, on = c("StudyID","mtime"),roll=T, rollends=c(T,T)]
Nga <- Ng[TestDate <= DischargeDate & 
            TestDate > AdmissionDate & 
            TestResult=="Positive" & 
            COVID19Diagnosis=="No" & 
            TestDate > (ScndVaccDate + 14)]
setorder(Nga,AdmissionDate,StudyID)
Nga <- Nga[,.SD[1],by = StudyID] 
nrow(Nga)

## Calculate difference between test and admission date 
Nga[, diff := as.numeric(difftime(TestDate,AdmissionDate, units = "days"))]

# labels events as in booster or non-booster phases
Nga[!is.na(ThrdVaccDate), phase := ifelse(TestDate > (ThrdVaccDate + 14),"booster","scnd dose")]
Nga[is.na(ThrdVaccDate), phase := "scnd dose"]
Nga[,event:=1]
Ngb <- Nga[,.(StudyID,TestDate,COVID19Diagnosis,event)]
setnames(Ngb, old = c("TestDate","event","COVID19Diagnosis"), new = c("TestDatePosAdmis","eventTestPosAdmis","COVID19DiagnosisPosAdmis"))

# merge into existing data set 
DT <- merge(DT,Ngb,by = c("StudyID"), all.x=T)

#DT[eventPosAdmis==1,StudyID][2]
Nga[StudyID == 7801074]
#DT[StudyID==15770]
#AH[StudyID==15770]
#Test[StudyID==15770]

################################################################################
########## COVID Death as the event ############################################
################ 2,819 from COVID ##############################################

DT[!is.na(DateOfDeath), table(CauseOfDeath)]
nrow(DT[!is.na(DateOfDeath), DateOfDeath < (ScndVaccDate + 14)])

## all had previous hospital admission or positive test
DT[CauseOfDeath == "COVID Confirmed",table(event)]
DT[CauseOfDeath == "COVID Confirmed", sum(event),by = COVID19Diagnosis]
DT[,COVIDDeath := ifelse(CauseOfDeath == "COVID Confirmed",1,0)]

## Define the event time 
DT[!is.na(AdmissionDate) & is.na(AdmissionDatePosTest),event.date := AdmissionDate]
DT[!is.na(AdmissionDatePosTest) & is.na(AdmissionDate),event.date := AdmissionDatePosTest]
DT[!is.na(AdmissionDatePosTest) & !is.na(AdmissionDate),event.date := pmin(AdmissionDatePosTest,AdmissionDate)]
DT[is.na(AdmissionDatePosTest) & is.na(AdmissionDate) & COVIDDeath==1,event.date := DateOfDeath]
nrow(DT[!is.na(event.date)])

DT[!is.na(TestDatePosAdmis) & is.na(event.date), event.date:= TestDatePosAdmis]
DT[!is.na(TestDatePosAdmis) & !is.na(AdmissionDate),event.date := pmin(TestDatePosAdmis,AdmissionDate)]

# Some deaths occur before hospitalisation
nrow(DT[CauseOfDeath=="COVID Confirmed" & DateOfDeath < event.date])
DT[CauseOfDeath=="COVID Confirmed" & DateOfDeath < event.date, event.date:= DateOfDeath]

#### Composite event indicator 
DT[,event := 0]

nrow(DT[!is.na(AdmissionDate)])
nrow(DT[!is.na(AdmissionDate) | !is.na(AdmissionDatePosTest)])
nrow(DT[!is.na(AdmissionDate) | !is.na(AdmissionDatePosTest)  | COVIDDeath==1])
nrow(DT[is.na(AdmissionDate) & is.na(AdmissionDatePosTest)  & COVIDDeath==1])
DT[!is.na(AdmissionDate) | !is.na(TestDatePosAdmis) | COVIDDeath==1, event:=1]
table(DT$event)

################################################################################
################### Exposure Definition ########################################
DT[,exposure := 0]

# If they have a third vaccination dose date then assign to exposure  = 1
DT[!is.na(ThrdVaccDate), exposure := 1]
DT[exposure == 1, tstart := ThrdVaccDate + 14]
DT[exposure == 0, tstart := ScndVaccDate + 14]

## these received third dose but died before the 14 day run in period - reassigned to the exposure group
length(DT[!is.na(ThrdVaccDate) & DateOfDeath > (ScndVaccDate + 14) & DateOfDeath <= tstart, StudyID])
DT[!is.na(ThrdVaccDate) & DateOfDeath > (ScndVaccDate + 14) & DateOfDeath <= tstart,`:=`(exposure= 0, tstart = ScndVaccDate + 14)]

## these had Covid hospitalisation before their third dose 
DT[!is.na(ThrdVaccDate) & event.date > (ScndVaccDate + 14) & event.date <= tstart, StudyID][1]
DT[!is.na(ThrdVaccDate) & event.date > (ScndVaccDate + 14) & event.date <= tstart,`:=`(exposure= 0, tstart = ScndVaccDate + 14)]

## these has Covid positive admission 
DT[!is.na(ThrdVaccDate) & event.date > (ScndVaccDate + 14) & event.date <= tstart,`:=`(exposure= 0, tstart = ScndVaccDate + 14)]

DT[StudyID==739,.(exposure,tstart,event.date,ScndVaccDate,AdmissionDate,AdmissionDatePosTest,event)]

# Sanity check  - Hospitlised before their second dose run in period began?  
DT[event.date <= (ScndVaccDate + 14)]

DT[exposure == 1 & event.date < tstart]

DT[DateOfDeath <= (ScndVaccDate + 14), StudyID]

################################################################################
################### time in study ##############################################

DT[event == 1, time := difftime(event.date, tstart,units="days")]
DT[exposure == 1 & event == 1 & is.na(event.date)]

## No event but died - make event date the last.day 
DT[event ==0 & !is.na(DateOfDeath) ,event.date := DateOfDeath]
DT[event==0 & is.na(DateOfDeath), event.date := last.day]

# no event but died, censored at time of death
nms <- c("exposure","tstart","SecondVaccinationDoseDate","ThirdVaccinationDoseDate","event","DateOfDeath","time")
DT[event == 0 & !is.na(DateOfDeath), time := difftime(event.date,tstart,units = "days")]

# no event and no death recorded - censored at end of study date
DT[event == 0 & is.na(DateOfDeath), time := difftime(event.date,tstart,units = "days")]

DT[is.na(time),table(COVIDDeath)]

# define the background infection rate variable - based on ONS infection survey. 
onsbreaks <- c(as.Date("2020-12-07"),sp)
DT[,onsdate := cut(tstart,breaks = onsbreaks)]
DT[,infrate:= idt$c.infrate[match(onsdate,as.character(idt$sp))]]

# Categories for age. 
Agbr <- c(18,50,55,60,65,70,75,80,111)
Aglbs <- c("18-49","50-54","55-59","60-64","65-69","70-74","75-79","80+")
DT[,c.Age := cut(Age,breaks = Agbr, labels = Aglbs, right=F)]
cut(c(50,55,60,65,70,75,80),breaks = Agbr, labels = Aglbs,right=F)

Agbr2 <- c(18,40,50,65,75,120)
Aglbs2 <- c("18-39","40-49","50-64","64-74","75+")
DT[,c.Age2 := cut(Age,breaks = Agbr2, labels = Aglbs2, right=F)]
cut(c(18,40,50,65,75),breaks = Agbr2, labels = Aglbs2,right=F)

# Gap between first and second vaccine doses - 
DT[,GapVaccDose := difftime(ScndVaccDate,FirstVaccinationDoseDate)]
Gbrs <- c(0,7*7,9*7,11*7,13*7,Inf)
Glbs <- c("0-6 weeks","7-8 weeks","9-10 weeks","11-12 weeks",">=13 weeks")
DT[,GapVaccDose := cut(as.numeric(GapVaccDose),breaks = Gbrs, labels = Glbs, right = F)]
table(DT$GapVaccDose)

# From date of second dose, time since previous positive test 
TSX <- merge(DT[,.(StudyID,ScndVaccDate,exposure)],TS,by = "StudyID")
TSX3 <- TSX[TestDate <= ScndVaccDate & TestResult=="Positive"]
tmp <- TSX3[,.(LastPosTest = tail(TestDate,1), SVD = ScndVaccDate[1]), by = StudyID]
tmp[,GapPrevHistC19 := as.numeric(difftime(SVD,LastPosTest,unit = "days"))]
tmp[,SVD := NULL]
DT <- merge(DT,tmp, by = "StudyID", all.x = T)
DT[is.na(GapPrevHistC19), GapPrevHistC19 := 0]
PHbrs <- c(0,0.001,3*30,6*30,9*30,2000)
PHlbs <- c("No prior infection","0-2 months","3-5 months","6-8 months",">=9 months")
DT[,GapPrevHistC19 := cut(GapPrevHistC19,breaks = PHbrs, labels = PHlbs, right = F)]


# create factor for IMD and make 5 (the least deprived) the reference level. 
DT[,IMDQF := factor(IMDQuintile,levels=1:5,labels = c("1-Most","2","3","4","5-Least"),ordered = TRUE)]

# Urban rural 
DT[,UrbanRural := factor(UrbanRural)]
DT[,UrbanRural := relevel(UrbanRural, ref = "Urban")]


# Number of PCR tests prior to second dose 
TS <- merge(Test,Symp,by = c("StudyID","TestDate"),all.x=T)
tmp <- merge(DT[,.(StudyID,ScndVaccDate,exposure)],TS,by = "StudyID")
tmp <- tmp[TestDate < ScndVaccDate]
tmp <- tmp[,.(NC19Tests= .N), by = StudyID]
DT <- merge(DT,tmp, by = "StudyID", all.x = T)
DT[is.na(NC19Tests), NC19Tests := 0]
PTbrs <- c(0,1,2,3,5,10,Inf)
PTlbs <- c("0","1","2","3-4","5","10+")
DT[,NC19Tests := cut(NC19Tests, breaks = PTbrs, labels = PTlbs, right = F,ordered_result=T)]
table(DT$NC19Tests, useNA="always")

# BMI 
bmibr <- c(0,18.5,25,30,35,40,Inf)
bmilbs <- c("<18.5","18.5-24.9","25.0-29.9","30.0-34.9","35.0-39.9",">=40")
DT[,c.BMI := cut(BMILatest,c(0,1),breaks = bmibr, labels = bmilbs, right = F)]
DT[is.na(c.BMI),c.BMI:="Missing"]
DT[,c.BMI := factor(c.BMI)]
DT[,c.BMI := relevel(c.BMI, ref = "18.5-24.9")]
table(DT$c.BMI)

# number of risk groups 
#nms <- names(DT)[grep("RiskGroup",names(DT))]
#tmp <- DT[,..nms]
#DT[,NRiskGrps := rowSums(tmp)]
#NRlbs <- c("0","1","2","3","4","5+")
#DT[,NumRiskGrps := cut(NRiskGrps,breaks= c(0,1,2,3,4,5,Inf),labels = NRlbs, right =F)]

# number of QCovid risk groups 
nms <- names(DT)[c(grep("QConditions",names(DT)),grep("QDrugs",names(DT)))]
DT[,(nms) := lapply(.SD,function(x){ifelse(!is.na(x),1,0)}), .SDcols = nms]
table(DT$QConditions_AsthmaLatestDate)
tmp <- DT[,..nms]
DT[,NRiskGrps := rowSums(tmp)]
NRlbs <- c("0","1","2","3","4","5+")
DT[,NumRiskGrps := cut(NRiskGrps,breaks= c(0,1,2,3,4,5,Inf),labels = NRlbs, right =F)]
DT[,table(NumRiskGrps)]

# Create missing level for Ethnicity and make white the reference level
DT[is.na(EthnicityCode), EthnicityCode := "X"] # code X for missing 
DT[,EthnicityCode := factor(EthnicityCode)]
DT[,EthnicityCode := relevel(EthnicityCode, ref = "W")]

# Smoking Status. Create missing level and make non-smoker the reference level 
DT[is.na(SmokingStatus_Category), SmokingStatus_Category:= "X"] # code X for missing 
DT[, SmokingStatus_Category := factor(SmokingStatus_Category)]
DT[, SmokingStatus_Category := relevel(SmokingStatus_Category, ref = "Non-Smoker")]

# Shielding : Flag if shielding earliest date precedes vacc dose date 
DT[,Shielding := 0]
DT[!is.na(ShieldingEarliestDate),Shielding := ifelse(ShieldingEarliestDate < (tstart-14),1,0)]

# Defining the brand specific dose. 
DT[ThirdVaccinationBrand=="AstraZeneca" & exposure==1,exposure.brsp:=1]
DT[ThirdVaccinationBrand=="Moderna" & exposure == 1, exposure.brsp:=2]
DT[ThirdVaccinationBrand=="Pfizer-BioNTech" & exposure == 1, exposure.brsp:=3]
DT[is.na(exposure.brsp) & exposure ==1, exposure.brsp:=4]
DT[is.na(exposure.brsp), exposure.brsp := 0]
DT[,exposure.brsp := factor(exposure.brsp, levels = 0:4,labels =c("No booster","AZ","MD","PB","UNK"))]

# Multi-morbidity score 
cut(1,c(-Inf,0,1,2,Inf), labels = c("<=0","0.01 - 1","1.01 - 2",">2"),right = T)
DT[,c.CMMS := cut(CMMSWithoutAgeAdjustment,c(-Inf,0,1,2,Inf), labels = c("<=0","0.01 - 1","1.01 - 2",">2"),right = T)]

# HIV 
DT[is.na(HIVEarliestDate), HIVStatus := 0]
DT[!is.na(HIVEarliestDate) | !is.na(HIVLatestDate), HIVStatus := 1]

# Background interval week 
#DT[,bg_interval := NULL]
DT[event==1 & event.date > as.Date("2021-10-14") & event.date <= as.Date("2021-12-31"), bg_interval := as.numeric(week(event.date))]
DT[event==1 & event.date > as.Date("2021-12-31"), bg_interval := as.numeric(week(event.date)+53)]
DT[event==0, bg_interval:=0]
DT[,bg_interval:=factor(bg_interval, levels=c(0,42:66),labels=c("Ref",paste0("wk",42:66)))]
table(DT$bg_interval)

# Study entry data and person-years calculation 
# To allow for at least a week of follow up - tstart cannot be after last.day -14
idx <- DT[tstart >= (last.day - 7), StudyID]
DT <- DT[!{StudyID %in% idx}]

DT[,ipy := as.numeric(time)/365.25]

## Checks
DT[ipy== 0,StudyID][1]
DT[time==0]
DT[time<0]

DT[,sum(event), by = exposure]

# Save the data for shortcuts 
saveRDS(object = DT, file = "data/analysis_cohort.rds")



