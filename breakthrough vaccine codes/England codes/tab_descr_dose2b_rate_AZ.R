# PB only 
#source("src/Define the cohort.R")
DT <- readRDS("data/analysis_cohort.rds")

## change thes e back to unordered factors (stop R fitting poly terms)
DT[,IMDQF := factor(IMDQF, ordered =F)]
DT[,IMDQF := relevel(IMDQF, ref = "5-Least")]
DT[,NC19Tests := factor(NC19Tests, ordered = F)]
DT[,NC19Tests := relevel(NC19Tests,ref = "0")]

DT <- DT[SecondVaccinationBrand=="AstraZeneca"]
table(DT$SecondVaccinationBrand)

delta.period <- as.Date("2021-06-01")
omicron.period <- as.Date("2021-12-20")

# exclude the events happending prior to the start of the delta period. 
DTx <- DT[event.date > delta.period & ScndVaccDate < omicron.period]

# if they didn't have the event and did not die - end date is start of omicron
DTx[event==0 & is.na(DateOfDeath),event.date := omicron.period]

# if they died after end of delta period 
DTx[DateOfDeath > omicron.period, `:=`(event.date = omicron.period, event = 0)]

# if they had an event after end of delta period 
DTx[event.date > omicron.period, `:=`(event.date = omicron.period, event = 0)]

# adjust the follow-up time 
DTx[event == 0 & is.na(DateOfDeath),ipy := as.numeric(difftime(omicron.period,delta.period,units = "days"))/365.25]
DTx[DateOfDeath < omicron.period, ipy := as.numeric(difftime(DateOfDeath,delta.period,units = "days"))/365.25]

# follow-up time is since ScndVaccdose in everyone  
DTx[ScndVaccDate > delta.period, `:=`(ipy = as.numeric(difftime(event.date,ScndVaccDate,units = "days"))/365.25)]
DTx[ScndVaccDate <= delta.period, `:=`(ipy = as.numeric(difftime(event.date,delta.period,units = "days"))/365.25)]

DTx[ThrdVaccDate >= omicron.period & 
      ScndVaccDate < delta.period, StudyID][1]
DTx[StudyID==2,.(event.date,event,exposure,DateOfDeath,tstart,ipy*365.241,ThrdVaccDate,ScndVaccDate,ipy)]

x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy))]
x[,rate := (events/personyrs)*1000]
x

## Omicron period 
DTy <- DT[event.date > omicron.period & exposure == 1]

# Should we recalculate person-years follow-up - note not correcting time variable in this
DTy[tstart < omicron.period, ipy := as.numeric(difftime(event.date,omicron.period,units = "days"))/365.25]
DTy[tstart >= omicron.period, ipy := as.numeric(difftime(event.date,tstart,units = "days"))/365.25]

DTy[,.(min(ipy),max(ipy))]

y <- DTy[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy))]
y[,rate := (events/personyrs)*1000]
y

##
A0 <- data.frame(Grp="bg_variant",levels = "Delta", N = x$N, events = x$events,personyrs = x$personyrs,rate=x$rate)
A0 <- rbind(A0,data.frame(Grp = "bg_variant",levels = "Omicron",y))
A0

################# Split follow up ##############################################

# expand the data set in the exposed. 
DT <- DT[!{exposure.brsp %in% c("AZ","UNK")}]
DTe1 <- DT[exposure==1]
DTe0 <- DT[exposure==0]

ts <- c(28-14,49-14,1000)
tlbs <- c("week 2-4 ","week 5-7","week 8+")
f0 <- as.formula(Surv(time,event) ~ exposure.brsp + GapPrevHistC19 +  GapVaccDose + 
                   Sex + c.Age + IMDQF + UrbanRural + NC19Tests + EthnicityCode + 
                   NumRiskGrps + c.BMI + NewNHSRegion + StudyID)
DTe1 <- setDT(survSplit(f0, data = DTe1,cut = ts, episode = "period",id="id"))
DTe1[,period := factor(period, levels = 1:3, labels = tlbs)]

# recalculate individual person-time 
DTe1[,ipy := (time - tstart)/365.25]

## the non-exposed 
tsc <- c(84 - 14,168 - 14,1000)
tlbs <- c("week 3 - 9","week 10 - 19","week 20+")
DTe0 <- setDT(survSplit(f0, data = DTe0,cut = ts, episode = "period",id="id"))
DTe0[,period := factor(period, levels = 1:3, labels = tlbs)]

# recalculate individual person-time 
DTe0[,ipy := (time - tstart)/365.25]

## merge back together 
DTs <- rbind(DTe1,DTe0)

## combine period and vacc dose categories - easier than doing interactions. 
DTs[,vaccInterv := paste0(exposure.brsp,".",period)]
DTs[,vaccInterv := factor(vaccInterv)]
DTs[,vaccInterv := relevel(vaccInterv, ref = "No booster.week 3 - 9")]
table(DTs$vaccInterv)

DTx <- DTs[,.(events = sum(event),
              N = length(unique(StudyID)),
              personyrs = sum(ipy)),by = vaccInterv]           

DTx[, rate := (events/personyrs)*1000]
setnames(DTx,old = "vaccInterv", new = "levels")
A1 <- data.frame(Grp="vacc_dose_interval",DTx)

# everything else 
nms <- c("Sex","c.Age","EthnicityCode","NumRiskGrps","c.BMI","IMDQF",
         "UrbanRural","NC19Tests","GapPrevHistC19","GapVaccDose","NewNHSRegion")
lbs <- c("Sex","Age","Ethnicity","#Risk groups","Body mass index",
         "Deprivation index","Urban/Rural","# PCR tests","Prior infection","Vaccine gap","NHS Region")
for (i in 1:length(lbs)){
  x <- DT[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy)), by = eval(nms[i])]
  x[,rate := (events/personyrs)*1000]
  setorder(x)
  setnames(x,nms[i],"levels")
  A1 <- rbind(A1,data.frame(Grp = lbs[i],x))
}  
A1

AT <- rbind(A0,A1)
View(AT)
#write.csv(AT,file = "results/FINAL/t_desc_dose2b_rate_AZ.csv")

setDT(AT)
AT[Grp=="# PCR test"]
AT
