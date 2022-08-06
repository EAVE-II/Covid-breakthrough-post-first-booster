## table descriptives omicron rate. 
rm(list=ls())
DT <- readRDS("data/analysis_cohort.rds")

# Overall 
omd <- as.Date("2021-12-20")
DTx <- DT[event.date > omd & exposure.brsp %in% c("MD","PB")]

# Recalculate person-years follow-up if there start date was before start of period
DTx[tstart < omd, ipy := as.numeric(difftime(event.date,omd,units = "days"))/365.25]

x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy))]
x[,rate := (events/personyrs)*1000]
A <- data.frame(Grp="bg_variant",levels = "Omicron", N = x$N, events = x$events,personyrs = x$personyrs,rate=x$rate)
nms <- c("Sex","c.Age","EthnicityCode","NumRiskGrps","c.BMI","IMDQF",
         "UrbanRural","NC19Tests","GapPrevHistC19","GapVaccDose","NewNHSRegion")
lbs <- c("Sex","Age","Ethnicity","#Risk groups","Body mass index",
         "Deprivation index","Urban/Rural","# PCR tests","Prior infection","Vaccine gap","NHS Region")

for (i in 1:length(lbs)){
  x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy)), by = eval(nms[i])]
  x[,rate := (events/personyrs)*1000]
  setorder(x)
  setnames(x,nms[i],"levels")
  A <- rbind(A,data.frame(Grp = lbs[i],x))
}  

View(A)
write.csv(A, file = "results/FINAL/t_desc_doseb_omicron_rate.csv")

# Pfizer-BioNTech 
omd <- as.Date("2021-12-20")
DTx <- DT[event.date > omd & exposure.brsp=="PB"]

# Recalculate person-years follow-up if there start date was before start of period
DTx[tstart < omd, ipy := as.numeric(difftime(event.date,omd,units = "days"))/365.25]

x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy))]
x[,rate := (events/personyrs)*1000]
APB <- data.frame(Grp="bg_variant",levels = "Omicron", N = x$N, events = x$events,personyrs = x$personyrs,rate=x$rate)
for (i in 1:length(lbs)){
  x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy)), by = eval(nms[i])]
  x[,rate := (events/personyrs)*1000]
  setorder(x)
  setnames(x,nms[i],"levels")
  APB <- rbind(APB,data.frame(Grp = lbs[i],x))
}  

View(APB)
write.csv(APB, file = "results/FINAL/t_desc_doseb_omicron_rate_PB.csv")

# Moderna booster. 
omd <- as.Date("2021-12-20")
DTx <- DT[event.date > omd & exposure.brsp=="MD"]

# Recalculate person-years follow-up if there start date was before start of period
DTx[tstart < omd, ipy := as.numeric(difftime(event.date,omd,units = "days"))/365.25]

x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy))]
x[,rate := (events/personyrs)*1000]
AMD <- data.frame(Grp="bg_variant",levels = "Omicron", N = x$N, events = x$events,personyrs = x$personyrs,rate=x$rate)
for (i in 1:length(lbs)){
  x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy)), by = eval(nms[i])]
  x[,rate := (events/personyrs)*1000]
  setorder(x)
  setnames(x,nms[i],"levels")
  AMD <- rbind(AMD,data.frame(Grp = lbs[i],x))
}  
View(AMD)
write.csv(AMD, file = "results/FINAL/t_desc_doseb_omicron_rate_MD.csv")


