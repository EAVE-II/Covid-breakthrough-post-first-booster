## This is the t_pois_coef_dose2 analysis 
#source("src/Define the cohort.R")
DT <- readRDS("data/analysis_cohort.rds")

## change thes e back to unordered factors (stop R fitting poly terms)
DT[,IMDQF := factor(IMDQF, ordered =F)]
DT[,IMDQF := relevel(IMDQF, ref = "5-Least")]
DT[,NC19Tests := factor(NC19Tests, ordered = F)]
DT[,NC19Tests := relevel(NC19Tests,ref = "0")]

DT <- DT[SecondVaccinationBrand=="Pfizer-BioNTech"]

# split, so that follow-up time can be different. 
# expand the data set in the exposed. 
DTe1 <- DT[exposure==1]
DTe1 <- DTe1[exposure.brsp %in% c("MD","PB")]
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
tlbs <- c("week 3-9","week 10-19","week 20+")
DTe0 <- setDT(survSplit(f0, data = DTe0,cut = ts, episode = "period",id="id"))
DTe0[,period := factor(period, levels = 1:3, labels = tlbs)]

# recalculate individual person-time 
DTe0[,ipy := (time - tstart)/365.25]

## merge back together 
DTs <- rbind(DTe1,DTe0)

## combine period and vacc dose categories - easier than doing interactions. 
DTs[,vaccInterv := paste0(exposure.brsp,".",period)]
DTs[,vaccInterv := factor(vaccInterv)]
DTs[,vaccInterv := relevel(vaccInterv, ref = "No booster.week 3-9")]
table(DTs$vaccInterv)

aggr <- c("vaccInterv","GapPrevHistC19","GapVaccDose","Sex","c.Age",
          "IMDQF","UrbanRural","NumRiskGrps","NC19Tests","EthnicityCode",
          "c.BMI","NewNHSRegion")

DTx <- DTs[,.(n.events = sum(event),
              n.ids = length(unique(StudyID)),
              pyrs = sum(ipy)),by = aggr]           

f <- as.formula(n.events ~ vaccInterv + c.Age + GapPrevHistC19 + GapVaccDose + 
                  Sex + IMDQF + UrbanRural + NumRiskGrps  + NC19Tests + EthnicityCode +
                  c.BMI + NewNHSRegion)

m <- glm(f, offset = log(pyrs), family = poisson, data = DTx)
su <- summary(m)
pvs <- su$coefficients[,'Pr(>|z|)']
su
mt <- tidy(m)
mte <- tidy(m, exponentiate = T)
View(mt)
#write.csv(mt,file = "results/FINAL/t_pois_coef_dose2_PB.csv" )



