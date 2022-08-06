## This is the t_pois_coef_doseb analysis 
#source("src/Define the cohort.R")
DT <- readRDS("data/analysis_cohort.rds")

## change thes e back to unordered factors (stop R fitting poly terms)
DT[,IMDQF := factor(IMDQF, ordered =F)]
DT[,IMDQF := relevel(IMDQF, ref = "5-Least")]
DT[,NC19Tests := factor(NC19Tests, ordered = F)]
DT[,NC19Tests := relevel(NC19Tests,ref = "0")]

DTx <- DT[exposure.brsp=="PB"]
table(DTx$exposure.brsp)

# split, so that follow-up time can be different. 
ts <- c(28-14,49-14,1000)
tlbs <- c("week 2 - 4","week 5-7","week 8+")
f0 <- as.formula(Surv(time,event) ~ exposure + GapPrevHistC19 +  GapVaccDose + 
                   Sex + c.Age + IMDQF + UrbanRural + NC19Tests + EthnicityCode + 
                   NumRiskGrps + c.BMI + NewNHSRegion + StudyID)
DTx1 <- setDT(survSplit(f0, data = DTx,cut = ts, episode = "period",id="id"))
DTx1[,period := factor(period, levels = 1:3, labels = tlbs)]

# recalculate individual person-time 
DTx1[,ipy := (time - tstart)/365.25]
aggr <- c("period","GapPrevHistC19","GapVaccDose","Sex","c.Age",
          "IMDQF","UrbanRural","NumRiskGrps","NC19Tests","EthnicityCode",
          "c.BMI","NewNHSRegion")

ADT <- DTx1[,.(n.events = sum(event),
               n.ids = length(unique(StudyID)),
               pyrs = sum(ipy)),by = aggr]           

f <- as.formula(n.events ~ period + c.Age + GapPrevHistC19 + GapVaccDose + 
                  Sex + c.Age + IMDQF + UrbanRural + NumRiskGrps  + NC19Tests + EthnicityCode +  
                  c.BMI + NewNHSRegion)
m <- glm(f, offset = log(pyrs), family = poisson, data = ADT)
su <- summary(m)
pvs <- su$coefficients[,'Pr(>|z|)']
su
mt <- tidy(m)
mte <- tidy(m, exponentiate = T)
View(mt)
write.csv(mt,file = "results/FINAL/t_pois_coef_doseb_PB.csv" )



