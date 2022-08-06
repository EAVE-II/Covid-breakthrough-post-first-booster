# OMICRON PERIOD ONLY 
# Rate ratios by clinical indicators by < 75 and over 75's 
DT <- readRDS("data/analysis_cohort.rds")

## change thes e back to unordered factors (stop R fitting poly terms)
DT[,IMDQF := factor(IMDQF, ordered =F)]
DT[,IMDQF := relevel(IMDQF, ref = "5-Least")]
DT[,NC19Tests := factor(NC19Tests, ordered = F)]
DT[,NC19Tests := relevel(NC19Tests,ref = "0")]

# Booster only 
DTx <- DT[exposure==1]

# expand the data set in the exposed. 
ts <- c(28-14,49-14,1000)
tlbs <- c("week 2 - 4","week 5-7","week 8+")
f0 <- as.formula(Surv(time,event) ~ Asplenia_RiskGroup + Asthma_RiskGroup + 
                   SevereAsthma_RiskGroup + CHD_RiskGroup + ChronicKidney_RiskGroup + 
                   ChronicLiver_RiskGroup + ChronicResp_RiskGroup + Diabetes_RiskGroup +  
                   Immuno_RiskGroup + Neuro_RiskGroup +  SevereMentalIllness_RiskGroup +   
                   FluRiskFactor + Shielding + WiderLearningDisability_RiskGroup + MorbidObesity_RiskGroup + 
                   HIVStatus + EthnicityCode + GapPrevHistC19 + Sex + IMDQF + UrbanRural + NewNHSRegion + StudyID)
DTx1 <- setDT(survSplit(f0, data = DTx,cut = ts, episode = "period",id="id"))
DTx1[,period := factor(period, levels = 1:3, labels = tlbs)]

# recalculate individual person-time 
DTx1[,ipy := (time - tstart)/365.25]
aggr <- c("period","Asplenia_RiskGroup","Asthma_RiskGroup", 
            "SevereAsthma_RiskGroup","CHD_RiskGroup","ChronicKidney_RiskGroup",
            "ChronicLiver_RiskGroup","ChronicResp_RiskGroup","Diabetes_RiskGroup",  
            "Immuno_RiskGroup","Neuro_RiskGroup","SevereMentalIllness_RiskGroup", "FluRiskFactor",   
            "Shielding","WiderLearningDisability_RiskGroup",
            "MorbidObesity_RiskGroup","HIVStatus","GapPrevHistC19","Sex","EthnicityCode",
            "IMDQF","UrbanRural",
              "NewNHSRegion")

ADT <- DTx1[,.(n.events = sum(event),
               n.ids = length(unique(StudyID)),
               pyrs = sum(ipy)),by = aggr]           

f <- as.formula(n.events ~ Asplenia_RiskGroup + Asthma_RiskGroup + 
                  SevereAsthma_RiskGroup + CHD_RiskGroup + ChronicKidney_RiskGroup + 
                  ChronicLiver_RiskGroup + ChronicResp_RiskGroup + Diabetes_RiskGroup +  
                  Immuno_RiskGroup + Neuro_RiskGroup +  SevereMentalIllness_RiskGroup +     
                  FluRiskFactor + Shielding + WiderLearningDisability_RiskGroup + 
                  MorbidObesity_RiskGroup + HIVStatus +  
                  EthnicityCode + GapPrevHistC19 + 
                  Sex + IMDQF + UrbanRural + NewNHSRegion)

m <- glm(f, offset = log(pyrs), family = poisson, data = ADT)
su <- summary(m)
pvs <- su$coefficients[,'Pr(>|z|)']
mt <- tidy(m)
mte <- tidy(m, exponentiate = T)
setDT(mt)
setDT(mte)
View(mte)
cis <- confint.default(m)
write.csv(mt,file = "results/t_pois_coef_doseb_grnbook.csv")
