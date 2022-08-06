# Poisson model coef booster dose cohort QCovid Risk factors.  
# 
DT <- readRDS("data/analysis_cohort.rds")
table(DT$NC19Tests)

## change thes e back to unordered factors (stop R fitting poly terms)
DT[,IMDQF := factor(IMDQF, ordered =F)]
DT[,IMDQF := relevel(IMDQF, ref = "5-Least")]
DT[,NC19Tests := factor(NC19Tests, ordered = F)]
DT[,NC19Tests := relevel(NC19Tests,ref = "0")]

# Booster only 
omd <- as.Date("2021-12-20")
DTx <- DT[event.date > omd & exposure.brsp %in% c("MD","PB")]

# expand the data set in the exposed. 
ts <- c(28-14,49-14,1000)
tlbs <- c("week 2 - 4","week 5-7","week 8+")
f0 <- as.formula(Surv(time,event) ~ QConditions_AsthmaLatestDate +                        
                 QConditions_AtrialFibrillationLatestDate +            
                 QConditions_BloodBoneCancerLatestDate +               
                 QConditions_BoneMarrowLatestDate +                    
                 QConditions_CerebralPalsyLatestDate +                 
                 QConditions_CF_BronchiectasisLatestDate +             
                 QConditions_ChemotherapyLatestDate +                  
                 QConditions_ChronicKidneyLatestDate +                 
                 QConditions_CirrhosisLatestDate +                     
                 QConditions_CongenitalHeartLatestDate +               
                 QConditions_COPDLatestDate +                          
                 QConditions_CoronaryHeartLatestDate +              
                 QConditions_DementiaLatestDate +                      
                 QConditions_DiabetesLatestDate +                      
                 QConditions_EpilepsyLatestDate +                      
                 QConditions_FractureLatestDate +                      
                 QConditions_HeartFailureLatestDate +                  
                 QConditions_LearningDisabilityCardiacLatestDate +     
                 QConditions_LungOralCancerLatestDate +                
                 QConditions_MND_MS_MYAS_HuntingtonsLatestDate +       
                 QConditions_ParkinsonsLatestDate +                    
                 QConditions_PeripheralVascularLatestDate +            
                 QConditions_PulmonaryHypeFibrosisLatestDate +         
                 QConditions_RA_SLE_SeronegativeArthiritisLatestDate + 
                 QConditions_RadiotherapyLatestDate +                  
                 QConditions_SevereMental_IllnessLatestDate +          
                 QConditions_SickleCell_ImmunodeficiencyLatestDate +   
                 QConditions_SolidOrganTransplantLatestDate +          
                 QConditions_StrokeTIALatestDate +                     
                 QConditions_ThrombosisPulmonaryEmbolusLatestDate +    
                 QDrugs_AntiLeukotriene_LABALatestDate +               
                 QDrugs_ImmunosuppressantsLatestDate +                 
                 QDrugs_OralSteroidsLatestDate + 
                 HIVStatus + 
                 Sex + 
                 GapPrevHistC19 + 
                 GapVaccDose + 
                 c.Age + 
                 IMDQF + 
                 UrbanRural + 
                 NC19Tests +
                 EthnicityCode + 
                 c.BMI + 
                 NewNHSRegion + 
                 StudyID)

DTx1 <- setDT(survSplit(f0, data = DTx,cut = ts, episode = "period",id="id"))
DTx1[,period := factor(period, levels = 1:3, labels = tlbs)]

# recalculate individual person-time 
DTx1[,ipy := (time - tstart)/365.25]

rskf <- c("QConditions_AsthmaLatestDate",                       
          "QConditions_AtrialFibrillationLatestDate",           
          "QConditions_BloodBoneCancerLatestDate",              
          "QConditions_BoneMarrowLatestDate",                   
          "QConditions_CerebralPalsyLatestDate",                
          "QConditions_CF_BronchiectasisLatestDate",            
          "QConditions_ChemotherapyLatestDate",                 
          "QConditions_ChronicKidneyLatestDate",                
          "QConditions_CirrhosisLatestDate",                    
          "QConditions_CongenitalHeartLatestDate",              
          "QConditions_COPDLatestDate",                         
          "QConditions_CoronaryHeartLatestDate",                
          "QConditions_DementiaLatestDate",                     
          "QConditions_DiabetesLatestDate",                     
          "QConditions_EpilepsyLatestDate",                     
          "QConditions_FractureLatestDate",                     
          "QConditions_HeartFailureLatestDate",                 
          "QConditions_LearningDisabilityCardiacLatestDate",    
          "QConditions_LungOralCancerLatestDate",               
          "QConditions_MND_MS_MYAS_HuntingtonsLatestDate",      
          "QConditions_ParkinsonsLatestDate",                   
          "QConditions_PeripheralVascularLatestDate",           
          "QConditions_PulmonaryHypeFibrosisLatestDate",        
          "QConditions_RA_SLE_SeronegativeArthiritisLatestDate",
          "QConditions_RadiotherapyLatestDate",                 
          "QConditions_SevereMental_IllnessLatestDate",         
          "QConditions_SickleCell_ImmunodeficiencyLatestDate",  
          "QConditions_SolidOrganTransplantLatestDate",         
          "QConditions_StrokeTIALatestDate",                   
          "QConditions_ThrombosisPulmonaryEmbolusLatestDate",
          "QDrugs_AntiLeukotriene_LABALatestDate",              
          "QDrugs_ImmunosuppressantsLatestDate",                
          "QDrugs_OralSteroidsLatestDate",
          "HIVStatus")

adj <- c("GapPrevHistC19","GapVaccDose","Sex","c.Age","IMDQF",
         "UrbanRural","NC19Tests","EthnicityCode","c.BMI","NewNHSRegion")

aggr <- c("period", rskf,adj)
ADT <- DTx1[,.(n.events = sum(event),
               n.ids = length(unique(StudyID)),
               pyrs = sum(ipy)),by = aggr] 

R <- data.frame(term=NULL,estimate=NULL,std.error=NULL,statistic=NULL,p.value=NULL)
j <- 18
for (j in 1:length(rskf)){
nms <- c("n.events","pyrs","period",rskf[j],adj)
tmp <- ADT[,..nms]  
m <- glm(n.events ~ .-pyrs, offset = log(pyrs), family = poisson, data = tmp)
mt <- setDT(tidy(m))
r <- mt[term==rskf[j]]
R <- rbind(R,r) 
cat("iteration # ",j,"\n")
}
R
write.csv(R,file = "results/FINAL/t_pois_coef_doseb_QCovidmodels.csv")


### Alternative way - one at a time
# j <- 10
# aggr <- c("period", rskf[j],adj)
# ADT <- DTx1[,.(n.events = sum(event),
#                n.ids = length(unique(StudyID)),
#                pyrs = sum(ipy)),by = aggr] 
# nms <- c("n.events","pyrs","period",rskf[j],adj)
# tmp <- ADT[,..nms]  
# m <- glm(n.events ~ . -pyrs, offset = log(pyrs), family = poisson, data = tmp)
# mt <- setDT(tidy(m))
# mt[term==rskf[j]] 
# summary(m)






