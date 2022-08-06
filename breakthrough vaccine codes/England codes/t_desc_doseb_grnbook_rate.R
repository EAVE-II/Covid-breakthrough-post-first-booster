# t_desc_doseb_GreenBook_rate

DT <- readRDS("data/analysis_cohort.rds")

## change thes e back to unordered factors (stop R fitting poly terms)
DT[,IMDQF := factor(IMDQF, ordered =F)]
DT[,IMDQF := relevel(IMDQF, ref = "5-Least")]
DT[,NC19Tests := factor(NC19Tests, ordered = F)]
DT[,NC19Tests := relevel(NC19Tests,ref = "0")]

DTx <- DT[exposure==1]

x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy))]
x[,rate := (events/personyrs)*1000]
A <- data.frame(Grp="Total",levels = "", N = x$N, events = x$events,personyrs = x$personyrs,rate=x$rate)
lbs <- c("Asplenia","Asthma","Severe Asthma","Coronary Heart disease","Chronic Kidney disease","Chronic Liver disease",
         "Chronic Respiratory disease","Diabetes","Immuno-compromised","Neurological disease","Severe mental illness",
         "Flu risk","Shielding","Wider learning disability","Morbid Obesity","HIV")
nms <- c("Asplenia_RiskGroup","Asthma_RiskGroup","SevereAsthma_RiskGroup","CHD_RiskGroup","ChronicKidney_RiskGroup",
         "ChronicLiver_RiskGroup","ChronicResp_RiskGroup","Diabetes_RiskGroup",  
         "Immuno_RiskGroup","Neuro_RiskGroup","SevereMentalIllness_RiskGroup", "FluRiskFactor",   
         "Shielding","WiderLearningDisability_RiskGroup",
         "MorbidObesity_RiskGroup","HIVStatus")
for (i in 1:length(lbs)){
  x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy)), by = eval(nms[i])]
  x[,rate := (events/personyrs)*1000]
  setorder(x)
  setnames(x,nms[i],"levels")
  A <- rbind(A,data.frame(Grp = lbs[i],x))
}  
View(A)

write.csv(A, file = "results/t_desc_doseb_grnbook_rate.csv")
