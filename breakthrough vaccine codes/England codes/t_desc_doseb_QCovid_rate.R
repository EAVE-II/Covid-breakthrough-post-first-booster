# t_desc_doseb_QCovid_rate
library(stringr)
DT <- readRDS("data/analysis_cohort.rds")

## change thes e back to unordered factors (stop R fitting poly terms)
DT[,IMDQF := factor(IMDQF, ordered =F)]
DT[,IMDQF := relevel(IMDQF, ref = "5-Least")]
DT[,NC19Tests := factor(NC19Tests, ordered = F)]
DT[,NC19Tests := relevel(NC19Tests,ref = "0")]

omd <- as.Date("2021-12-20")
DTx <- DT[event.date > omd & exposure.brsp %in% c("MD","PB")]

x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy))]
x[,rate := (events/personyrs)*1000]
A <- data.frame(Grp="Total",levels = "", N = x$N, events = x$events,personyrs = x$personyrs,rate=x$rate)
names(DT)
nms <- c(names(DT)[c(grep("QConditions",names(DT)),grep("QDrugs",names(DT)))],"HIVStatus")
lbs <- gsub(".*_","",nms)
lbs <- str_remove(lbs,"LatestDate")

for (i in 1:length(lbs)){
  x <- DTx[,.(N = length(StudyID),events = sum(event), personyrs = sum(ipy)), by = eval(nms[i])]
  x[,rate := (events/personyrs)*1000]
  setorder(x)
  setnames(x,nms[i],"levels")
  A <- rbind(A,data.frame(Grp = lbs[i],x))
}  
View(A)

write.csv(A, file = "results/FINAL/t_desc_doseb_QCovid_rate.csv")

