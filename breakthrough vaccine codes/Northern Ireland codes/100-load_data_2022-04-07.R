#source("scripts/00-libraries_and_functions.R")

## 10 Apr - added the prefix '../' to the file location to allow one version of this part 
## and us all to pull into our projects 

dt_demographics <- readRDS("../input/clean/demographics.rds")
dt_covid_vaccines_wide <- readRDS("../input/clean/covid_vaccines_wide.rds")
dt_admissions <- readRDS("../input/clean/admissions.rds")
dt_deaths <- readRDS("../input/clean/deaths.rds")
dt_tests <- readRDS("../input/clean/tests.rds")
dt_epd <- readRDS("../input/clean/epd.rds")
dt_care_home_residents_pillar_2 <- readRDS("../input/clean/care_home_residents_pillar2.rds")
dt_cog_uk <- readRDS("../input/Cog_UK_Feb22_1.rds")


message("****All files read in - move to script 101*****")
