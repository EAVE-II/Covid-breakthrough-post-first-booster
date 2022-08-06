#################################################################
#### Remove novavax and dose 2 moderna
#################################################################

## Align to the protocol and remove Moderna dose 2

dt_cohort = dt_cohort %>%
  filter(Dose2_name!="Novavax" & Dose2_name!="Moderna")

# Check
#z = dt_cohort %>%
 # group_by(Dose2_name) %>%
  #tally()


#############################################
#### Clean days between dose 1 and 2
#############################################

# Exclude those with dose 2 before dose 1 and within 21 days

dt_cohort <- dt_cohort %>% 
  filter(days_btn_1_2 >= min_days_between_first_and_second_dose) %>%
  mutate(wks_btn_1_2 = cut(days_btn_1_2, breaks=c(1,41,55,69,83,max(days_btn_1_2)), 
                            labels=c("3-6 weeks","7-8 weeks","9-10 weeks","11-12 weeks","13+ weeks")))

z = dt_cohort %>%
  dplyr::select(study_id, days_btn_1_2, wks_btn_1_2)

dt_cohort$wks_btn_1_2 <- factor(dt_cohort$wks_btn_1_2, ordered = TRUE,
                                   levels = c("3-6 weeks", "7-8 weeks", "9-10 weeks",
                                              "11-12 weeks", "13+ weeks"))
#z = dt_cohort %>%
 # group_by(wks_btn_1_2) %>%
#  summarise(min = min(days_btn_1_2), max = max(days_btn_1_2))


# Save =========================================================================
cat("Save\n")

qsave(
  dt_cohort,
  file = ("input/dt_cohort.qs")
)

message("***Now removed Novavax and Dose 2 Moderna and applied a 21 day rule between dose 1 and 2***")


