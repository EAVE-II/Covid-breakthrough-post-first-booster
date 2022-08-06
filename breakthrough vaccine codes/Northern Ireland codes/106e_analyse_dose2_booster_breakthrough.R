cat("Clearing workspace and loading stuff\n")

source("scripts/000_libraries_and_functions.R")

gc()
memory.limit(64000)
options(scipen = 100)

# Load =========================================================================
cat("Load\n")

d_analysis <- qread("results/d_analysis_stats.qs") 

d_analysis %>% count(vacc_dose_interval)

# Count person-years ===========================================================
cat("Count person-years\n")

# x = d_analysis %>%
#     count(bg_variant)

# total person-years by vaccine dose and variant
t_pyears_vacc_variant <-
    d_analysis %>%
    calc_pyears(
        vacc_dose_interval,
        bg_variant
    ) %>%
    mutate(bg_variant = factor(bg_variant) %>% forcats::fct_relevel("Delta", "Omicron")) %>%
    arrange(bg_variant, vacc_dose_interval)

# x = t_pyears_vacc_variant %>%
#     count(bg_variant)

print(t_pyears_vacc_variant)

# Clean ========================================================================
cat("Clean\n")

d_analysis <-
    d_analysis %>%
    # select analysis variables
    dplyr::select(
        study_id,
        # survival analysis things
        tstart,
        tstop,
        event_flg,
        # vaccine
        vacc_dose2_name,
        vacc_dose_interval,
        # covariates
        sex,
        age_gp,
        vacc_dose1_dose2_diff_cat,
        dose2_prior_infection_cat,
        nimdm,
        ur_combined,
        BNF_group,
        n_tests_gp,
        lgd,
        # background
        bg_interval,
        bg_avg_infection,
        bg_variant,
        # study design
        sample_weight
    ) %>%
    # final edits
    mutate(
        bg_interval = factor(bg_interval),
        # remove vaccine name from dose 2 interval
        vacc_dose_interval = vacc_dose_interval %>%
            as.character() %>%
            stringr::str_replace("(dose2)_[A-Z]+_(day[0-9]+)", "\\1_\\2") %>%
            factor(),
        # squish booster dose intervals to match analysis plan
        # vacc_dose_interval = vacc_dose_interval %>%
        #     forcats::fct_collapse(
        #         doseb_MD_day028 = c(
        #             "doseb_MD_day028",
        #             "doseb_MD_day042",
        #             "doseb_MD_day056",
        #             "doseb_MD_day070"
        #         ),
        #         doseb_PB_day028 = c(
        #             "doseb_PB_day028",
        #             "doseb_PB_day042",
        #             "doseb_PB_day056",
        #             "doseb_PB_day070"
        #         )
        #    ),
        # set reference categories
        sex = forcats::fct_relevel(sex, "Female"),
        age_gp = forcats::fct_relevel(age_gp, "18-49"),
        bg_variant = forcats::fct_relevel(bg_variant, "Delta"),
        vacc_dose_interval = forcats::fct_relevel(vacc_dose_interval, "dose2_day014"),
        dose2_prior_infection_cat = forcats::fct_relevel(dose2_prior_infection_cat, "No prior infection"),
        nimdm = forcats::fct_relevel(nimdm, "5 - least deprived")
    )

#d_analysis %>% count(sex)

# Aggregate to person-years ====================================================
cat("Aggregate to person-years\n")

# t = d_analysis %>% calc_pyears(dose2_prior_infection_cat)

d_pyears <- d_analysis %>%
    calc_pyears(
        vacc_dose2_name,
        vacc_dose_interval,
        sex,
        age_gp,  # age_3cat age_4cat age_5y_cat
        vacc_dose1_dose2_diff_cat,
        dose2_prior_infection_cat,
        nimdm,
        ur_combined,
        BNF_group,
        n_tests_gp,
        lgd,
        bg_interval,
        bg_avg_infection,
        bg_variant
    )

# Fit Poisson model ====================================================
cat("Fit overall Poisson model\n")

# interval_check = d_pyears %>% count(vacc_dose_interval, event)
# write.csv(interval_check, "help_documents/interval_event_check.csv", row.names=FALSE)

frml <- event ~ vacc_dose_interval +
                        dose2_prior_infection_cat +   
                        vacc_dose1_dose2_diff_cat +
                        sex +
                        age_gp +
                        nimdm +
                        ur_combined +
                        BNF_group +
                        n_tests_gp +
                        lgd + #not presented in tables
                        bg_interval + #not presented in tables
                        offset(log(pyears))


pois_dose2 <- glm(
    data    = d_pyears,
    family  = poisson,
    formula = frml
)

summary(pois_dose2)$coefficients

x<-as.data.frame(round(exp(pois_dose2$coefficients), digits=3))
x<-cbind(x,round(exp(pois_dose2$coefficients-coef(summary(pois_dose2))[,2]*1.96),3))
x<-cbind(x,round(exp(pois_dose2$coefficients+coef(summary(pois_dose2))[,2]*1.96),3))
x<-cbind(x, paste(x[,1]," (",x[,2],"-",x[,3],")", sep = ""))
colnames(x)[1]<-"Adjusted Rate Ratio"
colnames(x)[2]<-"LCI"
colnames(x)[3]<-"UCI"
colnames(x)[4]<-"ARR (LCI-UCI)"

z<-as.data.frame(pois_dose2$coefficients)
z<-cbind(z,round((pois_dose2$coefficients-coef(summary(pois_dose2))[,2]*1.96),3))
z<-cbind(z,round((pois_dose2$coefficients+coef(summary(pois_dose2))[,2])*1.96,3))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))

colnames(z)[1]<-"coefficients"
colnames(z)[2]<-"LCI - coef"
colnames(z)[3]<-"UCI - coef"
colnames(z)[4]<-"coef (LCI coef - UCI - coef)"

write.csv(x, file="output/national_study/Outputs for approval/Dose2/t_pois_irr_full.csv")
write.csv(z, file="output/national_study/Outputs for approval/Dose2/t_pois_coef_full.csv")

# Fit Dose2-AZ specific Poisson model ==========================================
cat("Fit Dose2-AZ specific Poisson model\n")

pois_dose2_az <- glm(
    data    = d_pyears,
    subset  = vacc_dose2_name == "AZ",
    family  = poisson,
    formula = frml
)

AZ<-as.data.frame(round(exp(pois_dose2_az$coefficients), digits=3))
AZ<-cbind(AZ,round(exp(pois_dose2_az$coefficients-coef(summary(pois_dose2_az))[,2]*1.96),3))
AZ<-cbind(AZ,round(exp(pois_dose2_az$coefficients+coef(summary(pois_dose2_az))[,2]*1.96),3))
AZ<-cbind(AZ, paste(AZ[,1]," (",AZ[,2],"-",AZ[,3],")", sep = ""))
colnames(AZ)[1]<-"Adjusted Rate Ratio"
colnames(AZ)[2]<-"LCI"
colnames(AZ)[3]<-"UCI"
colnames(AZ)[4]<-"ARR (LCI-UCI)"

AZc<-as.data.frame(pois_dose2_az$coefficients)
AZc<-cbind(AZc,round((pois_dose2_az$coefficients-coef(summary(pois_dose2_az))[,2]*1.96),3))
AZc<-cbind(AZc,round((pois_dose2_az$coefficients+coef(summary(pois_dose2_az))[,2])*1.96,3))
AZc<-cbind(AZc, paste(AZc[,1]," (",AZc[,2],"-",AZc[,3],")", sep = ""))

colnames(AZc)[1]<-"coefficients"
colnames(AZc)[2]<-"LCI - coef"
colnames(AZc)[3]<-"UCI - coef"
colnames(AZc)[4]<-"coef (LCI coef - UCI - coef)"

write.csv(AZ, file="output/national_study/Outputs for approval/Dose2/t_pois_irr_AZ.csv")
write.csv(AZc, file="output/national_study/Outputs for approval/Dose2/t_pois_coef_AZ.csv")

# Fit Dose2-PB specific Poisson model ==========================================
cat("Fit Dose2-PB specific Poisson model\n")

pois_dose2_pb <- glm(
    data    = d_pyears,
    subset  = vacc_dose2_name == "PB",
    family  = poisson,
    formula = frml
)

PB<-as.data.frame(round(exp(pois_dose2_pb$coefficients), digits=3))
PB<-cbind(PB,round(exp(pois_dose2_pb$coefficients-coef(summary(pois_dose2_pb))[,2]*1.96),3))
PB<-cbind(PB,round(exp(pois_dose2_pb$coefficients+coef(summary(pois_dose2_pb))[,2]*1.96),3))
PB<-cbind(PB, paste(PB[,1]," (",PB[,2],"-",PB[,3],")", sep = ""))
colnames(PB)[1]<-"Adjusted Rate Ratio"
colnames(PB)[2]<-"LCI"
colnames(PB)[3]<-"UCI"
colnames(PB)[4]<-"ARR (LCI-UCI)"

PBc<-as.data.frame(pois_dose2_pb$coefficients)
PBc<-cbind(PBc,round((pois_dose2_pb$coefficients-coef(summary(pois_dose2_pb))[,2]*1.96),3))
PBc<-cbind(PBc,round((pois_dose2_pb$coefficients+coef(summary(pois_dose2_pb))[,2])*1.96,3))
PBc<-cbind(PBc, paste(PBc[,1]," (",PBc[,2],"-",PBc[,3],")", sep = ""))

colnames(PBc)[1]<-"coefficients"
colnames(PBc)[2]<-"LCI - coef"
colnames(PBc)[3]<-"UCI - coef"
colnames(PBc)[4]<-"coef (LCI coef - UCI - coef)"

write.csv(PB, file="output/national_study/Outputs for approval/Dose2/t_pois_irr_PB.csv")
write.csv(PBc, file="output/national_study/Outputs for approval/Dose2/t_pois_coef_PB.csv")


# # Extract coef =================================================================
# cat("Extract coef\n")
# 
# expr_term <- attributes(terms(pois_dose2))$term.labels
# expr_term <- stringr::str_c(expr_term, collapse = "|")
# expr_term <- stringr::str_c("(", expr_term, ")(.+)")
# 
# make_tbl_coef <- function(model) {
#     dplyr::bind_cols(
#         tidy(model, exponentiate = TRUE),
#         # Wald confidence intervals
#         exp(confint_tidy(model, func = stats::confint.default))
#     ) %>%
#     mutate(
#         xvar = stringr::str_replace(term, expr_term, "\\1"),
#         xlbl = stringr::str_replace(term, expr_term, "\\2"),
#         .after = term
#     )%>%
#     dplyr::select(-term)
# }
# 
# t_coef_dose2    <- make_tbl_coef(pois_dose2)
# t_coef_dose2_az <- make_tbl_coef(pois_dose2_az)
# t_coef_dose2_pb <- make_tbl_coef(pois_dose2_pb)

# Save =========================================================================
cat("Save\n")

qsave(t_pyears_vacc_variant, file = "results/t_dose2_pyears_vacc_variant.qs")

qsave(pois_dose2,    file = "results/pois_dose2.qs")
qsave(pois_dose2_az, file = "results/pois_dose2_az.qs")
qsave(pois_dose2_pb, file = "results/pois_dose2_pb.qs")

# qsave(t_coef_dose2,    file = "results/t_coef_dose2.qs")
# qsave(t_coef_dose2_az, file = "results/t_coef_dose2_az.qs")
# qsave(t_coef_dose2_pb, file = "results/t_coef_dose2_pb.qs")

# # Print tidy coef table ========================================================
# cat("Print tidy coef table\n")
# 
# tidy_coef <- . %>%
#     mutate(
#         coef = stringr::str_glue("{est} ({low}, {high})",
#             est  = base::format(round(estimate,  2), nsmall = 2),
#             low  = base::format(round(conf.low,  2), nsmall = 2),
#             high = base::format(round(conf.high, 2), nsmall = 2)
#         )
#     ) %>%
#     dplyr::select(xvar, xlbl, coef)
# 
# table_lookup = tidy_coef(t_coef_dose2) %>%
#     dplyr::select(xvar, xlbl)
# 
# t_coef <- dplyr::bind_cols(
#    tidy_coef(t_coef_dose2)    %>% dplyr::select(xvar, xlbl, dose2 = coef),
#    tidy_coef(t_coef_dose2_az) %>% 
#        right_join(table_lookup, by=c("xvar","xlbl"), all.y = TRUE) 
#         %>% dplyr::select(dose2_az = coef),
#    tidy_coef(t_coef_dose2_pb) %>% dplyr::select(dose2_pb = coef)
# )
# 
# t_coef %>%
#     knitr::kable() %>%
#     kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
# print()

#beep()


