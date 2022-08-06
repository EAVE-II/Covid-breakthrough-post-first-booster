cat("Clearing workspace and loading stuff\n")

source("scripts/000_libraries_and_functions.R")

gc()
memory.limit(64000)
options(scipen = 100)

omicron_start_date <- ymd("2021-12-20")

# Load =========================================================================
cat("Load\n")

d_sample   <- qread("results/d_sample_b.qs")
d_analysis <- qread("results/d_analysis_stats_b.qs")
d_bg       <- qread("input/d_background.qs")


# Filter booster people only during Omicron period ============================
cat("Filter booster people only during Omicron period\n")

alf_w_booster <- d_sample %>%
    # has booster
    filter(
        assertr::not_na(vacc_doseb_date),
        assertr::not_na(vacc_doseb_name)
    ) %>%
    # person is still in the study as of 20th Decemeber 2021
    filter(
        (omicron_start_date < outcome_date        | is.na(outcome_date)),
        (omicron_start_date < death_noncovid_date | is.na(death_noncovid_date)),
        (omicron_start_date < Dose4_date          | is.na(Dose4_date)),
        (omicron_start_date < study_end_date      | is.na(study_end_date))
    ) %>%
    # person has not had an outcome prior to booster
    filter(
        (vacc_doseb_date < outcome_date        | is.na(outcome_date)),
        (vacc_doseb_date < death_noncovid_date | is.na(death_noncovid_date)),
        (vacc_doseb_date < Dose4_date          | is.na(Dose4_date)),
        (vacc_doseb_date < study_end_date      | is.na(study_end_date))
    ) %>%
    # keep only id column
    select(study_id) %>% collect()

msg_n1 <- format(nrow(d_sample) - nrow(alf_w_booster), big.mark = ",")
msg_n2 <- format(nrow(alf_w_booster), big.mark = ",")
cat("    Dropped ", msg_n1, " people\n    Keeping ", msg_n2, " people\n", sep = "")

d_analysis <- d_analysis %>% inner_join( alf_w_booster, by = "study_id" )

# Remake analysis data set =====================================================
cat("Remake analysis data set\n")

# omicron background intervals
omicron_bg_intervals <- d_bg %>%
    filter(bg_start >= omicron_start_date) %>%
    select(bg_interval) %>%
    unlist()

d_analysis <-
    d_analysis %>%
    # select analysis variables
    select(
        study_id,
        # survival analysis things
        tstart,
        tstop,
        event_flg,
        # vaccine
        vacc_doseb_name,
        vacc_dose_interval,
        # covariates
        sex,
        age_gp,
        BNF_group,
        nimdm,
        ur_combined,
        lgd,
        n_tests_gp,
        dose2_prior_infection_cat,
        vacc_dose1_dose2_diff_cat,
        # background
        bg_interval,
        bg_avg_infection,
        bg_variant,
        # study design
        sample_weight
    ) %>%
    # get only booster dose intervals over omicron background intervals
    filter(
        stringr::str_detect(vacc_dose_interval, "doseb"),
        bg_interval %in% omicron_bg_intervals
    ) %>%
    # reset the follow-up clock to zero at 20th Decemeber
    dtplyr::lazy_dt() %>%
    group_by(study_id) %>%
    mutate(
        tstop  = tstop  - min(tstart),
        tstart = tstart - min(tstart)
    ) %>%
    ungroup() %>%
    as_tibble() %>%
    # final edits
    mutate(
        bg_interval = factor(bg_interval),
        # remove vaccine name from booster interval
        vacc_dose_interval = vacc_dose_interval %>%
            as.character() %>%
            stringr::str_replace("(doseb)_[A-Z]+_(day[0-9]+)", "\\1_\\2") %>%
            factor(),
        # drop any unused levels
        vacc_dose1_dose2_diff_cat = forcats::fct_drop(vacc_dose1_dose2_diff_cat),
        # set reference categories
        sex = forcats::fct_relevel(sex, "Female"),
        age_gp = forcats::fct_relevel(age_gp, "18-49"),
        vacc_dose_interval = forcats::fct_relevel(vacc_dose_interval, "doseb_day014"),
        dose2_prior_infection_cat = forcats::fct_relevel(dose2_prior_infection_cat, "No prior infection"),
        nimdm = forcats::fct_relevel(nimdm, "5 - least deprived")
    )


# Aggregate to person-years ====================================================
cat("Aggregate to person-years\n")

d_pyears <-
    d_analysis %>%
    calc_pyears(
        vacc_doseb_name,
        vacc_dose_interval,
        sex,
        age_gp, # age_4cat age_5y_cat
        BNF_group,
        nimdm,
        ur_combined,
        n_tests_gp,
        dose2_prior_infection_cat,
        vacc_dose1_dose2_diff_cat,
        lgd,
        bg_interval,
        bg_avg_infection,
        bg_variant
    )

# Fit overall Poisson model ====================================================
cat("Fit overall booster Poisson models\n")

# frml_bgavg <- event ~ vacc_dose_interval +
#                 dose2_prior_infection_cat +
#                 vacc_dose1_dose2_diff_cat +              
#                 sex +
#                 age_gp +
#                 nimdm +              
#                 ur_combined +              
#                 BNF_group +
#                 n_tests_gp +
#                 lgd +
#                 bg_avg_infection +
#                 offset(log(pyears))

frml_full <- event ~ vacc_dose_interval +
                  dose2_prior_infection_cat +
                  vacc_dose1_dose2_diff_cat +
                  sex +
                  age_gp +
                  nimdm +
                  ur_combined +
                  BNF_group +
                  n_tests_gp +
                  lgd +
                  bg_interval +
                  offset(log(pyears))

# pois_doseb_bgavg <- glm(
#     data    = d_pyears,
#     family  = poisson,
#     formula = frml_bgavg
# )

pois_doseb <- glm(
    data    = d_pyears,
    family  = poisson,
    formula = frml_full
)

## Write out coefficients ===========

summary(pois_doseb)$coefficients

x<-as.data.frame(round(exp(pois_doseb$coefficients), digits=3))
x<-cbind(x,round(exp(pois_doseb$coefficients-coef(summary(pois_doseb))[,2]*1.96),3))
x<-cbind(x,round(exp(pois_doseb$coefficients+coef(summary(pois_doseb))[,2]*1.96),3))
x<-cbind(x, paste(x[,1]," (",x[,2],"-",x[,3],")", sep = ""))
colnames(x)[1]<-"Adjusted Rate Ratio"
colnames(x)[2]<-"LCI"
colnames(x)[3]<-"UCI"
colnames(x)[4]<-"ARR (LCI-UCI)"

z<-as.data.frame(pois_doseb$coefficients)
z<-cbind(z,round((pois_doseb$coefficients-coef(summary(pois_doseb))[,2]*1.96),3))
z<-cbind(z,round((pois_doseb$coefficients+coef(summary(pois_doseb))[,2])*1.96,3))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))

colnames(z)[1]<-"coefficients"
colnames(z)[2]<-"LCI - coef"
colnames(z)[3]<-"UCI - coef"
colnames(z)[4]<-"coef (LCI coef - UCI - coef)"

write.csv(x, file="output/national_study/Outputs for approval/Booster/t_pois_irr_b_all.csv")
write.csv(z, file="output/national_study/Outputs for approval/Booster/t_pois_coef_b_all.csv")

#pois_doseb_common <- glm(
 #   data    = d_pyears,
  #  family  = poisson,
   # formula = frml_common
#)


# Fit overall booster Cox PH model =============================================
#cat("Fit overall booster Cox PH model\n")

#frml_cph <- Surv(tstart, tstop, event_flg) ~ vacc_dose_interval +
                          #                   sex +
                           #                  age_5y_cat +
                            #                 qcovid_cat +
                             #                ethn_cat +
                              #               wimd2019_quintile +
                               #              urban_rural_class +
                                #             test_pre_dose2_cat +
                                 #            dose2_prior_infection_cat +
                                  #           vacc_dose1_dose2_diff_cat +
                                   #          bmi_cat +
                                    #         strata(health_board)

#cph_doseb_full <- coxph(
 #   data    = d_analysis,
  #  id      = alf_e,
   # formula = frml_cph,
    #weights = sample_weight
#)

# quick comparison of estimates ------------------------------------------------
#cat("\tquick comparison of estimates\n")

#too_big <- function(x) {ifelse(x > 10^6, Inf, x)}

#t_pois_bgavg_coef <-
 #   bind_cols(
  #      tidy(pois_doseb_bgavg, exponentiate = TRUE),
   #     exp(confint_tidy(pois_doseb_bgavg, func = stats::confint.default))
    #) %>%
    #mutate(
     #   pois_bgavg_coef = str_glue("{est} ({conf_low}, {conf_high})",
      #      est       = format(round(too_big(estimate), 2), nsmall = 2),
       #     conf_low  = format(round(too_big(conf.low), 2), nsmall = 2),
        #    conf_high = format(round(too_big(conf.high), 2), nsmall = 2),
        #)
    #) %>%
    #select(term, pois_bgavg_coef)

#t_pois_full_coef <-
 #   bind_cols(
  #      tidy(pois_doseb_full, exponentiate = TRUE),
   #     exp(confint_tidy(pois_doseb_full, func = stats::confint.default))
    #) %>%
    #mutate(
     #   pois_full_coef = str_glue("{est} ({conf_low}, {conf_high})",
      #      est       = format(round(too_big(estimate), 2), nsmall = 2),
       #     conf_low  = format(round(too_big(conf.low), 2), nsmall = 2),
        #    conf_high = format(round(too_big(conf.high), 2), nsmall = 2),
        #)
    #) %>%
    #select(term, pois_full_coef)

#t_cph_coef <-
 #   tidy(cph_doseb_full, exponentiate = TRUE, conf.int = TRUE) %>%
  #  mutate(
   #     cph_coef = str_glue("{est} ({conf_low}, {conf_high})",
    #        est       = format(round(estimate, 2), nsmall = 2),
     #       conf_low  = format(round(conf.low, 2), nsmall = 2),
      #      conf_high = format(round(conf.high, 2), nsmall = 2),
       # )
    #) %>%
    #select(term, cph_coef)

#t_pois_bgavg_coef %>%
 #   full_join(t_pois_full_coef, by = "term") %>%
  #  full_join(t_cph_coef, by = "term") %>%
   # kable() %>%
    #kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
    #print()

# Fit booster-MD specific Poisson model ==========================================
cat("Fit booster-MD specific Poisson models\n")

pois_doseb_md_full <- glm(
    data    = d_pyears,
    subset  = vacc_doseb_name == "MD",
    family  = poisson,
    formula = frml_full
)

#pois_doseb_md_common <- glm(
 #   data    = d_pyears,
  #  subset  = vacc_doseb_name == "MD",
   # family  = poisson,
  #  formula = frml_common
#)

# Write out coefficients

MD<-as.data.frame(round(exp(pois_doseb_md_full$coefficients), digits=3))
MD<-cbind(MD,round(exp(pois_doseb_md_full$coefficients-coef(summary(pois_doseb_md_full))[,2]*1.96),3))
MD<-cbind(MD,round(exp(pois_doseb_md_full$coefficients+coef(summary(pois_doseb_md_full))[,2]*1.96),3))
MD<-cbind(MD, paste(MD[,1]," (",MD[,2],"-",MD[,3],")", sep = ""))
colnames(MD)[1]<-"Adjusted Rate Ratio"
colnames(MD)[2]<-"LCI"
colnames(MD)[3]<-"UCI"
colnames(MD)[4]<-"ARR (LCI-UCI)"

MDc<-as.data.frame(pois_doseb_md_full$coefficients)
MDc<-cbind(MDc,round((pois_doseb_md_full$coefficients-coef(summary(pois_doseb_md_full))[,2]*1.96),3))
MDc<-cbind(MDc,round((pois_doseb_md_full$coefficients+coef(summary(pois_doseb_md_full))[,2])*1.96,3))
MDc<-cbind(MDc, paste(MDc[,1]," (",MDc[,2],"-",MDc[,3],")", sep = ""))

colnames(MDc)[1]<-"coefficients"
colnames(MDc)[2]<-"LCI - coef"
colnames(MDc)[3]<-"UCI - coef"
colnames(MDc)[4]<-"coef (LCI coef - UCI - coef)"

write.csv(MD, file="output/national_study/Outputs for approval/Booster/t_pois_irr_b_MD.csv")
write.csv(MDc, file="output/national_study/Outputs for approval/Booster/t_pois_coef_b_MD.csv")



# Fit booster-PB specific Poisson model ==========================================
cat("Fit booster-PB specific Poisson models\n")

pois_doseb_pb_full <- glm(
    data    = d_pyears,
    subset  = vacc_doseb_name == "PB",
    family  = poisson,
    formula = frml_full
)

#pois_doseb_pb_common <- glm(
 #   data    = d_pyears,
  #  subset  = vacc_doseb_name == "PB",
   # family  = poisson,
    #formula = frml_common
#)

# Write out coefficients =======================================

PB<-as.data.frame(round(exp(pois_doseb_pb_full$coefficients), digits=3))
PB<-cbind(PB,round(exp(pois_doseb_pb_full$coefficients-coef(summary(pois_doseb_pb_full))[,2]*1.96),3))
PB<-cbind(PB,round(exp(pois_doseb_pb_full$coefficients+coef(summary(pois_doseb_pb_full))[,2]*1.96),3))
PB<-cbind(PB, paste(PB[,1]," (",PB[,2],"-",PB[,3],")", sep = ""))
colnames(PB)[1]<-"Adjusted Rate Ratio"
colnames(PB)[2]<-"LCI"
colnames(PB)[3]<-"UCI"
colnames(PB)[4]<-"ARR (LCI-UCI)"

PBc<-as.data.frame(pois_doseb_pb_full$coefficients)
PBc<-cbind(PBc,round((pois_doseb_pb_full$coefficients-coef(summary(pois_doseb_pb_full))[,2]*1.96),3))
PBc<-cbind(PBc,round((pois_doseb_pb_full$coefficients+coef(summary(pois_doseb_pb_full))[,2])*1.96,3))
PBc<-cbind(PBc, paste(PBc[,1]," (",PBc[,2],"-",PBc[,3],")", sep = ""))

colnames(PBc)[1]<-"coefficients"
colnames(PBc)[2]<-"LCI - coef"
colnames(PBc)[3]<-"UCI - coef"
colnames(PBc)[4]<-"coef (LCI coef - UCI - coef)"

write.csv(PB, file="output/national_study/Outputs for approval/Booster/t_pois_irr_b_PB.csv")
write.csv(PBc, file="output/national_study/Outputs for approval/Booster/t_pois_coef_b_PB.csv")


# # Extract coef =================================================================
# cat("Extract coef\n")
# 
# expr_term <- attributes(terms(pois_doseb_full))$term.labels
# expr_term <- stringr::str_c(expr_term, collapse = "|")
# expr_term <- stringr::str_c("(", expr_term, ")(.+)")
# 
# make_tbl_coef <- function(model) {
#     bind_cols(
#         broom::tidy(model, exponentiate = TRUE),
#         # Wald confidence intervals
#         exp(confint.default(model))
#     ) %>%
#     janitor::clean_names() %>%
#     rename(
#         conf_low = x2_5_percent,
#         conf_high = x97_5_percent
#     ) %>%
#     mutate(
#         xvar = stringr::str_replace(term, expr_term, "\\1"),
#         xlbl = stringr::str_replace(term, expr_term, "\\2"),
#         .after = term
#     )%>%
#     select(-term)
# }
# 
# t_coef_doseb    <- make_tbl_coef(pois_doseb_full)
# t_coef_doseb_md <- make_tbl_coef(pois_doseb_md_full)
# t_coef_doseb_pb <- make_tbl_coef(pois_doseb_pb_full)
# 
# 
# Save =========================================================================
cat("Save\n")


qsavem(
    #pois_doseb_bgavg,
    pois_doseb,
    #pois_doseb_common,
    #cph_doseb_full,
    pois_doseb_md_full,
    #pois_doseb_md_common,
    pois_doseb_pb_full,
    #pois_doseb_pb_common,
    file = "results/booster_breakthrough_models.qsm"
)

# qsave(t_coef_doseb,    file = "results/t_coef_doseb.qs")
# qsave(t_coef_doseb_md, file = "results/t_coef_doseb_md.qs")
# qsave(t_coef_doseb_pb, file = "results/t_coef_doseb_pb.qs")
# 
# 
# # Print tidy coef table ========================================================
# cat("Print tidy coef table\n")
# 
# too_small <- function(x) {ifelse(x < 1/10^9, -Inf, x)}
# too_big <- function(x) {ifelse(x > 10^9, Inf, x)}
# 
# tidy_coef <- . %>%
#     mutate(
#         coef = stringr::str_glue("{est} ({low}, {high})",
#             est  = format(round(estimate,  2), nsmall = 2),
#             low  = format(round(too_small(conf_low),  2), nsmall = 2),
#             high = format(round(too_big(conf_high), 2), nsmall = 2)
#         )
#     ) %>%
#     select(xvar, xlbl, coef)
# 
# t_coef <-
#     tidy_coef(t_coef_doseb) %>% rename(dose_b = coef) %>%
#     left_join(
#         y = tidy_coef(t_coef_doseb_md) %>% rename(doseb_md = coef),
#         by = c("xvar", "xlbl")
#     ) %>%
#     left_join(
#         y = tidy_coef(t_coef_doseb_pb) %>% rename(doseb_pb = coef),
#         by = c("xvar", "xlbl")
#     )
# 
# t_coef %>%
# knitr::kable() %>%
# kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
# print()


