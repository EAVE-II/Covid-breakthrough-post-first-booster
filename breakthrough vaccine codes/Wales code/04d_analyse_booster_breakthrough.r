source("r_clear_and_load.r")

# Load =========================================================================
cat("Load\n")

d_analysis <- qread(s_drive("d_analysis.qs"))
d_bg       <- qread(s_drive("d_background.qs"))


# Analyse only booster people during Omicron period ============================
cat("Analyse only booster people during Omicron period\n")

d_analysis <- d_analysis %>% filter(in_study_20dec21)


# Remake analysis data set =====================================================
cat("Remake analysis data set\n")

# Omicron background intervals
bg_intervals <- d_bg %>%
    filter(bg_start >= ymd("2021-12-20")) %>%
    select(bg_interval) %>%
    unlist()

d_analysis <-
    d_analysis %>%
    # select analysis variables
    select(
        alf_e,
        # survival analysis things
        tstart,
        tstop,
        event_flg,
        # vaccine
        vacc_doseb_name,
        vacc_dose_interval,
        # covariates
        sex,
        age_3cat,
        age_4cat,
        age_5y_cat,
        ethn_cat,
        qcovid_cat,
        wimd2019_quintile,
        urban_rural_class,
        health_board,
        test_pre_dose2_cat,
        dose2_prior_infection_cat,
        vacc_dose1_dose2_diff_cat,
        bmi_cat,
        # background
        bg_interval,
        bg_avg_infection,
        bg_variant,
        # study design
        cc_group,
        sample_weight
    ) %>%
    # get only booster dose intervals over Omicron background intervals
    filter(
        str_detect(vacc_dose_interval, "doseb"),
        bg_interval %in% bg_intervals
    ) %>%
    # reset the follow-up clock to zero at 20th Decemeber
    lazy_dt() %>%
    group_by(alf_e) %>%
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
            str_replace("(doseb)_[A-Z]+_(day[0-9]+)", "\\1_\\2") %>%
            factor(),
        # drop any unused levels
        vacc_dose1_dose2_diff_cat = fct_drop(vacc_dose1_dose2_diff_cat)
    )


# Describe =====================================================================
cat("Describe\n")

calc_rate <- function(x) {
    control_weight <- d_analysis %>%
        filter(cc_group == "control") %>%
        select(sample_weight) %>%
        distinct() %>%
        unlist()

    case_weight <- d_analysis %>%
        filter(cc_group == "case") %>%
        select(sample_weight) %>%
        distinct() %>%
        unlist()

    d_analysis %>%
    lazy_dt() %>%
    group_by(vacc_doseb_name, !!sym(x)) %>%
    summarise(
        case_n    = n_distinct(ifelse(cc_group == "case",    alf_e, NA_real_), na.rm = TRUE),
        control_n = n_distinct(ifelse(cc_group == "control", alf_e, NA_real_), na.rm = TRUE),
        pyears    = sum((tstop - tstart) * sample_weight) / 365.25,
        event     = sum(event_flg)
    ) %>%
    ungroup() %>%
    as_tibble() %>%
    mutate(
        xvar = x,
        n = case_n*case_weight + control_n * control_weight,
        .before = case_n
    ) %>%
    select(-case_n, -control_n) %>%
    rename(xlbl = !!sym(x)) %>%
    select(xvar, xlbl, vacc_doseb_name, everything()) %>%
    arrange(xvar, xlbl, vacc_doseb_name) %>%
    mutate(rate = event / pyears * 1000)
}

x_vars <- c(
    "bg_variant",
    "vacc_dose_interval",
    "sex",
    "age_4cat",
    "qcovid_cat",
    "ethn_cat",
    "wimd2019_quintile",
    "urban_rural_class",
    "test_pre_dose2_cat",
    "dose2_prior_infection_cat",
    "vacc_dose1_dose2_diff_cat",
    "bmi_cat",
    "health_board"
)

t_desc_doseb_omicron_rate <- lapply(x_vars, calc_rate)
t_desc_doseb_omicron_rate <- bind_rows(t_desc_doseb_omicron_rate)


# Save descriptives ============================================================
cat("Save descriptives\n")

qsave(t_desc_doseb_omicron_rate, file = "results/t_desc_doseb_omicron_rate.qs")

t_desc_doseb_omicron_rate %>%
filter(xlbl != 0) %>%
kable(
    format.args = list(big.mark = ","),
    digits = c(0,0,0,0,1,0,1)
) %>%
kable_styling(
    bootstrap_options = "striped",
    full_width = FALSE
) %>%
print()


# Set reference categories =====================================================
cat("Set reference categories\n")

d_analysis <-
    d_analysis %>%
    mutate(
        bmi_cat = fct_relevel(bmi_cat, "18.5-24.9"),
        vacc_dose_interval = fct_relevel(vacc_dose_interval, "doseb_day014"),
        dose2_prior_infection_cat = fct_relevel(dose2_prior_infection_cat, "No prior infection"),
        wimd2019_quintile = fct_relevel(wimd2019_quintile, "5least")
    )


# Aggregate to person-years ====================================================
cat("Aggregate to person-years\n")

d_pyears <-
    d_analysis %>%
    calc_pyears(
        vacc_doseb_name,
        vacc_dose_interval,
        sex,
        age_5y_cat, # age_4cat age_5y_cat
        qcovid_cat,
        ethn_cat,
        wimd2019_quintile,
        urban_rural_class,
        test_pre_dose2_cat,
        dose2_prior_infection_cat,
        vacc_dose1_dose2_diff_cat,
        bmi_cat,
        health_board,
        bg_interval,
        bg_avg_infection,
        bg_variant
    )


# Fit overall Poisson model ====================================================
cat("Fit overall booster Poisson models\n")

frml_bgavg <- event ~ vacc_dose_interval +
                      sex +
                      age_5y_cat +
                      qcovid_cat +
                      ethn_cat +
                      wimd2019_quintile +
                      urban_rural_class +
                      test_pre_dose2_cat +
                      dose2_prior_infection_cat +
                      vacc_dose1_dose2_diff_cat +
                      bmi_cat +
                      health_board +
                      bg_avg_infection +
                      offset(log(pyears))

frml <- event ~ vacc_dose_interval +
                     sex +
                     age_5y_cat +
                     qcovid_cat +
                     ethn_cat +
                     wimd2019_quintile +
                     urban_rural_class +
                     test_pre_dose2_cat +
                     dose2_prior_infection_cat +
                     vacc_dose1_dose2_diff_cat +
                     bmi_cat +
                     health_board +
                     bg_interval +
                     offset(log(pyears))

pois_doseb_bgavg <- glm(
    data    = d_pyears,
    family  = poisson,
    formula = frml_bgavg
)

pois_doseb <- glm(
    data    = d_pyears,
    family  = poisson,
    formula = frml
)


# Fit overall booster Cox PH model =============================================
cat("Fit overall booster Cox PH model\n")

frml_cph <- Surv(tstart, tstop, event_flg) ~ vacc_dose_interval +
                                             sex +
                                             age_5y_cat +
                                             qcovid_cat +
                                             ethn_cat +
                                             wimd2019_quintile +
                                             urban_rural_class +
                                             test_pre_dose2_cat +
                                             dose2_prior_infection_cat +
                                             vacc_dose1_dose2_diff_cat +
                                             bmi_cat +
                                             strata(health_board)

cph_doseb <- coxph(
    data    = d_analysis,
    id      = alf_e,
    formula = frml_cph,
    weights = sample_weight
)

# quick comparison of estimates ------------------------------------------------
cat("\tquick comparison of estimates\n")

too_big <- function(x) {ifelse(x > 10^6, Inf, x)}

t_pois_bgavg_coef <-
    bind_cols(
        tidy(pois_doseb_bgavg, exponentiate = TRUE),
        exp(confint_tidy(pois_doseb_bgavg, func = stats::confint.default))
    ) %>%
    mutate(
        pois_bgavg_coef = str_glue("{est} ({conf_low}, {conf_high})",
            est       = format(round(too_big(estimate), 2), nsmall = 2),
            conf_low  = format(round(too_big(conf.low), 2), nsmall = 2),
            conf_high = format(round(too_big(conf.high), 2), nsmall = 2),
        )
    ) %>%
    select(term, pois_bgavg_coef)

t_pois_coef <-
    bind_cols(
        tidy(pois_doseb, exponentiate = TRUE),
        exp(confint_tidy(pois_doseb, func = stats::confint.default))
    ) %>%
    mutate(
        pois_coef = str_glue("{est} ({conf_low}, {conf_high})",
            est       = format(round(too_big(estimate), 2), nsmall = 2),
            conf_low  = format(round(too_big(conf.low), 2), nsmall = 2),
            conf_high = format(round(too_big(conf.high), 2), nsmall = 2),
        )
    ) %>%
    select(term, pois_coef)

t_cph_coef <-
    tidy(cph_doseb, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
        cph_coef = str_glue("{est} ({conf_low}, {conf_high})",
            est       = format(round(estimate, 2), nsmall = 2),
            conf_low  = format(round(conf.low, 2), nsmall = 2),
            conf_high = format(round(conf.high, 2), nsmall = 2),
        )
    ) %>%
    select(term, cph_coef)

t_pois_bgavg_coef %>%
    full_join(t_pois_coef, by = "term") %>%
    full_join(t_cph_coef, by = "term") %>%
    kable() %>%
    kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
    print()

# Fit booster-MD specific Poisson model ==========================================
cat("Fit booster-MD specific Poisson models\n")

pois_doseb_md <- glm(
    data    = d_pyears,
    subset  = vacc_doseb_name == "MD",
    family  = poisson,
    formula = frml
)


# Fit booster-PB specific Poisson model ==========================================
cat("Fit booster-PB specific Poisson models\n")

pois_doseb_pb <- glm(
    data    = d_pyears,
    subset  = vacc_doseb_name == "PB",
    family  = poisson,
    formula = frml
)


# Extract coef =================================================================
cat("Extract coef\n")

expr_term <- attributes(terms(pois_doseb))$term.labels
expr_term <- str_c(expr_term, collapse = "|")
expr_term <- str_c("(", expr_term, ")(.+)")

make_tbl_coef <- function(model) {
    bind_cols(
        tidy(model, exponentiate = TRUE),
        # Wald confidence intervals
        exp(confint.default(model))
    ) %>%
    clean_names() %>%
    rename(
        conf_low = x2_5_percent,
        conf_high = x97_5_percent
    ) %>%
    mutate(
        xvar = str_replace(term, expr_term, "\\1"),
        xlbl = str_replace(term, expr_term, "\\2"),
        .after = term
    )%>%
    select(-term)
}

t_coef_doseb    <- make_tbl_coef(pois_doseb)
t_coef_doseb_md <- make_tbl_coef(pois_doseb_md)
t_coef_doseb_pb <- make_tbl_coef(pois_doseb_pb)


# Save =========================================================================
cat("Save\n")

qsave(d_analysis, file = s_drive("d_analysis_doseb.qs"))

qsavem(
    pois_doseb_bgavg,
    pois_doseb,
    cph_doseb,
    pois_doseb_md,
    pois_doseb_pb,
    file = "results/booster_breakthrough_models.qsm"
)

qsave(t_coef_doseb,    file = "results/t_coef_doseb.qs")
qsave(t_coef_doseb_md, file = "results/t_coef_doseb_md.qs")
qsave(t_coef_doseb_pb, file = "results/t_coef_doseb_pb.qs")


# Print tidy coef table ========================================================
cat("Print tidy coef table\n")

too_small <- function(x) {ifelse(x < 1/10^9, -Inf, x)}
too_big <- function(x) {ifelse(x > 10^9, Inf, x)}

tidy_coef <- . %>%
    mutate(
        coef = str_glue("{est} ({low}, {high})",
            est  = format(round(estimate,  2), nsmall = 2),
            low  = format(round(too_small(conf_low),  2), nsmall = 2),
            high = format(round(too_big(conf_high), 2), nsmall = 2)
        )
    ) %>%
    select(xvar, xlbl, coef)

t_coef <-
    tidy_coef(t_coef_doseb) %>% rename(dose_b = coef) %>%
    left_join(
        y = tidy_coef(t_coef_doseb_md) %>% rename(doseb_md = coef),
        by = c("xvar", "xlbl")
    ) %>%
    left_join(
        y = tidy_coef(t_coef_doseb_pb) %>% rename(doseb_pb = coef),
        by = c("xvar", "xlbl")
    )

t_coef %>%
kable() %>%
kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
print()
