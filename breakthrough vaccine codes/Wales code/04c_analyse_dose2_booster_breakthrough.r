source("r_clear_and_load.r")

# Load =========================================================================
cat("Load\n")

d_analysis <- qread(s_drive("d_analysis.qs"))
d_bg       <- qread(s_drive("d_background.qs"))


# Analyse from Monday 13th Sept 2021 ===========================================
cat("Analyse from Monday 13th Sept 2021\n")

d_analysis <- d_analysis %>% filter(in_study_13sep21)


# Remake analysis data set =====================================================
cat("Remake analysis data set\n")

# relevant background intervals
bg_intervals <- d_bg %>%
    filter(bg_start >= ymd("2021-09-13")) %>%
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
        vacc_dose2_name,
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
    # get relevant intervals
    filter(
        bg_interval %in% bg_intervals
    ) %>%
    # reset the follow-up clock to zero at epoch
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
        # remove vaccine name from dose2 interval
        vacc_dose_interval = vacc_dose_interval %>%
            as.character() %>%
            str_replace("(dose2)_[A-Z]+_(day[0-9]+)", "\\1_\\2") %>%
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
    group_by(vacc_dose2_name, !!sym(x)) %>%
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
    select(xvar, xlbl, vacc_dose2_name, everything()) %>%
    arrange(xvar, xlbl, vacc_dose2_name) %>%
    mutate(rate = event / pyears * 1000)
}

x_vars <- c(
    "bg_variant",
    "vacc_dose_interval",
    "sex",
    "age_4cat",
    "ethn_cat",
    "qcovid_cat",
    "bmi_cat",
    "wimd2019_quintile",
    "urban_rural_class",
    "test_pre_dose2_cat",
    "dose2_prior_infection_cat",
    "vacc_dose1_dose2_diff_cat",
    "health_board"
)

t_desc_dose2b_rate <- lapply(x_vars, calc_rate)
t_desc_dose2b_rate <- bind_rows(t_desc_dose2b_rate)

# Save descriptives ============================================================
cat("Save descriptives\n")

qsave(t_desc_dose2b_rate, file = "results/t_desc_dose2b_rate.qs")

t_desc_dose2b_rate %>%
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
        vacc_dose_interval = fct_relevel(vacc_dose_interval, "dose2_day014"),
        dose2_prior_infection_cat = fct_relevel(dose2_prior_infection_cat, "No prior infection"),
        wimd2019_quintile = fct_relevel(wimd2019_quintile, "5least")
    )


# Aggregate to person-years ====================================================
cat("Aggregate to person-years\n")

d_pyears <-
    d_analysis %>%
    calc_pyears(
        vacc_dose2_name,
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
cat("Fit overall Poisson model\n")

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

pois_dose2 <- glm(
    data    = d_pyears,
    family  = poisson,
    formula = frml
)


# Fit Dose2-AZ specific Poisson model ==========================================
cat("Fit Dose2-AZ specific Poisson model\n")

pois_dose2_az <- glm(
    data    = d_pyears,
    subset  = vacc_dose2_name == "AZ",
    family  = poisson,
    formula = frml
)


# Fit Dose2-PB specific Poisson model ==========================================
cat("Fit Dose2-PB specific Poisson model\n")

pois_dose2_pb <- glm(
    data    = d_pyears,
    subset  = vacc_dose2_name == "PB",
    family  = poisson,
    formula = frml
)


# Extract coef =================================================================
cat("Extract coef\n")

expr_term <- attributes(terms(pois_dose2))$term.labels
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

t_coef_dose2    <- make_tbl_coef(pois_dose2)
t_coef_dose2_az <- make_tbl_coef(pois_dose2_az)
t_coef_dose2_pb <- make_tbl_coef(pois_dose2_pb)


# Save models ==================================================================
cat("Save models\n")

qsave(d_analysis, file = s_drive("d_analysis_dose2.qs"))

qsave(pois_dose2,    file = "results/pois_dose2.qs")
qsave(pois_dose2_az, file = "results/pois_dose2_az.qs")
qsave(pois_dose2_pb, file = "results/pois_dose2_pb.qs")

qsave(t_coef_dose2,    file = "results/t_coef_dose2.qs")
qsave(t_coef_dose2_az, file = "results/t_coef_dose2_az.qs")
qsave(t_coef_dose2_pb, file = "results/t_coef_dose2_pb.qs")


# Print tidy coef table ========================================================
cat("Print tidy coef table\n")

tidy_coef <- . %>%
    mutate(
        coef = str_glue("{est} ({low}, {high})",
            est  = format(round(estimate,  2), nsmall = 2),
            low  = format(round(conf_low,  2), nsmall = 2),
            high = format(round(conf_high, 2), nsmall = 2)
        )
    ) %>%
    select(xvar, xlbl, coef)

t_coef <- bind_cols(
    tidy_coef(t_coef_dose2)    %>% select(xvar, xlbl, dose2 = coef),
    tidy_coef(t_coef_dose2_az) %>% select(dose2_az = coef),
    tidy_coef(t_coef_dose2_pb) %>% select(dose2_pb = coef)
)

t_coef %>%
kable() %>%
kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
print()

beep()
