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
        # qcovid
        qc_b2_82,
        qc_b2_leukolaba,
        qc_b2_prednisolone,
        qc_b_af,
        qc_b_ccf,
        qc_b_asthma,
        qc_b_bloodcancer,
        qc_b_cerebralpalsy,
        qc_b_chd,
        qc_b_cirrhosis,
        qc_b_congenheart,
        qc_b_copd,
        qc_b_dementia,
        qc_b_epilepsy,
        qc_b_fracture4,
        qc_b_neurorare,
        qc_b_parkinsons,
        qc_b_pulmhyper,
        qc_b_pulmrare,
        qc_b_pvd,
        qc_b_ra_sle,
        qc_b_respcancer,
        qc_b_semi,
        qc_b_sicklecelldisease,
        qc_b_stroke,
        qc_diabetes1,
        qc_diabetes2,
        qc_b_vte,
        qc_chemo_cat,
        qc_carehome,
        qc_homeless,
        qc_learndisab,
        qc_downs,
        #qc_p_marrow6,          # <!> omitted because no one has this condition
        qc_p_radio6,
        qc_p_solidtransplant,
        qc_ckd3,
        qc_ckd4,
        qc_ckd5,
        renal_transplant,
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
    # reset the follow-up clock to zero
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

# List QCovid variables ========================================================
cat("List QCovid variables\n")

x_qcovid <- c(
    "qc_b2_82",
    "qc_b2_leukolaba",
    "qc_b2_prednisolone",
    "qc_b_af",
    "qc_b_ccf",
    "qc_b_asthma",
    "qc_b_bloodcancer",
    "qc_b_cerebralpalsy",
    "qc_b_chd",
    "qc_b_cirrhosis",
    "qc_b_congenheart",
    "qc_b_copd",
    "qc_b_dementia",
    "qc_b_epilepsy",
    "qc_b_fracture4",
    "qc_b_neurorare",
    "qc_b_parkinsons",
    "qc_b_pulmhyper",
    "qc_b_pulmrare",
    "qc_b_pvd",
    "qc_b_ra_sle",
    "qc_b_respcancer",
    "qc_b_semi",
    "qc_b_sicklecelldisease",
    "qc_b_stroke",
    "qc_diabetes1",
    "qc_diabetes2",
    "qc_b_vte",
    "qc_chemo_cat",
    "qc_carehome",
    "qc_homeless",
    "qc_learndisab",
    "qc_downs",
    "qc_p_radio6",
    "qc_p_solidtransplant",
    "qc_ckd3",
    "qc_ckd4",
    "qc_ckd5",
    "renal_transplant"
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
    group_by(!!sym(x)) %>%
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
    select(xvar, xlbl, everything()) %>%
    arrange(xvar, xlbl) %>%
    mutate(rate = event / pyears * 1000)
}

t_desc_doseb_qcovid_rate <- lapply(x_qcovid, calc_rate)
t_desc_doseb_qcovid_rate <- bind_rows(t_desc_doseb_qcovid_rate)


# Save descriptives ============================================================
cat("Save descriptives\n")

qsave(t_desc_doseb_qcovid_rate, file = "results/t_desc_doseb_qcovid_rate.qs")

t_desc_doseb_qcovid_rate %>%
filter(xlbl != 0) %>%
arrange(desc(rate)) %>%
kable(
    format.args = list(big.mark = ","),
    digits = c(0,0, 0,1,0,1)
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

# Estimate QCovid rate ratios ==================================================
cat("Estimate QCovid rate ratios\n")

get_rate_ratios <- function(xvar) {
    # xvar = "qc_b2_82"

    # (1) calc pyears using saturated forumla plus xvar
    # (2) fit poisson model
    # (3) extract estimate for qcovid measure with 95% CI

    cat(xvar)

    # vacc_doseb_name is needed below if we do subset modelling
    # can use any of the following age categories: age_3cat age_4cat age_5y_cat
    d_pyears <- d_analysis %>%
        calc_pyears(
            vacc_dose_interval,
            sex,
            age_5y_cat,
            wimd2019_quintile,
            urban_rural_class,
            health_board,
            test_pre_dose2_cat,
            dose2_prior_infection_cat,
            vacc_dose1_dose2_diff_cat,
            bmi_cat,
            health_board,
            bg_interval,
            !!sym(xvar)
        )

    cat(".")

    frml_pois <- event ~ vacc_dose_interval +
                         sex +
                         age_5y_cat +
                         wimd2019_quintile +
                         urban_rural_class +
                         health_board +
                         test_pre_dose2_cat +
                         dose2_prior_infection_cat +
                         vacc_dose1_dose2_diff_cat +
                         bmi_cat +
                         health_board +
                         bg_interval +
                         offset(log(pyears))

    frml_pois <- update(frml_pois, str_c(". ~ . + ", xvar))

    pois_doseb <- glm(
        data    = d_pyears,
        family  = poisson,
        formula = frml_pois
    )

    cat(".")

    expr_term <- attributes(terms(pois_doseb))$term.labels
    expr_term <- str_c(expr_term, collapse = "|")
    expr_term <- str_c("(", expr_term, ")(.+)")

    t_coef <-
        bind_cols(
            tidy(pois_doseb, exponentiate = TRUE),
            # Wald confidence intervals
            exp(confint.default(pois_doseb))
        ) %>%
        clean_names() %>%
        rename(
            conf_low = x2_5_percent,
            conf_high = x97_5_percent
        ) %>%
        filter(str_detect(term, xvar))

    cat(".")
    cat("\n")

    return(t_coef)
}

t_doseb_qcovid <- lapply(x_qcovid, get_rate_ratios)
t_doseb_qcovid <- bind_rows(t_doseb_qcovid)

# Save =========================================================================
cat("Save\n")

qsave(
    t_doseb_qcovid,
    file = "results/t_doseb_qcovid_rr.qs"
)

# Print coef ===================================================================
cat("Print coef\n")

t_doseb_qcovid %>%
select(term, estimate, conf_low, conf_high) %>%
mutate(conf_high = ifelse(conf_high > 10^6, Inf, conf_high)) %>%
arrange(desc(estimate)) %>%
kable(
    format.args = list(big.mark = ","),
    digits = c(0,2,2,2)
) %>%
kable_styling(
    bootstrap_options = "striped",
    full_width = FALSE
) %>%
print()
