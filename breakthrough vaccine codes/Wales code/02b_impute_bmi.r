source("r_clear_and_load.r")

# Load =========================================================================
cat("Load\n")

d_cohort <- qread(s_drive("d_cohort_clean.qs"))

# Impute BMI ===================================================================
cat("Impute BMI\n")

d_cohort <- d_cohort %>% mutate(qc_bmi_log = log(qc_bmi))

# create predictor matrix
pred.mat <- matrix(
    data = 0,
    nrow = ncol(d_cohort),
    ncol = ncol(d_cohort),
    dimnames = list(
        names(d_cohort),
        names(d_cohort)
    )
)

# set BMI to use the following measures
bmi_xvar <- c(
    # outcome
    "has_death_noncovid",
    "has_death_covid",
    "has_covid_hosp",
    "has_covid_infection",
    # exposure
    "has_vacc_dose1",
    "has_vacc_dose2",
    "has_vacc_doseb",
    "vacc_dose1_dose2_diff_cat",
    # demographics
    "sex",
    "age_5y_cat",
    "ethn_cat",
    # area
    "wimd2019_quintile",
    "health_board",
    "urban_rural_class",
    # qcovid
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
    #"qc_b_epilepsy",
    "qc_b_fracture4",
    "qc_b_neurorare",
    "qc_b_parkinsons",
    "qc_b_pulmhyper",
    "qc_b_pulmrare",
    "qc_b_pvd",
    #"qc_b_ra_sle",
    "qc_b_respcancer",
    "qc_b_semi",
    "qc_b_sicklecelldisease",
    #"qc_b_stroke",
    "qc_diabetes1",
    "qc_diabetes2",
    "qc_b_vte",
    "qc_chemo_cat",
    "qc_carehome",
    "qc_homeless",
    #"qc_learndisab",
    #"qc_downs",
    "qc_p_radio6",
    "qc_ckd3",
    "qc_ckd4",
    "qc_ckd5",
    # testing
    "test_pre08dec2020_cat"
)

pred.mat[
    rownames(pred.mat) == "qc_bmi_log",
    colnames(pred.mat) %in% bmi_xvar
] <- 1


# method, choose from: sample, norm, norm.boot, rf
mice.method <- rep("", ncol(d_cohort))
mice.method[rownames(pred.mat) == "qc_bmi_log"] <- "norm"

md_imp_bmi <- mice(
    data            = d_cohort,
    m               = 1,
    predictorMatrix = pred.mat,
    method          = mice.method,
    printFlag       = TRUE,
    maxit           = 5
)

qsave(md_imp_bmi, file = s_drive("md_imp_bmi.qs"))

md_imp_bmi <- qread(s_drive("md_imp_bmi.qs"))

# summary plot of BMI imputations
d_bmi_imp <-
    md_imp_bmi %>%
    complete(
        action = "long",
        include = TRUE
    )

p_bmi_imp <-
    d_bmi_imp %>%
    select(imp = .imp, alf_e, qc_bmi_log, sex, age_4cat, sex) %>%
    mutate(imp = factor(imp, 0:1, c("observed", "cmp1"))) %>%
    ggplot(aes(x = exp(qc_bmi_log), group = imp, colour = imp)) +
    geom_density() +
    facet_grid(age_4cat ~ sex) +
    labs(x = "qc_bmi_log") +
    xlim(10, 52)

ggsave(
    p_bmi_imp,
    file = "results/p_impute_bmi_density.png",
    width = p_width * 2,
    height = p_height
)

print(p_bmi_imp)

# Save =========================================================================
cat("Save\n")

d_bmi_imp <-
    md_imp_bmi %>%
    complete(
        action = "broad",
        include = TRUE
    ) %>%
    rename(
        alf_e   = alf_e.0,
        bmi_imp = qc_bmi_log.1
    ) %>%
    select(
        alf_e,
        bmi_imp
    ) %>%
    mutate(
        bmi_imp = exp(bmi_imp)
    )

qsave(
    x = d_bmi_imp,
    file = s_drive("d_bmi_imp.qs")
)

# Join to cohort ===============================================================
cat("Join to cohort\n")

d_cohort_clean <- qread(s_drive("d_cohort_clean.qs"))

if (any(names(d_cohort_clean) == "bmi_imp"))
    stop("d_cohort_clean already has bmi_imp")

d_bmi_imp <- qread(s_drive("d_bmi_imp.qs"))

d_cohort_clean <-
    d_cohort_clean %>%
    left_join(d_bmi_imp, by = "alf_e") %>%
    mutate(
        bmi_cat = case_when(
                              bmi_imp < 18.5 ~ "<18.5",     # underweight
            18.5 <= bmi_imp & bmi_imp < 25.0 ~ "18.5-24.9", # healthly weight
            25.0 <= bmi_imp & bmi_imp < 30.0 ~ "25.0-29.9", # overweight
            30.0 <= bmi_imp & bmi_imp < 35.0 ~ "30.0-34.9", # obese - part 1
            35.0 <= bmi_imp & bmi_imp < 40.0 ~ "35.0-39.9", # obese - part 2
            40.0 <= bmi_imp                  ~ "40.0+"
        ) %>%
        factor() %>%
        fct_relevel("<18.5")
    )

qsave(
    d_cohort_clean,
    file = s_drive("d_cohort_clean.qs")
)

d_cohort_clean %>% count(bmi_cat) %>% print()

beep()
