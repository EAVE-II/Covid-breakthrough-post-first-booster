source("r_clear_and_load.r")


# Load =========================================================================
cat("Load\n")

d_cohort_raw <- qread(s_drive("d_cohort_raw.qs"))


# Collapse dose 3 and booster ==================================================
cat("Collapse dose 3 and booster\n")

d_cohort_raw <- d_cohort_raw %>%
    mutate(
        vacc_dose3b_date = pmin(vacc_dose3_date, vacc_doseb_date, na.rm = TRUE),
        vacc_dose3b_name = if_else(
            condition = not_na(vacc_dose3_date) & vacc_dose3b_date == vacc_dose3_date,
            true = vacc_dose3_name,
            false = vacc_doseb_name
        )
    ) %>%
    select(
        -matches("vacc_dose3_.+"),
        -matches("vacc_doseb_.+")
    ) %>%
    rename(
        vacc_doseb_date = vacc_dose3b_date,
        vacc_doseb_name = vacc_dose3b_name
    )


# Select sample ================================================================
cat("Select cohort\n")

lkp_vacc_primary <- c(
    "Astrazeneca",
    "Pfizer Biontech",
    "Pfizer child"
)

lkp_vacc_booster <- c(
    "Pfizer Biontech",
    "Pfizer child",
    "Moderna"
)


d_cohort_raw <-
    d_cohort_raw %>%
    mutate(
        age = floor(interval(wob, vacc_start_date) / dyears(1))
    ) %>%
    mutate(
        has_wob_sex           = not_na(wob) &
                                not_na(gndr_cd),
        is_18_105             = age >= 18 & age <= 105,
        is_welsh_resident     = c20_start_date <= ymd('2020-01-01') &
                                c20_end_date >= vacc_start_date &
                                (is.na(death_date) | death_date >= ymd('2020-12-08')),
        has_ralf              = not_na(ralf_e) &
                                not_na(lsoa2011_cd),
        has_hh_lte10          = hh_all_n <= 10,
        has_sail_gp           = gp_end_date >= vacc_start_date,
        has_qcovid            = qc_flg == 1,
        has_good_vacc_record  = (is.na(has_bad_vacc_record) | has_bad_vacc_record == 0) &
                                (is.na(vacc_dose1_date) | vacc_dose1_date >= ymd('2020-12-08')),
        is_study_vacc_primary = (is.na(vacc_dose1_date) | vacc_dose1_name %in% lkp_vacc_primary) &
                                (is.na(vacc_dose2_date) | vacc_dose2_name %in% lkp_vacc_primary) &
                                (is.na(vacc_dose2_date) | vacc_dose2_name == vacc_dose1_name),
        is_study_booster      = (is.na(vacc_doseb_date) | vacc_doseb_name %in% lkp_vacc_booster)
    )


d_cohort_raw  %>% count(age)  %>% arrange(desc(age))


# sample selection summary -----------------------------------------------------

t_cohort_selection <- tribble(
    ~step, ~criteria, ~n,
        1, "In cohort from Jan 2020 to Dec 2020",
            d_cohort_raw %>% nrow(),
        2, "Has wob and sex recorded",
            d_cohort_raw %>% filter(has_wob_sex) %>% nrow(),
        3, "Is aged between 18 and 105 from Dec 2020",
            d_cohort_raw %>% filter(has_wob_sex, is_18_105) %>% nrow(),
        4, "Is Welsh resident from Jan 2020",
            d_cohort_raw %>% filter(has_wob_sex, is_18_105, is_welsh_resident) %>% nrow(),
        5, "Has RALF",
            d_cohort_raw %>% filter(has_wob_sex, is_18_105, is_welsh_resident, has_ralf) %>% nrow(),
        6, "Has household size <= 10",
            d_cohort_raw %>% filter(has_wob_sex, is_18_105, is_welsh_resident, has_ralf, has_hh_lte10) %>% nrow(),
        7, "Is registered with SAIL GP",
            d_cohort_raw %>% filter(has_wob_sex, is_18_105, is_welsh_resident, has_ralf, has_hh_lte10, has_sail_gp) %>% nrow(),
        8, "Has QCovid measures",
            d_cohort_raw %>% filter(has_wob_sex, is_18_105, is_welsh_resident, has_ralf, has_hh_lte10, has_sail_gp, has_qcovid) %>% nrow(),
        9, "Has good vacc records",
            d_cohort_raw %>% filter(has_wob_sex, is_18_105, is_welsh_resident, has_ralf, has_hh_lte10, has_sail_gp, has_qcovid, has_good_vacc_record) %>% nrow(),
       10, "If has first and second dose, then with AZ or PB",
            d_cohort_raw %>% filter(has_wob_sex, is_18_105, is_welsh_resident, has_ralf, has_hh_lte10, has_sail_gp, has_qcovid, has_good_vacc_record, is_study_vacc_primary) %>% nrow(),
       11, "If has booster dose, then with MD or PB",
            d_cohort_raw %>% filter(has_wob_sex, is_18_105, is_welsh_resident, has_ralf, has_hh_lte10, has_sail_gp, has_qcovid, has_good_vacc_record, is_study_vacc_primary, is_study_booster) %>% nrow(),
    ) %>%
    mutate(
        n_diff = n - lag(n),
        p_diff = round(n_diff / first(n), 3) * 100
    )

print(as.data.frame(t_cohort_selection))


# apply criteria ---------------------------------------------------------------

d_cohort_raw <-
    d_cohort_raw %>%
    filter(
        has_wob_sex,
        is_18_105,
        is_welsh_resident,
        has_ralf,
        has_hh_lte10,
        has_sail_gp,
        has_qcovid,
        has_good_vacc_record,
        is_study_vacc_primary,
        is_study_booster
    )


# Cleaning =====================================================================
cat("Cleaning\n")

rename_vacc <- function(x) {
    case_when(
        x == "Astrazeneca"     ~ "AZ",
        x == "Moderna"         ~ "MD",
        x == "Pfizer Biontech" ~ "PB",
        x == "Pfizer child"    ~ "PB"
    )
}

lkp_wimd <- c(
    "1most"  = "1",
    "2"      = "2",
    "3"      = "3",
    "4"      = "4",
    "5least" = "5"
)

d_cohort_clean <-
    d_cohort_raw %>%
    mutate(
        # demographics
        age_3cat = case_when(
                age <= 17                ~ "17_and_under",
                age %>% between(18,  64) ~ "18_64",
                age %>% between(65,  79) ~ "65_79",
                age >= 80                ~ "80+"
            ) %>%
            factor() %>%
            fct_explicit_na(),
        age_4cat = case_when(
                age <= 17                ~ "17_and_under",
                age %>% between(18,  49) ~ "18_49",
                age %>% between(50,  64) ~ "50_64",
                age %>% between(65,  79) ~ "65_79",
                age >= 80                ~ "80+"
            ) %>%
            factor() %>%
            fct_explicit_na(),
        age_5y_cat = case_when(
                age <= 17                ~ "17_and_under",
                age %>% between(18,  49) ~ "18_49",
                age %>% between(50,  54) ~ "50_54",
                age %>% between(55,  59) ~ "55_59",
                age %>% between(60,  64) ~ "60_64",
                age %>% between(65,  69) ~ "65_69",
                age %>% between(70,  74) ~ "70_74",
                age %>% between(75,  79) ~ "75_79",
                age >= 80                ~ "80+"
            ) %>%
            factor() %>%
            fct_explicit_na(),
        sex = factor(gndr_cd, 1:2, c("Male", "Female")) %>% fct_relevel("Female"),
        ethn_cat = fct_infreq(ethn_cat) %>% fct_explicit_na("Unknown")
    ) %>%
    mutate(
        # area
        wimd2019_quintile = fct_recode(as.character(wimd2019_quintile), !!!lkp_wimd),
        health_board = factor(health_board),
        urban_rural_class =
            urban_rural_class %>%
            factor() %>%
            fct_collapse(
                "Urban" = c(
                    "C1 Urban city and town",
                    "C2 Urban city and town in a sparse setting"
                ),
                "Rural" = c(
                    "D1 Rural town and fringe",
                    "D2 Rural town and fringe in a sparse setting",
                    "E1 Rural village and dispersed",
                    "E2 Rural village and dispersed in a sparse setting"
                )
            )
    ) %>%
    mutate(
        # vaccination
        has_vacc_dose1  = as.numeric(not_na(vacc_dose1_date)),
        has_vacc_dose2  = as.numeric(not_na(vacc_dose2_date)),
        has_vacc_doseb  = as.numeric(not_na(vacc_doseb_date)),
        vacc_dose1_name = rename_vacc(vacc_dose1_name),
        vacc_dose2_name = rename_vacc(vacc_dose2_name),
        vacc_doseb_name = rename_vacc(vacc_doseb_name)
    ) %>%
    mutate(
        # time between 1st and 2nd dose
        vacc_dose1_dose2_diff_week = floor(interval(vacc_dose1_date, vacc_dose2_date) / dweeks(1)),
        vacc_dose1_dose2_diff_cat = case_when(
                has_vacc_dose1 == 0 | has_vacc_dose2 == 0   ~ "No dose 1 or dose 2",
                between(vacc_dose1_dose2_diff_week,  0,  6) ~ "00-06wk",
                between(vacc_dose1_dose2_diff_week,  7,  8) ~ "07-08wk",
                between(vacc_dose1_dose2_diff_week,  9, 10) ~ "09-10wk",
                between(vacc_dose1_dose2_diff_week, 11, 12) ~ "11-12wk",
                vacc_dose1_dose2_diff_week >= 13            ~ "13+wk"
            ) %>% factor()
    ) %>%
    mutate(
        # prior PCR testing
        test_pre08dec2020_cat = case_when(
                is.na(test_pre08dec2020_n)         ~ "00",
                between(test_pre08dec2020_n, 0, 2) ~ str_pad(test_pre08dec2020_n, width = 2, pad = "0"),
                between(test_pre08dec2020_n, 3, 9) ~ "03-09",
                test_pre08dec2020_n >= 10          ~ "10+"
            ) %>%
            factor(),
        test_pre_dose2_cat = case_when(
                is.na(test_pre_dose2_n)         ~ "00",
                between(test_pre_dose2_n, 0, 2) ~ str_pad(test_pre_dose2_n, width = 2, pad = "0"),
                between(test_pre_dose2_n, 3, 4) ~ "03-04",
                between(test_pre_dose2_n, 5, 9) ~ "05-09",
                test_pre_dose2_n >= 10          ~ "10+"
            ) %>%
            factor()
    ) %>%
    mutate(
        has_covid_infection = as.numeric(not_na(infection1_test_date)),
        has_covid_hosp      = as.numeric(not_na(hosp_admis_date)),
        # time since most recent infection prior to dose 2
        dose2_prior_infection_date = case_when(
                !is.na(infection4_test_date) & infection4_test_date < vacc_dose2_date ~ infection4_test_date,
                !is.na(infection3_test_date) & infection3_test_date < vacc_dose2_date ~ infection3_test_date,
                !is.na(infection2_test_date) & infection2_test_date < vacc_dose2_date ~ infection2_test_date,
                !is.na(infection1_test_date) & infection1_test_date < vacc_dose2_date ~ infection1_test_date
            ),
        dose2_prior_infection_week = floor(interval(dose2_prior_infection_date, vacc_dose2_date) / dweeks()),
        dose2_prior_infection_cat = case_when(
            is.na(dose2_prior_infection_date)           ~ "No prior infection",
            between(dose2_prior_infection_week,  0, 12) ~ "00-12wk", # 0-2 months
            between(dose2_prior_infection_week, 13, 25) ~ "13-25wk", # 3-5 months
            between(dose2_prior_infection_week, 26, 38) ~ "26-38wk", # 6-8 months
            dose2_prior_infection_week >= 39            ~ "39+wk",   # 9+ months
        ) %>% factor()
    ) %>%
    mutate(
        # health care utilisation
        hcu_hosp_spell_n   = if_else(qc_flg == 1, as.integer(replace_na(hcu_hosp_spell_n, 0)), hcu_hosp_spell_n),
        hcu_hosp_spell_cat = case_when(
                hcu_hosp_spell_n <= 1 ~ as.character(hcu_hosp_spell_n),
                hcu_hosp_spell_n >= 2 ~ "2+"
            ) %>%
            factor(),
        hcu_gp_attendance_n   = if_else(qc_flg == 1, as.integer(replace_na(hcu_gp_attendance_n, 0)), hcu_gp_attendance_n),
        hcu_gp_attendance_cat = case_when(
                hcu_gp_attendance_n == 0             ~ "00",
                between(hcu_gp_attendance_n, 01, 19) ~ "01-19",
                between(hcu_gp_attendance_n, 20, 34) ~ "20-34",
                between(hcu_gp_attendance_n, 35, 59) ~ "35-59",
                hcu_gp_attendance_n >= 60            ~ "60plus"
            ) %>%
            factor(),
        hcu_gp_prescription_n   = if_else(qc_flg == 1, as.integer(replace_na(hcu_gp_prescription_n, 0)), hcu_gp_prescription_n),
        hcu_gp_prescription_cat = case_when(
                hcu_gp_prescription_n == 0             ~ "00",
                between(hcu_gp_prescription_n, 01, 06) ~ "01-06",
                between(hcu_gp_prescription_n, 07, 19) ~ "07-19",
                between(hcu_gp_prescription_n, 20, 49) ~ "20-49",
                hcu_gp_prescription_n >= 50            ~ "50plus"
            ) %>%
            factor()
    ) %>%
    mutate(
        # health
        hypt_cat = factor(hypertension_flg, 0:1, c("No", "Yes")),
        shielded_cat = factor(shielded_flg, 0:1, c("No", "Yes"))
    ) %>%
    mutate(
        # death
        death_noncovid_date = if_else(death_covid_flg == 0, death_date, NA_Date_),
        death_covid_date    = if_else(death_covid_flg == 1, death_date, NA_Date_),
        has_death_noncovid  = as.numeric(not_na(death_noncovid_date)),
        has_death_covid     = as.numeric(not_na(death_covid_date))
    ) %>%
    mutate(
        # tidy dates
        c20_end_date   = na_if(c20_end_date, max(c20_end_date)),
        study_end_date = study_end_date,
    )

# QCovid score -----------------------------------------------------------------

d_cohort_clean <-
    d_cohort_clean %>%
    # split out some of the categories into separate indicators
    mutate(
        qc_diabetes1     = as.numeric(qc_diabetes_cat == 1),
        qc_diabetes2     = as.numeric(qc_diabetes_cat == 2),
        qc_carehome      = as.numeric(qc_home_cat     == 1),
        qc_homeless      = as.numeric(qc_home_cat     == 2),
        qc_learndisab    = as.numeric(qc_learn_cat    == 1),
        qc_downs         = as.numeric(qc_learn_cat    == 2),
        qc_ckd3          = as.numeric(qc_renal_cat    == 2),
        qc_ckd4          = as.numeric(qc_renal_cat    == 3),
        qc_ckd5          = as.numeric(qc_renal_cat    >= 4),
        renal_transplant = as.numeric(qc_renal_cat    == 6)
    ) %>%
    # turn everything into 0/1 flag and count up how many items someone has
    mutate(
        qc_b2_82               = as.numeric(qc_b2_82               == 1),
        qc_b2_leukolaba        = as.numeric(qc_b2_leukolaba        == 1),
        qc_b2_prednisolone     = as.numeric(qc_b2_prednisolone     == 1),
        qc_b_af                = as.numeric(qc_b_af                == 1),
        qc_b_ccf               = as.numeric(qc_b_ccf               == 1),
        qc_b_asthma            = as.numeric(qc_b_asthma            == 1),
        qc_b_bloodcancer       = as.numeric(qc_b_bloodcancer       == 1),
        qc_b_cerebralpalsy     = as.numeric(qc_b_cerebralpalsy     == 1),
        qc_b_chd               = as.numeric(qc_b_chd               == 1),
        qc_b_cirrhosis         = as.numeric(qc_b_cirrhosis         == 1),
        qc_b_congenheart       = as.numeric(qc_b_congenheart       == 1),
        qc_b_copd              = as.numeric(qc_b_copd              == 1),
        qc_b_dementia          = as.numeric(qc_b_dementia          == 1),
        qc_b_epilepsy          = as.numeric(qc_b_epilepsy          == 1),
        qc_b_fracture4         = as.numeric(qc_b_fracture4         == 1),
        qc_b_neurorare         = as.numeric(qc_b_neurorare         == 1),
        qc_b_parkinsons        = as.numeric(qc_b_parkinsons        == 1),
        qc_b_pulmhyper         = as.numeric(qc_b_pulmhyper         == 1),
        qc_b_pulmrare          = as.numeric(qc_b_pulmrare          == 1),
        qc_b_pvd               = as.numeric(qc_b_pvd               == 1),
        qc_b_ra_sle            = as.numeric(qc_b_ra_sle            == 1),
        qc_b_respcancer        = as.numeric(qc_b_respcancer        == 1),
        qc_b_semi              = as.numeric(qc_b_semi              == 1),
        qc_b_sicklecelldisease = as.numeric(qc_b_sicklecelldisease == 1),
        qc_b_stroke            = as.numeric(qc_b_stroke            == 1),
        qc_diabetes_cat        = as.numeric(qc_diabetes_cat        != 0),
        qc_b_vte               = as.numeric(qc_b_vte               == 1),
        qc_chemo_cat           = as.numeric(qc_chemo_cat           != 0),
        qc_home_cat            = as.numeric(qc_home_cat            != 0),
        qc_learn_cat           = as.numeric(qc_learn_cat           != 0),
        qc_p_radio6            = as.numeric(qc_p_radio6            == 1),
        qc_p_solidtransplant   = as.numeric(qc_p_solidtransplant   == 1),
        qc_renal_cat           = as.numeric(qc_renal_cat           != 1),
        # simple sum score of how many QCovid risk factors does someone have
        qcovid_score =
            qc_b2_82 +
            qc_b2_leukolaba +
            qc_b2_prednisolone +
            qc_b_af +
            qc_b_ccf +
            qc_b_asthma +
            qc_b_bloodcancer +
            qc_b_cerebralpalsy +
            qc_b_chd +
            qc_b_cirrhosis +
            qc_b_congenheart +
            qc_b_copd +
            qc_b_dementia +
            qc_b_epilepsy +
            qc_b_fracture4 +
            qc_b_neurorare +
            qc_b_parkinsons +
            qc_b_pulmhyper +
            qc_b_pulmrare +
            qc_b_pvd +
            qc_b_ra_sle +
            qc_b_respcancer +
            qc_b_semi +
            qc_b_sicklecelldisease +
            qc_b_stroke +
            qc_diabetes_cat +
            qc_b_vte +
            qc_chemo_cat +
            qc_home_cat +
            qc_learn_cat +
            qc_p_radio6 +
            qc_p_solidtransplant +
            qc_renal_cat,
        qcovid_cat = case_when(
                qcovid_score >= 5 ~ "5+",
                TRUE              ~ as.character(qcovid_score)
            ) %>% factor()
    ) %>%
    select(
        -qc_diabetes_cat,
        -qc_home_cat,
        -qc_learn_cat,
        -qc_renal_cat
    )

# Save =========================================================================
cat("Save\n")

qsave(
    t_cohort_selection,
    file = "results/t_cohort_selection.qs"
)

qsave(
    d_cohort_clean,
    file = s_drive("d_cohort_clean.qs")
)

beep()
