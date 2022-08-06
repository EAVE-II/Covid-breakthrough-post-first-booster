source("r_clear_and_load.r")

# Get raw cohort ===============================================================
cat("Get raw cohort\n")

con <- sail_open()

q_cohort <- "
WITH
    dacvap_hosp AS
    (
        SELECT alf_e, MIN(epi_start_date) AS hosp_admis_date
        FROM sailw1151v.dacvap_hosp_long
        WHERE
            spell_admis_method = 'Emergency admission'
            AND (
                admis_covid19_cause_flg = 1
                OR admis_with_covid19_flg = 1
                OR covid19_during_admis_flg = 1
                OR within_14days_positive_pcr_flg = 1
            )
        GROUP BY alf_e
    )
SELECT
-- main id
    cohort.alf_e,
-- c19 cohort 20
    cohort.pers_id_e,
    cohort.wob,
    cohort.gndr_cd,
    cohort.wds_start_date,
    cohort.wds_end_date,
    cohort.gp_start_date,
    cohort.gp_end_date,
    cohort.c20_start_date,
    cohort.c20_end_date,
-- ethnicity
    cohort.ethn_cat,
-- health care worker
    cohort.hcw_dec2020_flg,
    cohort.hcw_sep2021_flg,
-- area and household info at 2020-12-07
    cohort.ralf_e,
    cohort.ralf_sts_cd,
    cohort.lsoa2011_cd,
    cohort.wimd2019_quintile,
    cohort.lad2011_nm,
    cohort.health_board,
    cohort.urban_rural_class,
    cohort.hh_all_n,
    cohort.hh_child_n,
    cohort.hh_adult_n,
-- shielded patient
    shielded_flg,
-- qcovid as at 2020-12-07
    cohort.qc_flg,
    cohort.qc_b2_82,
    cohort.qc_b2_leukolaba,
    cohort.qc_b2_prednisolone,
    cohort.qc_b_af,
    cohort.qc_b_ccf,
    cohort.qc_b_asthma,
    cohort.qc_b_bloodcancer,
    cohort.qc_b_cerebralpalsy,
    cohort.qc_b_chd,
    cohort.qc_b_cirrhosis,
    cohort.qc_b_congenheart,
    cohort.qc_b_copd,
    cohort.qc_b_dementia,
    cohort.qc_b_epilepsy,
    cohort.qc_b_fracture4,
    cohort.qc_b_neurorare,
    cohort.qc_b_parkinsons,
    cohort.qc_b_pulmhyper,
    cohort.qc_b_pulmrare,
    cohort.qc_b_pvd,
    cohort.qc_b_ra_sle,
    cohort.qc_b_respcancer,
    cohort.qc_b_semi,
    cohort.qc_b_sicklecelldisease,
    cohort.qc_b_stroke,
    cohort.qc_diabetes_cat,
    cohort.qc_b_vte,
    cohort.qc_bmi,
    cohort.qc_chemo_cat,
    cohort.qc_home_cat,
    cohort.qc_learn_cat,
    cohort.qc_p_marrow6,
    cohort.qc_p_radio6,
    cohort.qc_p_solidtransplant,
    cohort.qc_renal_cat,
-- hypertension over 5 years prior to 2020-12-07
    cohort.hypertension_flg,
-- health care utilisation over previous 2 years to 2020-12-07
    cohort.hcu_hosp_spell_n,
    cohort.hcu_gp_attendance_n,
    cohort.hcu_gp_prescription_n,
-- vaccine data quality
    vacc.has_bad_vacc_record,
-- first dose
    vacc.vacc_dose1_date,
    vacc.vacc_dose1_name,
    vacc.vacc_dose1_reaction_ind,
    vacc.vacc_dose1_reaction_cd,
-- second dose
    vacc.vacc_dose2_date,
    vacc.vacc_dose2_name,
    vacc.vacc_dose2_reaction_ind,
    vacc.vacc_dose2_reaction_cd,
-- third dose
    vacc.vacc_dose3_date,
    vacc.vacc_dose3_name,
    vacc.vacc_dose3_reaction_ind,
    vacc.vacc_dose3_reaction_cd,
-- booster dose
    vacc.vacc_doseb_date,
    vacc.vacc_doseb_name,
    vacc.vacc_doseb_reaction_ind,
    vacc.vacc_doseb_reaction_cd,
-- pcr test history
    pcr_test.test_ever_flg,
    pcr_test.test_pre08dec2020_n,
    pcr_test.test_pre16sep2021_n,
-- number of infections,
    pcr_test.infection_n,
-- positive tests 90 days apart
    pcr_test.infection1_test_date,
    pcr_test.infection2_test_date,
    pcr_test.infection3_test_date,
    pcr_test.infection4_test_date,
-- first hospital admission related to COVID-19
    dacvap_hosp.hosp_admis_date,
-- death
    death.death_date,
    death.death_covid_flg
FROM sailw1151v.dacvap_cohort AS cohort
LEFT JOIN sailw1151v.dacvap_vacc AS vacc
    ON cohort.alf_e = vacc.alf_e
LEFT JOIN sailw1151v.dacvap_pcr_test AS pcr_test
    ON cohort.alf_e = pcr_test.alf_e
LEFT JOIN dacvap_hosp
    ON cohort.alf_e = dacvap_hosp.alf_e
LEFT JOIN sailw1151v.dacvap_death AS death
    ON cohort.alf_e = death.alf_e
;"

d_cohort <- sail_run(con, q_cohort)


# Get number of pcr tests pre dose 2 ===========================================
cat("Get number of pcr tests pre dose 2\n")

q_test_pre_dose2 <- "
WITH
    pcr AS
    (
        -- if good, use spcm_collected_dt, else try spcm_recieved_dt
        SELECT
            cohort.alf_e,
            CASE
                WHEN test.spcm_collected_dt BETWEEN '2020-01-01' AND CURRENT DATE THEN test.spcm_collected_dt
                ELSE test.SPCM_RECEIVED_DT
            END AS test_date,
            test.covid19testresult AS test_result
        FROM
            sailw1151v.dacvap_cohort AS cohort
        INNER JOIN
            sail0911v.patd_df_covid_lims_testresults AS test
            ON cohort.alf_e = test.alf_e
        WHERE
            test.alf_e IS NOT NULL
            AND test.alf_sts_cd IN (1, 4, 39)
            AND test.spcm_collected_dt IS NOT NULL
            AND test.covid19testresult IS NOT NULL
    )
-- count number of tests before second dose
SELECT
    cohort.alf_e,
    COALESCE(SUM(pcr.test_date < vacc.vacc_dose2_date), 0) AS test_pre_dose2_n
FROM sailw1151v.dacvap_cohort AS cohort
LEFT JOIN
    sailw1151v.dacvap_vacc AS vacc
    ON cohort.alf_e = vacc.alf_e
LEFT JOIN pcr
    ON cohort.alf_e = pcr.alf_e
GROUP BY cohort.alf_e;"

d_test_pre_dose2 <- sail_run(con, q_test_pre_dose2)

# add on to cohort table
d_cohort <- d_cohort %>% left_join(d_test_pre_dose2, by = "alf_e")

# check alf_e is unique
if (nrow(d_cohort) != n_distinct(d_cohort$alf_e)) {
    stop("alf_e is not unique")
}


# Save =========================================================================
cat("Save\n")

qsave(
    d_cohort,
    file = s_drive("d_cohort_raw.qs")
)

# calculate study end date =====================================================

pcr_date <-
    d_cohort %>%
    filter(
        not_na(infection1_test_date) |
        not_na(infection2_test_date) |
        not_na(infection3_test_date) |
        not_na(infection4_test_date)
    ) %>%
    mutate(test_date = pmax(
        infection1_test_date,
        infection2_test_date,
        infection3_test_date,
        infection4_test_date,
        na.rm = TRUE
    )) %>%
    summarise(max = max(test_date))

vacc_date <-
    d_cohort %>%
    filter(vacc_doseb_date <= today()) %>%
    summarise(max = max(vacc_doseb_date))

if (study_end_date != min(pcr_date$max, vacc_date$max)) {
    cat(
        "You need to update study end date to be\n",
        as.character(min(pcr_date$max, vacc_date$max)), "\n",
        sep = ""
    )
}


# Goodbye ======================================================================

sail_close(con)
beep()
