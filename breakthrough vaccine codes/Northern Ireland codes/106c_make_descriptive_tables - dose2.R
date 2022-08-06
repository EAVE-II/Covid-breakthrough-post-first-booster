cat("Clearing workspace and loading stuff\n")

source("scripts/000_libraries_and_functions.R")

# Load =========================================================================
cat("Load\n")

d_sample <- qread("results/d_sample.qs")

d_analysis <- qread("results/d_analysis_stats.qs")

# Cleaning =====================================================================
cat("Clean\n")

# Need to make sure only MD and PB included
# Check vaccine doseb names
# d_sample %>% group_by(vacc_doseb_name) %>% tally
# remove NAs

d_sample = d_sample %>%
    filter(!is.na(vacc_doseb_name) & vacc_doseb_name!="AZ")
    

names(d_analysis)

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
        vacc_dose_interval, #dose and follow up time combined 
        # covariates
        sex,
        age_gp,
        vacc_dose1_dose2_diff_cat,
        dose2_prior_infection_cat,
        nimdm,
        ur_combined,
        lgd,
        BNF_group,
        n_tests_gp,
    # background
        bg_interval,
        bg_avg_infection,
        bg_variant,
        # study design
        sample_weight
    ) %>%
    mutate(
        # keep only dose number and vaccine name for the interval
        vacc_dose_name = vacc_dose_interval %>%
            as.character() %>%
            stringr::str_replace("_day[0-9]+", "") %>%
            factor()
    )

d_analysis %>% count(vacc_dose_name)

# Calculate sample counts and column percentages ===============================
cat("Calculate column percentages\n")

make_tbl_np <- function(vacc, xvar) {
    tbl_n <-
        d_sample %>%
        janitor::tabyl(!!sym(xvar), !!sym(vacc), na.rm = FALSE) %>%
        janitor::adorn_totals("col") %>%
        tibble::as_tibble() %>%
        mutate(xvar = xvar) %>%
        rename(xlvl = !!sym(xvar)) %>%
        dplyr::select(xvar, everything()) %>%
        rename_with(
            .cols = where(is.numeric),
            .fn = ~ stringr::str_c(.x, "_n") 
        )

    tbl_p <-
        tbl_n %>%
        mutate(across(
            .cols = where(is.numeric),
            .fns = ~ (.x / sum(.x))*100
        )) %>%
        rename_with(
            .cols = where(is.numeric),
            .fn = ~ stringr::str_replace(.x, "_n", "_p")
        )

    tbl <-
        left_join(tbl_n, tbl_p, by = c("xvar", "xlvl")) %>%
        dplyr::select(
            xvar,
            xlvl,
            matches("Total"),
            matches("AZ"),
            matches("MD"),
            matches("PB")
        )

    return(tbl)
}

t_dose2_np <- rbind(
    make_tbl_np("vacc_dose2_name", "sex"),
    make_tbl_np("vacc_dose2_name", "age_gp"),
    #make_tbl_np("vacc_dose2_name", "lgd"),
    make_tbl_np("vacc_dose2_name", "vacc_dose1_dose2_diff_cat"),
    make_tbl_np("vacc_dose2_name", "dose2_prior_infection_cat"),
    make_tbl_np("vacc_dose2_name", "nimdm"),
    make_tbl_np("vacc_dose2_name", "ur_combined"),
    make_tbl_np("vacc_dose2_name", "BNF_group"),
    make_tbl_np("vacc_dose2_name", "n_tests_gp")
)

t_dose2_np <- t_dose2_np %>%
    rename("Characteristic" = xvar, "Level" = xlvl, "Total (number)" = Total_n,
    "Total (percentage)" = Total_p, "AZ (number)" = AZ_n, "AZ (percentage)" = AZ_p,
    "Pfizer (number)" = PB_n, "Pfizer (percentage)" = PB_p)

t_dose2_np$xvar = plyr::revalue(t_dose2_np$Characteristic, c("age_gp" = "Age group (years)",
                                                   "vacc_dose1_dose2_diff_cat" = "Gap between vaccine doses",
                                                   "dose2_prior_infection_cat" = "Prior history of COVID-19",
                                                   "nimdm" = "Deprivation status",
                                                   "ur_combined" = "Urban/Rural index",
                                                   "BNF_group" = "Number of risk groups (BNF)",
                                                   "n_tests_gp" = "Number of previous tests"))


# Calculate person-years and rates =============================================
cat("Calculate rates\n")

tidy_pyears <- function(x) {
    d_analysis %>%
    calc_pyears(vacc_dose_name, !!sym(x)) %>%
    rename(xlvl = !!sym(x)) %>%
    mutate(xvar = x) %>%
    dplyr::select(vacc_dose_name, xvar, xlvl, pyears, event)
}

d_py <- rbind(
    tidy_pyears("sex"),
    tidy_pyears("age_gp"), # age_3cat age_4cat age_5y_cat
    tidy_pyears("lgd"),
    tidy_pyears("vacc_dose1_dose2_diff_cat"),
    tidy_pyears("dose2_prior_infection_cat"),
    tidy_pyears("nimdm"),
    tidy_pyears("ur_combined"),
    tidy_pyears("BNF_group"),
    tidy_pyears("n_tests_gp"),
    tidy_pyears("bg_variant")
)

calc_rate <- function(dose) {
    d_py %>%
    filter(stringr::str_detect(vacc_dose_name, dose)) %>%
    mutate(
        vacc_dose_name = stringr::str_replace(vacc_dose_name, "dose[2b]_", ""),
        rate = event / pyears * 1000
    ) %>%
    tidyr::pivot_wider(
        names_from = vacc_dose_name,
        names_glue = "{vacc_dose_name}_{.value}",
        values_from = c(pyears, event, rate),
    ) %>%
    rowwise() %>%
    mutate(
        Total_pyears = sum(c_across(matches("_pyears"))),
        Total_event  = sum(c_across(matches("_event"))),
        Total_rate   = Total_event / Total_pyears * 1000
    ) %>%
    ungroup() %>%
    dplyr::select(
        xvar,
        xlvl,
        matches("Total"),
        matches("AZ"),
        matches("MD"),
        matches("PB")
    )
}

t_dose2_rate <- calc_rate("dose2") %>%
    rename("Characteristic" = xvar, "Level" = xlvl, "Total person years" = Total_pyears)

t_dose2_rate$Characteristic = plyr::revalue(t_dose2_rate$Characteristic, c("age_gp" = "Age group (years)",
                                                                           "vacc_dose1_dose2_diff_cat" = "Gap between vaccine doses",
                                                                           "dose2_prior_infection_cat" = "Prior history of COVID-19",
                                                                           "nimdm" = "Deprivation status",
                                                                           "ur_combined" = "Urban/Rural index",
                                                                           "BNF_group" = "Number of risk groups (BNF)",
                                                                           "n_tests_gp" = "Number of previous tests"))

# Save =========================================================================
cat("Save\n")

qsave(t_dose2_np,   file = "results/t_dose2_np_sample.qs")
qsave(t_dose2_rate, file = "results/t_dose2_rate_sample.qs")

write.csv(t_dose2_np, "output/national_study/Outputs for approval/Dose2/dose2_descriptive_np.csv", row.names=FALSE)
write.csv(t_dose2_rate, "output/national_study/Outputs for approval/Dose2/dose2_descriptive_rates.csv", row.names=FALSE)
