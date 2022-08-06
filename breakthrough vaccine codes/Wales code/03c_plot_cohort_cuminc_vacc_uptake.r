source("r_clear_and_load.r")

# Load clean sample ============================================================
cat("Load clean sample\n")

d_cohort_clean <- qread(s_drive("d_cohort_clean.qs"))

# d_cohort_clean %>% glimpse()

# Make vacc state dataset ======================================================
cat("Make vacc state dataset\n")

# make a multi-state data set to calculate cumulative incidences for transitions
# between vaccination statues: dose 1, 2, 3, booster, and terminal events:
# death, move out of cohort, end of follow-up

d_vacc_state <-
    d_cohort_clean %>%
    mutate(
        stop_date = pmin(
            death_covid_date,
            death_noncovid_date,
            c20_end_date,
            study_end_date,
            na.rm = TRUE
        ),
        stop_cat = case_when(
            stop_date == death_noncovid_date ~ "death",
            stop_date == death_covid_date    ~ "death",
            stop_date == c20_end_date        ~ "move_out",
            stop_date == study_end_date      ~ "study_end"
        )
    ) %>%
    select(
        alf_e,
        age_3cat,
        vacc_dose1_date,
        vacc_dose1_name,
        vacc_dose2_date,
        vacc_dose2_name,
        vacc_doseb_date,
        vacc_doseb_name,
        stop_date,
        stop_cat
    ) %>%
    mutate(
        # if we have vacc records after follow-up ends, blank them
        vacc_dose1_date = if_else(vacc_dose1_date >= stop_date, NA_Date_,      vacc_dose1_date, vacc_dose1_date),
        vacc_dose1_name = if_else(vacc_dose1_date >= stop_date, NA_character_, vacc_dose1_name, vacc_dose1_name),
        vacc_dose2_date = if_else(vacc_dose2_date >= stop_date, NA_Date_,      vacc_dose2_date, vacc_dose2_date),
        vacc_dose2_name = if_else(vacc_dose2_date >= stop_date, NA_character_, vacc_dose2_name, vacc_dose2_name),
        vacc_doseb_date = if_else(vacc_doseb_date >= stop_date, NA_Date_,      vacc_doseb_date, vacc_doseb_date),
        vacc_doseb_name = if_else(vacc_doseb_date >= stop_date, NA_character_, vacc_doseb_name, vacc_doseb_name)
    ) %>%
    pivot_longer(
        cols      = matches("_date$"),
        names_to  = "event_type",
        values_to = "event_date"
    ) %>%
    filter(not_na(event_date)) %>%
    mutate(
        event_type = case_when(
            event_type == "stop_date" ~ stop_cat,
            event_type == "vacc_dose1_date"   ~ str_c("dose1_", vacc_dose1_name),
            event_type == "vacc_dose2_date"   ~ str_c("dose2_", vacc_dose2_name),
            event_type == "vacc_doseb_date"   ~ str_c("doseb_", vacc_doseb_name)
        ),
        event_type = factor(event_type) %>% fct_relevel("study_end"),
        event_type_simple = str_replace(event_type, "(dose[12b])_.*", "\\1"),
        event_type_simple = factor(event_type_simple) %>% fct_relevel("study_end")
    ) %>%
    select(
        -stop_cat,
        -vacc_dose1_name,
        -vacc_dose2_name,
        -vacc_doseb_name
    ) %>%
    # make survival analysis variables
    mutate(
        tstop = interval(vacc_start_date - ddays(1), event_date) / ddays()
    ) %>%
    lazy_dt() %>%
    group_by(alf_e) %>%
    mutate(tstart = lag(tstop, default = 0)) %>%
    ungroup() %>%
    as_tibble() %>%
    select(
        alf_e,
        age_3cat,
        tstart,
        tstop,
        event_type_simple,
        event_type,
        event_date,
        everything()
    ) %>%
    # basic checks to make sure nothing is broke
    verify(tstart >= 0) %>%
    verify(tstop > 0) %>%
    verify(tstart < tstop) %>%
    verify(not_na(event_type))


# Save =========================================================================
cat("Save\n")

qsave(d_vacc_state, file = s_drive("d_vacc_state.qs"))


# Load vacc state ==============================================================
cat("Load vacc state\n")

d_vacc_state <- qread(s_drive("d_vacc_state.qs"))


# Fit empirical survival curves ================================================
cat("Fit empirical survival curves\n")

sf_vacc_state <- survfit(
    Surv(tstart, tstop, event_type) ~ age_3cat,
    data = d_vacc_state,
    id = alf_e
)

# Plot survival curves =========================================================
cat("Plot survival curves\n")

expr_dose <- "(dose[123b])_(.+)"

epoch <- vacc_start_date - ddays(1)

lkp_state <- c(
    "Unvaccinated" = "(s0)",
    "Move out"     = "move_out",
    "Death"        = "death",
    "First dose"   = "dose1",
    "Second dose"  = "dose2",
    "Booster"      = "doseb"
)

lkp_vacc_name <- c(
    "Pfizer"      = "PB",
    "AstraZeneca" = "AZ",
    "Moderna"     = "MD"
)

lkp_fill <- c(
    "Pfizer"       = "#66c2a5",
    "AstraZeneca"  = "#fc8d62",
    "Moderna"      = "#8da0cb"
)

p_vacc_state <-
    sf_vacc_state %>%
    tidy() %>%
    mutate(
        strata    = str_replace(strata, "age_3cat=", ""),
        strata    = str_replace(strata, "_", "-"),
        vacc_name = str_replace(state, expr_dose, "\\2"),
        state     = str_replace(state, expr_dose, "\\1"),
        state     = factor(state, lkp_state, names(lkp_state)),
        vacc_name = factor(vacc_name, lkp_vacc_name, names(lkp_vacc_name)),
        time      = epoch + ddays(time)
    ) %>%
    ggplot(aes(x = time, y = estimate, group = vacc_name, fill = vacc_name)) +
    facet_grid(state ~ strata) +
    geom_area() +
    scale_x_date(
        name = "",
        date_labels = "%b\n%Y"
    ) +
    scale_y_continuous(
        name = "Cumulative incidence",
        breaks = pretty_breaks()
    ) +
    scale_fill_manual(
        breaks = names(lkp_fill),
        values = lkp_fill,
        na.value = "grey50"
    ) +
    theme(
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank()
    )

print(p_vacc_state)


# Save =========================================================================
cat("Save\n")

qsavem(
    sf_vacc_state,
    p_vacc_state,
    file = "results/cuminc_uptake.qsm"
)
