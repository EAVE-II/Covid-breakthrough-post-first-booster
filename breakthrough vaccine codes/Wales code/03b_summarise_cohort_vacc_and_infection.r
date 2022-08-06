source("r_clear_and_load.r")

# Load clean cohort ============================================================
cat("Load clean cohort\n")

d_cohort <- qread(s_drive("d_cohort_clean.qs"))

d_bg <- qread(s_drive("d_background.qs"))

# Vaccinations =================================================================
cat("Vaccinations\n")

d_vacc <-
    d_cohort %>%
    mutate(
        vacc_dose2_diff = interval(vacc_dose1_date, vacc_dose2_date) / dweeks(1),
        vacc_doseb_diff = interval(vacc_dose2_date, vacc_doseb_date) / dweeks(1)
    ) %>%
    select(
        alf_e,
        age_3cat,
        dose1_date = vacc_dose1_date,
        dose1_name = vacc_dose1_name,
        dose2_date = vacc_dose2_date,
        dose2_name = vacc_dose2_name,
        dose2_diff = vacc_dose2_diff,
        doseb_date = vacc_doseb_date,
        doseb_name = vacc_doseb_name,
        doseb_diff = vacc_doseb_diff
    ) %>%
    pivot_longer(
        cols = matches("^dose"),
        names_to = c("dose_seq", ".value"),
        names_sep = "_"
    ) %>%
    rename(
        dose_date = date,
        dose_name = name,
        dose_diff = diff
    )

# freq tables ------------------------------------------------------------------

t_vacc_dose_freq <-
    d_vacc %>%
    count(dose_name, dose_seq) %>%
    mutate(n = round(n, -1)) %>%
    pivot_wider(names_from = dose_seq, values_from = n) %>%
    mutate(dose_name = replace_na(dose_name, "Not received"))

t_age_vacc_dose_freq <-
    d_vacc %>%
    count(age_3cat, dose_name, dose_seq) %>%
    mutate(n = round(n, -1)) %>%
    pivot_wider(names_from = dose_seq, values_from = n) %>%
    mutate(dose_name = replace_na(dose_name, "Not received"))

t_vacc_dose_pattern <-
    d_cohort %>%
    count(
        vacc_dose1_name,
        vacc_dose2_name,
        vacc_doseb_name
    ) %>%
    mutate(
        n = round(n, -1),
        p = round(n / sum(n) * 100, 1)
    ) %>%
    arrange(desc(n))

# time between doses -----------------------------------------------------------

p_vacc_diff <-
    d_vacc %>%
    filter(not_na(dose_diff)) %>%
    mutate(
        dose_seq  = str_c(dose_seq, "-", dose_name),
        dose_diff = floor(dose_diff)
    ) %>%
    count(dose_seq, dose_diff) %>%
    mutate(
        n = if_else(between(n, 1, 9), as.integer(10), n)
    )  %>%
    filter(dose_seq %in% c(
        "dose2-AZ",
        "dose2-PB",
        "doseb-MD",
        "doseb-PB"
    )) %>%
    ggplot(aes(x = dose_diff, y = n)) +
    facet_wrap(~dose_seq, ncol = 1, scales = "free_y") +
    geom_line() +
    scale_y_continuous(labels = comma)

print(p_vacc_diff)

# number of vaccinations over time, by type and dose ---------------------------

p_vacc_week <-
    d_vacc %>%
    filter(not_na(dose_date)) %>%
    mutate(dose_date = floor_date(dose_date, "week")) %>%
    count(dose_name, dose_seq, dose_date, name = "n_per_week") %>%
    mutate(
        n_per_week = if_else(between(n_per_week, 1, 9), as.integer(10), n_per_week)
    ) %>%
    ggplot(aes(x = dose_date, y = n_per_week, colour = dose_seq)) +
    facet_wrap(~dose_name, ncol = 1) +
    geom_line() +
    scale_x_date(
        name = "",
        date_breaks = "2 months",
        date_labels = "%b\n%y"
    ) +
    scale_colour_brewer(palette = "Set1") +
    theme(legend.position = "bottom")

print(p_vacc_week)

# Infections ===================================================================
cat("Infections\n")

d_infect <-
    d_cohort %>%
    select(alf_e, age_3cat, matches("^infection[0-9]_")) %>%
    rename_with(
        .cols = matches("^infection"),
        .fn = str_replace,
        pattern = "infection(.+)_date",
        replacement = "\\1"
    ) %>%
    pivot_longer(
        cols = c(-alf_e, -age_3cat),
        names_to = c("infection_seq", "event_type"),
        names_sep = "_",
        values_to = "event_date"
    ) %>%
    filter(not_na(event_date)) %>%
    filter(event_date >= ymd("2020-03-01")) %>%
    select(-infection_seq)

# plot infections, admissions and deaths by age ================================

d_hosp_death <-
    d_cohort %>%
    mutate(event_date = pmin(hosp_admis_date, death_covid_date, na.rm = TRUE)) %>%
    select(alf_e, age_3cat, event_date) %>%
    filter(not_na(event_date)) %>%
    filter(event_date >= ymd("2020-03-01")) %>%
    mutate(event_type = "hosp_death")

d_infect_hosp_death <-
    bind_rows(d_infect, d_hosp_death)

p_infect_hosp_death <-
    d_infect_hosp_death %>%
    mutate(event_week = floor_date(event_date, "week")) %>%
    count(event_type, event_week) %>%
    mutate(event_type = factor(event_type, c("test", "hosp_death"))) %>%
    mutate(n = if_else(between(n, 1, 9), as.integer(10), n)) %>%
    ggplot(aes(x = event_week, y = n)) +
    facet_wrap(~event_type, ncol = 1, scales = "free_y") +
    geom_col() +
    # formatting
    scale_x_date(
        name = "",
        date_breaks = "2 months",
        date_labels = "%b\n%y"
    ) +
    scale_y_continuous(
        labels = comma
    )
print(p_infect_hosp_death)


# Infection relative to dose ===================================================
cat("Infection relative to dose\n")

d_s <-
    d_cohort %>%
    select(
        alf_e,
        vacc_dose1_date,
        vacc_dose2_date,
        vacc_doseb_date
    )

d_infect2dose <-
    d_infect %>%
    filter(event_type == "test") %>%
    left_join(d_s, by = "alf_e") %>%
    filter(not_na(vacc_dose1_date)) %>%
    mutate(
        infect2dose1 = case_when(
            is.na(vacc_dose2_date) ~ interval(vacc_dose1_date, event_date) / dweeks(),
            event_date < vacc_dose2_date ~ interval(vacc_dose1_date, event_date) / dweeks()
        ),
        infect2dose2 = case_when(
            is.na(vacc_dose2_date) ~ NA_real_,
            vacc_dose1_date < event_date & is.na(vacc_doseb_date) ~ interval(vacc_dose2_date, event_date) / dweeks(),
            vacc_dose1_date < event_date & event_date < vacc_doseb_date ~ interval(vacc_dose2_date, event_date) / dweeks()
        ),
        infect2doseb = case_when(
            is.na(vacc_doseb_date) ~ NA_real_,
            vacc_dose2_date < event_date ~ interval(vacc_doseb_date, event_date) / dweeks()
        )
    ) %>%
    select(
        alf_e,
        matches("infect2dose.")
    ) %>%
    pivot_longer(
        cols = matches("infect2dose."),
        names_to = "dose",
        values_to = "time"
    ) %>%
    mutate(
        dose = str_replace(dose, "infect2", "")
    ) %>%
    filter(not_na(time)) %>%
    mutate(time = floor(time))

p_infect2dose <-
    d_infect2dose %>%
    count(dose, time) %>%
    mutate(n = if_else(between(n, 1, 9), as.integer(10), n)) %>%
    ggplot(aes(x = time, y = n)) +
    facet_wrap(~dose, ncol = 1) +
    geom_col() +
    scale_x_continuous(
        limits = c(-50, 50),
        breaks = pretty_breaks()
    )

# Save =========================================================================
cat("Save\n")

qsavem(
    t_vacc_dose_freq,
    t_vacc_dose_pattern,
    t_age_vacc_dose_freq,
    p_vacc_week,
    p_vacc_diff,
    p_infect_hosp_death,
    p_infect2dose,
    file = "results/cohort_freq_plots.qsm"
)
