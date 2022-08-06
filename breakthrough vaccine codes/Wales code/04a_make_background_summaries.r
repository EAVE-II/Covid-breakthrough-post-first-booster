source("r_clear_and_load.r")

# load =========================================================================
cat("load\n")

d_cohort_raw <- qread(s_drive("d_cohort_raw.qs"))

d_variant_testing <- qread(s_drive("d_variant_testing.qs"))


# Count week infections ========================================================
cat("Count weekly infections\n")

d_week_infection <-
    d_cohort_raw %>%
    select(
        infection1_test_date,
        infection2_test_date,
        infection3_test_date,
        infection4_test_date
    ) %>%
    pivot_longer(cols = everything(), values_to = "infection_date") %>%
    filter(not_na(infection_date)) %>%
    mutate(infection_date = floor_date(infection_date, "week")) %>%
    count(infection_date, name = "infection_n") %>%
    filter(infection_date >= floor_date(vacc_start_date, "week"))


# Find changepoints ============================================================
cat("Find changepoints\n")

fit_cpt <- cpt.meanvar(
    data      = d_week_infection$infection_n,
    method    = "PELT", # "AMOC", "PELT", "SegNeigh", "BinSeg"
    penalty   = "BIC"   # "None", "SIC", "BIC", "MBIC", AIC", "Hannan-Quinn", "Asymptotic", "Manual", "CROPS"
)

plot(fit_cpt)

infection_changepoints <- d_week_infection$infection_date[cpts(fit_cpt)]

# Calculate dominant variant ===================================================
cat("Calculate dominant variant\n")

d_dominant_variant <-
    d_variant_testing %>%
    mutate(variant = as.character(variant)) %>%
    group_by(week_date) %>%
    mutate(
        prop = n / sum(n),
        is_max = prop == max(prop)
    ) %>%
    ungroup() %>%
    # pick dominant variant of the week
    filter(is_max) %>%
    select(
        week_date,
        variant,
        variant_pop = prop
    ) %>%
    # sanity check: we only have one dominant variant per week
    assert(is_uniq, week_date) %>%
    # pick first time variant was dominant
    group_by(variant) %>%
    summarise(first_week = min(week_date)) %>%
    ungroup() %>%
    arrange(first_week) %>%
    # pick variants one before vaccination started and all after
    mutate(
        first_week_diff = interval(vacc_start_date, first_week) / ddays(1),
        keep = if_else(first_week_diff > 0, NA_real_, first_week_diff),
        keep = max(keep, na.rm = TRUE),
        keep = first_week_diff >= keep
    ) %>%
    filter(keep) %>%
    select(-first_week_diff, -keep) %>%
    # reset min first_week to be week vaccinations started
    mutate(
        first_week = if_else(
            condition = first_week == min(first_week),
            true = floor_date(vacc_start_date, "week"),
            false = first_week
        )
    )


# Combine and calculate background rate within each interval ===================
cat("Combine and calculate rate\n")

x_weeks <- seq(
    from = min(d_week_infection$infection_date),
    to   = max(d_week_infection$infection_date),
    by   = "week"
)

d_background <-
    data.frame(infection_date = ymd(x_weeks)) %>%
    left_join(d_week_infection, by = "infection_date") %>%
    mutate(infection_n = replace_na(infection_n, 0)) %>%
    mutate(bg_start = infection_date) %>%
    left_join(d_dominant_variant, by = c("bg_start" = "first_week")) %>%
    rename(bg_variant = variant) %>%
    fill(bg_variant, .direction = "downup") %>%
    arrange(bg_start) %>%
    mutate(
        bg_interval      = str_c("bg_interval", row_number()),
        bg_end           = lead(bg_start) - ddays(1),
        bg_end           = replace_na(bg_end, max(d_week_infection$infection_date)),
        bg_length        = interval(bg_start, (bg_end + ddays(1)) ) / dweeks(1),
        bg_avg_infection = infection_n / bg_length,
        bg_end_week      = floor_date(bg_end, "week")
    ) %>%
    select(
        bg_interval,
        bg_start,
        bg_end,
        bg_end_week,
        bg_variant,
        bg_avg_infection
    )


# plot =========================================================================
cat("plot\n")

p_bg <-
    d_week_infection %>%
    ggplot(aes(x = infection_date, y = infection_n)) +
    geom_col() +
    # variant lines
    geom_vline(
        data = d_dominant_variant,
        mapping = aes(xintercept = first_week),
        linetype = 2,
        colour = "blue",
        size = 1
    ) +
    geom_label_repel(
        data = d_dominant_variant,
        mapping = aes(x = first_week, y = 70000, label = variant),
        direction = "y",
        hjust = 0,
        nudge_x = 5,
        min.segment.length = 9000
    ) +
    # background intervals
    geom_segment(
        data = d_background,
        mapping = aes(
            x = bg_start, xend = bg_end_week,
            y = bg_avg_infection, yend = bg_avg_infection
        ),
        colour ="red",
        size = 1
    ) +
    # formatting
    scale_x_date(
        date_breaks = "1 month",
        date_labels = "%b\n%Y"
    ) +
    scale_y_continuous(
        breaks = pretty_breaks(),
        labels = comma
    ) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
p_bg

# save =========================================================================
cat("save\n")

qsave(
    d_background,
    file = s_drive("d_background.qs")
)

qsave(
    p_bg,
    file = "results/p_bg.qs"
)

# print ========================================================================

print(d_background)
print(p_bg)
