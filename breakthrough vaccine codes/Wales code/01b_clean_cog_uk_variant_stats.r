source("r_clear_and_load.r")

# Load =========================================================================
cat("load\n")

d_variant_testing <-
    read_csv(
        file = "cog-uk-variant-stats/2022-01-13_lin_ew_wales_barplot.csv",
        col_types = cols()
    ) %>%
    clean_names() %>%
    rename(lineage = x1)

# Clean ========================================================================
cat("Clean\n")

# epi_week 9 is week starting Monday 24th Feb 2020

d_variant_testing <-
    d_variant_testing %>%
    pivot_longer(
        cols = -lineage,
        names_to = "week_epi",
        values_to = "n"
    ) %>%
    mutate(
        week_epi  = str_replace(week_epi, "x", ""),
        week_epi  = as.numeric(week_epi),
        week_date = ymd("2020-02-24") + dweeks(week_epi - min(week_epi))
    ) %>%
    mutate(
    variant = lineage %>% factor() %>%
        fct_collapse(
            "Alpha"       = "B.1.1.7",
            "Delta"       = "B.1.617.2",
            "Delta"       = "AY.4",
            "Delta"       = "AY",
            "Omicron"     = "B.1.1.529",
            "Omicron"     = "BA.1",
            "Omicron"     = "BA.2",
            "Omicron"     = "BA.3",
            "Summer_2020" = "B.1.177",
            "Other B"     = "B",
            "(Other)"     = c(
                "A",
                "AZ.5",
                "B.1.351",
                "B.1.621",
                "B.1.630",
                "B.1.640",
                "C.1.2",
                "C.37",
                "B.1.427",
                "B.1.429",
                "B.1.525",
                "B.1.526",
                "B.1.617.1",
                "P.1",
                "P.2",
                "P.3",
                "Other"
        )
        ) %>%
        fct_relevel(
            "Other B",
            "Summer_2020",
            "Alpha",
            "Delta",
            "Omicron",
            "(Other)"
        ) %>%
        fct_rev()
    ) %>%
    select(
        week_epi,
        week_date,
        lineage,
        variant,
        n
    )

# Save =========================================================================
cat("Save\n")

qsave(
    d_variant_testing,
    file = s_drive("d_variant_testing.qs")
)
