# Read libraries and fuctions  ==============================================================

source("scripts/000_libraries_and_functions.R")

# load =========================================================================
cat("load\n")

dt_cog_uk <- readRDS("input/Cog_UK_Feb22_1.rds")
 
## Need to check with Declan about lineage mapping 

d_variant_testing <- dt_cog_uk %>%
    rename(lineage = Lineage) %>%
    mutate(lineage = stringr::str_replace_all(lineage, " ", "")) %>% collect()

variant = read.csv("input/VOCVUI_Master_Lookup.csv") %>%
    rename(lineage = Lineage) %>%
    dplyr::select(lineage, variant) %>%
    mutate(lineage = stringr::str_replace_all(lineage, " ", ""))

# Clean ========================================================================

# epi_week 9 is week starting Monday 24th Feb 2020

d_variant_testing <- d_variant_testing %>%
    mutate(week_epi = isoweek(ymd(SpecimenDate))) %>%
        filter(!is.na(week_epi)) %>% # remove records with no spec date
    mutate(
        week_epi  = stringr::str_replace(week_epi, "x", ""),
        week_epi  = as.numeric(week_epi),
        week_date = lubridate::floor_date(as.Date(SpecimenDate), "week")
        ) %>%
    dplyr::select(-c(Study_ID, SpecimenDate, Anon_Specimen_ID)) %>%
    left_join(variant, by = c(lineage = "lineage")) %>%
    group_by(lineage,variant,week_epi, week_date) %>%
    tally() %>%
       dplyr::select(
        week_epi,
        week_date,
        lineage,
        variant,
        n
    ) %>% 
    filter(week_date>=booster_start_date-7) %>% #took the week before to capture the start of the vaccine programme
    arrange(week_date)

# updated this from vacc start date (dec 2020) to booster start date Sep 21 30 April


# Save =========================================================================
cat("Save\n")

qsave(
    d_variant_testing,
    file = ("input/d_variant_testing.qs")
)


# Plot =========================================================================
cat("Plot\n")

lkp_fill <- c(
    "#444444",
    "#377eb8",
    "#ff7f00",
    "#984ea3",
    "#4daf4a",
    "#e41a1c",
    "#4EEE94",
    "#FF82AB"
)

p_variant_week_freq <-
    d_variant_testing %>%
    ggplot(aes(x = week_date, y = n, fill = variant)) +
    geom_col() +
    scale_fill_manual(
        values = lkp_fill
    ) +
    theme(
        legend.position = c(0,1),
        legend.justification =  c(0, 1)
    )

p_variant_week_prop <-
    d_variant_testing %>%
    group_by(week_date) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = week_date, y = prop, fill = variant)) +
    geom_col() +
    scale_fill_manual(
        values = lkp_fill
    ) +
    theme(
        legend.position = c(0,1),
        legend.justification =  c(0, 1)
    )

print(p_variant_week_prop)


# Save =========================================================================
cat("Save\n")

qsavem(
    p_variant_week_freq,
    p_variant_week_prop,
    file = "results/variant.qsm"
)

message("***Variant information added and aggregated by week***")
