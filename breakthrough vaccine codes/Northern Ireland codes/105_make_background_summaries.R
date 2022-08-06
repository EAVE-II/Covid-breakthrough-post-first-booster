cat("Clearing workspace and loading stuff\n")

source("scripts/000_libraries_and_functions.R")

# load =========================================================================
cat("load\n")

d_cohort_raw <- qread("input/d_cohort_raw.qs") # from 102

d_variant_testing <- qread("input/d_variant_testing.qs") #from 104 COG data 

# De-duplicate the infections day using a 90 day window

## First create a wide dataset by infection number 
# d_cohort_raw_wide = d_cohort_raw %>%
#     #First create an indicator that identifies infection number
#     mutate(type = "infection") %>% collect() %>%
#     group_by(study_id) %>%
#     mutate(n =  1:n()) %>%
#     mutate(infection_n = paste(type, n, sep ="_")) %>%
#     #And then pivot wider to create infections per column 
#     tidyr::pivot_wider(
#         id_cols = study_id,
#         names_from = infection_n,
#         values_from = specimen_date
#     )
# 
# ## Now want to compare each infection to preceding and if its <90 days keep the previous 
# ## Can convert this to a function when I have time
# test = d_cohort_raw_wide %>%
#     ## Compare Infection 1 and 2 - and keep the first
#     mutate(btn1.2 = as.numeric(difftime(infection_2, infection_1, units = "days"))) %>%
#     mutate(clean1 = ifelse(btn1.2<90, as.character(infection_1), as.character(infection_2))) %>%
#     ## Compare infection 3 to previous and change if >=90 days 
#     mutate(btn3.c1 = as.numeric(difftime(infection_3, clean1, units = "days"))) %>%
#     mutate(clean2 = ifelse(btn3.c1<90, as.character(clean1), as.character(infection_3))) %>%
#     ## Compare infection 4 to previous and change if >=90 days 
#     mutate(btn4.c2 = as.numeric(difftime(infection_4, clean2, units = "days"))) %>%
#     mutate(clean3 = ifelse(btn4.c2<90, as.character(clean2), as.character(infection_4))) %>%
#     ## Compare infection 5 to previous and change if >=90 days 
#     mutate(btn5.c3 = as.numeric(difftime(infection_5, clean3, units = "days"))) %>%
#     mutate(clean4 = ifelse(btn5.c3<90, as.character(clean3), as.character(infection_5))) %>%
#     # Compare infection 6 to previous and change if >=90 days 
#     mutate(btn6.c4 = as.numeric(difftime(infection_6, clean4, units = "days"))) %>%
#     mutate(clean5 = ifelse(btn6.c4<90, as.character(clean4), as.character(infection_6))) %>%
#     # Compare infection 7 to previous and change if >=90 days 
#     mutate(btn7.c5 = as.numeric(difftime(infection_7, clean5, units = "days"))) %>%
#     mutate(clean6 = ifelse(btn7.c5<90, as.character(clean5), as.character(infection_7))) %>%
#     # Compare infection 8 to previous and change if >=90 days 
#     mutate(btn8.c6 = as.numeric(difftime(infection_8, clean6, units = "days"))) %>%
#     mutate(clean7 = ifelse(btn8.c6<90, as.character(clean6), as.character(infection_8))) %>%
#     # Compare infection 9 to previous and change if >=90 days 
#     mutate(btn9.c7 = as.numeric(difftime(infection_9, clean7, units = "days"))) %>%
#     mutate(clean8 = ifelse(btn9.c7<90, as.character(clean7), as.character(infection_9))) %>%
#     # Compare infection 10 to previous and change if >=90 days 
#     mutate(btn10.c8 = as.numeric(difftime(infection_10, clean8, units = "days"))) %>%
#     mutate(clean9 = ifelse(btn10.c8<90, as.character(clean8), as.character(infection_10))) %>%
#     # Compare infection 11 to previous and change if >=90 days 
#     mutate(btn11.c9 = as.numeric(difftime(infection_11, clean9, units = "days"))) %>%
#     mutate(clean10 = ifelse(btn11.c9<90, as.character(clean9), as.character(infection_11))) %>%
#     # Compare infection 12 to previous and change if >=90 days 
#     mutate(btn12.c10 = as.numeric(difftime(infection_12, clean10, units = "days"))) %>%
#     mutate(clean11 = ifelse(btn12.c10<90, as.character(clean10), as.character(infection_12))) %>%
#     # Compare infection 13 to previous and change if >=90 days 
#     mutate(btn13.c11 = as.numeric(difftime(infection_13, clean11, units = "days"))) %>%
#     mutate(clean12 = ifelse(btn13.c11<90, as.character(clean11), as.character(infection_13))) 
# 
# ## had to use as.character to get ifelse to work, now change back to date
# 
# date_cols = c("clean1","clean2", "clean3","clean4", "clean5", "clean6", "clean7",
#               "clean8", "clean9", "clean10", "clean11", "clean12")
# 
# test = test %>% mutate_at(vars(date_cols), as.Date, format="%Y-%m-%d")
# 
# ## Now prepare to pivot this back and clean
# d_cohort_dedup = test %>%
#     select(study_id, infection_1, clean1, clean2, clean3,clean4, clean5, clean6, clean7,
#            clean8, clean9, clean10, clean11, clean12) %>%
#     rename(clean0 = infection_1) %>%
#     tidyr::pivot_longer(
#         cols = starts_with('clean'),
#         names_to = "infection",
#         values_to = "specimen_date"
#     ) %>%
#     filter(!is.na(specimen_date)) %>%
#     select(-infection) %>%
#     distinct(study_id, specimen_date) %>%
#     mutate(type = "infection") %>%
#     group_by(study_id) %>%
#     mutate(n =  1:n()) %>%
#     mutate(infection_n = paste(type, n, sep ="_")) %>%
#     select(-c(type,n))
# 
# ## Check nothing outrageous
# #z = d_cohort_dedup %>%
#  #   group_by(infection_n) %>%
#   #  tally()
# # two people have had four infections
# 
# ## Save ==============================================================
# ## Save the de-dupl test dataset to save not running again
# cat("Save\n")
# 
# qsave(
#   d_cohort_dedup,
#   file = ("input/d_pos_dedupl.qs")
# )

### Load data ==============================================================

d_cohort_dedup <- qread("input/d_pos_dedupl.qs") %>%
  select(-infection_n)

# Count week infections ====================================================
cat("Count weekly infections\n")

d_week_infection <-
    d_cohort_dedup %>%
    mutate(infection_date = specimen_date) %>%
    mutate(infection_date = floor_date(infection_date, "week")) %>%
    group_by(infection_date) %>%
    tally() %>% rename (infection_n = n) %>%
  # changed from vacc start date to booster start date
  # reverted back on 3 May
    filter(infection_date >= floor_date(vacc_start_date, "week")) #counts from start of vacc start date

# Find changepoints ============================================================
cat("Find changepoints\n")

fit_cpt <- changepoint::cpt.meanvar(
    data      = d_week_infection$infection_n,
    method    = "PELT", # "AMOC", "PELT", "SegNeigh", "BinSeg"
    penalty   = "BIC"   # "None", "SIC", "BIC", "MBIC", AIC", "Hannan-Quinn", "Asymptotic", "Manual", "CROPS"
)

#plot(fit_cpt)

infection_changepoints <- d_week_infection$infection_date[changepoint::cpts(fit_cpt)]

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
    dplyr::select(
        week_date,
        variant,
        variant_pop = prop
    ) %>%
    # sanity check: we only have one dominant variant per week
    #testit::assert(is_uniq, week_date) %>%
    ## this wasn't working so test manually
    #z = d_dominant_variant %>%
    #group_by(week_date)
    # 102 observations - same as d_dominant_variant so must be unique
    #pick first time variant was dominant
    group_by(variant) %>%
    summarise(first_week = min(week_date)) %>%
    ungroup() %>%
    arrange(first_week) %>%
    # pick variants one before vaccination started and all after
    mutate(
      # changed from vacc start date to booster start date
      # reverted back on 3 May as it was causing issues further on 
        first_week_diff = interval(vacc_start_date, first_week) / ddays(1),
        keep = if_else(first_week_diff > 0, NA_real_, first_week_diff),
        keep = max(keep, na.rm = TRUE),
        keep = first_week_diff >= keep
    ) %>%
    filter(keep) %>%
    dplyr::select(-first_week_diff, -keep) %>%
    # reset min first_week to be week vaccinations started
    mutate(
        first_week = if_else(
            condition = first_week == min(first_week),
            # changed from vacc start date to booster start date
            # reverted back 3 May
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
    data.frame(infection_date = ymd(x_weeks)) %>% #create all the weeks
    left_join(d_week_infection, by = "infection_date") %>% #join to infection and fill missing weeks with 0
    mutate(infection_n = tidyr::replace_na(infection_n, 0)) %>%
    mutate(bg_start = infection_date) %>%
    left_join(d_dominant_variant, by = c("bg_start" = "first_week")) %>%
    rename(bg_variant = variant) %>%
    tidyr::fill(bg_variant, .direction = "downup") %>%
    arrange(bg_start) %>% #sort by week from start of study (7 Dec 2020)
    mutate(
        bg_interval      = stringr::str_c("bg_interval", row_number()), #concat bg_int with the row number which is the week number
        bg_end           = lead(bg_start) - ddays(1), #last date of the week 
        bg_end           = tidyr::replace_na(bg_end, max(d_week_infection$infection_date)), #replaces missing but we didn't have any
        bg_length        = interval(bg_start, (bg_end + ddays(1)) ) / dweeks(1),
        bg_avg_infection = infection_n / bg_length,
        bg_end_week      = floor_date(bg_end, "week")
    ) %>%
    dplyr::select(
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
    ggrepel::geom_label_repel(
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
        breaks = scales::pretty_breaks(),
        labels = scales::comma
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
    file = ("input/d_background.qs")
)

qsave(
    p_bg,
    file = "results/p_bg.qs"
)

# print ========================================================================

print(d_background)
print(p_bg)
