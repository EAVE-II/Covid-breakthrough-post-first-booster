cat("Clearing workspace and loading stuff\n")

# clear workspace ==============================================================

rm(list = ls())
gc()

# options ======================================================================

options(
  dplyr.summarise.inform = FALSE,
  readr.show_col_types = FALSE,
  lubridate.week.start = 1 # Monday
)

## Packages
# 
# library(dplyr)
# library(ggplot2)
# library(lubridate)
# library(survival)
# library(qs)
# library(table1)

# load packages ================================================================

pkgs <- c(
  # load these first to not highjack the other packages
  "purrr",
  "MASS",
  "mgcv",
  # alphabetical
  "assertr",
  "beepr",
  "broom",
  "changepoint.np",
  "dplyr",
  "dtplyr",
  "egg",
  "forcats",
  "ggfortify",
  "ggplot2",
  "ggthemes",
  "ggrepel",
  "gtsummary",
  "knitr",
  "kableExtra",
  "mice",
  "janitor",
  "lubridate",
  "patchwork",
  "qs",
  "RcppRoll",
  "rlang",
  "rmarkdown",
  #"sailr",
  "scales",
  "speedglm",
  "stringr",
  "readr",
  "survey",
  "survival",
  "table1",
  "tibble",
  "tidyquant",
  "tidyr"
)

for (pkg in pkgs) {
  suppressWarnings(
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  )
}


# useful dates =================================================================

vacc_start_date    <- ymd("2020-12-08")
booster_start_date <- ymd("2021-09-16")
study_end_date     <- ymd("2022-02-17")

cohort_age_date = as.Date("2020-11-01")  # one month earlier because of rounding down of dates of birth - i.e. we had month of birth only
## Qu - protocol is age @7 Dec 2020; just pointing out but doesn't likely  make
## much difference   ## we only have month of birth

pfizer_available_community_date <- as.Date("2021-02-01")

na_to_zero <- function(DT, cols) {
  for (j in cols)
    set(DT, which(is.na(DT[[j]])), j, 0)
}

calc_age <- function(dob, refdate) {
  period <- as.period(interval(dob, refdate), unit = "year")
  period$year
}

age.cat <- function(x, lower = 0, upper, by = agebandsize, sep = "-", above.char = "+") {
  labs <- c(paste(seq(lower, upper - by, by = by),
                  seq(lower + by - 1, upper - 1, by = by),
                  sep = sep),
            paste(upper, above.char, sep = ""))
  
  cut(floor(x), breaks = c(seq(lower, upper, by = by), Inf), 
      right = FALSE, labels = labs)
}

# custom functions =============================================================

KableOne <- function(x, ...) {
  k1 <- print(x, quote = TRUE, printToggle = FALSE, ...)
  rownames(k1) <- gsub(pattern = '\\"', replacement = '', rownames(k1))
  colnames(k1) <- gsub(pattern = '\\"', replacement = '', colnames(k1))
  return(k1)
}

qcut <- function(x, n) {
  quantiles = seq(0, 1, length.out = n+1)
  cutpoints = unname(quantile(x, quantiles, na.rm = TRUE))
  as.character(cut(x, cutpoints, include.lowest = TRUE))
}


calc_pyears <- function(.data, ...) {
  .data %>%
    dtplyr::lazy_dt() %>%
    group_by(...) %>%
    summarise(
      pyears = sum((tstop - tstart) * sample_weight) / 365.25,
      event = sum(event_flg)
    ) %>%
    ungroup() %>%
    as_tibble()
}

calc_pyears2 <- function(.data, ...) {
  .data %>%
    dtplyr::lazy_dt() %>%
    group_by(...) %>%
    summarise(
      pyears = sum(tstop - tstart) / 365.25,
      event = sum(event_flg)
    ) %>%
    ungroup() %>%
    as_tibble()
}

# plot dimensions ==============================================================

p_width  <- 5.2
p_height <- 8.75



