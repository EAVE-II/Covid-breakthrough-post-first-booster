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

library(dplyr)
library(data.table)
library(ggplot2)
library(lubridate)
library(survival)
library(qs)
library(table1)

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


summary_factorlist_wt <- function(data, dependent, explanatory) {
  summary_tbl_list <- list()
  
  for (i in 1:length(explanatory)) {
    if (is.numeric(data[, get(explanatory[i])])) {
      z_mean <- data[, .(mean = round(weighted.mean(get(explanatory[i]), w = weight, na.rm = TRUE), 1)), get(dependent)]
      z_mean_sd <- data[, round(sqrt(spatstat.geom::weighted.var(get(explanatory[i]), w = weight)), 1), get(dependent)]
      names(z_mean_sd) <- c("get", "sd")
      z_mean <- as.data.table(merge(z_mean, z_mean_sd, by = "get", all = TRUE))
      z_mean[, mean.sd := paste0(mean, " (", sd, ")")]
      z_mean <- z_mean[, .(get, mean.sd)]
      z_mean$characteristic <- explanatory[i]
      z_mean$measure <- "mean.sd"
      names(z_mean) <- c("dependent", "value", "characteristic", "measure")
      
      z_median <-
        data[, .(
          median = spatstat.geom::weighted.median(get(explanatory[i]), w = weight),
          q1 = spatstat.geom::weighted.quantile(get(explanatory[i]), w = weight, probs = 0.25),
          q3 = spatstat.geom::weighted.quantile(get(explanatory[i]), w = weight, probs = 0.75)
        ), get(dependent)
        ]
      z_median$characteristic <- explanatory[i]
      z_median[, median.iqr := paste0(round(median, 1), " (", round(q3 - q1, 1), ")")]
      z_median$measure <- "median.iqr"
      z_median <- z_median[, .(get, median.iqr, characteristic, measure)]
      names(z_median) <- c("dependent", "value", "characteristic", "measure")
      
      summary_tbl_list[[i]] <- rbind(z_mean, z_median)
      
    } else {
      
      count_by_group <- data[, .(n = sum(weight)), .(get(explanatory[i]), dependent = get(dependent))]
      count_by_group[, percent := sprintf("%.1f", round(100 * n / sum(n), 1)), dependent]
      count_by_group[, n_perc := paste0(round(n, 0), " (", percent, "%)")]
      count_by_group[, measure := "n (%)"]
      count_by_group <- count_by_group[, .(dependent, value = n_perc, characteristic = paste0(explanatory[i], ": ", get), measure)]
      summary_tbl_list[[i]] <- count_by_group
    }
  } 
  
  summary_tbl_wt <- rbindlist(summary_tbl_list)
  summary_tbl_wt
}

