cat("Clearing workspace and loading stuff\n")

source("scripts/000_libraries_and_functions.R")

x_seed <- 125460070

# Load =========================================================================
cat("Load\n")

d_cohort <- qread("input/dt_cohort_clean.qs") # from step 103

d_background <- qread("input/d_background.qs") #  from step 105

# Clearning======================================================================
# Remove the NAs and AZ for booster dose

# d_cohort %>% group_by(Booster_name) %>% tally()

d_cohort = d_cohort %>% filter(!is.na(Booster_name) & Booster_name!="AZ")

# Select analysis sample =======================================================
cat("Select analysis sample\n")
#z = d_cohort %>%
 # group_by(has_outcome)%>%
  #tally() %>%
  #collect()

d_cohort <-
  d_cohort %>%
  mutate(
    epoch = Dose2_date + ddays(14), #field with value 14 days after dose date
    # outcome: covid-19 hosp or death
    outcome_date = pmin(covid_admission_date, death_covid_date, na.rm = TRUE)) %>% collect() %>%
  mutate(
    # has outcome before any other event which would stop follow-up
    has_outcome = epoch <= outcome_date &
      (is.na(death_noncovid_date) | outcome_date <= death_noncovid_date) &
      (is.na(Dose4_date)       | outcome_date <= Dose4_date) &
      outcome_date <= study_end_date) %>% 
  mutate(
    has_outcome = if_else(is.na(outcome_date), FALSE, has_outcome, has_outcome)) %>%
  mutate(
    # sample selection criteria
    has_dose2 = assertr::not_na(Dose1_date) &
      assertr::not_na(Dose2_date) &
      (is.na(death_noncovid_date) | epoch <= death_noncovid_date) &
      (is.na(Dose4_date) | epoch <= Dose4_date) &
      epoch <= study_end_date &
      Dose2_name %in% c("AZ", "MD", "PB"),
    no_prv_outcome = is.na(outcome_date) | (epoch <= outcome_date),
    valid_doseb = is.na(Booster_date) | Booster_name %in% c("MD", "PB")
  ) 
  
summary(d_cohort)

d_cohort %>% count(has_dose2, has_outcome, no_prv_outcome)
## all of the has_dose2 == FALSE have epoch > than the end of the study
## period i.e. they received their second dose close to the end of the
## study period or after it

t_analysis_sample_selection <- tribble(
    ~step, ~description, ~total_n, ~outcome_n,
        0, "All cohort",
            d_cohort %>% nrow(),
            d_cohort %>% filter(has_outcome) %>% nrow(),
        1, "Has second dose",
            d_cohort %>% filter(has_dose2) %>% nrow(),
            d_cohort %>% filter(has_dose2) %>% filter(has_outcome) %>% nrow(),
        2, "No outcome prior to second dose",
            d_cohort %>% filter(has_dose2, no_prv_outcome) %>% nrow(),
            d_cohort %>% filter(has_dose2, no_prv_outcome) %>% filter(has_outcome) %>% nrow(),
        3, "If had a booster, booster is MD or PB",
            d_cohort %>% filter(has_dose2, no_prv_outcome, valid_doseb) %>% nrow(),
            d_cohort %>% filter(has_dose2, no_prv_outcome, valid_doseb) %>% filter(has_outcome) %>% nrow()
)

t_analysis_sample_selection <-
    t_analysis_sample_selection %>%
    mutate(total_n_diff = total_n - lag(total_n), .after = total_n) %>%
    mutate(outcome_n_diff = outcome_n - lag(outcome_n), .after = outcome_n)

print(t_analysis_sample_selection)

## drop all those without valid data 
## takes us to 1166639 
## using the dataset from 04a (d_pop_cohort rather than dt_cohort_clean) which is clean so can remove this step 
#d_cohort = d_cohort %>% filter(has_dose2, no_prv_outcome, valid_doseb)

# Make analysis sample =========================================================
cat("Make analysis sample\n")

d_case <-
    d_cohort %>%
    filter(has_dose2, no_prv_outcome, valid_doseb) %>%
    filter(has_outcome) %>%
    mutate(
        cc_group = "case",
        sample_weight = 1
    )

set.seed(x_seed)

control_total_n <-
    d_cohort %>%
    filter(has_dose2, no_prv_outcome, valid_doseb) %>%
    filter(!has_outcome) %>%
    count() %>%
    unlist()

d_control <-
    d_cohort %>%
    filter(has_dose2, no_prv_outcome, valid_doseb) %>%
    filter(!has_outcome) %>%
    slice_sample(n = 10 * nrow(d_case)) %>%
    mutate(
        cc_group = "control",
        sample_weight = control_total_n / (10 * nrow(d_case))
    )

d_sample <- bind_rows(d_case, d_control)

# quick count

t_cc_group_n <- d_sample %>% count(cc_group)

## rename variables to match Welsh script

d_sample = d_sample %>%
  rename(vacc_dose2_name = Dose2_name, vacc_doseb_name = Booster_name,
  vacc_dose2_date = Dose2_date, vacc_doseb_date = Booster_date)

# split out the columns into:
#   * baseline covariates: age, sex, qcovid, health board, ...
#   * events: outcome date, move out date, vaccination dates, ...
# we will then transform d_events into a start-stop dataset and
# eventually left join back on the baseline measures

d_baseline <-
    d_sample %>%
    dplyr::select(
    # person id
        study_id,
    # vaccine names
        vacc_dose2_name,
        vacc_doseb_name,
    # characteristics
        sex,
        age_gp,
        vacc_dose1_dose2_diff_cat,
        dose2_prior_infection_cat,
        nimdm,
        ur_combined,
        BNF_group,
        lgd,
        n_tests_gp,
      # study design
        cc_group,
        sample_weight
    )

d_events <-
    d_sample %>%
    dplyr::select(
    # person id
        study_id,
    # event dates
        outcome_date,
        death_noncovid_date,
        Dose4_date,
        study_end_date,
    # vaccination
        vacc_dose2_date,
        vacc_dose2_name,
        vacc_doseb_date,
        vacc_doseb_name,
    # study design
        cc_group
    ) %>% collect()


# convert interval dates to days ===============================================
cat("Convert interval dates to days\n")

# convert dates to days since start of follow-up and determine the event
# also make sure everything is valid with cc_group

d_events <-
    d_events %>%
    mutate(
        # start date and week of follow-up
        start_date   = vacc_dose2_date + ddays(14),
        start_day    = start_date,
        start_week   = floor_date(start_day, "week"),
        # dates for post dose 2 cut points
        dose2_day014 = vacc_dose2_date + ddays(14),
        dose2_day070 = vacc_dose2_date + ddays(70), # dose 2: week  2 - 9
        dose2_day140 = vacc_dose2_date + ddays(140), # dose 2: week 10 - 19
        #dose2_day140 = vacc_dose2_date + ddays(140), # dose 2: week 24+
        # dates for post booster cut points
        doseb_day014 = vacc_doseb_date + ddays( 14),
        doseb_day035 = vacc_doseb_date + ddays(35), # dose B: week  2 - 4
        doseb_day056 = vacc_doseb_date + ddays(56), # dose B: days 5 - 7
        #doseb_day056 = vacc_doseb_date + ddays(56), # dose B: days 35 - 55
        #doseb_day056 = vacc_doseb_date + ddays( 112) # dose B: week  6
        #doseb_day056 = vacc_doseb_date + ddays( 56), # dose B: week  8
        #doseb_day070 = vacc_doseb_date + ddays( 70), # dose B: week 10
        #doseb_day084 = vacc_doseb_date + ddays( 84), # dose B: week 12
        #doseb_day098 = vacc_doseb_date + ddays( 98), # dose B: week 14
        #doseb_day112 = vacc_doseb_date + ddays(112)  # dose B: week 16
    ) %>%
    # dose2_day*_date columns should only have a value if the person had
    # not had reached day 14 of their booster dose at that time point
    mutate(across(
        .cols = matches("dose2_day[0-9]+"),
        .fns  = ~ if_else(.x >= doseb_day014, NA_Date_, .x, .x)
    )) %>%
    # stop follow-up if they have the outcome, 4th dose, die, or we reach
    # the end of the study, if multiple events on the same day, prioritise
    # the outcome
    mutate(
        stop_date = pmin(
            outcome_date,
            death_noncovid_date,
            Dose4_date,
            study_end_date,
            na.rm = TRUE
        ),
        stop_day = stop_date,
        event_cat = case_when(
            stop_date == outcome_date        ~ "covid19_hosp_death",
            stop_date == death_noncovid_date ~ "death_noncovid",
            stop_date == Dose4_date       ~ "fourth dose",
            stop_date == study_end_date      ~ "study_end"
        ),
        event_cat = factor(event_cat) %>% forcats::fct_relevel("study_end"),
        event_flg = as.numeric(event_cat == "covid19_hosp_death")
    ) %>%
    # convert dates to days from start of follow-up
    mutate(across(
        .cols = matches("(start_day|dose[2b]_day[0-9]+|stop_day)"),
        .fns  = ~ interval(start_day, .x) / ddays(1)
    )) %>%
    # this is a bit of a hack, but if someone has an event on the same
    # day as day 14 of their second dose, we assume we associate it with
    # day 14 rather than 0-13 days i.e. the vaccine is effective now and it
    # failed for this person
    mutate(
        stop_day = if_else(stop_day == 0, 0.5, stop_day)
    ) %>%
    # keep values for vaccination intervals that are within follow-up,
    # otherwise set to NA
    mutate(across(
        .cols = matches("^dose[2b]_day[0-9]+"),
        .fns = ~ if_else(.x < stop_day, .x, NA_real_, NA_real_)
    ))

## Can't get the verify to work so check manually
# sanity checks
#verify(
 #   start_day < stop_day
# z = d_events %>%
#   dplyr::select(start_day, stop_day) %>%
# mutate(flag = ifelse(start_day >=stop_day, 1,0)) %>%
#   count(flag)

## Looks OK
 
#    verify(
  #      (event_flg == 1 & cc_group == "case") |
  #     (event_flg == 0 & cc_group == "control")

# z = d_events %>%
#  group_by(event_flg, cc_group) %>%
#  tally() %>% collect()


# Add in background dates ======================================================
cat("Add background dates\n")

d_background <-
    d_background %>%
    dplyr::select(
        bg_interval,
        bg_start,
        bg_avg_infection,
        bg_variant
    )

# make a wide dataset of the intervals
d_background_wide <-
    d_background %>%
    dplyr::select(bg_interval, bg_start) %>%
    tidyr::pivot_wider(
        names_from  = bg_interval,
        values_from = bg_start
    )

# make sure background intervals only have values if they overlap with
# follow-up

#class(d_events$stop_date)

d_events <-
    cbind(d_events, d_background_wide) %>%
    # pivot first so it is easier to manipulate one column of intervals,
    # rather than X columns for X intervals
    tidyr::pivot_longer(
        cols = matches("^bg_interval"),
        names_to = "bg_interval",
        values_to = "bg_date"
    ) %>%
    mutate(
        # make sure no background periods are after follow-up
        bg_date = if_else(bg_date < stop_date, bg_date, NA_Date_, NA_Date_),
        # convert bg_interval columns from date to days from start_date for
        bg_day = interval(start_date, bg_date) / ddays(1)
    ) %>%
    # keep the background period just before follow-up starts, and all after
    mutate(bg_rank = if_else(bg_day > 0, NA_real_, bg_day)) %>%
    dtplyr::lazy_dt() %>%
    group_by(study_id) %>%
    mutate(bg_rank = max(bg_rank, na.rm = TRUE)) %>%
    ungroup() %>%
    as_tibble() %>%
    mutate(bg_day = if_else(bg_day < bg_rank, NA_real_, bg_day, bg_day)) %>%
    # negative bg_day values are due to the background interval starting
    # before follow-up, so we can fix these negative values to zero
    mutate(bg_day = pmax(bg_day, 0)) %>%
    # reshape back to wide
    dplyr::select(
        -bg_date,
        -bg_rank
    ) %>%
    tidyr::pivot_wider(
        names_from = "bg_interval",
        values_from = bg_day
    )


# Chop it like it's hot! -------------------------------------------------------
cat("Chop it like it's hot")

# first, chop up follow-up time according to the vaccination intervals
d_events_dose_long <-
    d_events %>%
    dplyr::select(
        study_id,
        matches("^dose2"),
        matches("^doseb")
    ) %>%
    #dplyr::select(-c(vacc_dose2_date, Dose2_name)) %>%
    dtplyr::lazy_dt() %>%
    tidyr::pivot_longer(
        cols = matches("dose"),
        names_to = "dose_interval",
        values_to = "tstart"
    ) %>%
    filter(assertr::not_na(tstart)) %>%
    as_tibble()

cat("!")

# second, chop up follow-up time according to the background intervals
d_events_bg_long <-
    d_events %>%
    dplyr::select(
        study_id,
        matches("bg_interval")
    ) %>%
    dtplyr::lazy_dt() %>%
    tidyr::pivot_longer(
        cols = matches("bg_interval"),
        names_to = "bg_interval",
        values_to = "tstart"
    )%>%
    filter(assertr::not_na(tstart)) %>%
    as_tibble()

cat("!")

# combine the two long datasets and add stop day and event flag
class(d_events_dose_long$tstart)
class(d_events_bg_long$tstart)

d_events_long <-
    d_events_dose_long %>%
    mutate(tstart_t = as.numeric(tstart)) %>%
    full_join(
        y  = d_events_bg_long,
        by = c("study_id", "tstart")
    ) %>%
    left_join(
        y  = d_events %>% dplyr::select(study_id, start_day, stop_day, event_flg),
        by = "study_id"
    ) %>%
    dplyr::select(
        study_id,
        tstart,
        event_flg,
        stop_day,
        dose_interval,
        bg_interval
    ) %>%
    arrange(study_id, tstart) %>%
    tidyr::fill(
        dose_interval,
        bg_interval
    ) 

## Can't get verify to work so do manually
## verify(assertr::not_na(stop_day))
# verify(assertr::not_na(event_flg)) 
# verify(assertr::not_na(tstart))
#summary(d_events_long)

cat("!")

# make tstop and update the event_flg
d_events_long <-
    d_events_long %>%
    dtplyr::lazy_dt() %>%
    group_by(study_id) %>%
    mutate(
        tstop = lead(tstart),
        event_flg = if_else(row_number() == n(), event_flg, 0),
        .after = tstart
    ) %>%
    ungroup() %>%
    mutate(
        tstop = if_else(is.na(tstop), stop_day, tstop)
    ) %>%
    as_tibble()

cat("!")

# simple sanity checks
#d_events_long <-
 #   d_events_long %>%
  #  verify(tstart < tstop) %>%
# z = d_events_long %>%
# mutate(flag = ifelse(tstart<tstop, 1,0)) %>%
# group_by(flag) %>%
# tally()
# verify(not_na(event_flg)) %>%
    #verify(not_na(dose_interval)) %>%
    #verify(not_na(bg_interval))
#summary(d_events_long)

cat("!\n")

d_events_long  %>% filter(study_id == d_events_long$study_id[29])

# Checks =======================================================================
cat("Checks\n")

# do some stronger checks to make sure we've not screwed anything up:
#   (1) we've not lost any people along the way
#   (2) number of events matches those in original d_events dataframe
#   (3) per person, starts with tstart 0 and ends with tstop == stop_day
#   (4) per person, tstart is 0 or lag(tstop)
#   (5) per person, dose_interval starts at dose2_day14

check_person_count <- n_distinct(d_events_long$study_id) == n_distinct(d_sample)

check_event_count <- sum(d_events_long$event_flg) == sum(d_events$event_flg)

check_per_person <-
    d_events_long %>%
    dtplyr::lazy_dt() %>%
    group_by(study_id) %>%
    summarise(
        has_tstart_eq_0 = first(tstart) == 0,
        has_tstop_eq_stopday = last(tstop) == last(stop_day),
        is_tstart_sequential = all(tstart == 0 | tstart == lag(tstop)),
        has_dose2_day014 = first(dose_interval == "dose2_day014")
    ) %>%
    ungroup() %>%
    as_tibble() %>%
    summarise(
        has_tstart_eq_0      = all(has_tstart_eq_0),
        has_tstop_eq_stopday = all(has_tstop_eq_stopday),
        is_tstart_sequential = all(is_tstart_sequential),
        has_dose2_day014     = all(has_dose2_day014)
    ) %>%
    mutate(
        all = all(
            has_tstart_eq_0,
            has_tstop_eq_stopday,
            is_tstart_sequential
            #has_dose2_day014
        )
    )

if (!check_person_count)
    stop("Number of people in d_events_long does not match d_sample")

if (!check_event_count)
    stop("Number of events in d_events_long does not match d_events")

if (!check_per_person$all)
    stop("Person checks did not pass -- go look at check_per_person")


# Add back the analysis covariates =============================================
cat("Add covariates\n")

expr_vacc <- "([A-Z]{2})_(dose.)_(day.+)"

d_analysis <-
    d_events_long %>%
    left_join(d_baseline, by = "study_id") %>%
    left_join(
        y  = d_background %>% dplyr::select(-bg_start),
        by = "bg_interval"
    ) %>%
    # prefix dose interval with name of vaccine
    mutate(
        is_dose2 = stringr::str_detect(dose_interval, "^dose2"),
        vacc_dose_interval = if_else(
            condition = is_dose2,
            true      = stringr::str_c(vacc_dose2_name, "_", dose_interval),
            false     = stringr::str_c(vacc_doseb_name, "_", dose_interval)
        ),
        # switch pattern of vacc dose interval to be:
        # dose number - vaccine name - time period
        vacc_dose_interval = vacc_dose_interval %>%
            stringr::str_replace(expr_vacc, "\\2_\\1_\\3") %>%
            factor()
    )


# Save =========================================================================
cat("Save\n")

qsave(
    d_sample,
    file = "results/d_sample_b.qs" #person level dataset for descriptives
)

qsave(
    d_analysis,
    file = "results/d_analysis_stats_b.qs" #dataset for poisson - has time series variables
)

#beep()

## Some checks
# d_sample %>% count(has_outcome)
# d_sample %>% count(vacc_dose2_name, has_outcome)

