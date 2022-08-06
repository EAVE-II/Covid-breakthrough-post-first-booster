source("r_clear_and_load.r")

x_seed <- 125460070


# Load =========================================================================
cat("Load\n")

d_cohort <-
    qread(s_drive("d_cohort_clean.qs")) %>%
    rename(move_out_date = c20_end_date)

d_background <- qread(s_drive("d_background.qs"))


# Select analysis sample =======================================================
cat("Select analysis sample\n")

date_13sep21 <- ymd("2021-09-13")
date_20dec21 <- ymd("2021-12-20")

d_cohort <-
    d_cohort %>%
    mutate(
        epoch = vacc_dose2_date + ddays(14),
        # outcome: covid-19 hosp or death
        outcome_date = pmin(hosp_admis_date, death_covid_date, na.rm = TRUE),
        # has outcome before any other event which would stop follow-up
        has_outcome = epoch <= outcome_date &
                    (is.na(death_noncovid_date) | outcome_date <= death_noncovid_date) &
                    (is.na(move_out_date)       | outcome_date <= move_out_date) &
                    outcome_date <= study_end_date,
        has_outcome = if_else(is.na(outcome_date), FALSE, has_outcome, has_outcome)
    ) %>%
    mutate(
        # sample selection criteria
        has_dose2 = not_na(vacc_dose1_date) &
                    not_na(vacc_dose2_date) &
                    (is.na(death_noncovid_date) | epoch <= death_noncovid_date) &
                    (is.na(move_out_date)       | epoch <= move_out_date) &
                    epoch <= study_end_date &
                    vacc_dose2_name %in% c("AZ", "PB"),
        no_prv_outcome = is.na(outcome_date) | (epoch <= outcome_date),
        valid_doseb = is.na(vacc_doseb_date) | vacc_doseb_name %in% c("MD", "PB"),
        # person is still in the study as of 13th September 2021
        in_study_13sep21 =
            (date_13sep21 < outcome_date        | is.na(outcome_date)) &
            (date_13sep21 < death_noncovid_date | is.na(death_noncovid_date)) &
            (date_13sep21 < move_out_date       | is.na(move_out_date)) &
            (date_13sep21 < study_end_date      | is.na(study_end_date)),
        # person:
        #   - has had booster
        #   - no events prior to booster
        #   - is still in the study as of 20th Decemeber 2021
        in_study_20dec21 =
            # had booster
            not_na(vacc_doseb_date) & not_na(vacc_doseb_name) &
            # still in study from 20 Decemember 2021
            (date_20dec21 < outcome_date        | is.na(outcome_date)) &
            (date_20dec21 < death_noncovid_date | is.na(death_noncovid_date)) &
            (date_20dec21 < move_out_date       | is.na(move_out_date)) &
            (date_20dec21 < study_end_date      | is.na(study_end_date)) &
            # no events prior to booster
            (vacc_doseb_date < outcome_date        | is.na(outcome_date)) &
            (vacc_doseb_date < death_noncovid_date | is.na(death_noncovid_date)) &
            (vacc_doseb_date < move_out_date       | is.na(move_out_date)) &
            (vacc_doseb_date < study_end_date      | is.na(study_end_date))
    )

t_analysis_sample_selection <- tribble(
    ~step, ~description, ~total_n, ~outcome_n,
        0, "All cohort with at least two doses",
            d_cohort %>% nrow(),
            d_cohort %>% filter(has_outcome) %>% nrow(),
        1, "Second dose is AZ or PB",
            d_cohort %>% filter(has_dose2) %>% nrow(),
            d_cohort %>% filter(has_dose2) %>% filter(has_outcome) %>% nrow(),
        2, "No outcome prior to second dose",
            d_cohort %>% filter(has_dose2, no_prv_outcome) %>% nrow(),
            d_cohort %>% filter(has_dose2, no_prv_outcome) %>% filter(has_outcome) %>% nrow(),
        3, "If had a booster, booster is MD or PB",
            d_cohort %>% filter(has_dose2, no_prv_outcome, valid_doseb) %>% nrow(),
            d_cohort %>% filter(has_dose2, no_prv_outcome, valid_doseb) %>% filter(has_outcome) %>% nrow(),
        4, "In study as of 13th September 2021",
            d_cohort %>% filter(has_dose2, no_prv_outcome, valid_doseb, in_study_13sep21) %>% nrow(),
            d_cohort %>% filter(has_dose2, no_prv_outcome, valid_doseb, in_study_13sep21) %>% filter(has_outcome) %>% nrow(),
        5, "In study as of 20th December 2021, and has had booster",
            d_cohort %>% filter(has_dose2, no_prv_outcome, valid_doseb, in_study_20dec21) %>% nrow(),
            d_cohort %>% filter(has_dose2, no_prv_outcome, valid_doseb, in_study_20dec21) %>% filter(has_outcome) %>% nrow()
)

t_analysis_sample_selection <-
    t_analysis_sample_selection %>%
    mutate(total_n_diff = total_n - lag(total_n), .after = total_n) %>%
    mutate(outcome_n_diff = outcome_n - lag(outcome_n), .after = outcome_n)

print(t_analysis_sample_selection)


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

# split out the columns into:
#   * baseline covariates: age, sex, qcovid, health board, ...
#   * events: outcome date, move out date, vaccination dates, ...
# we will then transform d_events into a start-stop dataset and
# eventually left join back on the baseline measures

d_baseline <-
    d_sample %>%
    select(
    # person id
        alf_e,
    # vaccine names
        vacc_dose2_name,
        vacc_doseb_name,
    # characteristics
        sex,
        age_3cat,
        age_4cat,
        age_5y_cat,
        ethn_cat,
        wimd2019_quintile,
        urban_rural_class,
        qcovid_cat,
        health_board,
        test_pre_dose2_cat,
        dose2_prior_infection_cat,
        vacc_dose1_dose2_diff_cat,
        bmi_cat,
        hypt_cat,
    # qcovid
        qc_b2_82,
        qc_b2_leukolaba,
        qc_b2_prednisolone,
        qc_b_af,
        qc_b_ccf,
        qc_b_asthma,
        qc_b_bloodcancer,
        qc_b_cerebralpalsy,
        qc_b_chd,
        qc_b_cirrhosis,
        qc_b_congenheart,
        qc_b_copd,
        qc_b_dementia,
        qc_b_epilepsy,
        qc_b_fracture4,
        qc_b_neurorare,
        qc_b_parkinsons,
        qc_b_pulmhyper,
        qc_b_pulmrare,
        qc_b_pvd,
        qc_b_ra_sle,
        qc_b_respcancer,
        qc_b_semi,
        qc_b_sicklecelldisease,
        qc_b_stroke,
        qc_diabetes1,
        qc_diabetes2,
        qc_b_vte,
        qc_chemo_cat,
        qc_carehome,
        qc_homeless,
        qc_learndisab,
        qc_downs,
        qc_p_marrow6,
        qc_p_radio6,
        qc_p_solidtransplant,
        qc_ckd3,
        qc_ckd4,
        qc_ckd5,
        renal_transplant,
    # study design
        cc_group,
        sample_weight,
        in_study_13sep21,
        in_study_20dec21
    )

d_events <-
    d_sample %>%
    select(
    # person id
        alf_e,
    # event dates
        outcome_date,
        death_noncovid_date,
        move_out_date,
        study_end_date,
    # vaccination
        vacc_dose2_date,
        vacc_dose2_name,
        vacc_doseb_date,
        vacc_doseb_name,
    # study design
        cc_group
    )


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
        dose2_day014 = vacc_dose2_date + ddays( 14), # dose 2: week  2
        dose2_day084 = vacc_dose2_date + ddays( 91), # dose 2: week 13
        dose2_day168 = vacc_dose2_date + ddays(182), # dose 2: week 26
        # dates for post booster cut points
        doseb_day014 = vacc_doseb_date + ddays( 14), # dose B: week  2
        doseb_day028 = vacc_doseb_date + ddays( 35), # dose B: week  5
        doseb_day056 = vacc_doseb_date + ddays( 56)  # dose B: week  8
    ) %>%
    # dose2_day*_date columns should only have a value if the person had
    # not had reached day 14 of their booster dose at that time point
    mutate(across(
        .cols = matches("dose2_day[0-9]+"),
        .fns  = ~ if_else(.x >= doseb_day014, NA_Date_, .x, .x)
    )) %>%
    # stop follow-up if they have the outcome, move out, die, or we reach
    # the end of the study, if multiple events on the same day, prioritise
    # the outcome
    mutate(
        stop_date = pmin(
            outcome_date,
            death_noncovid_date,
            move_out_date,
            study_end_date,
            na.rm = TRUE
        ),
        stop_day = stop_date,
        event_cat = case_when(
            stop_date == outcome_date        ~ "covid19_hosp_death",
            stop_date == death_noncovid_date ~ "death_noncovid",
            stop_date == move_out_date       ~ "move_out",
            stop_date == study_end_date      ~ "study_end"
        ),
        event_cat = factor(event_cat) %>% fct_relevel("study_end"),
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
    )) %>%
    # sanity checks
    verify(
        start_day < stop_day
    ) %>%
    verify(
        (event_flg == 1 & cc_group == "case") |
        (event_flg == 0 & cc_group == "control")
    )


# Add in background dates ======================================================
cat("Add background dates\n")

d_background <-
    d_background %>%
    select(
        bg_interval,
        bg_start,
        bg_avg_infection,
        bg_variant
    )

# make a wide dataset of the intervals
d_background_wide <-
    d_background %>%
    select(bg_interval, bg_start) %>%
    pivot_wider(
        names_from  = bg_interval,
        values_from = bg_start
    )

# make sure background intervals only have values if they overlap with
# follow-up
d_events <-
    bind_cols(d_events, d_background_wide) %>%
    # pivot first so it is easier to manipulate one column of intervals,
    # rather than X columns for X intervals
    pivot_longer(
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
    lazy_dt() %>%
    group_by(alf_e) %>%
    mutate(bg_rank = max(bg_rank, na.rm = TRUE)) %>%
    ungroup() %>%
    as_tibble() %>%
    mutate(bg_day = if_else(bg_day < bg_rank, NA_real_, bg_day, bg_day)) %>%
    # negative bg_day values are due to the background interval starting
    # before follow-up, so we can fix these negative values to zero
    mutate(bg_day = pmax(bg_day, 0)) %>%
    # reshape back to wide
    select(
        -bg_date,
        -bg_rank
    ) %>%
    pivot_wider(
        names_from = "bg_interval",
        values_from = bg_day
    )


# Chop it like it's hot! -------------------------------------------------------
cat("Chop it like it's hot")

# first, chop up follow-up time according to the vaccination intervals
d_events_dose_long <-
    d_events %>%
    select(
        alf_e,
        matches("^dose2"),
        matches("^doseb")
    ) %>%
    lazy_dt() %>%
    pivot_longer(
        cols = matches("dose"),
        names_to = "dose_interval",
        values_to = "tstart"
    ) %>%
    filter(not_na(tstart)) %>%
    as_tibble()

cat("!")

# second, chop up follow-up time according to the background intervals
d_events_bg_long <-
    d_events %>%
    select(
        alf_e,
        matches("bg_interval")
    ) %>%
    lazy_dt() %>%
    pivot_longer(
        cols = matches("bg_interval"),
        names_to = "bg_interval",
        values_to = "tstart"
    ) %>%
    filter(not_na(tstart)) %>%
    as_tibble()

cat("!")

# combine the two long datasets and add stop day and event flag
d_events_long <-
    d_events_dose_long %>%
    full_join(
        y  = d_events_bg_long,
        by = c("alf_e", "tstart")
    ) %>%
    left_join(
        y  = d_events %>% select(alf_e, start_day, stop_day, event_flg),
        by = "alf_e"
    ) %>%
    select(
        alf_e,
        tstart,
        event_flg,
        stop_day,
        dose_interval,
        bg_interval
    ) %>%
    arrange(alf_e, tstart) %>%
    fill(
        dose_interval,
        bg_interval
    ) %>%
    verify(not_na(stop_day)) %>%
    verify(not_na(event_flg)) %>%
    verify(not_na(tstart))

cat("!")

# make tstop and update the event_flg
d_events_long <-
    d_events_long %>%
    lazy_dt() %>%
    group_by(alf_e) %>%
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
d_events_long <-
    d_events_long %>%
    verify(tstart < tstop) %>%
    verify(not_na(event_flg)) %>%
    verify(not_na(dose_interval)) %>%
    verify(not_na(bg_interval))

cat("!\n")


# Checks =======================================================================
cat("Checks\n")

# do some stronger checks to make sure we've not screwed anything up:
#   (1) we've not lost any people along the way
#   (2) number of events matches those in original d_events dataframe
#   (3) per person, starts with tstart 0 and ends with tstop == stop_day
#   (4) per person, tstart is 0 or lag(tstop)
#   (5) per person, dose_interval starts at dose2_day14

check_person_count <- n_distinct(d_events_long$alf_e) == n_distinct(d_sample)

check_event_count <- sum(d_events_long$event_flg) == sum(d_events$event_flg)

check_per_person <-
    d_events_long %>%
    lazy_dt() %>%
    group_by(alf_e) %>%
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
            is_tstart_sequential,
            has_dose2_day014
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
    left_join(d_baseline, by = "alf_e") %>%
    left_join(
        y  = d_background %>% select(-bg_start),
        by = "bg_interval"
    ) %>%
    # prefix dose interval with name of vaccine
    mutate(
        is_dose2 = str_detect(dose_interval, "^dose2"),
        vacc_dose_interval = if_else(
            condition = is_dose2,
            true      = str_c(vacc_dose2_name, "_", dose_interval),
            false     = str_c(vacc_doseb_name, "_", dose_interval)
        ),
        # switch pattern of vacc dose interval to be:
        # dose number - vaccine name - time period
        vacc_dose_interval = vacc_dose_interval %>%
            str_replace(expr_vacc, "\\2_\\1_\\3") %>%
            factor()
    )


# Save =========================================================================
cat("Save\n")

qsave(
    t_analysis_sample_selection,
    file = "results/t_analysis_sample_selection.qs"
)

qsave(
    t_cc_group_n,
    file = "results/t_analysis_sample_cc_group.qs"
)

qsave(
    d_sample,
    file = s_drive("d_sample.qs")
)

qsave(
    d_analysis,
    file = s_drive("d_analysis.qs")
)

beep()
