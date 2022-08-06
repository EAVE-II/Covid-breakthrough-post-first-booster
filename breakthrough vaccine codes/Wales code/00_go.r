stop("Don't actually source this script, that would be bad for business!")

# Info =========================================================================

# Tuesday  8th December  2021 - Vaccination programme started
# Monday  20th September 2021 - Booster programme started

# Q: What does vaccine uptake look like across age groups?
# Q: How effective is the booster vaccine against variant-specific infections?

# Date extraction ==============================================================

#source("01a_get_raw_cohort.r")
#source("01b_clean_cog_uk_variant_stats.r")
#source("01c_get_raw_hosp.r")

# Data cleaning ================================================================

source("02a_select_cohort.r")
source("02b_impute_bmi.r")

# Describe and explore cohort ==================================================

source("03a_plot_cog_uk_variant_stats.r")
source("03b_summarise_cohort_vacc_and_infection.r")
source("03c_plot_cohort_cuminc_vacc_uptake.r")
source("03d_plot_raw_hosp_counts.r")

# Analyse dose 2 and booster breakthrough ======================================

source("04a_make_background_summaries.r")
source("04b_create_analysis_sample.r")
source("04c_analyse_dose2_booster_breakthrough.r")
source("04d_analyse_booster_breakthrough.r")
source("04e_analyse_booster_breakthrough_conditions.r")

# Report =======================================================================

source("r_clear_and_load.r")

render(
    input = "99_report.rmd",
    output_format = html_document(toc = TRUE),
    quiet = TRUE
)
