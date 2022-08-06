source("r_clear_and_load.r")

# Get raw cohort ===============================================================
cat("Get raw cohort\n")

con <- sail_open()

q_hosp_raw <- "
SELECT
	DATE(DATE_TRUNC('WEEK', epi_start_date)) AS hosp_start_week,
	spell_admis_method,
	CASE
		WHEN admis_covid19_cause_flg = 1 THEN 'covid19_cause'
		WHEN admis_with_covid19_flg = 1 THEN 'with_covid19'
		WHEN covid19_during_admis_flg = 1 THEN 'covid19_during'
		WHEN within_14days_positive_pcr_flg = 1 THEN 'positive_pcr'
	END AS epi_covid19_status,
	COUNT(*) AS n
FROM
	sailw1151v.dacvap_hosp_long
GROUP BY
	DATE_TRUNC('WEEK', epi_start_date),
	spell_admis_method,
	CASE
		WHEN admis_covid19_cause_flg = 1 THEN 'covid19_cause'
		WHEN admis_with_covid19_flg = 1 THEN 'with_covid19'
		WHEN covid19_during_admis_flg = 1 THEN 'covid19_during'
		WHEN within_14days_positive_pcr_flg = 1 THEN 'positive_pcr'
	END
;"

d_hosp_raw <- sail_run(con, q_hosp_raw)


# Save =========================================================================
cat("Save\n")

qsave(
    d_hosp_raw,
    file = s_drive("d_hosp_week_raw.qs")
)


# Goodbye ======================================================================

sail_close(con)
beep()
