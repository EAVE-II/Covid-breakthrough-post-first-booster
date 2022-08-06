source("r_clear_and_load.r")

# Load =========================================================================
cat("Load\n")

d_hosp <- qread(s_drive("d_hosp_week_raw.qs"))

# Plot by admission method =====================================================
cat("Plot by admission method\n")

p_hosp_admis <-
	d_hosp %>%
	filter(hosp_start_week >= ymd("2020-01-01")) %>%
	mutate(spell_admis_method = replace_na(spell_admis_method, "Unknown")) %>%
	group_by(hosp_start_week, spell_admis_method) %>%
	summarise(n = sum(n)) %>%
	mutate(n = ifelse(1 <= n & n <= 9, 10, n)) %>%
	ungroup() %>%
	ggplot(aes(x = hosp_start_week, y = n)) +
	facet_wrap(~spell_admis_method, ncol = 1) +
	geom_col()
p_hosp_admis

# Plot emergency admission by covid status =====================================
cat("Plot emergency admission by covid status\n")

p_hosp_emg_c19 <-
	d_hosp %>%
	filter(hosp_start_week >= ymd("2020-01-01")) %>%
	filter(spell_admis_method == "Emergency admission") %>%
	group_by(hosp_start_week, epi_covid19_status) %>%
	summarise(n = sum(n)) %>%
	mutate(n = ifelse(1 <= n & n <= 9, 10, n)) %>%
	ggplot(aes(x = hosp_start_week, y = n)) +
	facet_wrap(~epi_covid19_status, ncol = 1) +
	geom_col()
p_hosp_emg_c19

# Save =========================================================================
cat("Save\n")

qsave(p_hosp_admis,   file = "results/p_hosp_admis.qs")
qsave(p_hosp_emg_c19, file = "results/p_hosp_emg_c19.qs")
