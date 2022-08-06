source("r_clear_and_load.r")

# Load =========================================================================
cat("Load\n")

d_variant_testing <- qread(s_drive("d_variant_testing.qs"))


# Plot =========================================================================
cat("Plot\n")

lkp_fill <- c(
  "#444444",
  "#377eb8",
  "#ff7f00",
  "#984ea3",
  "#4daf4a",
  "#e41a1c"
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
