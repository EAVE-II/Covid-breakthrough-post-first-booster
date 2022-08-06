library(tidyverse)
library(scales)
library(cowplot)

ARR_data <- read.csv("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/output/ARR_14-28_days_as_ref_b.csv")

expr_coef <- "([\\.0-9]+) \\(([\\.0-9]+), ([\\.0-9]+)\\)"

### dose interval
ARR_data1 <- ARR_data[c(2,3,6,7,9,10,12,13),]
z <- ARR_data1 %>% 
  mutate(
    xvar="vacc_dose_interval",
    xlbl=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
t_coef_overall <- z[,c(6:10)]

### age
ARR_data1 <- ARR_data[c(38:44),]
z <- ARR_data1 %>% 
  mutate(
    xvar="age_gp",
    xlbl=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
t_coef_overall <- rbind(t_coef_overall,z[,c(6:10)])

### BMI
ARR_data1 <- ARR_data[c(50:54),]
z <- ARR_data1 %>% 
  mutate(
    xvar="bmi_gp",
    xlbl=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
t_coef_overall <- rbind(t_coef_overall,z[,c(6:10)])

### Number of risk groups
ARR_data1 <- ARR_data[c(45:49),]
z <- ARR_data1 %>% 
  mutate(
    xvar="qcovid_cat",
    xlbl=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
t_coef_overall <- rbind(t_coef_overall,z[,c(6:10)])

### SIMD
ARR_data1 <- ARR_data[c(55:58),]
z <- ARR_data1 %>% 
  mutate(
    xvar="simd2019_quintile",
    xlbl=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
t_coef_overall <- rbind(t_coef_overall,z[,c(6:10)])

# tidy plot of demographics and health characteristics =======


t_coef_ref <- tribble(
  ~xvar,                ~xlbl, ~doseb_est, ~doseb_conf_low, ~doesb_conf_high,
  "vacc_dose_interval", "pv_period_fv2_2:9",  1.000, 1.000, 1.000,
  "age_gp",           "18-49",         1.000, 1.000, 1.000,
  "bmi_gp",            "18.5-24.9",     1.000, 1.000, 1.000,
  "qcovid_cat",         "0",             1.000, 1.000, 1.000,
  "simd2019_quintile",  "5-least",      1.000, 1.000, 1.000
)

lkp_xvar <- c(
  "Dose\ninterval"           = "vacc_dose_interval",
  "Age"                        = "age_gp",
  "BMI"                        = "bmi_gp",
  "QCovid\nscore"              = "qcovid_cat",
  "Area\ndeprivation\nquintile" = "simd2019_quintile"
)

lkp_xlbl <- c(
  # dose interval
  "Day 14-69" = "pv_period_fv2_2:9",
  "Day 70-139" = "pv_period_fv2_10:19",
  "Day 140+" = "pv_period_fv2_20+",
  "Day 14-34 Mo" = "pv_period_fv3_2:4_Mo",
  "Day 14-34 PB" = "pv_period_fv3_2:4_PB",
  "Day 35-55 Mo" = "pv_period_fv3_5:8_Mo",
  "Day 35-55 PB" = "pv_period_fv3_5:8_PB",
  "Day 56+ Mo" = "pv_period_fv3_9+_Mo",
  "Day 56+ PB"   = "pv_period_fv3_9+_PB",
  # age
  "18-49"  = "18-49",
  "50-54"  = "age_gp50-54",
  "55-59"  = "age_gp55-59",
  "60-64"  = "age_gp60-64",
  "65-69"  = "age_gp65-69",
  "70-74"  = "age_gp70-74",
  "75-79"  = "age_gp75-79",
  "80+"    = "age_gp80+",
  # bmi
  "<18.5"     = "bmi_gp<18.5",
  "18.5-24.9" = "18.5-24.9",
  "25.0-29.9" = "bmi_gp25-29.9",
  "30.0-34.9" = "bmi_gp30-34.9",
  "35.0-39.9" = "bmi_gp35-39.9",
  "40.0+"     = "bmi_gp40+",
  # qcovid
  "0"   = "0",
  "1"   = "n_risk_gps1",
  "2"   = "n_risk_gps2",
  "3"   = "n_risk_gps3",
  "4"   = "n_risk_gps4",
  "5+"  = "n_risk_gps5+",
  # simd
  "5 - Least deprived" = "5-least",
  "4"                  = "simd2020_sc_quintile4",
  "3"                  = "simd2020_sc_quintile3",
  "2"                  = "simd2020_sc_quintile2",
  "1"  = "simd2020_sc_quintile1 - High"
)

lkp_colour <- c(
  "reference" = "#000000",
  "estimate"  = "#0a6aa6"
)

p_coef_overall <-
  bind_rows(
    t_coef_overall %>% mutate(est_type = "estimate"),
    t_coef_ref     %>% mutate(est_type = "reference")
  ) 
p_coef_overall$doseb_conf_high<-coalesce(p_coef_overall$doseb_conf_high,
                                         p_coef_overall$doesb_conf_high)
p_coef_overall$doesb_conf_high<-NULL
p_coef_overall_o <-p_coef_overall%>% 
  filter(
    xvar %in% lkp_xvar,
    xlbl %in% lkp_xlbl
  ) %>% 
  mutate(
    xvar = factor(xvar, lkp_xvar, names(lkp_xvar)),
    xlbl = factor(xlbl, lkp_xlbl, names(lkp_xlbl))
  ) %>% 
  ggplot(aes(
    x = doseb_est,
    xmin = doseb_conf_low,
    xmax = doseb_conf_high,
    y = xlbl,
    colour = est_type
  )) +
  facet_grid(
    xvar ~ .,
    scales = "free",
    space = "free_y",
    switch = "y"
  ) +
  geom_vline(xintercept = 1) +
  geom_pointrange() +
  scale_x_continuous(
    name = "Adjusted risk ratio (95% CI)",
    breaks = pretty_breaks()
  ) +
  scale_colour_manual(values = lkp_colour) +
  theme(
    axis.title.y      = element_blank(),
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.x      = element_text(face = "bold"),
    strip.text.y.left = element_text(face = "bold"),
    legend.position   = "none",
    axis.text = element_text(face="bold")
  ) +
  coord_cartesian(
    xlim = c(0,7.5)
  ) +
  labs(
    title = "Demographic and health summaries"
  )

p_coef_overall_o

# plot QCovid items ==========================
ARR_data <- read.csv("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/output/ARR_conds_14-28_days_as_ref.csv")

expr_coef <- "([\\.0-9]+) \\(([\\.0-9]+), ([\\.0-9]+)\\)"

# conditions
z <- ARR_data %>% 
  mutate(
    term=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
t_coef_qcovid <- z[,c(6:9)]

lkp_qcovid <- c(
  "Blood or bone marrrow cancer" = "Q_DIAG_BLOOD_CANCER",
  "Rare neuron disease" = "Q_DIAG_NEURO",
  "Immunosuppressed" = "immuno",
  "Parkinson's disease" = "Q_DIAG_PARKINSONS",
  "Diabetes type 1" = "Q_DIAG_DIABETES_1",
  "Lung cancer" = "Q_DIAG_RESP_CANCER",
  "Cirrhosis of the liver" = "Q_DIAG_CIRRHOSIS",
  "Dementia" = "Q_DIAG_DEMENTIA",
  "Pulmonary hypertension" = "Q_DIAG_PULM_HYPER",
  "COPD" = "Q_DIAG_COPD",
  "Epilepsy" = "Q_DIAG_EPILEPSY",
  "Cerebral palsy" = "Q_DIAG_CEREBRALPALSY",
  "Kidney disease"        = "Q_DIAG_CKD",
  "Rheumatoid arthritis or SLE" = "Q_DIAG_RA_SLE",
  "Learning disability or Down's Syndrome" = "Q_LEARN_CAT",
  "Thrombosis or pulmonary embolus" = "Q_DIAG_VTE",
  "Care home or homeless" = "Q_HOME_CAT",
  "Stroke or TIA" = "Q_DIAG_STROKE",
  "Peripheral vascular disease" = "Q_DIAG_PVD",
  "Sickle cell disease" = "Q_DIAG_SICKLE_CELL",
  "Heart failure" = "Q_DIAG_CCF",
  "Diabetes type 2" = "Q_DIAG_DIABETES_2",
  "Atrial fibrillation" = "Q_DIAG_AF",
  "Prior fracture: hip, wrist, spine or humerus" = "Q_DIAG_FRACTURE",
  "Coronary heart disease" = "Q_DIAG_CHD",
  "Asthma" = "Q_DIAG_ASTHMA",
  "Severe mental illness" = "Q_DIAG_SEV_MENT_ILL",
  "Congenital heart disease" = "Q_DIAG_CONGEN_HD"
)

lkp_colour_qcovid <- "#fb8c61"

p_coef_qcovid <-
  t_coef_qcovid %>% 
  filter(term != "Q_DIAG_CEREBRALPALSY") %>% 
  filter(term %in% lkp_qcovid) %>% 
  mutate(
    term = factor(term, lkp_qcovid, names(lkp_qcovid)),
    term = fct_rev(term)
  ) %>% 
  ggplot(aes(
    x = doseb_est,
    xmin = doseb_conf_low,
    xmax = doseb_conf_high,
    y = term
  )) +
  geom_vline(xintercept = 1) +
  geom_pointrange(colour = lkp_colour_qcovid) +
  scale_x_continuous(
    name = "Adjusted risk ratio (95% CI)",
    breaks = pretty_breaks()
  ) +
  theme(
    axis.title.y      = element_blank(),
    strip.background  = element_blank(),
    axis.text = element_text(face="bold")
  ) +
  coord_cartesian(
    xlim = c(0.5, 3)
  ) +
  labs(
    title = "QCovid1 items"
  )

p_coef_qcovid

# combine =====

p_coef <- plot_grid(p_coef_overall_o, p_coef_qcovid)

ggsave(
  plot = p_coef,
  filename = "/conf/EAVE/GPanalysis/progs/UA/p_vacc_breakthrough_coef.png",
  width = 35,
  height = 20,
  units = "cm"
)

print(p_coef)
