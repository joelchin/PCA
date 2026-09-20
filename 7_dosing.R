# ============================================================
# 7 — Opioid Dose and Consumption
# Paediatric PCA Study
# Addenbrooke's Hospital
#
# Purpose:
# Characterise opioid dosing across the cohort
# Primary: background infusion rate and bolus dose by drug and age
# Secondary: total estimated background dose per episode (weighted)
# Comparative: oxycodone vs morphine dose in matched population
#
# Table 4 (the direct oxycodone-vs-morphine comparison in the matched
# population) reports background_mg_kg_hr_true as its primary
# background-rate metric — a duration-weighted per-hour rate
# (total_background_mg_kg divided by episode_duration_dose_hrs),
# consistent with the nominal-vs-delivered dose correction applied
# elsewhere in the pipeline (see 9_comsumption.R and
# 4_flags.R/10_comparative_PSM.R). background_rate_mg_kg_hr_median (a
# simple median of per-syringe rates, which doesn't account for how
# long each syringe was actually running) is retained alongside it as
# a legacy/reference figure, since it's the metric other descriptive
# tables in this script (Tables 1, 2, 3, 5, 6) still report and
# comparing the two side by side is informative.
#
# Data source: pca_mar_clean
# - background_rate_ml: extracted from PROGRAMME text (audit.R s2)
# - bolus_dose_ml: extracted from PROGRAMME text
# - oxycodone_conc_mg_per_ml: extracted from MED_ORDER
# - morphine_conc_mg_per_ml: extracted from MED_ORDER
# - TAKEN_TIME / Discontinue_Time: syringe start/end
#
# Dose derivation:
# background_rate_mg_hr = background_rate_ml * concentration
# background_rate_mg_kg_hr = background_rate_mg_hr / WEIGHT_KG
# bolus_dose_mg = bolus_dose_ml * concentration
# bolus_dose_mg_kg = bolus_dose_mg / WEIGHT_KG
# total_background_mg = sum(rate_mg_hr * syringe_duration_hrs)
# per episode across all syringes (weighted by syringe duration)
#
# Note: 181 oxycodone records missing Discontinue_Time (last syringe)
# → imputed with episode_end from pca_episodes_flagged
# ============================================================

library(tidyverse)
library(lubridate)

stopifnot(
  exists("pca_mar_clean"),
  exists("pca_episodes_flagged"),
  exists("episode_metrics")
)

cat("============================================================\n")
cat("ANALYSIS 10 — OPIOID DOSE AND CONSUMPTION\n")
cat("============================================================\n\n")

fmt_med_iqr <- function(x, digits = 2) {
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0) return("NA")
  sprintf("%.*f (%.*f–%.*f)",
          digits, median(x),
          digits, quantile(x, 0.25),
          digits, quantile(x, 0.75))
}

fmt_n_pct <- function(n, total) {
  sprintf("%d (%.1f%%)", n, 100 * n / total)
}

# ============================================================
# STEP 1 — BUILD SYRINGE-LEVEL DOSE DATASET
# ============================================================

cat("--- Building syringe-level dose dataset ---\n\n")

# ── CONCENTRATION FIX ──────────────────────────────────────────────────────
# 2,046 morphine records use MED_ORDER format:
# "morphine 50 mg in 50 mL paediatric PCA/NCA (> 50 kg)"
# This is the >50kg adult formulation — concentration = 1 mg/ml
# The audit.R regex captures "morphine sulfate X mg in Y mL" but not
# "morphine X mg in Y mL" (without sulfate) — fix here and flag for audit.R

cat("Fixing missing morphine concentration (>50kg formulation)...\n")

pca_mar_clean <- pca_mar_clean |>
  mutate(
    morphine_conc_mg_per_ml = case_when(
      !is.na(morphine_conc_mg_per_ml) ~ morphine_conc_mg_per_ml,
      drug == "morphine" &
        MED_ORDER == "morphine 50 mg in 50 mL paediatric PCA/NCA (> 50 kg)" ~ 1.0,
      TRUE ~ morphine_conc_mg_per_ml
    )
  )

n_missing_after <- sum(is.na(pca_mar_clean$morphine_conc_mg_per_ml[
  pca_mar_clean$drug == "morphine"]))
cat(sprintf("Missing morphine concentration after fix: %d\n\n", n_missing_after))

# TODO: Fix audit.R regex to capture this format in next pipeline rerun
# Add to morphine concentration extraction:
# str_extract(MED_ORDER, "morphine (?:sulfate )?([0-9.]+) mg in ([0-9.]+) mL")

# Bridge: assign each syringe record to its episode via time overlap
# A syringe belongs to an episode if TAKEN_TIME falls within
# episode_start → episode_end AND drug matches

episode_ends <- pca_episodes_flagged |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end,
         WEIGHT_KG, AGE, age_group, drug)

# Use findInterval approach — faster than many-to-many join
# Split by patient, match syringe times to episode intervals

cat("Linking syringe records to episodes via time overlap...\n")

syringe_episode_link <- pca_mar_clean |>
  filter(drug %in% c("oxycodone", "morphine")) |>
  select(PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME, drug) |>
  left_join(
    episode_ends,
    by = c("PAT_ENC_CSN_ID", "drug"),
    relationship = "many-to-many"
  ) |>
  filter(
    TAKEN_TIME >= episode_start,
    is.na(episode_end) | TAKEN_TIME <= episode_end
  ) |>
  # If a syringe matches multiple episodes take the closest start
  group_by(PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME, drug) |>
  slice_min(abs(as.numeric(TAKEN_TIME - episode_start)), n = 1) |>
  ungroup() |>
  select(PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME, drug, episode_id,
         episode_end, WEIGHT_KG, AGE, age_group)

cat(sprintf("Syringe records linked: %d / %d\n",
            nrow(syringe_episode_link),
            sum(pca_mar_clean$drug %in% c("oxycodone", "morphine"))))
cat("\n")

dose_data <- pca_mar_clean |>
  filter(drug %in% c("oxycodone", "morphine")) |>
  left_join(syringe_episode_link,
            by = c("PAT_ENC_CSN_ID", "ORDER_MED_ID", "TAKEN_TIME", "drug"),
            relationship = "many-to-many") |>
  filter(!is.na(episode_id)) |>
  mutate(
    # Concentration per drug
    concentration_mg_per_ml = case_when(
      drug == "oxycodone" ~ oxycodone_conc_mg_per_ml,
      drug == "morphine" ~ morphine_conc_mg_per_ml,
      TRUE ~ NA_real_
    ),
    # Impute missing Discontinue_Time with episode_end
    syringe_end = if_else(
      is.na(Discontinue_Time),
      episode_end,
      Discontinue_Time
    ),
    # Syringe duration in hours
    syringe_duration_hrs = as.numeric(
      difftime(syringe_end, TAKEN_TIME, units = "hours")
    ),
    # Cap implausible syringe durations (>7 days = data error)
    syringe_duration_hrs = if_else(
      syringe_duration_hrs > 168 | syringe_duration_hrs < 0,
      NA_real_,
      syringe_duration_hrs
    ),
    # Background rate in mg/hr
    background_rate_mg_hr = background_rate_ml * concentration_mg_per_ml,
    # Background rate in mg/kg/hr
    background_rate_mg_kg_hr = background_rate_mg_hr / WEIGHT_KG,
    # Bolus dose in mg
    bolus_dose_mg = bolus_dose_ml * concentration_mg_per_ml,
    # Bolus dose in mg/kg
    bolus_dose_mg_kg = bolus_dose_mg / WEIGHT_KG,
    # Loading dose in mg/kg
    loading_dose_mg = loading_dose_ml * concentration_mg_per_ml,
    loading_dose_mg_kg = loading_dose_mg / WEIGHT_KG,
    # Background dose delivered this syringe (mg)
    syringe_background_mg = background_rate_mg_hr * syringe_duration_hrs,
    # Background dose per kg this syringe
    syringe_background_mg_kg = syringe_background_mg / WEIGHT_KG
  )

cat("Syringe records built:\n")
cat(sprintf(" Oxycodone: %d\n", sum(dose_data$drug == "oxycodone")))
cat(sprintf(" Morphine: %d\n", sum(dose_data$drug == "morphine")))
cat(sprintf(" Missing concentration: %d\n",
            sum(is.na(dose_data$concentration_mg_per_ml))))
cat(sprintf(" Missing syringe duration after imputation: %d\n",
            sum(is.na(dose_data$syringe_duration_hrs))))
cat("\n")

# ============================================================
# STEP 2 — EPISODE-LEVEL DOSE SUMMARY
# Aggregate syringe records to episode level
# ============================================================

cat("--- Aggregating to episode level ---\n\n")

episode_dose <- dose_data |>
  group_by(PAT_ENC_CSN_ID, episode_id, drug, AGE, age_group, WEIGHT_KG) |>
  summarise(
    n_syringes = n(),
    # Background rate — median across syringes (typical prescribed rate).
    # Retained for reference — see header note. Not the primary
    # comparison metric, since it doesn't account for how long each
    # syringe actually ran.
    background_rate_mg_kg_hr_median = median(background_rate_mg_kg_hr,
                                             na.rm = TRUE),
    # Bolus dose — first syringe (prescribed at initiation)
    bolus_dose_mg_kg_first = first(bolus_dose_mg_kg),
    # Loading dose — first syringe only
    had_loading_dose = any(!is.na(loading_dose_mg_kg)),
    loading_dose_mg_kg = first(loading_dose_mg_kg[!is.na(loading_dose_mg_kg)]),
    # Total background dose delivered (weighted sum across syringes)
    total_background_mg = sum(syringe_background_mg, na.rm = TRUE),
    total_background_mg_kg = sum(syringe_background_mg_kg, na.rm = TRUE),
    # Episode duration from dose data (sum of syringe durations)
    episode_duration_dose_hrs = sum(syringe_duration_hrs, na.rm = TRUE),
    # Daily dose (mg/kg/day) — normalised to 24hrs
    daily_background_mg_kg = if_else(
      episode_duration_dose_hrs > 0,
      total_background_mg_kg / episode_duration_dose_hrs * 24,
      NA_real_
    ),
    # Concentration used
    concentration_mg_per_ml = first(concentration_mg_per_ml),
    .groups = "drop"
  ) |>
  mutate(
    # ============================================================
    # background_mg_kg_hr_true
    #
    # The duration-weighted per-hour delivered-dose rate —
    # total_background_mg_kg divided by episode_duration_dose_hrs
    # (mathematically identical to daily_background_mg_kg / 24, just
    # computed directly at the per-hour scale this pipeline reports
    # consumption in elsewhere). This is the metric that should be
    # used for any oxycodone-vs-morphine RATE comparison — it
    # properly weights each syringe's contribution by how long it
    # actually ran, unlike background_rate_mg_kg_hr_median above.
    # ============================================================
    background_mg_kg_hr_true = if_else(
      episode_duration_dose_hrs > 0,
      total_background_mg_kg / episode_duration_dose_hrs,
      NA_real_
    ),
    # Flag implausible total doses (>10x expected max)
    flag_implausible_dose = daily_background_mg_kg > 2 |
      daily_background_mg_kg < 0
  )

cat(sprintf("Episode dose summaries: %d\n", nrow(episode_dose)))
cat(sprintf(" Implausible doses flagged: %d\n",
            sum(episode_dose$flag_implausible_dose, na.rm = TRUE)))
cat("\n")

# Remove implausible doses for analysis
episode_dose_clean <- episode_dose |>
  filter(!flag_implausible_dose | is.na(flag_implausible_dose))

# ============================================================
# TABLE 1 — PRESCRIBED DOSING REGIMEN BY DRUG
# Background rate, bolus dose, loading dose
# ============================================================

cat("============================================================\n")
cat("TABLE 1 — PRESCRIBED DOSING REGIMEN BY DRUG\n")
cat("============================================================\n\n")

dose_by_drug <- episode_dose_clean |>
  group_by(drug) |>
  summarise(
    n = n(),
    # Background rate
    background_mg_kg_hr_median = round(median(background_rate_mg_kg_hr_median,
                                              na.rm = TRUE), 4),
    background_mg_kg_hr_q25 = round(quantile(background_rate_mg_kg_hr_median,
                                             0.25, na.rm = TRUE), 4),
    background_mg_kg_hr_q75 = round(quantile(background_rate_mg_kg_hr_median,
                                             0.75, na.rm = TRUE), 4),
    # Bolus dose
    bolus_mg_kg_median = round(median(bolus_dose_mg_kg_first,
                                      na.rm = TRUE), 4),
    bolus_mg_kg_q25 = round(quantile(bolus_dose_mg_kg_first,
                                     0.25, na.rm = TRUE), 4),
    bolus_mg_kg_q75 = round(quantile(bolus_dose_mg_kg_first,
                                     0.75, na.rm = TRUE), 4),
    # Loading dose
    pct_had_loading_dose = round(100 * mean(had_loading_dose,
                                            na.rm = TRUE), 1),
    loading_mg_kg_median = round(median(loading_dose_mg_kg,
                                        na.rm = TRUE), 4),
    # Total and daily dose
    total_bg_mg_kg_median = round(median(total_background_mg_kg,
                                         na.rm = TRUE), 2),
    daily_bg_mg_kg_median = round(median(daily_background_mg_kg,
                                         na.rm = TRUE), 4),
    .groups = "drop"
  )

cat("Background infusion rate (mg/kg/hr) — median (IQR) [legacy,\n")
cat("median-of-per-syringe-rates — see header note]:\n")
for (i in seq_len(nrow(dose_by_drug))) {
  cat(sprintf(" %-12s %s\n", dose_by_drug$drug[i],
              fmt_med_iqr(episode_dose_clean$background_rate_mg_kg_hr_median[
                episode_dose_clean$drug == dose_by_drug$drug[i]])))
}
cat("\nBolus dose (mg/kg) — median (IQR):\n")
for (i in seq_len(nrow(dose_by_drug))) {
  cat(sprintf(" %-12s %s\n", dose_by_drug$drug[i],
              fmt_med_iqr(episode_dose_clean$bolus_dose_mg_kg_first[
                episode_dose_clean$drug == dose_by_drug$drug[i]])))
}
cat("\nLoading dose prescribed — n (%):\n")
for (i in seq_len(nrow(dose_by_drug))) {
  n_load <- sum(episode_dose_clean$had_loading_dose[
    episode_dose_clean$drug == dose_by_drug$drug[i]], na.rm = TRUE)
  n_tot <- sum(episode_dose_clean$drug == dose_by_drug$drug[i])
  cat(sprintf(" %-12s %s\n", dose_by_drug$drug[i], fmt_n_pct(n_load, n_tot)))
}
cat("\nEstimated daily background dose (mg/kg/day) — median (IQR):\n")
for (i in seq_len(nrow(dose_by_drug))) {
  cat(sprintf(" %-12s %s\n", dose_by_drug$drug[i],
              fmt_med_iqr(episode_dose_clean$daily_background_mg_kg[
                episode_dose_clean$drug == dose_by_drug$drug[i]])))
}
cat("\n")

# ============================================================
# TABLE 2 — DOSING BY AGE GROUP (oxycodone)
# Key for dosing guidance — what doses are actually used?
# ============================================================

cat("============================================================\n")
cat("TABLE 2 — OXYCODONE DOSING BY AGE GROUP\n")
cat("============================================================\n\n")

oxy_dose_by_age <- episode_dose_clean |>
  filter(drug == "oxycodone") |>
  group_by(age_group) |>
  summarise(
    n = n(),
    background_mg_kg_hr = fmt_med_iqr(background_rate_mg_kg_hr_median),
    bolus_mg_kg = fmt_med_iqr(bolus_dose_mg_kg_first),
    pct_loading = round(100 * mean(had_loading_dose, na.rm = TRUE), 1),
    daily_bg_mg_kg = fmt_med_iqr(daily_background_mg_kg),
    conc_mg_per_ml_median = round(median(concentration_mg_per_ml, na.rm = TRUE), 3),
    .groups = "drop"
  ) |>
  arrange(age_group)

cat("Oxycodone dosing by age group:\n\n")
print(oxy_dose_by_age)
cat("\n")

# Concentration used by age group
cat("--- Oxycodone concentration (mg/ml) by age group ---\n")
episode_dose_clean |>
  filter(drug == "oxycodone") |>
  group_by(age_group) |>
  summarise(
    n = n(),
    conc_median = round(median(concentration_mg_per_ml, na.rm = TRUE), 3),
    conc_min = round(min(concentration_mg_per_ml, na.rm = TRUE), 3),
    conc_max = round(max(concentration_mg_per_ml, na.rm = TRUE), 3),
    .groups = "drop"
  ) |>
  arrange(age_group) |>
  print()
cat("\n")

# ============================================================
# TABLE 3 — DOSING BY AGE GROUP (morphine)
# ============================================================

cat("============================================================\n")
cat("TABLE 3 — MORPHINE DOSING BY AGE GROUP\n")
cat("============================================================\n\n")

mor_dose_by_age <- episode_dose_clean |>
  filter(drug == "morphine") |>
  group_by(age_group) |>
  summarise(
    n = n(),
    background_mg_kg_hr = fmt_med_iqr(background_rate_mg_kg_hr_median),
    bolus_mg_kg = fmt_med_iqr(bolus_dose_mg_kg_first),
    pct_loading = round(100 * mean(had_loading_dose, na.rm = TRUE), 1),
    daily_bg_mg_kg = fmt_med_iqr(daily_background_mg_kg),
    .groups = "drop"
  ) |>
  arrange(age_group)

print(mor_dose_by_age)
cat("\n")

# ============================================================
# TABLE 4 — DOSE COMPARISON IN MATCHED POPULATION
# Does oxycodone require more or less opioid than morphine?
#
# background_mg_kg_hr_true (duration-weighted rate) is reported and
# tested as the primary background-rate comparison here, computed on
# the same matched population (psm_matched_keys.rds) as the
# manuscript's actual Table 4 (10_comparative_PSM.R /
# 16_publication_tables.R). Note this is not the same figure as that
# table's background-rate reference row, which is sourced from
# episode_true_total_dose.rds (8_true_total_dose.R) and restricted to
# the flowsheet-coverage subset — this table's background_mg_kg_hr_true
# uses the full MAR-linked population instead. Both are legitimate,
# differently-scoped background-rate metrics; they should not be
# expected to match exactly. See header note for full rationale.
# background_rate_mg_kg_hr_median is retained alongside, explicitly
# labelled legacy/reference, not primary.
# ============================================================

cat("============================================================\n")
cat("TABLE 4 — DOSE IN MATCHED POPULATION (oxycodone vs morphine)\n")
cat("============================================================\n\n")

# Load matched keys
matched_keys <- tryCatch(
  readRDS("1 - data/3 - processed_data/psm_matched_keys.rds"),
  error = function(e) {
    cat("psm_matched_keys.rds not found — run 10_comparative_PSM.R first\n")
    return(NULL)
  }
)

if (!is.null(matched_keys)) {
  matched_dose <- episode_dose_clean |>
    semi_join(matched_keys, by = c("PAT_ENC_CSN_ID", "episode_id"))
  
  cat("Matched episodes with dose data:",
      nrow(matched_dose), "\n\n")
  
  cat("--- PRIMARY: duration-weighted background rate (mg/kg/hr) ---\n")
  matched_dose |>
    group_by(drug) |>
    summarise(
      n = n(),
      background_mg_kg_hr_true = fmt_med_iqr(background_mg_kg_hr_true),
      .groups = "drop"
    ) |>
    print()
  cat("\n")
  
  cat("--- Other dose metrics (background rate shown as legacy\n")
  cat("reference only — see header note) ---\n")
  matched_dose |>
    group_by(drug) |>
    summarise(
      n = n(),
      background_mg_kg_hr_median_LEGACY = fmt_med_iqr(background_rate_mg_kg_hr_median),
      bolus_mg_kg = fmt_med_iqr(bolus_dose_mg_kg_first),
      daily_bg_mg_kg = fmt_med_iqr(daily_background_mg_kg),
      total_bg_mg_kg = fmt_med_iqr(total_background_mg_kg),
      pct_loading = round(100 * mean(had_loading_dose,
                                     na.rm = TRUE), 1),
      .groups = "drop"
    ) |>
    print()
  cat("\n")
  
  # Statistical comparison of doses in matched population
  cat("--- Wilcoxon tests (matched population) ---\n")
  mor_matched <- matched_dose |> filter(drug == "morphine")
  oxy_matched <- matched_dose |> filter(drug == "oxycodone")
  
  cat("PRIMARY:\n")
  wt_primary <- wilcox.test(oxy_matched$background_mg_kg_hr_true,
                            mor_matched$background_mg_kg_hr_true, na.rm = TRUE)
  cat(sprintf(" %-40s p=%.4f\n", "background_mg_kg_hr_true", wt_primary$p.value))
  
  cat("\nOther metrics (legacy/secondary):\n")
  for (var in c("background_rate_mg_kg_hr_median",
                "bolus_dose_mg_kg_first",
                "daily_background_mg_kg",
                "total_background_mg_kg")) {
    wt <- wilcox.test(oxy_matched[[var]], mor_matched[[var]], na.rm = TRUE)
    cat(sprintf(" %-40s p=%.4f\n", var, wt$p.value))
  }
  cat("\n")
  
  cat("If PRIMARY's p-value here diverges materially from the\n")
  cat("manuscript's actual locked Table 4 consumption figures\n")
  cat("(10_comparative_PSM.R / 16_publication_tables.R), investigate\n")
  cat("before treating either as final. Note the two aren't expected to\n")
  cat("match exactly regardless: this table's background rate is\n")
  cat("MAR-linked (episode_dose.rds, full matched population), while\n")
  cat("Table 4's background-rate reference row is flowsheet-validated\n")
  cat("(episode_true_total_dose.rds, restricted to >=80% flowsheet\n")
  cat("coverage) — a materially different, smaller population.\n\n")
}

# ============================================================
# TABLE 5 — DOSE BY ERA
# Did dosing change as experience grew?
# ============================================================

cat("============================================================\n")
cat("TABLE 5 — OXYCODONE DOSING BY ERA\n")
cat("============================================================\n\n")

era_data <- pca_episodes_flagged |>
  select(PAT_ENC_CSN_ID, episode_id, era)

episode_dose_clean |>
  filter(drug == "oxycodone") |>
  left_join(era_data, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  group_by(era) |>
  summarise(
    n = n(),
    background_mg_kg_hr = fmt_med_iqr(background_rate_mg_kg_hr_median),
    bolus_mg_kg = fmt_med_iqr(bolus_dose_mg_kg_first),
    daily_bg_mg_kg = fmt_med_iqr(daily_background_mg_kg),
    .groups = "drop"
  ) |>
  print()
cat("\n")

# ============================================================
# TABLE 6 — BOLUS DOSE RANGE CHECK
# Are doses within expected safe ranges?
# For manuscript methods: what concentration and dose ranges used
# ============================================================

cat("============================================================\n")
cat("TABLE 6 — DOSE RANGE SUMMARY (safety check)\n")
cat("============================================================\n\n")

cat("--- Oxycodone ---\n")
oxy_range <- episode_dose_clean |>
  filter(drug == "oxycodone") |>
  summarise(
    n = n(),
    # Background rate
    bg_mg_kg_hr_min = round(min(background_rate_mg_kg_hr_median, na.rm = TRUE), 4),
    bg_mg_kg_hr_max = round(max(background_rate_mg_kg_hr_median, na.rm = TRUE), 4),
    bg_mg_kg_hr_median = round(median(background_rate_mg_kg_hr_median, na.rm = TRUE), 4),
    # Bolus
    bolus_mg_kg_min = round(min(bolus_dose_mg_kg_first, na.rm = TRUE), 4),
    bolus_mg_kg_max = round(max(bolus_dose_mg_kg_first, na.rm = TRUE), 4),
    bolus_mg_kg_median = round(median(bolus_dose_mg_kg_first, na.rm = TRUE), 4),
    # Concentration
    conc_min = round(min(concentration_mg_per_ml, na.rm = TRUE), 3),
    conc_max = round(max(concentration_mg_per_ml, na.rm = TRUE), 3)
  )
print(t(oxy_range))
cat("\n")

cat("--- Morphine ---\n")
mor_range <- episode_dose_clean |>
  filter(drug == "morphine") |>
  summarise(
    n = n(),
    bg_mg_kg_hr_min = round(min(background_rate_mg_kg_hr_median, na.rm = TRUE), 4),
    bg_mg_kg_hr_max = round(max(background_rate_mg_kg_hr_median, na.rm = TRUE), 4),
    bg_mg_kg_hr_median = round(median(background_rate_mg_kg_hr_median, na.rm = TRUE), 4),
    bolus_mg_kg_min = round(min(bolus_dose_mg_kg_first, na.rm = TRUE), 4),
    bolus_mg_kg_max = round(max(bolus_dose_mg_kg_first, na.rm = TRUE), 4),
    bolus_mg_kg_median = round(median(bolus_dose_mg_kg_first, na.rm = TRUE), 4),
    conc_min = round(min(concentration_mg_per_ml, na.rm = TRUE), 3),
    conc_max = round(max(concentration_mg_per_ml, na.rm = TRUE), 3)
  )
print(t(mor_range))
cat("\n")

# ============================================================
# SAVE EPISODE DOSE DATASET
# ============================================================

saveRDS(episode_dose_clean,
        "1 - data/3 - processed_data/episode_dose.rds")

# This script builds/saves the object under the name
# episode_dose_clean, but every downstream script (16_publication_tables.R,
# 9_comsumption.R, 8_true_total_dose.R) checks for and uses an in-session object
# literally called episode_dose. Aliased here so the session object
# name matches what everything downstream expects.
episode_dose <- episode_dose_clean

cat("episode_dose.rds saved\n\n")
cat("episode_dose (alias of episode_dose_clean) available in session\n\n")

cat("============================================================\n")
cat("ANALYSIS 10 COMPLETE\n")
cat("Key outputs:\n")
cat(" - Background infusion rates by drug and age group\n")
cat(" - Bolus doses by drug and age group\n")
cat(" - Total estimated background dose per episode\n")
cat(" - Matched comparison: oxycodone vs morphine dose (now using\n")
cat("   background_mg_kg_hr_true as PRIMARY, see header note)\n")
cat(" - Dosing range for manuscript methods section\n")
cat(" - Dose trends by era\n")
cat("\n")
cat("LIMITATIONS:\n")
cat(" - Background dose estimated from prescribed rate × syringe duration\n")
cat(" - Does not capture bolus activations (flowsheet data needed)\n")
cat(" - 181 missing syringe end times imputed with episode_end\n")
cat(" - Daily dose calculation assumes continuous infusion at prescribed rate\n")
cat(" - Tables 1, 2, 3, 5, 6 still report background_rate_mg_kg_hr_median\n")
cat("   (the legacy, non-duration-weighted metric) as their primary rate\n")
cat("   figure — only this script's Table 4 (the direct matched-\n")
cat("   population comparison) uses the duration-weighted metric,\n")
cat("   since it carries the highest risk of producing a conflicting\n")
cat("   result if left on the legacy metric. Worth a follow-up pass to\n")
cat("   update the remaining descriptive tables for full consistency.\n")
cat("============================================================\n")