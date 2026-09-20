# ============================================================
# 9 — OXYCODONE CONSUMPTION (mcg/kg/hr)
# Global reporting + PSM comparison vs morphine
# Benchmarked against Howard et al. 2010 (GOSH morphine NCA)
# ============================================================
#
# Consumption is computed as total_background_mg /
# episode_duration_dose_hrs from episode_dose.rds — the same
# delivered-dose quantity that 7_dosing.R, 8_true_total_dose.R Section 1,
# and the Table 4/6 background-only figures all use. This is not the
# same as oxycodone_mg from pca_episodes_flagged, which is the
# nominal drug content of the syringe/bag as prescribed (regex-
# extracted from MED_ORDER text, e.g. "Oxycodone 50mg in 50mL" -> 50)
# rather than the amount actually delivered to the patient over the
# episode — a syringe can be filled with 50mg and only have a
# fraction infused before being changed or the episode ending.
# Using the delivered-dose quantity also makes this figure genuinely
# comparable to Howard et al. 2010, since Howard's rate is derived
# from infusion pump records (delivered dose), not prescribed/nominal
# syringe content.
#
# A diagnostic comparison (nominal-strength method vs delivered-dose
# method, side by side) is printed near the top of this script's
# output, so the scale of the difference between the two approaches
# is visible before the PSM comparison further down.
#
# A >=2hr minimum-duration exclusion is applied consistently across
# every table below (global, surgical classification, age group, and
# the PSM comparison). This threshold is empirically justified: a
# sensitivity sweep (0-24hrs) on the full unfiltered cohort shows
# median consumption is essentially flat across every threshold
# (morphine exactly 4.0 mcg/kg/hr at every cutoff; oxycodone 3.93 at
# 0hrs vs 3.89 at 24hrs — about 1% drift, not a meaningful shift),
# and short episodes are not disproportionately noisy or extreme
# (max/IQR-width are, if anything, higher in the >=24hr bin than in
# any short bin). 2hrs matches the existing PSM exclusion criterion
# (exclude_short_pca) used elsewhere in the pipeline, and the sweep
# confirms this isn't costing anything by being too loose or too
# strict.
# ============================================================

MIN_DURATION_HRS <- 2   # see sensitivity note above

library(tidyverse)
library(lubridate)

cat("============================================================\n")
cat("ANALYSIS 12 — OXYCODONE CONSUMPTION (mcg/kg/hr)\n")
cat("============================================================\n\n")

# Require episode_shift_table, weight data, AND episode_dose
# (episode_dose.rds from 7_dosing.R -- NEW requirement for the fix)
stopifnot(exists("episode_shift_table"))
stopifnot(exists("pca_episodes_flagged"))
stopifnot(exists("s4_surgical_keys"))
stopifnot(exists("episode_dose"))

# ============================================================
# HELPER FUNCTIONS
# ============================================================

fmt_median_iqr <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  sprintf("%.2f (%.2f–%.2f)",
          round(median(x), digits),
          round(quantile(x, 0.25), digits),
          round(quantile(x, 0.75), digits))
}

# ============================================================
# DIAGNOSTIC — OLD (nominal syringe strength) vs NEW (delivered
# background dose), same episodes, so the scale of the fix is
# visible before anything downstream is trusted
#
# oxycodone_mg is not carried to episode level in
# pca_episodes_flagged (it's a per-syringe field, collapsed away
# during Section 3's episode aggregation), so it's rebuilt here as
# an episode-level sum of oxycodone_mg across all syringes in the
# episode, using dose_data (built in 7_dosing.R, which already links
# each syringe to its correct episode_id).
# ============================================================

cat("--- DIAGNOSTIC: OLD (nominal mg) vs NEW (delivered mg) method ---\n\n")

stopifnot(exists("dose_data"))  # from 7_dosing.R

nominal_mg_by_episode <- dose_data |>
  filter(drug == "oxycodone") |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  summarise(
    oxycodone_mg_nominal_total = sum(oxycodone_mg, na.rm = TRUE),
    .groups = "drop"
  )

old_vs_new <- pca_episodes_flagged |>
  filter(drug == "oxycodone") |>
  select(PAT_ENC_CSN_ID, episode_id, WEIGHT_KG, episode_duration_hrs) |>
  inner_join(nominal_mg_by_episode, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  inner_join(
    episode_dose |>
      filter(drug == "oxycodone") |>
      select(PAT_ENC_CSN_ID, episode_id, total_background_mg, episode_duration_dose_hrs),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  ) |>
  mutate(
    old_mcg_per_kg_hr = (oxycodone_mg_nominal_total * 1000) / (WEIGHT_KG * episode_duration_hrs),
    new_mcg_per_kg_hr = (total_background_mg * 1000) / (WEIGHT_KG * episode_duration_dose_hrs),
    ratio_old_to_new  = old_mcg_per_kg_hr / new_mcg_per_kg_hr
  )

cat(sprintf("Episodes with both fields available: %d\n\n", nrow(old_vs_new)))

old_vs_new |>
  summarise(
    n = n(),
    old_median = round(median(old_mcg_per_kg_hr, na.rm = TRUE), 2),
    new_median = round(median(new_mcg_per_kg_hr, na.rm = TRUE), 2),
    median_ratio_old_to_new = round(median(ratio_old_to_new, na.rm = TRUE), 2)
  ) |>
  print()

cat("\nIf median_ratio_old_to_new is well above 1 (e.g. several-fold),\n")
cat("this confirms the old figure reflects unused syringe content,\n")
cat("not delivered dose. Proceeding with the corrected (NEW) method\n")
cat("for all figures below.\n\n")
# ============================================================
# CALCULATE CONSUMPTION (delivered background dose)
# ============================================================

cat("Calculating oxycodone consumption (mcg/kg/hr)...\n\n")

daily_consumption_prefloor <- episode_dose |>
  filter(drug == "oxycodone") |>
  select(PAT_ENC_CSN_ID, episode_id, WEIGHT_KG,
         total_background_mg, episode_duration_dose_hrs) |>
  mutate(
    mcg_per_kg_hr = (total_background_mg * 1000) / (WEIGHT_KG * episode_duration_dose_hrs)
  )

daily_consumption <- daily_consumption_prefloor |>
  filter(episode_duration_dose_hrs >= MIN_DURATION_HRS)

cat(sprintf("Duration floor applied (>=%.0fhrs): %d -> %d episodes (%d dropped)\n\n",
            MIN_DURATION_HRS, nrow(daily_consumption_prefloor), nrow(daily_consumption),
            nrow(daily_consumption_prefloor) - nrow(daily_consumption)))

# ============================================================
# SUMMARY — EPISODE LEVEL
# ============================================================

cat("OXYCODONE CONSUMPTION — EPISODE LEVEL\n")
cat("(mcg/kg/hr, delivered background dose — comparable to Howard et al. 2010)\n\n")

episode_summary <- daily_consumption |>
  summarise(
    n_episodes = n(),
    median_mcg_per_kg_hr = round(median(mcg_per_kg_hr, na.rm=TRUE), 2),
    q25 = round(quantile(mcg_per_kg_hr, 0.25, na.rm=TRUE), 2),
    q75 = round(quantile(mcg_per_kg_hr, 0.75, na.rm=TRUE), 2),
    .groups = "drop"
  )

print(episode_summary)

cat("\nHoward et al. 2010 benchmarks (morphine, D0 = day of NCA start):\n")
cat("  Overall D0->D1: 18.48 -> 11.75 mcg/kg/hr (mean, 36%% daily fall)\n")
cat("  (NOTE: earlier versions of this cat() incorrectly cited '8.4 D0 /\n")
cat("  27.23 D1' -- those are actually neonate-vs-16yr+ age comparison\n")
cat("  figures, BOTH measured on D0, not a D0-vs-D1 comparison. Fixed.)\n\n")

# ============================================================
# CONSUMPTION BY SURGICAL CLASSIFICATION
# ============================================================

cat("OXYCODONE CONSUMPTION BY SURGICAL CLASSIFICATION\n\n")

consumption_by_class <- daily_consumption |>
  mutate(
    is_surgical = (paste(PAT_ENC_CSN_ID, episode_id, sep="_") %in%
                     paste(s4_surgical_keys$PAT_ENC_CSN_ID, s4_surgical_keys$episode_id, sep="_")),
    classification = if_else(is_surgical, "Surgical", "Non-surgical")
  ) |>
  group_by(classification) |>
  summarise(
    n_episodes = n(),
    median_mcg = round(median(mcg_per_kg_hr, na.rm=TRUE), 2),
    q25 = round(quantile(mcg_per_kg_hr, 0.25, na.rm=TRUE), 2),
    q75 = round(quantile(mcg_per_kg_hr, 0.75, na.rm=TRUE), 2),
    .groups = "drop"
  )

print(consumption_by_class)
cat("\n")

# ============================================================
# CONSUMPTION BY AGE GROUP
# ============================================================

cat("OXYCODONE CONSUMPTION BY AGE GROUP\n\n")

consumption_by_age <- daily_consumption |>
  left_join(
    pca_episodes_flagged |>
      filter(drug == "oxycodone") |>
      select(PAT_ENC_CSN_ID, episode_id, age_group),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  ) |>
  group_by(age_group) |>
  summarise(
    n_episodes = n(),
    median_mcg = round(median(mcg_per_kg_hr, na.rm=TRUE), 2),
    q25 = round(quantile(mcg_per_kg_hr, 0.25, na.rm=TRUE), 2),
    q75 = round(quantile(mcg_per_kg_hr, 0.75, na.rm=TRUE), 2),
    .groups = "drop"
  ) |>
  arrange(age_group)

print(consumption_by_age)
cat("\n")

# ============================================================
# PSM COMPARISON: OXYCODONE vs MORPHINE
# ============================================================

cat("============================================================\n")
cat("PSM COMPARISON: OXYCODONE vs MORPHINE CONSUMPTION\n")
cat("============================================================\n\n")

# Reads psm_matched_keys.rds directly from disk (rather than
# depending on matched_outcomes already existing in memory), so this
# script is self-contained regardless of run order, as long as
# 10_comparative_PSM.R has been run at some point.
matched_outcomes <- tryCatch(
  readRDS("1 - data/3 - processed_data/psm_matched_keys.rds") |>
    select(PAT_ENC_CSN_ID, episode_id, drug),
  error = function(e) {
    stop("psm_matched_keys.rds not found — run 10_comparative_PSM.R first (at least once, on disk).")
  }
)

cat("Matched keys loaded:", nrow(matched_outcomes), "episodes\n\n")

# episode_dose already carries drug + WEIGHT_KG + total_background_mg
# + episode_duration_dose_hrs for BOTH drugs, so no separate
# oxy/mor join-with-suffix dance is needed here (that was a source
# of fragility in the old version).
morphine_oxy_psm_prefloor <- matched_outcomes |>
  select(PAT_ENC_CSN_ID, episode_id, drug) |>
  left_join(
    episode_dose |>
      select(PAT_ENC_CSN_ID, episode_id, drug,
             WEIGHT_KG, total_background_mg, episode_duration_dose_hrs),
    by = c("PAT_ENC_CSN_ID", "episode_id", "drug")
  ) |>
  mutate(
    mcg_per_kg_hr_actual = (total_background_mg * 1000) /
      (WEIGHT_KG * episode_duration_dose_hrs)
  )

# NOTE: this floor should be redundant here in practice, since
# matched_outcomes (from 10_comparative_PSM.R) is already built from the
# PSM population, which has the exclude_short_pca flag (>=2hrs)
# applied upstream. Applying it explicitly again here is a
# belt-and-braces check, not expected to change n -- if
# it DOES drop episodes, that's worth investigating (would mean
# the upstream exclusion isn't being carried through as assumed).
morphine_oxy_psm <- morphine_oxy_psm_prefloor |>
  filter(is.na(episode_duration_dose_hrs) | episode_duration_dose_hrs >= MIN_DURATION_HRS)

cat(sprintf("Duration floor check (>=%.0fhrs): %d -> %d episodes (%d dropped)\n",
            MIN_DURATION_HRS, nrow(morphine_oxy_psm_prefloor), nrow(morphine_oxy_psm),
            nrow(morphine_oxy_psm_prefloor) - nrow(morphine_oxy_psm)))
cat("(expected: 0 dropped, since PSM population already excludes <2hr episodes upstream)\n\n")

cat(sprintf("Matched episodes with dose data: %d / %d\n\n",
            sum(!is.na(morphine_oxy_psm$mcg_per_kg_hr_actual)),
            nrow(morphine_oxy_psm)))

psm_comparison <- morphine_oxy_psm |>
  group_by(drug) |>
  summarise(
    n_episodes = n(),
    median_mcg = round(median(mcg_per_kg_hr_actual, na.rm=TRUE), 2),
    q25 = round(quantile(mcg_per_kg_hr_actual, 0.25, na.rm=TRUE), 2),
    q75 = round(quantile(mcg_per_kg_hr_actual, 0.75, na.rm=TRUE), 2),
    .groups = "drop"
  )

print(psm_comparison)

cat("\nFor reference — PREVIOUSLY LOCKED (nominal-strength method,\n")
cat("now known incorrect) PSM figures:\n")
cat("  Morphine:  20.8 (14.6-41.5) mcg/kg/hr\n")
cat("  Oxycodone: 16.1 (12.6-25.4) mcg/kg/hr\n\n")

cat("Check: does the DIRECTION and RELATIVE MAGNITUDE of the\n")
cat("morphine-vs-oxycodone gap survive the correction? The\n")
cat("manuscript's Discussion argues oxycodone shows lower\n")
cat("consumption despite 1:1 dosing (interpreted as a potency/\n")
cat("efficacy signal, consistent with Lenz 2009 / Silvasti 1998 /\n")
cat("Kalso 1991 and POPCORN's 1:1 approach). That interpretation\n")
cat("depends on the RATIO between drugs, not the absolute numbers\n")
cat("-- if the ratio changes materially, the framing needs\n")
cat("re-examining, not just the reported figures.\n\n")

cat("============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n")
cat("Oxycodone consumption now calculated as delivered background\n")
cat("dose (rate x syringe duration, summed per episode) / weight /\n")
cat("episode duration -- SAME methodology as episode_dose.rds\n")
cat("background-only figures (7_dosing.R) and 8_true_total_dose.R\n")
cat("Section 1's background_mcg_kg_hr. Comparable to Howard et al.\n")
cat("2010 (delivered-dose-based, not nominal-syringe-strength-based).\n\n")
cat("NOTE: this is BACKGROUND dose only, consistent with existing\n")
cat("Table 1/4 methodology. It does NOT include bolus doses --\n")
cat("see 8_true_total_dose.R (episode_true_total_dose.rds) for a bolus-\n")
cat("inclusive TOTAL delivered dose figure, available for the\n")
cat("flowsheet-linked subset only (~73-79%% of episodes).\n")
cat("============================================================\n\n")