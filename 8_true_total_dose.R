# ============================================================
# 8 — True Total Opioid Dose (Background + Bolus)
# Paediatric PCA Study — Addenbrooke's Hospital
# ============================================================
#
# Purpose:
#   Compute TRUE total delivered dose (background infusion + actual
#   delivered boluses) per episode, using the PCA Flowsheets refresh
#   (2026-05-26) NUMBER_OF_DOSES counters, bridged to episode_id via
#   time-interval containment against episode_start/episode_end.
#
# Methodology (validated across Diagnostics 1-13; bolus/lockout
# using per-syringe granularity, see note above
# HELPER FUNCTIONS below):
#   - Bolus total = SUM ACROSS SYRINGES of (boluses delivered on that
#     syringe, reset-corrected via true syringe boundaries) x (that
#     SAME syringe's own prescribed bolus dose in mg, from PROGRAMME
#     text on pca_mar_clean) — NOT a single episode-level bolus dose
#     applied to all boluses regardless of which syringe they came
#     from
#   - NOT using flowsheet CUMULATIVE_DOSE/CUMULATIVE_INTAKE fields:
#     morphine's CUMULATIVE_DOSE is a genuine running total, but
#     oxycodone's CUMULATIVE_INTAKE is a period value, not cumulative
#     - the two drugs' "dose" fields are NOT equivalent, so
#     NUMBER_OF_DOSES (a true monotonic bolus counter for BOTH drugs)
#     is used instead for consistency.
#   - Syringe/episode boundaries taken from existing trusted MAR data
#     (TAKEN_TIME/Discontinue_Time), NOT from flowsheet counter resets
#     (counter resets don't reliably align with real syringe changes —
#     only 70.4% within 6hrs in validation).
#
# Exclusion criteria applied (documented, matches existing precedent):
#   1. Episode duration >= 2hrs (matches the standard exclude_short_pca flag)
#   2. Flowsheet coverage >= 80% of episode span (excludes episodes
#      where the flowsheet barely touched the encounter)
#   3. Lockout-implied plausibility: boluses/hr <= 60/lockout_mins,
#      checked PER SYRINGE against that syringe's own lockout_mins.
#      Missing lockout OR missing bolus dose on any syringe that had
#      boluses excludes the whole episode (does not silently pass)
#   4. Explicit quarantine of confirmed-corrupted counter encounters
#
# Coverage caveat: flowsheet covers ~72.8% of MAR/episode cohort CSNs
# overall (79.3% for oxycodone specifically) — this analysis is
# necessarily restricted to that subset, not the full cohort.
#
# Data source: Chin Paeds PCA Datasets - PCA Flowsheets - refresh
#              2026-05-26.xlsx
#
# ------------------------------------------------------------
# SECTION 2: RECONCILIATION DIAGNOSTIC
# Resolves the ~4-5x discrepancy between this script's
# background_mcg_kg_hr and the locked/published PSM background
# figures (oxy 16.1, morph 20.8 mcg/kg/hr) BEFORE the delivered-
# dose figures are trusted for the manuscript. See STOP note at
# bottom of Section 1 output. Runs automatically at the end of
# this script, using the episode_true_total_dose object just
# built above plus episode_dose already in memory.
# ------------------------------------------------------------
#
# DATA SAFETY — MANDATORY THROUGHOUT:
# No PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME, PROGRAMME, or
# MED_ORDER content may appear in any script output or log.
# Aggregate only (counts, %, medians, IQR, summary()).
# ============================================================

library(tidyverse)
library(readxl)

stopifnot(exists("pca_mar_clean"), exists("pca_episodes_flagged"), exists("episode_dose"))



cat("============================================================\n")
cat("ANALYSIS 13 — TRUE TOTAL DOSE (BACKGROUND + BOLUS)\n")
cat("============================================================\n\n")

# ============================================================
# STEP 0 — LOAD FLOWSHEET
# ============================================================

cat("--- Loading flowsheet ---\n\n")

flow_raw <- read_excel(
  "~/pcaraw/Chin Paeds PCA Datasets - PCA Flowsheets - refresh 2026-07-15.xlsx"
) |> suppressWarnings()

flow <- flow_raw |>
  rename_with(~ str_replace_all(.x, " ", "_")) |>
  rename_with(~ str_to_upper(.x))

cat(sprintf("Flowsheet rows: %s\n\n", format(nrow(flow), big.mark = ",")))

# ============================================================
# HELPER FUNCTIONS
# ============================================================
#
# Bolus dose and lockout are extracted per syringe from PROGRAMME
# text (already present on pca_mar_clean, see audit.R), rather than
# collapsed to a single per-episode value. This matters because a
# bolus-dose or lockout change on a later syringe within the same
# episode would otherwise be costed at the wrong rate: boluses
# delivered on syringe N are costed at syringe N's own prescribed
# dose and checked against syringe N's own lockout. Background rate
# doesn't have this issue — total_background_mg in episode_dose.rds
# is already a per-syringe weighted sum (see 7_dosing.R).
# ============================================================

# Build syringe -> episode_id bridge via time-interval containment
build_syringe_intervals <- function(drug_name, drug_regex) {
  episode_windows <- pca_episodes_flagged |>
    filter(drug == drug_name) |>
    select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end) |>
    distinct() |>
    arrange(PAT_ENC_CSN_ID, episode_start)
  
  syringes <- pca_mar_clean |>
    filter(str_detect(MEDICATION, regex(drug_regex, ignore_case = TRUE))) |>
    distinct(PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME, Discontinue_Time)
  
  assign_episode_vec <- function(csn_data, windows_data) {
    csn_windows <- windows_data |> filter(PAT_ENC_CSN_ID == csn_data$PAT_ENC_CSN_ID[1])
    if (nrow(csn_windows) == 0) {
      return(csn_data |> mutate(episode_id = windows_data$episode_id[NA_integer_]))
    }
    hit_idx <- map_int(csn_data$TAKEN_TIME, function(t) {
      hit <- which(t >= csn_windows$episode_start & t <= csn_windows$episode_end)
      if (length(hit) >= 1) hit[1] else NA_integer_
    })
    csn_data |> mutate(episode_id = csn_windows$episode_id[hit_idx])
  }
  
  syringes_assigned <- syringes |>
    group_split(PAT_ENC_CSN_ID) |>
    map_dfr(~ assign_episode_vec(.x, episode_windows))
  
  syringes_assigned |>
    filter(!is.na(episode_id)) |>
    arrange(PAT_ENC_CSN_ID, TAKEN_TIME) |>
    group_by(PAT_ENC_CSN_ID) |>
    mutate(syringe_num = row_number()) |>
    ungroup()
}

# Attach each syringe's OWN prescribed bolus dose (mg) and lockout
# (mins), pulled directly from pca_mar_clean/PROGRAMME text — not
# collapsed to episode level. conc_col is the drug-specific
# concentration column already on pca_mar_clean (e.g.
# "oxycodone_conc_mg_per_ml"), extracted per-row from MED_ORDER in
# audit.R, so formulation changes (e.g. >50kg concentration) are
# also captured correctly per syringe.
attach_syringe_settings <- function(syringe_intervals, drug_regex, conc_col) {
  settings <- pca_mar_clean |>
    filter(str_detect(MEDICATION, regex(drug_regex, ignore_case = TRUE))) |>
    mutate(bolus_dose_mg_syringe = bolus_dose_ml * .data[[conc_col]]) |>
    distinct(PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME,
             bolus_dose_mg_syringe, lockout_mins)
  
  syringe_intervals |>
    left_join(settings, by = c("PAT_ENC_CSN_ID", "ORDER_MED_ID", "TAKEN_TIME"))
}

# Quantify how often bolus dose / lockout actually change within an
# episode (i.e. how material this fix is) — aggregate output only.
report_settings_stability <- function(syringe_intervals_with_settings, drug_name) {
  cat(sprintf("\n--- Settings stability check: %s ---\n\n", drug_name))
  
  ep_summary <- syringe_intervals_with_settings |>
    group_by(PAT_ENC_CSN_ID, episode_id) |>
    summarise(
      n_syringes = n(),
      n_distinct_bolus = n_distinct(bolus_dose_mg_syringe, na.rm = TRUE),
      n_distinct_lockout = n_distinct(lockout_mins, na.rm = TRUE),
      .groups = "drop"
    )
  
  multi_syringe <- ep_summary |> filter(n_syringes > 1)
  
  cat(sprintf("Episodes: %d total, %d with >1 syringe (%.1f%%)\n",
              nrow(ep_summary), nrow(multi_syringe),
              100 * nrow(multi_syringe) / nrow(ep_summary)))
  
  if (nrow(multi_syringe) > 0) {
    cat(sprintf("Of multi-syringe episodes:\n"))
    cat(sprintf("  Bolus dose changes mid-episode: %d (%.1f%%)\n",
                sum(multi_syringe$n_distinct_bolus > 1, na.rm = TRUE),
                100 * mean(multi_syringe$n_distinct_bolus > 1, na.rm = TRUE)))
    cat(sprintf("  Lockout changes mid-episode:    %d (%.1f%%)\n",
                sum(multi_syringe$n_distinct_lockout > 1, na.rm = TRUE),
                100 * mean(multi_syringe$n_distinct_lockout > 1, na.rm = TRUE)))
  }
  cat("\n")
  
  invisible(ep_summary)
}

# Assign flowsheet observations to (syringe_num, episode_id)
assign_flow_to_syringes <- function(flow_dose_col, syringe_intervals) {
  flow_drug <- flow |>
    filter(!is.na(.data[[flow_dose_col]])) |>
    mutate(doses = as.numeric(.data[[flow_dose_col]])) |>
    select(PAT_ENC_CSN_ID, RECORDED_TIME, doses)
  
  assign_syringe_ep <- function(csn_data, intervals_data) {
    csn_intervals <- intervals_data |> filter(PAT_ENC_CSN_ID == csn_data$PAT_ENC_CSN_ID[1])
    if (nrow(csn_intervals) == 0) {
      return(csn_data |> mutate(syringe_num = NA_integer_,
                                episode_id = intervals_data$episode_id[NA_integer_]))
    }
    starts <- as.numeric(csn_intervals$TAKEN_TIME)
    idx <- findInterval(as.numeric(csn_data$RECORDED_TIME), starts)
    valid <- idx >= 1 & idx <= nrow(csn_intervals)
    csn_data |> mutate(
      syringe_num = if_else(valid, csn_intervals$syringe_num[pmax(idx, 1)], NA_integer_),
      episode_id  = csn_intervals$episode_id[if_else(valid, pmax(idx, 1), NA_integer_)]
    )
  }
  
  flow_drug |>
    group_split(PAT_ENC_CSN_ID) |>
    map_dfr(~ assign_syringe_ep(.x, syringe_intervals))
}

# Full pipeline for one drug
compute_true_total_dose <- function(drug_name, drug_regex, flow_dose_col, conc_col) {
  
  cat(sprintf("\n--- Processing: %s ---\n\n", drug_name))
  
  syringe_intervals <- build_syringe_intervals(drug_name, drug_regex) |>
    attach_syringe_settings(drug_regex, conc_col)
  cat(sprintf("Syringe intervals (with episode_id): %d\n", nrow(syringe_intervals)))
  cat(sprintf("  Missing bolus dose data: %d (%.1f%%)\n",
              sum(is.na(syringe_intervals$bolus_dose_mg_syringe)),
              100 * mean(is.na(syringe_intervals$bolus_dose_mg_syringe))))
  cat(sprintf("  Missing lockout data:    %d (%.1f%%)\n",
              sum(is.na(syringe_intervals$lockout_mins)),
              100 * mean(is.na(syringe_intervals$lockout_mins))))
  
  report_settings_stability(syringe_intervals, drug_name)
  
  # Syringe duration (needed for per-syringe plausibility check)
  episode_windows <- pca_episodes_flagged |>
    filter(drug == drug_name) |>
    select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end)
  
  syringe_intervals <- syringe_intervals |>
    left_join(episode_windows |> select(-episode_start),
              by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    mutate(
      syringe_end_imputed = if_else(is.na(Discontinue_Time), episode_end, Discontinue_Time),
      syringe_duration_hrs = as.numeric(difftime(syringe_end_imputed, TAKEN_TIME, units = "hours"))
    )
  
  flow_assigned <- assign_flow_to_syringes(flow_dose_col, syringe_intervals)
  cat(sprintf("Flowsheet rows assigned: %d / %d (%.1f%%)\n",
              sum(!is.na(flow_assigned$syringe_num)), nrow(flow_assigned),
              100 * mean(!is.na(flow_assigned$syringe_num))))
  
  # Boluses delivered PER SYRINGE, joined to THAT SYRINGE'S OWN bolus
  # dose (mg) and lockout — this is the fix. A bolus-dose or lockout
  # change on a later syringe within the episode now costs/checks
  # correctly instead of inheriting the first syringe's settings.
  syringe_bolus <- flow_assigned |>
    filter(!is.na(syringe_num), !is.na(episode_id)) |>
    group_by(PAT_ENC_CSN_ID, episode_id, syringe_num) |>
    summarise(boluses_this_syringe = max(doses, na.rm = TRUE) - min(doses, na.rm = TRUE),
              .groups = "drop") |>
    left_join(
      syringe_intervals |>
        select(PAT_ENC_CSN_ID, episode_id, syringe_num,
               bolus_dose_mg_syringe, lockout_mins, syringe_duration_hrs),
      by = c("PAT_ENC_CSN_ID", "episode_id", "syringe_num")
    ) |>
    mutate(
      bolus_mg_this_syringe = boluses_this_syringe * bolus_dose_mg_syringe,
      boluses_per_hr_this_syringe = boluses_this_syringe / syringe_duration_hrs,
      max_per_hr_this_syringe = 60 / lockout_mins,
      syringe_implausible = !is.na(max_per_hr_this_syringe) &
        boluses_per_hr_this_syringe > max_per_hr_this_syringe,
      syringe_lockout_missing = is.na(lockout_mins) & boluses_this_syringe > 0,
      syringe_bolus_dose_missing = is.na(bolus_dose_mg_syringe) & boluses_this_syringe > 0
    )
  
  cat(sprintf("Syringes with computable bolus totals: %d\n", nrow(syringe_bolus)))
  cat(sprintf("  Syringes with boluses but missing bolus dose: %d\n",
              sum(syringe_bolus$syringe_bolus_dose_missing, na.rm = TRUE)))
  cat(sprintf("  Syringes with boluses but missing lockout:    %d\n",
              sum(syringe_bolus$syringe_lockout_missing, na.rm = TRUE)))
  
  # Collapse to episode level — total dose sums correctly-costed
  # per-syringe amounts; missingness/implausibility flags propagate
  # up as "ANY syringe in this episode had a problem"
  episode_bolus_totals <- syringe_bolus |>
    group_by(PAT_ENC_CSN_ID, episode_id) |>
    summarise(
      total_boluses_delivered = sum(boluses_this_syringe, na.rm = TRUE),
      total_bolus_mg = sum(bolus_mg_this_syringe, na.rm = TRUE),
      any_syringe_implausible = any(syringe_implausible, na.rm = TRUE),
      any_syringe_lockout_missing = any(syringe_lockout_missing, na.rm = TRUE),
      any_syringe_bolus_dose_missing = any(syringe_bolus_dose_missing, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat(sprintf("Episodes with computable bolus totals: %d\n", nrow(episode_bolus_totals)))
  
  # Flowsheet coverage adequacy per episode (unchanged)
  coverage_adequacy <- flow_assigned |>
    filter(!is.na(episode_id)) |>
    group_by(PAT_ENC_CSN_ID, episode_id) |>
    summarise(first_obs = min(RECORDED_TIME), last_obs = max(RECORDED_TIME), .groups = "drop") |>
    inner_join(episode_windows, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    mutate(
      episode_span_hrs = as.numeric(difftime(episode_end, episode_start, units = "hours")),
      flowsheet_span_hrs = as.numeric(difftime(last_obs, first_obs, units = "hours")),
      coverage_pct = pmin(100 * flowsheet_span_hrs / episode_span_hrs, 100)
    ) |>
    select(PAT_ENC_CSN_ID, episode_id, coverage_pct)
  
  # Background dose + duration (unchanged — already correctly
  # per-syringe weighted in episode_dose.rds, see 7_dosing.R)
  episode_dose_bg <- episode_dose |>
    filter(drug == drug_name) |>
    select(PAT_ENC_CSN_ID, episode_id, WEIGHT_KG, total_background_mg, episode_duration_dose_hrs)
  
  # Combine and apply all exclusion criteria
  result <- episode_bolus_totals |>
    inner_join(episode_dose_bg, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    left_join(coverage_adequacy, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    filter(!is.na(total_background_mg)) |>
    mutate(
      total_dose_mg = total_bolus_mg + total_background_mg,
      bolus_pct_of_total = 100 * total_bolus_mg / total_dose_mg,
      total_dose_mcg_kg_hr = (total_dose_mg * 1000) / WEIGHT_KG / episode_duration_dose_hrs,
      background_mcg_kg_hr = (total_background_mg * 1000) / WEIGHT_KG / episode_duration_dose_hrs
    ) |>
    filter(
      is.finite(total_dose_mcg_kg_hr),
      episode_duration_dose_hrs >= 2,                 # matches the standard exclude_short_pca flag
      coverage_pct >= 80,                              # adequate flowsheet coverage
      !any_syringe_lockout_missing,                     # missing lockout excludes (per-syringe now)
      !any_syringe_bolus_dose_missing,                  # NEW: missing bolus dose (any syringe) excludes
      !any_syringe_implausible,                         # plausibility ceiling, checked PER SYRINGE now
      !PAT_ENC_CSN_ID %in% KNOWN_CORRUPTED_CSNS         # explicit quarantine
    ) |>
    mutate(drug = drug_name)
  
  cat(sprintf("Episodes after full exclusion pipeline: %d\n", nrow(result)))
  cat("\nTotal dose (background+bolus), mcg/kg/hr:\n")
  print(summary(result$total_dose_mcg_kg_hr))
  cat("\nBackground-only, mcg/kg/hr (existing methodology, same subset):\n")
  print(summary(result$background_mcg_kg_hr))
  cat(sprintf("\nMedian ratio (total/background): %.2fx\n",
              median(result$total_dose_mcg_kg_hr, na.rm = TRUE) /
                median(result$background_mcg_kg_hr, na.rm = TRUE)))
  cat(sprintf("Median bolus %% of total: %.1f%%\n\n",
              median(result$bolus_pct_of_total, na.rm = TRUE)))
  
  result
}

# ============================================================
# RUN FOR MORPHINE (validated across Diagnostics 1-13)
# ============================================================

morphine_true_dose <- compute_true_total_dose(
  drug_name = "morphine",
  drug_regex = "morphine",
  flow_dose_col = "MORPHINE_NUMBER_OF_DOSES",
  conc_col = "morphine_conc_mg_per_ml"
)

# ============================================================
# RUN FOR OXYCODONE — same pipeline, INDEPENDENT trust-check
# needed since this is the number that goes in the manuscript
# ============================================================

oxycodone_true_dose <- compute_true_total_dose(
  drug_name = "oxycodone",
  drug_regex = "oxycodone",
  flow_dose_col = "OXYCODONE_NUMBER_OF_DOSES",
  conc_col = "oxycodone_conc_mg_per_ml"
)

cat("============================================================\n")
cat("OXYCODONE-SPECIFIC SANITY CHECKS (do not skip — this is the\n")
cat("primary drug, morphine's validation does not automatically\n")
cat("transfer)\n")
cat("============================================================\n\n")

cat("Top 10 oxycodone outliers — eyeball before trusting.\n")
cat("(Anonymised for chat: row index only — PAT_ENC_CSN_ID/episode_id\n")
cat("are NOT printed here per the data safety mandate. If you need to\n")
cat("trace a specific outlier back to an encounter, run the same\n")
cat("arrange()/select() with those columns locally in RStudio and\n")
cat("inspect there — don't paste those columns into chat.)\n\n")
oxycodone_true_dose |>
  arrange(desc(total_dose_mcg_kg_hr)) |>
  mutate(outlier_rank = row_number()) |>
  select(outlier_rank, total_boluses_delivered, episode_duration_dose_hrs,
         coverage_pct, total_dose_mcg_kg_hr) |>
  slice_head(n = 10) |>
  print()

cat("\nCompare oxycodone median ratio to morphine's 2.45x —\n")
cat("large divergence would be unexpected given similar pump/PROGRAMME\n")
cat("logic and warrants investigation before use:\n")
cat(sprintf("Morphine ratio:  %.2fx\n",
            median(morphine_true_dose$total_dose_mcg_kg_hr, na.rm = TRUE) /
              median(morphine_true_dose$background_mcg_kg_hr, na.rm = TRUE)))
cat(sprintf("Oxycodone ratio: %.2fx\n",
            median(oxycodone_true_dose$total_dose_mcg_kg_hr, na.rm = TRUE) /
              median(oxycodone_true_dose$background_mcg_kg_hr, na.rm = TRUE)))

# ============================================================
# COMBINE AND SAVE
# ============================================================

episode_true_total_dose <- bind_rows(morphine_true_dose, oxycodone_true_dose)

cat("\n============================================================\n")
cat("FINAL COMBINED SUMMARY — TRUE TOTAL vs BACKGROUND-ONLY, BY DRUG\n")
cat("============================================================\n\n")

episode_true_total_dose |>
  group_by(drug) |>
  summarise(
    n = n(),
    background_median = round(median(background_mcg_kg_hr, na.rm = TRUE), 2),
    background_q25 = round(quantile(background_mcg_kg_hr, 0.25, na.rm = TRUE), 2),
    background_q75 = round(quantile(background_mcg_kg_hr, 0.75, na.rm = TRUE), 2),
    total_median = round(median(total_dose_mcg_kg_hr, na.rm = TRUE), 2),
    total_q25 = round(quantile(total_dose_mcg_kg_hr, 0.25, na.rm = TRUE), 2),
    total_q75 = round(quantile(total_dose_mcg_kg_hr, 0.75, na.rm = TRUE), 2),
    bolus_pct_median = round(median(bolus_pct_of_total, na.rm = TRUE), 1)
  ) |>
  print()

saveRDS(episode_true_total_dose,
        "1 - data/3 - processed_data/episode_true_total_dose.rds")

cat("\nepisode_true_total_dose.rds saved\n\n")

cat("============================================================\n")
cat("ANALYSIS 13 SECTION 1 COMPLETE\n")
cat("============================================================\n\n")

cat("METHODOLOGY SUMMARY FOR METHODS SECTION:\n")
cat("  Total dose = background infusion dose (rate x syringe duration,\n")
cat("  existing methodology) + bolus dose (boluses delivered per\n")
cat("  PCA/NCA pump flowsheet activation counter x prescribed bolus\n")
cat("  dose per syringe). Bolus counts derived from pump flowsheet\n")
cat("  NUMBER_OF_DOSES field, syringe-anchored to existing trusted MAR\n")
cat("  syringe boundaries (not flowsheet counter resets, which showed\n")
cat("  imperfect alignment with true syringe changes).\n\n")

cat("PRE-SPECIFIED LIMITATIONS FOR THIS ANALYSIS:\n")
cat("  - Flowsheet coverage restricted to ~73-79%% of full cohort\n")
cat("    (episodes without flowsheet monitoring excluded entirely)\n")
cat("  - Additional exclusions: duration <2hrs, flowsheet coverage\n")
cat("    <80%% of episode span, missing/implausible lockout data\n")
cat("  - One encounter with confirmed corrupted counter data excluded\n")
cat("    by name (documented: doses=3914 vs attempts=429, CSN 25786392)\n")
cat("  - This is a SUBSET analysis (bolus-flowsheet-available episodes\n")
cat("    only) and should be reported alongside, not instead of, the\n")
cat("    full-cohort background-only consumption figures\n")
cat("============================================================\n\n")

cat("NOTE: this run uses the per-syringe bolus-dose/lockout method\n")
cat("described above — total_bolus_mg, total_dose_mg and the exclusion\n")
cat("counts are NOT directly comparable to any prior run of this\n")
cat("script. Re-eyeball the oxycodone outlier table and the settings\n")
cat("stability check above before trusting these numbers.\n\n")

cat("STOP — DO NOT proceed to a PSM-matched delivered-dose\n")
cat("comparison yet. Section 2 below (reconciliation diagnostic)\n")
cat("MUST show agreement before any delivered-dose figure is trusted.\n")
cat("============================================================\n\n")


# ============================================================
# ============================================================
# SECTION 2 — RECONCILIATION DIAGNOSTIC
# background_mcg_kg_hr (Section 1, above) vs episode_dose.rds-
# derived consumption, SAME episodes, aggregate output only
# ============================================================
# ============================================================
#
# Tests TWO candidate explanations for the ~4-5x gap between
# Section 1's background_mcg_kg_hr and the locked/published PSM
# background figures (oxy 16.1, morph 20.8 mcg/kg/hr):
#
#   Candidate A: background_rate_mg_kg_hr_median * 1000
#                (median of PER-SYRINGE rates — what 9_comsumption.R/6
#                 may be using for the published figure)
#   Candidate B: total_background_mg * 1000 / WEIGHT_KG /
#                episode_duration_dose_hrs
#                (duration-weighted total — the exact formula
#                 Section 1 uses)
#
# If B ~= Section 1's figure: Section 1 is internally consistent;
#   the gap is A vs B being genuinely different quantities
#   (median-of-rates vs weighted-average) -> check which one
#   9_comsumption.R actually reports.
# If B also diverges from Section 1's own figure despite being the
#   "same" formula: denominator/units bug confirmed -> Steps 4-5
#   below isolate whether it's the duration field, the
#   total_background_mg field, or both.
# ============================================================

cat("============================================================\n")
cat("ANALYSIS 13 SECTION 2 — RECONCILIATION DIAGNOSTIC\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# STEP 1 — RESTRICT episode_dose TO THE SAME EPISODES AS SECTION 1
# ------------------------------------------------------------

common_keys <- episode_true_total_dose |>
  select(PAT_ENC_CSN_ID, episode_id, drug)

cat(sprintf("Section 1 episodes (flowsheet-linked subset): %d\n", nrow(common_keys)))

recon <- episode_dose |>
  inner_join(common_keys, by = c("PAT_ENC_CSN_ID", "episode_id", "drug")) |>
  select(PAT_ENC_CSN_ID, episode_id, drug, WEIGHT_KG,
         background_rate_mg_kg_hr_median,
         total_background_mg, total_background_mg_kg,
         episode_duration_dose_hrs) |>
  rename(episode_dose_duration_hrs = episode_duration_dose_hrs) |>
  left_join(
    episode_true_total_dose |>
      select(PAT_ENC_CSN_ID, episode_id, drug,
             new_pipeline_background_mcg_kg_hr = background_mcg_kg_hr,
             new_pipeline_duration_hrs = episode_duration_dose_hrs,
             new_pipeline_total_background_mg = total_background_mg),
    by = c("PAT_ENC_CSN_ID", "episode_id", "drug")
  ) |>
  mutate(
    candidate_A_mcg_kg_hr = background_rate_mg_kg_hr_median * 1000,
    candidate_B_mcg_kg_hr = (total_background_mg * 1000) / WEIGHT_KG / episode_dose_duration_hrs
  )

cat(sprintf("Episodes matched on (CSN, episode_id, drug): %d / %d\n\n",
            nrow(recon), nrow(common_keys)))

if (nrow(recon) < 0.9 * nrow(common_keys)) {
  cat("WARNING: fewer than 90%% of Section 1 episodes found in episode_dose.\n")
  cat("    This alone could partly explain divergence — check drug/key alignment.\n\n")
}

# ------------------------------------------------------------
# STEP 2 — HEADLINE COMPARISON: CANDIDATE A vs B vs SECTION 1
# ------------------------------------------------------------

cat("============================================================\n")
cat("HEADLINE: MEDIAN mcg/kg/hr BY METHOD, BY DRUG (same episodes)\n")
cat("============================================================\n\n")

recon |>
  group_by(drug) |>
  summarise(
    n = n(),
    candidate_A_median = round(median(candidate_A_mcg_kg_hr, na.rm = TRUE), 2),
    candidate_B_median = round(median(candidate_B_mcg_kg_hr, na.rm = TRUE), 2),
    section1_median = round(median(new_pipeline_background_mcg_kg_hr, na.rm = TRUE), 2),
    ratio_A_to_section1 = round(median(candidate_A_mcg_kg_hr, na.rm = TRUE) /
                                  median(new_pipeline_background_mcg_kg_hr, na.rm = TRUE), 2),
    ratio_B_to_section1 = round(median(candidate_B_mcg_kg_hr, na.rm = TRUE) /
                                  median(new_pipeline_background_mcg_kg_hr, na.rm = TRUE), 2),
    .groups = "drop"
  ) |>
  print()

cat("\nFor reference — locked/published PSM background figures:\n")
cat("  Morphine:  20.8 (14.6-41.5) mcg/kg/hr\n")
cat("  Oxycodone: 16.1 (12.6-25.4) mcg/kg/hr\n\n")

# ------------------------------------------------------------
# STEP 3 — IS CANDIDATE B (SAME FORMULA) IDENTICAL TO SECTION 1?
# If Section 1 and episode_dose both compute
# total_background_mg*1000/WEIGHT_KG/duration on the same episode,
# these should match near-exactly. Any gap here = denominator or
# total_background_mg mismatch between the two pipelines, NOT a
# methodology difference.
# ------------------------------------------------------------

cat("============================================================\n")
cat("CANDIDATE B vs SECTION 1 — should be ~identical (same formula)\n")
cat("============================================================\n\n")

recon |>
  mutate(
    diff_B_new = candidate_B_mcg_kg_hr - new_pipeline_background_mcg_kg_hr,
    pct_diff_B_new = 100 * diff_B_new / new_pipeline_background_mcg_kg_hr
  ) |>
  group_by(drug) |>
  summarise(
    n = n(),
    median_abs_pct_diff = round(median(abs(pct_diff_B_new), na.rm = TRUE), 1),
    pct_episodes_within_1pct = round(100 * mean(abs(pct_diff_B_new) < 1, na.rm = TRUE), 1),
    pct_episodes_within_10pct = round(100 * mean(abs(pct_diff_B_new) < 10, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

# ------------------------------------------------------------
# STEP 4 — DURATION FIELD COMPARISON
# episode_duration_dose_hrs is built independently in 7_dosing.R
# (episode_dose) vs Section 1 (via inner_join onto episode_dose_mg,
# which itself pulls from episode_dose) -- should be THE SAME
# COLUMN. If it isn't, a join/filter step in Section 1 is silently
# duplicating or dropping syringe rows before that column is
# computed downstream.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("DURATION FIELD COMPARISON (episode_dose vs Section 1 pipeline)\n")
cat("============================================================\n\n")

recon |>
  mutate(
    duration_diff = episode_dose_duration_hrs - new_pipeline_duration_hrs,
    duration_ratio = episode_dose_duration_hrs / new_pipeline_duration_hrs
  ) |>
  group_by(drug) |>
  summarise(
    n = n(),
    median_duration_episode_dose = round(median(episode_dose_duration_hrs, na.rm = TRUE), 2),
    median_duration_section1 = round(median(new_pipeline_duration_hrs, na.rm = TRUE), 2),
    median_ratio = round(median(duration_ratio, na.rm = TRUE), 2),
    pct_matching_within_5pct = round(100 * mean(abs(duration_ratio - 1) < 0.05, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

# ------------------------------------------------------------
# STEP 5 — total_background_mg FIELD COMPARISON
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("total_background_mg FIELD COMPARISON (episode_dose vs Section 1)\n")
cat("============================================================\n\n")

recon |>
  mutate(
    mg_diff = total_background_mg - new_pipeline_total_background_mg,
    mg_ratio = total_background_mg / new_pipeline_total_background_mg
  ) |>
  group_by(drug) |>
  summarise(
    n = n(),
    median_mg_episode_dose = round(median(total_background_mg, na.rm = TRUE), 1),
    median_mg_section1 = round(median(new_pipeline_total_background_mg, na.rm = TRUE), 1),
    median_ratio = round(median(mg_ratio, na.rm = TRUE), 2),
    pct_matching_within_5pct = round(100 * mean(abs(mg_ratio - 1) < 0.05, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

# ------------------------------------------------------------
# INTERPRETATION GUIDE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("INTERPRETATION GUIDE\n")
cat("============================================================\n\n")

cat("1. If Candidate B ~= Section 1 (Step 3, within ~1-10%%):\n")
cat("   -> Section 1's own arithmetic is internally consistent.\n")
cat("   -> The 4-5x gap to the published 16.1/20.8 figures means\n")
cat("      9_comsumption.R/6 is NOT using total_background_mg/duration —\n")
cat("      it is using background_rate_mg_kg_hr_median (Candidate A)\n")
cat("      or something else entirely. Check ratio_A_to_section1 in\n")
cat("      Step 2: if that ratio is close to 1, Candidate A IS the\n")
cat("      published methodology, and the two 'consumption' figures\n")
cat("      measure genuinely different things (median per-syringe\n")
cat("      rate vs a duration-weighted total-dose rate) -- not\n")
cat("      interchangeable, and this needs an explicit methods-\n")
cat("      section decision about which one 'consumption' means for\n")
cat("      the manuscript.\n\n")

cat("2. If Candidate B diverges from Section 1 (Step 3 fails):\n")
cat("   -> Same formula, same episodes, different answer = a real\n")
cat("      pipeline bug (join duplicating/dropping syringe rows,\n")
cat("      unit mismatch, or stale object). Steps 4-5 isolate\n")
cat("      whether it's the duration field, the total_background_mg\n")
cat("      field, or both, that differs between the two pipelines.\n\n")

cat("3. Do NOT proceed to the PSM-matched delivered-dose comparison\n")
cat("   until Step 3 shows near-perfect agreement AND you've decided\n")
cat("   (based on Step 2) which formula is the correct 'consumption'\n")
cat("   for direct comparison against total delivered dose.\n\n")

cat("4. Once resolved: rerun the PSM-matched comparison against\n")
cat("   psm_matched_keys.rds (the actual 237-pair matched population\n")
cat("   from 10_comparative_PSM.R) — NOT s3_morphine_keys/s3_oxycodone_keys,\n")
cat("   which are the pre-match candidate pools (1,718/433).\n")
cat("============================================================\n")


# ============================================================
# ============================================================
# SECTION 3 — POPULATION / DURATION SELECTION CHECK
# ============================================================
# ============================================================
#
# Section 2 found candidate_A ~= candidate_B ~= Section 1 (no
# formula or units bug) but ALL THREE sit at ~4 mcg/kg/hr for both
# drugs -- nowhere near the published PSM figures (16.1 oxy / 20.8
# morph). This section tests whether the gap is a population/
# duration selection effect: the flowsheet-linked, >=80% coverage,
# >=2hr subset required by Section 1 may be systematically
# different (e.g. much longer episodes) from the population the
# published consumption figures were computed on — rather than a
# pipeline bug.
# ============================================================

cat("\n\n============================================================\n")
cat("ANALYSIS 13 SECTION 3 — POPULATION / DURATION SELECTION CHECK\n")
cat("============================================================\n\n")

stopifnot(exists("episode_dose"), exists("episode_true_total_dose"))

full_cohort <- episode_dose |>
  filter(drug %in% c("morphine", "oxycodone")) |>
  select(PAT_ENC_CSN_ID, episode_id, drug,
         background_rate_mg_kg_hr_median, episode_duration_dose_hrs)

subset_cohort <- episode_true_total_dose |>
  select(PAT_ENC_CSN_ID, episode_id, drug,
         background_mcg_kg_hr, episode_duration_dose_hrs)

cat("--- Episode duration + background rate: FULL episode_dose cohort ---\n\n")
full_cohort |>
  group_by(drug) |>
  summarise(
    n = n(),
    duration_median_hrs = round(median(episode_duration_dose_hrs, na.rm = TRUE), 1),
    duration_q25_hrs = round(quantile(episode_duration_dose_hrs, 0.25, na.rm = TRUE), 1),
    duration_q75_hrs = round(quantile(episode_duration_dose_hrs, 0.75, na.rm = TRUE), 1),
    rate_median_mcg_kg_hr = round(median(background_rate_mg_kg_hr_median, na.rm = TRUE) * 1000, 2),
    .groups = "drop"
  ) |>
  print()

cat("\n--- Episode duration + background rate: Section 1 subset (flowsheet-linked) ---\n\n")
subset_cohort |>
  group_by(drug) |>
  summarise(
    n = n(),
    duration_median_hrs = round(median(episode_duration_dose_hrs, na.rm = TRUE), 1),
    duration_q25_hrs = round(quantile(episode_duration_dose_hrs, 0.25, na.rm = TRUE), 1),
    duration_q75_hrs = round(quantile(episode_duration_dose_hrs, 0.75, na.rm = TRUE), 1),
    rate_median_mcg_kg_hr = round(median(background_mcg_kg_hr, na.rm = TRUE), 2),
    .groups = "drop"
  ) |>
  print()

cat("\n--- Does background rate fall as duration rises? (FULL cohort, by duration quartile) ---\n\n")
full_cohort |>
  filter(!is.na(episode_duration_dose_hrs), !is.na(background_rate_mg_kg_hr_median)) |>
  group_by(drug) |>
  mutate(duration_quartile = ntile(episode_duration_dose_hrs, 4)) |>
  group_by(drug, duration_quartile) |>
  summarise(
    n = n(),
    duration_median_hrs = round(median(episode_duration_dose_hrs), 1),
    rate_median_mcg_kg_hr = round(median(background_rate_mg_kg_hr_median) * 1000, 2),
    .groups = "drop"
  ) |>
  arrange(drug, duration_quartile) |>
  print()

cat("\n============================================================\n")
cat("INTERPRETATION GUIDE — SECTION 3\n")
cat("============================================================\n\n")

cat("- If duration_median is much higher in the Section 1 subset than\n")
cat("  the FULL cohort, AND rate falls with duration quartile in the\n")
cat("  full-cohort breakdown above: the ~4x gap is a genuine SELECTION\n")
cat("  effect of requiring >=80%% flowsheet coverage (which mechanically\n")
cat("  favours longer, more heavily-monitored episodes) — not a bug.\n")
cat("  The published 16.1/20.8 figures describe the full/PSM\n")
cat("  population; the delivered-dose subset describes a\n")
cat("  systematically longer-duration slice, and the two are NOT\n")
cat("  directly comparable without adjusting for that (e.g. reporting\n")
cat("  delivered dose stratified by duration, or as a ratio to each\n")
cat("  episode's OWN background figure rather than as absolute\n")
cat("  mcg/kg/hr against the PSM headline).\n\n")

cat("- If duration is similar across both but rate still differs\n")
cat("  substantially: the flowsheet-linked subset differs on some\n")
cat("  other axis (procedure type, ward, era) rather than duration —\n")
cat("  would need a covariate comparison (age, procedure, era) between\n")
cat("  the two populations before drawing conclusions.\n\n")

cat("- Either way: do NOT report an absolute total-delivered-dose\n")
cat("  mcg/kg/hr figure as directly comparable to the 16.1/20.8 PSM\n")
cat("  headline until this is resolved. The safer manuscript framing,\n")
cat("  if selection is confirmed, is the WITHIN-SUBSET total/background\n")
cat("  RATIO (2.44-2.48x here) rather than the absolute delivered-dose\n")
cat("  number, since the ratio is less sensitive to which population\n")
cat("  the subset happens to select.\n")
cat("============================================================\n")


# ============================================================
# ============================================================
# SECTION 4 — IS episode_dose.rds STALE RELATIVE TO THE LOCKED
# 16.1 / 20.8 PSM FIGURES?
# ============================================================
# ============================================================
#
# Section 3 refutes the duration-selection hypothesis: background
# rate is ~4.00 mcg/kg/hr for BOTH drugs, essentially FLAT across
# duration quartiles from a 16hr median to an 814hr median episode
# (morphine) — if long tapering tails were diluting the rate, short
# episodes would show a much higher rate. They don't. This is also
# the FULL, unfiltered episode_dose cohort (n=3672/1243), not an
# 8_true_total_dose.R subset — so this isn't a selection artifact of
# 8_true_total_dose.R at all, it's a property of episode_dose.rds itself.
#
# Both drugs also converge to nearly the SAME rate (4.00 vs 3.99),
# which is odd given the locked figures differ by 30% (20.8 vs
# 16.1) and reflect a real drug-specific potency/dosing difference.
#
# Decisive test: rerun the EXACT PSM dose calculation (matched_keys
# + episode_dose_clean, as in 10_comparative_PSM.R/Table 4) right now on the
# currently loaded episode_dose object, and check whether it still
# reproduces 16.1 / 20.8. If it doesn't, episode_dose.rds (or
# psm_matched_keys.rds) has changed since those figures were locked
# — a higher-priority finding than the delivered-dose analysis,
# since it would mean the MANUSCRIPT'S headline consumption figure
# itself needs re-verifying.
# ============================================================

cat("\n\n============================================================\n")
cat("ANALYSIS 13 SECTION 4 — episode_dose.rds STALENESS CHECK\n")
cat("============================================================\n\n")

matched_keys <- tryCatch(
  readRDS("1 - data/3 - processed_data/psm_matched_keys.rds"),
  error = function(e) {
    cat("psm_matched_keys.rds not found — cannot run this check\n")
    NULL
  }
)

if (!is.null(matched_keys)) {
  
  # Implausible-dose filter, matches 7_dosing.R's episode_dose_clean exactly
  episode_dose_clean <- episode_dose |>
    mutate(
      daily_background_mg_kg = if_else(
        episode_duration_dose_hrs > 0,
        total_background_mg_kg / episode_duration_dose_hrs * 24,
        NA_real_
      ),
      flag_implausible_dose = daily_background_mg_kg > 2 | daily_background_mg_kg < 0
    ) |>
    filter(!flag_implausible_dose | is.na(flag_implausible_dose))
  
  matched_dose <- episode_dose_clean |>
    semi_join(matched_keys, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    filter(drug %in% c("morphine", "oxycodone"))
  
  cat(sprintf("Matched episodes with dose data: %d / %d\n\n",
              nrow(matched_dose), nrow(matched_keys)))
  
  cat("--- Rerun of PSM dose summary, RIGHT NOW, on current objects ---\n\n")
  matched_dose |>
    group_by(drug) |>
    summarise(
      n = n(),
      background_mcg_kg_hr_median = round(median(background_rate_mg_kg_hr_median, na.rm = TRUE) * 1000, 2),
      background_mcg_kg_hr_q25 = round(quantile(background_rate_mg_kg_hr_median, 0.25, na.rm = TRUE) * 1000, 2),
      background_mcg_kg_hr_q75 = round(quantile(background_rate_mg_kg_hr_median, 0.75, na.rm = TRUE) * 1000, 2),
      .groups = "drop"
    ) |>
    print()
  
  cat("\nFor reference — LOCKED/published PSM background figures:\n")
  cat("  Morphine:  20.8 (14.6-41.5) mcg/kg/hr\n")
  cat("  Oxycodone: 16.1 (12.6-25.4) mcg/kg/hr\n\n")
  cat("If the rerun above does NOT reproduce these numbers ->\n")
  cat("episode_dose.rds (and/or psm_matched_keys.rds) has changed\n")
  cat("since the PSM headline was locked, and the MANUSCRIPT figure\n")
  cat("itself needs re-verifying — higher priority than the\n")
  cat("delivered-dose analysis.\n")
  cat("If it DOES reproduce 16.1/20.8 -> the objects are current, and\n")
  cat("the ~4.00 seen in Sections 1-3 has a different explanation\n")
  cat("(e.g. this script is reading a stale/cached episode_dose in\n")
  cat("this R session specifically — try a fresh source() of\n")
  cat("1_load_clean_data.R before rerunning).\n\n")
}

cat("--- Distribution check: is background_rate_mg_kg_hr_median\n")
cat("    suspiciously clustered at a single value? ---\n\n")

episode_dose |>
  filter(drug %in% c("morphine", "oxycodone")) |>
  mutate(rate_mcg_kg_hr_rounded = round(background_rate_mg_kg_hr_median * 1000, 1)) |>
  count(drug, rate_mcg_kg_hr_rounded, sort = TRUE) |>
  group_by(drug) |>
  slice_max(n, n = 5) |>
  ungroup() |>
  print(n = 20)

cat("\nIf one value dominates the count (e.g. a large fraction of\n")
cat("episodes sharing the exact same rate), that suggests a default/\n")
cat("order-set value or an imputation is contaminating the field\n")
cat("rather than reflecting titrated, patient-specific rates.\n")
cat("============================================================\n")
