# ============================================================
# 8 — True total delivered dose (background + bolus)
# Paediatric PCA/NCA study — Addenbrooke's Hospital (CUH NHS FT)
# Author: J Chin
# ============================================================
#
# Purpose
#   Computes the total delivered opioid dose per episode: background infusion
#   plus the boluses actually delivered, using the pump flowsheet
#   NUMBER_OF_DOSES counter linked to episodes by time-interval containment.
#
# Inputs
#   In session: pca_mar_clean, pca_episodes_flagged, episode_dose (from 7_dosing.R).
#   Raw file (project root): "Chin Paeds PCA Datasets - PCA Flowsheets - refresh
#     2026-07-15.xlsx" (read directly; change FLOWSHEET_FILE for a newer refresh).
#
# Output
#   episode_true_total_dose.rds
#
# Method
#   - Bolus total = sum across syringes of (boluses delivered on that syringe,
#     corrected for counter resets using true syringe boundaries) x (that same
#     syringe's own prescribed bolus dose in mg, from the PROGRAMME text on
#     pca_mar_clean). Bolus dose and lockout are therefore per syringe, not one
#     episode-level value.
#   - The flowsheet CUMULATIVE_DOSE / CUMULATIVE_INTAKE fields are NOT used: for
#     morphine CUMULATIVE_DOSE is a running total, but for oxycodone
#     CUMULATIVE_INTAKE is a period value, so the two are not equivalent.
#     NUMBER_OF_DOSES is a monotonic bolus counter for both drugs and is used
#     instead.
#   - Syringe and episode boundaries come from the MAR (TAKEN_TIME /
#     Discontinue_Time), not from flowsheet counter resets, which do not
#     reliably align with real syringe changes.
#
# Exclusions
#   1. Episode duration < 2 h (matches exclude_short_pca in the PSM population).
#   2. Flowsheet coverage < 80% of the episode span.
#   3. Lockout-implied plausibility: boluses per hour must not exceed
#      60 / lockout_mins, checked per syringe against that syringe's own
#      lockout. Missing lockout or bolus dose on any syringe with boluses
#      excludes the whole episode (it does not silently pass).
#   4. Corrupted pump counters: any encounter with a NUMBER_OF_DOSES reading more
#      than COUNTER_RATIO_LIMIT (2) times the matching attempts counter. The rule uses
#      the data, so no encounter ID is stored anywhere; counts are printed for k = 1,
#      1.5, 2 and 5.
#
# Data protection
#   No PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME, PROGRAMME or MED_ORDER content is
#   printed; output is aggregate only (counts, percentages, medians, IQR).
#
# Coverage
#   Only episodes present in the flowsheet can be analysed, so this is a SUBSET
#   analysis: report it alongside, not instead of, the full-cohort
#   background-only consumption figures. Coverage is printed at run time.
# ============================================================

library(tidyverse)
library(readxl)

stopifnot(exists("pca_mar_clean"), exists("pca_episodes_flagged"), exists("episode_dose"))


cat("============================================================\n")
cat("SCRIPT 8 — TRUE TOTAL DOSE (BACKGROUND + BOLUS)\n")
cat("============================================================\n\n")

# ============================================================
# STEP 0 — LOAD FLOWSHEET
# ============================================================

cat("--- Loading flowsheet ---\n\n")

# Raw flowsheet export (project root). Change this filename for a newer refresh.
FLOWSHEET_FILE <- "Chin Paeds PCA Datasets - PCA Flowsheets - refresh 2026-07-15.xlsx"

flow_raw <- read_excel(FLOWSHEET_FILE) |> suppressWarnings()

flow <- flow_raw |>
  rename_with(~ str_replace_all(.x, " ", "_")) |>
  rename_with(~ str_to_upper(.x))

cat(sprintf("Flowsheet rows: %s\n\n", format(nrow(flow), big.mark = ",")))

# ------------------------------------------------------------
# COUNTER INTEGRITY CHECK
# A PCA/NCA pump cannot deliver more boluses than were demanded, so a
# NUMBER_OF_DOSES reading far above the matching attempts counter means the
# counter is corrupted. Encounters with such a reading are excluded below. The
# limit is deliberately generous (a legitimate counter has doses <= attempts), and
# only counts are printed, never identifiers.
# ------------------------------------------------------------
COUNTER_RATIO_LIMIT <- 2   # exclude if NUMBER_OF_DOSES > COUNTER_RATIO_LIMIT x attempts

find_attempts_col <- function(dose_col) {
  prefix <- sub("_NUMBER_OF_DOSES$", "", dose_col)
  attempt_cols <- names(flow)[str_detect(names(flow), "ATTEMPT")]
  drug_specific <- attempt_cols[str_detect(attempt_cols, fixed(prefix))]
  generic <- attempt_cols[!str_detect(attempt_cols, "MORPHINE|OXYCODONE")]
  chosen <- if (length(drug_specific) == 1) drug_specific
            else if (length(drug_specific) == 0 && length(generic) == 1) generic
            else character(0)
  if (length(chosen) != 1) {
    stop("Cannot identify a single attempts column for ", dose_col,
         ". Flowsheet columns containing 'ATTEMPT': ",
         if (length(attempt_cols)) paste(attempt_cols, collapse = ", ") else "(none)",
         ". Edit find_attempts_col() to name the right one.")
  }
  chosen
}

counter_readings <- map_dfr(
  c("MORPHINE_NUMBER_OF_DOSES", "OXYCODONE_NUMBER_OF_DOSES"),
  function(dose_col) {
    attempts_col <- find_attempts_col(dose_col)
    tibble(
      PAT_ENC_CSN_ID = flow$PAT_ENC_CSN_ID,
      doses    = as.numeric(flow[[dose_col]]),
      attempts = as.numeric(flow[[attempts_col]])
    )
  }
) |>
  filter(!is.na(doses), !is.na(attempts))

cat("Encounters with a NUMBER_OF_DOSES reading above k x the attempts counter:\n")
for (k in c(1, 1.5, 2, 5)) {
  cat(sprintf("  k = %-3g : %d encounter(s)\n", k,
              n_distinct(counter_readings$PAT_ENC_CSN_ID[counter_readings$doses > k * counter_readings$attempts])))
}

counter_corrupt_csns <- counter_readings |>
  filter(doses > COUNTER_RATIO_LIMIT * attempts) |>
  distinct(PAT_ENC_CSN_ID) |>
  pull(PAT_ENC_CSN_ID)
cat(sprintf("Excluded as corrupted (k = %g): %d encounter(s)\n\n",
            COUNTER_RATIO_LIMIT, length(counter_corrupt_csns)))

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
      !PAT_ENC_CSN_ID %in% counter_corrupt_csns          # corrupted pump counter (see COUNTER INTEGRITY CHECK)
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
cat("(Row index only — PAT_ENC_CSN_ID/episode_id are deliberately not printed.\n")
cat("To trace an outlier back to an encounter, run the same arrange()/select()\n")
cat("with those columns interactively on the secure desktop.)\n\n")
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
cat("SCRIPT 8 COMPLETE\n")
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
cat("  - Flowsheet coverage restricted to ~73-79% of full cohort\n")
cat("    (episodes without flowsheet monitoring excluded entirely)\n")
cat("  - Additional exclusions: duration <2hrs, flowsheet coverage\n")
cat("    <80% of episode span, missing/implausible lockout data\n")
cat(sprintf("  - %d encounter(s) with a corrupted pump counter excluded\n",
            length(counter_corrupt_csns)))
cat(sprintf("    (NUMBER_OF_DOSES more than %g x the attempts counter)\n", COUNTER_RATIO_LIMIT))
cat("  - This is a SUBSET analysis (bolus-flowsheet-available episodes\n")
cat("    only) and should be reported alongside, not instead of, the\n")
cat("    full-cohort background-only consumption figures\n")
cat("============================================================\n\n")
