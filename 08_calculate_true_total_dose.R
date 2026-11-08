# ============================================================
# 08 - Total delivered dose (background + bolus)
# ============================================================
# Total dose per episode = background infusion dose (07_calculate_doses.R)
# + delivered boluses, from the PCA flowsheet extract:
#
#   - Boluses delivered come from the NUMBER_OF_DOSES counter, which is a
#     monotonic bolus count for both drugs. The flowsheet cumulative dose
#     fields are not used: for morphine the field is a running total,
#     for oxycodone it is a period value, so the two are not equivalent.
#   - Bolus total = sum over syringes of (boluses delivered on that
#     syringe) x (that syringe's own prescribed bolus dose in mg).
#   - Syringe and episode boundaries come from the MAR start and stop
#     times, not from counter resets, which do not reliably coincide with
#     syringe changes.
#
# Episodes are excluded if:
#   1. duration is under 2 hours
#   2. flowsheet coverage is under 80% of the episode span
#   3. boluses per hour exceed 60 / lockout minutes on any syringe
#   4. a syringe with boluses has no lockout or no bolus dose recorded
#   5. the pump counter is known to be corrupted (encounter list read
#      from a local file, see below)
#
# The analysis is therefore restricted to episodes with flowsheet data.
# Output: episode_true_total_dose.rds (aggregate output only is printed).
# ============================================================

library(tidyverse)
library(readxl)

stopifnot(exists("pca_mar_clean"), exists("pca_episodes_flagged"), exists("episode_dose"))

# Encounters with a confirmed corrupted pump dose counter (doses delivered
# far exceeding attempts) are excluded by encounter ID. The IDs are not
# stored in this repository: they are read from a local file (one
# PAT_ENC_CSN_ID per line) in PCA_DATA_DIR/manual_corrections/. If the
# file is absent, none are excluded.
corrupted_csn_file <- file.path(Sys.getenv("PCA_DATA_DIR"), "manual_corrections",
                                "corrupted_counter_csns.txt")
KNOWN_CORRUPTED_CSNS <- if (file.exists(corrupted_csn_file)) {
  suppressWarnings(as.numeric(readLines(corrupted_csn_file, warn = FALSE)))
} else numeric(0)
KNOWN_CORRUPTED_CSNS <- KNOWN_CORRUPTED_CSNS[!is.na(KNOWN_CORRUPTED_CSNS)]

# ============================================================
# Load the flowsheet extract
# ============================================================

flow_file <- file.path(
  Sys.getenv("PCA_DATA_DIR"),
  paste0(Sys.getenv("PCA_FILE_PREFIX", unset = "PCA Datasets - "),
         "PCA Flowsheets",
         Sys.getenv("PCA_FILE_SUFFIX", unset = " - refresh.xlsx"))
)
stopifnot("PCA flowsheet extract not found - check PCA_DATA_DIR / PCA_FILE_PREFIX / PCA_FILE_SUFFIX" =
            file.exists(flow_file))
flow_raw <- read_excel(flow_file) |> suppressWarnings()

flow <- flow_raw |>
  rename_with(~ str_replace_all(.x, " ", "_")) |>
  rename_with(~ str_to_upper(.x))

cat(sprintf("Flowsheet rows: %s\n\n", format(nrow(flow), big.mark = ",")))

# ============================================================
# Helper functions
# ============================================================
# Bolus dose and lockout are taken per syringe from the PROGRAMME text on
# pca_mar_clean, not collapsed to one value per episode. Boluses
# delivered on a syringe are costed at that syringe's own prescribed dose
# and checked against its own lockout, so a mid-episode change in
# settings is handled correctly. The background dose already comes from
# a per-syringe sum (total_background_mg in episode_dose.rds).

# Link each syringe to an episode by time-interval containment
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

# Attach each syringe's own prescribed bolus dose (mg) and lockout
# (minutes). conc_col is the drug-specific concentration column on
# pca_mar_clean, extracted per row from MED_ORDER.
attach_syringe_settings <- function(syringe_intervals, drug_regex, conc_col) {
  settings <- pca_mar_clean |>
    filter(str_detect(MEDICATION, regex(drug_regex, ignore_case = TRUE))) |>
    mutate(bolus_dose_mg_syringe = bolus_dose_ml * .data[[conc_col]]) |>
    distinct(PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME,
             bolus_dose_mg_syringe, lockout_mins)
  
  syringe_intervals |>
    left_join(settings, by = c("PAT_ENC_CSN_ID", "ORDER_MED_ID", "TAKEN_TIME"))
}

# How often bolus dose and lockout change within an episode (aggregate only)
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

# Assign flowsheet observations to a syringe and episode
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

# Total delivered dose for one drug
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
  
  # Syringe duration, for the per-syringe plausibility check
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
  
  # Boluses delivered per syringe (change in the pump dose counter),
  # joined to that syringe's own bolus dose and lockout
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
  
  # Collapse to episode level; a missing-data or plausibility flag on any
  # syringe flags the episode
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
  
  # Flowsheet coverage per episode
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
  
  # Background dose and duration from episode_dose.rds (07_calculate_doses.R)
  episode_dose_bg <- episode_dose |>
    filter(drug == drug_name) |>
    select(PAT_ENC_CSN_ID, episode_id, WEIGHT_KG, total_background_mg, episode_duration_dose_hrs)
  
  # Combine and apply the exclusion criteria
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
      episode_duration_dose_hrs >= 2,                  # minimum duration, as exclude_short_pca
      coverage_pct >= 80,                              # flowsheet covers at least 80% of the episode
      !any_syringe_lockout_missing,                    # bolus counted but lockout missing
      !any_syringe_bolus_dose_missing,                 # bolus counted but bolus dose missing
      !any_syringe_implausible,                        # boluses per hour above the lockout limit
      !PAT_ENC_CSN_ID %in% KNOWN_CORRUPTED_CSNS        # corrupted pump counter
    ) |>
    mutate(drug = drug_name)
  
  cat(sprintf("Episodes after full exclusion pipeline: %d\n", nrow(result)))
  cat("\nTotal dose (background+bolus), mcg/kg/hr:\n")
  print(summary(result$total_dose_mcg_kg_hr))
  cat("\nBackground only, mcg/kg/hr (same episodes):\n")
  print(summary(result$background_mcg_kg_hr))
  cat(sprintf("\nMedian ratio (total/background): %.2fx\n",
              median(result$total_dose_mcg_kg_hr, na.rm = TRUE) /
                median(result$background_mcg_kg_hr, na.rm = TRUE)))
  cat(sprintf("Median bolus %% of total: %.1f%%\n\n",
              median(result$bolus_pct_of_total, na.rm = TRUE)))
  
  result
}


morphine_true_dose <- compute_true_total_dose(
  drug_name = "morphine",
  drug_regex = "morphine",
  flow_dose_col = "MORPHINE_NUMBER_OF_DOSES",
  conc_col = "morphine_conc_mg_per_ml"
)

oxycodone_true_dose <- compute_true_total_dose(
  drug_name = "oxycodone",
  drug_regex = "oxycodone",
  flow_dose_col = "OXYCODONE_NUMBER_OF_DOSES",
  conc_col = "oxycodone_conc_mg_per_ml"
)

# ============================================================
# COMBINE AND SAVE
# ============================================================

episode_true_total_dose <- bind_rows(morphine_true_dose, oxycodone_true_dose)

cat("\nTotal vs background-only dose by drug:\n")
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

cat("\nepisode_true_total_dose.rds saved\n")
