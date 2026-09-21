# ============================================================
# 3 — Episode-by-shift table and episode-level metrics
# Paediatric PCA/NCA study — Addenbrooke's Hospital (CUH NHS FT)
# Author: J Chin
# ============================================================
#
# Purpose
#   Expands each PCA/NCA episode into one row per 12-hour shift and derives
#   episode-level summaries (pain, analgesic and side-effect markers by shift,
#   demographics).
#
# Inputs (files, 1 - data/3 - processed_data/)
#   admissions_clean, pca_episodes, non_pca_mar_clean, pain_scores_clean,
#   operations_clean, anaesthesia_events_clean, 1_shifts.
#
# Outputs (.rds)
#   episode_shift_table  one row per shift per active episode
#                        (temporal analyses, pain sequencing)
#   episode_metrics      one row per episode (derived summaries)
#   Also reads drug_lists.R, only to check its antiemetic list against the one used here.
#
# Notes
#   - A shift is a clock-anchored 12-hour block, used for temporal sequencing
#     only (not day/night classification). shift_in_episode counts shifts from
#     the episode start; shift_index_from_pca_start (after D. Brooks) is 0 for
#     the shift in which PCA started, negative before, positive after.
#   - The shift-based side-effect fields built here (e.g. shift 1 = prophylaxis
#     at induction; shift 2 onwards) are kept for descriptive and temporal
#     analyses. The side-effect OUTCOMES used in the comparative analysis are
#     the anaesthesia-stop-anchored flags built in 4_flags.R (Flag 9).
#   - Runtime: about 10 minutes.
# ============================================================

library(tidyverse)
library(lubridate)

# Antiemetics counted in this script's shift-level descriptive fields. The outcome
# flags (4_flags.R, Flag 9) use antiemetic_drugs from drug_lists.R; any difference
# between the two lists is reported here rather than left to drift silently.
shift_table_antiemetics <- c("ondansetron", "cyclizine", "metoclopramide",
                             "droperidol", "haloperidol")
source("drug_lists.R")
if (exists("antiemetic_drugs")) {
  list_diff <- c(setdiff(shift_table_antiemetics, tolower(antiemetic_drugs)),
                 setdiff(tolower(antiemetic_drugs), shift_table_antiemetics))
  if (length(list_diff) > 0) {
    warning("shift-table antiemetic list differs from drug_lists.R antiemetic_drugs: ",
            paste(list_diff, collapse = ", "))
  }
}

# ============================================================
# LOAD CLEANED DATA
# ============================================================

admissions_clean         <- readRDS("1 - data/3 - processed_data/admissions_clean.rds")
pca_episodes             <- readRDS("1 - data/3 - processed_data/pca_episodes.rds")
non_pca_mar_clean        <- readRDS("1 - data/3 - processed_data/non_pca_mar_clean.rds")
pain_scores_clean        <- readRDS("1 - data/3 - processed_data/pain_scores_clean.rds")
operations_clean         <- readRDS("1 - data/3 - processed_data/operations_clean.rds")
anaesthesia_events_clean <- readRDS("1 - data/3 - processed_data/anaesthesia_events_clean.rds")
shifts                   <- readRDS("1 - data/3 - processed_data/1_shifts.rds")

cat("All cleaned data loaded\n")

# ============================================================
# STEP 1 — EXPAND EPISODES TO ONE ROW PER ACTIVE SHIFT
#
# shift_in_episode = temporal index within episode
#   1 = first shift, 2 = second shift etc.
#
# shift_index_from_pca_start = index anchored to PCA start
#   0 = shift PCA started
#   negative = pre-PCA shifts (not used in main analysis)
#   positive = post-PCA start shifts
#
# Note: shift is the unit of clinical time — 12 hours
# Used for temporal sequencing NOT day/night classification
# ============================================================

cat("\nBuilding episode-shift expansion...\n")

episode_shifts <- pca_episodes |>
  filter(
    !is.na(episode_start_shift),
    !is.na(episode_end_shift)
  ) |>
  select(
    PAT_ENC_CSN_ID, episode_id, drug, mode,
    episode_start, episode_end, episode_duration_hrs,
    episode_start_shift, episode_end_shift,
    bolus_dose_ml, lockout_mins,
    n_orders, n_syringes, has_missing_stop, drug_switch
  ) |>
  rowwise() |>
  mutate(
    shift_number = list(
      seq(episode_start_shift, episode_end_shift, by = 1)
    )
  ) |>
  unnest(shift_number) |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  mutate(
    # Position within episode (1 = first shift)
    shift_in_episode = row_number(),
    # Index anchored to PCA start
    # 0 = shift PCA started
    shift_index_from_pca_start = shift_number - episode_start_shift,
    shift_phase = case_when(
      shift_index_from_pca_start == 0 ~ "PCA_start_shift",
      shift_index_from_pca_start > 0  ~ "post_PCA_start",
      TRUE                            ~ "pre_PCA"
    ),
    # Flag first shift — may reflect prophylactic prescribing
    # rather than treatment of established opioid side effects
    is_first_shift = shift_in_episode == 1
  ) |>
  ungroup()

cat("Episode-shift rows:", nrow(episode_shifts), "\n")
cat("Unique episodes:", n_distinct(
  paste(episode_shifts$PAT_ENC_CSN_ID,
        episode_shifts$episode_id)), "\n")
cat("Median shifts per episode:",
    median(table(paste(episode_shifts$PAT_ENC_CSN_ID,
                       episode_shifts$episode_id))), "\n")

# ============================================================
# STEP 2 — JOIN SHIFT METADATA AND ANAESTHESIA STOP TIME
# Anaesthesia stop retained for descriptive context only
# Not used as primary analysis window filter
# ============================================================

cat("\nJoining shift metadata and anaesthesia stop times...\n")

anaesthesia_stop_times <- anaesthesia_events_clean |>
  group_by(PAT_ENC_CSN_ID) |>
  summarise(
    anaesthesia_stop          = min(anaesthesia_stop, na.rm = TRUE),
    anaesthesia_duration_mins = first(anaesthesia_duration_mins),
    .groups = "drop"
  )

episode_shifts <- episode_shifts |>
  left_join(
    shifts |> select(shift_number, start_time, end_time),
    by = "shift_number"
  ) |>
  rename(
    shift_start = start_time,
    shift_end   = end_time
  ) |>
  left_join(anaesthesia_stop_times, by = "PAT_ENC_CSN_ID") |>
  mutate(
    # Hours from anaesthesia stop to start of this shift
    # Retained for descriptive/sensitivity analyses
    hrs_from_anaesthesia_stop = as.numeric(
      difftime(shift_start, anaesthesia_stop, units = "hours")
    )
  )

cat("Shift metadata and anaesthesia stop times joined\n")

# ============================================================
# STEP 3 — JOIN PAIN SCORES BY SHIFT
#
# Pain scores linked by shift number — temporal attribution
# without requiring exact timestamp matching
#
# Data quality note:
# Pain scores show 4-hourly temporal clustering (6am, 10am,
# 2pm, 6pm, 10pm) confirming real-time documentation.
# This was confirmed from raw timestamp distribution
# independently of the shift system.
# ============================================================

cat("\nJoining pain scores...\n")

pain_by_shift <- pain_scores_clean |>
  filter(!is.na(pain_score_shift)) |>
  group_by(PAT_ENC_CSN_ID, pain_score_shift) |>
  summarise(
    pain_max          = max(score_value, na.rm = TRUE),
    pain_median       = median(score_value, na.rm = TRUE),
    pain_n_scores     = n(),
    had_severe_pain   = any(score_value >= 7, na.rm = TRUE),
    had_moderate_pain = any(score_value >= 4, na.rm = TRUE),
    had_mild_pain     = any(
      score_value >= 1 & score_value <= 3, na.rm = TRUE
    ),
    had_any_pain      = any(score_value > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  rename(shift_number = pain_score_shift)

episode_shifts <- episode_shifts |>
  left_join(pain_by_shift,
            by = c("PAT_ENC_CSN_ID", "shift_number")) |>
  mutate(
    no_pain_assessment = is.na(pain_n_scores),
    had_severe_pain    = replace_na(had_severe_pain, FALSE),
    had_moderate_pain  = replace_na(had_moderate_pain, FALSE),
    had_mild_pain      = replace_na(had_mild_pain, FALSE),
    had_any_pain       = replace_na(had_any_pain, FALSE)
  )

cat("Pain scores joined\n")
cat("Shifts with pain assessment:",
    sum(!episode_shifts$no_pain_assessment), "\n")
cat("Shifts without pain assessment:",
    sum(episode_shifts$no_pain_assessment), "\n")

# ============================================================
# STEP 4 — JOIN SIDE EFFECT MARKERS BY SHIFT
#
# Attribution:
# Drug given=TRUE during active episode shift only
# inner_join ensures only episode-active shifts included
# Preop, intraop, post-episode drugs automatically excluded
#
# Naloxone distinguished by formulation:
# - low_dose_pruritus (40mcg) = pruritus proxy
# - reversal_respiratory (400mcg) = key safety outcome
#
# Dexamethasone excluded:
# - Repeated scheduled dosing (every 6-8hrs over days)
# - Predominantly intraop prophylaxis or anti-inflammatory
# - Inconsistent with PONV attribution
# ============================================================

cat("\nJoining side effect markers...\n")

side_effects_by_shift <- non_pca_mar_clean |>
  filter(
    drug_category == "side_effect_marker",
    given == TRUE,
    !is.na(med_admin_shift)
  ) |>
  # inner_join restricts to episode-active shifts only
  # automatically excludes preop/intraop/post-episode drugs
  inner_join(
    episode_shifts |>
      select(PAT_ENC_CSN_ID, shift_number) |>
      distinct(),
    by = c("PAT_ENC_CSN_ID",
           "med_admin_shift" = "shift_number")
  ) |>
  group_by(PAT_ENC_CSN_ID, med_admin_shift) |>
  summarise(
    antiemetic_given   = any(tolower(GENERIC_NAME) %in% shift_table_antiemetics),
    antipruritic_given = any(
      tolower(GENERIC_NAME) == "chlorphenamine" |
        naloxone_type == "low_dose_pruritus"
    ),
    naloxone_any       = any(
      tolower(GENERIC_NAME) == "naloxone"
    ),
    naloxone_reversal  = any(
      naloxone_type == "reversal_respiratory",
      na.rm = TRUE
    ),
    naloxone_pruritus  = any(
      naloxone_type == "low_dose_pruritus",
      na.rm = TRUE
    ),
    any_side_effect    = TRUE,
    .groups = "drop"
  ) |>
  rename(shift_number = med_admin_shift)

episode_shifts <- episode_shifts |>
  left_join(side_effects_by_shift,
            by = c("PAT_ENC_CSN_ID", "shift_number")) |>
  mutate(
    antiemetic_given   = replace_na(antiemetic_given, FALSE),
    antipruritic_given = replace_na(antipruritic_given, FALSE),
    naloxone_any       = replace_na(naloxone_any, FALSE),
    naloxone_reversal  = replace_na(naloxone_reversal, FALSE),
    naloxone_pruritus  = replace_na(naloxone_pruritus, FALSE),
    any_side_effect    = replace_na(any_side_effect, FALSE)
  )

cat("Side effect markers joined\n")
cat("Shifts with antiemetic:",
    sum(episode_shifts$antiemetic_given), "\n")
cat("  of which shift 1:",
    sum(episode_shifts$antiemetic_given &
          episode_shifts$is_first_shift), "\n")
cat("  of which shift 2+:",
    sum(episode_shifts$antiemetic_given &
          !episode_shifts$is_first_shift), "\n")
cat("Shifts with antipruritic:",
    sum(episode_shifts$antipruritic_given), "\n")
cat("Shifts with naloxone reversal:",
    sum(episode_shifts$naloxone_reversal), "\n")

# ============================================================
# STEP 5 — JOIN LAXATIVES BY SHIFT
#
# Opioid-induced constipation proxy
# Expected temporal pattern: increases with shift_in_episode
# Early = pre-existing or non-opioid
# Late = more likely opioid-induced
# ============================================================

cat("\nJoining laxatives...\n")

laxatives_by_shift <- non_pca_mar_clean |>
  filter(
    drug_category == "laxative",
    given == TRUE,
    !is.na(med_admin_shift)
  ) |>
  inner_join(
    episode_shifts |>
      select(PAT_ENC_CSN_ID, shift_number) |>
      distinct(),
    by = c("PAT_ENC_CSN_ID",
           "med_admin_shift" = "shift_number")
  ) |>
  group_by(PAT_ENC_CSN_ID, med_admin_shift) |>
  summarise(
    laxative_given = TRUE,
    .groups = "drop"
  ) |>
  rename(shift_number = med_admin_shift)

episode_shifts <- episode_shifts |>
  left_join(laxatives_by_shift,
            by = c("PAT_ENC_CSN_ID", "shift_number")) |>
  mutate(
    laxative_given = replace_na(laxative_given, FALSE)
  )

cat("Laxatives joined\n")

# ============================================================
# STEP 6 — JOIN ANALGESICS BY SHIFT
# ============================================================

cat("\nJoining analgesics...\n")

analgesics_by_shift <- non_pca_mar_clean |>
  filter(
    drug_category == "analgesic",
    given == TRUE,
    !is.na(med_admin_shift)
  ) |>
  inner_join(
    episode_shifts |>
      select(PAT_ENC_CSN_ID, shift_number) |>
      distinct(),
    by = c("PAT_ENC_CSN_ID",
           "med_admin_shift" = "shift_number")
  ) |>
  group_by(PAT_ENC_CSN_ID, med_admin_shift) |>
  summarise(
    paracetamol_given = any(
      tolower(GENERIC_NAME) == "paracetamol"
    ),
    nsaid_given       = any(analgesic_class == "NSAID"),
    opioid_given      = any(analgesic_class == "opioid"),
    adjuvant_given    = any(analgesic_class == "adjuvant"),
    neuropathic_given = any(analgesic_class == "neuropathic"),
    .groups = "drop"
  ) |>
  rename(shift_number = med_admin_shift)

episode_shifts <- episode_shifts |>
  left_join(analgesics_by_shift,
            by = c("PAT_ENC_CSN_ID", "shift_number")) |>
  mutate(
    paracetamol_given = replace_na(paracetamol_given, FALSE),
    nsaid_given       = replace_na(nsaid_given, FALSE),
    opioid_given      = replace_na(opioid_given, FALSE),
    adjuvant_given    = replace_na(adjuvant_given, FALSE),
    neuropathic_given = replace_na(neuropathic_given, FALSE)
  )

cat("Analgesics joined\n")

# ============================================================
# STEP 7 — JOIN DEMOGRAPHICS
# ============================================================

cat("\nJoining demographics...\n")

episode_shifts <- episode_shifts |>
  left_join(
    admissions_clean |>
      select(PAT_ENC_CSN_ID, AGE, age_days, age_group,
             WEIGHT_KG, GENDER, ADMITTING_SPECIALTY,
             case_type),
    by = "PAT_ENC_CSN_ID"
  ) |>
  left_join(
    operations_clean |>
      filter(primary_procedure == TRUE) |>
      select(PAT_ENC_CSN_ID, proc_category) |>
      distinct(PAT_ENC_CSN_ID, .keep_all = TRUE),
    by = "PAT_ENC_CSN_ID"
  )

cat("Demographics joined\n")

# ============================================================
# STEP 8 — DERIVE EPISODE-LEVEL METRICS
#
# Three side effect windows:
#
# Window A — Full episode (shift 1 onwards)
#   Broadest — includes prophylactic shift 1 administrations
#   Use for: sensitivity analysis, population burden
#
# Window B — Post-shift-1 (shift 2 onwards) ← PRIMARY
#   Excludes shift 1 prophylactic prescribing at induction
#   More likely to represent established opioid side effects
#   Use for: primary outcome, main comparison
#
# Window C — Shift 1 only
#   Captures prophylactic effect
#   Use for: demonstrating prophylaxis inflation
#   Difference between A and B shows prophylaxis contribution
#
# Three metric types per window:
# 1. Binary: any administration during window
# 2. Time to first: shift_in_episode at first administration
# 3. Burden: proportion of shifts with administration
# ============================================================

cat("\nDeriving episode metrics...\n")

episode_metrics <- episode_shifts |>
  group_by(PAT_ENC_CSN_ID, episode_id, drug, mode) |>
  summarise(

    # ---- Episode duration ----
    total_shifts         = n(),
    episode_duration_hrs = first(episode_duration_hrs),

    # ---- Pain metrics ----
    shifts_with_assessment    = sum(!no_pain_assessment),
    shifts_no_assessment      = sum(no_pain_assessment),
    pct_shifts_no_assessment  = mean(no_pain_assessment) * 100,
    shifts_severe_pain        = sum(had_severe_pain),
    shifts_moderate_pain      = sum(had_moderate_pain),
    pct_shifts_severe         = mean(had_severe_pain) * 100,
    pct_shifts_moderate       = mean(had_moderate_pain) * 100,

    # Time to pain control
    # First shift where pain mild/absent AND assessment done
    first_controlled_shift = min(
      shift_in_episode[
        !had_moderate_pain &
          !had_severe_pain &
          !no_pain_assessment
      ],
      na.rm = TRUE
    ),

    # ---- WINDOW A: Full episode (shift 1 onwards) ----
    # Includes prophylactic shift 1 — use for sensitivity only

    any_antiemetic_full      = any(antiemetic_given),
    any_antipruritic_full    = any(antipruritic_given),
    any_naloxone_full        = any(naloxone_any),
    any_naloxone_reversal    = any(naloxone_reversal),
    any_naloxone_pruritus    = any(naloxone_pruritus),
    any_side_effect_full     = any(any_side_effect),
    any_laxative_full        = any(laxative_given),

    pct_shifts_antiemetic_full    = mean(antiemetic_given) * 100,
    pct_shifts_antipruritic_full  = mean(antipruritic_given) * 100,
    pct_shifts_side_effect_full   = mean(any_side_effect) * 100,

    shift_first_antiemetic_full   = min(
      shift_in_episode[antiemetic_given], na.rm = TRUE
    ),
    shift_first_antipruritic_full = min(
      shift_in_episode[antipruritic_given], na.rm = TRUE
    ),
    shift_first_naloxone_full     = min(
      shift_in_episode[naloxone_any], na.rm = TRUE
    ),

    # ---- WINDOW B: Post-shift-1 (shift 2 onwards) ← PRIMARY ----
    # Excludes shift 1 prophylactic prescribing at induction
    # More likely to represent established opioid side effects

    any_antiemetic           = any(
      antiemetic_given[shift_in_episode > 1],
      na.rm = TRUE
    ),
    any_antipruritic         = any(
      antipruritic_given[shift_in_episode > 1],
      na.rm = TRUE
    ),
    any_side_effect          = any(
      any_side_effect[shift_in_episode > 1],
      na.rm = TRUE
    ),
    any_laxative             = any(
      laxative_given[shift_in_episode > 1],
      na.rm = TRUE
    ),

    pct_shifts_antiemetic    = mean(
      antiemetic_given[shift_in_episode > 1],
      na.rm = TRUE
    ) * 100,
    pct_shifts_antipruritic  = mean(
      antipruritic_given[shift_in_episode > 1],
      na.rm = TRUE
    ) * 100,
    pct_shifts_side_effect   = mean(
      any_side_effect[shift_in_episode > 1],
      na.rm = TRUE
    ) * 100,

    shift_first_antiemetic   = min(
      shift_in_episode[
        antiemetic_given & shift_in_episode > 1
      ],
      na.rm = TRUE
    ),
    shift_first_antipruritic = min(
      shift_in_episode[
        antipruritic_given & shift_in_episode > 1
      ],
      na.rm = TRUE
    ),
    shift_first_naloxone     = min(
      shift_in_episode[
        naloxone_any & shift_in_episode > 1
      ],
      na.rm = TRUE
    ),
    shift_first_laxative     = min(
      shift_in_episode[
        laxative_given & shift_in_episode > 1
      ],
      na.rm = TRUE
    ),

    # ---- WINDOW C: Shift 1 only ----
    # Captures prophylactic effect at induction
    # Difference between Window A and B = prophylaxis contribution

    any_antiemetic_shift1    = any(
      antiemetic_given[shift_in_episode == 1],
      na.rm = TRUE
    ),
    any_antipruritic_shift1  = any(
      antipruritic_given[shift_in_episode == 1],
      na.rm = TRUE
    ),
    any_side_effect_shift1   = any(
      any_side_effect[shift_in_episode == 1],
      na.rm = TRUE
    ),

    # ---- Co-analgesia metrics ----
    pct_shifts_paracetamol   = mean(paracetamol_given) * 100,
    pct_shifts_nsaid         = mean(nsaid_given) * 100,
    pct_shifts_adjuvant      = mean(adjuvant_given) * 100,
    pct_shifts_neuropathic   = mean(neuropathic_given) * 100,
    pct_shifts_laxative      = mean(laxative_given) * 100,

    .groups = "drop"

  ) |>

  # Clean up Inf values where events never occurred
  mutate(
    across(
      c(first_controlled_shift,
        shift_first_antiemetic_full,
        shift_first_antipruritic_full,
        shift_first_naloxone_full,
        shift_first_antiemetic,
        shift_first_antipruritic,
        shift_first_naloxone,
        shift_first_laxative),
      ~ if_else(is.infinite(.x), NA_real_, .x)
    )
  )

cat("Episode metrics derived\n")
cat("Total episodes with metrics:", nrow(episode_metrics), "\n")

# ---- Sanity check — compare windows ----
cat("\nSide effect rates by drug — Window A (full episode):\n")
episode_metrics |>
  group_by(drug) |>
  summarise(
    n                        = n(),
    pct_antiemetic_full      = round(
      mean(any_antiemetic_full) * 100, 1
    ),
    pct_antipruritic_full    = round(
      mean(any_antipruritic_full) * 100, 1
    ),
    pct_reversal             = round(
      mean(any_naloxone_reversal) * 100, 2
    )
  ) |>
  print()

cat("\nSide effect rates by drug — Window B (shift 2+ PRIMARY):\n")
episode_metrics |>
  group_by(drug) |>
  summarise(
    n                        = n(),
    pct_antiemetic           = round(
      mean(any_antiemetic, na.rm = TRUE) * 100, 1
    ),
    pct_antipruritic         = round(
      mean(any_antipruritic, na.rm = TRUE) * 100, 1
    )
  ) |>
  print()

cat("\nShift 1 prophylaxis rates by drug — Window C:\n")
episode_metrics |>
  group_by(drug) |>
  summarise(
    n                        = n(),
    pct_antiemetic_shift1    = round(
      mean(any_antiemetic_shift1, na.rm = TRUE) * 100, 1
    ),
    pct_antipruritic_shift1  = round(
      mean(any_antipruritic_shift1, na.rm = TRUE) * 100, 1
    )
  ) |>
  print()

# ============================================================
# STEP 9 — SAVE BOTH TABLES
# ============================================================

dir.create("1 - data/3 - processed_data",
           showWarnings = FALSE, recursive = TRUE)

saveRDS(episode_shifts,
        "1 - data/3 - processed_data/episode_shift_table.rds")

saveRDS(episode_metrics,
        "1 - data/3 - processed_data/episode_metrics.rds")

cat("\nBoth tables saved:\n")
cat("  episode_shift_table.rds —",
    nrow(episode_shifts), "rows\n")
cat("  episode_metrics.rds —",
    nrow(episode_metrics), "rows\n")
cat("\nFiles in processed_data:\n")
list.files("1 - data/3 - processed_data/")
