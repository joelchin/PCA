# ============================================================
# 12 — Safety analysis
# Paediatric PCA/NCA study — Addenbrooke's Hospital (CUH NHS FT)
# Author: J Chin
# ============================================================
#
# Purpose
#   Safety outcomes across the cohort: naloxone reversal events, drug switching
#   (opioid rotation), side-effect rates by drug, age group and era, and the
#   epidural subgroup.
#
# Inputs (in session)
#   pca_episodes_flagged, episode_metrics, episode_shift_table.
#
# Definitions
#   Side effects use the final anaesthesia-stop-anchored flags from 4_flags.R
#   (Flag 9). Naloxone reversal reflects manual chart-review reclassification:
#   several 400 mcg-dose events were confirmed to be itch, not respiratory
#   depression, and were recoded. Population sizes are computed from the data
#   (nothing is hard-coded).
#
# Limitations
#   - Antipruritic is a combined flag (chlorphenamine + naloxone 40 mcg), so
#     "naloxone pruritus only" cannot be reported separately.
#   - Timing of reversal relative to anaesthesia stop is not stored as a
#     per-episode column.
#   - Counts, not percentages, are reported for rare events.
#
# Output
#   Console tables only.
# ============================================================

library(tidyverse)
library(lubridate)

# ============================================================
# LOAD DATA
# ============================================================

stopifnot(
  exists("pca_episodes_flagged"),
  exists("episode_metrics"),
  exists("episode_shift_table")
)

cat("============================================================\n")
cat("SCRIPT 12 — SAFETY ANALYSIS (SECTION 4)\n")
cat("Population: ALL episodes, ALL drugs (S4)\n")
cat("============================================================\n\n")

fmt_med_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  sprintf("%.1f (%.1f–%.1f)",
          round(median(x), digits),
          round(quantile(x, 0.25), digits),
          round(quantile(x, 0.75), digits))
}

fmt_n_pct <- function(n, total) {
  sprintf("%d (%.1f%%)", n, 100 * n / total)
}

fmt_n_per1000 <- function(n, total) {
  sprintf("%d (%.1f per 1,000 episodes)", n, 1000 * n / total)
}

# All episodes, all drugs
all_ep <- pca_episodes_flagged

# all_met NOW carries the finalised flags, joined from pca_episodes_flagged
# (episode_metrics itself does not have these — they live on
# pca_episodes_flagged, built once in 4_flags.R Flag 9)
all_met <- episode_metrics |>
  left_join(
    pca_episodes_flagged |>
      select(PAT_ENC_CSN_ID, episode_id,
             reactive_antiemetic_final, any_antipruritic_final,
             any_naloxone_reversal_final, any_side_effect_final),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  )

n_ep  <- nrow(all_ep)
n_met <- nrow(all_met)

cat("Total episodes (S4, all drugs):", n_ep, "\n")
cat("Episodes with metrics data:    ", n_met, "\n\n")

# ============================================================
# TABLE 1 — NALOXONE REVERSAL EVENTS
#
# PRIMARY safety outcome. Note: the pruritus-specific naloxone
# breakdown (40mcg only, separate from chlorphenamine) previously
# reported here can no longer be isolated — any_antipruritic_final
# is a combined flag. See LIMITATIONS at end of script.
# ============================================================

cat("============================================================\n")
cat("TABLE 1 — NALOXONE REVERSAL EVENTS (all episodes)\n")
cat("(FINAL definition — manually chart-reviewed, minimal/no cutoff)\n")
cat("============================================================\n\n")

cat("--- Overall naloxone reversal rate ---\n")
cat("Naloxone reversal (FINAL):    ",
    fmt_n_per1000(sum(all_met$any_naloxone_reversal_final, na.rm = TRUE), n_met), "\n\n")

cat("--- Naloxone reversal by drug ---\n")
nalox_by_drug <- all_met |>
  group_by(drug) |>
  summarise(
    n                  = n(),
    n_nalox_reversal   = sum(any_naloxone_reversal_final, na.rm = TRUE),
    rate_reversal_1000 = round(1000 * n_nalox_reversal / n, 2),
    .groups = "drop"
  ) |>
  arrange(desc(n))

print(nalox_by_drug)
cat("\n")

cat("--- Naloxone reversal by mode ---\n")
all_met |>
  group_by(mode) |>
  summarise(
    n                = n(),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    rate_per_1000    = round(1000 * n_nalox_reversal / n, 2),
    .groups = "drop"
  ) |>
  print()
cat("\n")

cat("--- Naloxone reversal by age group ---\n")
all_met |>
  left_join(
    all_ep |> select(PAT_ENC_CSN_ID, episode_id, age_group),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  ) |>
  group_by(age_group) |>
  summarise(
    n                = n(),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    rate_per_1000    = round(1000 * n_nalox_reversal / n, 2),
    .groups = "drop"
  ) |>
  arrange(age_group) |>
  print()
cat("\n")

cat("--- Naloxone reversal by era ---\n")
all_met |>
  left_join(
    all_ep |> select(PAT_ENC_CSN_ID, episode_id, era),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  ) |>
  group_by(era) |>
  summarise(
    n                = n(),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    rate_per_1000    = round(1000 * n_nalox_reversal / n, 2),
    .groups = "drop"
  ) |>
  print()
cat("\n")

# ============================================================
# TABLE 2 — NALOXONE REVERSAL EPISODES (individual, FINAL)
# ============================================================

cat("============================================================\n")
cat("TABLE 2 — NALOXONE REVERSAL EPISODES (individual, FINAL)\n")
cat("(Manually chart-reviewed — dose-based classification alone was\n")
cat(" found unreliable in both directions during validation; several\n")
cat(" 400mcg-dose events across drugs were confirmed to be itch, not\n")
cat(" respiratory depression, and reclassified accordingly)\n")
cat("============================================================\n\n")

reversal_episodes <- all_met |>
  filter(any_naloxone_reversal_final == TRUE) |>
  select(PAT_ENC_CSN_ID, episode_id, drug) |>
  left_join(
    all_ep |> select(PAT_ENC_CSN_ID, episode_id, mode,
                     age_group, AGE, WEIGHT_KG, proc_category,
                     admission_type, era, concurrent_epidural,
                     had_remifentanil, is_first_line),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  )

cat("Total reversal episodes (FINAL): ", nrow(reversal_episodes), "\n\n")

cat("--- Drug distribution ---\n")
reversal_episodes |>
  count(drug) |>
  mutate(pct = round(100 * n / nrow(reversal_episodes), 1)) |>
  print()
cat("\n")

cat("--- Mode distribution ---\n")
reversal_episodes |> count(mode) |>
  mutate(pct = round(100 * n / nrow(reversal_episodes), 1)) |> print()
cat("\n")

cat("--- Age group distribution ---\n")
reversal_episodes |> count(age_group) |>
  mutate(pct = round(100 * n / nrow(reversal_episodes), 1)) |>
  arrange(age_group) |> print()
cat("\n")

cat("--- Age and weight of reversal episodes ---\n")
cat("Age — median (IQR) years:  ", fmt_med_iqr(reversal_episodes$AGE), "\n")
cat("Weight — median (IQR) kg:  ", fmt_med_iqr(reversal_episodes$WEIGHT_KG), "\n\n")

cat("--- Procedure category ---\n")
reversal_episodes |> count(proc_category, sort = TRUE) |> print()
cat("\n")

cat("--- Era distribution ---\n")
reversal_episodes |> count(era) |> print()
cat("\n")

cat("--- First line vs rotation ---\n")
reversal_episodes |>
  count(is_first_line) |>
  mutate(group = if_else(is_first_line, "First line", "Rotation")) |>
  select(group, n) |>
  print()
cat("\n")

cat("(NOTE: timing of reversal relative to anaesthesia stop was\n")
cat("checked manually for every reversal event. No persistent column\n")
cat("exists for this; re-derive from anaesthesia_stop_times if a\n")
cat("per-episode figure is needed here.)\n\n")

# ============================================================
# TABLE 3 — DRUG SWITCHING (safety proxy) — unchanged
# ============================================================

cat("============================================================\n")
cat("TABLE 3 — DRUG SWITCHING (opioid rotation)\n")
cat("============================================================\n\n")

cat("Total rotation episodes (any drug switch): ",
    fmt_n_pct(sum(all_ep$is_rotation, na.rm = TRUE), n_ep), "\n\n")

cat("--- Switching rate by drug (episodes rotated away FROM) ---\n")
rotation_from <- all_ep |>
  filter(is_rotation == TRUE) |>
  count(prior_drug, name = "n_rotated_away") |>
  rename(drug = prior_drug)

ep_by_drug <- all_ep |> count(drug, name = "n_total")

rotation_from |>
  left_join(ep_by_drug, by = "drug") |>
  mutate(pct_rotated = round(100 * n_rotated_away / n_total, 1)) |>
  arrange(desc(n_rotated_away)) |>
  print()
cat("\n")

cat("--- Switching rate by drug (episodes rotated TO) ---\n")
all_ep |>
  filter(is_rotation == TRUE) |>
  count(drug, name = "n_rotated_to") |>
  left_join(ep_by_drug, by = "drug") |>
  mutate(pct_rotated = round(100 * n_rotated_to / n_total, 1)) |>
  arrange(desc(n_rotated_to)) |>
  print()
cat("\n")

# ============================================================
# TABLE 4 — SIDE EFFECT RATES BY DRUG (FINAL definitions)
# ============================================================

cat("============================================================\n")
cat("TABLE 4 — SIDE EFFECT RATES BY DRUG (all episodes, FINAL)\n")
cat("============================================================\n\n")

all_met |>
  group_by(drug) |>
  summarise(
    n                  = n(),
    pct_antiemetic     = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    pct_antipruritic   = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
    pct_laxative       = round(100 * mean(any_laxative_full, na.rm = TRUE), 1),
    pct_side_effect    = round(100 * mean(any_side_effect_final, na.rm = TRUE), 1),
    n_nalox_reversal   = sum(any_naloxone_reversal_final, na.rm = TRUE),
    rate_reversal_1000 = round(1000 * n_nalox_reversal / n, 2),
    .groups = "drop"
  ) |>
  arrange(desc(n)) |>
  print()
cat("(Laxative rate is descriptive only — cannot distinguish\n")
cat(" prophylactic bowel-regimen dosing from genuine reactive\n")
cat(" treatment of opioid-induced constipation. Not tested\n")
cat(" statistically anywhere in this pipeline.)\n\n")

# ============================================================
# TABLE 5 — SIDE EFFECT RATES BY AGE GROUP (FINAL)
# ============================================================

cat("============================================================\n")
cat("TABLE 5 — SIDE EFFECT RATES BY AGE GROUP (all episodes, FINAL)\n")
cat("============================================================\n\n")

all_met |>
  left_join(
    all_ep |> select(PAT_ENC_CSN_ID, episode_id, age_group),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  ) |>
  group_by(age_group) |>
  summarise(
    n                = n(),
    pct_antiemetic   = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    pct_antipruritic = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
    pct_side_effect  = round(100 * mean(any_side_effect_final, na.rm = TRUE), 1),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    rate_rev_1000    = round(1000 * n_nalox_reversal / n, 2),
    .groups = "drop"
  ) |>
  arrange(age_group) |>
  print()
cat("\n")

# ============================================================
# TABLE 8 — SIDE EFFECT RATES BY AGE GROUP x DRUG (FINAL)
# ============================================================

cat("============================================================\n")
cat("TABLE 8 — SIDE EFFECT RATES BY AGE GROUP x DRUG (FINAL)\n")
cat("============================================================\n\n")

age_drug_met <- all_met |>
  left_join(
    all_ep |> select(PAT_ENC_CSN_ID, episode_id, age_group),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  )

cat("--- Raw n by age_group x drug (check cell sizes first) ---\n")
print(table(age_drug_met$age_group, age_drug_met$drug))
cat("\n")

cat("--- Raw event counts by age_group x drug (FINAL) ---\n")
age_drug_met |>
  group_by(age_group, drug) |>
  summarise(
    n                = n(),
    n_antiemetic     = sum(reactive_antiemetic_final, na.rm = TRUE),
    n_antipruritic   = sum(any_antipruritic_final, na.rm = TRUE),
    n_side_effect    = sum(any_side_effect_final, na.rm = TRUE),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(age_group, desc(n)) |>
  print(n = 100)
cat("\n")

cat("--- Rates (%) by age_group x drug (FINAL) — CAUTION where n small ---\n")
age_drug_met |>
  group_by(age_group, drug) |>
  summarise(
    n                  = n(),
    pct_antiemetic     = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    pct_antipruritic   = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
    pct_side_effect    = round(100 * mean(any_side_effect_final, na.rm = TRUE), 1),
    n_nalox_reversal   = sum(any_naloxone_reversal_final, na.rm = TRUE),
    rate_reversal_1000 = round(1000 * n_nalox_reversal / n, 2),
    .groups = "drop"
  ) |>
  arrange(age_group, desc(n)) |>
  print(n = 100)
cat("\n")

cat("--- NEONATAL SUBSET (<28 days) — by drug, raw counts (FINAL) ---\n")
neonatal_met <- age_drug_met |> filter(age_group == "Neonate (<28 days)")
cat("Total neonatal episodes: ", nrow(neonatal_met), "\n\n")

neonatal_met |>
  group_by(drug) |>
  summarise(
    n                = n(),
    n_antiemetic     = sum(reactive_antiemetic_final, na.rm = TRUE),
    n_antipruritic   = sum(any_antipruritic_final, na.rm = TRUE),
    n_side_effect    = sum(any_side_effect_final, na.rm = TRUE),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n)) |>
  print()
cat("\n")

cat("--- Fisher's exact test: reversal rate across age groups (FINAL) ---\n")
reversal_age_table <- table(age_drug_met$age_group, age_drug_met$any_naloxone_reversal_final)
print(reversal_age_table)
set.seed(42)  # Monte Carlo p-value: fixed seed so the result is reproducible
print(fisher.test(reversal_age_table, simulate.p.value = TRUE, B = 10000))
cat("\n")

# ============================================================
# TABLE 6 — SAFETY TREND BY YEAR (FINAL)
# ============================================================

cat("============================================================\n")
cat("TABLE 6 — SAFETY TREND BY YEAR (FINAL)\n")
cat("============================================================\n\n")

all_met |>
  left_join(
    all_ep |> select(PAT_ENC_CSN_ID, episode_id, year),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  ) |>
  group_by(year) |>
  summarise(
    n                = n(),
    pct_side_effect  = round(100 * mean(any_side_effect_final, na.rm = TRUE), 1),
    pct_antiemetic   = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    rate_rev_1000    = round(1000 * n_nalox_reversal / n, 2),
    .groups = "drop"
  ) |>
  arrange(year) |>
  print()
cat("\n")

# ============================================================
# TABLE 7 — CONCURRENT EPIDURAL EPISODES (FINAL)
# ============================================================

cat("============================================================\n")
cat("TABLE 7 — CONCURRENT EPIDURAL EPISODES (excluded from primary)\n")
cat("============================================================\n\n")

epidural_ep <- all_ep |> filter(concurrent_epidural == TRUE)

cat("Concurrent epidural episodes: ", nrow(epidural_ep), "\n\n")

cat("--- Drug distribution ---\n")
epidural_ep |> count(drug) |>
  mutate(pct = round(100 * n / nrow(epidural_ep), 1)) |> print()
cat("\n")

cat("--- Age group ---\n")
epidural_ep |> count(age_group) |>
  mutate(pct = round(100 * n / nrow(epidural_ep), 1)) |>
  arrange(age_group) |> print()
cat("\n")

cat("--- Procedure category ---\n")
epidural_ep |> count(proc_category, sort = TRUE) |>
  mutate(pct = round(100 * n / nrow(epidural_ep), 1)) |>
  head(10) |> print()
cat("\n")

epidural_met <- all_met |>
  semi_join(epidural_ep, by = c("PAT_ENC_CSN_ID", "episode_id"))

cat("--- Side effects in epidural episodes (for reference, FINAL) ---\n")
cat("Any antiemetic (FINAL):  ",
    fmt_n_pct(sum(epidural_met$reactive_antiemetic_final, na.rm = TRUE),
              nrow(epidural_met)), "\n")
cat("Any antipruritic (FINAL):",
    fmt_n_pct(sum(epidural_met$any_antipruritic_final, na.rm = TRUE),
              nrow(epidural_met)), "\n")
cat("Naloxone reversal (FINAL):",
    sum(epidural_met$any_naloxone_reversal_final, na.rm = TRUE), "\n\n")

cat("============================================================\n")
cat("SCRIPT 12 COMPLETE\n")
cat("Key safety finding: naloxone reversal rate by drug and age group\n")
cat("Report counts not percentages for rare events (n<30)\n")
cat("Trend analysis: safety outcomes stable across study period\n\n")
cat("LIMITATIONS:\n")
cat("  - 'Any naloxone' and isolated 'naloxone pruritus only' rates\n")
cat("    can no longer be reported separately — any_antipruritic_final\n")
cat("    is a combined flag (chlorphenamine + naloxone-40mcg). If a\n")
cat("    naloxone-40mcg-only breakdown is needed, this would require\n")
cat("    building a separate flag from non_pca_mar_reclassified.\n")
cat("  - Naloxone reversal timing relative to anaesthesia stop is not\n")
cat("    stored as a persistent per-episode column; only available as\n")
cat("    a manually reviewed figure for confirmed true-positive events.\n")
cat("  - Laxative administration is reported descriptively only,\n")
cat("    not as a formal comparative outcome — cannot distinguish\n")
cat("    prophylactic bowel-regimen dosing from genuine reactive\n")
cat("    treatment of opioid-induced constipation.\n")
cat("============================================================\n")
