# ============================================================
# 11 — Section 2: Oxycodone Descriptive
# Paediatric PCA Study
# Addenbrooke's Hospital
#
# Purpose:
# Describes oxycodone PCA/NCA episodes in detail
# Primary novel contribution of the paper
#
# Population:
#   Table 1 (Demographics): unique oxycodone patients (PAT_ID dedup)
#   All other tables: S4 — all oxycodone episodes (n=1,237)
#
# Side-effect definitions (see 4_flags.R Flag 9):
#   reactive_antiemetic_final, any_antipruritic_final,
#   any_naloxone_reversal_final, any_side_effect_final
#
# Pain/co-analgesia columns (pct_shifts_severe, first_controlled_
# shift, pct_shifts_paracetamol, etc.) are joined onto s4_oxy once at
# build time, so Tables 5, 6, 7, and 8 all read from the same base
# object with these columns present.
#
# Table 1's case type line shows both the original, unrefined
# case_type ("does any OPCS-coded operation exist for this
# admission") and case_type_refined side by side — the original for
# an audit trail, the refined version for actual reporting.
# ============================================================

library(tidyverse)
library(lubridate)

# ============================================================
# LOAD DATA
# ============================================================

stopifnot(
  exists("pca_episodes_flagged"),
  exists("episode_metrics"),
  exists("episode_shift_table"),
  exists("admissions_clean")
)

stopifnot("case_type_refined" %in% names(pca_episodes_flagged))

cat("============================================================\n")
cat("ANALYSIS 5 — OXYCODONE DESCRIPTIVE (SECTION 2)\n")
cat("============================================================\n\n")

# ============================================================
# LOAD CONSUMPTION DATA
# ============================================================

cat("Calculating oxycodone consumption (mcg/kg/hr) directly from episode_dose...\n\n")

episode_dose_for_consumption <- tryCatch(
  readRDS("1 - data/3 - processed_data/episode_dose.rds"),
  error = function(e) {
    cat("episode_dose.rds not found — run 7_dosing.R first\n")
    return(NULL)
  }
)

if (!is.null(episode_dose_for_consumption)) {
  MIN_DURATION_HRS <- 2
  daily_consumption_prefloor <- episode_dose_for_consumption |>
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
} else {
  daily_consumption <- NULL
  cat("Consumption section will be skipped — episode_dose.rds not available\n\n")
}

# ============================================================
# DEFINE POPULATIONS — S4 (all oxycodone episodes)
# BUGFIX: pain/co-analgesia columns joined here, ONCE, so every
# downstream table can reference them via plain s4_oxy$column.
# ============================================================

s4_oxy <- pca_episodes_flagged |>
  filter(drug == "oxycodone") |>
  semi_join(episode_metrics |> filter(drug == "oxycodone"),
            by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  left_join(
    episode_metrics |>
      select(PAT_ENC_CSN_ID, episode_id,
             pct_shifts_severe, pct_shifts_moderate, pct_shifts_no_assessment,
             first_controlled_shift, pct_shifts_paracetamol,
             pct_shifts_nsaid, pct_shifts_adjuvant, pct_shifts_neuropathic),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  )

# ---- Table 1 population: TRUE PATIENT-LEVEL dedup ----
oxycodone_patients <- s4_oxy |>
  left_join(admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID),
            by = "PAT_ENC_CSN_ID") |>
  arrange(PAT_ID, episode_start) |>
  distinct(PAT_ID, .keep_all = TRUE)

n_patients <- nrow(oxycodone_patients)
n_episodes <- nrow(s4_oxy)

cat("Population sizes:\n")
cat(sprintf("  Unique oxycodone PATIENTS (Table 1 population):  n = %d\n", n_patients))
cat(sprintf("  All oxycodone EPISODES, S4 (all other tables):  n = %d\n", n_episodes))
cat("\n")

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

cat(sprintf("Reactive antiemetic (FINAL): n=%d, %.1f%%\n",
            sum(s4_oxy$reactive_antiemetic_final, na.rm=TRUE),
            100*mean(s4_oxy$reactive_antiemetic_final, na.rm=TRUE)))
cat(sprintf("Antipruritic (FINAL): n=%d, %.1f%%\n",
            sum(s4_oxy$any_antipruritic_final, na.rm=TRUE),
            100*mean(s4_oxy$any_antipruritic_final, na.rm=TRUE)))
cat(sprintf("Naloxone reversal (FINAL): n=%d\n\n",
            sum(s4_oxy$any_naloxone_reversal_final, na.rm=TRUE)))

# ============================================================
# TABLE 1 — DEMOGRAPHICS (patient-level, PAT_ID-deduplicated)
# ============================================================

cat("============================================================\n")
cat("TABLE 1 — OXYCODONE PATIENT DEMOGRAPHICS\n")
cat("(Unique patients, PAT_ID-deduplicated, earliest episode kept)\n")
cat("============================================================\n\n")

cat("Unique patients:", n_patients, "\n\n")

cat("Age — median (IQR) years:  ", fmt_med_iqr(oxycodone_patients$AGE), "\n")
cat("Weight — median (IQR) kg:  ", fmt_med_iqr(oxycodone_patients$WEIGHT_KG), "\n")
cat("Male — n (%):              ",
    fmt_n_pct(sum(oxycodone_patients$GENDER == "Male", na.rm = TRUE), n_patients), "\n\n")

cat("--- Age group ---\n")
oxycodone_patients |> count(age_group) |> mutate(pct = round(100 * n / n_patients, 1)) |> print()
cat("\n")

# ============================================================
# CASE TYPE
#
# Both the original, unrefined case_type (does any OPCS-coded
# operation exist for this admission — a blunt presence check) and
# case_type_refined are shown explicitly: the original for an audit
# trail, the refined version for actual use.
# ============================================================

cat("--- Case type (ORIGINAL — presence-of-any-operation only,\n")
cat("does not account for incidental/non-explanatory procedures) ---\n")
oxycodone_patients |> count(case_type) |> mutate(pct = round(100 * n / n_patients, 1)) |> print()
cat("\n")

cat("--- Case type (REFINED — USE THIS FOR REPORTING. Excludes\n")
cat("episodes where the labelled procedure is entirely minor/\n")
cat("incidental, has no preceding operation at all, or PCA\n")
cat("substantially predates the procedure by >2h — see 4_flags.R\n")
cat("case_type_refined for the full 5-mechanism derivation) ---\n")
oxycodone_patients |> count(case_type_refined) |> mutate(pct = round(100 * n / n_patients, 1)) |> print()
cat("\n")

cat("--- Admission type ---\n")
oxycodone_patients |> count(admission_type) |> mutate(pct = round(100 * n / n_patients, 1)) |>
  arrange(desc(n)) |> print()
cat("\n")

cat("--- Era ---\n")
oxycodone_patients |> count(era) |> mutate(pct = round(100 * n / n_patients, 1)) |> print()
cat("\n")

# ============================================================
# TABLE 2 — DOSING AND MODE
# ============================================================

cat("============================================================\n")
cat("TABLE 2 — DOSING AND MODE\n")
cat("============================================================\n\n")

cat("--- Mode distribution ---\n")
s4_oxy |> count(mode) |> mutate(pct = round(100 * n / n_episodes, 1)) |> print()
cat("\n")

cat("--- Episode duration ---\n")
cat("All S4:    ", fmt_med_iqr(s4_oxy$episode_duration_hrs), "hours\n")
cat("PCA only:  ", fmt_med_iqr(s4_oxy$episode_duration_hrs[s4_oxy$mode == "PCA"]), "hours\n")
cat("NCA only:  ", fmt_med_iqr(s4_oxy$episode_duration_hrs[s4_oxy$mode == "NCA"]), "hours\n\n")

cat("--- Duration by age group ---\n")
s4_oxy |>
  group_by(age_group) |>
  summarise(
    n = n(),
    median_hrs = round(median(episode_duration_hrs, na.rm = TRUE), 1),
    q25_hrs = round(quantile(episode_duration_hrs, 0.25, na.rm = TRUE), 1),
    q75_hrs = round(quantile(episode_duration_hrs, 0.75, na.rm = TRUE), 1),
    .groups = "drop"
  ) |> print()
cat("\n")

cat("--- First line vs rotation ---\n")
cat("First line (oxycodone as initial drug): ", fmt_n_pct(sum(s4_oxy$is_first_line), n_episodes), "\n")
cat("Rotation to oxycodone:                  ", fmt_n_pct(sum(!s4_oxy$is_first_line), n_episodes), "\n\n")

cat("--- Intraoperative drug flags ---\n")
cat("Remifentanil:       ", fmt_n_pct(sum(s4_oxy$had_remifentanil, na.rm = TRUE), n_episodes), "\n")
cat("Intraop ketamine:   ", fmt_n_pct(sum(s4_oxy$had_intraop_ketamine, na.rm = TRUE), n_episodes), "\n")
cat("Intraop clonidine:  ", fmt_n_pct(sum(s4_oxy$had_intraop_clonidine, na.rm = TRUE), n_episodes), "\n")
cat("Intrathecal opioid (spinal): ", fmt_n_pct(sum(s4_oxy$neuraxial_intrathecal_mar, na.rm = TRUE), n_episodes), "\n")
cat("Peripheral nerve block: ", fmt_n_pct(sum(s4_oxy$had_peripheral, na.rm = TRUE), n_episodes), "\n\n")

# ============================================================
# TABLE 3 — PROCEDURE CATEGORY
# ============================================================

cat("============================================================\n")
cat("TABLE 3 — PROCEDURE CATEGORY (S4, all oxycodone episodes)\n")
cat("============================================================\n\n")

s4_oxy |> count(proc_category) |> mutate(pct = round(100 * n / n_episodes, 1)) |>
  arrange(desc(n)) |> print()
cat("\n")

cat("Note: episodes with proc_category == NA (no preceding operation\n")
cat("found at all) are excluded from this breakdown by count()'s\n")
cat("default NA handling — check separately if this table should\n")
cat("explicitly show a 'No procedure found' row:\n")
cat("  sum(is.na(s4_oxy$proc_category)), \"/\", n_episodes\n\n")

# ============================================================
# TABLE 4 — SIDE EFFECTS BY MODE (FINAL definitions)
# ============================================================

cat("============================================================\n")
cat("TABLE 4 — SIDE EFFECTS BY MODE\n")
cat("(FINAL anaesthesia-stop-anchored definitions — see header)\n")
cat("============================================================\n\n")

se_by_mode <- s4_oxy |>
  group_by(mode) |>
  summarise(
    n                  = n(),
    pct_antiemetic     = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    pct_antipruritic   = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
    pct_side_effect    = round(100 * mean(any_side_effect_final, na.rm = TRUE), 1),
    n_nalox_reversal   = sum(any_naloxone_reversal_final, na.rm = TRUE),
    pct_nalox_reversal = round(100 * mean(any_naloxone_reversal_final, na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("Side effect rates by mode (S4, FINAL definitions):\n")
print(se_by_mode)
cat("\n")

cat("--- Overall (all modes) ---\n")
cat("Reactive antiemetic (FINAL): ", fmt_n_pct(sum(s4_oxy$reactive_antiemetic_final, na.rm = TRUE), n_episodes), "\n")
cat("Any antipruritic (FINAL):    ", fmt_n_pct(sum(s4_oxy$any_antipruritic_final, na.rm = TRUE), n_episodes), "\n")
cat("Any side effect (FINAL):     ", fmt_n_pct(sum(s4_oxy$any_side_effect_final, na.rm = TRUE), n_episodes), "\n")
cat("Naloxone reversal (FINAL):   ", fmt_n_pct(sum(s4_oxy$any_naloxone_reversal_final, na.rm = TRUE), n_episodes), "\n\n")

# ============================================================
# TABLE 5 — SIDE EFFECTS BY AGE GROUP
# ============================================================

cat("============================================================\n")
cat("TABLE 5 — SIDE EFFECTS BY AGE GROUP\n")
cat("============================================================\n\n")

se_by_age <- s4_oxy |>
  group_by(age_group) |>
  summarise(
    n                = n(),
    pct_antiemetic   = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    pct_antipruritic = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
    pct_side_effect  = round(100 * mean(any_side_effect_final,  na.rm = TRUE), 1),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(age_group)

print(se_by_age)
cat("\n")

# ============================================================
# TABLE 6 — SIDE EFFECTS BY ERA
# ============================================================

cat("============================================================\n")
cat("TABLE 6 — SIDE EFFECTS BY ERA\n")
cat("(Tests whether outcomes changed as oxycodone experience grew)\n")
cat("============================================================\n\n")

se_by_era <- s4_oxy |>
  group_by(era) |>
  summarise(
    n                = n(),
    pct_antiemetic   = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    pct_antipruritic = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
    pct_side_effect  = round(100 * mean(any_side_effect_final,  na.rm = TRUE), 1),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    .groups = "drop"
  )

print(se_by_era)
cat("\n")

# ============================================================
# TABLE 7 — PAIN OUTCOMES (S4) — now uses plain s4_oxy throughout
# ============================================================

cat("============================================================\n")
cat("TABLE 7 — PAIN OUTCOMES (S4, all ages, stratified by age group)\n")
cat("============================================================\n\n")

s4_shifts <- episode_shift_table |>
  semi_join(s4_oxy, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  filter(shift_in_episode >= 2)

cat("Total shifts analysed (shift 2+):", nrow(s4_shifts), "\n\n")

cat("--- Overall pain outcomes ---\n")
cat("Episodes:", n_episodes, "\n")
cat("Pct shifts severe (≥7) — median (IQR):   ", fmt_med_iqr(s4_oxy$pct_shifts_severe), "\n")
cat("Pct shifts moderate (≥4) — median (IQR): ", fmt_med_iqr(s4_oxy$pct_shifts_moderate), "\n")
cat("First controlled shift — median (IQR):    ", fmt_med_iqr(s4_oxy$first_controlled_shift), "\n\n")

cat("--- Pain outcomes by age group ---\n")
cat("(Note: zero rates in younger children reflect scoring unreliability,\n")
cat(" not necessarily absence of pain — interpret with caution)\n\n")

s4_oxy |>
  group_by(age_group) |>
  summarise(
    n                       = n(),
    pct_severe_median       = round(median(pct_shifts_severe,        na.rm = TRUE), 1),
    pct_severe_iqr_25       = round(quantile(pct_shifts_severe, 0.25, na.rm = TRUE), 1),
    pct_severe_iqr_75       = round(quantile(pct_shifts_severe, 0.75, na.rm = TRUE), 1),
    pct_moderate_median     = round(median(pct_shifts_moderate,      na.rm = TRUE), 1),
    pct_no_assess_median    = round(median(pct_shifts_no_assessment, na.rm = TRUE), 1),
    first_controlled_median = round(median(first_controlled_shift,   na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  arrange(age_group) |>
  print()
cat("\n")

cat("--- Pain outcomes by mode ---\n")
s4_oxy |>
  group_by(mode) |>
  summarise(
    n                    = n(),
    pct_severe_median    = round(median(pct_shifts_severe,  na.rm = TRUE), 1),
    pct_severe_q25       = round(quantile(pct_shifts_severe, 0.25, na.rm = TRUE), 1),
    pct_severe_q75       = round(quantile(pct_shifts_severe, 0.75, na.rm = TRUE), 1),
    pct_no_assessment    = round(median(pct_shifts_no_assessment, na.rm = TRUE), 1),
    .groups = "drop"
  ) |> print()
cat("\n")

# ============================================================
# TABLE 8 — CO-ANALGESIA (S4) — now uses plain s4_oxy
# ============================================================

cat("============================================================\n")
cat("TABLE 8 — CO-ANALGESIA (S4)\n")
cat("============================================================\n\n")

cat("Paracetamol — pct shifts given, median (IQR): ", fmt_med_iqr(s4_oxy$pct_shifts_paracetamol), "\n")
cat("NSAID — pct shifts given, median (IQR):        ", fmt_med_iqr(s4_oxy$pct_shifts_nsaid), "\n")
cat("Adjuvant — pct shifts given, median (IQR):     ", fmt_med_iqr(s4_oxy$pct_shifts_adjuvant), "\n")
cat("Neuropathic — pct shifts given, median (IQR):  ", fmt_med_iqr(s4_oxy$pct_shifts_neuropathic), "\n\n")

cat("--- Co-analgesia by mode ---\n")
s4_oxy |>
  group_by(mode) |>
  summarise(
    n                      = n(),
    paracetamol_median_pct = round(median(pct_shifts_paracetamol, na.rm = TRUE), 1),
    nsaid_median_pct       = round(median(pct_shifts_nsaid,        na.rm = TRUE), 1),
    adjuvant_median_pct    = round(median(pct_shifts_adjuvant,     na.rm = TRUE), 1),
    .groups = "drop"
  ) |> print()
cat("\n")

# ============================================================
# TABLE 9 — ROTATION INTO OXYCODONE
# ============================================================

cat("============================================================\n")
cat("TABLE 9 — ROTATION INTO OXYCODONE\n")
cat("============================================================\n\n")

oxy_rotations <- pca_episodes_flagged |>
  filter(drug == "oxycodone", is_rotation == TRUE)

cat("Total oxycodone rotation episodes:", nrow(oxy_rotations), "\n\n")

cat("--- Prior drug ---\n")
oxy_rotations |>
  count(prior_drug, sort = TRUE) |>
  mutate(pct = round(100 * n / nrow(oxy_rotations), 1)) |>
  print()
cat("\n")

cat("--- Rotation episode side effects vs first-line oxycodone (FINAL definitions, S4) ---\n")

s4_oxy |>
  group_by(is_first_line) |>
  summarise(
    n                = n(),
    pct_antiemetic   = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    pct_antipruritic = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
    pct_side_effect  = round(100 * mean(any_side_effect_final,  na.rm = TRUE), 1),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(group = if_else(is_first_line, "First line", "Rotation")) |>
  select(group, everything(), -is_first_line) |>
  print()
cat("\n")

# ============================================================
# TABLE 10 — MULTIPLE ROTATIONS
# ============================================================

cat("============================================================\n")
cat("TABLE 10 — MULTIPLE OPIOID ROTATIONS\n")
cat("============================================================\n\n")

rotations_per_patient <- pca_episodes_flagged |>
  filter(is_rotation == TRUE) |>
  count(PAT_ENC_CSN_ID, name = "n_rotations")

cat("Patients with at least one rotation:", nrow(rotations_per_patient), "\n\n")

rotation_dist <- rotations_per_patient |>
  count(n_rotations) |>
  mutate(pct = round(100 * n / nrow(rotations_per_patient), 1))

cat("Distribution of rotation counts per patient:\n")
print(rotation_dist)
cat("\n")

cat("Max rotations in one patient:", max(rotations_per_patient$n_rotations), "\n")
cat("Patients with 3+ rotations:", sum(rotations_per_patient$n_rotations >= 3), "\n\n")

cat("--- Most complex rotation sequences (patients with 3+ switches) ---\n")
complex_patients <- rotations_per_patient |>
  filter(n_rotations >= 3) |>
  pull(PAT_ENC_CSN_ID)

pca_episodes_flagged |>
  filter(PAT_ENC_CSN_ID %in% complex_patients) |>
  arrange(PAT_ENC_CSN_ID, episode_start) |>
  group_by(PAT_ENC_CSN_ID) |>
  summarise(
    sequence = paste(drug, collapse = " → "),
    n_drugs  = n_distinct(drug),
    .groups  = "drop"
  ) |>
  count(sequence, sort = TRUE) |>
  head(10) |>
  print()
cat("\n")

# ============================================================
# SUPPLEMENTARY — S4 CONSUMPTION AND SURGICAL CLASSIFICATION
# (already correctly using s4_surgical_keys/s4_medical_keys, which
# are built from case_type_refined in 6_cohort_classification.R —
# this section was already consistent; Table 1 above was the one
# needing the fix)
# ============================================================

cat("============================================================\n")
cat("SUPPLEMENTARY — OUTCOMES BY SURGICAL CLASSIFICATION (S4)\n")
cat("============================================================\n\n")

cat("By reason for PCA:\n")
n_surgical <- sum(s4_oxy$PAT_ENC_CSN_ID %in% s4_surgical_keys$PAT_ENC_CSN_ID &
                    s4_oxy$episode_id %in% s4_surgical_keys$episode_id)
n_nonsurgical <- n_episodes - n_surgical

cat(sprintf("  Surgical (major OPCS):     %d (%.1f%%)\n", n_surgical, 100*n_surgical/n_episodes))
cat(sprintf("  Non-surgical (minor/none): %d (%.1f%%)\n\n", n_nonsurgical, 100*n_nonsurgical/n_episodes))

if (!is.null(daily_consumption)) {
  cat("Oxycodone consumption (mcg/kg/hr):\n")
  consume_overall <- daily_consumption |>
    summarise(median = round(median(mcg_per_kg_hr, na.rm=TRUE), 2),
              q25 = round(quantile(mcg_per_kg_hr, 0.25, na.rm=TRUE), 2),
              q75 = round(quantile(mcg_per_kg_hr, 0.75, na.rm=TRUE), 2))
  cat(sprintf("  Overall: %.2f (%.2f–%.2f) mcg/kg/hr\n",
              consume_overall$median, consume_overall$q25, consume_overall$q75))
  
  consume_surg <- daily_consumption |>
    mutate(is_surgical = (paste(PAT_ENC_CSN_ID, episode_id, sep="_") %in%
                            paste(s4_surgical_keys$PAT_ENC_CSN_ID, s4_surgical_keys$episode_id, sep="_"))) |>
    filter(is_surgical == TRUE) |>
    summarise(median = round(median(mcg_per_kg_hr, na.rm=TRUE), 2),
              q25 = round(quantile(mcg_per_kg_hr, 0.25, na.rm=TRUE), 2),
              q75 = round(quantile(mcg_per_kg_hr, 0.75, na.rm=TRUE), 2))
  cat(sprintf("  Surgical: %.2f (%.2f–%.2f) mcg/kg/hr\n",
              consume_surg$median, consume_surg$q25, consume_surg$q75))
  
  consume_nonsurg <- daily_consumption |>
    mutate(is_surgical = (paste(PAT_ENC_CSN_ID, episode_id, sep="_") %in%
                            paste(s4_surgical_keys$PAT_ENC_CSN_ID, s4_surgical_keys$episode_id, sep="_"))) |>
    filter(is_surgical == FALSE) |>
    summarise(median = round(median(mcg_per_kg_hr, na.rm=TRUE), 2),
              q25 = round(quantile(mcg_per_kg_hr, 0.25, na.rm=TRUE), 2),
              q75 = round(quantile(mcg_per_kg_hr, 0.75, na.rm=TRUE), 2))
  cat(sprintf("  Non-surgical: %.2f (%.2f–%.2f) mcg/kg/hr\n\n",
              consume_nonsurg$median, consume_nonsurg$q25, consume_nonsurg$q75))
} else {
  cat("Consumption statistics skipped — daily_consumption not available\n\n")
}

s4_surgical_metrics <- s4_oxy |>
  semi_join(s4_surgical_keys, by = c("PAT_ENC_CSN_ID", "episode_id"))
s4_nonsurgical_metrics <- s4_oxy |>
  anti_join(s4_surgical_keys, by = c("PAT_ENC_CSN_ID", "episode_id"))

cat("S4 Surgical (major OPCS, n =", nrow(s4_surgical_metrics), "):\n")
cat("  Reactive antiemetic (FINAL): ",
    fmt_n_pct(sum(s4_surgical_metrics$reactive_antiemetic_final, na.rm = TRUE),
              nrow(s4_surgical_metrics)), "\n")
cat("  Any side effect (FINAL):     ",
    fmt_n_pct(sum(s4_surgical_metrics$any_side_effect_final, na.rm = TRUE),
              nrow(s4_surgical_metrics)), "\n\n")

cat("S4 Non-surgical (minor/none, n =", nrow(s4_nonsurgical_metrics), "):\n")
cat("  Reactive antiemetic (FINAL): ",
    fmt_n_pct(sum(s4_nonsurgical_metrics$reactive_antiemetic_final, na.rm = TRUE),
              nrow(s4_nonsurgical_metrics)), "\n")
cat("  Any side effect (FINAL):     ",
    fmt_n_pct(sum(s4_nonsurgical_metrics$any_side_effect_final, na.rm = TRUE),
              nrow(s4_nonsurgical_metrics)), "\n\n")

cat("(Note: genuine propensity-matched oxycodone-vs-morphine comparison\n")
cat("is reported separately in 10_comparative_PSM.R — this section shows\n")
cat("unadjusted S4 descriptive rates only, not a matched comparison.)\n\n")

# ============================================================
# TABLE 11 — OXYCODONE DOSING (S4)
# ============================================================

cat("============================================================\n")
cat("TABLE 11 — OXYCODONE DOSING (S4, from 7_dosing.R)\n")
cat("Requires episode_dose.rds from 7_dosing.R\n")
cat("============================================================\n\n")

episode_dose <- tryCatch(
  readRDS("1 - data/3 - processed_data/episode_dose.rds"),
  error = function(e) {
    cat("episode_dose.rds not found — run 7_dosing.R first\n")
    return(NULL)
  }
)

if (!is.null(episode_dose)) {
  
  oxy_dose <- episode_dose |>
    filter(drug == "oxycodone") |>
    semi_join(s4_oxy, by = c("PAT_ENC_CSN_ID", "episode_id"))
  
  cat(sprintf("Oxycodone episodes with dose data: %d / %d\n\n", nrow(oxy_dose), n_episodes))
  
  cat("--- Overall dosing (S4) ---\n")
  cat("Background rate (mg/kg/hr) — median (IQR): ", fmt_med_iqr(oxy_dose$background_rate_mg_kg_hr_median, 4), "\n")
  cat("Bolus dose (mg/kg) — median (IQR):          ", fmt_med_iqr(oxy_dose$bolus_dose_mg_kg_first, 4), "\n")
  cat("Daily background dose (mg/kg/day) — median (IQR): ", fmt_med_iqr(oxy_dose$daily_background_mg_kg, 3), "\n")
  cat("Loading dose prescribed — n (%):            ",
      fmt_n_pct(sum(oxy_dose$had_loading_dose, na.rm = TRUE), nrow(oxy_dose)), "\n\n")
  
  cat("--- Dosing by age group ---\n")
  oxy_dose |>
    group_by(age_group) |>
    summarise(
      n = n(),
      background_mg_kg_hr = fmt_med_iqr(background_rate_mg_kg_hr_median, 4),
      bolus_mg_kg = fmt_med_iqr(bolus_dose_mg_kg_first, 4),
      daily_bg_mg_kg = fmt_med_iqr(daily_background_mg_kg, 3),
      pct_loading = round(100 * mean(had_loading_dose, na.rm = TRUE), 1),
      conc_mg_per_ml = round(median(concentration_mg_per_ml, na.rm = TRUE), 3),
      .groups = "drop"
    ) |> arrange(age_group) |> print()
  cat("\n")
  
  cat("--- Dosing by mode ---\n")
  oxy_dose |>
    left_join(s4_oxy |> select(PAT_ENC_CSN_ID, episode_id, mode), by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    group_by(mode) |>
    summarise(
      n = n(),
      background_mg_kg_hr = fmt_med_iqr(background_rate_mg_kg_hr_median, 4),
      bolus_mg_kg = fmt_med_iqr(bolus_dose_mg_kg_first, 4),
      daily_bg_mg_kg = fmt_med_iqr(daily_background_mg_kg, 3),
      pct_loading = round(100 * mean(had_loading_dose, na.rm = TRUE), 1),
      .groups = "drop"
    ) |> print()
  cat("\n")
  
  cat("--- Dosing by era ---\n")
  oxy_dose |>
    left_join(s4_oxy |> select(PAT_ENC_CSN_ID, episode_id, era), by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    group_by(era) |>
    summarise(
      n = n(),
      background_mg_kg_hr = fmt_med_iqr(background_rate_mg_kg_hr_median, 4),
      bolus_mg_kg = fmt_med_iqr(bolus_dose_mg_kg_first, 4),
      daily_bg_mg_kg = fmt_med_iqr(daily_background_mg_kg, 3),
      .groups = "drop"
    ) |> print()
  cat("\n")
  
  cat("--- Concentration used by age group ---\n")
  oxy_dose |>
    group_by(age_group) |>
    summarise(
      n = n(),
      conc_median = round(median(concentration_mg_per_ml, na.rm = TRUE), 3),
      conc_min = round(min(concentration_mg_per_ml, na.rm = TRUE), 3),
      conc_max = round(max(concentration_mg_per_ml, na.rm = TRUE), 3),
      .groups = "drop"
    ) |> arrange(age_group) |> print()
  cat("\n")
  
} else {
  cat("Dosing section skipped — run 7_dosing.R first\n\n")
}

cat("============================================================\n")
cat("ANALYSIS 5 COMPLETE\n")
cat("Table 1: unique patients (n=", n_patients, "). All other tables:\n")
cat("S4, all episodes (n=", n_episodes, "). FINAL side-effect definitions\n")
cat("throughout, shared with 10_comparative_PSM.R and global reporting.\n")
cat("PSM-matched comparison: 10_comparative_PSM.R (not duplicated here).\n")
cat("============================================================\n")