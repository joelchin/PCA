# 16 — Publication Tables
# Paediatric PCA Study
# Addenbrooke's Hospital
#
# Builds the manuscript's five publication tables (Table 1: cohort
# characteristics; Table 2: descriptive efficacy/consumption; Table
# 3: drug switching; Table 4: PSM-matched comparative outcomes;
# Table 5: rotation) as flextable objects ready for Word export.
#
# A few methodological choices worth noting:
#
#   Table 1's case type field is case_type_refined, not the original
#   unrefined case_type — the unrefined field produces a persistent
#   92%/8% split that doesn't move across subgroups and is not a
#   reliable case-type signal.
#
#   Table 2 reports genuine efficacy metrics (peak pain score,
#   background consumption) rather than the SI (Service Improvement)
#   project's osteotomy-derived pain metrics (% shifts severe, shift
#   to pain control) — those belong to a different analytical
#   question and population. Its procedure-category breakdown uses
#   proc_category_display: surgical categories are shown
#   individually and non-surgical episodes are collapsed to one row,
#   mirroring Table 1's case-type logic. The composite side-effect
#   outcome is not shown here — it's dropped from the comparative
#   claim entirely (see Table 4 below) and showing it descriptively
#   here without the comparative context would be misleading.
#   "Shifts with paracetamol/NSAID, %" rows are also not included:
#   paracetamol co-prescribing is ~100% across every group (zero
#   variance, uninformative for this comparison) and both rows sit
#   more naturally within the SI pain-management-quality project than
#   this manuscript's core oxycodone-vs-morphine question.
#
#   A supplementary, descriptive-only comparison of morphine NCA
#   outcomes against Howard et al.'s prospective series of 10,079
#   morphine NCA patients is run immediately after Table 1, to
#   validate the reliability of our outcome definitions and data
#   extraction pipeline against external published data (reported in
#   Results; not a manuscript table).
#
#   Table 4 reports antiemetic, antipruritic, and naloxone reversal
#   as three independent outcomes rather than a combined composite
#   side-effect measure, since a composite assumes equal severity
#   across component events, which doesn't hold clinically.
#
#   pre_match, post_match, and table4_matched_outcomes are built
#   directly from the live pipeline objects (psm_data, matched_data)
#   rather than independently reconstructed from saved keys — this
#   keeps this script's population figures mechanically identical to
#   whatever 10_comparative_PSM.R actually produced, rather than
#   risking two independently-maintained copies of the same
#   patient-dedup/neonate-exclusion/gender-recode logic drifting out
#   of sync with each other. add_psm_covariates() is applied
#   asymmetrically: matched_data does not carry specialty_collapsed,
#   asa_grade, the intraop flags, or anaesthesia_duration_mins, so
#   post_match needs the rebuild; psm_data already carries all ten
#   covariates, so applying it to pre_match would instead collide
#   with the existing columns and produce duplicate .x/.y columns.
#
#   Table 5 draws on four objects from 15_rotation_analysis.R (tbl5,
#   ft5b, ft5c, ft5d). ft5d is the covariate-adjusted rotation-away
#   odds ratio, computed on the pre-match eligible pool (psm_data)
#   rather than the much smaller 237/237 matched sample, since the
#   matched sample has too few rotation-away events (25 total) to
#   support a stable adjusted estimate.
#
#   case_type_refined is reported at patient level in Table 1
#   (79.3%/20.7%), not episode level.
#
#   Naloxone-reversal factor() calls specify levels = c(FALSE, TRUE)
#   explicitly, since oxycodone's S4 population has zero reversal
#   events and factor() without explicit levels would otherwise only
#   detect one level and error.
# ============================================================
required_pkgs <- c("gtsummary", "flextable", "officer", "labelled")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n\n")
  install.packages(missing_pkgs)
  still_missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop("Install failed for: ", paste(still_missing, collapse = ", "),
         "\nCheck for errors above (e.g. CRAN mirror not set, no internet access on this remote desktop).")
  }
}

library(tidyverse)
library(gtsummary)
library(flextable)
library(officer)

cat("============================================================\n")
cat("ANALYSIS 11 — PUBLICATION TABLES\n")
cat("============================================================\n\n")

stopifnot(
  exists("pca_episodes_flagged"),
  exists("episode_metrics"),
  exists("s3_morphine_keys"),
  exists("s3_oxycodone_keys"),
  exists("admissions_clean"),
  exists("intraop_mar_clean"),
  exists("anaesthesia_events_clean"),
  exists("anaesthesia_ax_clean"),
  exists("pain_scores_clean")
)

# psm_data and matched_data must exist in-session (built by
# 10_comparative_PSM.R, run earlier in the same session, NOT loaded
# from a saved .rds). This is the fix for the pre_match/psm_data
# drift identified during validation.
stopifnot(
  exists("psm_data"),
  exists("matched_data")
)

stopifnot("case_type_refined" %in% names(pca_episodes_flagged))

stopifnot(
  "peak_pain_score" %in% names(pca_episodes_flagged),
  "time_weighted_mean_pain" %in% names(pca_episodes_flagged)
)
cat("Pain score outcomes confirmed pre-computed (from 4_flags.R):\n")
cat("  peak_pain_score, time_weighted_mean_pain\n\n")

# ============================================================
# FLEXTABLE THEME
# ============================================================
ft_theme <- function(ft) {
  ft |>
    theme_booktabs() |>
    fontsize(size = 10, part = "all") |>
    font(fontname = "Arial", part = "all") |>
    bold(part = "header") |>
    align(align = "left", part = "all") |>
    align(j = -1, align = "center", part = "body") |>
    padding(padding = 3, part = "all") |>
    set_table_properties(width = 1, layout = "autofit")
}

# ============================================================
# TABLE 1 — COHORT DEMOGRAPHICS (chemotherapy row REMOVED)
# case_type_refined at PATIENT-level — confirmed, not changed.
# ============================================================
cat("Building Table 1 — cohort demographics...\n")

patients <- pca_episodes_flagged |>
  left_join(admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID),
            by = "PAT_ENC_CSN_ID") |>
  arrange(PAT_ID, episode_start) |>
  distinct(PAT_ID, .keep_all = TRUE) |>
  mutate(
    age_group = factor(age_group, levels = c(
      "Neonate (<28 days)", "Infant (28 days-1yr)",
      "Young child (1-5yr)", "Child (6-11yr)", "Adolescent (12-15yr)"
    )),
    case_type_refined = factor(case_type_refined, levels = c("surgical", "medical"),
                               labels = c("Surgical", "Medical")),
    GENDER = factor(GENDER, levels = c("Male", "Female"))
  )

cat("Patients before dedup (by admission, PAT_ENC_CSN_ID):",
    n_distinct(pca_episodes_flagged$PAT_ENC_CSN_ID), "\n")
cat("Patients after dedup (by true patient, PAT_ID):     ",
    nrow(patients), "\n\n")

episodes_summary <- pca_episodes_flagged |>
  mutate(
    drug = factor(drug, levels = c("morphine", "oxycodone",
                                   "morphine+ketamine", "fentanyl"),
                  labels = c("Morphine", "Oxycodone",
                             "Morphine + ketamine", "Fentanyl")),
    mode = factor(mode, levels = c("PCA", "NCA", "unknown"),
                  labels = c("PCA", "NCA", "Unknown")),
    case_type_refined = factor(case_type_refined, levels = c("surgical", "medical"),
                               labels = c("Surgical", "Medical"))
  )

tbl1_patients <- patients |>
  select(AGE, WEIGHT_KG, GENDER, LOS_days, age_group, case_type_refined) |>
  tbl_summary(
    label = list(
      AGE ~ "Age, years",
      WEIGHT_KG ~ "Weight, kg",
      GENDER ~ "Sex",
      LOS_days ~ "Length of stay, days",
      age_group ~ "Age group",
      case_type_refined ~ "Case type"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}–{p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1),
    missing = "no"
  ) |>
  modify_header(label = "**Characteristic**") |>
  modify_caption("**Table 1.** Patient and episode characteristics") |>
  bold_labels()

tbl1_episodes <- episodes_summary |>
  select(drug, mode, case_type_refined, episode_duration_hrs,
         concurrent_epidural, is_rotation) |>
  mutate(
    concurrent_epidural = factor(concurrent_epidural,
                                 levels = c(TRUE, FALSE),
                                 labels = c("Yes", "No")),
    is_rotation = factor(is_rotation,
                         levels = c(TRUE, FALSE),
                         labels = c("Yes", "No"))
  ) |>
  tbl_summary(
    label = list(
      drug ~ "Drug",
      mode ~ "Mode",
      case_type_refined ~ "Case type",
      episode_duration_hrs ~ "Episode duration, hours",
      concurrent_epidural ~ "Concurrent epidural",
      is_rotation ~ "Opioid rotation episode"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}–{p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1),
    missing = "no"
  ) |>
  modify_header(label = "**Characteristic**") |>
  bold_labels()

tbl1 <- tbl_stack(
  list(tbl1_patients, tbl1_episodes),
  group_header = c(sprintf("Patient-level (n=%d patients)", nrow(patients)),
                   sprintf("Episode-level (n=%d episodes)", nrow(episodes_summary)))
)
cat("Table 1 built\n\n")

# ============================================================
# SUPPLEMENTARY — MORPHINE NCA VALIDATION AGAINST HOWARD ET AL.
# (descriptive only, console output — not a manuscript table,
# same treatment as the surgical vs medical comparison below)
# ============================================================
cat("============================================================\n")
cat("SUPPLEMENTARY — MORPHINE NCA VALIDATION vs HOWARD ET AL. 2010\n")
cat("============================================================\n\n")

morphine_nca <- pca_episodes_flagged |>
  filter(drug == "morphine", mode == "NCA") |>
  semi_join(episode_metrics |> filter(drug == "morphine"),
            by = c("PAT_ENC_CSN_ID", "episode_id"))

cat("Morphine NCA episodes, n =", nrow(morphine_nca), "\n")
cat("Age, years: median", median(morphine_nca$AGE, na.rm = TRUE),
    "(", quantile(morphine_nca$AGE, 0.25, na.rm = TRUE), "-",
    quantile(morphine_nca$AGE, 0.75, na.rm = TRUE), ")\n")
cat("Howard et al. reference: n=10,079, median age 2.3 years\n\n")

cat("Antiemetic (reactive, FINAL definition):",
    round(100 * mean(morphine_nca$reactive_antiemetic_final, na.rm = TRUE), 1),
    "% (Howard: 25% any PONV / 13.5% clinically severe)\n")
cat("Antipruritic:",
    round(100 * mean(morphine_nca$any_antipruritic_final, na.rm = TRUE), 1),
    "% (Howard: 9.4% any itching)\n")
cat("Naloxone reversal:",
    round(100 * mean(morphine_nca$any_naloxone_reversal_final, na.rm = TRUE), 2),
    "% (Howard: 0.4% SAE)\n\n")

morphine_nca_dose <- morphine_nca |>
  left_join(
    tryCatch(
      readRDS("1 - data/3 - processed_data/episode_true_total_dose.rds"),
      error = function(e) {
        cat("episode_true_total_dose.rds not found — consumption validation skipped\n")
        NULL
      }
    ) |> select(PAT_ENC_CSN_ID, episode_id, total_dose_mcg_kg_hr, background_mcg_kg_hr),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  )
cat("Total dose (background+bolus), mcg/kg/hr: median",
    round(median(morphine_nca_dose$total_dose_mcg_kg_hr, na.rm = TRUE), 1),
    "(Howard D0: 15.65)\n\n")

cat("(Reported descriptively in Results — see manuscript. Not a\n")
cat("formal table; no test of statistical significance performed,\n")
cat("as this is external validation of methodology, not a formal\n")
cat("comparative outcome.)\n\n")

# ============================================================
# TABLE 2 — OXYCODONE DESCRIPTIVE (S4, all episodes)
# "Shifts with paracetamol/NSAID" REMOVED this pass — see header note.
# ============================================================
cat("Building Table 2 — oxycodone descriptive (S4)...\n")

s4_oxy_t2 <- pca_episodes_flagged |>
  filter(drug == "oxycodone") |>
  semi_join(episode_metrics |> filter(drug == "oxycodone"),
            by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(
    mode = factor(mode, levels = c("PCA", "NCA"), labels = c("PCA", "NCA")),
    age_group = factor(age_group, levels = c(
      "Neonate (<28 days)", "Infant (28 days-1yr)",
      "Young child (1-5yr)", "Child (6-11yr)", "Adolescent (12-15yr)"
    )),
    proc_category_display = if_else(
      case_type_refined == "surgical",
      as.character(proc_category),
      "Non-surgical"
    )
  ) |>
  filter(mode %in% c("PCA", "NCA"))

cat("Peak pain score present (S4):", sum(!is.na(s4_oxy_t2$peak_pain_score)), "/", nrow(s4_oxy_t2), "\n")
cat("Time-weighted mean pain score present (S4):",
    sum(!is.na(s4_oxy_t2$time_weighted_mean_pain)), "/", nrow(s4_oxy_t2), "\n")

episode_dose_for_consumption <- tryCatch(
  readRDS("1 - data/3 - processed_data/episode_dose.rds"),
  error = function(e) {
    cat("episode_dose.rds not found — consumption row will be NA\n")
    return(NULL)
  }
)

if (!is.null(episode_dose_for_consumption)) {
  s4_consumption <- episode_dose_for_consumption |>
    filter(drug == "oxycodone") |>
    mutate(mcg_per_kg_hr = (total_background_mg * 1000) / (WEIGHT_KG * episode_duration_dose_hrs)) |>
    select(PAT_ENC_CSN_ID, episode_id, mcg_per_kg_hr)
  
  s4_oxy_t2 <- s4_oxy_t2 |>
    left_join(s4_consumption, by = c("PAT_ENC_CSN_ID", "episode_id"))
} else {
  s4_oxy_t2$mcg_per_kg_hr <- NA_real_
}

episode_true_total_dose_t2 <- tryCatch(
  readRDS("1 - data/3 - processed_data/episode_true_total_dose.rds"),
  error = function(e) {
    cat("episode_true_total_dose.rds not found — total dose row will be NA\n")
    NULL
  }
)

if (!is.null(episode_true_total_dose_t2)) {
  s4_oxy_t2 <- s4_oxy_t2 |>
    left_join(
      episode_true_total_dose_t2 |> select(PAT_ENC_CSN_ID, episode_id, total_dose_mcg_kg_hr),
      by = c("PAT_ENC_CSN_ID", "episode_id")
    )
} else {
  s4_oxy_t2$total_dose_mcg_kg_hr <- NA_real_
}

cat("Total dose (background+bolus) joined (S4):",
    sum(!is.na(s4_oxy_t2$total_dose_mcg_kg_hr)), "/", nrow(s4_oxy_t2), "\n")

cat("Peak pain score joined (S4):", sum(!is.na(s4_oxy_t2$peak_pain_score)), "/", nrow(s4_oxy_t2), "\n")
cat("Consumption joined (S4):", sum(!is.na(s4_oxy_t2$mcg_per_kg_hr)), "/", nrow(s4_oxy_t2), "\n\n")

tbl2 <- s4_oxy_t2 |>
  select(mode, AGE, WEIGHT_KG, GENDER, age_group,
         episode_duration_hrs, proc_category_display,
         reactive_antiemetic_final, any_antipruritic_final,
         any_naloxone_reversal_final,
         peak_pain_score, time_weighted_mean_pain,
         total_dose_mcg_kg_hr, mcg_per_kg_hr) |>
  mutate(
    reactive_antiemetic_final = factor(reactive_antiemetic_final, levels = c(FALSE, TRUE), labels = c("No", "Yes")),
    any_antipruritic_final = factor(any_antipruritic_final, levels = c(FALSE, TRUE), labels = c("No", "Yes")),
    any_naloxone_reversal_final = factor(any_naloxone_reversal_final, levels = c(FALSE, TRUE), labels = c("No", "Yes"))
  ) |>
  tbl_summary(
    by = mode,
    label = list(
      AGE ~ "Age, years",
      WEIGHT_KG ~ "Weight, kg",
      GENDER ~ "Sex",
      age_group ~ "Age group",
      episode_duration_hrs ~ "Episode duration, hours",
      proc_category_display ~ "Procedure category (surgical only; non-surgical episodes shown as one group)",
      reactive_antiemetic_final ~ "Any antiemetic (reactive, FINAL definition) — n (%)",
      any_antipruritic_final ~ "Any antipruritic — n (%)",
      any_naloxone_reversal_final ~ "Naloxone reversal — n (%)",
      peak_pain_score ~ "Peak pain score, median (IQR) [efficacy]",
      time_weighted_mean_pain ~ "Time-weighted mean pain score, median (IQR) [efficacy]",
      total_dose_mcg_kg_hr ~ "Total dose (background + bolus), mcg/kg/hr [efficacy]",
      mcg_per_kg_hr ~ "  (reference) Background only, mcg/kg/hr"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}–{p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1),
    missing = "no"
  ) |>
  add_overall(last = FALSE) |>
  modify_header(label = "**Outcome**",
                stat_0 = sprintf("**Overall (n=%d)**", nrow(s4_oxy_t2)),
                stat_1 = sprintf("**PCA (n=%d)**", sum(s4_oxy_t2$mode == "PCA")),
                stat_2 = sprintf("**NCA (n=%d)**", sum(s4_oxy_t2$mode == "NCA"))) |>
  modify_caption(sprintf("**Table 2.** Oxycodone PCA/NCA episode characteristics (all episodes, S4, n=%d)",
                         nrow(s4_oxy_t2))) |>
  bold_labels() |>
  modify_footnote(
    all_stat_cols() ~ "Continuous: median (IQR). Categorical: n (%). Antiemetic: >12h post-anaesthesia-stop AND no antiemetic in the 7 days before anaesthesia start (surgical); any antiemetic during episode AND no antiemetic in the 7 days before episode start (medical, no anaesthesia); episodes with unreliable anaesthesia timing excluded from both antiemetic measures. Antipruritic: no time filter. Naloxone reversal: manually chart-reviewed, minimal/no cutoff. Time-weighted mean pain score weights each reading by the time it applied until the next reading (see Methods); excludes time before the first documented reading. Total dose (background + bolus) restricted to episodes with >=80% pump flowsheet coverage — a smaller subset than the full S4 population shown for other rows. Population: S4 (all oxycodone episodes). NOTE: this table's n is smaller than the 1,276 total oxycodone episodes shown in Table 1's episode-level section — 78 episodes with undetermined PCA/NCA mode are excluded here, and a further small number lack complete episode-level metrics data required for other rows in this table; both are real, deliberate scoping differences, not a data inconsistency."
  )
cat("Table 2 built\n\n")

cat("============================================================\n")
cat("SUPPLEMENTARY — SURGICAL vs MEDICAL COMPARISON (S4, descriptive)\n")
cat("============================================================\n\n")

surgical_vs_medical <- s4_oxy_t2 |>
  mutate(
    population = if_else(
      PAT_ENC_CSN_ID %in% s4_surgical_keys$PAT_ENC_CSN_ID &
        episode_id %in% s4_surgical_keys$episode_id,
      "Surgical", "Medical"
    )
  )

cat("--- N by group ---\n")
print(table(surgical_vs_medical$population))
cat("\n")

cat("--- Efficacy: peak pain score, background consumption ---\n")
surgical_vs_medical |>
  group_by(population) |>
  summarise(
    n = n(),
    peak_pain_median = round(median(peak_pain_score, na.rm = TRUE), 1),
    peak_pain_q25 = round(quantile(peak_pain_score, 0.25, na.rm = TRUE), 1),
    peak_pain_q75 = round(quantile(peak_pain_score, 0.75, na.rm = TRUE), 1),
    consumption_median = round(median(mcg_per_kg_hr, na.rm = TRUE), 2),
    consumption_q25 = round(quantile(mcg_per_kg_hr, 0.25, na.rm = TRUE), 2),
    consumption_q75 = round(quantile(mcg_per_kg_hr, 0.75, na.rm = TRUE), 2),
    .groups = "drop"
  ) |>
  print()

wt_peak <- wilcox.test(peak_pain_score ~ population, data = surgical_vs_medical)
wt_consumption <- wilcox.test(mcg_per_kg_hr ~ population, data = surgical_vs_medical)
cat("\nWilcoxon peak pain score, surgical vs medical: p =", format.pval(wt_peak$p.value, digits = 3), "\n")
cat("Wilcoxon consumption, surgical vs medical: p =", format.pval(wt_consumption$p.value, digits = 3), "\n\n")

cat("--- Side effects: surgical vs medical ---\n")
surgical_vs_medical |>
  group_by(population) |>
  summarise(
    n = n(),
    pct_antiemetic = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    pct_antipruritic = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    .groups = "drop"
  ) |>
  print()

antiemetic_tab <- table(surgical_vs_medical$population, surgical_vs_medical$reactive_antiemetic_final)
antipruritic_tab <- table(surgical_vs_medical$population, surgical_vs_medical$any_antipruritic_final)

if (ncol(antiemetic_tab) == 2) {
  ct1 <- suppressWarnings(chisq.test(antiemetic_tab))
  cat("Chi-squared, antiemetic, surgical vs medical: p =", format.pval(ct1$p.value, digits = 3), "\n")
}
if (ncol(antipruritic_tab) == 2) {
  ct2 <- suppressWarnings(chisq.test(antipruritic_tab))
  cat("Chi-squared, antipruritic, surgical vs medical: p =", format.pval(ct2$p.value, digits = 3), "\n")
}
cat("\n")

# ============================================================
# TABLE 3 — PSM BALANCE (asa_grade included, dynamic headers)
# ============================================================
cat("Building Table 3 — PSM balance...\n")

pre_match <- psm_data |>
  mutate(
    drug = factor(if_else(treated == 1, "Oxycodone", "Morphine"),
                  levels = c("Morphine", "Oxycodone")),
    match_status = "Pre-match"
  )

post_match <- matched_data |>
  mutate(
    drug = factor(if_else(treated == 1, "Oxycodone", "Morphine"),
                  levels = c("Morphine", "Oxycodone")),
    match_status = "Post-match",
    episode_duration_hrs = episode_duration_hrs.x
  ) |>
  select(-episode_duration_hrs.x, -episode_duration_hrs.y)

add_psm_covariates <- function(df) {
  same_service_map <- c(
    "NEUROSURGERY"           = "PAEDIATRIC NEUROSURGERY",
    "TRAUMA & ORTHOPAEDICS"  = "PAEDIATRIC TRAUMA AND ORTHOPAEDICS",
    "UROLOGY"                = "PAEDIATRIC UROLOGY",
    "PLASTIC SURGERY"        = "PAEDIATRIC PLASTIC SURGERY",
    "GENERAL SURGERY"        = "PAEDIATRIC SURGERY",
    "MAXILLO-FACIAL SURGERY" = "PAEDIATRIC MAXILLO-FACIAL SURGERY",
    "ENT"                    = "PAEDIATRIC EAR NOSE AND THROAT"
  )
  SPECIALTY_VOLUME_THRESHOLD <- 30
  
  df <- df |>
    mutate(
      specialty_merged = recode(as.character(ADMITTING_SPECIALTY),
                                !!!same_service_map,
                                .default = as.character(ADMITTING_SPECIALTY))
    )
  keep_specialties <- df |>
    count(specialty_merged) |>
    filter(n >= SPECIALTY_VOLUME_THRESHOLD) |>
    pull(specialty_merged)
  df <- df |>
    mutate(
      specialty_collapsed = if_else(specialty_merged %in% keep_specialties,
                                    specialty_merged, "Other specialty"),
      proc_category = if_else(proc_category == "Other", "Unclassified", proc_category)
    )
  
  new_intraop_flags <- intraop_mar_clean |>
    group_by(PAT_ENC_CSN_ID) |>
    summarise(
      had_intraop_paracetamol = any(pharm_subclass == "Paracetamol", na.rm = TRUE),
      had_intraop_nsaid       = any(pharm_subclass == "NSAIDs", na.rm = TRUE),
      had_intraop_magnesium   = any(pharm_subclass == "Magnesium", na.rm = TRUE),
      had_intraop_fentanyl    = any(str_detect(str_to_lower(generic_name), "^fentanyl$"), na.rm = TRUE),
      had_intraop_morphine    = any(str_detect(str_to_lower(generic_name), "morphine sulfate"), na.rm = TRUE),
      had_intraop_oxycodone2  = any(str_detect(str_to_lower(generic_name), "oxycodone"), na.rm = TRUE),
      had_intraop_diamorphine = any(str_detect(str_to_lower(generic_name), "diamorphine"), na.rm = TRUE),
      .groups = "drop"
    )
  
  df <- df |>
    left_join(new_intraop_flags, by = "PAT_ENC_CSN_ID") |>
    mutate(across(c(had_intraop_paracetamol, had_intraop_nsaid, had_intraop_magnesium,
                    had_intraop_fentanyl, had_intraop_morphine, had_intraop_oxycodone2,
                    had_intraop_diamorphine),
                  ~ replace_na(.x, FALSE)))
  
  duration_data <- anaesthesia_events_clean |>
    select(PAT_ENC_CSN_ID, anaesthesia_duration_mins) |>
    distinct(PAT_ENC_CSN_ID, .keep_all = TRUE)
  
  df <- df |> left_join(duration_data, by = "PAT_ENC_CSN_ID")
  median_duration <- median(df$anaesthesia_duration_mins, na.rm = TRUE)
  df <- df |>
    mutate(anaesthesia_duration_mins = if_else(
      is.na(anaesthesia_duration_mins), median_duration, anaesthesia_duration_mins
    ))
  
  asa_data <- anaesthesia_ax_clean |>
    filter(!is.na(asa_grade)) |>
    distinct(PAT_ENC_CSN_ID, .keep_all = TRUE) |>
    select(PAT_ENC_CSN_ID, asa_grade)
  
  df <- df |> left_join(asa_data, by = "PAT_ENC_CSN_ID")
  if (sum(is.na(df$asa_grade)) > 0) {
    asa_mode <- names(sort(table(df$asa_grade), decreasing = TRUE))[1]
    df <- df |> mutate(asa_grade = if_else(is.na(asa_grade), asa_mode, asa_grade))
  }
  df <- df |> mutate(asa_grade = factor(asa_grade, ordered = TRUE))
  
  df
}

post_match <- add_psm_covariates(post_match)

shared_proc_category_levels <- union(
  levels(factor(pre_match$proc_category)), levels(factor(post_match$proc_category))
)
shared_specialty_levels <- union(
  levels(factor(pre_match$specialty_collapsed)), levels(factor(post_match$specialty_collapsed))
)

pre_match <- pre_match |>
  mutate(
    proc_category = factor(proc_category, levels = shared_proc_category_levels),
    specialty_collapsed = factor(specialty_collapsed, levels = shared_specialty_levels)
  )
post_match <- post_match |>
  mutate(
    proc_category = factor(proc_category, levels = shared_proc_category_levels),
    specialty_collapsed = factor(specialty_collapsed, levels = shared_specialty_levels)
  )

balance_vars <- c("AGE", "WEIGHT_KG", "GENDER", "mode",
                  "proc_category", "specialty_collapsed", "asa_grade",
                  "had_peripheral", "year",
                  "had_remifentanil", "had_intraop_ketamine",
                  "had_intraop_clonidine", "neuraxial_intrathecal_mar",
                  "had_intraop_paracetamol", "had_intraop_nsaid",
                  "had_intraop_magnesium", "had_intraop_fentanyl",
                  "had_intraop_morphine", "had_intraop_oxycodone2",
                  "had_intraop_diamorphine", "anaesthesia_duration_mins")

tbl3_labels <- list(
  AGE ~ "Age, years", WEIGHT_KG ~ "Weight, kg", GENDER ~ "Sex",
  mode ~ "Mode", proc_category ~ "Procedure category",
  specialty_collapsed ~ "Specialty", asa_grade ~ "ASA grade",
  had_peripheral ~ "Peripheral nerve block", year ~ "Year",
  had_remifentanil ~ "Remifentanil", had_intraop_ketamine ~ "Intraop ketamine",
  had_intraop_clonidine ~ "Intraop clonidine",
  neuraxial_intrathecal_mar ~ "Intrathecal opioid",
  had_intraop_paracetamol ~ "Intraop paracetamol",
  had_intraop_nsaid ~ "Intraop NSAID", had_intraop_magnesium ~ "Intraop magnesium",
  had_intraop_fentanyl ~ "Intraop fentanyl", had_intraop_morphine ~ "Intraop morphine",
  had_intraop_oxycodone2 ~ "Intraop oxycodone", had_intraop_diamorphine ~ "Intraop diamorphine",
  anaesthesia_duration_mins ~ "Anaesthesia duration, min"
)

tbl3_pre <- pre_match |>
  select(drug, all_of(balance_vars)) |>
  mutate(across(where(is.logical), ~ factor(., labels = c("No", "Yes")))) |>
  tbl_summary(
    by = drug, label = tbl3_labels,
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1), missing = "no"
  ) |>
  modify_header(label = "**Variable**",
                stat_1 = sprintf("**Morphine (n=%d)**", sum(pre_match$drug == "Morphine")),
                stat_2 = sprintf("**Oxycodone (n=%d)**", sum(pre_match$drug == "Oxycodone"))) |>
  bold_labels()

tbl3_post <- post_match |>
  select(drug, all_of(balance_vars)) |>
  mutate(across(where(is.logical), ~ factor(., labels = c("No", "Yes")))) |>
  tbl_summary(
    by = drug, label = tbl3_labels,
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1), missing = "no"
  ) |>
  modify_header(label = "**Variable**",
                stat_1 = sprintf("**Morphine (n=%d)**", sum(post_match$drug == "Morphine")),
                stat_2 = sprintf("**Oxycodone (n=%d)**", sum(post_match$drug == "Oxycodone"))) |>
  bold_labels()

tbl3 <- tbl_merge(
  list(tbl3_pre, tbl3_post),
  tab_spanner = c("**Pre-match**", "**Post-match**")
) |>
  modify_caption("**Table 3.** Propensity score matching — covariate balance") |>
  modify_footnote(
    all_stat_cols() ~ "Continuous: median (IQR). Categorical: n (%). All post-match standardised mean differences <0.10. Concurrent epidural analgesia is not shown as a balance covariate: it is an eligibility criterion for this population (no episode with a concurrent epidural is included), not a variable requiring matching."
  )
cat("Table 3 built\n\n")

# ============================================================
# TABLE 4 — PRIMARY OUTCOMES
# (Composite side effect REMOVED — see header note. Reports
# antiemetic, antipruritic, and naloxone reversal independently.)
# ============================================================
cat("Building Table 4 — matched outcomes...\n")

table4_base <- matched_data |>
  mutate(
    drug = factor(if_else(treated == 1, "Oxycodone", "Morphine"),
                  levels = c("Morphine", "Oxycodone")),
    episode_duration_hrs = episode_duration_hrs.x
  ) |>
  select(-episode_duration_hrs.x, -episode_duration_hrs.y)

table4_matched_outcomes <- table4_base |>
  left_join(
    episode_metrics |>
      select(-any_of(setdiff(names(table4_base),
                             c("PAT_ENC_CSN_ID", "episode_id")))),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  ) |>
  mutate(
    across(c(reactive_antiemetic_final,
             any_antipruritic_final, any_naloxone_reversal_final),
           ~ factor(., levels = c(FALSE, TRUE), labels = c("No", "Yes")))
  )

episode_true_total_dose <- tryCatch(
  readRDS("1 - data/3 - processed_data/episode_true_total_dose.rds"),
  error = function(e) {
    cat("episode_true_total_dose.rds not found — Total dose rows will be skipped\n")
    NULL
  }
)

if (!is.null(episode_true_total_dose)) {
  table4_matched_outcomes <- table4_matched_outcomes |>
    left_join(
      episode_true_total_dose |> select(PAT_ENC_CSN_ID, episode_id, total_dose_mcg_kg_hr, background_mcg_kg_hr),
      by = c("PAT_ENC_CSN_ID", "episode_id")
    )
}

cat(sprintf("Peak pain score present: %d / %d matched episodes\n",
            sum(!is.na(table4_matched_outcomes$peak_pain_score)),
            nrow(table4_matched_outcomes)))
cat(sprintf("Time-weighted mean pain score present: %d / %d matched episodes\n\n",
            sum(!is.na(table4_matched_outcomes$time_weighted_mean_pain)),
            nrow(table4_matched_outcomes)))

tbl4 <- table4_matched_outcomes |>
  select(drug,
         reactive_antiemetic_final,
         any_antipruritic_final, any_naloxone_reversal_final,
         peak_pain_score, time_weighted_mean_pain,
         episode_duration_hrs,
         any_of(c("total_dose_mcg_kg_hr", "background_mcg_kg_hr"))) |>
  tbl_summary(
    by = drug,
    label = list(
      reactive_antiemetic_final ~ "Antiemetic (reactive, FINAL definition)",
      any_antipruritic_final ~ "Antipruritic",
      any_naloxone_reversal_final ~ "Naloxone reversal (respiratory)",
      peak_pain_score ~ "Peak pain score, median (IQR)",
      time_weighted_mean_pain ~ "Time-weighted mean pain score, median (IQR)",
      episode_duration_hrs ~ "Episode duration, hours [not a pre-defined outcome domain — see Methods]",
      total_dose_mcg_kg_hr ~ "Total dose (background + bolus), mcg/kg/hr [primary consumption outcome]",
      background_mcg_kg_hr ~ "  (reference) Background only, same linked subset, mcg/kg/hr"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}–{p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1),
    missing = "no"
  ) |>
  add_p(
    test = list(
      any_naloxone_reversal_final ~ "fisher.test",
      all_continuous() ~ "wilcox.test",
      all_categorical() ~ "chisq.test"
    ),
    pvalue_fun = ~ style_pvalue(., digits = 3)
  ) |>
  add_overall(last = FALSE) |>
  modify_header(
    label = "**Outcome**",
    stat_0 = sprintf("**Overall (n=%d)**", nrow(table4_matched_outcomes)),
    stat_1 = sprintf("**Morphine (n=%d)**", sum(table4_matched_outcomes$drug == "Morphine")),
    stat_2 = sprintf("**Oxycodone (n=%d)**", sum(table4_matched_outcomes$drug == "Oxycodone")),
    p.value = "**p-value**"
  ) |>
  modify_caption("**Table 4.** Primary outcomes — propensity-matched comparison") |>
  bold_labels() |>
  bold_p(t = 0.05) |>
  modify_footnote(
    all_stat_cols() ~ "Continuous: median (IQR). Categorical: n (%). p-values: chi-squared (categorical) unless noted, Wilcoxon rank-sum (continuous), Fisher's exact test for naloxone reversal (sparse event count). Three co-equal outcome domains are reported (side effects, efficacy, consumption); side-effect components (antiemetic, antipruritic, naloxone reversal) are reported individually rather than as a combined composite. Side-effect outcomes use FINAL anaesthesia-stop-anchored definitions (see Methods). Peak pain score and time-weighted mean pain score are the efficacy domain; the latter weights each reading by the time it applied until the next reading, excluding time before the first documented reading (see Methods). Episode duration is reported as a descriptive characteristic, not one of the three pre-defined outcome domains. Total dose (background + bolus) restricted to episodes with >=80% pump flowsheet coverage; linkage not random by drug (p=0.013) but a separate IPW sensitivity check confirmed the naive comparison is robust to this. NOTE: naloxone reversal counts shown here (matched subset, n=474) differ from the full-cohort chart-reviewed figure (0 oxycodone vs 3 morphine, all episodes) reported elsewhere in the manuscript — both are correct, describing different populations (matched vs full cohort); do not present them as conflicting."
  )

cat("--- Laxative (descriptive only, not part of tbl4's formal p-value table) ---\n")
laxative_descriptive <- table4_matched_outcomes |>
  group_by(drug) |>
  summarise(n = n(), pct_laxative = round(100 * mean(any_laxative_full, na.rm = TRUE), 1),
            .groups = "drop")
print(laxative_descriptive)
cat("(Not tested statistically. Report as descriptive text, not in Table 4's formal rows.)\n\n")

cat("Table 4 built\n\n")

# ============================================================
# TABLE 5 — ROTATION (descriptive panels + adjusted OR)
#
# Built by 15_rotation_analysis.R, run BEFORE this script — produces
# FOUR objects, not three: tbl5 (5a, descriptive rate), ft5b
# (destinations), ft5c (multi-rotation pattern), and ft5d (NEW —
# covariate-adjusted rotation-away OR, computed on the pre-match
# eligible pool due to too few events — 25 — in the 237/237 matched
# sample alone; see script 15 for full rationale).
# ============================================================

stopifnot(
  exists("tbl5"), exists("ft5b"), exists("ft5c"), exists("ft5d")
)
cat("Table 5 objects (tbl5, ft5b, ft5c, ft5d) confirmed present — built by\n")
cat("15_rotation_analysis.R, not constructed here.\n\n")



# ============================================================
# EXPORT TO WORD
# ============================================================
cat("Exporting to Word...\n")
output_path <- sprintf("manuscript_tables_R_%s.docx",
                       format(Sys.time(), "%Y%m%d_%H%M%S"))
cat("Output file:", output_path, "\n")
ft1 <- as_flex_table(tbl1) |> ft_theme()
ft2 <- as_flex_table(tbl2) |> ft_theme()
ft3 <- as_flex_table(tbl3) |> ft_theme()
ft4 <- as_flex_table(tbl4) |> ft_theme()
ft5a <- as_flex_table(tbl5) |> ft_theme()

# NOTE — the study end date in the manuscript title below still
# needs confirming (April vs July 2026 — see Methods discussion).
# Update once max(admissions_clean$ADM_DATE) is checked.
doc <- read_docx() |>
  body_add_par("Oxycodone PCA/NCA in Children: Manuscript Tables",
               style = "heading 1") |>
  body_add_par("Addenbrooke's Hospital · October 2014 – [CONFIRM DATE]",
               style = "Normal") |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft1) |>
  body_add_par("", style = "Normal") |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft2) |>
  body_add_par("", style = "Normal") |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft3) |>
  body_add_par("", style = "Normal") |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft4) |>
  body_add_par("", style = "Normal") |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft5a) |>
  body_add_par("", style = "Normal") |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft5b) |>
  body_add_par("", style = "Normal") |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft5c) |>
  body_add_par("", style = "Normal") |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft5d)

print(doc, target = output_path)

cat("\n============================================================\n")
cat("ANALYSIS 11 COMPLETE\n")
cat("Output:", output_path, "(Tables 1-5, including Panel 5d)\n")
cat("============================================================\n")