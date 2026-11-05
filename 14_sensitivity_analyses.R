# ============================================================
# 14 — Sensitivity Analyses
# Paediatric PCA Study
# Addenbrooke's Hospital
#
# run_psm()'s formula matches 09_comparative_psm.R's production
# formula exactly (specialty_collapsed, asa_grade, expanded intraop
# covariates), so every sensitivity check tests the same core model
# as the primary analysis rather than a simpler stand-in. Outcomes
# use the finalised anaesthesia-stop-anchored definitions (04_derive_episode_flags.R
# Flag 9). The primary row is read live from psm_matched_keys.rds
# rather than hardcoded.
#
# THIS PASS: every pool is built from the PRIMARY definition
# (06_classify_cohorts.R) and varies ONE thing at a time. Before,
# the regional-free, strict-exclusion and all-cases pools were built
# from separate filters that (a) did not apply the first-exposure rule
# (the child's first-ever recorded PCA/NCA episode) or the medical
# oncology exclusion, (b) used the coarser case_type instead of
# case_type_refined, and (c) did not require mode PCA/NCA.
# make_pool() below carries the shared definition, and a check
# confirms it reproduces the primary S3 pool exactly.
#
# CHANGED (Oct 2026): neonates are NO LONGER excluded, in this
# script or in 09_comparative_psm.R (a neonate without a comparator
# within the caliper simply drops out at matching). run_psm() also
# now applies the same patient-level deduplication (one episode per
# true patient per drug, earliest kept) that script 10 applies before
# matching, so the 1:2, 0.3 SD and other checks match on the same
# kind of pool as the primary.
#
# Requires, in the same session: 06_classify_cohorts.R (builds
# s3_*_keys and first_exposure_keys) and 09_comparative_psm.R
# (builds psm_matched_keys.rds).
#
# Regional-technique sensitivity is a single "regional-free" pool
# excluding all regional techniques (epidural, caudal, spinal,
# peripheral nerve block) together, rather than testing caudal and
# spinal exclusion separately — the clinically meaningful question is
# whether the result holds when PCA/NCA is the only analgesic
# technique. A separate "regional-inclusive" pool isn't a distinct
# check, since that's the existing primary population, which already
# excludes only definitive epidural and treats everything else as a
# covariate.
#
# "All cases (incl. medical)" is an exploratory comparison, not a
# same-population sensitivity check: S3/primary is surgical-only by
# construction, so this is a genuinely different population rather
# than a variation on one analytic choice.
#
# A patient-level deduplication check on the primary matched result
# specifically (not a re-match, but a check of whether the existing
# matched-pairs OR changes when repeat-patient episodes are
# deduplicated to one per patient per arm) is included because peak
# pain score was previously found to be a clustering artefact under
# this kind of check, and the composite side-effect outcome hadn't
# been checked for the same issue.
#
# Episode gap threshold (6h/24h vs the locked 12h) is not tested
# here — it requires a full pipeline rerun from 02_clean_data.R
# with a different gap value, not a filter on already-built episodes.
#
# Primary outcome throughout: any_side_effect_final
# ============================================================

library(tidyverse)
library(lubridate)
library(MatchIt)

# ============================================================
# LOAD DATA
# ============================================================

stopifnot(
  exists("pca_episodes_flagged"),
  exists("episode_metrics"),
  exists("intraop_mar_clean"),
  exists("anaesthesia_ax_clean"),
  exists("anaesthesia_events_clean"),
  exists("admissions_clean"),   # needed for PAT_ID, dedup
  exists("s3_morphine_keys"),
  exists("s3_oxycodone_keys")
)

if (!exists("first_exposure_keys")) {
  stop("first_exposure_keys not found in the session. Run 06_classify_cohorts.R ",
       "(it builds first_exposure_keys) in this session first.")
}

cat("============================================================\n")
cat("ANALYSIS 9 — SENSITIVITY ANALYSES\n")
cat("============================================================\n\n")

fmt_n_pct <- function(n, total) {
  sprintf("%d (%.1f%%)", n, 100 * n / total)
}

# ============================================================
# SHARED PRE-PROCESSING — same covariates as 09_comparative_psm.R
# ============================================================

build_matching_covariates <- function(data) {
  
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
  
  data <- data |>
    mutate(
      specialty_merged = recode(as.character(ADMITTING_SPECIALTY),
                                !!!same_service_map,
                                .default = as.character(ADMITTING_SPECIALTY))
    )
  
  keep_specialties <- data |>
    count(specialty_merged) |>
    filter(n >= SPECIALTY_VOLUME_THRESHOLD) |>
    pull(specialty_merged)
  
  data <- data |>
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
  
  data <- data |>
    left_join(new_intraop_flags, by = "PAT_ENC_CSN_ID") |>
    mutate(across(c(had_intraop_paracetamol, had_intraop_nsaid, had_intraop_magnesium,
                    had_intraop_fentanyl, had_intraop_morphine, had_intraop_oxycodone2,
                    had_intraop_diamorphine),
                  ~ replace_na(.x, FALSE)))
  
  duration_data <- anaesthesia_events_clean |>
    select(PAT_ENC_CSN_ID, anaesthesia_duration_mins) |>
    distinct(PAT_ENC_CSN_ID, .keep_all = TRUE)
  
  data <- data |> left_join(duration_data, by = "PAT_ENC_CSN_ID")
  median_duration <- median(data$anaesthesia_duration_mins, na.rm = TRUE)
  data <- data |>
    mutate(anaesthesia_duration_mins = if_else(
      is.na(anaesthesia_duration_mins), median_duration, anaesthesia_duration_mins
    ))
  
  asa_data <- anaesthesia_ax_clean |>
    filter(!is.na(asa_grade)) |>
    distinct(PAT_ENC_CSN_ID, .keep_all = TRUE) |>
    select(PAT_ENC_CSN_ID, asa_grade)
  
  data <- data |> left_join(asa_data, by = "PAT_ENC_CSN_ID")
  if (sum(is.na(data$asa_grade)) > 0) {
    asa_mode <- names(sort(table(data$asa_grade), decreasing = TRUE))[1]
    data <- data |> mutate(asa_grade = if_else(is.na(asa_grade), asa_mode, asa_grade))
  }
  data <- data |> mutate(asa_grade = factor(asa_grade, ordered = TRUE))
  
  data <- data |> filter(GENDER != "Not specified")
  
  data
}

PSM_FORMULA <- treated ~ AGE + WEIGHT_KG + GENDER + proc_category +
  mode + had_peripheral + year + specialty_collapsed + asa_grade +
  had_remifentanil + had_intraop_ketamine + had_intraop_clonidine +
  neuraxial_intrathecal_mar +
  had_intraop_paracetamol + had_intraop_nsaid + had_intraop_magnesium +
  had_intraop_fentanyl + had_intraop_morphine + had_intraop_oxycodone2 +
  had_intraop_diamorphine + anaesthesia_duration_mins

# ============================================================
# SHARED POOL DEFINITION — the primary S3 definition from
# 06_classify_cohorts.R, varying one choice at a time.
# Common to every pool: first line, first episode, no concurrent
# epidural, mode PCA/NCA, medical oncology excluded, and the child's
# first-ever recorded PCA/NCA episode (first_exposure_keys, built in
# script 6). Options change ONE thing:
#   surgical_only = FALSE : all cases incl. medical (exploratory)
#   strict        = TRUE  : strict instead of conservative minor-
#                           procedure exclusion set
#   regional_free = TRUE  : also exclude epidural, spinal, peripheral
#                           nerve block
# ============================================================

make_pool <- function(drug_name, treated_val,
                      surgical_only = TRUE, strict = FALSE, regional_free = FALSE) {
  d <- pca_episodes_flagged |>
    filter(
      drug == drug_name,
      is_first_line == TRUE,
      is_first_episode == TRUE,
      concurrent_epidural == FALSE,
      mode %in% c("PCA", "NCA"),
      !coalesce(str_detect(toupper(as.character(ADMITTING_SPECIALTY)), "MEDICAL ONCOLOGY"), FALSE)
    ) |>
    semi_join(first_exposure_keys, by = c("PAT_ENC_CSN_ID", "episode_id"))
  
  if (surgical_only) d <- d |> filter(case_type_refined == "surgical")
  
  if (strict) {
    d <- d |> filter(exclude_minor_procedure_any_strict == FALSE)
  } else {
    d <- d |> filter(exclude_minor_procedure_any == FALSE)
  }
  
  if (regional_free) {
    d <- d |> filter(
      had_lda_epidural            == FALSE,
      neuraxial_epidural_mar_only == FALSE,
      neuraxial_intrathecal_mar   == FALSE,
      had_peripheral              == FALSE
    )
  }
  
  d |> mutate(treated = treated_val, year = as.numeric(year(episode_start)))
}

# ============================================================
# Core PSM function
# ============================================================

run_psm <- function(data, label, ratio = 1, caliper = 0.2) {
  
  cat(sprintf("\n--- %s ---\n", label))
  
  # Neonates are NOT excluded (as in 09_comparative_psm.R); any without
  # a comparator within the caliper simply drop out at matching.
  cat(sprintf("Neonates in pool (not excluded): %d\n",
              sum(data$age_group == "Neonate (<28 days)", na.rm = TRUE)))
  
  # Patient-level deduplication, as in 09_comparative_psm.R: one episode
  # per true patient (PAT_ID) per drug, earliest kept, BEFORE matching.
  n_before_dedup <- nrow(data)
  data <- data |>
    left_join(admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID),
              by = "PAT_ENC_CSN_ID") |>
    arrange(PAT_ID, episode_start) |>
    distinct(PAT_ID, drug, .keep_all = TRUE)
  cat(sprintf("Patient-level dedup: %d -> %d episodes (%d duplicate episodes removed)\n",
              n_before_dedup, nrow(data), n_before_dedup - nrow(data)))
  
  cat(sprintf("Pre-match: morphine n=%d, oxycodone n=%d\n",
              sum(data$treated == 0), sum(data$treated == 1)))
  
  data <- data |>
    mutate(WEIGHT_KG = if_else(is.na(WEIGHT_KG), median(WEIGHT_KG, na.rm = TRUE), WEIGHT_KG))
  
  data <- build_matching_covariates(data)
  
  data <- data |>
    mutate(
      proc_category        = factor(proc_category),
      specialty_collapsed  = factor(specialty_collapsed),
      GENDER               = factor(GENDER),
      mode                 = factor(mode)
    )
  
  cat(sprintf("Matching pool after GENDER fix: morphine n=%d, oxycodone n=%d\n",
              sum(data$treated == 0), sum(data$treated == 1)))
  
  set.seed(42)
  match_out <- tryCatch({
    matchit(
      formula     = PSM_FORMULA,
      data        = data,
      method      = "nearest",
      distance    = "glm",
      link        = "logit",
      ratio       = ratio,
      caliper     = caliper,
      std.caliper = TRUE
    )
  }, error = function(e) {
    cat("MatchIt error:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(match_out)) return(NULL)
  
  matched <- match.data(match_out)
  
  cat(sprintf("Matched: oxycodone n=%d, morphine n=%d\n",
              sum(matched$treated == 1), sum(matched$treated == 0)))
  cat(sprintf("Unmatched oxycodone: %d\n",
              sum(data$treated == 1) - sum(matched$treated == 1)))
  
  if ("age_group" %in% names(matched)) {
    n_neo_matched <- sum(matched$age_group == "Neonate (<28 days)", na.rm = TRUE)
    cat(sprintf("Neonate episodes in matched set: %d\n", n_neo_matched))
  }
  
  outcomes <- matched |>
    select(PAT_ENC_CSN_ID, episode_id, drug, treated) |>
    left_join(
      pca_episodes_flagged |>
        select(PAT_ENC_CSN_ID, episode_id,
               reactive_antiemetic_final, any_antipruritic_final,
               any_naloxone_reversal_final, any_side_effect_final),
      by = c("PAT_ENC_CSN_ID", "episode_id")
    )
  
  result <- outcomes |>
    group_by(drug) |>
    summarise(
      n               = n(),
      pct_side_effect = round(100 * mean(any_side_effect_final, na.rm = TRUE), 1),
      pct_antiemetic  = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
      pct_antiprurit  = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
      n_nalox_rev     = sum(any_naloxone_reversal_final, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(result)
  
  tab <- table(outcomes$drug, outcomes$any_side_effect_final)
  if (ncol(tab) == 2) {
    ct <- suppressWarnings(chisq.test(tab))
    or <- (tab["oxycodone", 2] / tab["oxycodone", 1]) /
      (tab["morphine",  2] / tab["morphine",  1])
    cat(sprintf("Composite SE: OR=%.2f, p=%.4f\n", or, ct$p.value))
  }
  
  cat("\n")
  invisible(outcomes)
}

summarise_sensitivity <- function(outcomes, label) {
  if (is.null(outcomes)) {
    cat(sprintf("%-45s  FAILED\n", label))
    return(invisible(NULL))
  }
  mor <- outcomes |> filter(drug == "morphine")
  oxy <- outcomes |> filter(drug == "oxycodone")
  n_mor <- nrow(mor)
  n_oxy <- nrow(oxy)
  pct_mor <- round(100 * mean(mor$any_side_effect_final, na.rm = TRUE), 1)
  pct_oxy <- round(100 * mean(oxy$any_side_effect_final, na.rm = TRUE), 1)
  tab <- table(outcomes$drug, outcomes$any_side_effect_final)
  or  <- NA_real_
  pv  <- NA_real_
  if (ncol(tab) == 2) {
    ct <- suppressWarnings(chisq.test(tab))
    or <- round((tab["oxycodone",2]/tab["oxycodone",1]) /
                  (tab["morphine", 2]/tab["morphine", 1]), 2)
    pv <- round(ct$p.value, 4)
  }
  cat(sprintf("%-45s  n=%d, %s%%  n=%d, %s%%  OR=%.2f  p=%.4f\n",
              label, n_mor, pct_mor, n_oxy, pct_oxy,
              ifelse(is.na(or), 0, or),
              ifelse(is.na(pv), 1, pv)))
}

# ============================================================
# BASE POOL — primary population (same as 09_comparative_psm.R)
# ============================================================

base_morphine <- pca_episodes_flagged |>
  semi_join(s3_morphine_keys, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(treated = 0L, year = as.numeric(year(episode_start)))

base_oxycodone <- pca_episodes_flagged |>
  semi_join(s3_oxycodone_keys, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(treated = 1L, year = as.numeric(year(episode_start)))

base_pool <- bind_rows(base_morphine, base_oxycodone)

cat("Base pool (primary analysis, from the S3 keys, before patient-level dedup):\n")
cat(sprintf("  Morphine:  n=%d\n", nrow(base_morphine)))
cat(sprintf("  Oxycodone: n=%d\n", nrow(base_oxycodone)))

# Check: make_pool() with default options must reproduce the primary S3
# pool, otherwise the sensitivity pools below are not "primary + one
# change". Any difference is flagged loudly.
chk_m <- make_pool("morphine", 0L)
chk_o <- make_pool("oxycodone", 1L)
same_m <- nrow(chk_m) == nrow(base_morphine) &&
  nrow(anti_join(chk_m, base_morphine, by = c("PAT_ENC_CSN_ID", "episode_id"))) == 0
same_o <- nrow(chk_o) == nrow(base_oxycodone) &&
  nrow(anti_join(chk_o, base_oxycodone, by = c("PAT_ENC_CSN_ID", "episode_id"))) == 0
cat(sprintf("make_pool() default vs primary S3 pool — morphine: %d vs %d | oxycodone: %d vs %d\n",
            nrow(chk_m), nrow(base_morphine), nrow(chk_o), nrow(base_oxycodone)))
if (same_m && same_o) {
  cat("  Identical: sensitivity pools are primary + one change.\n\n")
} else {
  warning("make_pool() does NOT reproduce the primary S3 pool. The filters in ",
          "make_pool() have drifted from 06_classify_cohorts.R; compare them ",
          "before trusting S3, S4 and the exploratory rows.")
  cat("  WARNING: not identical — see the warning at the end of the run.\n\n")
}

# ============================================================
# S1 — PSM RATIO 1:2
# ============================================================

cat("============================================================\n")
cat("S1 — PSM RATIO 1:2 (primary uses 1:1)\n")
cat("============================================================\n")
s1_outcomes <- run_psm(base_pool, "PSM 1:2", ratio = 2)

# ============================================================
# S2 — BROADER CALIPER (0.3 SD)
# ============================================================

cat("============================================================\n")
cat("S2 — BROADER CALIPER 0.3 SD (primary uses 0.2 SD)\n")
cat("============================================================\n")
s2_outcomes <- run_psm(base_pool, "PSM caliper 0.3 SD", ratio = 1, caliper = 0.3)

# ============================================================
# S3 — REGIONAL-FREE POOL
# Excludes ALL regional techniques together (epidural, caudal,
# spinal, peripheral nerve block) — replaces old separate
# caudal-only / spinal-only checks with one coherent "PCA/NCA was
# the ONLY analgesic technique" comparison.
# ============================================================

cat("============================================================\n")
cat("S3 — REGIONAL-FREE POOL\n")
cat("(Excludes epidural, caudal, spinal AND peripheral nerve block —\n")
cat(" primary only excludes definitive epidural, keeps others as\n")
cat(" covariates)\n")
cat("============================================================\n")

morphine_regional_free  <- make_pool("morphine",  0L, regional_free = TRUE)
oxycodone_regional_free <- make_pool("oxycodone", 1L, regional_free = TRUE)

pool_regional_free <- bind_rows(morphine_regional_free, oxycodone_regional_free)

cat(sprintf("Regional-free pool: morphine n=%d, oxycodone n=%d\n",
            nrow(morphine_regional_free), nrow(oxycodone_regional_free)))

s3_outcomes <- run_psm(pool_regional_free, "PSM regional-free", ratio = 1)

# ============================================================
# S4 — PREDEFINED EXCLUSIONS: STRICT vs CONSERVATIVE
# Both already computed in 04_derive_episode_flags.R; never previously compared here.
# ============================================================

cat("============================================================\n")
cat("S4 — PREDEFINED EXCLUSIONS: STRICT (vs primary's CONSERVATIVE)\n")
cat("============================================================\n")

morphine_strict  <- make_pool("morphine",  0L, strict = TRUE)
oxycodone_strict <- make_pool("oxycodone", 1L, strict = TRUE)

pool_strict <- bind_rows(morphine_strict, oxycodone_strict)

cat(sprintf("Strict-exclusion pool: morphine n=%d, oxycodone n=%d\n",
            nrow(morphine_strict), nrow(oxycodone_strict)))

s4_outcomes <- run_psm(pool_strict, "PSM strict exclusion set", ratio = 1)

# ============================================================
# EXPLORATORY (not a same-population sensitivity check) —
# ALL CASES INCLUDING MEDICAL. S3/primary is surgical-only BY
# CONSTRUCTION, so this tests a genuinely different population,
# not one varied analytic choice.
# ============================================================

cat("============================================================\n")
cat("EXPLORATORY — ALL CASES INCLUDING MEDICAL\n")
cat("(NOT a same-population sensitivity check — primary population\n")
cat(" is surgical-only by construction, this tests a genuinely\n")
cat(" different, broader population)\n")
cat("============================================================\n")

morphine_all <- make_pool("morphine", 0L, surgical_only = FALSE) |>
  mutate(proc_category = replace_na(proc_category, "Medical/none"))

oxycodone_all <- make_pool("oxycodone", 1L, surgical_only = FALSE) |>
  mutate(proc_category = replace_na(proc_category, "Medical/none"))

pool_all_cases <- bind_rows(morphine_all, oxycodone_all)

cat(sprintf("All cases pool: morphine n=%d, oxycodone n=%d\n",
            nrow(morphine_all), nrow(oxycodone_all)))

exploratory_outcomes <- run_psm(pool_all_cases, "EXPLORATORY: all cases (incl medical)", ratio = 1)

# ============================================================
# S5 — PATIENT-LEVEL DEDUPLICATION CHECK ON THE PRIMARY RESULT
# Not a re-match — checks whether the PRIMARY matched-pairs OR
# changes when repeat-patient episodes are deduplicated to one per
# patient per arm. Peak pain score was already found to be a
# clustering artefact once during validation; composite side effect has
# never been checked for the same issue.
# ============================================================

cat("============================================================\n")
cat("S5 — PATIENT-LEVEL DEDUPLICATION CHECK (on PRIMARY matched pairs)\n")
cat("============================================================\n")

primary_matched_keys <- tryCatch(
  readRDS("1 - data/3 - processed_data/psm_matched_keys.rds"),
  error = function(e) {
    cat("psm_matched_keys.rds not found — run 09_comparative_psm.R first.\n")
    return(NULL)
  }
)

if (!is.null(primary_matched_keys)) {
  
  primary_matched_full <- primary_matched_keys |>
    select(PAT_ENC_CSN_ID, episode_id, drug) |>
    left_join(admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID),
              by = "PAT_ENC_CSN_ID") |>
    left_join(
      pca_episodes_flagged |> select(PAT_ENC_CSN_ID, episode_id, episode_start, any_side_effect_final),
      by = c("PAT_ENC_CSN_ID", "episode_id")
    )
  
  n_dup_patients <- primary_matched_full |>
    count(PAT_ID) |>
    filter(n > 1) |>
    nrow()
  
  cat(sprintf("Patients appearing more than once in matched pairs: %d\n", n_dup_patients))
  
  # Deduplicate: one episode per patient PER ARM (a patient could
  # legitimately appear once in each arm if they had both drugs
  # matched at different times — only collapse WITHIN each arm)
  primary_matched_dedup <- primary_matched_full |>
    arrange(PAT_ID, drug, episode_start) |>
    distinct(PAT_ID, drug, .keep_all = TRUE)
  
  cat(sprintf("Episode-level n: morphine=%d, oxycodone=%d\n",
              sum(primary_matched_full$drug == "morphine"),
              sum(primary_matched_full$drug == "oxycodone")))
  cat(sprintf("Patient-deduplicated n: morphine=%d, oxycodone=%d\n\n",
              sum(primary_matched_dedup$drug == "morphine"),
              sum(primary_matched_dedup$drug == "oxycodone")))
  
  cat("--- Episode-level (primary, as-is) ---\n")
  summarise_sensitivity(primary_matched_full, "Episode-level (primary)")
  
  cat("\n--- Patient-deduplicated ---\n")
  summarise_sensitivity(primary_matched_dedup, "Patient-deduplicated")
  
} else {
  cat("Dedup check skipped — psm_matched_keys.rds not available.\n\n")
}

# ============================================================
# SENSITIVITY SUMMARY TABLE
# ============================================================

cat("\n============================================================\n")
cat("SENSITIVITY SUMMARY\n")
cat("============================================================\n\n")

cat(sprintf("%-45s  %-20s  %-20s  %-8s  %s\n",
            "Analysis", "Morphine", "Oxycodone", "OR", "p"))
cat(strrep("-", 105), "\n")

if (!is.null(primary_matched_keys)) {
  summarise_sensitivity(primary_matched_full, "PRIMARY (1:1, surgical, 0.2SD, live)")
} else {
  cat(sprintf("%-45s  SKIPPED (see warning above)\n", "PRIMARY"))
}

summarise_sensitivity(s1_outcomes, "S1: PSM 1:2")
summarise_sensitivity(s2_outcomes, "S2: Caliper 0.3 SD")
summarise_sensitivity(s3_outcomes, "S3: Regional-free (no epi/caudal/spinal/PNB)")
summarise_sensitivity(s4_outcomes, "S4: strict exclusion set")
summarise_sensitivity(exploratory_outcomes, "EXPLORATORY: all cases (incl medical)")
if (!is.null(primary_matched_keys)) {
  summarise_sensitivity(primary_matched_dedup, "S5: Patient-deduplicated (primary pairs)")
}

cat("\n")
cat("============================================================\n")
cat("ANALYSIS 9 COMPLETE\n")
cat("All sensitivity checks now use the SAME core matching formula\n")
cat("as 09_comparative_psm.R (specialty_collapsed, asa_grade, expanded\n")
cat("intraop covariates) and the FINALISED side-effect definitions, and\n")
cat("the same pool definition as the primary (first-ever PCA/NCA episode,\n")
cat("medical oncology excluded, neonates NOT excluded, patient-level\n")
cat("dedup applied) varied one choice at a time.\n")
cat("Compare each row against the LIVE primary row above — do not\n")
cat("assume consistency with any previously-reported figure.\n\n")
cat("STILL DEFERRED: episode gap threshold sensitivity (6h/24h vs\n")
cat("locked 12h) — requires a full pipeline rerun from\n")
cat("02_clean_data.R with a different gap value, not just a filter\n")
cat("on already-built episodes.\n")
cat("============================================================\n")