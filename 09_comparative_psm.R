# ============================================================
# 09 - Propensity score matching: oxycodone vs morphine
# ============================================================
# Builds the matching pool from the primary comparative population
# (s3_*_keys from 06_classify_cohorts.R), runs 1:1 nearest-neighbour
# propensity score matching with MatchIt, and saves the matched keys
# (psm_matched_keys.rds). psm_data, matched_data, match_out and
# psm_formula_expanded are left in the session for the scripts that
# follow; matched outcomes are analysed in those scripts.
#
# Matching
#   Method:    nearest neighbour, 1:1, without replacement
#   Distance:  logistic regression propensity score (glm, logit link)
#   Caliper:   0.2 standard deviations of the logit propensity score
#   Seed:      42
#
# Covariates (21): age, weight, sex, procedure category, admitting
# specialty, delivery mode, peripheral nerve block, year, ASA grade,
# anaesthesia duration, and intraoperative remifentanil, ketamine,
# clonidine, intrathecal route (spinal), paracetamol, NSAID, magnesium,
# fentanyl, morphine, oxycodone and diamorphine. The intraoperative drugs
# are proxies for case complexity (confounding by indication), not
# treatment effects.
#
#   - Procedure urgency is not a covariate: documentation is sparse and
#     declines over the study period.
#   - Admitting specialty is merged where adult and paediatric labels
#     name the same service, and specialties with fewer than 30
#     episodes in the pool are grouped as "Other specialty", to avoid
#     sparse cells against procedure category.
#   - Sex "Not specified" is dropped, and procedure category "Other" is
#     merged into "Unclassified"; both are near-zero-variance levels
#     that cause separation in the propensity model.
#   - Missing weight, anaesthesia duration and ASA grade are imputed with
#     the median, median and most common grade respectively (counts are
#     printed below).
#   - Neonates are not excluded; a neonate with no comparator within the
#     caliper is simply left unmatched.
#   - The pool is deduplicated to one episode per child (earliest) before
#     matching. After the first-exposure restriction in
#     06_classify_cohorts.R this should remove nothing; it is kept as a
#     safeguard.
# ============================================================

library(tidyverse)
library(lubridate)
library(MatchIt)

stopifnot(
  exists("pca_episodes_flagged"),
  exists("episode_metrics"),
  exists("s3_morphine_keys"),
  exists("s3_oxycodone_keys"),
  exists("admissions_clean"),
  exists("intraop_mar_clean"),
  exists("anaesthesia_events_clean"),
  exists("anaesthesia_ax_clean")
)

PS_CALIPER <- 0.2

# ============================================================
# BUILD PRE-MATCH DATASET
# ============================================================

cat("--- Building pre-match dataset ---\n\n")

morphine_pool <- pca_episodes_flagged |>
  semi_join(s3_morphine_keys, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(treated = 0L)

oxycodone_pool <- pca_episodes_flagged |>
  semi_join(s3_oxycodone_keys, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(treated = 1L)

cat("Eligible pool, before deduplication and the sex fix:\n")
cat(sprintf("  Morphine:   n = %d\n", nrow(morphine_pool)))
cat(sprintf("  Oxycodone:  n = %d\n", nrow(oxycodone_pool)))
cat("\n")

# ============================================================
# Neonates (not excluded; counted only)
# ============================================================
# The number matched is printed after matching.

if ("age_group" %in% names(morphine_pool)) {
  n_neo_mor <- sum(morphine_pool$age_group == "Neonate (<28 days)", na.rm = TRUE)
  n_neo_oxy <- sum(oxycodone_pool$age_group == "Neonate (<28 days)", na.rm = TRUE)
  cat("--- Neonates (not excluded) ---\n")
  cat("Neonates in eligible pool: morphine =", n_neo_mor,
      "| oxycodone =", n_neo_oxy, "\n\n")
} else {
  cat("age_group column not found; cannot count neonates.\n\n")
}

# ============================================================
# One episode per child
# ============================================================
# Keeps each child's earliest episode within each drug pool, so that a
# child cannot contribute more than one row to the matching pool.

morphine_pool_dedup <- morphine_pool |>
  left_join(admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID), by = "PAT_ENC_CSN_ID") |>
  arrange(PAT_ID, episode_start) |>
  distinct(PAT_ID, .keep_all = TRUE)

oxycodone_pool_dedup <- oxycodone_pool |>
  left_join(admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID), by = "PAT_ENC_CSN_ID") |>
  arrange(PAT_ID, episode_start) |>
  distinct(PAT_ID, .keep_all = TRUE)

cat("--- Patient-level deduplication ---\n")
cat(sprintf("Morphine:  %d episodes -> %d unique patients (%d duplicate episodes removed)\n",
            nrow(morphine_pool), nrow(morphine_pool_dedup),
            nrow(morphine_pool) - nrow(morphine_pool_dedup)))
cat(sprintf("Oxycodone: %d episodes -> %d unique patients (%d duplicate episodes removed)\n\n",
            nrow(oxycodone_pool), nrow(oxycodone_pool_dedup),
            nrow(oxycodone_pool) - nrow(oxycodone_pool_dedup)))

morphine_pool  <- morphine_pool_dedup
oxycodone_pool <- oxycodone_pool_dedup

psm_data <- bind_rows(morphine_pool, oxycodone_pool) |>
  mutate(year = as.numeric(year(episode_start)))

cat("Total pre-match pool (post-dedup):", nrow(psm_data), "\n\n")

# Medical oncology is excluded in 06_classify_cohorts.R; stop if any
# remain (06 was not re-run)
if (any(str_detect(toupper(as.character(psm_data$ADMITTING_SPECIALTY)), "MEDICAL ONCOLOGY"),
        na.rm = TRUE)) {
  stop("Medical oncology episodes are still in the matching pool. Re-run ",
       "06_classify_cohorts.R (it now excludes them), then re-run this script.")
}

# ============================================================
# Drop sex "Not specified" (near-zero-variance level)
# ============================================================

n_not_specified <- sum(psm_data$GENDER == "Not specified", na.rm = TRUE)
cat("Episodes with GENDER == 'Not specified':", n_not_specified,
    "out of", nrow(psm_data), "\n")

if (n_not_specified > 0) {
  psm_data <- psm_data |> filter(GENDER != "Not specified")
}

cat("Pre-match pool after the sex fix:", nrow(psm_data), "\n")
cat(sprintf("  Morphine:   n = %d\n", sum(psm_data$treated == 0)))
cat(sprintf("  Oxycodone:  n = %d\n\n", sum(psm_data$treated == 1)))

# ============================================================
# ADMITTING_SPECIALTY -> specialty_collapsed
# ============================================================


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

psm_data <- psm_data |>
  mutate(
    specialty_merged = recode(as.character(ADMITTING_SPECIALTY),
                              !!!same_service_map,
                              .default = as.character(ADMITTING_SPECIALTY))
  )

keep_specialties <- psm_data |>
  count(specialty_merged) |>
  filter(n >= SPECIALTY_VOLUME_THRESHOLD) |>
  pull(specialty_merged)

psm_data <- psm_data |>
  mutate(
    specialty_collapsed = if_else(specialty_merged %in% keep_specialties,
                                  specialty_merged, "Other specialty")
  )

cat("ADMITTING_SPECIALTY levels: raw =", n_distinct(psm_data$ADMITTING_SPECIALTY),
    "| collapsed =", n_distinct(psm_data$specialty_collapsed), "\n\n")

# ============================================================
# Intraoperative covariate flags
# ============================================================

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

psm_data <- psm_data |>
  left_join(new_intraop_flags, by = "PAT_ENC_CSN_ID") |>
  mutate(across(c(had_intraop_paracetamol, had_intraop_nsaid, had_intraop_magnesium,
                  had_intraop_fentanyl, had_intraop_morphine, had_intraop_oxycodone2,
                  had_intraop_diamorphine),
                ~ replace_na(.x, FALSE)))


# ============================================================
# Anaesthesia duration
# ============================================================

duration_data <- anaesthesia_events_clean |>
  select(PAT_ENC_CSN_ID, anaesthesia_duration_mins) |>
  distinct(PAT_ENC_CSN_ID, .keep_all = TRUE)

psm_data <- psm_data |> left_join(duration_data, by = "PAT_ENC_CSN_ID")

n_missing_duration <- sum(is.na(psm_data$anaesthesia_duration_mins))
cat("anaesthesia_duration_mins missing:", n_missing_duration, "/", nrow(psm_data), "\n")

if (n_missing_duration > 0) {
  median_duration <- median(psm_data$anaesthesia_duration_mins, na.rm = TRUE)
  psm_data <- psm_data |>
    mutate(anaesthesia_duration_mins = if_else(
      is.na(anaesthesia_duration_mins), median_duration, anaesthesia_duration_mins
    ))
  cat("Imputed", n_missing_duration, "missing values with the median (",
      round(median_duration, 1), "min)\n\n")
}

# ============================================================
# ASA grade
# ============================================================
# From anaesthesia_ax_clean; the first non-missing grade is used if an
# admission has more than one assessment row.

asa_data <- anaesthesia_ax_clean |>
  filter(!is.na(asa_grade)) |>
  distinct(PAT_ENC_CSN_ID, .keep_all = TRUE) |>
  select(PAT_ENC_CSN_ID, asa_grade)

psm_data <- psm_data |> left_join(asa_data, by = "PAT_ENC_CSN_ID")

n_missing_asa <- sum(is.na(psm_data$asa_grade))
cat("asa_grade missing:", n_missing_asa, "/", nrow(psm_data),
    sprintf("(%.1f%%)\n", 100 * n_missing_asa / nrow(psm_data)))

if (n_missing_asa > 0) {
  cat("Imputing missing ASA grade with the most common grade\n\n")
  asa_mode <- names(sort(table(psm_data$asa_grade), decreasing = TRUE))[1]
  psm_data <- psm_data |>
    mutate(asa_grade = if_else(is.na(asa_grade), asa_mode, asa_grade))
}

psm_data <- psm_data |> mutate(asa_grade = factor(asa_grade, ordered = TRUE))

cat("ASA grade distribution (pre-match):\n")
print(table(psm_data$drug, psm_data$asa_grade))
cat("\n")

# ============================================================
# Check the matching variables
# ============================================================

matching_vars <- c(
  "AGE", "WEIGHT_KG", "GENDER", "proc_category", "mode",
  "had_peripheral", "year", "specialty_collapsed", "asa_grade",
  "had_remifentanil", "had_intraop_ketamine", "had_intraop_clonidine",
  "neuraxial_intrathecal_mar",
  "had_intraop_paracetamol", "had_intraop_nsaid", "had_intraop_magnesium",
  "had_intraop_fentanyl", "had_intraop_morphine", "had_intraop_oxycodone2",
  "had_intraop_diamorphine", "anaesthesia_duration_mins"
)

missing_vars <- matching_vars[!matching_vars %in% names(psm_data)]
if (length(missing_vars) > 0) {
  cat("Missing variables:", paste(missing_vars, collapse = ", "), "\n")
} else {
  cat("All matching variables present\n")
}
cat("\n")

cat("--- Missingness in matching variables ---\n")
psm_data |>
  summarise(across(all_of(matching_vars), ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0) |>
  print()
cat("\n")

# ============================================================
# Merge procedure category "Other" into "Unclassified"
# ============================================================

n_other <- sum(psm_data$proc_category == "Other", na.rm = TRUE)
cat("Episodes with proc_category == 'Other':", n_other, "\n")

psm_data <- psm_data |>
  mutate(
    proc_category = if_else(proc_category == "Other", "Unclassified", proc_category)
  )

cat("Procedure category levels after merge:",
    n_distinct(psm_data$proc_category), "\n\n")

psm_data <- psm_data |>
  mutate(
    proc_category        = factor(proc_category),
    ADMITTING_SPECIALTY  = factor(ADMITTING_SPECIALTY),
    specialty_collapsed  = factor(specialty_collapsed),
    GENDER               = factor(GENDER),
    mode                 = factor(mode)
  )

psm_data <- psm_data |>
  mutate(WEIGHT_KG = if_else(is.na(WEIGHT_KG), median(WEIGHT_KG, na.rm = TRUE), WEIGHT_KG))

cat("WEIGHT_KG missing after imputation:", sum(is.na(psm_data$WEIGHT_KG)), "\n\n")

psm_formula_expanded <- treated ~ AGE + WEIGHT_KG + GENDER + proc_category +
  mode + had_peripheral + year + specialty_collapsed + asa_grade +
  had_remifentanil + had_intraop_ketamine + had_intraop_clonidine +
  neuraxial_intrathecal_mar +
  had_intraop_paracetamol + had_intraop_nsaid + had_intraop_magnesium +
  had_intraop_fentanyl + had_intraop_morphine + had_intraop_oxycodone2 +
  had_intraop_diamorphine + anaesthesia_duration_mins

# ============================================================
# Propensity model convergence and separation check
# ============================================================

glm_warnings <- character(0)
glm_check <- withCallingHandlers(
  glm(psm_formula_expanded, data = psm_data, family = binomial),
  warning = function(w) {
    glm_warnings <<- c(glm_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(glm_warnings) > 0) {
  cat("GLM warnings raised during propensity model fit:\n")
  for (w in unique(glm_warnings)) cat("   -", w, "\n")
  cat("\n")
} else {
  cat("No GLM warnings; model converged without apparent separation.\n\n")
}

coef_table <- summary(glm_check)$coefficients
large_coefs <- coef_table[abs(coef_table[, "Estimate"]) > 5 | coef_table[, "Std. Error"] > 5, , drop = FALSE]

if (nrow(large_coefs) > 0) {
  cat("Terms with |estimate| > 5 or SE > 5 (possible separation indicators):\n")
  print(large_coefs)
  cat("\n")
} else {
  cat("No individual coefficients with |estimate| or SE > 5.\n\n")
}


# ============================================================
# Matching
# ============================================================

psm_formula <- psm_formula_expanded

cat("Running MatchIt...\n")
set.seed(42)

match_out <- tryCatch({
  matchit(
    formula  = psm_formula,
    data     = psm_data,
    method   = "nearest",
    distance = "glm",
    link     = "logit",
    ratio    = 1,
    caliper  = PS_CALIPER,
    std.caliper = TRUE
  )
}, error = function(e) {
  cat("MatchIt error:", conditionMessage(e), "\n")
  NULL
})

if (is.null(match_out)) {
  cat("Matching failed: check variable missingness and factor levels above\n")
  stop("PSM failed")
}

cat("Matching complete\n\n")
print(summary(match_out))

# ============================================================
# EXTRACT MATCHED DATASET
# ============================================================

matched_data <- match.data(match_out)

# Denominators are taken from psm_data, the actual matching pool
n_pool_oxy     <- sum(psm_data$treated == 1)
n_pool_mor     <- sum(psm_data$treated == 0)
n_matched_oxy  <- sum(matched_data$treated == 1)
n_matched_mor  <- sum(matched_data$treated == 0)
n_unmatched    <- n_pool_oxy - n_matched_oxy

cat(sprintf("\nMatched pairs:          %d\n", n_matched_oxy))
cat(sprintf("Oxycodone matched:      %d / %d\n", n_matched_oxy, n_pool_oxy))
cat(sprintf("Morphine matched:       %d / %d\n", n_matched_mor, n_pool_mor))
cat(sprintf("Oxycodone unmatched:    %d\n\n", n_unmatched))

# Neonates in the pool and matched
cat("--- Neonates in matching ---\n")
if ("age_group" %in% names(psm_data)) {
  neo_pre <- psm_data |>
    filter(age_group == "Neonate (<28 days)") |>
    count(drug, name = "n_in_pool")
  neo_post <- matched_data |>
    filter(age_group == "Neonate (<28 days)") |>
    count(drug, name = "n_matched")
  print(full_join(neo_pre, neo_post, by = "drug"))
  cat("\n")
} else {
  cat("age_group not available; cannot report neonates.\n\n")
}

matched_keys <- matched_data |>
  select(PAT_ENC_CSN_ID, episode_id, drug, treated, subclass, distance, weights)

saveRDS(matched_keys, "1 - data/3 - processed_data/psm_matched_keys.rds")
cat("Matched keys saved to psm_matched_keys.rds\n\n")

# ============================================================
# POST-MATCH BALANCE
# ============================================================

smd_summary <- summary(match_out, standardize = TRUE)
print(smd_summary$sum.matched)
cat("\nMaximum post-match |SMD|:",
    round(max(abs(smd_summary$sum.matched[, "Std. Mean Diff."]), na.rm = TRUE), 3), "\n")
cat("Terms with |SMD| > 0.10:",
    sum(abs(smd_summary$sum.matched[, "Std. Mean Diff."]) > 0.10, na.rm = TRUE),
    "of", nrow(smd_summary$sum.matched), "\n")
