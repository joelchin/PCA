# ============================================================
# 10 — Section 3: PSM Comparative Analysis
# Paediatric PCA Study
# Addenbrooke's Hospital
#
# Propensity-matched comparison of oxycodone vs morphine, addressing
# indication bias in the unadjusted comparison from Section 2.
#
# Matching covariates
#
#   The matching formula includes ASA grade as a pre-treatment
#   covariate. This population (S3) is 100% surgical by construction,
#   so ASA grade needs no "not applicable" coding for medical
#   patients the way it would in the full S4 cohort. It's sourced
#   from anaesthesia_ax_clean (built in 2_cleaningscript.R, Section
#   7B) and joined by PAT_ENC_CSN_ID.
#
#   Procedure urgency is deliberately not included as a covariate.
#   Coverage is too sparse (~17.6% of assessment-SDE rows, and
#   declining over the dataset's history) to support inclusion in a
#   matching formula, so it's retained only as a descriptive field.
#
#   ADMITTING_SPECIALTY (28 levels) is collapsed to specialty_collapsed
#   (9 levels): apparent same-service adult/paediatric label pairs are
#   merged, and remaining specialties with n<30 are grouped into
#   "Other specialty". This avoids a sparse-cell problem against
#   proc_category that would otherwise risk unstable or extreme
#   propensity scores.
#
#   Intraoperative covariates (had_remifentanil/ketamine/clonidine/
#   magnesium/nsaid/paracetamol/fentanyl/morphine/oxycodone/
#   diamorphine) are included based on pre-match imbalance checks:
#   several show large imbalances between drug groups (up to +41.6pp,
#   p<2e-16), consistent with confounding by indication rather than
#   purely theoretical justification. had_intraop_alfentanil is the
#   one exception, excluded for showing no imbalance (+0.7pp, p=0.41).
#
#   GENDER == "Not specified" (a data-entry gap) and
#   proc_category == "Other" (collapsed into "Unclassified") are
#   dropped/regrouped as near-zero-variance levels that otherwise
#   produce separation in the propensity model (standard errors in
#   the hundreds for both).
#
#   A drug x era interaction term on the primary outcome is not
#   significant (p=0.55), so era is not treated as an effect modifier
#   and the pooled model with year as a covariate is used throughout.
#
# Side-effect outcome definitions
#
#   Side-effect outcomes (composite, reactive antiemetic,
#   antipruritic, naloxone reversal) use the same
#   anaesthesia-stop-anchored definitions as global (S4) reporting in
#   4_flags.R (Flag 9), computed once there and read directly from
#   pca_episodes_flagged rather than recomputed independently in this
#   script:
#     - Naloxone reversal: minimal/no cutoff — true reversal events
#       sit well outside any plausible washout window
#     - Antipruritic: no time filter (chlorphenamine + naloxone at
#       the 40mcg pruritus-consistent dose) — sensitivity-curve
#       testing found no natural break anywhere in the 0-48h range,
#       consistent with itch tracking cumulative opioid exposure
#       rather than a clearable perioperative confound
#     - Reactive antiemetic: >12h post-anaesthesia-stop (not episode
#       start, to separate the general anaesthetic/volatile agent
#       effect from the postoperative opioid-PCA effect), excluding
#       nonopioid_nausea_confound_final (a patient-level flag
#       combining specialty/OPCS-category risk with confirmed
#       pre-op symptomatic status — see 4_flags.R Flag 8/9). No
#       burden-percentage threshold is applied: a discretisation
#       cliff at exactly 50% burden affects around 8% of episodes in
#       both directions, and high burden is more consistent with a
#       probable non-opioid cause than a persistent opioid-
#       attributable pattern once checked against known confounds.
#
#   Since S3 is 100% surgical, every episode here follows the
#   "surgical" branch of each rule; the medical-episode logic that
#   applies in S4/global reporting is not relevant in this script.
#   any_laxative uses no time filter, consistent with the rest of the
#   project.
#
# Age-stratified / peak-pain diagnostics
#
#   The exploratory diagnostic sections further down in this script
#   read peak_pain_score and time_weighted_mean_pain directly from
#   pca_episodes_flagged (computed once in 4_flags.R) rather than
#   recomputing a separate per-episode median pain score, so there is
#   one pain-score methodology in use across the project rather than
#   two independently-derived ones.
#
# Population:
#   Pre-match: morphine n=1,718 / oxycodone n=433
#   (First line, first episode, no concurrent epidural,
#    surgical only, the predefined exclusion set applied)
#
# Method:
#   Nearest neighbour 1:1 matching
#   Caliper: 0.2 SD of logit propensity score
#   Package: MatchIt
#
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
  exists("s3_morphine_keys"),
  exists("s3_oxycodone_keys"),
  exists("intraop_mar_clean"),
  exists("anaesthesia_ax_clean")   # needed for asa_grade
)

cat("============================================================\n")
cat("ANALYSIS 6 — PSM COMPARATIVE ANALYSIS (SECTION 3)\n")
cat("============================================================\n\n")

cat("Not sourcing 9_comsumption.R here — it requires\n")
cat("matched_outcomes, which doesn't exist until much later in this\n")
cat("script (that would be a circular dependency). Not needed anyway:\n")
cat("the Dose Outcomes section below reads episode_true_total_dose.rds\n")
cat("/ episode_dose.rds directly.\n\n")

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

cat("Pre-match pool (before neonate exclusion):\n")
cat(sprintf("  Morphine:   n = %d\n", nrow(morphine_pool)))
cat(sprintf("  Oxycodone:  n = %d\n", nrow(oxycodone_pool)))
cat("\n")

# ============================================================
# EXCLUDE NEONATES
# ============================================================

if ("age_group" %in% names(morphine_pool)) {
  n_neo_mor <- sum(morphine_pool$age_group == "Neonate (<28 days)", na.rm = TRUE)
  n_neo_oxy <- sum(oxycodone_pool$age_group == "Neonate (<28 days)", na.rm = TRUE)
  cat("--- Excluding neonates ---\n")
  cat("Neonates found: morphine =", n_neo_mor, "| oxycodone =", n_neo_oxy, "\n")
  
  morphine_pool  <- morphine_pool  |> filter(age_group != "Neonate (<28 days)")
  oxycodone_pool <- oxycodone_pool |> filter(age_group != "Neonate (<28 days)")
  
  cat("Post-exclusion: morphine n =", nrow(morphine_pool),
      "| oxycodone n =", nrow(oxycodone_pool), "\n\n")
} else {
  cat("age_group column not found — cannot confirm neonate exclusion.\n")
  cat("   Check manually before treating this population as final.\n\n")
}

# ============================================================
# PATIENT-LEVEL DEDUPLICATION (July 2026) — BUILT INTO PSM ITSELF
#
# S3 already restricts to is_first_episode==TRUE (first episode WITHIN
# an admission), but this doesn't prevent the SAME real patient
# appearing more than once if they have multiple admissions (e.g.
# readmitted for a separate procedure, or rotated back to the same
# drug in a later admission). Without this fix, such a patient could
# contribute two rows to the same drug's matching pool, effectively
# double-counting their outcome and inflating apparent significance
# for continuous outcomes especially (their consistent dose/pain
# tendency counted twice) — the same mechanism that made peak pain
# score a clustering artefact in the old, undeduplicated version of
# this analysis (p=0.046 -> p=0.475 once checked).
#
# Fix: deduplicate each drug's pool to ONE episode per true patient
# (PAT_ID), keeping the earliest episode — BEFORE matching runs, so
# the matching pool itself, not just the final reported outcome, is
# clean. This is now the DEFAULT behaviour of this script, not an
# optional sensitivity check.
# ============================================================

stopifnot(exists("admissions_clean"))

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

# ============================================================
# [SEPARATION FIX 1] Drop GENDER == "Not specified"
# ============================================================

cat("--- GENDER == 'Not specified' check ---\n\n")

n_not_specified <- sum(psm_data$GENDER == "Not specified", na.rm = TRUE)
cat("Episodes with GENDER == 'Not specified':", n_not_specified,
    "out of", nrow(psm_data), "\n")

if (n_not_specified > 0) {
  cat("Dropping these from the matching population (data-entry gap,\n")
  cat("not a genuine third category — retaining would reintroduce\n")
  cat("separation in the propensity model).\n\n")
  psm_data <- psm_data |> filter(GENDER != "Not specified")
} else {
  cat("None found in this pull — no action needed.\n\n")
}

cat("Pre-match pool after GENDER fix:", nrow(psm_data), "\n\n")

# ============================================================
# CLEAN ADMITTING_SPECIALTY -> specialty_collapsed
# ============================================================

cat("--- Cleaning ADMITTING_SPECIALTY -> specialty_collapsed ---\n\n")

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
# BUILD EXPANDED INTRAOP COVARIATE FLAGS
# ============================================================

cat("--- Building expanded intraop covariate flags ---\n\n")

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

cat("New intraop flags added: had_intraop_paracetamol, had_intraop_nsaid,\n")
cat("  had_intraop_magnesium, had_intraop_fentanyl, had_intraop_morphine,\n")
cat("  had_intraop_oxycodone2, had_intraop_diamorphine\n")
cat("(had_intraop_alfentanil excluded — no imbalance signal, p=0.41)\n\n")

# ============================================================
# ADD anaesthesia_duration_mins
# ============================================================

cat("--- Adding anaesthesia_duration_mins covariate ---\n\n")

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
  cat("Imputed", n_missing_duration, "missing values with median (",
      round(median_duration, 1), "mins ) — small enough proportion",
      "not to warrant a more complex imputation strategy.\n\n")
}

# ============================================================
# ADD asa_grade COVARIATE
# ============================================================
# Clean addition — S3 is 100% surgical by construction, so unlike S4
# there is no "not applicable for medical patients" coding problem to
# solve. Sourced from anaesthesia_ax_clean (2_cleaningscript.R
# Section 7B), joined by PAT_ENC_CSN_ID. If an admission has more than
# one anaesthesia assessment row, the first non-missing grade is used
# — worth re-checking this is a reasonable choice if any admission
# genuinely has two different ASA grades recorded (unlikely within a
# single admission, but not yet explicitly verified).
#
# PROCEDURE URGENCY: deliberately not added — see header note. Too
# sparse (~17.6% coverage, declining) to support inclusion in a
# matching formula.
# ============================================================

cat("--- Adding asa_grade covariate ---\n\n")

asa_data <- anaesthesia_ax_clean |>
  filter(!is.na(asa_grade)) |>
  distinct(PAT_ENC_CSN_ID, .keep_all = TRUE) |>
  select(PAT_ENC_CSN_ID, asa_grade)

psm_data <- psm_data |> left_join(asa_data, by = "PAT_ENC_CSN_ID")

n_missing_asa <- sum(is.na(psm_data$asa_grade))
cat("asa_grade missing:", n_missing_asa, "/", nrow(psm_data),
    sprintf("(%.1f%%)\n", 100 * n_missing_asa / nrow(psm_data)))

if (n_missing_asa > 0) {
  cat("Non-zero missingness in asa_grade for a 100%% surgical\n")
  cat("population is unexpected — worth checking whether this reflects\n")
  cat("genuine documentation gaps or a linkage problem before imputing.\n")
  cat("Using mode imputation (most common grade) as a placeholder —\n")
  cat("revisit if missingness is substantial (>5%%).\n\n")
  asa_mode <- names(sort(table(psm_data$asa_grade), decreasing = TRUE))[1]
  psm_data <- psm_data |>
    mutate(asa_grade = if_else(is.na(asa_grade), asa_mode, asa_grade))
} else {
  cat("No missingness — no imputation needed.\n\n")
}

psm_data <- psm_data |> mutate(asa_grade = factor(asa_grade, ordered = TRUE))

cat("ASA grade distribution (pre-match):\n")
print(table(psm_data$drug, psm_data$asa_grade))
cat("\n")

# ============================================================
# CHECK MATCHING VARIABLES
# ============================================================

cat("--- Checking matching variable availability ---\n")
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
# PRE-MATCH BALANCE TABLE
# ============================================================

cat("============================================================\n")
cat("TABLE 1 — PRE-MATCH BALANCE\n")
cat("============================================================\n\n")

pre_match_balance <- psm_data |>
  group_by(drug) |>
  summarise(
    n                   = n(),
    age_median          = round(median(AGE, na.rm = TRUE), 1),
    age_q25             = round(quantile(AGE, 0.25, na.rm = TRUE), 1),
    age_q75             = round(quantile(AGE, 0.75, na.rm = TRUE), 1),
    weight_median       = round(median(WEIGHT_KG, na.rm = TRUE), 1),
    pct_male            = round(100 * mean(GENDER == "Male", na.rm = TRUE), 1),
    pct_pca             = round(100 * mean(mode == "PCA", na.rm = TRUE), 1),
    pct_remifentanil    = round(100 * mean(had_remifentanil, na.rm = TRUE), 1),
    pct_ketamine        = round(100 * mean(had_intraop_ketamine, na.rm = TRUE), 1),
    pct_clonidine       = round(100 * mean(had_intraop_clonidine, na.rm = TRUE), 1),
    pct_magnesium       = round(100 * mean(had_intraop_magnesium, na.rm = TRUE), 1),
    pct_nsaid           = round(100 * mean(had_intraop_nsaid, na.rm = TRUE), 1),
    pct_paracetamol     = round(100 * mean(had_intraop_paracetamol, na.rm = TRUE), 1),
    pct_intraop_oxy     = round(100 * mean(had_intraop_oxycodone2, na.rm = TRUE), 1),
    pct_intraop_mor     = round(100 * mean(had_intraop_morphine, na.rm = TRUE), 1),
    pct_spinal          = round(100 * mean(neuraxial_intrathecal_mar, na.rm = TRUE), 1),
    pct_peripheral      = round(100 * mean(had_peripheral, na.rm = TRUE), 1),
    asa_median          = round(median(as.numeric(as.character(asa_grade)), na.rm = TRUE), 1),
    .groups = "drop"
  )
print(pre_match_balance)
cat("\n")

# ============================================================
# [SEPARATION FIX 2] Collapse proc_category "Other" into "Unclassified"
# ============================================================

cat("--- proc_category 'Other' -> 'Unclassified' ---\n\n")

n_other <- sum(psm_data$proc_category == "Other", na.rm = TRUE)
cat("Episodes with proc_category == 'Other':", n_other, "\n")

psm_data <- psm_data |>
  mutate(
    proc_category = if_else(proc_category == "Other", "Unclassified", proc_category)
  )

cat("Collapsed into 'Unclassified'. New proc_category level count:",
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

psm_formula_previous <- treated ~ AGE + WEIGHT_KG + GENDER + proc_category +
  mode + had_peripheral + year + specialty_collapsed +
  had_remifentanil + had_intraop_ketamine + had_intraop_clonidine +
  neuraxial_intrathecal_mar

psm_formula_expanded <- treated ~ AGE + WEIGHT_KG + GENDER + proc_category +
  mode + had_peripheral + year + specialty_collapsed + asa_grade +
  had_remifentanil + had_intraop_ketamine + had_intraop_clonidine +
  neuraxial_intrathecal_mar +
  had_intraop_paracetamol + had_intraop_nsaid + had_intraop_magnesium +
  had_intraop_fentanyl + had_intraop_morphine + had_intraop_oxycodone2 +
  had_intraop_diamorphine + anaesthesia_duration_mins

# ============================================================
# SEPARATION / CONVERGENCE CHECK
# (Re-run after adding asa_grade — a sparse grade, e.g. very few
# grade 4/5 concentrated in one arm, could reintroduce the same
# instability GENDER/proc_category caused before.)
# ============================================================

cat("============================================================\n")
cat("[CHECK] Separation/convergence check — expanded formula (incl. asa_grade)\n")
cat("============================================================\n\n")

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
  cat("\nCheck whether asa_grade is the new source of instability —\n")
  cat("e.g. a sparse ASA grade x proc_category combination — before\n")
  cat("assuming it's the same GENDER/proc_category issue as before.\n\n")
} else {
  cat("No GLM warnings — model converged without apparent separation.\n\n")
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

cat("============================================================\n\n")

# ============================================================
# DIAGNOSTIC COMPARISON — PREVIOUS vs EXPANDED
# ============================================================

cat("============================================================\n")
cat("[DIAGNOSTIC] PREVIOUS (remi/ketamine/clonidine only) vs\n")
cat("             EXPANDED (+ magnesium/nsaid/paracetamol/opioids/asa_grade)\n")
cat("============================================================\n\n")

set.seed(42)
match_diag_previous <- tryCatch(
  matchit(psm_formula_previous, data = psm_data, method = "nearest",
          distance = "glm", link = "logit", ratio = 1,
          caliper = 0.2, std.caliper = TRUE),
  error = function(e) { cat("PREVIOUS match error:", conditionMessage(e), "\n"); NULL }
)

set.seed(42)
match_diag_expanded <- tryCatch(
  matchit(psm_formula_expanded, data = psm_data, method = "nearest",
          distance = "glm", link = "logit", ratio = 1,
          caliper = 0.2, std.caliper = TRUE),
  error = function(e) { cat("EXPANDED match error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(match_diag_previous) && !is.null(match_diag_expanded)) {
  n_prev <- sum(match_diag_previous$weights == 1) / 2
  n_exp  <- sum(match_diag_expanded$weights == 1) / 2
  
  smd_prev <- summary(match_diag_previous, standardize = TRUE)$sum.matched
  smd_exp  <- summary(match_diag_expanded, standardize = TRUE)$sum.matched
  max_smd_prev <- max(abs(smd_prev[, "Std. Mean Diff."]), na.rm = TRUE)
  max_smd_exp  <- max(abs(smd_exp[, "Std. Mean Diff."]), na.rm = TRUE)
  
  cat(sprintf("Matched pairs — PREVIOUS: %d | EXPANDED: %d\n", n_prev, n_exp))
  cat(sprintf("Max |SMD| post-match — PREVIOUS: %.3f | EXPANDED: %.3f\n", max_smd_prev, max_smd_exp))
  cat("(Target <0.1 for every covariate — check full tables if either exceeds this)\n\n")
  cat("--- Full balance table, PREVIOUS ---\n"); print(smd_prev)
  cat("\n--- Full balance table, EXPANDED ---\n"); print(smd_exp)
  
  cat("\n--- Balance on asa_grade specifically (EXPANDED match only) ---\n")
  asa_rows <- grep("asa_grade", rownames(smd_exp), value = TRUE)
  if (length(asa_rows) > 0) {
    print(smd_exp[asa_rows, , drop = FALSE])
  }
} else {
  cat("Comparison skipped — one or both diagnostic matches failed.\n")
}
cat("\n============================================================\n\n")

# ------------------------------------------------------------
# PRODUCTION MATCH — uses EXPANDED formula (now incl. asa_grade)
# ------------------------------------------------------------

psm_formula <- psm_formula_expanded

cat("Running MatchIt (production match, expanded formula incl. asa_grade)...\n")
set.seed(42)

match_out <- tryCatch({
  matchit(
    formula  = psm_formula,
    data     = psm_data,
    method   = "nearest",
    distance = "glm",
    link     = "logit",
    ratio    = 1,
    caliper  = 0.2,
    std.caliper = TRUE
  )
}, error = function(e) {
  cat("MatchIt error:", conditionMessage(e), "\n")
  NULL
})

if (is.null(match_out)) {
  cat("Matching failed — check variable missingness and factor levels above\n")
  stop("PSM failed")
}

cat("Matching complete\n\n")
print(summary(match_out))

# ============================================================
# EXTRACT MATCHED DATASET
# ============================================================

matched_data <- match.data(match_out)

n_matched_oxy  <- sum(matched_data$treated == 1)
n_matched_mor  <- sum(matched_data$treated == 0)
n_unmatched    <- nrow(oxycodone_pool) - n_matched_oxy

cat(sprintf("\nMatched pairs:          %d\n", n_matched_oxy))
cat(sprintf("Oxycodone matched:      %d / %d\n", n_matched_oxy, nrow(oxycodone_pool)))
cat(sprintf("Morphine matched:       %d / %d\n", n_matched_mor, nrow(morphine_pool)))
cat(sprintf("Oxycodone unmatched:    %d\n\n", n_unmatched))

matched_keys <- matched_data |>
  select(PAT_ENC_CSN_ID, episode_id, drug, treated, subclass, distance, weights)

saveRDS(matched_keys, "1 - data/3 - processed_data/psm_matched_keys.rds")
cat("Matched keys saved to psm_matched_keys.rds\n\n")

# ============================================================
# POST-MATCH BALANCE TABLE
# ============================================================

cat("============================================================\n")
cat("TABLE 2 — POST-MATCH BALANCE\n")
cat("============================================================\n\n")

post_match_balance <- matched_data |>
  group_by(drug) |>
  summarise(
    n                = n(),
    age_median       = round(median(AGE, na.rm = TRUE), 1),
    age_q25          = round(quantile(AGE, 0.25, na.rm = TRUE), 1),
    age_q75          = round(quantile(AGE, 0.75, na.rm = TRUE), 1),
    weight_median    = round(median(WEIGHT_KG, na.rm = TRUE), 1),
    pct_male         = round(100 * mean(GENDER == "Male", na.rm = TRUE), 1),
    pct_pca          = round(100 * mean(mode == "PCA", na.rm = TRUE), 1),
    pct_remifentanil = round(100 * mean(had_remifentanil, na.rm = TRUE), 1),
    pct_ketamine     = round(100 * mean(had_intraop_ketamine, na.rm = TRUE), 1),
    pct_clonidine    = round(100 * mean(had_intraop_clonidine, na.rm = TRUE), 1),
    pct_magnesium    = round(100 * mean(had_intraop_magnesium, na.rm = TRUE), 1),
    pct_nsaid        = round(100 * mean(had_intraop_nsaid, na.rm = TRUE), 1),
    pct_paracetamol  = round(100 * mean(had_intraop_paracetamol, na.rm = TRUE), 1),
    pct_intraop_oxy  = round(100 * mean(had_intraop_oxycodone2, na.rm = TRUE), 1),
    pct_intraop_mor  = round(100 * mean(had_intraop_morphine, na.rm = TRUE), 1),
    pct_spinal       = round(100 * mean(neuraxial_intrathecal_mar, na.rm = TRUE), 1),
    pct_peripheral   = round(100 * mean(had_peripheral, na.rm = TRUE), 1),
    asa_median       = round(median(as.numeric(as.character(asa_grade)), na.rm = TRUE), 1),
    .groups = "drop"
  )
print(post_match_balance)
cat("\n")

cat("--- Standardised Mean Differences (post-match) ---\n")
cat("(Values <0.1 indicate good balance)\n\n")
smd_summary <- summary(match_out, standardize = TRUE)
print(smd_summary$sum.matched)
cat("\n")

# ============================================================
# JOIN OUTCOMES — NOW USING FINALISED SIDE-EFFECT FLAGS
# ============================================================
# Reactive antiemetic, antipruritic, naloxone reversal, and composite
# side effect are read DIRECTLY from pca_episodes_flagged
# (reactive_antiemetic_final, any_antipruritic_final,
# any_naloxone_reversal_final, any_side_effect_final) — built ONCE in
# 4_flags.R Flag 9, shared with global/S4 reporting. This replaces the
# old inline shift-count-based recomputation entirely (previously ~50
# lines here duplicating logic that now lives in one place).
# any_laxative is unchanged, still sourced from episode_metrics.
# ============================================================

cat("--- Joining episode metrics + finalised side-effect flags to matched dataset ---\n\n")

matched_outcomes <- matched_data |>
  select(PAT_ENC_CSN_ID, episode_id, drug, treated, subclass) |>
  left_join(
    episode_metrics |> select(-drug, -mode),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  ) |>
  left_join(
    pca_episodes_flagged |>
      select(PAT_ENC_CSN_ID, episode_id,
             reactive_antiemetic_final, any_antipruritic_final,
             any_naloxone_reversal_final, any_side_effect_final),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  )

n_missing_metrics <- sum(is.na(matched_outcomes$any_side_effect_final))
cat("Episodes missing finalised outcome data:", n_missing_metrics, "\n\n")

cat(sprintf("Reactive antiemetic (final): n=%d\n",
            sum(matched_outcomes$reactive_antiemetic_final, na.rm = TRUE)))
cat(sprintf("Antipruritic (final): n=%d\n",
            sum(matched_outcomes$any_antipruritic_final, na.rm = TRUE)))
cat(sprintf("Naloxone reversal (final): n=%d\n",
            sum(matched_outcomes$any_naloxone_reversal_final, na.rm = TRUE)))
cat(sprintf("Composite side effect (final): n=%d\n\n",
            sum(matched_outcomes$any_side_effect_final, na.rm = TRUE)))

# ============================================================
# TABLE 3 — PRIMARY OUTCOME
# ============================================================

cat("============================================================\n")
cat("TABLE 3 — PRIMARY OUTCOME: COMPOSITE SIDE EFFECT\n")
cat("(Anaesthesia-stop-anchored definitions — see header note)\n")
cat("============================================================\n\n")

primary_outcome <- matched_outcomes |>
  group_by(drug) |>
  summarise(
    n               = n(),
    n_side_effect   = sum(any_side_effect_final, na.rm = TRUE),
    pct_side_effect = round(100 * mean(any_side_effect_final, na.rm = TRUE), 1),
    n_antiemetic    = sum(reactive_antiemetic_final, na.rm = TRUE),
    pct_antiemetic  = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    n_antipruritic  = sum(any_antipruritic_final, na.rm = TRUE),
    pct_antiprurit  = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
    n_nalox_rev     = sum(any_naloxone_reversal_final, na.rm = TRUE),
    .groups = "drop"
  )
print(primary_outcome)
cat("\n")

cat("--- Binary outcomes (chi-squared) ---\n\n")
for (var in c("any_side_effect_final", "reactive_antiemetic_final",
              "any_antipruritic_final")) {
  tab <- table(matched_outcomes$drug, matched_outcomes[[var]])
  if (ncol(tab) < 2) {
    cat(sprintf("%-30s — insufficient variation\n", var))
    next
  }
  ct  <- suppressWarnings(chisq.test(tab))
  or  <- (tab["oxycodone", 2] / tab["oxycodone", 1]) /
    (tab["morphine",  2] / tab["morphine",  1])
  cat(sprintf("%-30s OR=%.2f  p=%.4f\n", var, or, ct$p.value))
}
cat("\n")

cat("============================================================\n")
cat("TABLE 4 — SECONDARY OUTCOMES (matched pairs)\n")
cat("============================================================\n\n")

secondary_outcomes <- matched_outcomes |>
  group_by(drug) |>
  summarise(
    n                         = n(),
    pct_antiemetic            = round(100 * mean(reactive_antiemetic_final,   na.rm = TRUE), 1),
    pct_antipruritic          = round(100 * mean(any_antipruritic_final,      na.rm = TRUE), 1),
    pct_laxative              = round(100 * mean(any_laxative,                na.rm = TRUE), 1),
    n_nalox_reversal          = sum(any_naloxone_reversal_final,              na.rm = TRUE),
    pct_severe_median         = round(median(pct_shifts_severe,              na.rm = TRUE), 1),
    pct_severe_q25            = round(quantile(pct_shifts_severe, 0.25,      na.rm = TRUE), 1),
    pct_severe_q75            = round(quantile(pct_shifts_severe, 0.75,      na.rm = TRUE), 1),
    first_controlled_median   = round(median(first_controlled_shift,         na.rm = TRUE), 1),
    duration_median           = round(median(episode_duration_hrs,           na.rm = TRUE), 1),
    duration_q25              = round(quantile(episode_duration_hrs, 0.25,   na.rm = TRUE), 1),
    duration_q75              = round(quantile(episode_duration_hrs, 0.75,   na.rm = TRUE), 1),
    paracetamol_pct_median    = round(median(pct_shifts_paracetamol,         na.rm = TRUE), 1),
    nsaid_pct_median          = round(median(pct_shifts_nsaid,               na.rm = TRUE), 1),
    adjuvant_pct_median       = round(median(pct_shifts_adjuvant,            na.rm = TRUE), 1),
    .groups = "drop"
  )
print(t(secondary_outcomes))
cat("\n")

# ============================================================
# PAIN OUTCOMES BY AGE GROUP (matched pairs)
#
# This is an internal, exploratory diagnostic section, not one of
# the manuscript's numbered tables.
# ============================================================

cat("============================================================\n")
cat("PAIN OUTCOMES BY AGE GROUP (matched pairs)\n")
cat("============================================================\n\n")

pain_matched <- matched_outcomes |>
  left_join(
    pca_episodes_flagged |>
      select(PAT_ENC_CSN_ID, episode_id, age_group),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  )

pain_matched |>
  group_by(drug, age_group) |>
  summarise(
    n                    = n(),
    pct_severe_median    = round(median(pct_shifts_severe, na.rm = TRUE), 1),
    pct_severe_q25       = round(quantile(pct_shifts_severe, 0.25, na.rm = TRUE), 1),
    pct_severe_q75       = round(quantile(pct_shifts_severe, 0.75, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  arrange(age_group, drug) |>
  print()
cat("\n")

cat("NOTE ON THE ABOVE (% shifts severe, >=7 Serlin cutpoint):\n")
cat("This threshold is the SI (osteotomy) project's core\n")
cat("stratification variable, not independently validated as a drug-\n")
cat("comparison efficacy endpoint for this analysis. Retained above for\n")
cat("continuity, but the peak/time-weighted pain score section below\n")
cat("is the recommended primary pain comparison.\n\n")

# ============================================================
# PEAK & TIME-WEIGHTED MEAN PAIN SCORE
#
# Reads peak_pain_score and time_weighted_mean_pain directly from
# pca_episodes_flagged (computed once in 4_flags.R, specifically to
# avoid the floor-effect problem a plain median suffers from — see
# that script's header) rather than recomputing an independent
# pain-score comparison from pain_scores_clean, so there is one
# pain-score methodology in use across the project rather than two.
# ============================================================

cat("============================================================\n")
cat("PEAK & TIME-WEIGHTED MEAN PAIN SCORE (matched pairs)\n")
cat("============================================================\n\n")

stopifnot(
  "peak_pain_score" %in% names(pca_episodes_flagged),
  "time_weighted_mean_pain" %in% names(pca_episodes_flagged)
)

episode_pain_summary <- pca_episodes_flagged |>
  filter(PAT_ENC_CSN_ID %in% matched_data$PAT_ENC_CSN_ID,
         episode_id %in% matched_data$episode_id) |>
  select(PAT_ENC_CSN_ID, episode_id, age_group,
         peak_pain_score, time_weighted_mean_pain) |>
  left_join(matched_data |> select(PAT_ENC_CSN_ID, episode_id, drug),
            by = c("PAT_ENC_CSN_ID", "episode_id"))

cat("Episodes with a peak_pain_score value:",
    sum(!is.na(episode_pain_summary$peak_pain_score)),
    "/", nrow(episode_pain_summary), "\n")
cat("Episodes with a time_weighted_mean_pain value:",
    sum(!is.na(episode_pain_summary$time_weighted_mean_pain)),
    "/", nrow(episode_pain_summary), "\n\n")

cat("--- OVERALL: peak vs time-weighted mean pain score, by drug ---\n\n")
print(episode_pain_summary |> group_by(drug) |>
        summarise(n = n(),
                  peak_median          = median(peak_pain_score, na.rm=TRUE),
                  peak_q25             = quantile(peak_pain_score, 0.25, na.rm=TRUE),
                  peak_q75             = quantile(peak_pain_score, 0.75, na.rm=TRUE),
                  time_weighted_median = median(time_weighted_mean_pain, na.rm=TRUE),
                  time_weighted_q75    = quantile(time_weighted_mean_pain, 0.75, na.rm=TRUE),
                  .groups = "drop"))

wt_peak          <- wilcox.test(peak_pain_score ~ drug, data = episode_pain_summary)
wt_time_weighted <- suppressWarnings(wilcox.test(time_weighted_mean_pain ~ drug, data = episode_pain_summary))
cat("\nWilcoxon peak score p =          ", format.pval(wt_peak$p.value, digits = 3), "\n")
cat("Wilcoxon time-weighted mean p =  ", format.pval(wt_time_weighted$p.value, digits = 3), "\n\n")

cat("--- BY AGE GROUP: peak pain score ---\n\n")
for (ag in unique(episode_pain_summary$age_group)) {
  sub <- episode_pain_summary |> filter(age_group == ag)
  if (n_distinct(sub$drug) < 2 || nrow(sub) < 10) next
  cat("---", ag, "(n =", nrow(sub), ") ---\n")
  print(sub |> group_by(drug) |>
          summarise(n = n(),
                    peak_median = median(peak_pain_score, na.rm=TRUE),
                    peak_q25    = quantile(peak_pain_score, 0.25, na.rm=TRUE),
                    peak_q75    = quantile(peak_pain_score, 0.75, na.rm=TRUE),
                    .groups = "drop"))
  wt_sub <- tryCatch(wilcox.test(peak_pain_score ~ drug, data = sub), error = function(e) NULL)
  if (!is.null(wt_sub)) cat("Wilcoxon (peak) p =", format.pval(wt_sub$p.value, digits = 3), "\n")
  cat("\n")
}

cat("6 comparisons run above (1 overall + up to 5 strata) — NOT\n")
cat("multiple-comparison corrected. Treat as exploratory/hypothesis-\n")
cat("generating, not confirmatory.\n\n")

cat("============================================================\n")
cat("TABLE 6 — STATISTICAL TESTS (matched pairs)\n")
cat("Note: chi-squared used here; use conditional logistic for paper\n")
cat("============================================================\n\n")

cat("--- Binary outcomes (chi-squared) ---\n\n")
for (var in c("any_side_effect_final", "reactive_antiemetic_final",
              "any_antipruritic_final",
              "any_naloxone_reversal_final")) {
  tab <- table(matched_outcomes$drug, matched_outcomes[[var]])
  if (ncol(tab) < 2) {
    cat(sprintf("%-35s — insufficient variation\n", var))
    next
  }
  ct <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
  if (is.null(ct)) next
  or <- (tab["oxycodone", 2] / tab["oxycodone", 1]) /
    (tab["morphine",  2] / tab["morphine",  1])
  cat(sprintf("%-35s OR=%.2f  p=%.4f\n", var, or, ct$p.value))
}
cat("\n")

cat("--- Continuous outcomes (Wilcoxon rank-sum) ---\n\n")
mor <- matched_outcomes |> filter(drug == "morphine")
oxy <- matched_outcomes |> filter(drug == "oxycodone")
for (var in c("pct_shifts_severe", "first_controlled_shift",
              "episode_duration_hrs", "pct_shifts_nsaid", "pct_shifts_paracetamol")) {
  x  <- oxy[[var]][!is.na(oxy[[var]])]
  y  <- mor[[var]][!is.na(mor[[var]])]
  wt <- wilcox.test(x, y)
  cat(sprintf("%-35s p=%.4f\n", var, wt$p.value))
}
cat("\n")

cat("\n--- Laxative (any_laxative_full) — DESCRIPTIVE ONLY, not a formal outcome ---\n")
cat("(Cannot distinguish prophylactic bowel-regimen dosing from a\n")
cat(" genuine reactive response to opioid-induced constipation; no\n")
cat(" order-type field exists to separate these. Reported as a\n")
cat(" descriptive rate only, not tested statistically.)\n\n")

matched_outcomes |>
  group_by(drug) |>
  summarise(n = n(), pct_laxative = round(100 * mean(any_laxative_full, na.rm = TRUE), 1),
            .groups = "drop") |>
  print()
cat("\n")

# ============================================================
# DOSE OUTCOMES — EFFICACY SECTION (unchanged from before)
# ============================================================

cat("============================================================\n")
cat("DOSE OUTCOMES (EFFICACY) — MATCHED PAIRS [bolus-inclusive]\n")
cat("Requires episode_true_total_dose.rds from 8_true_total_dose.R\n")
cat("============================================================\n\n")

episode_true_total_dose <- tryCatch(
  readRDS("1 - data/3 - processed_data/episode_true_total_dose.rds"),
  error = function(e) {
    cat("episode_true_total_dose.rds not found — run 8_true_total_dose.R first\n")
    return(NULL)
  }
)

episode_dose <- tryCatch(
  readRDS("1 - data/3 - processed_data/episode_dose.rds"),
  error = function(e) {
    cat("episode_dose.rds not found — run 7_dosing.R first\n")
    return(NULL)
  }
)

if (!is.null(episode_true_total_dose)) {
  
  matched_total_dose <- episode_true_total_dose |>
    semi_join(matched_keys, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    filter(drug %in% c("morphine", "oxycodone"))
  
  linkage_by_drug <- matched_keys |>
    select(PAT_ENC_CSN_ID, episode_id, drug) |>
    left_join(
      episode_true_total_dose |> select(PAT_ENC_CSN_ID, episode_id) |> mutate(is_linked = TRUE),
      by = c("PAT_ENC_CSN_ID", "episode_id")
    ) |>
    mutate(is_linked = replace_na(is_linked, FALSE))
  
  cat("--- Flowsheet linkage rate by drug (this matched population) ---\n")
  print(linkage_by_drug |> group_by(drug) |>
          summarise(n = n(), n_linked = sum(is_linked),
                    pct_linked = round(100*mean(is_linked), 1), .groups = "drop"))
  link_test <- suppressWarnings(chisq.test(table(linkage_by_drug$drug, linkage_by_drug$is_linked)))
  cat("Chi-squared (linkage x drug) p =", format.pval(link_test$p.value, digits = 3), "\n")
  cat("(Non-random linkage confirmed separately not to change the conclusion —\n")
  cat(" see ipw_correction_total_dose.R for the full sensitivity check.)\n\n")
  
  cat(sprintf("Matched episodes with bolus-inclusive dose data: %d / %d (%.1f%%)\n\n",
              nrow(matched_total_dose), nrow(matched_keys),
              100*nrow(matched_total_dose)/nrow(matched_keys)))
  
  cat("--- PRIMARY: Total dose (background + bolus), mcg/kg/hr ---\n")
  total_dose_summary <- matched_total_dose |>
    group_by(drug) |>
    summarise(n = n(), total_dose = fmt_med_iqr(total_dose_mcg_kg_hr, 2), .groups = "drop")
  print(total_dose_summary)
  
  wt_total <- wilcox.test(total_dose_mcg_kg_hr ~ drug, data = matched_total_dose)
  cat("Wilcoxon (TOTAL dose) p =", format.pval(wt_total$p.value, digits = 3), "\n\n")
  
  cat("--- SECONDARY/REFERENCE: Background only, SAME linked subset ---\n")
  bg_summary <- matched_total_dose |>
    group_by(drug) |>
    summarise(n = n(), background = fmt_med_iqr(background_mcg_kg_hr, 2), .groups = "drop")
  print(bg_summary)
  
  wt_bg <- wilcox.test(background_mcg_kg_hr ~ drug, data = matched_total_dose)
  cat("Wilcoxon (background only, same subset) p =", format.pval(wt_bg$p.value, digits = 3), "\n\n")
  
  cat("--- Bolus as % of total dose, by drug ---\n")
  bolus_pct_summary <- matched_total_dose |>
    group_by(drug) |>
    summarise(n = n(), bolus_pct = fmt_med_iqr(bolus_pct_of_total, 1), .groups = "drop")
  print(bolus_pct_summary)
  
  wt_bolus <- wilcox.test(bolus_pct_of_total ~ drug, data = matched_total_dose)
  cat("Wilcoxon (bolus % of total) p =", format.pval(wt_bolus$p.value, digits = 3), "\n\n")
  
} else {
  cat("Bolus-inclusive dose outcomes skipped — run 8_true_total_dose.R first\n\n")
  matched_total_dose <- NULL
}

if (!is.null(episode_dose)) {
  matched_loading <- episode_dose |>
    semi_join(matched_keys, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    filter(drug %in% c("morphine", "oxycodone"))
  
  cat("--- Loading dose prescribed, n (%) [from episode_dose.rds] ---\n")
  load_summary <- matched_loading |>
    group_by(drug) |>
    summarise(n = n(), pct_loading = round(100 * mean(had_loading_dose, na.rm = TRUE), 1), .groups = "drop")
  print(load_summary)
  
  load_tab <- table(matched_loading$drug, matched_loading$had_loading_dose)
  if (ncol(load_tab) == 2) {
    ct <- suppressWarnings(chisq.test(load_tab))
    cat("had_loading_dose p =", format.pval(ct$p.value, digits = 3), "\n")
  }
  cat("\n")
} else {
  matched_loading <- NULL
}

# ============================================================
# SUMMARY TABLE FOR MANUSCRIPT
# ============================================================

cat("============================================================\n")
cat("MANUSCRIPT SUMMARY TABLE\n")
cat("(Morphine vs Oxycodone — matched pairs)\n")
cat("Population: Surgical cases (major OPCS only)\n")
cat("Side effects: anaesthesia-stop-anchored definitions (see header)\n")
cat("============================================================\n\n")

mor <- matched_outcomes |> filter(drug == "morphine")
oxy <- matched_outcomes |> filter(drug == "oxycodone")

n_mor <- nrow(mor)
n_oxy <- nrow(oxy)

cat(sprintf("%-40s  Morphine (n=%d)   Oxycodone (n=%d)\n", "", n_mor, n_oxy))
cat(strrep("-", 75), "\n")

print_row <- function(label, mor_val, oxy_val) {
  cat(sprintf("%-40s  %-20s %s\n", label, mor_val, oxy_val))
}

print_row("Composite side effect — n (%)",
          fmt_n_pct(sum(mor$any_side_effect_final, na.rm=TRUE), n_mor),
          fmt_n_pct(sum(oxy$any_side_effect_final, na.rm=TRUE), n_oxy))

print_row("Antiemetic (reactive) — n (%)",
          fmt_n_pct(sum(mor$reactive_antiemetic_final, na.rm=TRUE), n_mor),
          fmt_n_pct(sum(oxy$reactive_antiemetic_final, na.rm=TRUE), n_oxy))

print_row("Antipruritic — n (%)",
          fmt_n_pct(sum(mor$any_antipruritic_final, na.rm=TRUE), n_mor),
          fmt_n_pct(sum(oxy$any_antipruritic_final, na.rm=TRUE), n_oxy))

print_row("Naloxone reversal — n",
          as.character(sum(mor$any_naloxone_reversal_final, na.rm=TRUE)),
          as.character(sum(oxy$any_naloxone_reversal_final, na.rm=TRUE)))

print_row("Laxative — n (%)",
          fmt_n_pct(sum(mor$any_laxative, na.rm=TRUE), n_mor),
          fmt_n_pct(sum(oxy$any_laxative, na.rm=TRUE), n_oxy))

if (exists("episode_pain_summary")) {
  mor_pain <- episode_pain_summary |> filter(drug == "morphine")
  oxy_pain <- episode_pain_summary |> filter(drug == "oxycodone")
  
  print_row("Peak pain score, median (IQR) [primary]",
            fmt_med_iqr(mor_pain$peak_pain_score),
            fmt_med_iqr(oxy_pain$peak_pain_score))
  
  print_row("Time-weighted mean pain score, median (IQR)",
            fmt_med_iqr(mor_pain$time_weighted_mean_pain),
            fmt_med_iqr(oxy_pain$time_weighted_mean_pain))
}

print_row("  (secondary, SI-derived cutpoint) % shifts severe pain",
          fmt_med_iqr(mor$pct_shifts_severe),
          fmt_med_iqr(oxy$pct_shifts_severe))

print_row("First controlled shift, median (IQR)",
          fmt_med_iqr(mor$first_controlled_shift),
          fmt_med_iqr(oxy$first_controlled_shift))

print_row("Episode duration (hrs), median (IQR)",
          fmt_med_iqr(mor$episode_duration_hrs),
          fmt_med_iqr(oxy$episode_duration_hrs))

print_row("% shifts paracetamol, median (IQR)",
          fmt_med_iqr(mor$pct_shifts_paracetamol),
          fmt_med_iqr(oxy$pct_shifts_paracetamol))

print_row("% shifts NSAID, median (IQR)",
          fmt_med_iqr(mor$pct_shifts_nsaid),
          fmt_med_iqr(oxy$pct_shifts_nsaid))

if (!is.null(matched_total_dose)) {
  mor_total <- matched_total_dose |> filter(drug == "morphine")
  oxy_total <- matched_total_dose |> filter(drug == "oxycodone")
  
  print_row(sprintf("Total dose (mcg/kg/hr), median (IQR) [n=%d/%d linked]",
                    nrow(mor_total), nrow(oxy_total)),
            fmt_med_iqr(mor_total$total_dose_mcg_kg_hr, 2),
            fmt_med_iqr(oxy_total$total_dose_mcg_kg_hr, 2))
  
  print_row("  (reference) Background only, same subset",
            fmt_med_iqr(mor_total$background_mcg_kg_hr, 2),
            fmt_med_iqr(oxy_total$background_mcg_kg_hr, 2))
} else {
  cat("Total dose row skipped — episode_true_total_dose.rds unavailable\n")
}

cat("\n")
cat("============================================================\n")
cat("ANALYSIS 6 COMPLETE\n")
cat("Matched keys saved: psm_matched_keys.rds\n")
cat("Side-effect outcomes use anaesthesia-stop-anchored definitions,\n")
cat("shared with global/S4 reporting (4_flags.R Flag 9) — NOT the old\n")
cat("shift-count-based definitions. Compare this run's primary outcome\n")
cat("against the old locked OR 1.12 (p=0.58) explicitly rather than\n")
cat("assuming the conclusion carries over unchanged.\n")
cat("============================================================\n")