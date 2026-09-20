# ============================================================
# 13 — Section 5: Subgroup Analyses
# Paediatric PCA Study
# Addenbrooke's Hospital
#
# Uses the finalised, anaesthesia-stop-anchored side-effect flags
# from 4_flags.R (Flag 9) throughout: any_side_effect_final,
# reactive_antiemetic_final, any_antipruritic_final,
# any_naloxone_reversal_final.
#
# any_laxative_full is used rather than any_laxative: the two are
# genuinely different columns in episode_metrics (1200 vs 1147
# events, only 5321/5374 episodes agree). any_laxative_full is the
# true full-episode, no-time-filter flag; what any_laxative (plain)
# represents is unclear and is not a substitute for it.
#
# matched_keys should be re-derived from a fresh run of
# 10_comparative_PSM.R before this script is trusted, so the
# antiemetic definition it relies on is current.
# ============================================================

library(tidyverse)
library(lubridate)

# ============================================================
# LOAD DATA
# ============================================================

stopifnot(
  exists("pca_episodes_flagged"),
  exists("episode_metrics")
)

matched_keys <- tryCatch(
  readRDS("1 - data/3 - processed_data/psm_matched_keys.rds"),
  error = function(e) {
    stop("psm_matched_keys.rds not found — run 10_comparative_PSM.R first")
  }
)

cat("============================================================\n")
cat("ANALYSIS 8 — SUBGROUP ANALYSES (SECTION 5)\n")
cat("Population: PSM matched pairs from 10_comparative_PSM.R\n")
cat("============================================================\n\n")

fmt_n_pct <- function(n, total) {
  sprintf("%d (%.1f%%)", n, 100 * n / total)
}

fmt_med_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  sprintf("%.1f (%.1f–%.1f)",
          round(median(x), digits),
          round(quantile(x, 0.25), digits),
          round(quantile(x, 0.75), digits))
}

# ============================================================
# BUILD MATCHED ANALYSIS DATASET — INCLUDES FINALISED FLAGS
# (via pca_episodes_flagged, first join — no separate join needed)
# ============================================================

cat("--- Building matched subgroup dataset ---\n\n")

matched_data <- pca_episodes_flagged |>
  semi_join(matched_keys, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  left_join(
    matched_keys |> select(PAT_ENC_CSN_ID, episode_id,
                           treated, subclass, distance),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  ) |>
  left_join(
    episode_metrics |> select(-drug, -mode),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  )

stopifnot(
  "reactive_antiemetic_final" %in% names(matched_data),
  "any_antipruritic_final" %in% names(matched_data),
  "any_naloxone_reversal_final" %in% names(matched_data),
  "any_side_effect_final" %in% names(matched_data),
  "any_laxative_full" %in% names(matched_data)
)

cat("Matched episodes loaded:", nrow(matched_data), "\n")
cat("Oxycodone:", sum(matched_data$drug == "oxycodone"), "\n")
cat("Morphine: ", sum(matched_data$drug == "morphine"), "\n\n")

# ============================================================
# SUBGROUP COMPARISON FUNCTION (FINAL definitions)
# ============================================================

subgroup_compare <- function(data, subgroup_var, label) {
  
  cat(sprintf("--- %s ---\n", label))
  
  levels <- unique(data[[subgroup_var]])
  levels <- levels[!is.na(levels)]
  
  results <- map_dfr(levels, function(lv) {
    sub <- data |> filter(.data[[subgroup_var]] == lv)
    
    mor <- sub |> filter(drug == "morphine")
    oxy <- sub |> filter(drug == "oxycodone")
    
    n_mor <- nrow(mor)
    n_oxy <- nrow(oxy)
    
    if (n_mor < 5 || n_oxy < 5) {
      return(tibble(
        subgroup      = as.character(lv),
        n_morphine    = n_mor,
        n_oxycodone   = n_oxy,
        pct_se_mor    = NA_real_,
        pct_se_oxy    = NA_real_,
        OR            = NA_real_,
        p_value       = NA_real_,
        note          = "insufficient n"
      ))
    }
    
    pct_mor <- round(100 * mean(mor$any_side_effect_final, na.rm = TRUE), 1)
    pct_oxy <- round(100 * mean(oxy$any_side_effect_final, na.rm = TRUE), 1)
    
    tab <- table(sub$drug, sub$any_side_effect_final)
    
    if (ncol(tab) < 2) {
      return(tibble(
        subgroup    = as.character(lv),
        n_morphine  = n_mor,
        n_oxycodone = n_oxy,
        pct_se_mor  = pct_mor,
        pct_se_oxy  = pct_oxy,
        OR          = NA_real_,
        p_value     = NA_real_,
        note        = "no variation"
      ))
    }
    
    ct <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
    or <- tryCatch(
      (tab["oxycodone", 2] / tab["oxycodone", 1]) /
        (tab["morphine",  2] / tab["morphine",  1]),
      error = function(e) NA_real_
    )
    
    tibble(
      subgroup    = as.character(lv),
      n_morphine  = n_mor,
      n_oxycodone = n_oxy,
      pct_se_mor  = pct_mor,
      pct_se_oxy  = pct_oxy,
      OR          = round(or, 2),
      p_value     = if (!is.null(ct)) round(ct$p.value, 4) else NA_real_,
      note        = ""
    )
  })
  
  print(results)
  cat("\n")
  invisible(results)
}

# ============================================================
# TABLE 1 — SUBGROUP: AGE GROUP
# ============================================================

cat("============================================================\n")
cat("TABLE 1 — SUBGROUP: AGE GROUP\n")
cat("Primary outcome: composite side effect (FINAL)\n")
cat("============================================================\n\n")

subgroup_compare(matched_data, "age_group", "Age group")

cat("--- Antiemetic rate by age group (FINAL) ---\n")
matched_data |>
  group_by(age_group, drug) |>
  summarise(
    n              = n(),
    pct_antiemetic = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = drug, values_from = c(n, pct_antiemetic)) |>
  print()
cat("\n")

# ============================================================
# TABLE 2 — SUBGROUP: MODE (PCA vs NCA)
# ============================================================

cat("============================================================\n")
cat("TABLE 2 — SUBGROUP: MODE (PCA vs NCA)\n")
cat("============================================================\n\n")

subgroup_compare(matched_data, "mode", "Mode")

cat("--- Pain outcomes by mode (matched, all ages) ---\n")
matched_data |>
  group_by(mode, drug) |>
  summarise(
    n                 = n(),
    pct_severe_median = round(median(pct_shifts_severe, na.rm = TRUE), 1),
    pct_severe_q25    = round(quantile(pct_shifts_severe, 0.25, na.rm = TRUE), 1),
    pct_severe_q75    = round(quantile(pct_shifts_severe, 0.75, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()
cat("\n")

# ============================================================
# TABLE 3 — SUBGROUP: PROCEDURE CATEGORY
# ============================================================

cat("============================================================\n")
cat("TABLE 3 — SUBGROUP: PROCEDURE CATEGORY\n")
cat("(Restricted to categories with n≥10 per drug)\n")
cat("============================================================\n\n")

subgroup_compare(matched_data, "proc_category", "Procedure category")

# ============================================================
# TABLE 4 — SUBGROUP: ERA
# ============================================================

cat("============================================================\n")
cat("TABLE 4 — SUBGROUP: ERA\n")
cat("(Early 2014-2018 vs Recent 2019-2026)\n")
cat("============================================================\n\n")

subgroup_compare(matched_data, "era", "Era")

cat("--- Secondary outcomes by era (matched, FINAL) ---\n")
matched_data |>
  group_by(era, drug) |>
  summarise(
    n                = n(),
    pct_antiemetic   = round(100 * mean(reactive_antiemetic_final, na.rm = TRUE), 1),
    pct_antipruritic = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
    pct_laxative     = round(100 * mean(any_laxative_full, na.rm = TRUE), 1),
    n_nalox_reversal = sum(any_naloxone_reversal_final, na.rm = TRUE),
    .groups = "drop"
  ) |>
  print()
cat("\n")

# ============================================================
# TABLE 5 — SUBGROUP: FIRST-LINE VS ROTATION
# ============================================================

cat("============================================================\n")
cat("TABLE 5 — SUBGROUP: FIRST-LINE VS ROTATION\n")
cat("============================================================\n\n")

subgroup_compare(matched_data, "is_first_line", "First-line vs rotation")

cat("Note: PSM population is all first-line by inclusion criterion\n")
cat("(is_first_line == TRUE required for S3 matching pool)\n")
cat("Rotation subgroup comparison not applicable in matched population\n\n")

# ============================================================
# TABLE 6 — SUBGROUP: PERIPHERAL NERVE BLOCK STATUS
# ============================================================

cat("============================================================\n")
cat("TABLE 6 — SUBGROUP: PERIPHERAL NERVE BLOCK STATUS\n")
cat("============================================================\n\n")

subgroup_compare(matched_data, "had_peripheral", "Peripheral nerve block")

# ============================================================
# TABLE 7 — SUBGROUP: INTRAOP SPINAL
# ============================================================

cat("============================================================\n")
cat("TABLE 7 — SUBGROUP: INTRAOPERATIVE SPINAL\n")
cat("(neuraxial_intrathecal_mar — PSM covariate)\n")
cat("============================================================\n\n")

subgroup_compare(matched_data, "neuraxial_intrathecal_mar", "Intraop spinal")

# ============================================================
# TABLE 8 — LAXATIVE RATES, DESCRIPTIVE ONLY (any_laxative_full)
# Not a formal comparative outcome — see caveat below
# ============================================================

cat("============================================================\n")
cat("TABLE 8 — LAXATIVE RATES (matched pairs)\n")
cat("(Opioid-induced constipation proxy)\n")
cat("============================================================\n\n")

cat("--- Overall laxative rates (matched) ---\n")
matched_data |>
  group_by(drug) |>
  summarise(
    n            = n(),
    pct_laxative = round(100 * mean(any_laxative_full, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()
cat("\n")

cat("(Laxative rate DESCRIPTIVE ONLY — not tested statistically.\n")
cat(" Cannot distinguish prophylactic bowel-regimen dosing from\n")
cat(" genuine reactive treatment of opioid-induced constipation;\n")
cat(" no order-type field exists to separate these.)\n\n")

cat("--- Laxative rate by age group (matched) ---\n")
matched_data |>
  group_by(age_group, drug) |>
  summarise(
    n            = n(),
    pct_laxative = round(100 * mean(any_laxative_full, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = drug, values_from = c(n, pct_laxative)) |>
  print()
cat("\n")

cat("--- Laxative rate by mode (matched) ---\n")
matched_data |>
  group_by(mode, drug) |>
  summarise(
    n            = n(),
    pct_laxative = round(100 * mean(any_laxative_full, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = drug, values_from = c(n, pct_laxative)) |>
  print()
cat("\n")

# ============================================================
# SUMMARY — FOREST PLOT DATA (FINAL definitions)
# ============================================================

cat("============================================================\n")
cat("FOREST PLOT DATA — SUBGROUP ORs (FINAL)\n")
cat("(For manuscript Figure 2)\n")
cat("============================================================\n\n")

get_or_ci <- function(data, subgroup_label) {
  tab <- table(data$drug, data$any_side_effect_final)
  if (ncol(tab) < 2 || any(tab == 0)) return(NULL)
  ct <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  or <- (tab["oxycodone", 2] / tab["oxycodone", 1]) /
    (tab["morphine",  2] / tab["morphine",  1])
  log_or <- log(or)
  se_log_or <- sqrt(sum(1 / tab))
  ci_lo <- exp(log_or - 1.96 * se_log_or)
  ci_hi <- exp(log_or + 1.96 * se_log_or)
  tibble(
    subgroup = subgroup_label,
    n_oxy    = sum(data$drug == "oxycodone"),
    n_mor    = sum(data$drug == "morphine"),
    OR       = round(or, 2),
    CI_lo    = round(ci_lo, 2),
    CI_hi    = round(ci_hi, 2),
    p        = round(ct$p.value, 4)
  )
}

forest_data <- bind_rows(
  get_or_ci(matched_data, "Overall"),
  get_or_ci(matched_data |> filter(mode == "PCA"),  "PCA"),
  get_or_ci(matched_data |> filter(mode == "NCA"),  "NCA"),
  get_or_ci(matched_data |> filter(age_group == "Adolescent (12-15yr)"), "Adolescent"),
  get_or_ci(matched_data |> filter(age_group == "Child (6-11yr)"),       "Child 6-11yr"),
  get_or_ci(matched_data |> filter(age_group == "Young child (1-5yr)"),  "Young child 1-5yr"),
  get_or_ci(matched_data |> filter(is_first_line == TRUE),  "First line"),
  get_or_ci(matched_data |> filter(is_first_line == FALSE), "Rotation"),
  get_or_ci(matched_data |> filter(era == "Early (2014-2018)"),   "Early era"),
  get_or_ci(matched_data |> filter(era == "Recent (2019-2026)"),  "Recent era"),
  get_or_ci(matched_data |> filter(proc_category == "Orthopaedic/spine"),    "Orthopaedic/spine"),
  get_or_ci(matched_data |> filter(proc_category == "Neurosurgery"),          "Neurosurgery"),
  get_or_ci(matched_data |> filter(proc_category == "Colorectal"),            "Colorectal"),
  get_or_ci(matched_data |> filter(proc_category == "Skin/trauma/general"),   "Skin/trauma/general")
)

print(forest_data)
cat("\n")

saveRDS(forest_data, "1 - data/3 - processed_data/forest_plot_data.rds")
cat("Forest plot data saved\n\n")

cat("============================================================\n")
cat("ANALYSIS 8 COMPLETE\n")
cat("Pre-specified subgroup analyses complete (FINAL definitions)\n")
cat("Forest plot data saved for Figure 2\n")
cat("Note: subgroup analyses hypothesis-generating only\n")
cat("Not powered for formal subgroup comparisons\n")
cat("============================================================\n")