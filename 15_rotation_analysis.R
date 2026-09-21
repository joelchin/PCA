# ============================================================
# 15 — Rotation analysis (Table 5)
# Paediatric PCA/NCA study — Addenbrooke's Hospital (CUH NHS FT)
# Author: J Chin
# ============================================================
#
# Purpose
#   Builds the rotation panels of Table 5:
#   5a  rotation-away rate by drug (global: all episodes, all drugs)
#   5b  rotation destination (origin drug -> destination drug)
#   5c  multi-rotation summary (>= 2 genuine transitions within one admission),
#       using run-length encoding to distinguish genuine leave-and-return
#       sequences from repeated episodes on the same drug
#   5d  covariate-adjusted odds ratio for rotation away, oxycodone vs morphine,
#       from a Firth-corrected logistic model fitted on the pre-match eligible
#       pool (psm_data) with the same covariate set as the matching model. The
#       matched sample has too few events to support a model of this size.
#
# Inputs (in session)
#   pca_episodes_flagged in session; psm_data (built by 10_comparative_PSM.R, read
#   from psm_data.rds if it is not already in the session).
#
# Outputs (in session; consumed by 16_publication_tables.R)
#   tbl5 (gtsummary, 5a), ft5b, ft5c, ft5d (flextable, 5b-5d).
#
# Scope
#   5a-5c are descriptive and full-cohort. The chart-reviewed rotation reasons
#   (matched cases: side effect vs inadequate analgesia) are reported in the
#   Results text, not here. Run this script BEFORE 16_publication_tables.R.
# ============================================================

library(tidyverse)
library(gtsummary)
library(flextable)
library(logistf)

stopifnot(exists("pca_episodes_flagged"))

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

cat("============================================================\n")
cat("SCRIPT 15 - ROTATION ANALYSIS\n")
cat("============================================================\n\n")

#
# Purely descriptive — full cohort, no PSM/matched content, no
# chart-reviewed reasons (kept separate deliberately, and reported as
# prose in Results/Discussion instead). Two panels: rotation-away
# rate by drug, and destination (what drug people rotate to from
# each origin drug).
# ============================================================
cat("Building Table 5 — rotation, global descriptive...\n")

rotation_data <- pca_episodes_flagged |>
  arrange(PAT_ENC_CSN_ID, episode_start) |>
  group_by(PAT_ENC_CSN_ID) |>
  mutate(
    next_drug_rot = lead(drug),
    rotated_away_from = !is.na(next_drug_rot) & next_drug_rot != drug
  ) |>
  ungroup() |>
  mutate(
    drug = factor(drug, levels = c("morphine", "oxycodone",
                                   "morphine+ketamine", "fentanyl"),
                  labels = c("Morphine", "Oxycodone",
                             "Morphine + ketamine", "Fentanyl")),
    next_drug_rot = factor(next_drug_rot, levels = c("morphine", "oxycodone",
                                                     "morphine+ketamine", "fentanyl"),
                           labels = c("Morphine", "Oxycodone",
                                      "Morphine + ketamine", "Fentanyl"))
  )

# ---- Panel A: rotation-away rate by drug ----
tbl5_rate <- rotation_data |>
  select(drug, rotated_away_from) |>
  mutate(rotated_away_from = factor(rotated_away_from, levels = c(FALSE, TRUE),
                                    labels = c("No", "Yes"))) |>
  tbl_summary(
    by = drug,
    label = list(rotated_away_from ~ "Rotated away from this drug"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    missing = "no"
  ) |>
  add_overall(last = FALSE) |>
  modify_header(
    label = "**Outcome**",
    stat_0 = sprintf("**Overall (n=%d)**", nrow(rotation_data)),
    stat_1 = sprintf("**Morphine (n=%d)**", sum(rotation_data$drug == "Morphine")),
    stat_2 = sprintf("**Oxycodone (n=%d)**", sum(rotation_data$drug == "Oxycodone")),
    stat_3 = sprintf("**Morphine + ketamine (n=%d)**", sum(rotation_data$drug == "Morphine + ketamine")),
    stat_4 = sprintf("**Fentanyl (n=%d)**", sum(rotation_data$drug == "Fentanyl"))
  ) |>
  bold_labels()

# ---- Panel B: destination cross-tab (rotations only) ----
destination_summary <- rotation_data |>
  filter(rotated_away_from == TRUE) |>
  count(drug, next_drug_rot) |>
  group_by(drug) |>
  mutate(pct = round(100 * n / sum(n), 1),
         cell = sprintf("%d (%.1f%%)", n, pct)) |>
  ungroup() |>
  select(drug, next_drug_rot, cell) |>
  pivot_wider(names_from = next_drug_rot, values_from = cell, values_fill = "—")

cat("\nTable 5, Panel B (destination) — raw data for manual insertion:\n")
print(destination_summary)
cat("\n(Panel B built as a plain data frame above, not a gtsummary\n")
cat("object — cross-tab structure doesn't map cleanly to tbl_summary.\n")
cat("Format as a simple flextable below.)\n\n")

ft5b <- destination_summary |>
  rename(`Origin drug` = drug) |>
  flextable() |>
  ft_theme() |>
  set_caption("Table 5b. Rotation destination (origin drug -> destination drug), n (%)")

tbl5 <- tbl5_rate |>
  modify_caption(sprintf("**Table 5a.** Rotation-away rate by drug, all episodes (n=%d)",
                         nrow(rotation_data))) |>
  modify_footnote(
    all_stat_cols() ~ "n (%). 'Rotated away from this drug' indicates a later episode in the same admission used a different drug. Morphine and oxycodone rotation rates reflect switching between first-line single-agent opioids and are directly comparable. Morphine+ketamine and fentanyl involve substantially smaller episode volumes (n=229, n=210) and are not directly comparable to the two primary drugs: morphine+ketamine to morphine transitions occurred after notably longer episodes (median 142h vs 38.8h for other rotation types, p<0.001), consistent with ketamine augmentation during a prolonged, complex pain course rather than routine first-line rotation."
  )

# ---- Panel C: multi-rotation summary (corrected, RLE-based) ----
# Uses run-length encoding to collapse consecutive same-drug episodes
# before counting transitions/returns, avoiding a bug in an earlier
# version of this check that wrongly flagged simple same-drug
# repetition (e.g. "morphine, morphine, morphine") as a "return"
# pattern with no rotation ever having occurred.

multirotation_sequences <- pca_episodes_flagged |>
  arrange(PAT_ENC_CSN_ID, episode_start) |>
  group_by(PAT_ENC_CSN_ID) |>
  summarise(raw_sequence = list(drug), .groups = "drop") |>
  mutate(
    collapsed_sequence = map(raw_sequence, ~ rle(.x)$values),
    n_genuine_transitions = map_int(collapsed_sequence, length) - 1,
    has_genuine_return = map_lgl(collapsed_sequence, ~ any(duplicated(.x)))
  )

n_multi <- sum(multirotation_sequences$n_genuine_transitions >= 2)
n_return <- sum(multirotation_sequences$n_genuine_transitions >= 2 &
                  multirotation_sequences$has_genuine_return)

top_return_sequence <- multirotation_sequences |>
  filter(n_genuine_transitions >= 2, has_genuine_return) |>
  mutate(seq_str = map_chr(collapsed_sequence, ~paste(.x, collapse = " -> "))) |>
  count(seq_str, sort = TRUE) |>
  slice(1)

multirotation_summary <- tibble(
  Metric = c(
    "Admissions with >=2 genuine drug transitions",
    "  Of which, involved a return to a previously-used drug",
    sprintf("  Most common return sequence: \"%s\"", top_return_sequence$seq_str)
  ),
  Value = c(
    sprintf("%d (%.1f%% of all admissions)", n_multi, 100*n_multi/nrow(multirotation_sequences)),
    sprintf("%d (%.1f%% of multi-transition admissions)", n_return, 100*n_return/n_multi),
    sprintf("n=%d", top_return_sequence$n)
  )
)

ft5c <- multirotation_summary |>
  flextable() |>
  ft_theme() |>
  set_caption("Table 5c. Multiple rotations within a single admission") |>
  add_footer_lines(
    "Sequences collapsed by run-length encoding (consecutive episodes on the same drug counted once) before assessing transitions and returns, so repeated documentation of the same drug is not counted as a rotation. The dominant return pattern (morphine -> morphine+ketamine -> morphine, n=17) is consistent with temporary ketamine augmentation during a complex pain course followed by de-escalation back to morphine alone, rather than treatment-failure rotation (see Table 5a footnote)."
  )

cat("Table 5 (including Panel C) built\n\n")

# ============================================================
# HEADLINE SUMMARY — key top-line numbers for Results prose
# ============================================================

cat("\n")
cat("############################################################\n")
cat("# HEADLINE ROTATION NUMBERS\n")
cat("############################################################\n\n")

n_total_episodes <- nrow(pca_episodes_flagged)
n_total_rotations <- sum(rotation_data$rotated_away_from)
pct_total_rotations <- round(100 * n_total_rotations / n_total_episodes, 1)

cat(sprintf("Total episodes in cohort:                         %d\n", n_total_episodes))
cat(sprintf("Total rotation-away events:                       %d (%.1f%%)\n",
            n_total_rotations, pct_total_rotations))
cat(sprintf("Total admissions:                                 %d\n",
            nrow(multirotation_sequences)))
cat(sprintf("Admissions with >=1 genuine rotation:              %d (%.1f%%)\n",
            sum(multirotation_sequences$n_genuine_transitions >= 1),
            100*mean(multirotation_sequences$n_genuine_transitions >= 1)))
cat(sprintf("Admissions with >=2 genuine rotations:             %d (%.1f%%)\n",
            n_multi, 100*n_multi/nrow(multirotation_sequences)))
cat(sprintf("  Of those, showing a return-to-prior-drug pattern: %d (%.1f%%)\n",
            n_return, 100*n_return/n_multi))
cat(sprintf("  Most common return sequence: \"%s\" (n=%d)\n\n",
            top_return_sequence$seq_str, top_return_sequence$n))

cat("--- Ready-to-use Results sentence ---\n\n")
cat(sprintf(
  paste0(
    "\"Rotation away from the index opioid occurred in %d of %d episodes ",
    "(%.1f%%) across the full cohort. Of %d admissions, %d (%.1f%%) involved ",
    "more than one rotation; of these, %d (%.1f%%) showed a return to a ",
    "previously-used drug, most commonly %s (n=%d).\"\n\n"
  ),
  n_total_rotations, n_total_episodes, pct_total_rotations,
  nrow(multirotation_sequences),
  n_multi, 100*n_multi/nrow(multirotation_sequences),
  n_return, 100*n_return/n_multi,
  top_return_sequence$seq_str, top_return_sequence$n
))

# ============================================================
# PANEL D — COVARIATE-ADJUSTED ROTATION-AWAY OR (PRE-MATCH POOL)
# ============================================================
# Runs on the PRE-MATCH ELIGIBLE POOL (psm_data from
# 10_comparative_PSM.R), NOT the matched sample. Reasoning: the matched
# sample has too few rotation-away events in total — too sparse to trust a
# ~20-covariate model even with Firth correction.
# The pre-match pool is ~5x larger, still eligibility-restricted
# (first-line, first-episode, surgical, non-epidural — same criteria
# as matching itself) and uses the identical covariate formula, so
# confounding control is preserved without starving the model of
# events. This is the standard fallback when a matched sample's
# outcome is too sparse for within-match regression: adjust via
# regression on the full eligible pool instead of restricting further.
#
# treated == 1 is oxycodone, treated == 0 is morphine (confirmed
# from psm_data/matched_data construction in script 10) — so the OR
# below is expressed as oxycodone vs morphine (reference).

if (!exists("psm_data")) {
  psm_file <- "1 - data/3 - processed_data/psm_data.rds"
  if (!file.exists(psm_file)) stop("psm_data not found: run 10_comparative_PSM.R first.")
  psm_data <- readRDS(psm_file)
}
stopifnot(exists("rotation_data"))  # built earlier in this script

# Join the rotation-away flag onto the pre-match eligible pool.
eligible_rotation <- psm_data |>
  left_join(
    rotation_data |> select(PAT_ENC_CSN_ID, episode_id, rotated_away_from),
    by = c("PAT_ENC_CSN_ID", "episode_id")
  )

cat("Pre-match eligible episodes with rotation flag joined:", nrow(eligible_rotation), "\n")
cat("Missing rotation flag after join:", sum(is.na(eligible_rotation$rotated_away_from)), "\n")
cat("Rotation-away events in eligible pool:",
    sum(eligible_rotation$rotated_away_from, na.rm = TRUE), "\n")
cat("  — by arm:\n")
print(table(eligible_rotation$treated, eligible_rotation$rotated_away_from, useNA = "ifany"))

# Same covariate set as psm_formula_expanded, with treated (drug) as
# the predictor of interest rather than the matching target.
rotation_formula <- rotated_away_from ~ treated + AGE + WEIGHT_KG + GENDER +
  proc_category + mode + had_peripheral + year + specialty_collapsed +
  asa_grade + had_remifentanil + had_intraop_ketamine + had_intraop_clonidine +
  neuraxial_intrathecal_mar + had_intraop_paracetamol + had_intraop_nsaid +
  had_intraop_magnesium + had_intraop_fentanyl + had_intraop_morphine +
  had_intraop_oxycodone2 + had_intraop_diamorphine + anaesthesia_duration_mins

rotation_model <- logistf(
  rotation_formula,
  data = eligible_rotation,
  control = logistf.control(maxit = 1000, maxstep = 5)
)

cat("\n============================================================\n")
cat("FIRTH-CORRECTED ROTATION-AWAY MODEL — PRE-MATCH ELIGIBLE POOL\n")
cat("============================================================\n")
print(summary(rotation_model))

# Pull the treated coefficient specifically and exponentiate for the
# OR + 95% CI — treated==1 (oxycodone) vs treated==0 (morphine, ref).
treated_idx <- which(names(coef(rotation_model)) == "treated")
or_estimate <- exp(coef(rotation_model)[treated_idx])
or_ci <- exp(confint(rotation_model)[treated_idx, ])

cat(sprintf(
  "\nAdjusted OR for rotation-away, oxycodone vs morphine (pre-match eligible pool, n=%d):\n",
  nrow(eligible_rotation)
))
cat(sprintf("  OR = %.3f (95%% CI %.3f-%.3f)\n", or_estimate, or_ci[1], or_ci[2]))
cat(sprintf("  p  = %.4f\n", rotation_model$prob[treated_idx]))

# ============================================================
# BUILD ft5d — flextable object for 16_publication_tables.R
# ============================================================
# A single-row summary table, formatted consistently with the other
# Table 5 panels (as_flex_table |> ft_theme() pattern used throughout
# 16_publication_tables.R).

tbl5d_data <- tibble::tibble(
  Comparison = "Oxycodone vs Morphine (reference)",
  `Adjusted OR` = sprintf("%.3f", or_estimate),
  `95% CI` = sprintf("%.3f\u2013%.3f", or_ci[1], or_ci[2]),
  `p-value` = ifelse(rotation_model$prob[treated_idx] < 0.001,
                     "<0.001",
                     sprintf("%.3f", rotation_model$prob[treated_idx])),
  Population = sprintf("Pre-match eligible pool (n=%d)", nrow(eligible_rotation))
)

ft5d <- flextable::flextable(tbl5d_data) |>
  ft_theme() |>
  flextable::set_caption(
    caption = paste0(
      "Table 5d. Covariate-adjusted rotation-away odds ratio ",
      "(Firth-corrected logistic regression, same covariates as ",
      "propensity-score matching). Estimated on the pre-match ",
      "eligible pool rather than the matched sample, due to too few ",
      "rotation-away events within the matched sample alone to ",
      "support this covariate set."
    )
  )

cat("\nft5d built — ready for 16_publication_tables.R\n\n")

# The adjusted OR describes the pre-match eligible pool, not the matched
# pairs. State this explicitly in Methods/Results so readers know which
# population it relates to relative to Tables 3-4.

cat("Script 15 complete - tbl5, ft5b, ft5c, ft5d ready for 16_publication_tables.R\n\n")
