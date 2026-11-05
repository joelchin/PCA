# 15 — Rotation Analysis
# Paediatric PCA Study
# Addenbrooke's Hospital
#
# Extracted from 16_publication_tables.R into its own script, since
# rotation analysis is conceptually distinct from Tables 1-4 (all
# PSM/descriptive-population content) and had grown large enough to
# deserve separation.
#
# Builds Table 5 in four panels:
#   5a — rotation-away rate by drug, global (all episodes, all drugs)
#   5b — rotation destination (origin drug -> destination drug)
#   5c — multi-rotation summary (>=2 genuine transitions within a
#        single admission), using run-length encoding to correctly
#        distinguish genuine leave-and-return sequences from trivial
#        same-drug repetition (see script comments below for detail)
#   5d — covariate-adjusted rotation-away odds ratio (oxycodone vs
#        morphine), computed on the pre-match eligible pool
#        (psm_data) rather than the matched pairs, due to too few
#        events in the matched sample alone to support a ~20-covariate
#        model (the matched pair count and event count are computed
#        from matched_data and printed in the caption, not hardcoded)
#
# Produces: tbl5 (gtsummary object, panel 5a), ft5b, ft5c, ft5d
# (flextable objects, panels 5b/5c/5d) — consumed by
# 16_publication_tables.R for the Word export. Run this BEFORE
# 16_publication_tables.R, and AFTER 09_comparative_psm.R in the same
# session (it needs psm_data and matched_data).
#
# CHANGED (Oct 2026): every number quoted in the table footnotes and
# captions (episode counts, ketamine-to-morphine durations, the top
# return sequence and its n, matched pair and event counts) is now
# computed live from the data in this script, not typed in. Previously
# several were hardcoded and could drift from the data. Neonates are
# no longer excluded upstream (script 10), so the eligible pool used
# for Panel D includes any eligible neonates.
#
# Deliberately excludes: the chart-reviewed rotation REASONS, which
# are reported in the Results prose, not here.
# ============================================================

library(tidyverse)
library(gtsummary)
library(flextable)

stopifnot(exists("pca_episodes_flagged"))

ft_theme <- function(ft) {
  ft |>
    theme_booktabs() |>
    fontsize(size = 10, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |>   # qualified: ggpubr/survminer also define font()
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
# prose in Results/Discussion instead). Panels: rotation-away rate by
# drug, and destination (what drug people rotate to from each origin
# drug).
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

# ---- Live numbers for the Table 5a footnote ----------------------
# (previously typed into the footnote; now computed from the data)
n_mk  <- sum(rotation_data$drug == "Morphine + ketamine")
n_fen <- sum(rotation_data$drug == "Fentanyl")

# Duration of the origin episode, morphine+ketamine -> morphine
# transitions vs all other rotation types. Verify that this
# reproduces the figures quoted in the manuscript (142 h vs 38.8 h,
# p<0.001); if it does not, the manuscript must follow this output.
rot_dur <- rotation_data |>
  filter(rotated_away_from, !is.na(episode_duration_hrs)) |>
  mutate(mk_to_morphine = drug == "Morphine + ketamine" & next_drug_rot == "Morphine")

cat("--- Origin-episode duration: morphine+ketamine -> morphine vs other rotations ---\n")
rot_dur_summary <- rot_dur |>
  group_by(mk_to_morphine) |>
  summarise(n = n(),
            median_hrs = round(median(episode_duration_hrs), 1),
            q25 = round(quantile(episode_duration_hrs, 0.25), 1),
            q75 = round(quantile(episode_duration_hrs, 0.75), 1),
            .groups = "drop")
print(rot_dur_summary)

if (any(rot_dur$mk_to_morphine) && any(!rot_dur$mk_to_morphine)) {
  mk_med  <- median(rot_dur$episode_duration_hrs[rot_dur$mk_to_morphine])
  oth_med <- median(rot_dur$episode_duration_hrs[!rot_dur$mk_to_morphine])
  mk_p    <- wilcox.test(episode_duration_hrs ~ mk_to_morphine, data = rot_dur)$p.value
  mk_p_txt <- if (mk_p < 0.001) "p<0.001" else sprintf("p=%.3f", mk_p)
  cat(sprintf("Median %s h vs %s h, Wilcoxon %s\n\n",
              format(round(mk_med, 1)), format(round(oth_med, 1)), mk_p_txt))
  mk_footnote_txt <- sprintf(
    "morphine+ketamine to morphine transitions occurred after notably longer episodes (median %s h vs %s h for other rotation types, %s), consistent with ketamine augmentation during a prolonged, complex pain course rather than routine first-line rotation.",
    format(round(mk_med, 1)), format(round(oth_med, 1)), mk_p_txt)
} else {
  cat("No morphine+ketamine -> morphine transitions (or no comparators) found.\n\n")
  mk_footnote_txt <- "(morphine+ketamine to morphine duration comparison not available in this run)."
}

tbl5 <- tbl5_rate |>
  modify_caption(sprintf("**Table 5a.** Rotation-away rate by drug, all episodes (n=%d)",
                         nrow(rotation_data))) |>
  modify_footnote(
    all_stat_cols() ~ paste0(
      "n (%). 'Rotated away from this drug' indicates a later episode in the same admission used a different drug. ",
      "Morphine and oxycodone rotation rates reflect switching between first-line single-agent opioids and are directly comparable. ",
      sprintf("Morphine+ketamine and fentanyl involve substantially smaller episode volumes (n=%d, n=%d) and are not directly comparable to the two primary drugs: ",
              n_mk, n_fen),
      mk_footnote_txt
    )
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
    sprintf(
      "Sequences collapsed by run-length encoding (consecutive episodes on the same drug counted once) before assessing transitions and returns, so repeated documentation of the same drug is not counted as a rotation. The most common return pattern in this run was %s (n=%d); if this is morphine -> morphine+ketamine -> morphine it is consistent with temporary ketamine augmentation during a complex pain course followed by de-escalation back to morphine alone, rather than treatment-failure rotation (see Table 5a footnote).",
      top_return_sequence$seq_str, top_return_sequence$n)
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
# Added Oct 2026. Runs on the PRE-MATCH ELIGIBLE POOL (psm_data from
# 09_comparative_psm.R; pool sizes are printed by script 10 and below),
# NOT the matched sample. Reasoning: the matched sample has too few
# rotation-away events to trust a ~20-covariate model even with Firth
# correction (roughly 1 event per covariate). The pre-match pool is
# ~5x larger, still eligibility-restricted (same criteria as matching
# itself) and uses the identical covariate formula, so confounding
# control is preserved without starving the model of events. This is
# the standard fallback when a matched sample's outcome is too sparse
# for within-match regression: adjust via regression on the full
# eligible pool instead of restricting further.
#
# Neonates are no longer excluded upstream, so this pool includes any
# eligible neonates.
#
# treated == 1 is oxycodone, treated == 0 is morphine (confirmed
# from psm_data/matched_data construction in script 10) — so the OR
# below is expressed as oxycodone vs morphine (reference).

library(logistf)

stopifnot(exists("psm_data"))
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

t_start <- Sys.time()
rotation_model <- logistf(
  rotation_formula,
  data = eligible_rotation,
  control = logistf.control(maxit = 1000, maxstep = 5),
  plcontrol = logistpl.control(maxit = 1000),
  plconf = 2          # profile-likelihood CI for coefficient 2 only (treated); 1 is the intercept
)
stopifnot(names(coef(rotation_model))[2] == "treated")
cat("logistf (full model) took", round(difftime(Sys.time(), t_start, units = "mins"), 1), "min\n")

cat("\n============================================================\n")
cat("FIRTH-CORRECTED ROTATION-AWAY MODEL — PRE-MATCH ELIGIBLE POOL\n")
cat("============================================================\n")
print(summary(rotation_model))

# Pull the treated coefficient specifically and exponentiate for the
# OR + 95% CI — treated==1 (oxycodone) vs treated==0 (morphine, ref).
treated_idx <- 2
or_estimate <- exp(coef(rotation_model)[treated_idx])
or_ci <- exp(c(rotation_model$ci.lower[treated_idx], rotation_model$ci.upper[treated_idx]))

cat(sprintf(
  "\nAdjusted OR for rotation-away, oxycodone vs morphine (pre-match eligible pool, n=%d):\n",
  nrow(eligible_rotation)
))
cat(sprintf("  OR = %.3f (95%% CI %.3f-%.3f)\n", or_estimate, or_ci[1], or_ci[2]))
cat(sprintf("  p  = %.4f\n", rotation_model$prob[treated_idx]))

# ============================================================
# MATCHED-SAMPLE PAIR AND EVENT COUNTS FOR THE TABLE 5d CAPTION
# Computed from matched_data (script 10), not hardcoded.
# ============================================================
if (exists("matched_data")) {
  n_pairs_cap <- sum(matched_data$treated == 1)
  events_matched_cap <- rotation_data |>
    semi_join(matched_data, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    summarise(n = sum(rotated_away_from, na.rm = TRUE)) |>
    pull(n)
  cat(sprintf("\nMatched sample: %d pairs (%d episodes), %d rotation-away events\n",
              n_pairs_cap, 2 * n_pairs_cap, events_matched_cap))
  cap_matched_txt <- sprintf(
    "the matched sample (%d pairs), which had too few rotation-away events (n=%d) to support this covariate set.",
    n_pairs_cap, events_matched_cap
  )
} else {
  cat("\nmatched_data not found in the session: caption will not quote pair/event counts.\n")
  cap_matched_txt <- paste0(
    "the matched sample, which had too few rotation-away events ",
    "to support this covariate set."
  )
}

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
      "eligible pool rather than ", cap_matched_txt
    )
  )

cat("\nft5d built — ready for 16_publication_tables.R\n\n")

# The figure printed by this run is the one to cite. It supersedes the
# earlier candidate figures that appeared in draft material (0.124 and
# ~0.19, and 0.084 from the run before neonates were no longer
# excluded), which were computed under different, no-longer-current
# formulations and should not be cited. State explicitly in
# Methods/Results that this OR describes the pre-match eligible pool,
# not the matched pairs, so a reader isn't confused about which
# population this OR relates to relative to Table 3's matched sample.

cat("Script 15 complete - tbl5, ft5b, ft5c, ft5d ready for 16_publication_tables.R\n\n")