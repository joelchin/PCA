# ============================================================
# 15 - Opioid rotation analysis (Table 4)
# ============================================================
# Builds the rotation table of the manuscript:
#   tbl4a     (A) rotation-away rate by drug, all episodes
#   ft4b      (B) destination of rotation (origin drug -> destination drug)
#   ft4c      (C) covariate-adjusted odds ratio of rotation-away,
#             oxycodone vs morphine, from Firth-corrected logistic
#             regression on the pre-match eligible pool (psm_data): the
#             matched sample has too few events to support the full
#             covariate set
#   ft4_multi multiple rotations within an admission (not a manuscript
#             table; its figures are quoted in the Results text)
#
# Rotation-away means a later episode in the same admission used a
# different drug. The chart-reviewed reasons for rotation are reported
# in the Results text, not here.
#
# Run after 09_comparative_psm.R in the same R session (needs psm_data
# and matched_data) and before 12_publication_tables.R. Every number in
# the captions and footnotes is computed from the data.
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

# ============================================================
# Panels A and B: all episodes, descriptive
# ============================================================

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

# Panel A: rotation-away rate by drug
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

# Panel B: destination (rotations only)
destination_summary <- rotation_data |>
  filter(rotated_away_from == TRUE) |>
  count(drug, next_drug_rot) |>
  group_by(drug) |>
  mutate(pct = round(100 * n / sum(n), 1),
         cell = sprintf("%d (%.1f%%)", n, pct)) |>
  ungroup() |>
  select(drug, next_drug_rot, cell) |>
  pivot_wider(names_from = next_drug_rot, values_from = cell, values_fill = "—")

cat("\nPanel B (destination):\n")
print(destination_summary)
cat("\n")

ft4b <- destination_summary |>
  rename(`Origin drug` = drug) |>
  flextable() |>
  ft_theme() |>
  set_caption("Table 4B. Rotation destination (origin drug -> destination drug), n (%)")

# Numbers for the Panel A footnote
n_mk  <- sum(rotation_data$drug == "Morphine + ketamine")
n_fen <- sum(rotation_data$drug == "Fentanyl")

# Duration of the origin episode: morphine+ketamine -> morphine
# transitions vs all other rotation types
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

tbl4a <- tbl5_rate |>
  modify_caption(sprintf("**Table 4A.** Rotation-away rate by drug, all episodes (n=%d)",
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

# Multiple rotations within an admission. Run-length encoding collapses
# consecutive episodes on the same drug before counting transitions and
# returns, so repeated documentation of one drug is not a rotation.

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

ft4_multi <- multirotation_summary |>
  flextable() |>
  ft_theme() |>
  set_caption("Multiple rotations within a single admission") |>
  add_footer_lines(
    sprintf(
      "Sequences collapsed by run-length encoding (consecutive episodes on the same drug counted once) before assessing transitions and returns. The most common return pattern was %s (n=%d).",
      top_return_sequence$seq_str, top_return_sequence$n)
  )

cat("Panels A and B and the multiple-rotation summary built\n\n")

# ============================================================
# Headline rotation numbers (quoted in the Results text)
# ============================================================

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
cat(sprintf("  Of those, with a return to a prior drug:          %d (%.1f%%)\n",
            n_return, 100*n_return/n_multi))
cat(sprintf("  Most common return sequence: \"%s\" (n=%d)\n\n",
            top_return_sequence$seq_str, top_return_sequence$n))

# ============================================================
# Panel C: covariate-adjusted rotation-away odds ratio
# ============================================================
# Runs on the pre-match eligible pool (psm_data from 09_comparative_psm.R),
# not the matched sample: the matched sample has roughly one rotation
# event per covariate, too few for a model with this covariate set even
# with Firth correction. The pool uses the same eligibility criteria and
# the same covariates as the matching model. treated == 1 is oxycodone,
# treated == 0 is morphine (reference).

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
cat("  by arm:\n")
print(table(eligible_rotation$treated, eligible_rotation$rotated_away_from, useNA = "ifany"))

# Same covariates as the matching model, with treated (drug) as the
# predictor of interest
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
cat("FIRTH-CORRECTED ROTATION-AWAY MODEL: PRE-MATCH ELIGIBLE POOL\n")
cat("============================================================\n")
print(summary(rotation_model))

# Odds ratio and 95% CI for treated (oxycodone vs morphine)
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
# Matched pair and event counts for the Panel C caption
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
# Panel C table
# ============================================================

tbl5d_data <- tibble::tibble(
  Comparison = "Oxycodone vs Morphine (reference)",
  `Adjusted OR` = sprintf("%.3f", or_estimate),
  `95% CI` = sprintf("%.3f\u2013%.3f", or_ci[1], or_ci[2]),
  `p-value` = ifelse(rotation_model$prob[treated_idx] < 0.001,
                     "<0.001",
                     sprintf("%.3f", rotation_model$prob[treated_idx])),
  Population = sprintf("Pre-match eligible pool (n=%d)", nrow(eligible_rotation))
)

ft4c <- flextable::flextable(tbl5d_data) |>
  ft_theme() |>
  flextable::set_caption(
    caption = paste0(
      "Table 4C. Covariate-adjusted rotation-away odds ratio ",
      "(Firth-corrected logistic regression, same covariates as ",
      "propensity-score matching). Estimated on the pre-match ",
      "eligible pool rather than ", cap_matched_txt
    )
  )

cat("\nPanel C table built\n\n")

cat("Script 11 complete: tbl4a, ft4b, ft4_multi and ft4c are ready for 12_publication_tables.R\n")
