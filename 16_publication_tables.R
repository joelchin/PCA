# 16 — Publication Tables
# Paediatric PCA Study
# Addenbrooke's Hospital
#
# Builds the manuscript's publication tables as flextable objects ready
# for Word export (Table 1: cohort characteristics; Table 2: oxycodone
# descriptive; Table 3: PSM balance; Table 4: PSM-matched comparative
# outcomes; Table 5: rotation).
#
# REVISION NOTES (this version)
#   Table 1  now reports line of therapy (1st/2nd/3rd/4th+, counted as
#            distinct drug changes WITHIN an admission, matching the
#            is_first_line flag), episodes per patient in the patient-level
#            section, and a footnote stating that the global tables are
#            descriptive and that comparative inference is on the matched
#            (first-ever PCA/NCA episode, one-per-patient) cohort.
#   Table 2  is now stacked like Table 1: patient-level section (one row per
#            patient: age, weight, sex, age group) and episode-level section
#            (mode, line of therapy, duration, procedure category, side
#            effects, pain, consumption). The original by-PCA/NCA episode
#            table is still built as tbl2_by_mode but is NOT exported; add
#            it to the Word export if you want it.
#   Howard   comparison applies ONLY Howard et al.'s stated criteria
#            (morphine NCA after surgery, not a switch-in). Console summary
#            only (informal anchor, not validation, not exported).
#            Sensitivity pools: HOWARD_FAITHFUL.R.
#   Table 4  footnote now states the population explicitly.
#   Export   publication-style formatting (Arial 9 pt, three-line tables,
#            8 pt footnotes, A4 narrow margins, page break between tables),
#            stat-column headers that state the statistic instead of "N =",
#            the stacked-table "Group" header blanked.
#            Footnotes that cannot be attached to a stacked gtsummary table
#            are written as text under the table instead.
#
#   EARLIER PASS (first-exposure cohort, medical oncology excluded)
#   - Table S1 (tbl3): post-match specialty is now collapsed using the
#     PRE-match specialty list (as in the matching model). Before, the
#     n>=30 collapse was re-run inside the matched rows, which folded
#     Urology and PAEDIATRICS into "Other specialty" and showed them as
#     0 (0%).
#   - Table S1 footnote: the balance statement is computed from match_out
#     (maximum SMD, number of terms above 0.10, which terms) instead of a
#     hardcoded "<0.10". match_out must be in the session (script 10).
#   - Year is shown without thousands separator or decimals in Table S1.
#   - Table 1, Table 3 (tbl4) and Notes wording now describe the matched
#     cohort as each child's first-ever recorded PCA/NCA episode at this
#     centre, medical oncology excluded.
#   - Table 3 (tbl4) footnote: matched n and the flowsheet-linkage p-value
#     are computed from the data, not hardcoded.
#   - Study period in the title set to October 2014 to July 2026; the
#     data range is printed to the console as a check.
#
#   THIS PASS (Oct 2026)
#   - Footnotes shortened to match the manuscript (Tables 1, 2, 3 and S1):
#     each now says what is needed to read the table and nothing more;
#     rationale lives in the Methods.
#   - Table 2 footnote: the number of excluded oxycodone episodes and the
#     reasons (no recorded episode end time; undetermined PCA/NCA mode) are
#     computed from the data. Previously "78 episodes with undetermined
#     mode" and "1,276" were typed in; 78 is the whole-cohort unknown-mode
#     count, not the oxycodone one.
#   - Table 3 footnote: the sentence comparing matched naloxone counts with
#     a full-cohort chart-reviewed figure is removed. The duplicate
#     gtsummary footnote on the p-value column (test names) is removed,
#     since the tests are named in the main footnote.
#   - s4_surgical_keys (script 6) is checked at the top, so a missing
#     object fails at the start rather than after the long build.
#   - Neonates are no longer excluded upstream (script 10); nothing in
#     this script depended on that.
#
# Methodological choices carried over from earlier versions:
#
#   Table 1's case type field is case_type_refined, not the original
#   unrefined case_type.
#
#   Table 2 reports genuine efficacy metrics (peak pain score, time-weighted
#   mean pain, consumption) rather than the SI project's osteotomy-derived
#   metrics. Its procedure-category breakdown uses proc_category_display
#   (surgical categories individually, non-surgical collapsed to one row).
#   The composite side-effect outcome is not shown; paracetamol/NSAID
#   shift rows are not included (paracetamol co-prescribing ~100%).
#
#   Table 4 reports antiemetic, antipruritic and naloxone reversal as
#   independent outcomes (no composite).
#
#   pre_match, post_match and table4_matched_outcomes are built directly
#   from the live pipeline objects (psm_data, matched_data) rather than
#   reconstructed from saved keys, keeping population figures mechanically
#   identical to 10_comparative_PSM.R. add_psm_covariates() is applied only
#   to post_match (psm_data already carries those covariates).
#
#   Table 5 draws on four objects from 15_rotation_analysis.R (tbl5, ft5b,
#   ft5c, ft5d). ft5d is the covariate-adjusted rotation-away odds ratio,
#   computed on the pre-match eligible pool.
#
#   case_type_refined is reported at patient level in Table 1.
#
#   Naloxone-reversal factor() calls specify levels = c(FALSE, TRUE)
#   explicitly because oxycodone's S4 population has zero reversal events.
# ============================================================
required_pkgs <- c("gtsummary", "flextable", "officer", "labelled")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "), ".\n",
       "Install them in a FRESH R session (Session > Restart R) BEFORE running the pipeline, ",
       "not mid-session: installing or updating a package that other loaded packages depend on ",
       "(e.g. officer, used by flextable) fails with 'cannot be unloaded'.")
}

library(tidyverse)
library(gtsummary)
library(flextable)
library(officer)
library(MatchIt)   # needed so summary() works on match_out (S3 method for class "matchit")

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
  exists("pain_scores_clean"),
  exists("s4_surgical_keys")   # from 6_cohort_classification.R; used in the surgical-vs-medical check
)

# psm_data, matched_data and match_out must exist in-session (built by
# 10_comparative_PSM.R, run earlier in the same session, NOT loaded
# from a saved .rds).
stopifnot(
  exists("psm_data"),
  exists("matched_data"),
  exists("match_out")
)
if (!inherits(match_out, "matchit")) {
  stop("match_out is not a MatchIt object (class: ", paste(class(match_out), collapse = "/"),
       "). Re-run 10_comparative_PSM.R in THIS session; do not rely on objects restored from a saved workspace.")
}

stopifnot("matched_data must be the version rebuilt by 13_subgroups.R (needs episode_duration_hrs.x)" =
            "episode_duration_hrs.x" %in% names(matched_data))

stopifnot("case_type_refined" %in% names(pca_episodes_flagged))

stopifnot(
  "peak_pain_score" %in% names(pca_episodes_flagged),
  "time_weighted_mean_pain" %in% names(pca_episodes_flagged)
)
cat("Pain score outcomes confirmed pre-computed (from 4_flags.R):\n")
cat("  peak_pain_score, time_weighted_mean_pain\n\n")

# ============================================================
# FLEXTABLE THEME + FOOTNOTE HELPER
# ============================================================
ft_theme <- function(ft) {
  ft <- ft |>
    theme_booktabs() |>
    flextable::font(fontname = "Arial", part = "all") |>   # qualified: ggpubr/survminer also define font()
    fontsize(size = 9, part = "all") |>
    bold(part = "header") |>
    align(align = "left", j = 1, part = "all") |>
    align(align = "center", j = -1, part = "all") |>
    valign(valign = "center", part = "header") |>
    padding(padding.top = 2, padding.bottom = 2,
            padding.left = 4, padding.right = 4, part = "all") |>
    line_spacing(space = 1, part = "all") |>
    set_table_properties(width = 1, layout = "autofit")
  # footnotes smaller than the table body
  tryCatch(fontsize(ft, size = 8, part = "footer"), error = function(e) ft)
}

# Stacked gtsummary tables carry a row-group column whose header reads
# "Group"; blank it so the table reads cleanly.
fix_group_col <- function(ft) {
  gk <- grep("group", ft$col_keys, value = TRUE, ignore.case = TRUE)
  if (length(gk) == 1) {
    ft <- tryCatch(set_header_labels(ft, values = setNames(list(""), gk)),
                   error = function(e) ft)
  }
  ft
}

# Caption "Table N." in bold followed by the title (falls back to plain text
# if this flextable version does not accept formatted captions).
set_cap <- function(ft, num, title) {
  tryCatch(
    set_caption(ft, caption = as_paragraph(as_b(paste0(num, " ")), title)),
    error = function(e) set_caption(ft, caption = paste(num, title))
  )
}

# Try to attach a footnote to a gtsummary table (incl. stacked tables). If
# that fails in this gtsummary version, return the table unchanged and flag
# it so the footnote is written as a text paragraph in the Word export.
attach_footnote <- function(tbl, text) {
  out <- tryCatch(tbl |> modify_footnote(all_stat_cols() ~ text),
                  error = function(e) NULL)
  if (is.null(out)) {
    cat("  (footnote could not be attached to this table object; it will be\n")
    cat("   written as text under the table in the Word export)\n")
    list(tbl = tbl, ok = FALSE)
  } else {
    list(tbl = out, ok = TRUE)
  }
}

# ---- Denominator helpers: footnotes carry computed n, never typed n ----
n_ev   <- function(x) sum(!is.na(x))
fmt_n  <- function(x) formatC(x, big.mark = ",", format = "d")
and_list <- function(x) {
  if (length(x) <= 2) paste(x, collapse = " and ")
  else paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
}
# vars: named vector, name = wording used in the footnote, value = column
denom_sentence <- function(df, vars, total) {
  ev    <- vapply(vars, function(v) n_ev(df[[v]]), integer(1))
  short <- ev[ev < total]
  if (length(short) == 0) return("")
  sprintf("Evaluable episodes differ from %s for %s; percentages and medians use these denominators. ",
          fmt_n(total), and_list(sprintf("%s (n=%s)", names(short), fmt_n(short))))
}

# ============================================================
# SHARED BUILDING BLOCKS: line of therapy and true patient ID
# ============================================================
# drug_line_number counts DISTINCT DRUG CHANGES within an admission
# (PAT_ENC_CSN_ID), so three consecutive morphine episodes are all
# "1st line" and a switch to oxycodone afterwards is "2nd line". Grouped
# by admission to match is_first_line / is_rotation in 4_flags.R, so
# "1st line" is identical to is_first_line == TRUE.
pca_line_lookup <- pca_episodes_flagged |>
  arrange(PAT_ENC_CSN_ID, episode_start) |>
  group_by(PAT_ENC_CSN_ID) |>
  mutate(
    drug_changed     = is.na(lag(drug)) | drug != lag(drug),
    drug_line_number = cumsum(drug_changed)
  ) |>
  ungroup() |>
  select(PAT_ENC_CSN_ID, episode_id, drug_line_number)

admissions_patid <- admissions_clean |>
  select(PAT_ENC_CSN_ID, PAT_ID) |>
  distinct(PAT_ENC_CSN_ID, .keep_all = TRUE)

episodes_per_patient <- pca_episodes_flagged |>
  left_join(admissions_patid, by = "PAT_ENC_CSN_ID") |>
  count(PAT_ID, name = "n_episodes_patient")

line_label <- function(x) {
  factor(
    case_when(x == 1 ~ "1st line",
              x == 2 ~ "2nd line",
              x == 3 ~ "3rd line",
              x >= 4 ~ "4th line or later"),
    levels = c("1st line", "2nd line", "3rd line", "4th line or later")
  )
}

# ============================================================
# TABLE 1 — COHORT DEMOGRAPHICS
# case_type_refined at PATIENT-level.
# ============================================================
cat("Building Table 1 — cohort demographics...\n")

t1_footnote <- paste0(
  "Descriptive summary of all episodes, including rotations and later episodes; not used for between-drug comparison (Table 3 uses the matched cohort). ",
  "Patient-level: one row per patient (earliest episode); episodes per patient counts all admissions. ",
  "First line = first PCA/NCA episode in the admission (no prior opioid that admission); line of therapy counts distinct drug changes within the admission. ",
  "Abbreviations: IQR, interquartile range; PCA, patient-controlled analgesia; NCA, nurse-controlled analgesia.")

patients <- pca_episodes_flagged |>
  left_join(admissions_patid, by = "PAT_ENC_CSN_ID") |>
  arrange(PAT_ID, episode_start) |>
  distinct(PAT_ID, .keep_all = TRUE) |>
  left_join(episodes_per_patient, by = "PAT_ID") |>
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
  left_join(pca_line_lookup, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(
    drug = factor(drug, levels = c("morphine", "oxycodone",
                                   "morphine+ketamine", "fentanyl"),
                  labels = c("Morphine", "Oxycodone",
                             "Morphine + ketamine", "Fentanyl")),
    mode = factor(mode, levels = c("PCA", "NCA", "unknown"),
                  labels = c("PCA", "NCA", "Unknown")),
    case_type_refined = factor(case_type_refined, levels = c("surgical", "medical"),
                               labels = c("Surgical", "Medical")),
    drug_line_label = line_label(drug_line_number)
  )

tbl1_patients <- patients |>
  select(AGE, WEIGHT_KG, GENDER, LOS_days, n_episodes_patient, age_group, case_type_refined) |>
  tbl_summary(
    label = list(
      AGE ~ "Age, years",
      WEIGHT_KG ~ "Weight, kg",
      GENDER ~ "Sex",
      LOS_days ~ "Length of stay, days",
      n_episodes_patient ~ "PCA/NCA episodes per patient (all admissions)",
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
  modify_header(label = "Characteristic", stat_0 = "Median (IQR) or n (%)") |>
  modify_caption("Table 1. Patient and episode characteristics") |>
  bold_labels()

tbl1_episodes <- episodes_summary |>
  select(drug, mode, case_type_refined, drug_line_label, episode_duration_hrs,
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
      drug_line_label ~ "Line of therapy (within admission)",
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
  modify_header(label = "Characteristic") |>
  bold_labels()

tbl1 <- tbl_stack(
  list(tbl1_patients, tbl1_episodes),
  group_header = c(sprintf("Patient-level (n=%d patients)", nrow(patients)),
                   sprintf("Episode-level (n=%d episodes)", nrow(episodes_summary)))
)
r1 <- attach_footnote(tbl1, t1_footnote)
tbl1 <- r1$tbl
t1_fn_ok <- r1$ok
cat("Table 1 built\n\n")

# ============================================================
# SUPPLEMENTARY — MORPHINE NCA vs HOWARD ET AL. 2010
# Applies ONLY Howard's stated criteria: morphine NCA after (major)
# surgery, excluding patients changed to NCA from another technique.
# Console summary only (informal anchor for outcome definitions); not
# exported to Word and not a manuscript table.
# Full comparison, sensitivity pools and age bands: HOWARD_FAITHFUL.R
# ============================================================
cat("============================================================\n")
cat("SUPPLEMENTARY — MORPHINE NCA vs HOWARD ET AL. 2010 (Howard criteria only)\n")
cat("(descriptive context only; definitions differ — see Discussion)\n")
cat("============================================================\n\n")

morphine_nca <- pca_episodes_flagged |>
  filter(drug == "morphine", mode == "NCA", is_first_episode,
         case_type_refined == "surgical")

# Unfiltered antiemetic flags from the raw administrations (closest
# administration analogue of Howard's "nausea or vomiting on one or more
# occasions"); 2+ calendar days is a proxy for the "more than one
# postoperative day" arm of their severe PONV definition.
have_ae <- exists("antiemetic_events")
if (have_ae) {
  ae_flags <- morphine_nca |>
    select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end) |>
    left_join(antiemetic_events, by = "PAT_ENC_CSN_ID", relationship = "many-to-many") |>
    filter(!is.na(ADMIN_TIME), ADMIN_TIME >= episode_start, ADMIN_TIME <= episode_end) |>
    group_by(PAT_ENC_CSN_ID, episode_id) |>
    summarise(n_ae_days = n_distinct(as_date(ADMIN_TIME)), .groups = "drop")
  morphine_nca <- morphine_nca |>
    left_join(ae_flags, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    mutate(any_antiemetic_during = !is.na(n_ae_days),
           antiemetic_multiday = coalesce(n_ae_days, 0L) >= 2L)
} else {
  cat("antiemetic_events not found in session: antiemetic rows will show n/a\n")
  morphine_nca <- morphine_nca |>
    mutate(any_antiemetic_during = NA, antiemetic_multiday = NA)
}

n_adm <- nrow(morphine_nca)
n_pat <- n_distinct(left_join(morphine_nca, admissions_patid, by = "PAT_ENC_CSN_ID")$PAT_ID)
dur_h <- as.numeric(difftime(morphine_nca$episode_end, morphine_nca$episode_start, units = "hours"))
nal   <- morphine_nca$any_naloxone_reversal_final
neo_i <- which(morphine_nca$age_group == "Neonate (<28 days)")
old_i <- setdiff(seq_len(n_adm), neo_i)

fmt_ci <- function(x, n) {
  ci <- binom.test(x, n)$conf.int
  sprintf("%d/%s (%.2f%%; 95%% CI %.2f\u2013%.2f%%)", x, format(n, big.mark = ","),
          100 * x / n, 100 * ci[1], 100 * ci[2])
}
pct1 <- function(x) if (all(is.na(x))) "n/a" else sprintf("%.1f%%", 100 * mean(x, na.rm = TRUE))

howard_a <- data.frame(
  Measure = c("Patients", "Aged under 6 years", "Median age", "Mean NCA duration",
              "Naloxone reversal / serious adverse events",
              "   beyond the neonatal period", "   neonates (under 1 month)",
              "Any antiemetic during NCA (administration-based)",
              "Antiemetic on 2 or more days (proxy)*",
              "Antipruritic (administration-based)"),
  `This study` = c(
    sprintf("%s admissions (%s patients)", format(n_adm, big.mark = ","), format(n_pat, big.mark = ",")),
    sprintf("%.0f%%", 100 * mean(morphine_nca$AGE < 6, na.rm = TRUE)),
    sprintf("%.1f years", median(morphine_nca$AGE, na.rm = TRUE)),
    sprintf("%.1f h", mean(dur_h, na.rm = TRUE)),
    fmt_ci(sum(nal, na.rm = TRUE), sum(!is.na(nal))),
    fmt_ci(sum(nal[old_i], na.rm = TRUE), length(old_i)),
    sprintf("%d/%d", sum(nal[neo_i], na.rm = TRUE), length(neo_i)),
    pct1(morphine_nca$any_antiemetic_during),
    pct1(morphine_nca$antiemetic_multiday),
    pct1(morphine_nca$any_antipruritic_final)),
  `Howard et al. 2010` = c("10,079", "80% aged 6 or under", "2.3 years", "43.7 h",
                           "39 (0.39%)", "26/9,569 (0.27%)", "13/510 (2.5%)",
                           "25% any PONV (symptom-based)",
                           "13.5% clinically severe PONV",
                           "9.4% any itch (symptom-based)"),
  check.names = FALSE, stringsAsFactors = FALSE)

band_levels <- c("<1 month", "1-36 months", "3-6 years", "7-10 years", "11-15 years", "16+ years")
howard_b <- morphine_nca |>
  mutate(howard_age = case_when(
    age_group == "Neonate (<28 days)" ~ "<1 month",
    is.na(AGE) ~ NA_character_,
    AGE < 3  ~ "1-36 months",
    AGE < 7  ~ "3-6 years",
    AGE < 11 ~ "7-10 years",
    AGE < 16 ~ "11-15 years",
    TRUE     ~ "16+ years"),
    howard_age = factor(howard_age, levels = band_levels)) |>
  filter(!is.na(howard_age)) |>
  group_by(howard_age) |>
  summarise(n = n(),
            `Any antiemetic, %` = round(100 * mean(any_antiemetic_during), 1),
            `Antiemetic on 2+ days, %` = round(100 * mean(antiemetic_multiday), 1),
            `Antipruritic, %` = round(100 * mean(any_antipruritic_final, na.rm = TRUE), 1),
            `Naloxone, n` = sum(any_naloxone_reversal_final, na.rm = TRUE),
            .groups = "drop") |>
  rename(`Age band` = howard_age, `Episodes, n` = n) |>
  mutate(`Age band` = as.character(`Age band`))

cat("Howard comparison (console summary only, not exported to Word). Key figures:\n")
print(howard_a, row.names = FALSE)
cat("\nBy Howard age band:\n")
print(howard_b, width = Inf)
cat("\n(Antiemetic and antipruritic are administration-based; Howard recorded symptoms.\n")
cat("Use as an informal anchor for the outcome definitions, not as validation.)\n\n")

# ============================================================
# TABLE 2 — OXYCODONE DESCRIPTIVE (S4, all episodes)
# Stacked: patient-level + episode-level, mirroring Table 1.
# ============================================================
cat("Building Table 2 — oxycodone descriptive (S4)...\n")

s4_oxy_t2 <- pca_episodes_flagged |>
  filter(drug == "oxycodone") |>
  semi_join(episode_metrics |> filter(drug == "oxycodone"),
            by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  left_join(pca_line_lookup, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  left_join(admissions_patid, by = "PAT_ENC_CSN_ID") |>
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
    ),
    drug_line_label = line_label(drug_line_number)
  ) |>
  filter(mode %in% c("PCA", "NCA"))

cat("Line of therapy check — is_first_line TRUE must appear ONLY under '1st line':\n")
print(table(is_first_line = s4_oxy_t2$is_first_line,
            line = s4_oxy_t2$drug_line_label, useNA = "ifany"))
cat("\n")

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

# ---- Which oxycodone episodes are NOT in Table 2, and why (computed) ----
# An oxycodone episode is excluded if (a) it has no row in episode_metrics
# (in practice: no recorded episode end time, so no duration/shift data) or
# (b) its PCA/NCA mode is undetermined.
oxy_all_t2   <- pca_episodes_flagged |> filter(drug == "oxycodone")
oxy_met_t2   <- episode_metrics |> filter(drug == "oxycodone")
k2           <- c("PAT_ENC_CSN_ID", "episode_id")
oxy_in_met   <- semi_join(oxy_all_t2, oxy_met_t2, by = k2)
oxy_no_met   <- anti_join(oxy_all_t2, oxy_met_t2, by = k2)
n_no_metrics <- nrow(oxy_no_met)
n_no_end     <- sum(is.na(oxy_no_met$episode_end))
n_unk_mode   <- sum(!(oxy_in_met$mode %in% c("PCA", "NCA")))

cat(sprintf("Table 2 exclusions: %d oxycodone episodes in total; %d with no episode_metrics row (%d of them with no recorded end time); %d with undetermined mode; %d included\n",
            nrow(oxy_all_t2), n_no_metrics, n_no_end, n_unk_mode, nrow(s4_oxy_t2)))
if (nrow(oxy_all_t2) - n_no_metrics - n_unk_mode != nrow(s4_oxy_t2)) {
  cat("  WARNING: the exclusion counts do not add up to the Table 2 population;\n")
  cat("  check the filters before using the footnote.\n")
}

t2_excl_txt <- sprintf(
  "Of %s oxycodone episodes, %s are included; %d with %s and %d with undetermined PCA/NCA mode are excluded.",
  format(nrow(oxy_all_t2), big.mark = ","), format(nrow(s4_oxy_t2), big.mark = ","),
  n_no_metrics,
  if (n_no_end == n_no_metrics) "no recorded episode end time" else "no episode-level metrics",
  n_unk_mode)

if (n_ev(s4_oxy_t2$peak_pain_score) != n_ev(s4_oxy_t2$time_weighted_mean_pain))
  warning("Peak and time-weighted pain rows have different n: split them in the footnote")
t2_denom_txt <- denom_sentence(
  s4_oxy_t2,
  c("antiemetic"          = "reactive_antiemetic_final",
    "pain scores"         = "peak_pain_score",
    "total dose"          = "total_dose_mcg_kg_hr",
    "background infusion" = "mcg_per_kg_hr"),
  total = nrow(s4_oxy_t2))

t2_footnote <- paste0(
  "Continuous: median (IQR); categorical: n (%). Patient-level: one row per patient (earliest oxycodone episode); episode-level: every oxycodone episode, including rotations. ",
  "Line of therapy counts distinct drug changes within an admission. ",
  "Antiemetic and time-weighted pain definitions as in Methods (episodes with unreliable anaesthesia timing excluded from antiemetic measures); antipruritic, no time filter; naloxone reversal manually chart-reviewed. ",
  t2_denom_txt,
   "Total dose is restricted to episodes with \u226580% pump flowsheet coverage. ",
  t2_excl_txt, " ",
  "Abbreviations: IQR, interquartile range; PCA, patient-controlled analgesia; NCA, nurse-controlled analgesia.")

# ---- Original by-PCA/NCA episode-level table, retained as supplementary ----
tbl2_by_mode <- s4_oxy_t2 |>
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
      reactive_antiemetic_final ~ "Antiemetic administration (reactive definition), n (%)",
      any_antipruritic_final ~ "Any antipruritic — n (%)",
      any_naloxone_reversal_final ~ "Naloxone reversal — n (%)",
      peak_pain_score ~ "Peak pain score",
      time_weighted_mean_pain ~ "Time-weighted mean pain score",
      total_dose_mcg_kg_hr ~ "Total dose (background + bolus), µg/kg/h",
      mcg_per_kg_hr ~ "Background infusion only, µg/kg/h"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}–{p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1),
    missing = "no"
  ) |>
  add_overall(last = FALSE) |>
  modify_header(label = "Outcome",
                stat_0 = sprintf("Overall (n=%d)", nrow(s4_oxy_t2)),
                stat_1 = sprintf("PCA (n=%d)", sum(s4_oxy_t2$mode == "PCA")),
                stat_2 = sprintf("NCA (n=%d)", sum(s4_oxy_t2$mode == "NCA"))) |>
  modify_caption(sprintf("Oxycodone PCA/NCA episode characteristics by mode (all episodes, n=%d) — not exported",
                         nrow(s4_oxy_t2))) |>
  bold_labels() |>
  modify_footnote(all_stat_cols() ~ t2_footnote)

# ---- Main Table 2: patient-level + episode-level ----
t2_patients <- s4_oxy_t2 |>
  arrange(PAT_ID, episode_start) |>
  distinct(PAT_ID, .keep_all = TRUE)

tbl2_patients <- t2_patients |>
  select(AGE, WEIGHT_KG, GENDER, age_group) |>
  tbl_summary(
    label = list(AGE ~ "Age, years", WEIGHT_KG ~ "Weight, kg",
                 GENDER ~ "Sex", age_group ~ "Age group"),
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1), missing = "no"
  ) |>
  modify_header(label = "Characteristic", stat_0 = "Median (IQR) or n (%)") |>
  modify_caption(sprintf("Table 2. Oxycodone PCA/NCA patient and episode characteristics (%d patients / %d episodes)",
                         nrow(t2_patients), nrow(s4_oxy_t2))) |>
  bold_labels()

tbl2_episodes <- s4_oxy_t2 |>
  select(mode, drug_line_label, episode_duration_hrs, proc_category_display,
         reactive_antiemetic_final, any_antipruritic_final, any_naloxone_reversal_final,
         peak_pain_score, time_weighted_mean_pain, total_dose_mcg_kg_hr, mcg_per_kg_hr) |>
  mutate(
    reactive_antiemetic_final = factor(reactive_antiemetic_final, levels = c(FALSE, TRUE), labels = c("No", "Yes")),
    any_antipruritic_final = factor(any_antipruritic_final, levels = c(FALSE, TRUE), labels = c("No", "Yes")),
    any_naloxone_reversal_final = factor(any_naloxone_reversal_final, levels = c(FALSE, TRUE), labels = c("No", "Yes"))
  ) |>
  tbl_summary(
    label = list(
      mode ~ "Mode",
      drug_line_label ~ "Line of therapy (oxycodone as 1st/2nd/3rd+ opioid in the admission)",
      episode_duration_hrs ~ "Episode duration, hours",
      proc_category_display ~ "Procedure category (surgical only; non-surgical episodes shown as one group)",
      reactive_antiemetic_final ~ "Antiemetic administration (reactive definition), n (%)",
      any_antipruritic_final ~ "Any antipruritic — n (%)",
      any_naloxone_reversal_final ~ "Naloxone reversal — n (%)",
      peak_pain_score ~ "Peak pain score",
      time_weighted_mean_pain ~ "Time-weighted mean pain score",
      total_dose_mcg_kg_hr ~ "Total dose (background + bolus), µg/kg/h",
      mcg_per_kg_hr ~ "Background infusion only, µg/kg/h"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1), missing = "no"
  ) |>
  modify_header(label = "Characteristic") |>
  bold_labels()

tbl2 <- tbl_stack(
  list(tbl2_patients, tbl2_episodes),
  group_header = c(sprintf("Patient-level (n=%d patients)", nrow(t2_patients)),
                   sprintf("Episode-level (n=%d episodes)", nrow(s4_oxy_t2)))
)
r2 <- attach_footnote(tbl2, t2_footnote)
tbl2 <- r2$tbl
t2_fn_ok <- r2$ok
cat("Table 2 built (stacked patient/episode; by-mode version kept as tbl2_by_mode)\n\n")

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

# keep_specialties: the specialties kept as their own level. When NULL the
# n>=30 rule is applied within df itself; for post_match it is passed the
# PRE-match list so the matched table uses the same collapse as the matching.
add_psm_covariates <- function(df, keep_specialties = NULL) {
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
  if (is.null(keep_specialties)) {
    keep_specialties <- df |>
      count(specialty_merged) |>
      filter(n >= SPECIALTY_VOLUME_THRESHOLD) |>
      pull(specialty_merged)
  }
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

# Specialties kept as their own level in the matching model (pre-match list)
keep_from_pre <- setdiff(unique(as.character(pre_match$specialty_collapsed)), "Other specialty")
cat("Specialties kept as their own level (from the pre-match pool):\n")
print(keep_from_pre)
cat("\n")

post_match <- add_psm_covariates(post_match, keep_specialties = keep_from_pre)

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

# Year shown as a plain integer (no thousands separator or decimals)
tbl3_digits <- list(all_continuous() ~ 1,
                    year ~ function(x) sprintf("%.0f", x))

tbl3_pre <- pre_match |>
  select(drug, all_of(balance_vars)) |>
  mutate(across(where(is.logical), ~ factor(., labels = c("No", "Yes")))) |>
  tbl_summary(
    by = drug, label = tbl3_labels,
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = tbl3_digits, missing = "no"
  ) |>
  modify_header(label = "Variable",
                stat_1 = sprintf("Morphine (n=%d)", sum(pre_match$drug == "Morphine")),
                stat_2 = sprintf("Oxycodone (n=%d)", sum(pre_match$drug == "Oxycodone"))) |>
  bold_labels()

tbl3_post <- post_match |>
  select(drug, all_of(balance_vars)) |>
  mutate(across(where(is.logical), ~ factor(., labels = c("No", "Yes")))) |>
  tbl_summary(
    by = drug, label = tbl3_labels,
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = tbl3_digits, missing = "no"
  ) |>
  modify_header(label = "Variable",
                stat_1 = sprintf("Morphine (n=%d)", sum(post_match$drug == "Morphine")),
                stat_2 = sprintf("Oxycodone (n=%d)", sum(post_match$drug == "Oxycodone"))) |>
  bold_labels()

# Balance statement computed from the matching model (match_out, script 10)
smd_post <- summary(match_out, standardize = TRUE)$sum.matched[, "Std. Mean Diff."]
smd_post <- abs(smd_post[!is.na(smd_post) & names(smd_post) != "distance"])
smd_over <- sort(smd_post[smd_post > 0.10], decreasing = TRUE)

clean_term <- function(x) {
  x <- sub("^proc_category", "procedure: ", x)
  x <- sub("^specialty_collapsed", "specialty: ", x)
  x <- sub("^had_intraop_", "intraop ", x)
  x <- sub("^had_", "", x)
  sub("TRUE$", "", x)
}

smd_sentence <- if (length(smd_over) == 0) {
  sprintf("All post-match standardised mean differences <0.10 (maximum %.3f).", max(smd_post))
} else {
  sprintf("Maximum post-match standardised mean difference %.3f; %d of %d covariate terms exceeded 0.10 (%s).",
          max(smd_post), length(smd_over), length(smd_post),
          paste(clean_term(names(smd_over)), collapse = ", "))
}
cat("Balance statement for the Table S1 footnote:\n  ", smd_sentence, "\n\n")

tbl3 <- tbl_merge(
  list(tbl3_pre, tbl3_post),
  tab_spanner = c("Pre-match", "Post-match")
) |>
  modify_caption("Supplementary Table S1. Propensity score matching — covariate balance") |>
  modify_footnote(
    all_stat_cols() ~ paste0(
      "Continuous: median (IQR); categorical: n (%). ", smd_sentence, " ",
      "Concurrent epidural analgesia was an eligibility criterion (excluded), so is not a balance covariate. ",
      "Abbreviations: ASA, American Society of Anesthesiologists physical status; IQR, interquartile range; NCA, nurse-controlled analgesia; NSAID, nonsteroidal anti-inflammatory drug; PCA, patient-controlled analgesia.")
  )
cat("Table 3 built\n\n")

# ============================================================
# TABLE 4 — PRIMARY OUTCOMES
# (Composite side effect removed. Reports antiemetic, antipruritic and
# naloxone reversal independently.)
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

# Pain-score rows are shown twice: all ages, and the pre-specified age >= 6
# subset (the matching itself used all ages, so the >=6 rows are a subset of
# the matched cohort and some pairs lose a member).
table4_matched_outcomes <- table4_matched_outcomes |>
  mutate(
    peak_pain_score_6plus = if_else(!is.na(AGE) & AGE >= 6, peak_pain_score, NA_real_),
    time_weighted_mean_pain_6plus = if_else(!is.na(AGE) & AGE >= 6, time_weighted_mean_pain, NA_real_)
  )
n_peak_all <- sum(!is.na(table4_matched_outcomes$peak_pain_score))
n_tw_all   <- sum(!is.na(table4_matched_outcomes$time_weighted_mean_pain))
n_peak_6   <- sum(!is.na(table4_matched_outcomes$peak_pain_score_6plus))
n_tw_6     <- sum(!is.na(table4_matched_outcomes$time_weighted_mean_pain_6plus))
cat(sprintf("Pain rows: all ages n=%d (peak) / %d (time-weighted); age >=6 n=%d / %d\n\n",
            n_peak_all, n_tw_all, n_peak_6, n_tw_6))

n_matched <- nrow(table4_matched_outcomes)
by_arm <- function(var) {
  tab <- table4_matched_outcomes |> group_by(drug) |>
    summarise(n = n_ev(.data[[var]]), .groups = "drop")
  setNames(tab$n, as.character(tab$drug))
}
fmt_arm <- function(x) sprintf("%d; morphine %d, oxycodone %d", sum(x), x[["Morphine"]], x[["Oxycodone"]])

pain_arm <- by_arm("peak_pain_score")
t4_denom <- paste0("Evaluable episodes differ from ", n_matched,
                   " for pain scores (", fmt_arm(pain_arm), ")")
have_dose <- all(c("total_dose_mcg_kg_hr", "background_mcg_kg_hr") %in% names(table4_matched_outcomes))
if (have_dose) {
  dose_arm <- by_arm("total_dose_mcg_kg_hr"); bg_arm <- by_arm("background_mcg_kg_hr")
  t4_denom <- if (identical(dose_arm, bg_arm)) {
    paste0(t4_denom, " and for total dose and background infusion (", fmt_arm(dose_arm), ").")
  } else {
    paste0(t4_denom, ", total dose (", fmt_arm(dose_arm), "), and background infusion (", fmt_arm(bg_arm), ").")
  }
} else t4_denom <- paste0(t4_denom, ".")

# Footnote piece that depends on the data: whether flowsheet linkage
# differed by drug in THIS matched population.
n_t4 <- nrow(table4_matched_outcomes)

dose_sentence <- "Total dose is restricted to episodes with \u226580% pump flowsheet coverage."
if (!is.null(episode_true_total_dose)) {
  link_df <- table4_base |>
    select(PAT_ENC_CSN_ID, episode_id, drug) |>
    left_join(
      episode_true_total_dose |> distinct(PAT_ENC_CSN_ID, episode_id) |> mutate(linked = TRUE),
      by = c("PAT_ENC_CSN_ID", "episode_id")
    ) |>
    mutate(linked = replace_na(linked, FALSE))
  link_p <- tryCatch(
    suppressWarnings(chisq.test(table(link_df$drug, link_df$linked))$p.value),
    error = function(e) NA_real_)
  cat(sprintf("Flowsheet linkage by drug in the matched population: chi-squared p = %.3f\n\n", link_p))
  if (!is.na(link_p) && link_p < 0.05) {
    dose_sentence <- paste0(dose_sentence, sprintf(
      " Linkage differed by drug (p=%.3f); a separate inverse-probability-weighting sensitivity check confirmed the unweighted comparison is robust to this.", link_p))
  } else if (!is.na(link_p)) {
    dose_sentence <- paste0(dose_sentence, sprintf(
      " Linkage did not differ by drug (p=%.2f).", link_p))
  }
}

t4_footnote <- paste0(
  "Propensity-matched cohort (Methods). Continuous: median (IQR); categorical: n (%). ",
  "Tests: chi-squared (categorical), Wilcoxon rank-sum (continuous), Fisher's exact (naloxone reversal). ",
  "Pain scores are shown for all ages and for the pre-specified subset aged \u22656 years (some pairs lose a member). ",
  "Time-weighted mean pain weights each reading by the time until the next reading. ",
  "Episode duration is descriptive, not a pre-defined outcome. ",
  dose_sentence, " ", t4_denom, " ",
  "Abbreviations: IQR, interquartile range; NCA, nurse-controlled analgesia; PCA, patient-controlled analgesia.")

tbl4 <- table4_matched_outcomes |>
  select(drug,
         reactive_antiemetic_final,
         any_antipruritic_final, any_naloxone_reversal_final,
         peak_pain_score, time_weighted_mean_pain,
         peak_pain_score_6plus, time_weighted_mean_pain_6plus,
         episode_duration_hrs,
         any_of(c("total_dose_mcg_kg_hr", "background_mcg_kg_hr"))) |>
  tbl_summary(
    by = drug,
    label = list(
      reactive_antiemetic_final ~ "Antiemetic administration (reactive definition)",
      any_antipruritic_final ~ "Antipruritic",
      any_naloxone_reversal_final ~ "Naloxone reversal (respiratory)",
      peak_pain_score ~ sprintf("Peak pain score, all ages (n=%d), median (IQR)", n_peak_all),
      time_weighted_mean_pain ~ sprintf("Time-weighted mean pain score, all ages (n=%d), median (IQR)", n_tw_all),
      peak_pain_score_6plus ~ sprintf("Peak pain score, age \u22656 years (n=%d), median (IQR)", n_peak_6),
      time_weighted_mean_pain_6plus ~ sprintf("Time-weighted mean pain score, age \u22656 years (n=%d), median (IQR)", n_tw_6),
      episode_duration_hrs ~ "Episode duration, hours",
      total_dose_mcg_kg_hr ~ "Total dose (background + bolus), µg/kg/h",
      background_mcg_kg_hr ~ "Background infusion only (same episodes), µg/kg/h"
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
    label = "Outcome",
    stat_0 = sprintf("Overall (n=%d)", nrow(table4_matched_outcomes)),
    stat_1 = sprintf("Morphine (n=%d)", sum(table4_matched_outcomes$drug == "Morphine")),
    stat_2 = sprintf("Oxycodone (n=%d)", sum(table4_matched_outcomes$drug == "Oxycodone")),
    p.value = "p-value"
  ) |>
  modify_caption("Table 3. Primary outcomes — propensity-matched comparison") |>
  bold_labels() |>
  bold_p(t = 0.05) |>
  modify_footnote(all_stat_cols() ~ t4_footnote)

# Remove gtsummary's automatic footnote on the p-value column (the test
# names): the tests are already named in the main footnote above. If this
# gtsummary version rejects the call, the table is left as it was.
tbl4 <- tryCatch(tbl4 |> modify_footnote(p.value ~ NA_character_),
                 error = function(e) tbl4)

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
# FOUR objects: tbl5 (5a, descriptive rate), ft5b (destinations), ft5c
# (multi-rotation pattern) and ft5d (covariate-adjusted rotation-away OR,
# computed on the pre-match eligible pool due to too few events in the
# matched sample alone; see script 15 — its caption quotes the matched
# pair and event counts, computed live there).
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

# Data range, as a check on the study period in the title line below
cat("Data range for the title line: first episode_start =",
    format(min(pca_episodes_flagged$episode_start, na.rm = TRUE), "%d %b %Y"),
    "| last episode_end =",
    format(max(pca_episodes_flagged$episode_end, na.rm = TRUE), "%d %b %Y"), "\n")

pub <- function(tbl) as_flex_table(tbl) |> fix_group_col() |> ft_theme()
ft1  <- pub(tbl1)
ft2  <- pub(tbl2)
ft3  <- pub(tbl3)
ft4  <- pub(tbl4)
ft5a <- pub(tbl5)
ft5b <- ft_theme(ft5b)
ft5c <- ft_theme(ft5c)
ft5d <- ft_theme(ft5d)

fp_title <- fp_text(font.family = "Arial", font.size = 13, bold = TRUE)
fp_sub   <- fp_text(font.family = "Arial", font.size = 9, italic = TRUE, color = "555555")
fp_head  <- fp_text(font.family = "Arial", font.size = 11, bold = TRUE)
fp_note  <- fp_text(font.family = "Arial", font.size = 8)

# Safety net: remove literal markdown asterisks ("**bold**") that survive into the
# Word file (e.g. from captions/headers built in script 15). Works on a copy and
# replaces the output only if every step succeeds; otherwise the file is untouched.
strip_markdown_asterisks <- function(path) {
  tryCatch({
    if (!requireNamespace("zip", quietly = TRUE)) stop("package 'zip' not available")
    tmp <- file.path(tempdir(), paste0("docx_fix_", format(Sys.time(), "%H%M%S")))
    dir.create(tmp)
    utils::unzip(path, exdir = tmp)
    f <- file.path(tmp, "word", "document.xml")
    x <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    n <- length(gregexpr("**", x, fixed = TRUE)[[1]])
    if (identical(as.integer(gregexpr("**", x, fixed = TRUE)[[1]][1]), -1L)) return(0L)
    x <- gsub("**", "", x, fixed = TRUE)
    writeLines(x, f, useBytes = TRUE)
    rel <- list.files(tmp, recursive = TRUE, all.files = TRUE, include.dirs = FALSE)
    newzip <- tempfile(fileext = ".docx")
    zip::zip(zipfile = newzip, files = rel, root = tmp, mode = "mirror", include_directories = FALSE)
    file.copy(newzip, path, overwrite = TRUE)
    n
  }, error = function(e) {
    message("Asterisk clean-up skipped: ", conditionMessage(e))
    NA_integer_
  })
}

add_note <- function(doc, text) body_add_fpar(doc, fpar(ftext(text, fp_note)))
add_head <- function(doc, text) body_add_fpar(doc, fpar(ftext(text, fp_head)))

doc <- read_docx()
# A4 portrait, narrow margins (skipped silently if this officer version lacks it)
doc <- tryCatch(
  body_set_default_section(doc, prop_section(
    page_size = page_size(orient = "portrait", width = 8.27, height = 11.69),
    page_margins = page_mar(top = 0.8, bottom = 0.8, left = 0.8, right = 0.8))),
  error = function(e) doc)

# Study period in the title below: October 2014 to July 2026 (end date
# confirmed). The data range printed above is a check against it.
doc <- doc |>
  body_add_fpar(fpar(ftext("Oxycodone PCA/NCA in Children: Manuscript Tables", fp_title))) |>
  body_add_fpar(fpar(ftext("Addenbrooke's Hospital \u00b7 October 2014 \u2013 July 2026", fp_sub))) |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft1)
if (!t1_fn_ok) doc <- doc |> add_note(t1_footnote)
doc <- doc |>
  body_add_break() |>
  body_add_flextable(ft2)
if (!t2_fn_ok) doc <- doc |> add_note(t2_footnote)
doc <- doc |>
  body_add_break() |>
  body_add_flextable(ft3) |>
  body_add_break() |>
  body_add_flextable(ft4) |>
  body_add_break() |>
  body_add_flextable(ft5a) |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft5b) |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft5c) |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft5d) |>
  body_add_break() |>
  add_head("Notes and definitions") |>
  add_note("Matched cohort = primary comparative population: each child's first-ever recorded PCA/NCA episode at this centre (therefore first line and first episode), one episode per patient, surgical, no concurrent epidural, medical oncology excluded. Descriptive population = all oxycodone episodes, including rotations.") |>
  add_note("First line (as used in Tables 1 and 2) = first PCA/NCA episode in the admission (no prior opioid in that admission). The matched cohort is stricter: it requires no earlier recorded PCA/NCA episode in any admission at this centre. PCA/NCA at other hospitals cannot be excluded, so this does not establish first opioid exposure ever.") |>
  add_note("Time-weighted mean pain score: each reading is weighted by the time until the next reading (last observation carried forward); time before the first reading is excluded.") |>
  add_note("Side-effect measures are administration-based proxies (antiemetic for PONV, antipruritic for itch); naloxone reversal is manually chart-reviewed.")

print(doc, target = output_path)
n_ast <- strip_markdown_asterisks(output_path)
if (!is.na(n_ast) && n_ast > 0) cat("Removed", n_ast, "literal ** markers from the Word file\n")

cat("\n============================================================\n")
cat("ANALYSIS 11 COMPLETE\n")
cat("Output:", output_path, "(Tables 1-5 and notes)\n")
cat("============================================================\n")