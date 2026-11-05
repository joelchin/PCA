# ============================================================
# 18 — FINAL OUTCOMES — final comparison in the matched cohort
# Oxycodone vs morphine PCA/NCA
#
# Run AFTER 09_comparative_psm.R in the same R session (needs matched_data,
# pca_episodes_flagged, pain_scores_clean). Aggregate output only; saves
# outcomes_final_summary_<primary|sensitivity>.csv.
#
# What it does differently from the quick printouts in script 10:
#   * respects the matched pairs: cluster-robust SEs by pair (subclass)
#   * gives EFFECT SIZES WITH 95% CIs, not just p-values
#   * adjusts for residual imbalance after matching (ADJ below) and also
#     shows the UNADJUSTED estimate next to it
#   * side effects are also expressed per episode-hour (oxycodone episodes
#     are ~12 h longer, so "any event during the episode" favours more events)
#   * pain in the whole episode and in the FIRST 24 / 48 h (age >= 6)
#   * consumption as a ratio of means (oxycodone / morphine)
#   * loading dose described, and key outcomes re-estimated within
#     loading-dose yes / no strata
#
# NOT included (needs data this script does not have): fixed-window CONSUMPTION
# and event-time side-effect windows. The dose file is episode-level, so these
# would need the time-stamped flowsheet / administration records re-extracted.
#
# Reading the effect columns (oxycodone vs morphine):
#   pp     = risk difference in percentage points (+ = more with oxycodone)
#   ratio  = ratio of means or rate ratio (1 = no difference)
#   diff   = difference in means (+ = higher with oxycodone)
# ============================================================
library(tidyverse)
library(lubridate)
if (!requireNamespace("sandwich", quietly = TRUE)) install.packages("sandwich")
stopifnot(requireNamespace("sandwich", quietly = TRUE))

stopifnot(exists("matched_data"), exists("pca_episodes_flagged"), exists("pain_scores_clean"))
stopifnot("subclass" %in% names(matched_data))

KEYS    <- c("PAT_ENC_CSN_ID", "episode_id")
AGE_MIN <- 6
Z       <- qnorm(0.975)

# Covariates that stayed above SMD 0.10 after matching (strict cohort, caliper 0.10):
# year, Respiratory/ENT procedures, medical oncology specialty. Edit if the
# matching changes (check the "Covariates above 0.10" lines printed by script 10).
ADJ <- c("year", "resp_ent", "med_onc")

tag <- if (exists("STRICT_FIRST_EVER") && !STRICT_FIRST_EVER) "sensitivity" else "primary"
cat("============================================================\n")
cat("OUTCOMES_FINAL —", toupper(tag), "cohort\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. ANALYSIS FRAME (one row per matched episode)
# ------------------------------------------------------------
outcome_cols <- c("episode_start", "episode_end", "episode_duration_hrs",
                  "peak_pain_score", "time_weighted_mean_pain",
                  "reactive_antiemetic_final", "any_antipruritic_final",
                  "any_naloxone_reversal_final")

A <- matched_data |>
  select(all_of(KEYS), drug, treated, subclass,
         any_of(c("AGE", "WEIGHT_KG", "year", "proc_category", "specialty_collapsed")),
         any_of(outcome_cols))
need <- setdiff(outcome_cols, names(A))
if (length(need) > 0) {
  A <- A |> left_join(pca_episodes_flagged |> select(all_of(KEYS), any_of(need)), by = KEYS)
}

dose_tot <- tryCatch(
  readRDS("1 - data/3 - processed_data/episode_true_total_dose.rds") |>
    select(all_of(KEYS), total_dose_mcg_kg_hr, background_mcg_kg_hr) |>
    distinct(across(all_of(KEYS)), .keep_all = TRUE),
  error = function(e) { cat("episode_true_total_dose.rds not found: consumption rows skipped\n"); NULL })
dose_ld <- tryCatch(
  readRDS("1 - data/3 - processed_data/episode_dose.rds") |>
    select(all_of(KEYS), had_loading_dose) |>
    distinct(across(all_of(KEYS)), .keep_all = TRUE),
  error = function(e) { cat("episode_dose.rds not found: loading-dose rows skipped\n"); NULL })

if (!is.null(dose_tot)) A <- A |> left_join(dose_tot, by = KEYS) else {
  A$total_dose_mcg_kg_hr <- NA_real_; A$background_mcg_kg_hr <- NA_real_ }
if (!is.null(dose_ld)) A <- A |> left_join(dose_ld, by = KEYS) else A$had_loading_dose <- NA

A <- A |>
  mutate(
    arm = factor(if_else(treated == 1, "Oxycodone", "Morphine"), levels = c("Morphine", "Oxycodone")),
    ae  = as.integer(reactive_antiemetic_final),
    pr  = as.integer(any_antipruritic_final),
    nal = as.integer(any_naloxone_reversal_final),
    ld  = as.integer(coalesce(had_loading_dose, FALSE)),
    had_loading_dose = coalesce(had_loading_dose, FALSE),
    resp_ent = as.integer(proc_category == "Respiratory/ENT"),
    med_onc  = as.integer(specialty_collapsed == "PAEDIATRIC MEDICAL ONCOLOGY")
  )

stopifnot(!anyDuplicated(A[KEYS]))
n_pairs <- A |> count(subclass) |> filter(n == 2) |> nrow()
cat(sprintf("Matched episodes: %d (morphine %d, oxycodone %d) in %d complete pairs\n",
            nrow(A), sum(A$arm == "Morphine"), sum(A$arm == "Oxycodone"), n_pairs))
cat("Adjusting for:", paste(ADJ, collapse = ", "), "\n\n")

# ------------------------------------------------------------
# 2. HELPERS
# ------------------------------------------------------------
fmt_med <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("n/a")
  sprintf("%.1f (%.1f–%.1f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
}
fmt_np <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("n/a")
  sprintf("%d/%d (%.1f%%)", sum(x == 1), length(x), 100 * mean(x == 1))
}

# Effect of oxycodone vs morphine with cluster-robust (by matched pair) 95% CI.
#   type "rd"    : linear-probability model -> risk difference
#   type "rdiff" : linear model             -> difference in means
#   type "rr"    : quasi-Poisson, log link  -> rate/mean ratio
#   rate = TRUE  : adds offset(log(episode duration)) -> events per episode-hour
est_effect <- function(d, y, type, adj = ADJ, rate = FALSE) {
  tryCatch({
    d <- d |> drop_na(all_of(c(y, "arm", "subclass", adj, if (rate) "episode_duration_hrs")))
    if (rate) d <- d |> filter(episode_duration_hrs > 0)
    if (n_distinct(d$arm) < 2 || nrow(d) < 20) return(NULL)
    adj <- adj[vapply(adj, function(v) n_distinct(d[[v]]) > 1, logical(1))]
    rhs <- paste(c("arm", adj), collapse = " + ")
    if (rate) rhs <- paste0(rhs, " + offset(log(episode_duration_hrs))")
    f   <- as.formula(paste(y, "~", rhs))
    fit <- if (type == "rr") glm(f, data = d, family = quasipoisson(link = "log")) else lm(f, data = d)
    vc  <- sandwich::vcovCL(fit, cluster = d$subclass, type = "HC1")
    b   <- unname(coef(fit)["armOxycodone"])
    se  <- sqrt(vc["armOxycodone", "armOxycodone"])
    ex  <- type == "rr"
    list(est = if (ex) exp(b) else b,
         lo  = if (ex) exp(b - Z * se) else b - Z * se,
         hi  = if (ex) exp(b + Z * se) else b + Z * se,
         p   = 2 * pnorm(-abs(b / se)), n = nrow(d))
  }, error = function(e) NULL)
}

fmt_eff <- function(e, kind) {
  if (is.null(e)) return("n/a")
  switch(kind,
         pp    = sprintf("%+.1f pp (%.1f to %.1f)", 100 * e$est, 100 * e$lo, 100 * e$hi),
         ratio = sprintf("%.2f (%.2f to %.2f)", e$est, e$lo, e$hi),
         diff  = sprintf("%+.2f (%.2f to %.2f)", e$est, e$lo, e$hi))
}

# Matched-pair odds ratio from discordant pairs (McNemar-type, exact CI)
pair_or <- function(d, y) {
  tryCatch({
    w <- d |> select(subclass, arm, all_of(y)) |> drop_na() |>
      pivot_wider(names_from = arm, values_from = all_of(y)) |>
      drop_na(Morphine, Oxycodone)
    b  <- sum(w$Oxycodone == 1 & w$Morphine == 0)
    cc <- sum(w$Oxycodone == 0 & w$Morphine == 1)
    if (b + cc == 0) return("pair OR n/a (no discordant pairs)")
    bt <- binom.test(b, b + cc)
    sprintf("pair OR %.2f (%.2f–%.2f; discordant pairs %d vs %d)",
            ifelse(cc == 0, Inf, b / cc),
            bt$conf.int[1] / (1 - bt$conf.int[1]), bt$conf.int[2] / (1 - bt$conf.int[2]), b, cc)
  }, error = function(e) "pair OR n/a")
}

paired_wilcox_p <- function(d, y) {
  tryCatch({
    w <- d |> select(subclass, arm, all_of(y)) |> drop_na() |>
      pivot_wider(names_from = arm, values_from = all_of(y)) |> drop_na(Morphine, Oxycodone)
    if (nrow(w) < 10) return("paired Wilcoxon n/a")
    sprintf("paired Wilcoxon p=%.3f (%d pairs)",
            suppressWarnings(wilcox.test(w$Oxycodone, w$Morphine, paired = TRUE, exact = FALSE)$p.value), nrow(w))
  }, error = function(e) "paired Wilcoxon n/a")
}

RES <- list()
add_res <- function(domain, outcome, analysis, n, m_txt, o_txt, kind = NA, e_adj = NULL, e_unadj = NULL, note = "") {
  RES[[length(RES) + 1]] <<- tibble(
    domain = domain, outcome = outcome, analysis = analysis, n = n,
    morphine = m_txt, oxycodone = o_txt,
    effect_adjusted   = if (is.na(kind)) "" else fmt_eff(e_adj, kind),
    effect_unadjusted = if (is.na(kind)) "" else fmt_eff(e_unadj, kind),
    p_adjusted = if (is.null(e_adj)) NA_real_ else round(e_adj$p, 3),
    note = note)
}
by_arm <- function(d, y, f) c(f(d[[y]][d$arm == "Morphine"]), f(d[[y]][d$arm == "Oxycodone"]))

# ------------------------------------------------------------
# 3. SIDE EFFECTS
# ------------------------------------------------------------
for (s in list(c("ae", "Antiemetic (reactive)"), c("pr", "Antipruritic"))) {
  y <- s[1]; lab <- s[2]; v <- by_arm(A, y, fmt_np)
  add_res("Side effects", lab, "episode: risk difference", sum(!is.na(A[[y]])), v[1], v[2], "pp",
          est_effect(A, y, "rd"), est_effect(A, y, "rd", adj = character(0)), pair_or(A, y))
  add_res("Side effects", lab, "per episode-hour: rate ratio", sum(!is.na(A[[y]])), "", "", "ratio",
          est_effect(A, y, "rr", rate = TRUE), est_effect(A, y, "rr", adj = character(0), rate = TRUE),
          "events per hour of episode time; >1 = more with oxycodone")
}
add_res("Side effects", "Naloxone reversal", "counts only (too few events to model)", sum(!is.na(A$nal)),
        sprintf("%d/%d", sum(A$nal[A$arm == "Morphine"], na.rm = TRUE), sum(A$arm == "Morphine")),
        sprintf("%d/%d", sum(A$nal[A$arm == "Oxycodone"], na.rm = TRUE), sum(A$arm == "Oxycodone")))

# ------------------------------------------------------------
# 4. PAIN — whole episode (all ages, and age >= 6)
# ------------------------------------------------------------
for (age_lab in c("all ages", paste0("age >= ", AGE_MIN))) {
  d <- if (age_lab == "all ages") A else A |> filter(AGE >= AGE_MIN)
  for (s in list(c("time_weighted_mean_pain", "Time-weighted mean pain"), c("peak_pain_score", "Peak pain"))) {
    v <- by_arm(d, s[1], fmt_med)
    add_res("Pain", s[2], paste0("whole episode, ", age_lab, ": mean difference"), sum(!is.na(d[[s[1]]])), v[1], v[2], "diff",
            est_effect(d, s[1], "rdiff"), est_effect(d, s[1], "rdiff", adj = character(0)), paired_wilcox_p(d, s[1]))
  }
}

# ------------------------------------------------------------
# 5. PAIN — fixed windows from episode start (age >= 6)
# Same convention as 04_derive_episode_flags.R / TIME_IN_COMFORT.R: each reading applies until
# the next; time before the first reading is excluded; highest score at a
# timestamp. Windows make the comparison independent of episode length.
# ------------------------------------------------------------
win_stats <- function(frame, W) {
  frame |>
    select(all_of(KEYS), episode_start, episode_end) |>
    mutate(win_end = if_else(episode_end < episode_start + hours(W), episode_end, episode_start + hours(W))) |>
    inner_join(pain_scores_clean |> select(PAT_ENC_CSN_ID, score_time, score_value),
               by = "PAT_ENC_CSN_ID", relationship = "many-to-many") |>
    filter(!is.na(score_time), !is.na(score_value), score_time >= episode_start, score_time <= win_end) |>
    group_by(PAT_ENC_CSN_ID, episode_id, win_end, score_time) |>
    summarise(score_value = max(score_value), .groups = "drop") |>
    arrange(PAT_ENC_CSN_ID, episode_id, score_time) |>
    group_by(PAT_ENC_CSN_ID, episode_id) |>
    mutate(next_time = coalesce(lead(score_time), win_end),
           seg_hrs = as.numeric(difftime(next_time, score_time, units = "hours"))) |>
    summarise(tw_win = ifelse(sum(seg_hrs) > 0, sum(score_value * seg_hrs) / sum(seg_hrs), NA_real_),
              peak_win = max(score_value), .groups = "drop")
}
A6 <- A |> filter(AGE >= AGE_MIN)
for (W in c(24, 48)) {
  d <- A6 |> left_join(win_stats(A6, W), by = KEYS)
  for (s in list(c("tw_win", "Time-weighted mean pain"), c("peak_win", "Peak pain"))) {
    v <- by_arm(d, s[1], fmt_med)
    add_res("Pain", s[2], sprintf("first %d h, age >= %d: mean difference", W, AGE_MIN), sum(!is.na(d[[s[1]]])), v[1], v[2], "diff",
            est_effect(d, s[1], "rdiff"), est_effect(d, s[1], "rdiff", adj = character(0)), paired_wilcox_p(d, s[1]))
  }
}

# ------------------------------------------------------------
# 6. CONSUMPTION (episodes linked to the pump flowsheet, >= 80% coverage)
# ------------------------------------------------------------
if (!is.null(dose_tot)) {
  for (s in list(c("total_dose_mcg_kg_hr", "Total dose (background + boluses), ug/kg/h"),
                 c("background_mcg_kg_hr", "Background infusion only, ug/kg/h"))) {
    v <- by_arm(A, s[1], fmt_med)
    add_res("Consumption", s[2], "ratio of means (oxycodone / morphine)", sum(!is.na(A[[s[1]]])), v[1], v[2], "ratio",
            est_effect(A, s[1], "rr"), est_effect(A, s[1], "rr", adj = character(0)),
            "linked episodes only; loading dose is NOT included")
  }
}
v <- by_arm(A, "episode_duration_hrs", fmt_med)
add_res("Exposure", "Episode duration, h", "ratio of means (oxycodone / morphine)", sum(!is.na(A$episode_duration_hrs)), v[1], v[2], "ratio",
        est_effect(A, "episode_duration_hrs", "rr"), est_effect(A, "episode_duration_hrs", "rr", adj = character(0)),
        paired_wilcox_p(A, "episode_duration_hrs"))

# ------------------------------------------------------------
# 7. LOADING DOSE — described, then key outcomes within strata
# (loading is decided at the START of the episode, so it is a co-intervention,
#  not a matching covariate; strata are a diagnostic, not the primary estimate)
# ------------------------------------------------------------
if (!is.null(dose_ld)) {
  v <- by_arm(A, "ld", fmt_np)
  add_res("Loading dose", "Loading dose prescribed", "proportion", sum(!is.na(A$ld)), v[1], v[2], note = pair_or(A, "ld"))
  for (lv in c(FALSE, TRUE)) {
    d <- A |> filter(had_loading_dose == lv)
    lab <- paste0("loading dose = ", lv)
    if (n_distinct(d$arm) < 2) next
    v <- by_arm(d, "time_weighted_mean_pain", fmt_med)
    add_res("Loading strata", "Time-weighted mean pain", paste0(lab, ": mean difference"), sum(!is.na(d$time_weighted_mean_pain)), v[1], v[2], "diff",
            est_effect(d, "time_weighted_mean_pain", "rdiff"), est_effect(d, "time_weighted_mean_pain", "rdiff", adj = character(0)))
    if (!is.null(dose_tot)) {
      v <- by_arm(d, "total_dose_mcg_kg_hr", fmt_med)
      add_res("Loading strata", "Total dose, ug/kg/h", paste0(lab, ": ratio of means"), sum(!is.na(d$total_dose_mcg_kg_hr)), v[1], v[2], "ratio",
              est_effect(d, "total_dose_mcg_kg_hr", "rr"), est_effect(d, "total_dose_mcg_kg_hr", "rr", adj = character(0)))
    }
    v <- by_arm(d, "ae", fmt_np)
    add_res("Loading strata", "Antiemetic (reactive)", paste0(lab, ": risk difference"), sum(!is.na(d$ae)), v[1], v[2], "pp",
            est_effect(d, "ae", "rd"), est_effect(d, "ae", "rd", adj = character(0)))
  }
}

# ------------------------------------------------------------
# 8. OUTPUT
# ------------------------------------------------------------
res <- bind_rows(RES)
options(width = 220)
print(res |> select(domain, outcome, analysis, n, morphine, oxycodone, effect_adjusted, effect_unadjusted, p_adjusted),
      n = Inf, width = Inf)
cat("\nNotes column (pair OR / paired Wilcoxon / caveats):\n")
print(res |> filter(note != "") |> select(outcome, analysis, note), n = Inf, width = Inf)

out_file <- sprintf("outcomes_final_summary_%s.csv", tag)
write_csv(res, out_file)
cat(sprintf("\nSaved: %s\n", out_file))
cat("\nREAD: effects are oxycodone vs morphine. For 'no difference' claims quote the 95% CI,\n")
cat("not the p-value. Run once with STRICT_FIRST_EVER <- TRUE (primary) and once with FALSE\n")
cat("(sensitivity; re-run scripts 6 and 10 first) and compare the two CSVs.\n")