# ============================================================
# TIME TO EFFECT: how long until PCA/NCA achieves adequate analgesia?
# ============================================================
# Motivation: PCA/NCA is titrated to effect, which is why raw opioid
# consumption doesn't need an OME conversion to be interpretable
# (Comment 7). But "titrated to effect" implies there IS a moment
# effect is achieved -- and consumption only makes full sense once
# you know how quickly that happened. This defines a time-to-event
# endpoint: time from episode start to the FIRST pain score at or
# below a threshold, compared between oxycodone and morphine.
#
# Right-censoring: an episode that never reaches the threshold before
# ending is censored at episode end (effect not achieved within the
# observed window), handled properly via Kaplan-Meier / log-rank
# rather than just dropping those episodes or treating them as missing.
#
# As with every other arbitrary cutoff in this analysis, we test a
# RANGE of thresholds rather than committing to one.
#
# REVISION (this pass)
#   - Loading-dose adjustment and loading-dose strata removed: the
#     prescribed loading dose on the order is a "just in case" contingency,
#     not a dose known to have been given, so it is not used as an
#     exposure or covariate anywhere in the analysis.
#   - Section 3 (titration-phase consumption) removed: it multiplied the
#     whole-episode mean rate by the hours to threshold, which understates
#     front-loaded early dosing. Consumption is reported from the matched
#     total-dose analysis instead.
#   - The PRIMARY definition is now explicit (threshold <=3, sustained,
#     episodes starting above threshold) and the summary lists it first;
#     all other combinations are sensitivity. The Cox hazard ratio with
#     pair-clustered robust SEs is the primary inference (the log-rank test
#     does not account for the matching).
#   - The number of complete matched pairs among the analysed episodes is
#     printed for every definition (the age >= 6 rule and the at-risk
#     restriction both break some pairs).
#   - FIGURE EXPORT FIXED (Section 4). Two problems with the previous
#     export, neither affecting the statistics:
#       1. Blank first page in the PDF: print() on a survminer plot starts
#          a new page by default, leaving the device's own first page empty.
#          Now print(..., newpage = FALSE).
#       2. The "<=" sign came out as "=" in the PDF: pdf() uses a Latin-1
#          font with no such glyph. The PDF is now written with
#          cairo_pdf(), and the PNG with ragg::agg_png() when available
#          (falls back to png(type = "cairo")), both of which draw Unicode
#          correctly.
# ============================================================

if (!requireNamespace("survminer", quietly = TRUE)) {
  stop("survminer is not installed. Install it in a fresh R session (Session > Restart R) ",
       "BEFORE running the pipeline, then run the pipeline again.")
}

library(tidyverse)
library(survival)
# survminer is NOT attached (it loads ggpubr, whose font() masks flextable::font() and
# breaks the Word export in script 16 if both are in the same session). The KM plot
# below calls survminer::ggsurvplot() directly instead; the stats work without it.

cat("============================================================\n")
cat("TIME TO EFFECT ANALYSIS\n")
cat("============================================================\n\n")

stopifnot(exists("pca_episodes_flagged"))
stopifnot(exists("matched_data"))
stopifnot(exists("pain_scores_clean"))
stopifnot(exists("s4_surgical_keys"))

THRESHOLDS <- c(3, 4, 5)  # test a range -- no single "correct" cutoff

# PRIMARY definition (pre-specified; everything else is sensitivity):
# pain score <= 3, sustained (next reading also <= 3), restricted to
# episodes whose first charted score was above 3.
PRIMARY_THRESHOLD <- 3
PRIMARY_DEF <- "SUSTAINED | AT-RISK ONLY"

# Minimum age (years) for pain-score efficacy analyses. 6 follows the locked
# study rule. Set to 0 to include all ages (not recommended for pain scores).
# To compare rules without editing this file, set AGE_MIN_OVERRIDE (e.g. 0 for
# all ages) in the console before sourcing; rm(AGE_MIN_OVERRIDE) to go back.
AGE_MIN <- if (exists("AGE_MIN_OVERRIDE")) AGE_MIN_OVERRIDE else 6
cat(sprintf("Age rule in use: pain-score analyses include age >= %g years\n\n", AGE_MIN))

# Threshold rationale: Serlin et al. 1995 (Pain 61:277-84) established
# NRS 1-4 = mild pain on a 0-10 scale, the standard cited cutpoint
# system. A comparable active paediatric PCA trial (adolescent
# scoliosis, NCT07229495) uses NRS<=3 as its "adequate analgesia"
# efficacy endpoint and NRS>=4 as the threshold defining inadequate
# control requiring rescue analgesia -- i.e. the literature treats
# <=3 and >=4 as two framings of the same boundary. 3 is used as
# PRIMARY below; 4-5 are tested as sensitivity, reflecting Zernikow
# et al. 2013's finding that "optimal" cutpoints vary across
# populations (no single universally correct value).

# ============================================================
# SECTION 0: BUILD EPISODE + READING FRAME (matched cohort primary)
# ============================================================

has_subclass <- "subclass" %in% names(matched_data)
cat(if (has_subclass) "Matched-pair ID (subclass) found: Cox models use robust SEs clustered on pair.\n\n"
    else "No 'subclass' column on matched_data: Cox models will treat arms as independent.\n\n")
s3_keys <- matched_data |>
  select(PAT_ENC_CSN_ID, episode_id, drug, any_of("subclass")) |>
  distinct()

episode_frame <- pca_episodes_flagged |>
  inner_join(s3_keys, by = c("PAT_ENC_CSN_ID", "episode_id", "drug")) |>
  select(PAT_ENC_CSN_ID, episode_id, drug, AGE, episode_start, episode_end, any_of("subclass")) |>
  mutate(episode_duration_hrs = as.numeric(difftime(episode_end, episode_start, units = "hours"))) |>
  filter(episode_duration_hrs > 0)

# Age restriction (locked rule: pain-score efficacy analyses use age >= 6,
# because younger children have implausibly high zero-score rates, which
# would inflate "comfortable" time and the share starting below threshold).
# The matching itself used all ages (neonates included), so this is a
# subset of the matched cohort and can leave some pairs with one member.
cat(sprintf("Matched episodes with a valid window: %d | aged >= %g: %d | younger or age missing (excluded): %d\n\n",
            nrow(episode_frame), AGE_MIN,
            sum(!is.na(episode_frame$AGE) & episode_frame$AGE >= AGE_MIN),
            sum(is.na(episode_frame$AGE) | episode_frame$AGE < AGE_MIN)))
episode_frame <- episode_frame |> filter(!is.na(AGE), AGE >= AGE_MIN)

if (has_subclass) {
  n_pairs_age <- episode_frame |> count(subclass) |> filter(n == 2) |> nrow()
  cat(sprintf("After the age rule: %d episodes, of which %d are in complete matched pairs (%d episodes have lost their partner)\n\n",
              nrow(episode_frame), n_pairs_age, nrow(episode_frame) - 2 * n_pairs_age))
}

readings <- episode_frame |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end) |>
  left_join(
    pain_scores_clean |> select(PAT_ENC_CSN_ID, score_time, score_value),
    by = "PAT_ENC_CSN_ID", relationship = "many-to-many"
  ) |>
  filter(!is.na(score_time), !is.na(score_value),
         score_time >= episode_start, score_time <= episode_end) |>
  arrange(PAT_ENC_CSN_ID, episode_id, score_time)

cat(sprintf("Matched-cohort episodes with pain score data: %d\n\n",
            n_distinct(paste(readings$PAT_ENC_CSN_ID, readings$episode_id))))

# ---- FIX 1: drop episodes with NO pain readings at all -------------------
# Previously these were kept and treated as "never reached threshold",
# censored at full episode length -- but a missing score is not evidence
# that comfort was never achieved. They carry no information, so they
# are excluded here.
evaluable_keys <- readings |> distinct(PAT_ENC_CSN_ID, episode_id)
cat(sprintf("Matched episodes: %d | with >=1 reading (analysed): %d | no readings (excluded): %d\n",
            nrow(episode_frame), nrow(evaluable_keys),
            nrow(episode_frame) - nrow(evaluable_keys)))
episode_frame <- episode_frame |>
  semi_join(evaluable_keys, by = c("PAT_ENC_CSN_ID", "episode_id"))

# ---- FIX 2: collapse readings sharing an identical timestamp -------------
# Many flowsheet entries are charted at the SAME timestamp (earlier
# output: median time-to-next-reading was 0 h). With ties, (a) the
# 'first reading <= threshold' step returned several rows per episode,
# silently duplicating episodes in the survival data (946 rows analysed
# from ~468 matched episodes), and (b) the 'sustained' definition's
# "next reading" was often just a simultaneous duplicate. We now keep ONE
# value per distinct timestamp per episode -- the HIGHEST score charted
# at that moment (conservative: a concurrent score of 6 means the child
# was not yet comfortable).
readings_unique <- readings |>
  group_by(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end, score_time) |>
  summarise(score_value = max(score_value), .groups = "drop") |>
  arrange(PAT_ENC_CSN_ID, episode_id, score_time)

cat(sprintf("Readings before collapsing identical timestamps: %d | after: %d\n\n",
            nrow(readings), nrow(readings_unique)))

# Pain score at the FIRST charted moment of each episode (used to define
# who actually needed titrating -- see at_risk_only below)
initial_scores <- readings_unique |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  slice_min(score_time, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(PAT_ENC_CSN_ID, episode_id, initial_score = score_value)

cat("How many episodes START above each threshold (i.e. genuinely needed\n")
cat("titrating) -- if most start at/below threshold, 'time to threshold'\n")
cat("mostly measures time to first charted score, not analgesic effect:\n\n")
for (th0 in THRESHOLDS) {
  tmp <- episode_frame |>
    inner_join(initial_scores, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    group_by(drug) |>
    summarise(N = n(), n_start_above = sum(initial_score > th0),
              pct_start_above = round(100 * mean(initial_score > th0), 1),
              .groups = "drop")
  cat(sprintf("Threshold <= %d:\n", th0)); print(tmp, width = Inf); cat("\n")
}

# ============================================================
# SECTION 1: COMPUTE TIME-TO-THRESHOLD PER EPISODE, PER THRESHOLD
# ============================================================

compute_time_to_threshold <- function(threshold, sustained = FALSE, at_risk_only = FALSE) {
  # sustained = FALSE: "first touch" -- the first single reading at or
  #   below threshold counts as the event, even if pain climbs back up
  #   immediately after. Vulnerable to a single noisy/transient dip
  #   (e.g. a brief lull between boluses) being miscounted as "effect
  #   achieved" when the pain trajectory isn't actually settling.
  # sustained = TRUE: requires the NEXT reading (if one exists within
  #   the episode) to ALSO be at or below threshold -- a much better
  #   proxy for genuinely achieving and holding comfort, given pain
  #   trajectories are often non-monotonic (physio, position changes,
  #   breakthrough pain) rather than cleanly declining.
  # at_risk_only = TRUE: restrict to episodes whose FIRST charted score was
  #   ABOVE the threshold -- i.e. children who actually had something to
  #   titrate. Without this, an episode that starts at 0 trivially
  #   "reaches" the threshold at its first reading.
  frame <- episode_frame |>
    inner_join(initial_scores, by = c("PAT_ENC_CSN_ID", "episode_id"))
  if (at_risk_only) {
    frame <- frame |> filter(initial_score > threshold)
  }
  
  ordered_readings <- readings_unique |>
    group_by(PAT_ENC_CSN_ID, episode_id) |>
    mutate(next_score_value = lead(score_value)) |>
    ungroup()
  
  if (sustained) {
    candidates <- ordered_readings |>
      filter(score_value <= threshold,
             is.na(next_score_value) | next_score_value <= threshold)
  } else {
    candidates <- ordered_readings |>
      filter(score_value <= threshold)
  }
  
  first_below <- candidates |>
    group_by(PAT_ENC_CSN_ID, episode_id) |>
    slice_min(score_time, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(time_to_threshold_hrs = as.numeric(difftime(score_time, episode_start, units = "hours"))) |>
    select(PAT_ENC_CSN_ID, episode_id, time_to_threshold_hrs)
  
  out <- frame |>
    left_join(first_below, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
    mutate(
      reached_threshold = !is.na(time_to_threshold_hrs),
      # Censored at episode end if threshold never (sustainably) reached
      time_hrs = if_else(reached_threshold, time_to_threshold_hrs, episode_duration_hrs),
      event = as.numeric(reached_threshold)
    )
  
  # Guard: exactly ONE row per episode (the bug above would trip this)
  stopifnot(!anyDuplicated(out[c("PAT_ENC_CSN_ID", "episode_id")]))
  out
}

# ============================================================
# SECTION 2: KAPLAN-MEIER + LOG-RANK, ACROSS THRESHOLDS
# ============================================================

cat("SECTION 2: TIME TO FIRST REACHING THRESHOLD, BY DRUG\n")
cat("============================================================\n\n")

km_results <- list()

# Matched samples are not independent observations: each oxycodone episode
# has a partner. Robust (sandwich) SEs clustered on the matched pair are the
# standard way to respect that for time-to-event outcomes (Austin 2014).
f_unadj <- if (has_subclass) Surv(time_hrs, event) ~ drug + cluster(subclass) else Surv(time_hrs, event) ~ drug
f_adj   <- if (has_subclass) Surv(time_hrs, event) ~ drug + initial_score + cluster(subclass) else Surv(time_hrs, event) ~ drug + initial_score

# Hazard ratio helper. drug is alphabetical, so MORPHINE is the reference:
# HR > 1 = oxycodone reaches the threshold FASTER; HR < 1 = more slowly.
hr_text <- function(m) {
  tryCatch({
    ci <- summary(m)$conf.int
    r <- grep("^drug", rownames(ci))[1]
    sprintf("HR %.2f (95%% CI %.2f-%.2f)", ci[r, "exp(coef)"], ci[r, "lower .95"], ci[r, "upper .95"])
  }, error = function(e) "HR n/a")
}

# Each threshold is run under BOTH definitions:
#   "first touch" = first single reading at/below threshold
#   "sustained"   = first reading at/below threshold whose NEXT reading
#                   also stays at/below it (guards against transient dips)
definitions <- expand_grid(th = THRESHOLDS, at_risk = c(TRUE, FALSE), sustained = c(FALSE, TRUE))

for (i in seq_len(nrow(definitions))) {
  th <- definitions$th[i]
  sus <- definitions$sustained[i]
  ar <- definitions$at_risk[i]
  def_label <- paste0(if (sus) "SUSTAINED" else "FIRST TOUCH",
                      if (ar) " | AT-RISK ONLY" else " | ALL EPISODES")
  key <- paste0(th, "_", def_label)
  is_primary <- (th == PRIMARY_THRESHOLD && def_label == PRIMARY_DEF)
  
  cat(sprintf("--- Threshold: pain score <= %d | Definition: %s%s ---\n\n",
              th, def_label, if (is_primary) "  [PRIMARY]" else ""))
  
  tte_data <- compute_time_to_threshold(th, sustained = sus, at_risk_only = ar)
  
  n_complete_pairs <- if (has_subclass) tte_data |> count(subclass) |> filter(n == 2) |> nrow() else NA_integer_
  cat(sprintf("Episodes analysed: %d | complete matched pairs among them: %s\n",
              nrow(tte_data), ifelse(is.na(n_complete_pairs), "n/a", as.character(n_complete_pairs))))
  
  cat(sprintf("Episodes reaching threshold before episode end: %d / %d (%.1f%%)\n",
              sum(tte_data$event), nrow(tte_data),
              100 * sum(tte_data$event) / nrow(tte_data)))
  
  by_drug_reach <- tte_data |>
    group_by(drug) |>
    summarise(N = n(), n_reached = sum(event),
              pct_reached = round(100 * mean(event), 1), .groups = "drop")
  print(by_drug_reach, width = Inf)
  cat("\n")
  
  fit <- survfit(Surv(time_hrs, event) ~ drug, data = tte_data)
  cat("Median time to threshold (from KM estimate), by drug:\n")
  print(summary(fit)$table[, "median"])
  cat("\n")
  
  lr <- survdiff(Surv(time_hrs, event) ~ drug, data = tte_data)
  p_val <- 1 - pchisq(lr$chisq, df = 1)
  cat(sprintf("Log-rank test (does not account for the matching; the Cox HR below is the primary inference): chi-sq = %.3f, p = %.4f\n", lr$chisq, p_val))
  cat(if_else(p_val > 0.05,
              "-> No significant difference in TIME TO EFFECT between drugs\n\n",
              "-> SIGNIFICANT difference in time to effect between drugs -- investigate direction\n\n"))
  
  cox_unadj <- tryCatch(coxph(f_unadj, data = tte_data), error = function(e) NULL)
  cox_adj   <- tryCatch(coxph(f_adj, data = tte_data), error = function(e) NULL)
  ph_p <- tryCatch(cox.zph(cox_adj)$table["drug", "p"], error = function(e) NA_real_)
  hr_unadj <- if (is.null(cox_unadj)) "HR n/a" else hr_text(cox_unadj)
  hr_adj   <- if (is.null(cox_adj))   "HR n/a" else hr_text(cox_adj)
  cat(sprintf("Cox, oxycodone vs morphine, robust SE clustered on matched pair (HR>1 = oxycodone reaches threshold faster):\n  unadjusted: %s\n  adjusted for starting pain score: %s\n  proportional-hazards check for drug: p = %s (small p = curves cross / HR not constant)\n\n",
              hr_unadj, hr_adj, ifelse(is.na(ph_p), "n/a", sprintf("%.3f", ph_p))))
  
  km_results[[key]] <- list(fit = fit, data = tte_data, p = p_val, th = th, def = def_label,
                            hr_unadj = hr_unadj, hr_adj = hr_adj,
                            n_eps = nrow(tte_data), n_pairs = n_complete_pairs)
}

# ============================================================
# SECTION 3: SUMMARY ACROSS THRESHOLDS
# ============================================================

cat("SECTION 3: SUMMARY ACROSS THRESHOLDS\n")
cat("============================================================\n\n")

prim_key <- paste0(PRIMARY_THRESHOLD, "_", PRIMARY_DEF)
r <- km_results[[prim_key]]
cat("PRIMARY (pre-specified):\n")
cat(sprintf("  Threshold <=%d | %s | episodes %d, complete pairs %s\n    adjusted %s | unadjusted %s | log-rank p = %.4f\n\n",
            r$th, r$def, r$n_eps,
            ifelse(is.na(r$n_pairs), "n/a", as.character(r$n_pairs)),
            r$hr_adj, r$hr_unadj, r$p))

cat("Sensitivity definitions:\n")
for (k in setdiff(names(km_results), prim_key)) {
  r <- km_results[[k]]
  cat(sprintf("  Threshold <=%d | %-34s: episodes %d | log-rank p = %.4f | adjusted %s\n",
              r$th, r$def, r$n_eps, r$p, r$hr_adj))
}

cat("\nIf time-to-effect does not differ between drugs across thresholds,\n")
cat("this directly supports interpreting whole-episode consumption (already\n")
cat("shown comparable, Comment 7) as reflecting a genuinely comparable\n")
cat("titration process -- not two different speeds of reaching comfort\n")
cat("that happen to average out to the same total dose. Quote the hazard\n")
cat("ratio with its 95% CI, not just the p-value.\n")

cat("\n============================================================\n")

# ============================================================
# SECTION 4: MANUSCRIPT FIGURE (primary definition)
# ============================================================
# Publication-style Kaplan-Meier figure for the PRIMARY definition
# (pain score <= 3, sustained, episodes starting above 3):
#   * colour-blind-safe palette, light confidence bands, capitalised labels
#   * numbers at risk beneath the curves
#   * annotated with the adjusted Cox hazard ratio (robust SE clustered on
#     the matched pair) instead of the log-rank p-value, which ignores the
#     matching
#   * x-axis limited to 48 h (display only: the statistics use all follow-up)
# Saves time_to_effect_KM_figure.png (300 dpi) and .pdf, and prints a
# suggested figure legend with the actual n.
#
# KM_AS_CUMULATIVE = TRUE plots the proportion who HAVE reached the
# threshold (1 - survival; curves rise, higher = more comfortable).
# FALSE plots the proportion NOT yet reaching it (curves fall).

KM_AS_CUMULATIVE <- TRUE
KM_XMAX <- 48

tryCatch({
  r_prim <- km_results[[prim_key]]
  primary_fit  <- r_prim$fit
  primary_data <- r_prim$data
  
  hr_label <- paste0(
    "Oxycodone vs morphine\n",
    "Adjusted ", r_prim$hr_adj, "\n",
    "(HR > 1: oxycodone reaches the threshold sooner;\n",
    "SE clustered on matched pair)"
  )
  
  # "\u2264" is the less-than-or-equal sign. It is drawn correctly by the
  # PNG and PDF devices used below (ragg / cairo), but NOT by base pdf().
  y_lab <- if (KM_AS_CUMULATIVE) {
    sprintf("Proportion reaching sustained pain score \u2264 %d", PRIMARY_THRESHOLD)
  } else {
    sprintf("Proportion not yet reaching pain score \u2264 %d", PRIMARY_THRESHOLD)
  }
  
  p_km <- survminer::ggsurvplot(
    primary_fit, data = primary_data,
    fun = if (KM_AS_CUMULATIVE) "event" else NULL,
    conf.int = TRUE, conf.int.alpha = 0.12,
    pval = FALSE,
    risk.table = TRUE, risk.table.height = 0.24,
    risk.table.title = "Number at risk",
    risk.table.y.text = FALSE,
    palette = c("#D55E00", "#0072B2"),          # morphine, oxycodone
    xlim = c(0, KM_XMAX), break.time.by = 6,
    xlab = "Hours since PCA/NCA start",
    ylab = y_lab,
    legend.title = "",
    legend.labs = c("Morphine", "Oxycodone"),
    title = NULL,
    ggtheme = ggplot2::theme_classic(base_size = 12),
    tables.theme = survminer::theme_cleantable(),
    fontsize = 3.5
  )
  
  annot_y <- if (KM_AS_CUMULATIVE) 0.28 else 0.85
  p_km$plot <- p_km$plot +
    ggplot2::annotate("text", x = KM_XMAX - 1, y = annot_y, label = hr_label,
                      hjust = 1, vjust = 1, size = 3.4, lineheight = 1.05) +
    ggplot2::theme(legend.position = "top",
                   legend.text = ggplot2::element_text(size = 11))
  
  # ---- Export ---------------------------------------------------------
  # print(..., newpage = FALSE): survminer's print method otherwise calls
  # grid.newpage() on a freshly opened device, which leaves the first page
  # blank (the blank first page seen in the old PDF).
  
  # PNG, 300 dpi. ragg handles Unicode glyphs; fall back to Cairo if ragg
  # is not installed.
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png("time_to_effect_KM_figure.png", width = 7, height = 7, units = "in", res = 300)
  } else {
    png("time_to_effect_KM_figure.png", width = 7, height = 7, units = "in", res = 300, type = "cairo")
  }
  print(p_km, newpage = FALSE)
  dev.off()
  
  # PDF. cairo_pdf() embeds the real "<=" character; base pdf() cannot, which
  # is why the old PDF showed "=" instead.
  pdf_ok <- tryCatch({
    cairo_pdf("time_to_effect_KM_figure.pdf", width = 7, height = 7)
    print(p_km, newpage = FALSE)
    dev.off()
    TRUE
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    FALSE
  })
  
  if (pdf_ok) {
    cat("Saved: time_to_effect_KM_figure.png and time_to_effect_KM_figure.pdf\n\n")
  } else {
    cat("Saved: time_to_effect_KM_figure.png\n",
        "PDF NOT written: cairo_pdf() is unavailable in this R session. Use the PNG\n",
        "(it draws the <= sign correctly) or tell me and I will give a different PDF route.\n\n", sep = "")
  }
  
  n_by_arm <- table(primary_data$drug)
  cat("Suggested figure legend:\n")
  cat(sprintf(paste0(
    "Figure X. Time to sustained pain score \u2264 %d in propensity-matched children aged %g years or over ",
    "whose first charted pain score was above %d (n = %d episodes: morphine %d, oxycodone %d; ",
    "%s complete matched pairs). 'Sustained' = a reading \u2264 %d followed by a further reading \u2264 %d. ",
    "Shaded bands are 95%% confidence intervals; the x-axis is truncated at %d hours for display. ",
    "The hazard ratio is from a Cox model adjusted for starting pain score with robust standard errors ",
    "clustered on the matched pair.\n"),
    PRIMARY_THRESHOLD, AGE_MIN, PRIMARY_THRESHOLD, nrow(primary_data),
    n_by_arm[["morphine"]], n_by_arm[["oxycodone"]],
    ifelse(is.na(r_prim$n_pairs), "n/a", as.character(r_prim$n_pairs)),
    PRIMARY_THRESHOLD, PRIMARY_THRESHOLD, KM_XMAX))
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  cat("Figure generation skipped: ", conditionMessage(e), "\n",
      "(survminer must be installed; the statistics above are unaffected.)\n", sep = "")
})