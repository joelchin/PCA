# ============================================================
# 10 - Time to effect (post hoc; Figure 2)
# ============================================================
# Checks that comparable consumption reflects comparable analgesia: the
# time from episode start to a sustained pain score of 3 or below in the
# propensity-matched cohort.
#
# Definition
#   Event:      first pain score <= 3 followed by a second score <= 3 (a
#               final charted score <= 3 also counts)
#   Population: matched episodes in children aged 6 years or over whose
#               first charted score was above 3 (younger children have
#               many zero scores)
#   Censoring:  at episode end if the threshold is not reached
#   Exclusions: episodes with no pain readings; readings sharing a
#               timestamp are collapsed to the highest score
#   Analysis:   Kaplan-Meier curves; Cox regression for oxycodone vs
#               morphine adjusted for starting pain score, with robust
#               standard errors clustered on the matched pair. Morphine is
#               the reference, so HR > 1 means oxycodone reaches the
#               threshold sooner.
#
# The threshold of 3 is the upper end of the mild range used for
# paediatric PCA endpoints (Serlin et al. 1995 define mild as 1 to 4).
#
# Run after 09_comparative_psm.R in the same R session (needs
# matched_data and pain_scores_clean). survminer must be installed but is
# not attached, because it loads ggpubr, whose font() masks
# flextable::font() and would break the Word export in
# 12_publication_tables.R.
# Writes time_to_effect_KM_figure.png and .pdf to the working directory.
# ============================================================

if (!requireNamespace("survminer", quietly = TRUE)) {
  stop("survminer is not installed. Install it in a fresh R session ",
       "(Session > Restart R), then run the pipeline again.")
}

library(tidyverse)
library(survival)

stopifnot(
  exists("pca_episodes_flagged"),
  exists("matched_data"),
  exists("pain_scores_clean"),
  "subclass" %in% names(matched_data)
)

THRESHOLD <- 3
AGE_MIN   <- 6      # pain-score efficacy analyses use age >= 6
KM_XMAX   <- 48     # display limit for the figure (hours); statistics use all follow-up

# ============================================================
# Episode frame: matched episodes aged >= 6 with a valid window
# ============================================================

s3_keys <- matched_data |>
  select(PAT_ENC_CSN_ID, episode_id, drug, subclass) |>
  distinct()

episode_frame <- pca_episodes_flagged |>
  inner_join(s3_keys, by = c("PAT_ENC_CSN_ID", "episode_id", "drug")) |>
  select(PAT_ENC_CSN_ID, episode_id, drug, AGE, episode_start, episode_end, subclass) |>
  mutate(episode_duration_hrs = as.numeric(difftime(episode_end, episode_start, units = "hours"))) |>
  filter(episode_duration_hrs > 0)

# Matching used all ages, so the age rule can leave a pair with one member
cat(sprintf("Matched episodes with a valid window: %d | aged >= %g: %d | younger or age missing: %d\n",
            nrow(episode_frame), AGE_MIN,
            sum(!is.na(episode_frame$AGE) & episode_frame$AGE >= AGE_MIN),
            sum(is.na(episode_frame$AGE) | episode_frame$AGE < AGE_MIN)))
episode_frame <- episode_frame |> filter(!is.na(AGE), AGE >= AGE_MIN)

# ============================================================
# Pain readings within each episode
# ============================================================

readings <- episode_frame |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end) |>
  left_join(
    pain_scores_clean |> select(PAT_ENC_CSN_ID, score_time, score_value),
    by = "PAT_ENC_CSN_ID", relationship = "many-to-many"
  ) |>
  filter(!is.na(score_time), !is.na(score_value),
         score_time >= episode_start, score_time <= episode_end) |>
  arrange(PAT_ENC_CSN_ID, episode_id, score_time)

# Episodes with no pain readings carry no information and are excluded
evaluable_keys <- readings |> distinct(PAT_ENC_CSN_ID, episode_id)
cat(sprintf("Episodes: %d | with at least one reading: %d | without (excluded): %d\n",
            nrow(episode_frame), nrow(evaluable_keys),
            nrow(episode_frame) - nrow(evaluable_keys)))
episode_frame <- episode_frame |>
  semi_join(evaluable_keys, by = c("PAT_ENC_CSN_ID", "episode_id"))

# Readings charted at the same timestamp are collapsed to the highest
# score (a concurrent score of 6 means the child was not yet comfortable)
readings_unique <- readings |>
  group_by(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end, score_time) |>
  summarise(score_value = max(score_value), .groups = "drop") |>
  arrange(PAT_ENC_CSN_ID, episode_id, score_time)

cat(sprintf("Readings before collapsing identical timestamps: %d | after: %d\n\n",
            nrow(readings), nrow(readings_unique)))

# First charted score of each episode
initial_scores <- readings_unique |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  slice_min(score_time, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(PAT_ENC_CSN_ID, episode_id, initial_score = score_value)

# ============================================================
# Time to a sustained score <= threshold
# ============================================================

# Only episodes whose first score was above the threshold are analysed:
# an episode starting at or below it would trivially "reach" it.
tte_data <- episode_frame |>
  inner_join(initial_scores, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  filter(initial_score > THRESHOLD)

# Sustained: the next reading in the episode, if there is one, is also <= threshold
first_sustained <- readings_unique |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  mutate(next_score_value = lead(score_value)) |>
  ungroup() |>
  filter(score_value <= THRESHOLD,
         is.na(next_score_value) | next_score_value <= THRESHOLD) |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  slice_min(score_time, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(time_to_threshold_hrs = as.numeric(difftime(score_time, episode_start, units = "hours"))) |>
  select(PAT_ENC_CSN_ID, episode_id, time_to_threshold_hrs)

tte_data <- tte_data |>
  left_join(first_sustained, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(
    drug = factor(drug, levels = c("morphine", "oxycodone")),   # morphine = reference
    event = as.numeric(!is.na(time_to_threshold_hrs)),
    time_hrs = if_else(event == 1, time_to_threshold_hrs, episode_duration_hrs)
  )

stopifnot(!anyDuplicated(tte_data[c("PAT_ENC_CSN_ID", "episode_id")]))

n_complete_pairs <- tte_data |> count(subclass) |> filter(n == 2) |> nrow()
cat(sprintf("Analysed: %d episodes, %d in complete matched pairs\n",
            nrow(tte_data), n_complete_pairs))
cat(sprintf("Reaching the threshold before episode end: %d / %d (%.1f%%)\n",
            sum(tte_data$event), nrow(tte_data), 100 * mean(tte_data$event)))
print(tte_data |>
        group_by(drug) |>
        summarise(N = n(), n_reached = sum(event),
                  pct_reached = round(100 * mean(event), 1), .groups = "drop"))

# ============================================================
# Kaplan-Meier and Cox models
# ============================================================

fit <- survfit(Surv(time_hrs, event) ~ drug, data = tte_data)
cat("\nMedian time to threshold (hours, Kaplan-Meier), by drug:\n")
print(summary(fit)$table[, "median"])

# Robust standard errors clustered on the matched pair respect the matching
cox_unadj <- coxph(Surv(time_hrs, event) ~ drug + cluster(subclass), data = tte_data)
cox_adj   <- coxph(Surv(time_hrs, event) ~ drug + initial_score + cluster(subclass),
                   data = tte_data)

hr_text <- function(m) {
  ci <- summary(m)$conf.int
  r <- grep("^drug", rownames(ci))[1]
  sprintf("HR %.2f (95%% CI %.2f-%.2f)", ci[r, "exp(coef)"], ci[r, "lower .95"], ci[r, "upper .95"])
}
hr_adj <- hr_text(cox_adj)

cat("\nOxycodone vs morphine (HR > 1: oxycodone reaches the threshold sooner):\n")
cat("  unadjusted:                      ", hr_text(cox_unadj), "\n")
cat("  adjusted for starting pain score:", hr_adj, "\n")
ph_p <- tryCatch(cox.zph(cox_adj)$table["drug", "p"], error = function(e) NA_real_)
cat(sprintf("  proportional-hazards check for drug: p = %s\n\n",
            ifelse(is.na(ph_p), "n/a", sprintf("%.3f", ph_p))))

# ============================================================
# Figure 2
# ============================================================
# Cumulative proportion reaching the threshold (curves rise) with 95%
# confidence bands and numbers at risk, annotated with the adjusted HR.
# The PNG is written with ragg when available, and the PDF with
# cairo_pdf(), because base pdf() cannot draw the <= sign.

tryCatch({
  hr_label <- paste0(
    "Oxycodone vs morphine\n",
    "Adjusted ", hr_adj, "\n",
    "(HR > 1: oxycodone reaches the threshold sooner;\n",
    "SE clustered on matched pair)"
  )
  y_lab <- sprintf("Proportion reaching sustained pain score \u2264 %d", THRESHOLD)
  
  p_km <- survminer::ggsurvplot(
    fit, data = tte_data,
    fun = "event",
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
  
  p_km$plot <- p_km$plot +
    ggplot2::annotate("text", x = KM_XMAX - 1, y = 0.28, label = hr_label,
                      hjust = 1, vjust = 1, size = 3.4, lineheight = 1.05) +
    ggplot2::theme(legend.position = "top",
                   legend.text = ggplot2::element_text(size = 11))
  
  # print(..., newpage = FALSE) avoids a blank first page
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png("time_to_effect_KM_figure.png", width = 7, height = 7, units = "in", res = 300)
  } else {
    png("time_to_effect_KM_figure.png", width = 7, height = 7, units = "in", res = 300, type = "cairo")
  }
  print(p_km, newpage = FALSE)
  dev.off()
  
  pdf_ok <- tryCatch({
    cairo_pdf("time_to_effect_KM_figure.pdf", width = 7, height = 7)
    print(p_km, newpage = FALSE)
    dev.off()
    TRUE
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    FALSE
  })
  cat(if (pdf_ok) "Saved time_to_effect_KM_figure.png and .pdf\n\n"
      else "Saved time_to_effect_KM_figure.png (cairo_pdf() unavailable; PDF not written)\n\n")
  
  n_by_arm <- table(tte_data$drug)
  cat("Figure legend numbers:\n")
  cat(sprintf(paste0(
    "n = %d episodes (morphine %d, oxycodone %d), %d complete matched pairs; ",
    "children aged %g years or over whose first charted pain score was above %d.\n"),
    nrow(tte_data), n_by_arm[["morphine"]], n_by_arm[["oxycodone"]],
    n_complete_pairs, AGE_MIN, THRESHOLD))
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  cat("Figure generation skipped: ", conditionMessage(e), "\n",
      "(survminer must be installed; the statistics above are unaffected.)\n", sep = "")
})
