# ============================================================
# 5 — Load cleaned data into the R session
# Paediatric PCA/NCA study — Addenbrooke's Hospital (CUH NHS FT)
# Author: J Chin
# ============================================================
#
# Purpose
#   Start-of-session loader: reads every cleaned table and episode-level
#   object produced by scripts 2-4 from disk into the R session. Run this
#   first when returning to the analysis without re-running the whole
#   pipeline.
#
# Inputs
#   drug_lists.R and the .rds files in "1 - data/3 - processed_data/"
#   written by 2_cleaningscript.R, 3_shift_table.R and 4_flags.R.
#
# Objects created
#   admissions_clean, pca_mar_clean, non_pca_mar_clean, intraop_mar_clean,
#   nerve_blocks_clean, operations_clean, anaesthesia_events_clean,
#   anaesthesia_ax_clean, pain_scores_clean, lda_urethral_clean, lda_epidural_clean,
#   lda_epidural_assess_clean, shifts, episode_shift_table, episode_metrics,
#   pca_episodes_flagged (single source of truth for episodes) and
#   pca_episodes (an alias of pca_episodes_flagged kept for older scripts).
#
# Not loaded
#   Objects built in-session by later scripts (psm_data, matched_data,
#   match_out, tbl5/ft5*) are not saved to disk. Re-run 10_comparative_PSM.R in
#   the same R session before 15, and 15 before 16.
# ============================================================

library(tidyverse)
library(lubridate)
library(furrr)
library(future)

# ---- Project root (see 1_load_raw_data.R) ---------------------------------
proj_dir <- Sys.getenv("PCA_PROJECT_DIR", unset = "")
if (nzchar(proj_dir)) setwd(proj_dir)
if (!dir.exists("1 - data")) {
  stop("Project root not found. Set the working directory (or PCA_PROJECT_DIR) ",
       "to the folder that contains '1 - data/'.")
}

# ---- Drug lists (kept in one place; edit drug_lists.R, not the scripts) ---
source("drug_lists.R")

# ---- Cleaned tables ----------------------------------------------------------
admissions_clean          <- readRDS("1 - data/3 - processed_data/admissions_clean.rds")
pca_mar_clean             <- readRDS("1 - data/3 - processed_data/pca_mar_clean.rds")
non_pca_mar_clean         <- readRDS("1 - data/3 - processed_data/non_pca_mar_clean.rds")
intraop_mar_clean         <- readRDS("1 - data/3 - processed_data/intraop_mar_clean.rds")
nerve_blocks_clean        <- readRDS("1 - data/3 - processed_data/nerve_blocks_clean.rds")
operations_clean          <- readRDS("1 - data/3 - processed_data/operations_clean.rds")
anaesthesia_events_clean  <- readRDS("1 - data/3 - processed_data/anaesthesia_events_clean.rds")
anaesthesia_ax_clean      <- readRDS("1 - data/3 - processed_data/anaesthesia_ax_clean.rds")
pain_scores_clean         <- readRDS("1 - data/3 - processed_data/pain_scores_clean.rds")
shifts                    <- readRDS("1 - data/3 - processed_data/1_shifts.rds")
episode_shift_table       <- readRDS("1 - data/3 - processed_data/episode_shift_table.rds")
episode_metrics           <- readRDS("1 - data/3 - processed_data/episode_metrics.rds")

# pca_episodes_flagged is the single source of truth on disk.
# pca_episodes is kept as an alias because older scripts still use that name.
pca_episodes_flagged      <- readRDS("1 - data/3 - processed_data/pca_episodes_flagged.rds")
pca_episodes              <- pca_episodes_flagged

lda_urethral_clean        <- readRDS("1 - data/3 - processed_data/lda_urethral_clean.rds")
lda_epidural_clean        <- readRDS("1 - data/3 - processed_data/lda_epidural_clean.rds")
lda_epidural_assess_clean <- readRDS("1 - data/3 - processed_data/lda_epidural_assess_clean.rds")

cat("✅ All cleaned data loaded\n")
cat("Episodes:", nrow(pca_episodes), "\n")
cat("Episode-shift rows:", nrow(episode_shift_table), "\n")
cat("Episode metrics:", nrow(episode_metrics), "\n")
cat("LDA urethral catheters:", nrow(lda_urethral_clean), "\n")
cat("LDA epidural catheters:", nrow(lda_epidural_clean), "\n")
cat("LDA epidural assessments:", nrow(lda_epidural_assess_clean), "\n")
