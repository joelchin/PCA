# ============================================================
# 05 - Load processed data
# ============================================================
# Loads the processed tables saved by scripts 02-04 and the shared drug
# lists. Run at the start of an analysis session, in a fresh R session,
# before any of scripts 06 onwards.
# ============================================================

library(tidyverse)
library(lubridate)

code_dir <- normalizePath(getwd())
options(pca_code_dir = code_dir)
source(file.path(code_dir, "drug_lists.R"))

data_dir <- Sys.getenv("PCA_DATA_DIR")
if (!nzchar(data_dir)) stop("PCA_DATA_DIR is not set (see README).")
setwd(data_dir)

processed <- "1 - data/3 - processed_data/"

admissions_clean          <- readRDS(paste0(processed, "admissions_clean.rds"))
pca_mar_clean             <- readRDS(paste0(processed, "pca_mar_clean.rds"))
non_pca_mar_clean         <- readRDS(paste0(processed, "non_pca_mar_clean.rds"))
intraop_mar_clean         <- readRDS(paste0(processed, "intraop_mar_clean.rds"))
nerve_blocks_clean        <- readRDS(paste0(processed, "nerve_blocks_clean.rds"))
operations_clean          <- readRDS(paste0(processed, "operations_clean.rds"))
anaesthesia_events_clean  <- readRDS(paste0(processed, "anaesthesia_events_clean.rds"))
anaesthesia_ax_clean      <- readRDS(paste0(processed, "anaesthesia_ax_clean.rds"))
pain_scores_clean         <- readRDS(paste0(processed, "pain_scores_clean.rds"))
shifts                    <- readRDS(paste0(processed, "1_shifts.rds"))
episode_shift_table       <- readRDS(paste0(processed, "episode_shift_table.rds"))
episode_metrics           <- readRDS(paste0(processed, "episode_metrics.rds"))
lda_urethral_clean        <- readRDS(paste0(processed, "lda_urethral_clean.rds"))
lda_epidural_clean        <- readRDS(paste0(processed, "lda_epidural_clean.rds"))
lda_epidural_assess_clean <- readRDS(paste0(processed, "lda_epidural_assess_clean.rds"))

# pca_episodes_flagged is the saved episode table; pca_episodes is an
# alias because several scripts use that name
pca_episodes_flagged      <- readRDS(paste0(processed, "pca_episodes_flagged.rds"))
pca_episodes              <- pca_episodes_flagged

cat("All processed data loaded\n")
cat("Episodes:", nrow(pca_episodes), "\n")
cat("Episode-shift rows:", nrow(episode_shift_table), "\n")
cat("Episode metrics:", nrow(episode_metrics), "\n")
cat("LDA urethral catheters:", nrow(lda_urethral_clean), "\n")
cat("LDA epidural catheters:", nrow(lda_epidural_clean), "\n")
cat("LDA epidural assessments:", nrow(lda_epidural_assess_clean), "\n")
