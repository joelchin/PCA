# ============================================================
# 1 — Load raw EPIC extracts
# Paediatric PCA/NCA study — Addenbrooke's Hospital (CUH NHS FT)
# Author: J Chin
# ============================================================
#
# Purpose
#   Reads every raw EPIC export (July 2026 refresh) into the R session.
#   Nothing is cleaned or saved here; cleaning is done in 2_cleaningscript.R.
#
# Inputs
#   Excel files named "Chin Paeds PCA Datasets - <category> - refresh
#   2026-07-15.xlsx" in the project root. To load a newer refresh, change
#   `file_suffix` below (this is the only place the refresh date is set).
#
# Objects created
#   admissions, anaesthesia_events, anaesthesia_ax_sdes, intraop_mar,
#   nerve_blocks, non_pca_mar, operations, pca_mar, pca_flowsheets,
#   problem_list, lda_urethral, lda_epidural, lda_epidural_assess,
#   pain_avpu (combined), pain_scores, avpu_scores
#
# Format notes for this refresh
#   - LDAs is one workbook with three sheets (urethral catheters, epidural
#     catheters, epidural catheter assessments), read as three objects.
#   - Pain/AVPU is two workbooks ("pt 1", "pt 2"), stacked and then split by
#     score_type into pain scores and ACVPU scores.
#   - "Anaesthesia Ax SDEs" is a category that first appears in the July 2026
#     refresh.
#
# Working directory
#   Run from the project root (the folder that contains "1 - data/"), or set the
#   environment variable PCA_PROJECT_DIR (e.g. in ~/.Renviron) to that folder.
#   No path is hard-coded in the scripts.
# ============================================================

# ---- Project root ----------------------------------------------------------
# Run from the project root (the folder containing "1 - data/"), or set the
# environment variable PCA_PROJECT_DIR (e.g. in ~/.Renviron). No path is hard-coded.
proj_dir <- Sys.getenv("PCA_PROJECT_DIR", unset = "")
if (nzchar(proj_dir)) setwd(proj_dir)
if (!dir.exists("1 - data")) {
  stop("Project root not found. Set the working directory (or PCA_PROJECT_DIR) ",
       "to the folder that contains '1 - data/'.")
}

# ---- Packages (installed once, outside the pipeline) ----------------------
required_pkgs <- c("tidyverse", "lubridate", "readxl", "furrr", "future",
                   "tictoc", "MatchIt", "logistf", "gtsummary", "flextable",
                   "officer", "labelled")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install the missing packages first: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}

library(tidyverse)
library(lubridate)
library(readxl)

file_prefix <- "Chin Paeds PCA Datasets - "
file_suffix <- " - refresh 2026-07-15.xlsx"

cat("============================================================\n")
cat("LOADING RAW DATA\n")
cat("============================================================\n\n")

# --- Core extracts ---
admissions          <- read_excel(paste0(file_prefix, "Admissions", file_suffix))
anaesthesia_events  <- read_excel(paste0(file_prefix, "Anaesthesia Events", file_suffix))
anaesthesia_ax_sdes <- read_excel(paste0(file_prefix, "Anaesthesia Ax SDEs", file_suffix))
intraop_mar         <- read_excel(paste0(file_prefix, "Intraoperative MAR", file_suffix))
nerve_blocks        <- read_excel(paste0(file_prefix, "Nerve Blocks", file_suffix))
non_pca_mar         <- read_excel(paste0(file_prefix, "Non PCA MAR", file_suffix))
operations          <- read_excel(paste0(file_prefix, "Operations", file_suffix))
pca_mar             <- read_excel(paste0(file_prefix, "PCA MAR", file_suffix))
pca_flowsheets      <- read_excel(paste0(file_prefix, "PCA Flowsheets", file_suffix))
problem_list        <- read_excel(paste0(file_prefix, "Problem List", file_suffix))

# --- LDA workbook: three sheets, one file ---
lda_file <- paste0(file_prefix, "LDAs", file_suffix)
lda_urethral        <- read_excel(lda_file, sheet = "Urethral Catheters")
lda_epidural        <- read_excel(lda_file, sheet = "Epidural Catheters")
lda_epidural_assess <- read_excel(lda_file, sheet = "Epidural Catheter Assessments")

# --- Pain AVPU: exported as two parts, combined then split by score type ---
pain_avpu_pt1 <- read_excel(paste0(file_prefix, "Pain AVPU - pt 1", file_suffix))
pain_avpu_pt2 <- read_excel(paste0(file_prefix, "Pain AVPU - pt 2", file_suffix))

pain_avpu <- bind_rows(pain_avpu_pt1, pain_avpu_pt2) |>
  rename(
    score_type  = DISP_NAME,
    score_time  = RECORDED_TIME,
    score_value = MEAS_VALUE
  )

pain_scores <- pain_avpu |>
  filter(score_type == "Pain Score") |>
  mutate(score_value = as.numeric(score_value))

avpu_scores <- pain_avpu |>
  filter(score_type == "NEWS2/MEOWS: ACVPU")

# --- Load summary ---
cat("============================================================\n")
cat("LOAD SUMMARY\n")
cat("============================================================\n")
cat("Admissions:", nrow(admissions), "rows\n")
cat("PCA MAR:", nrow(pca_mar), "rows\n")
cat("Non-PCA MAR:", nrow(non_pca_mar), "rows\n")
cat("Intraop MAR:", nrow(intraop_mar), "rows\n")
cat("Nerve blocks:", nrow(nerve_blocks), "rows\n")
cat("Operations:", nrow(operations), "rows\n")
cat("Anaesthesia events:", nrow(anaesthesia_events), "rows\n")
cat("Anaesthesia Ax SDEs:", nrow(anaesthesia_ax_sdes), "rows\n")
cat("LDA urethral catheters:", nrow(lda_urethral), "rows\n")
cat("LDA epidural catheters:", nrow(lda_epidural), "rows\n")
cat("LDA epidural assessments:", nrow(lda_epidural_assess), "rows\n")
cat("PCA flowsheets:", nrow(pca_flowsheets), "rows\n")
cat("Problem list:", nrow(problem_list), "rows\n")
cat("Pain scores:", nrow(pain_scores), "rows\n")
cat("AVPU scores:", nrow(avpu_scores), "rows\n")
