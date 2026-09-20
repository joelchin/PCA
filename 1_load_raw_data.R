# ============================================================
# 1 — LOAD RAW DATA
#
# Reads the raw EPIC extract (one Excel file per data category, plus
# a three-sheet workbook for line/drain/airway - LDA - data) and does
# the minimal tidying needed before cleaning proper starts in
# 2_cleaningscript.R.
#
# Two aspects of the extract's structure are handled explicitly here:
#   - LDAs arrive as a single workbook with three sheets (urethral
#     catheters, epidural catheters, epidural catheter assessments),
#     read separately by sheet name rather than as one flat table.
#   - Pain AVPU is one combined feed containing both pain scores and
#     NEWS2/MEOWS ACVPU readings, distinguished by a score_type
#     column. It arrives in two parts, which are combined and then
#     split into pain_scores and avpu_scores below.
# ============================================================

setwd("//net.addenbrookes.nhs.uk/root/Users2-4/chinj/pcaraw")

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
