# ============================================================
# 01 — LOAD RAW DATA
#
# Reads the raw EPIC extract (one Excel file per data category, plus
# a three-sheet workbook for line/drain/airway - LDA - data) and does
# the minimal tidying needed before cleaning proper starts in
# 02_clean_data.R.
#
# Two aspects of the extract's structure are handled explicitly here:
#   - LDAs arrive as a single workbook with three sheets (urethral
#     catheters, epidural catheters, epidural catheter assessments),
#     read separately by sheet name rather than as one flat table.
#   - Pain AVPU is one combined feed containing both pain scores and
#     NEWS2/MEOWS ACVPU readings, distinguished by a score_type
#     column. It arrives in two parts, which are combined and then
#     split into pain_scores and avpu_scores below.
#
# Housekeeping built into this script (Oct 2026):
#   - Packages: only packages that are MISSING are installed, and only
#     when none of the pipeline packages is already attached in the
#     session. Re-running this script never reinstalls anything, and it
#     can no longer trigger the "cannot be unloaded" error (installing
#     over a package, such as officer, that is already loaded).
#   - Binary packages only: no Rtools is needed or used.
#   - The working directory and every input file are checked up front,
#     so a wrong path or a missing file stops the script with a clear
#     message instead of failing part-way through.
#   - A warning is printed if the session already contains objects
#     (a sign of a restored or stale workspace). The pipeline should
#     always be run from a fresh R session.
# ============================================================

# ------------------------------------------------------------
# 0. FRESH-SESSION CHECK
# ------------------------------------------------------------
stale_objects <- ls(envir = globalenv())
if (length(stale_objects) > 0) {
  warning("The R session is not empty (", length(stale_objects), " objects, e.g. ",
          paste(head(stale_objects, 5), collapse = ", "), "). ",
          "For a clean run, restart R (Session > Restart R) with workspace restore OFF, ",
          "then run this script again.", call. = FALSE)
  cat("\n*** WARNING: the session already contains objects. Restart R for a clean run. ***\n\n")
}

# ------------------------------------------------------------
# 1. WORKING DIRECTORY
# ------------------------------------------------------------
# Set PCA_DATA_DIR in your user-level .Renviron (see README), e.g.
#   PCA_DATA_DIR=//your-server/share/path/to/project_data
# The folder holds the raw extract and all derived data. It must live
# OUTSIDE this repository and is never committed.
data_dir <- Sys.getenv("PCA_DATA_DIR")

if (!nzchar(data_dir)) {
  stop("PCA_DATA_DIR is not set. Add it to your .Renviron (see README) ",
       "and restart R.")
}
if (!dir.exists(data_dir)) {
  stop("Data folder not found: ", data_dir,
       "\nCheck that you are connected to the network drive.")
}
setwd(data_dir)
cat("Working directory:", getwd(), "\n\n")

# ------------------------------------------------------------
# 2. PACKAGES
# Install only what is missing, and only in a session where none of
# the pipeline packages is loaded yet.
# ------------------------------------------------------------
pipeline_pkgs <- c("tidyverse", "lubridate", "readxl",
                   "flextable", "officer", "gtsummary", "labelled",
                   "furrr", "future", "tictoc",
                   "MatchIt", "logistf", "zip")

# Check "installed?" WITHOUT loading anything (requireNamespace() would load
# the package, which is exactly what we must avoid before an install).
is_installed <- function(pkgs) vapply(pkgs, function(p) nzchar(system.file(package = p)),
                                      logical(1))

missing_pkgs <- pipeline_pkgs[!is_installed(pipeline_pkgs)]

if (length(missing_pkgs) > 0) {
  attached_now <- pipeline_pkgs[paste0("package:", pipeline_pkgs) %in% search()]
  
  if (length(attached_now) > 0) {
    stop("These packages are missing: ", paste(missing_pkgs, collapse = ", "), ".\n",
         "Other pipeline packages (", paste(attached_now, collapse = ", "),
         ") are already loaded in this session, so installing now could fail with ",
         "'cannot be unloaded'.\n",
         "Restart R (Session > Restart R), then run this script again.")
  }
  
  cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  
  # CRAN mirror, in case none is set
  repos <- getOption("repos")
  if (is.null(repos) || is.na(repos["CRAN"]) || isTRUE(unname(repos["CRAN"]) == "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  # Windows binaries only: never try to build from source (no Rtools needed)
  options(install.packages.compile.from.source = "never")
  
  install.packages(missing_pkgs, type = "binary")
  
  still_missing <- missing_pkgs[!is_installed(missing_pkgs)]
  if (length(still_missing) > 0) {
    stop("Could not install: ", paste(still_missing, collapse = ", "), ".\n",
         "Check the messages above (network access, CRAN mirror, library permissions). ",
         "If a package has no binary for your R version yet, install it by hand ",
         "(or ask whoever supports your R installation) and run this script again.")
  }
  cat("Installation complete.\n\n")
}

cat("All pipeline packages are installed.\n")

library(tidyverse)
library(lubridate)
library(readxl)

cat("R version:", R.version.string, "\n")
cat("Key package versions: ",
    paste0(c("MatchIt", "logistf", "gtsummary", "flextable", "officer"), " ",
           vapply(c("MatchIt", "logistf", "gtsummary", "flextable", "officer"),
                  function(p) as.character(packageVersion(p)), character(1)),
           collapse = " | "), "\n\n", sep = "")

# ------------------------------------------------------------
# 3. INPUT FILES — check everything exists before reading anything
# ------------------------------------------------------------
# Extract file names are site-specific; set them in .Renviron (see README).
# Files are expected as: <prefix><table name><suffix>
file_prefix <- Sys.getenv("PCA_FILE_PREFIX", unset = "PCA Datasets - ")
file_suffix <- Sys.getenv("PCA_FILE_SUFFIX", unset = " - refresh.xlsx")

input_names <- c("Admissions", "Anaesthesia Events", "Anaesthesia Ax SDEs",
                 "Intraoperative MAR", "Nerve Blocks", "Non PCA MAR",
                 "Operations", "PCA MAR", "PCA Flowsheets", "Problem List",
                 "LDAs", "Pain AVPU - pt 1", "Pain AVPU - pt 2")
input_files <- paste0(file_prefix, input_names, file_suffix)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Input file(s) not found in ", getwd(), ":\n  ",
       paste(missing_files, collapse = "\n  "),
       "\nCheck the file names and the refresh date in 'file_suffix'.")
}

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