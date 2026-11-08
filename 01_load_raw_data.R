# ============================================================
# 01 - Load raw data
# ============================================================
# Reads the raw EPIC extract (one Excel file per data category plus a
# three-sheet workbook of lines/drains/airways, LDA) and sets up the
# session for 02_clean_data.R.
#
#   - LDAs arrive as one workbook with three sheets (urethral
#     catheters, epidural catheters, epidural catheter assessments).
#   - Pain and ACVPU observations arrive as one feed in two parts,
#     distinguished by a score type column; they are combined and then
#     split into pain_scores and avpu_scores.
#
# Only missing packages are installed (binaries only, no Rtools), and
# only if none of the pipeline packages is already loaded, which avoids
# the "cannot be unloaded" error. Input paths and files are checked
# before anything is read. Run from a fresh R session.
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
# Set PCA_DATA_DIR in your user-level .Renviron (see README). The folder
# holds the raw extract and all derived data and must live outside this
# repository.
data_dir <- Sys.getenv("PCA_DATA_DIR")

if (!nzchar(data_dir)) {
  stop("PCA_DATA_DIR is not set. Add it to your .Renviron (see README) ",
       "and restart R.")
}
if (!dir.exists(data_dir)) {
  stop("Data folder not found: ", data_dir,
       "\nCheck that you are connected to the network drive.")
}
# Remember where the code lives so drug_lists.R can be found after setwd()
options(pca_code_dir = normalizePath(getwd()))
setwd(data_dir)
cat("Working directory:", getwd(), "\n\n")

# ------------------------------------------------------------
# 2. PACKAGES
# ------------------------------------------------------------
pipeline_pkgs <- c("tidyverse", "lubridate", "readxl",
                   "flextable", "officer", "gtsummary", "labelled",
                   "furrr", "future", "tictoc",
                   "MatchIt", "logistf", "zip")

# Check for installation without loading (requireNamespace() would load it)
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
  
  # Set a CRAN mirror if none is set
  repos <- getOption("repos")
  if (is.null(repos) || is.na(repos["CRAN"]) || isTRUE(unname(repos["CRAN"]) == "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  # Binaries only
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
# 3. INPUT FILES: check everything exists before reading anything
# ------------------------------------------------------------
# File names are expected as <prefix><table name><suffix>; set the prefix
# and suffix in .Renviron (see README).
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

# Core extracts
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

# LDA workbook (three sheets)
lda_file <- paste0(file_prefix, "LDAs", file_suffix)
lda_urethral        <- read_excel(lda_file, sheet = "Urethral Catheters")
lda_epidural        <- read_excel(lda_file, sheet = "Epidural Catheters")
lda_epidural_assess <- read_excel(lda_file, sheet = "Epidural Catheter Assessments")

# Pain/ACVPU feed (two parts), combined then split by score type
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