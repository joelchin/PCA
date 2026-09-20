# ============================================================
# POST-PIPELINE VERIFICATION CHECKLIST
#
# A set of smoke tests to run after the full pipeline
# (2_cleaningscript.R through 16_publication_tables.R) completes,
# checking that key derived fields and headline figures look as
# expected before trusting output for a manuscript table. The main
# LOAD CLEAN DATA script (the pipeline's actual step 5) follows below
# this checklist.
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. BASIC SANITY — did the new fields actually get created?
# ------------------------------------------------------------

stopifnot(exists("pca_episodes_flagged"), exists("operations_clean"))

cat("primary_opcs_substantive present:",
    "primary_opcs_substantive" %in% names(operations_clean), "\n")
cat("proc_category present on episodes:",
    "proc_category" %in% names(pca_episodes_flagged), "\n")
cat("exclude_pca_before_procedure present:",
    "exclude_pca_before_procedure" %in% names(pca_episodes_flagged), "\n\n")

# ------------------------------------------------------------
# 2. PSM POPULATION — should match the current locked figure.
# ------------------------------------------------------------

cat("============================================================\n")
cat("CHECK: PSM matched pairs (expect 237)\n")
cat("============================================================\n")
if (exists("psm_matched_keys")) {
  cat("Matched pairs:", nrow(psm_matched_keys) / 2, "\n")
} else {
  cat("psm_matched_keys not loaded — check after running 10_comparative_PSM.R\n")
}
cat("If this has moved noticeably from 237, investigate before\n")
cat("trusting any other number — that population is supposed to be\n")
cat("robust to everything fixed during validation.\n\n")

# ------------------------------------------------------------
# 3. proc_category DISTRIBUTION — should show the predicted shape:
# no more invented categories, DDH/foot cluster folded into
# Orthopaedic/spine, amputation/mass-excision cluster folded into
# Skin/trauma/general, Unclassified narrower than before.
# ------------------------------------------------------------

cat("============================================================\n")
cat("CHECK: proc_category distribution (S4)\n")
cat("============================================================\n")
print(table(pca_episodes_flagged$proc_category, useNA = "ifany"))
cat("\nExpect: only the original 16 category names + NA. No\n")
cat("'(surgical)', '(major)', or 'Maxillofacial' names should appear —\n")
cat("if they do, an old cached operations_clean.rds is being read\n")
cat("instead of the freshly rebuilt one.\n\n")

# ------------------------------------------------------------
# 4. all_opcs_excluded — should have moved (fixed bug), compare
# scale against pre-fix count
# ------------------------------------------------------------

cat("============================================================\n")
cat("CHECK: all_opcs_excluded / exclude_pca_before_procedure counts\n")
cat("============================================================\n")
cat("all_opcs_excluded (S4):",
    sum(pca_episodes_flagged$all_opcs_excluded, na.rm = TRUE), "\n")
cat("exclude_pca_before_procedure (S4):",
    sum(pca_episodes_flagged$exclude_pca_before_procedure, na.rm = TRUE), "\n\n")

# ------------------------------------------------------------
# 5. N13.4 SPOT CHECK — the bundled N13.4/T30.9 row handled in
# 2_cleaningscript.R (Section 6A) should correctly show a
# T30.9-driven category, not N13.4/Unclassified
# ------------------------------------------------------------

cat("============================================================\n")
cat("CHECK: N13.4 + T30.9 bundled row — should now resolve via T30.9\n")
cat("============================================================\n")
n13_check <- operations_clean |>
  filter(str_detect(coalesce(PFMD_PROC_OPCS_CODES, ""), "N13.4")) |>
  select(PFMD_PROC_NAMES, PFMD_PROC_OPCS_CODES, primary_opcs,
         primary_opcs_substantive, proc_category)
print(n13_check, width = Inf)
cat("\nThe 'EXPLORATORY LAPAROTOMY' row should show\n")
cat("primary_opcs_substantive = T30.9 and proc_category =\n")
cat("Skin/trauma/general — NOT N13.4/Unclassified.\n\n")

# ------------------------------------------------------------
# 6. TABLE 2 — antiemetic rate, composite OR — confirm these
# haven't silently moved from the locked figures without you knowing
# ------------------------------------------------------------

cat("============================================================\n")
cat("CHECK: locked headline numbers, confirm still consistent\n")
cat("============================================================\n")
cat("Reactive antiemetic rate (S4 oxycodone, expect ~34.1%):\n")
oxy_s4 <- pca_episodes_flagged |> filter(drug == "oxycodone")
cat(round(100 * mean(oxy_s4$reactive_antiemetic_final, na.rm = TRUE), 1), "%\n")
cat("(N =", nrow(oxy_s4), ")\n\n")

cat("---\n")
cat("If ANYTHING above looks unexpected, stop before re-running Table\n")
cat("2/3 exports — better to catch a problem here than in a published\n")
cat("table.\n")

# ============================================================
# DRUG LISTS
# ============================================================
source("drug_lists.R")

# ============================================================
# LOAD CLEAN DATA — Paediatric PCA Study
# Run this at the start of every session
# ============================================================
library(tidyverse)
library(lubridate)
library(furrr)
library(future)

setwd("//net.addenbrookes.nhs.uk/root/Users2-4/chinj/pcaraw")

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
# pca_episodes is kept as an alias since older scripts still
# reference that name — both point to the same object.
pca_episodes_flagged      <- readRDS("1 - data/3 - processed_data/pca_episodes_flagged.rds")
pca_episodes              <- pca_episodes_flagged

lda_urethral_clean        <- readRDS("1 - data/3 - processed_data/lda_urethral_clean.rds")
lda_epidural_clean        <- readRDS("1 - data/3 - processed_data/lda_epidural_clean.rds")
lda_epidural_assess_clean <- readRDS("1 - data/3 - processed_data/lda_epidural_assess_clean.rds")

# pca_flowsheets_clean.rds is not produced by any script in this
# pipeline (1_load_raw_data.R loads the raw pca_flowsheets table, but
# nothing downstream cleans/saves it under this name — 8_true_total_
# dose.R reads the raw flowsheet Excel file directly instead). Loaded
# here as optional, since nothing else in the pipeline depends on it.
pca_flowsheets_clean <- tryCatch(
  readRDS("1 - data/3 - processed_data/pca_flowsheets_clean.rds"),
  error = function(e) {
    cat("pca_flowsheets_clean.rds not found — not produced by this\n")
    cat("pipeline; skipping (nothing downstream depends on it).\n")
    NULL
  }
)

cat("All cleaned data loaded\n")
cat("Episodes:", nrow(pca_episodes), "\n")
cat("Episode-shift rows:", nrow(episode_shift_table), "\n")
cat("Episode metrics:", nrow(episode_metrics), "\n")
cat("LDA urethral catheters:", nrow(lda_urethral_clean), "\n")
cat("LDA epidural catheters:", nrow(lda_epidural_clean), "\n")
cat("LDA epidural assessments:", nrow(lda_epidural_assess_clean), "\n")
if (!is.null(pca_flowsheets_clean)) {
  cat("PCA flowsheets:", nrow(pca_flowsheets_clean), "\n")
}