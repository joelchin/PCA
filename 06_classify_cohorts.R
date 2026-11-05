# 06 — GLOBAL POPULATION CLASSIFICATION
# Surgical (major OPCS) vs Non-surgical (minor/no OPCS)
# ============================================================
library(tidyverse)
library(lubridate)

# ---- Manual correction (site-specific) -------------------------------------
# A patient with GENDER "Not specified" in the extract had sex confirmed as
# female from the clinical record. The correction is applied to ALL of that
# patient's episodes. The patient ID is NOT stored in this repository: it is
# read from a local text file (one PAT_ID per line) kept alongside the data.
# If the file is absent, no correction is applied.
gender_fix_file <- file.path(Sys.getenv("PCA_DATA_DIR"), "manual_corrections",
                             "gender_fix_pat_ids.txt")
gender_fix_pat <- if (file.exists(gender_fix_file)) {
  trimws(readLines(gender_fix_file, warn = FALSE))
} else character(0)
gender_fix_pat <- gender_fix_pat[nzchar(gender_fix_pat)]

if (length(gender_fix_pat) > 0) {
  fix_csn <- admissions_clean |>
    filter(PAT_ID %in% gender_fix_pat) |>
    pull(PAT_ENC_CSN_ID)

  cat("Encounters corrected:", length(fix_csn), "\n")
  stopifnot(length(fix_csn) >= 1)                    # stops if an ID was mistyped

  pca_episodes_flagged <- pca_episodes_flagged |>
    mutate(GENDER = if_else(PAT_ENC_CSN_ID %in% fix_csn & GENDER == "Not specified",
                            "Female", GENDER))
} else {
  cat("No manual gender correction file found; no correction applied.\n")
}
# -----------------------------------------------------------------------------

cat("============================================================\n")
cat("ANALYSIS 4 — GLOBAL CLASSIFICATION\n")
cat("============================================================\n\n")
# Require cleaned data
stopifnot(exists("pca_episodes_flagged"))
stopifnot(exists("admissions_clean"))
stopifnot(exists("operations_clean"))
stopifnot("case_type_refined" %in% names(pca_episodes_flagged))
# admissions_clean/operations_clean are used here rather than the raw
# admissions/operations objects: operations_clean is a strict superset
# of operations (built via mutate()) and still carries
# PFMD_PROC_OPCS_CODES, the column used below, while the raw objects
# only exist transiently during initial loading and aren't reloaded
# in later pipeline stages.
#
# S3 and S4 both classify surgical/medical status consistently via
# case_type_refined (built in 04_derive_episode_flags.R), rather than each population
# tier maintaining its own independent classification logic — this
# avoids the surgical/medical split silently meaning something
# different depending on which population tier is being described.
# S2 (first oxycodone episode, no concurrent epidural) is not part of
# this pipeline's active population tiers.
# ============================================================
# S4 GLOBAL CLASSIFICATION (ALL OXYCODONE EPISODES)
# ============================================================
cat("Creating S4 global classification (all oxycodone)...\n\n")

all_oxy <- pca_episodes_flagged |>
  filter(
    drug == "oxycodone",
    !is.na(episode_end)   # excludes the 23 episodes with no end time: duration and
    # shift-based outcomes cannot be computed for them, so S4
    # stays aligned with episode_metrics (episodes with a
    # partial missing stop but an observed end are kept)
  )

cat("S4 excluded for missing end time:",
    sum(pca_episodes_flagged$drug == "oxycodone" & is.na(pca_episodes_flagged$episode_end)),
    "\n")
cat("S4 total (complete episodes only):", nrow(all_oxy), "\n\n")
all_oxy_classified <- all_oxy |>
  mutate(
    case_type = case_when(
      case_type_refined == "surgical" ~ "Surgical",
      case_type_refined == "medical"  ~ "Medical",
      TRUE                             ~ NA_character_
    )
  )
s4_medical_keys <- all_oxy_classified |>
  filter(case_type == "Medical") |>
  select(PAT_ENC_CSN_ID, episode_id)
s4_surgical_keys <- all_oxy_classified |>
  filter(case_type == "Surgical") |>
  select(PAT_ENC_CSN_ID, episode_id)
# ============================================================
# S3 — PSM POPULATION (morphine vs oxycodone comparative)
#
# Population-tier definitions live here in the classification layer,
# consistent with this project's separation of concerns between
# classification, flagging, and comparative analysis.
# ============================================================
cat("Creating S3 (PSM population: morphine vs oxycodone)...\n\n")

# S3's surgical filter uses case_type_refined rather than the coarser
# case_type ("does any OPCS-coded operation exist for this
# admission" — a blunt presence check). case_type_refined
# additionally excludes episodes where the operation present is
# genuinely non-explanatory for PCA need (entirely minor procedure
# codes, no operation actually precedes the episode, or the episode
# substantially predates the operation's anaesthesia timing by more
# than a 2h buffer — see 04_derive_episode_flags.R for the full derivation). These
# cases are not treated as surgical indications for PCA.

stopifnot("case_type_refined" %in% names(pca_episodes_flagged))

# ============================================================
# mode %in% c("PCA", "NCA") filter
#
# Without this filter, episodes with mode == "unknown" (a genuine
# data-quality gap — delivery mode couldn't be determined, not a
# meaningful clinical category) could enter the matched population
# and be treated by MatchIt as a legitimate third group to balance
# on. Table 2 already excludes these episodes ("Episodes with
# undetermined PCA/NCA mode excluded from this table"), so S3 is
# aligned with that here.
#
# Medical oncology (ADMITTING_SPECIALTY "PAEDIATRIC MEDICAL ONCOLOGY" /
# "MEDICAL ONCOLOGY") is excluded from S3 here, so every script that
# reads the S3 keys sees the same population. coalesce() keeps episodes
# with a missing specialty in, rather than letting NA silently drop them.
# ============================================================

s3_morphine_keys <- pca_episodes_flagged |>
  filter(
    drug == "morphine",
    is_first_line == TRUE,
    is_first_episode == TRUE,
    concurrent_epidural == FALSE,
    case_type_refined == "surgical",
    mode %in% c("PCA", "NCA"),
    !is.na(episode_end),
    exclude_minor_procedure_any == FALSE,
    !coalesce(str_detect(toupper(as.character(ADMITTING_SPECIALTY)), "MEDICAL ONCOLOGY"), FALSE)
  ) |>
  select(PAT_ENC_CSN_ID, episode_id)
s3_oxycodone_keys <- pca_episodes_flagged |>
  filter(
    drug == "oxycodone",
    is_first_line == TRUE,
    is_first_episode == TRUE,
    concurrent_epidural == FALSE,
    case_type_refined == "surgical",
    mode %in% c("PCA", "NCA"),
    !is.na(episode_end),
    exclude_minor_procedure_any == FALSE,
    !coalesce(str_detect(toupper(as.character(ADMITTING_SPECIALTY)), "MEDICAL ONCOLOGY"), FALSE)
  ) |>
  select(PAT_ENC_CSN_ID, episode_id)

# ============================================================
# FIRST EXPOSURE, FIRST EPISODE, FIRST LINE (patient level)  <-- NEW
#
# is_first_line / is_first_episode are judged WITHIN an admission, so a
# child with an earlier admission could re-enter as "first line" again.
# The study episode must now be the child's FIRST-EVER recorded PCA/NCA
# episode (of any opioid) at this centre. That makes it their first exposure
# to the study opioid, their first episode, and first line; an episode that
# follows an earlier PCA/NCA episode in any admission is excluded. A child's
# first-ever episode has only one drug, so no child can be in both arms.
# ============================================================
fe_pid <- admissions_clean |>
  select(PAT_ENC_CSN_ID, PAT_ID) |>
  distinct(PAT_ENC_CSN_ID, .keep_all = TRUE)

first_exposure_keys <- pca_episodes_flagged |>
  left_join(fe_pid, by = "PAT_ENC_CSN_ID") |>
  mutate(patient_key = coalesce(as.character(PAT_ID), paste0("csn_", PAT_ENC_CSN_ID))) |>
  arrange(patient_key, episode_start, PAT_ENC_CSN_ID, episode_id) |>
  group_by(patient_key) |>
  filter(row_number() == 1) |>
  ungroup() |>
  select(PAT_ENC_CSN_ID, episode_id)

s3_morphine_keys <- s3_morphine_keys |>
  semi_join(first_exposure_keys, by = c("PAT_ENC_CSN_ID", "episode_id"))
s3_oxycodone_keys <- s3_oxycodone_keys |>
  semi_join(first_exposure_keys, by = c("PAT_ENC_CSN_ID", "episode_id"))

cat("S3 pools after first-exposure restriction — morphine:", nrow(s3_morphine_keys),
    "| oxycodone:", nrow(s3_oxycodone_keys), "\n\n")

# ============================================================
# CROSS-ARM PATIENT EXCLUSION
#
# Found via a repeat-patient/independence check on the matched PSM
# output: 5 real children appeared in BOTH drug arms (morphine on one
# admission, oxycodone on a different, later admission). The existing
# dedup logic deduplicates EACH drug's pool separately (one episode
# per patient WITHIN morphine, one per patient WITHIN oxycodone) —
# nothing prevented the same real patient from legitimately passing
# both drugs' individual checks and being matched twice, once against
# a genuinely different morphine patient and once against a
# genuinely different oxycodone patient. Violates the independence
# assumption standard matching relies on.
#
# Direct decision: exclude these patients from S3 entirely (option 1
# of three considered — cluster-robust SEs or reporting-as-limitation
# were the alternatives, not used). Applied here, upstream of
# matching, so the propensity model itself never sees these
# ambiguous cross-arm patients.
# ============================================================

stopifnot(exists("admissions_clean"))

morphine_patients <- s3_morphine_keys |>
  left_join(admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID), by = "PAT_ENC_CSN_ID") |>
  distinct(PAT_ID)
oxycodone_patients <- s3_oxycodone_keys |>
  left_join(admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID), by = "PAT_ENC_CSN_ID") |>
  distinct(PAT_ID)

cross_arm_patients <- intersect(morphine_patients$PAT_ID, oxycodone_patients$PAT_ID)

cat("Cross-arm patients found (appear in BOTH morphine and oxycodone\n")
cat("S3 pools, via different admissions):", length(cross_arm_patients), "\n")

if (length(cross_arm_patients) > 0) {
  s3_morphine_keys <- s3_morphine_keys |>
    left_join(admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID), by = "PAT_ENC_CSN_ID") |>
    filter(!(PAT_ID %in% cross_arm_patients)) |>
    select(PAT_ENC_CSN_ID, episode_id)
  
  s3_oxycodone_keys <- s3_oxycodone_keys |>
    left_join(admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID), by = "PAT_ENC_CSN_ID") |>
    filter(!(PAT_ID %in% cross_arm_patients)) |>
    select(PAT_ENC_CSN_ID, episode_id)
  
  cat("S3 pools after cross-arm exclusion — morphine:", nrow(s3_morphine_keys),
      "| oxycodone:", nrow(s3_oxycodone_keys), "\n")
}

# ============================================================
# SAVE ALL KEYS
# ============================================================
cat("Saving classification keys...\n")
saveRDS(s4_medical_keys, "1 - data/3 - processed_data/s4_medical_keys.rds")
saveRDS(s4_surgical_keys, "1 - data/3 - processed_data/s4_surgical_keys.rds")
saveRDS(s3_morphine_keys, "1 - data/3 - processed_data/s3_morphine_keys.rds")
saveRDS(s3_oxycodone_keys, "1 - data/3 - processed_data/s3_oxycodone_keys.rds")
cat("s4_medical_keys.rds (all oxy, medical)\n")
cat("s4_surgical_keys.rds (all oxy, surgical)\n")
cat("s3_morphine_keys.rds (PSM population, morphine — NOW using case_type_refined)\n")
cat("s3_oxycodone_keys.rds (PSM population, oxycodone — NOW using case_type_refined)\n\n")
# ============================================================
# SUMMARY
# ============================================================
cat("============================================================\n")
cat("CLASSIFICATION SUMMARY\n")
cat("============================================================\n\n")
cat("S4 (All oxycodone):\n")
cat(sprintf(" Medical (no/minor OPCS): %d (%.1f%%)\n",
            nrow(s4_medical_keys), 100*nrow(s4_medical_keys)/nrow(all_oxy)))
cat(sprintf(" Surgical (major OPCS): %d (%.1f%%)\n\n",
            nrow(s4_surgical_keys), 100*nrow(s4_surgical_keys)/nrow(all_oxy)))
cat("S3 (PSM population, case_type_refined):\n")
cat(sprintf(" Morphine: n = %d\n", nrow(s3_morphine_keys)))
cat(sprintf(" Oxycodone: n = %d\n\n", nrow(s3_oxycodone_keys)))
cat("ACTION NEEDED: side_effect_temporal_profile.R still defines\n")
cat("its own inline copy of s3_morphine_keys/s3_oxycodone_keys — check\n")
cat("whether it also needs the case_type -> case_type_refined switch\n")
cat("applied, or it will silently diverge from this version.\n\n")
cat("S3/S4 consistently classify\n")
cat("surgical vs medical using case_type_refined — the opcs_data/\n")
cat("minor_procedure_opcs path (admission-level, first-operation-only,\n")
cat("hardcoded duplicate minor-procedure list) is not used\n")
cat("anywhere in this script. S2 is not part of this pipeline (confirmed\n")
cat("dead code, see note near top of script).\n\n")
cat("============================================================\n")