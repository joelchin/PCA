# ============================================================
# 06 - Classify cohorts
# ============================================================
# Builds the episode keys for the two analysis populations:
#
#   descriptive oxycodone population (s4_*_keys): all oxycodone episodes
#     with a recorded end time, split into surgical and medical using
#     case_type_refined (from 04_derive_episode_flags.R)
#   primary comparative population (s3_*_keys): each child's first-ever
#     recorded PCA/NCA episode, one episode per child, without
#     concurrent epidural analgesia, surgical, with the predefined
#     exclusions applied, a determined delivery mode and a recorded end
#     time, excluding medical oncology
#
# Keys are saved as .rds files in 1 - data/3 - processed_data/.
# ============================================================

library(tidyverse)
library(lubridate)

stopifnot(
  exists("pca_episodes_flagged"),
  exists("admissions_clean"),
  "case_type_refined" %in% names(pca_episodes_flagged)
)

# ------------------------------------------------------------
# Manual sex correction (site-specific)
# ------------------------------------------------------------
# One patient recorded as "Not specified" had sex confirmed as female
# from the clinical record. LOCAL VERSION: put the PAT_ID here (exactly as
# stored in admissions_clean). Do not publish this file with the ID in it.
gender_fix_pat <- c("PUT_PAT_ID_HERE")

if (length(gender_fix_pat) > 0) {
  fix_csn <- admissions_clean |>
    filter(PAT_ID %in% gender_fix_pat) |>
    pull(PAT_ENC_CSN_ID)
  
  cat("Encounters corrected:", length(fix_csn), "\n")
  stopifnot(length(fix_csn) >= 1)   # stops if an ID was mistyped
  
  pca_episodes_flagged <- pca_episodes_flagged |>
    mutate(GENDER = if_else(PAT_ENC_CSN_ID %in% fix_csn & GENDER == "Not specified",
                            "Female", GENDER))
} else {
  cat("No manual sex correction file found; no correction applied.\n")
}

# ============================================================
# DESCRIPTIVE OXYCODONE POPULATION
# ============================================================
cat("\nBuilding the descriptive oxycodone population...\n")

# Episodes with no recorded end time are excluded: duration and
# shift-based outcomes cannot be computed for them
all_oxy <- pca_episodes_flagged |>
  filter(
    drug == "oxycodone",
    !is.na(episode_end)
  )

cat("Excluded for missing end time:",
    sum(pca_episodes_flagged$drug == "oxycodone" & is.na(pca_episodes_flagged$episode_end)),
    "\n")
cat("Oxycodone episodes with an end time:", nrow(all_oxy), "\n\n")

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
# PRIMARY COMPARATIVE POPULATION (morphine vs oxycodone)
# ============================================================
cat("Building the primary comparative population...\n\n")

# Surgical status uses case_type_refined rather than case_type: it also
# excludes episodes whose recorded operation cannot explain PCA use
# (entirely minor procedure codes, no preceding operation, or PCA
# starting more than 2 hours before anaesthesia stop).
#
# Episodes with an undetermined delivery mode (mode == "unknown") are
# excluded so that they are not treated by MatchIt as a third group.
#
# Medical oncology is excluded here so that every script reading these
# keys sees the same population. coalesce() keeps episodes with a
# missing specialty in.

s3_filter <- function(df, drug_name, require_end_time = TRUE) {
  df |>
    filter(
      drug == drug_name,
      is_first_line == TRUE,
      is_first_episode == TRUE,
      concurrent_epidural == FALSE,
      case_type_refined == "surgical",
      mode %in% c("PCA", "NCA"),
      !require_end_time | !is.na(episode_end),
      exclude_minor_procedure_any == FALSE,
      !coalesce(str_detect(toupper(as.character(ADMITTING_SPECIALTY)),
                           "MEDICAL ONCOLOGY"), FALSE)
    ) |>
    select(PAT_ENC_CSN_ID, episode_id)
}

s3_morphine_keys  <- s3_filter(pca_episodes_flagged, "morphine")
s3_oxycodone_keys <- s3_filter(pca_episodes_flagged, "oxycodone")

# ------------------------------------------------------------
# First-ever PCA/NCA episode (patient level)
# ------------------------------------------------------------
# is_first_line and is_first_episode are judged within an admission, so a
# child with an earlier admission could re-enter as first line. The study
# episode must be the child's first-ever recorded PCA/NCA episode (of any
# opioid) at this centre. An episode that follows an earlier PCA/NCA
# episode in any admission is excluded. A first-ever episode has only
# one drug, so no child can be in both arms.
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

cat("Pools after first-exposure restriction - morphine:", nrow(s3_morphine_keys),
    "| oxycodone:", nrow(s3_oxycodone_keys), "\n\n")

# ------------------------------------------------------------
# Counts for the population flow diagram (Figure 1)
# ------------------------------------------------------------
# Eligible episodes before the end-time requirement, and how many that
# requirement removes (all other criteria, including first exposure,
# applied).
eligible_no_end_mor <- s3_filter(pca_episodes_flagged, "morphine", require_end_time = FALSE) |>
  semi_join(first_exposure_keys, by = c("PAT_ENC_CSN_ID", "episode_id"))
eligible_no_end_oxy <- s3_filter(pca_episodes_flagged, "oxycodone", require_end_time = FALSE) |>
  semi_join(first_exposure_keys, by = c("PAT_ENC_CSN_ID", "episode_id"))

cat("--- Flow diagram counts ---\n")
cat("All episodes:", nrow(pca_episodes_flagged),
    "| patients:", n_distinct(coalesce(
      as.character(fe_pid$PAT_ID[match(pca_episodes_flagged$PAT_ENC_CSN_ID, fe_pid$PAT_ENC_CSN_ID)]),
      paste0("csn_", pca_episodes_flagged$PAT_ENC_CSN_ID))), "\n")
cat("Morphine or oxycodone episodes:",
    sum(pca_episodes_flagged$drug %in% c("morphine", "oxycodone")), "\n")
cat("Eligible before the end-time requirement - morphine:", nrow(eligible_no_end_mor),
    "| oxycodone:", nrow(eligible_no_end_oxy),
    "| total:", nrow(eligible_no_end_mor) + nrow(eligible_no_end_oxy), "\n")
cat("Excluded for no recorded end time - morphine:",
    nrow(eligible_no_end_mor) - nrow(s3_morphine_keys),
    "| oxycodone:", nrow(eligible_no_end_oxy) - nrow(s3_oxycodone_keys), "\n\n")

# ------------------------------------------------------------
# Check that no child is in both arms
# ------------------------------------------------------------
# Expected to be zero after the first-exposure restriction. Any child
# found in both pools is removed from both, so that matching never sees
# the same child twice.
pat_lookup <- admissions_clean |> select(PAT_ENC_CSN_ID, PAT_ID)

morphine_patients <- s3_morphine_keys |>
  left_join(pat_lookup, by = "PAT_ENC_CSN_ID") |>
  distinct(PAT_ID)
oxycodone_patients <- s3_oxycodone_keys |>
  left_join(pat_lookup, by = "PAT_ENC_CSN_ID") |>
  distinct(PAT_ID)

cross_arm_patients <- intersect(morphine_patients$PAT_ID, oxycodone_patients$PAT_ID)
cat("Children in both arms:", length(cross_arm_patients), "\n")

if (length(cross_arm_patients) > 0) {
  s3_morphine_keys <- s3_morphine_keys |>
    left_join(pat_lookup, by = "PAT_ENC_CSN_ID") |>
    filter(!(PAT_ID %in% cross_arm_patients)) |>
    select(PAT_ENC_CSN_ID, episode_id)
  
  s3_oxycodone_keys <- s3_oxycodone_keys |>
    left_join(pat_lookup, by = "PAT_ENC_CSN_ID") |>
    filter(!(PAT_ID %in% cross_arm_patients)) |>
    select(PAT_ENC_CSN_ID, episode_id)
  
  cat("Pools after cross-arm exclusion - morphine:", nrow(s3_morphine_keys),
      "| oxycodone:", nrow(s3_oxycodone_keys), "\n")
}

# ============================================================
# SAVE KEYS
# ============================================================
processed <- "1 - data/3 - processed_data/"
saveRDS(s4_medical_keys,   paste0(processed, "s4_medical_keys.rds"))
saveRDS(s4_surgical_keys,  paste0(processed, "s4_surgical_keys.rds"))
saveRDS(s3_morphine_keys,  paste0(processed, "s3_morphine_keys.rds"))
saveRDS(s3_oxycodone_keys, paste0(processed, "s3_oxycodone_keys.rds"))

# ============================================================
# SUMMARY
# ============================================================
cat("\n============================================================\n")
cat("COHORT SUMMARY\n")
cat("============================================================\n")
cat("Descriptive oxycodone population (episodes with an end time):\n")
cat(sprintf("  Medical:  %d (%.1f%%)\n",
            nrow(s4_medical_keys), 100 * nrow(s4_medical_keys) / nrow(all_oxy)))
cat(sprintf("  Surgical: %d (%.1f%%)\n\n",
            nrow(s4_surgical_keys), 100 * nrow(s4_surgical_keys) / nrow(all_oxy)))
cat("Primary comparative population, before matching:\n")
cat(sprintf("  Morphine:  %d\n", nrow(s3_morphine_keys)))
cat(sprintf("  Oxycodone: %d\n", nrow(s3_oxycodone_keys)))