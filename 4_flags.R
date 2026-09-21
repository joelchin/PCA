# ============================================================
# 4 — Episode flags, exclusions and final outcome classifications
# Paediatric PCA/NCA study — Addenbrooke's Hospital (CUH NHS FT)
# Author: J Chin
# ============================================================
#
# Purpose
#   Adds every episode-level flag used downstream to pca_episodes and saves the
#   result as pca_episodes_flagged (the single source of truth for episodes).
#
# Inputs
#   Files (1 - data/3 - processed_data/): admissions_clean, pca_episodes,
#     non_pca_mar_clean, nerve_blocks_clean, operations_clean, intraop_mar_clean,
#     lda_epidural_clean, anaesthesia_events_clean,
#     naloxone_reclassifications (manual chart-review recoding of naloxone
#     events; created outside this pipeline, patient-level, not shareable;
#     must contain PAT_ENC_CSN_ID, ADMIN_TIME, reclassified_to — checked on load).
#   In session: pain_scores_clean. drug_lists.R is sourced below.
#
# Output
#   pca_episodes_flagged.rds
#
# Flags
#   1  Episode sequence (first line, first episode)
#   2  Concurrent nerve block / epidural
#   3  Chemotherapy
#   4  Intraoperative drug exposures
#   5  Exclusion flags: exclude_minor_procedure_any (conservative) and
#      exclude_minor_procedure_any_strict, incl. exclude_pca_before_procedure.
#      These apply to the S3/PSM population only.
#   6  Demographics. proc_category and PFMD_PROC_OPCS_CODES arrive already
#      attached at episode level from 2_cleaningscript.R (Section 6D);
#      operations_clean is deliberately NOT re-joined here.
#   7  Era
#   8  Admission type. admission_type keeps its original categories
#      (backward compatibility); admission_type_detailed and
#      likely_nonopioid_nausea_confound are added alongside it.
#   9  Final side-effect classifications (reactive antiemetic, antipruritic,
#      naloxone reversal, composite) — definitions documented at Flag 9.
#      These are computed once here and read by scripts 10-14; they are not
#      recomputed elsewhere.
#
# Concurrent epidural (concurrent_epidural)
#   Defined from two definitive sources: the LDA epidural table (continuous
#   catheter) and the intrathecal MAR route (spinal). Epidural-route MAR entries
#   without an LDA record are not classified as concurrent epidural because
#   they most likely represent single-shot caudal blocks (neuraxial_epidural_mar_only
#   is retained as a sensitivity flag).
#   Methods text: "Epidural route entries in the intraoperative MAR without a
#   corresponding LDA record were not classified as concurrent epidural, as
#   these likely represent caudal blocks in the paediatric population."
# ============================================================

library(tidyverse)
library(lubridate)
source("drug_lists.R")

# ============================================================
# LOAD DATA
# ============================================================

pca_episodes             <- readRDS("1 - data/3 - processed_data/pca_episodes.rds")
admissions_clean         <- readRDS("1 - data/3 - processed_data/admissions_clean.rds")
operations_clean         <- readRDS("1 - data/3 - processed_data/operations_clean.rds")
nerve_blocks_clean       <- readRDS("1 - data/3 - processed_data/nerve_blocks_clean.rds")
anaesthesia_events_clean <- readRDS("1 - data/3 - processed_data/anaesthesia_events_clean.rds")
non_pca_mar_clean        <- readRDS("1 - data/3 - processed_data/non_pca_mar_clean.rds")
intraop_mar_clean        <- readRDS("1 - data/3 - processed_data/intraop_mar_clean.rds")
lda_epidural_clean       <- readRDS("1 - data/3 - processed_data/lda_epidural_clean.rds")

cat("All data loaded\n")
cat("Episodes to flag:", nrow(pca_episodes), "\n")
cat("proc_category already attached (episode-level, from 2_cleaningscript.R):",
    "proc_category" %in% names(pca_episodes), "\n")

# ============================================================
# FLAG 1 - Episode sequence flags
# ============================================================

cat("\nAdding episode sequence flags...\n")

pca_episodes <- pca_episodes |>
  arrange(PAT_ENC_CSN_ID, episode_start) |>
  group_by(PAT_ENC_CSN_ID) |>
  mutate(
    episode_number   = row_number(),
    is_first_episode = episode_number == 1
  ) |>
  ungroup()

pca_episodes <- pca_episodes |>
  arrange(PAT_ENC_CSN_ID, drug, episode_start) |>
  group_by(PAT_ENC_CSN_ID, drug) |>
  mutate(
    drug_episode_number   = row_number(),
    is_first_drug_episode = drug_episode_number == 1
  ) |>
  ungroup()

pca_episodes <- pca_episodes |>
  arrange(PAT_ENC_CSN_ID, episode_start) |>
  group_by(PAT_ENC_CSN_ID) |>
  mutate(
    prior_drug    = lag(drug),
    is_first_line = is.na(prior_drug),
    is_rotation   = !is_first_line & drug != prior_drug
  ) |>
  ungroup()

cat("First episodes (any drug):",    sum(pca_episodes$is_first_episode), "\n")
cat("First drug episodes:",          sum(pca_episodes$is_first_drug_episode), "\n")
cat("First line episodes:",          sum(pca_episodes$is_first_line), "\n")
cat("Rotation episodes:",            sum(pca_episodes$is_rotation), "\n")

# ============================================================
# FLAG 2 - Concurrent nerve block flags
#
# Drawn from three sources, two of them definitive:
#
# Source 1: nerve_blocks table
#   Note: no neuraxial entries in current dataset —
#   epidurals documented via LDA system, not nerve blocks
#   flowsheet. Peripheral blocks captured correctly.
#
# Source 2: Intraop MAR route — split into two flags:
#   neuraxial_intrathecal_mar = spinal (definitive)
#   neuraxial_epidural_mar    = epidural route (uncertain)
#     epidural route without LDA = likely caudal
#
# Source 3: LDA epidural table (definitive)
#   had_lda_epidural = continuous epidural catheter placed
#
# concurrent_epidural (PRIMARY) = had_lda_epidural OR
#                                  neuraxial_intrathecal_mar
#
# neuraxial_epidural_mar_only (SENSITIVITY) = epidural MAR,
#   no LDA, no intrathecal — uncertain/likely caudal group
# ============================================================

cat("\nAdding nerve block flags...\n")

# Source 1: nerve_blocks table
nerve_block_episodes <- nerve_blocks_clean |>
  filter(!is.na(block_start_time)) |>
  left_join(
    pca_episodes |>
      select(PAT_ENC_CSN_ID, episode_id,
             episode_start, episode_end),
    by = "PAT_ENC_CSN_ID",
    relationship = "many-to-many"
  ) |>
  filter(
    block_start_time <= episode_end,
    block_end_time >= episode_start | is.na(block_end_time)
  ) |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  summarise(
    had_any_block         = TRUE,
    had_neuraxial         = any(block_category == "neuraxial",
                                na.rm = TRUE),
    had_peripheral        = any(block_category %in%
                                  c("lower_limb", "upper_limb", "truncal"),
                                na.rm = TRUE),
    had_unspecified_block = any(block_category == "unspecified",
                                na.rm = TRUE),
    block_categories      = paste(unique(block_category),
                                  collapse = ", "),
    .groups = "drop"
  )

pca_episodes <- pca_episodes |>
  left_join(nerve_block_episodes,
            by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(
    had_any_block         = replace_na(had_any_block, FALSE),
    had_neuraxial         = replace_na(had_neuraxial, FALSE),
    had_peripheral        = replace_na(had_peripheral, FALSE),
    had_unspecified_block = replace_na(had_unspecified_block, FALSE)
  )

# Add anaesthesia_csn_id if not already present
if (!"anaesthesia_csn_id" %in% names(pca_episodes)) {
  pca_episodes <- pca_episodes |>
    left_join(
      anaesthesia_events_clean |>
        select(PAT_ENC_CSN_ID, anaesthesia_csn_id),
      by = "PAT_ENC_CSN_ID"
    )
  cat("anaesthesia_csn_id added to pca_episodes\n")
}

# Source 2: intraop MAR route — split intrathecal from epidural
neuraxial_mar_flags <- intraop_mar_clean |>
  mutate(route_clean = coalesce(mar_route, "") |> str_to_lower()) |>
  group_by(PAT_ENC_CSN_ID, anaesthesia_csn_id) |>
  summarise(
    neuraxial_intrathecal_mar = any(
      str_detect(route_clean, "intrathecal"), na.rm = TRUE
    ),
    neuraxial_epidural_mar = any(
      str_detect(route_clean, "epidural"), na.rm = TRUE
    ),
    neuraxial_from_mar = any(
      str_detect(route_clean, "intrathecal|epidural"), na.rm = TRUE
    ),
    .groups = "drop"
  )

pca_episodes <- pca_episodes |>
  left_join(neuraxial_mar_flags,
            by = c("PAT_ENC_CSN_ID", "anaesthesia_csn_id")) |>
  mutate(
    neuraxial_intrathecal_mar = replace_na(neuraxial_intrathecal_mar, FALSE),
    neuraxial_epidural_mar    = replace_na(neuraxial_epidural_mar, FALSE),
    neuraxial_from_mar        = replace_na(neuraxial_from_mar, FALSE)
  )

# Source 3: LDA epidural table — continuous catheter (definitive)
lda_epidural_csns <- lda_epidural_clean |>
  distinct(PAT_ENC_CSN_ID) |>
  mutate(had_lda_epidural = TRUE)

pca_episodes <- pca_episodes |>
  left_join(lda_epidural_csns, by = "PAT_ENC_CSN_ID") |>
  mutate(
    had_lda_epidural = replace_na(had_lda_epidural, FALSE),

    # PRIMARY: definitive neuraxial sources only
    concurrent_epidural = had_lda_epidural,

    # concurrent_epidural = LDA epidural table only
    # Rationale: LDA captures continuous running catheters (definitive)
    # Intrathecal MAR (spinals) excluded — single intraop dose,
    # typically wears off within 12-18hrs; PCA remains primary
    # postoperative analgesic for majority of episode
    # neuraxial_intrathecal_mar retained as PSM covariate

    # SENSITIVITY: epidural MAR without LDA — likely caudals
    neuraxial_epidural_mar_only = neuraxial_epidural_mar &
      !had_lda_epidural &
      !neuraxial_intrathecal_mar,

    # RETAINED for compatibility
    neuraxial_mar_only = !had_neuraxial & neuraxial_from_mar
  )

cat("Episodes with any peripheral block:",
    sum(pca_episodes$had_any_block), "\n")
cat("Episodes with neuraxial (nerve blocks table):",
    sum(pca_episodes$had_neuraxial), "\n")
cat("Episodes with intrathecal MAR (spinal):",
    sum(pca_episodes$neuraxial_intrathecal_mar), "\n")
cat("Episodes with epidural MAR (includes caudals):",
    sum(pca_episodes$neuraxial_epidural_mar), "\n")
cat("Episodes with LDA epidural (continuous catheter):",
    sum(pca_episodes$had_lda_epidural), "\n")
cat("Episodes with concurrent epidural (LDA + intrathecal MAR):",
    sum(pca_episodes$concurrent_epidural), "\n")
cat("Episodes with epidural MAR only / no LDA (likely caudals):",
    sum(pca_episodes$neuraxial_epidural_mar_only), "\n")
cat("Episodes with peripheral block only:",
    sum(pca_episodes$had_peripheral &
          !pca_episodes$concurrent_epidural), "\n")
cat("\nConcurrent epidural by drug:\n")
print(table(pca_episodes$drug, pca_episodes$concurrent_epidural))

# ============================================================
# FLAG 3 - Chemotherapy flag
# ============================================================

cat("\nAdding chemotherapy flag...\n")

chemo_patients <- non_pca_mar_clean |>
  filter(drug_category == "chemotherapy") |>
  distinct(PAT_ENC_CSN_ID) |>
  mutate(had_chemotherapy = TRUE)

pca_episodes <- pca_episodes |>
  left_join(chemo_patients, by = "PAT_ENC_CSN_ID") |>
  mutate(had_chemotherapy = replace_na(had_chemotherapy, FALSE))

cat("Episodes with chemotherapy:",
    sum(pca_episodes$had_chemotherapy), "\n")
cat("\nChemotherapy by drug:\n")
print(table(pca_episodes$drug, pca_episodes$had_chemotherapy))

# ============================================================
# FLAG 4 - Intraoperative drug flags
# ============================================================

cat("\nAdding intraoperative drug flags...\n")

intraop_drug_flags <- intraop_mar_clean |>
  mutate(
    drug_name_clean = coalesce(generic_name, "") |>
      str_squish() |>
      str_to_lower(),
    route_clean = coalesce(mar_route, "") |>
      str_to_lower(),

    had_remifentanil_drug = as.integer(
      str_detect(drug_name_clean, "remifentanil")),

    had_opioid_drug = as.integer(
      str_detect(drug_name_clean,
                 "morphine|fentanyl|alfentanil|oxycodone") &
        !str_detect(drug_name_clean, "diamorphine") &
        !str_detect(route_clean, "intrathecal")
    ),

    had_intrathecal_opioid_drug = as.integer(
      str_detect(drug_name_clean, "diamorphine") |
        (str_detect(drug_name_clean,
                    "morphine|fentanyl|oxycodone") &
           str_detect(route_clean, "intrathecal"))
    ),

    had_ketamine_drug  = as.integer(
      str_detect(drug_name_clean, "ketamine")),
    had_clonidine_drug = as.integer(
      str_detect(drug_name_clean, "clonidine"))
  ) |>
  group_by(PAT_ENC_CSN_ID, anaesthesia_csn_id) |>
  summarise(
    had_remifentanil       = max(had_remifentanil_drug) == 1L,
    had_intraop_opioid     = max(had_opioid_drug) == 1L,
    had_intrathecal_opioid = max(had_intrathecal_opioid_drug) == 1L,
    had_intraop_ketamine   = max(had_ketamine_drug) == 1L,
    had_intraop_clonidine  = max(had_clonidine_drug) == 1L,
    .groups = "drop"
  )

pca_episodes <- pca_episodes |>
  left_join(
    intraop_drug_flags,
    by = c("PAT_ENC_CSN_ID", "anaesthesia_csn_id")
  ) |>
  mutate(
    had_remifentanil       = replace_na(had_remifentanil, FALSE),
    had_intraop_opioid     = replace_na(had_intraop_opioid, FALSE),
    had_intrathecal_opioid = replace_na(had_intrathecal_opioid, FALSE),
    had_intraop_ketamine   = replace_na(had_intraop_ketamine, FALSE),
    had_intraop_clonidine  = replace_na(had_intraop_clonidine, FALSE)
  )

cat("Episodes with remifentanil:",
    sum(pca_episodes$had_remifentanil), "\n")
cat("Episodes with intraop opioid (systemic only):",
    sum(pca_episodes$had_intraop_opioid), "\n")
cat("Episodes with intrathecal opioid:",
    sum(pca_episodes$had_intrathecal_opioid), "\n")
cat("Episodes with intraop ketamine:",
    sum(pca_episodes$had_intraop_ketamine), "\n")
cat("Episodes with intraop clonidine:",
    sum(pca_episodes$had_intraop_clonidine), "\n")

# ============================================================
# FLAG 5 - the predefined exclusion flags
# ============================================================

cat("\nAdding the predefined exclusion flags...\n")

# opcs_exclude now sourced from drug_lists.R (single source of truth,
# shared with 2_cleaningscript.R's primary_opcs_substantive logic) —
# local definition REMOVED, was previously duplicated here.

# ============================================================
# all_opcs_excluded — FIXED
#
# Previously checked only primary_opcs per operation — the same
# first-code-only bug found in proc_category's derivation. A row
# with a substantive companion code (e.g. N13.4 + T30.9, where T30.9
# is exploratory laparotomy) would have been wrongly counted as
# "excluded" if primary_opcs happened to be the minor code, even
# though a real companion procedure was also coded on that row.
#
# FIX: checks ALL codes in ALL operations for the admission — an
# admission only counts as all_opcs_excluded if EVERY code in EVERY
# operation row is on opcs_exclude, using all_opcs_in_row (the full
# per-row code list from 2_cleaningscript.R), not just primary_opcs.
# ============================================================

opcs_flags <- operations_clean |>
  mutate(
    row_fully_minor = map_lgl(
      all_opcs_in_row,
      ~ all(.x[.x != ""] %in% opcs_exclude)
    )
  ) |>
  group_by(PAT_ENC_CSN_ID) |>
  summarise(
    total_operations        = n_distinct(OP_ORDER),
    n_excluded_opcs         = sum(row_fully_minor, na.rm = TRUE),
    all_opcs_excluded       = all(row_fully_minor),
    has_multiple_operations = max(OP_ORDER, na.rm = TRUE) > 1,
    .groups = "drop"
  )

pca_episodes <- pca_episodes |>
  left_join(opcs_flags, by = "PAT_ENC_CSN_ID") |>
  mutate(
    all_opcs_excluded       = replace_na(all_opcs_excluded, FALSE),
    has_multiple_operations = replace_na(has_multiple_operations, FALSE)
  )

anaesthesia_flags <- anaesthesia_events_clean |>
  group_by(PAT_ENC_CSN_ID) |>
  summarise(
    anaesthesia_duration_mins = first(anaesthesia_duration_mins),
    anaesthesia_stop = suppressWarnings(
      min(anaesthesia_stop[is.finite(anaesthesia_stop)],
          na.rm = TRUE)
    ),
    .groups = "drop"
  ) |>
  mutate(
    anaesthesia_stop = if_else(
      is.infinite(anaesthesia_stop), NA_POSIXct_, anaesthesia_stop
    ),
    exclude_short_anaesthesia = !is.na(anaesthesia_duration_mins) &
      anaesthesia_duration_mins < 20,
    exclude_long_anaesthesia  = !is.na(anaesthesia_duration_mins) &
      anaesthesia_duration_mins > 960
  )

pca_episodes <- pca_episodes |>
  left_join(
    anaesthesia_flags |>
      select(PAT_ENC_CSN_ID, anaesthesia_stop,
             exclude_short_anaesthesia,
             exclude_long_anaesthesia),
    by = "PAT_ENC_CSN_ID"
  ) |>
  mutate(
    exclude_short_anaesthesia = replace_na(exclude_short_anaesthesia, FALSE),
    exclude_long_anaesthesia  = replace_na(exclude_long_anaesthesia, FALSE)
  )

pca_episodes <- pca_episodes |>
  mutate(
    hrs_anaesthesia_to_pca = if_else(
      !is.na(anaesthesia_stop),
      as.numeric(difftime(episode_start, anaesthesia_stop,
                          units = "hours")),
      NA_real_
    ),
    exclude_late_pca_start = case_when(
      is.na(hrs_anaesthesia_to_pca) ~ NA,
      hrs_anaesthesia_to_pca > 24   ~ TRUE,
      TRUE                          ~ FALSE
    ),
    exclude_short_pca = !is.na(episode_duration_hrs) &
      episode_duration_hrs < 2
  )

# ============================================================
# exclude_pca_before_procedure
#
# Formalizes the pre-procedure-timing diagnostic into a genuine
# exclusion criterion. Catches PCA starting BEFORE the admission's
# anaesthesia_stop at all — the opposite direction from
# exclude_late_pca_start (which catches PCA starting >24h AFTER).
#
# RESTRICTED TO SINGLE-OPERATION ADMISSIONS ONLY. anaesthesia_events
# has exactly one row per admission (confirmed on review, 0/3,727
# admissions have >1 row), but operations are frequently on
# genuinely DIFFERENT calendar days. This means the captured
# anaesthesia_stop cannot be reliably attributed to a specific
# operation when multiple operations occurred (e.g. a child has a
# major first operation, PCA legitimately starts right after it,
# then returns to theatre later — if the captured anaesthesia_stop
# belongs to that LATER visit rather than the real index operation,
# this flag would wrongly fire). Restricting to single-op admissions
# avoids acting on that ambiguity.
#
# NOTE: this restriction does not change exclude_minor_procedure_any's final
# output in practice — has_multiple_operations already excludes
# every multi-op admission via its own criterion regardless, so this
# flag only ever contributes new information for single-op
# admissions where has_multiple_operations doesn't already apply.
# ============================================================

# ============================================================
# exclude_pca_before_procedure — WITH 2-HOUR BUFFER
#
# Buffer added after validating the actual gap distribution: of 404
# episodes originally caught by a bare "any negative gap" rule, 60
# (14.9%) had gaps under 30 minutes — plausibly PCA being programmed
# in theatre before extubation, ordinary workflow, not genuine
# incidental timing. Only 16 more (4.0%) fell in the 30min-2hr
# ambiguous range. The remaining 328 (81.2%) showed gaps of 2+ hours
# (median ~14.6h, matching the scale already validated in Male
# genitalia/Vascular-lines), confirming the flag genuinely detects a
# real, broader-than-expected phenomenon — NOT primarily charting
# noise. The 2-hour threshold below excludes only that validated
# majority, leaving the ambiguous/likely-noise 76 episodes
# unexcluded rather than guessing either way on them.
# ============================================================

PCA_BEFORE_PROCEDURE_BUFFER_HOURS <- 2

pca_episodes <- pca_episodes |>
  mutate(
    gap_hrs_procedure_to_pca = if_else(
      !is.na(anaesthesia_stop),
      as.numeric(difftime(anaesthesia_stop, episode_start, units = "hours")),
      NA_real_
    ),
    exclude_pca_before_procedure =
      !has_multiple_operations &
      !is.na(anaesthesia_stop) &
      gap_hrs_procedure_to_pca > PCA_BEFORE_PROCEDURE_BUFFER_HOURS
  )

cat("PCA before procedure, single-op only, >2h buffer applied:",
    sum(pca_episodes$exclude_pca_before_procedure), "\n")

cat("\nexclude_pca_before_procedure_opdate_fallback and\n")
cat("exclude_missing_proc_category MOVED to after Flag 6 — both need\n")
cat("case_type, which doesn't exist on pca_episodes until the\n")
cat("demographics join runs. exclude_minor_procedure_any/exclude_minor_procedure_any_strict\n")
cat("are correspondingly finalized after Flag 6 too, not here.\n\n")

cat("\nPredefined exclusion flag summary:\n")
cat("Non-substantive OPCS:",
    sum(pca_episodes$all_opcs_excluded), "\n")
cat("Multiple operations:",
    sum(pca_episodes$has_multiple_operations), "\n")
cat("Short anaesthesia (<20 mins):",
    sum(pca_episodes$exclude_short_anaesthesia), "\n")
cat("Long anaesthesia (>16 hrs):",
    sum(pca_episodes$exclude_long_anaesthesia), "\n")
cat("Late PCA start (>24hrs):",
    sum(pca_episodes$exclude_late_pca_start, na.rm = TRUE),
    "(NA =", sum(is.na(pca_episodes$exclude_late_pca_start)), ")\n")
cat("(Final exclude_minor_procedure_any, incl. case_type-dependent criteria,\n")
cat("computed after Flag 6 — see below)\n")

# ============================================================
# has_reliable_anaesthesia_timing
#
# Closes a gap deliberately deferred earlier in development. For
# multi-operation admissions, anaesthesia_events only ever captures
# ONE anaesthetic record per admission (confirmed: 0/3,727 admissions
# have >1 row) — the OLD admission-level anaesthesia_stop join gave
# that single record to EVERY episode in the admission regardless of
# which operation it actually followed, and that logic (in
# anaesthesia_flags above) is unchanged by the fixes above.
#
# This flag identifies exactly which episodes' anaesthesia_stop is
# genuinely trustworthy: TRUE for single-operation admissions
# (trivially unambiguous), and for multi-op admissions, TRUE only
# when the episode's 6D-attributed operation (attributed_op_order)
# IS the specific operation with genuine ID-matched anaesthesia data
# (exact ANAESTHESIA_CSN_ID match, not date-proximity guessing).
# FALSE otherwise — those episodes' anaesthesia_stop is confirmed
# borrowed from a DIFFERENT operation.
#
# Validated (global, all drugs) prior to the proc_category
# fixes: 4,765/5,181 (92.0%) reliable among episodes with any
# anaesthesia data at all; 416 episodes definitively mis-anchored.
# Re-check this split on the actual oxycodone S4 population before
# deciding how to treat the FALSE group for antiemetic windowing —
# see Flag 9 below for the two options.
# ============================================================

op_with_real_anaesthesia <- operations_clean |>
  select(PAT_ENC_CSN_ID, OP_ORDER, ANAESTHESIA_CSN_ID) |>
  inner_join(
    anaesthesia_events_clean |>
      select(anaesthesia_csn_id, anaesthesia_stop) |>
      rename(ANAESTHESIA_CSN_ID = anaesthesia_csn_id) |>
      filter(!is.na(anaesthesia_stop)),
    by = "ANAESTHESIA_CSN_ID"
  ) |>
  distinct(PAT_ENC_CSN_ID, .keep_all = TRUE) |>
  select(PAT_ENC_CSN_ID, reliable_op_order = OP_ORDER)

pca_episodes <- pca_episodes |>
  left_join(op_with_real_anaesthesia, by = "PAT_ENC_CSN_ID") |>
  mutate(
    has_reliable_anaesthesia_timing = case_when(
      !has_multiple_operations ~ TRUE,
      has_multiple_operations & !is.na(reliable_op_order) &
        attributed_op_order == reliable_op_order ~ TRUE,
      TRUE ~ FALSE
    )
  )

cat("\nhas_reliable_anaesthesia_timing (oxycodone S4):\n")
print(table(
  pca_episodes$drug == "oxycodone",
  pca_episodes$has_reliable_anaesthesia_timing
))

# ============================================================
# FLAG 6 - Demographics
#
# operations_clean join REMOVED.
# proc_category and PFMD_PROC_OPCS_CODES now arrive on pca_episodes
# ALREADY ATTACHED at episode level from 2_cleaningscript.R Section
# 6B (matched to the operation nearest in time before each specific
# episode, via OpDate — not the old admission-level OP_ORDER==1
# join). Re-joining operations_clean here would either overwrite
# those corrected values with the old, less accurate admission-level
# ones, or create duplicate .x/.y columns — both wrong. Only the
# admissions_clean demographics join remains.
# ============================================================

cat("\nAdding demographics...\n")

pca_episodes <- pca_episodes |>
  left_join(
    admissions_clean |>
      select(PAT_ENC_CSN_ID, AGE, age_days, age_group,
             WEIGHT_KG, GENDER, ADMITTING_SPECIALTY,
             case_type, LOS_days),
    by = "PAT_ENC_CSN_ID"
  )

cat("Demographics joined (proc_category already present from source)\n")
cat("proc_category present:", "proc_category" %in% names(pca_episodes), "\n")
cat("PFMD_PROC_OPCS_CODES present:", "PFMD_PROC_OPCS_CODES" %in% names(pca_episodes), "\n")

# ============================================================
# case_type-DEPENDENT EXCLUSION FLAGS
#
# MOVED HERE from Flag 5, where case_type did not yet exist on
# pca_episodes (it's only attached by the demographics join directly
# above). Referencing case_type inside Flag 5 either errored or —
# more dangerously — silently picked up a stale case_type object left
# over in the R session from earlier work, misaligned with the actual
# row order. Confirmed via post-rebuild check: exactly one S3 episode
# still had NA proc_category despite exclude_missing_proc_category
# supposedly covering that exact pattern — this placement bug is why.
# ============================================================

# ---- exclude_pca_before_procedure_opdate_fallback ----
# Catches PCA-predates-only-recorded-operation cases where
# anaesthesia_stop is NA (so exclude_pca_before_procedure's own
# !is.na(anaesthesia_stop) precondition can't fire) — confirmed via
# OpDate instead (e.g. PCA started 2016-02-07 22:55, sole operation
# dated 2016-02-08).
pca_episodes <- pca_episodes |>
  mutate(
    exclude_pca_before_procedure_opdate_fallback =
      !has_multiple_operations &
      is.na(anaesthesia_stop) &
      case_type == "surgical" &
      is.na(proc_category)
  )

cat("\nPCA before procedure, OpDate fallback (anaesthesia_stop NA cases):",
    sum(pca_episodes$exclude_pca_before_procedure_opdate_fallback, na.rm = TRUE), "\n")

# ---- exclude_missing_proc_category ----
# GENERAL CATCH-ALL: any surgical-labelled episode with no procedure
# category at all cannot meaningfully enter a procedure-matched
# comparison, regardless of which specific mechanism left the
# category unset.
pca_episodes <- pca_episodes |>
  mutate(
    exclude_missing_proc_category = case_type == "surgical" & is.na(proc_category)
  )

cat("Missing proc_category, general catch-all:",
    sum(pca_episodes$exclude_missing_proc_category, na.rm = TRUE), "\n")

# ---- FINAL exclude_minor_procedure_any / exclude_minor_procedure_any_strict ----
# Computed HERE, not in Flag 5 — this is the first point where every
# contributing flag (including the two case_type-dependent ones just
# above) actually exists.
pca_episodes <- pca_episodes |>
  mutate(
    exclude_minor_procedure_any = all_opcs_excluded |
      has_multiple_operations |
      exclude_short_anaesthesia |
      exclude_long_anaesthesia |
      coalesce(exclude_late_pca_start, FALSE) |
      exclude_short_pca |
      exclude_pca_before_procedure |
      exclude_pca_before_procedure_opdate_fallback |
      exclude_missing_proc_category,
    exclude_minor_procedure_any_strict = all_opcs_excluded |
      has_multiple_operations |
      exclude_short_anaesthesia |
      exclude_long_anaesthesia |
      replace_na(exclude_late_pca_start, FALSE) |
      exclude_short_pca |
      exclude_pca_before_procedure |
      exclude_pca_before_procedure_opdate_fallback |
      exclude_missing_proc_category
  )

cat("\nAny Brooks exclusion (conservative, FINAL):",
    sum(pca_episodes$exclude_minor_procedure_any), "\n")
cat("Any Brooks exclusion (strict/sensitivity, FINAL):",
    sum(pca_episodes$exclude_minor_procedure_any_strict), "\n")

cat("\nVerification — surgical episodes with NA proc_category NOT\n")
cat("caught by exclude_minor_procedure_any (should be 0):\n")
cat(sum(pca_episodes$case_type == "surgical" &
          is.na(pca_episodes$proc_category) &
          !pca_episodes$exclude_minor_procedure_any, na.rm = TRUE), "\n")

# ============================================================
# case_type_refined
#
# case_type (above, from admissions_clean) is a blunt admission-level
# presence check: "does ANY OPCS-coded operation exist for this
# admission" → surgical, else medical. This doesn't ask whether the
# operation was substantive enough to plausibly explain PCA need.
#
# case_type_refined reclassifies an episode to "medical" when ANY of:
#   (a) its 6D-attributed operation is entirely minor (row_fully_minor
#       — every code on that specific operation record is on the
#       minor list, incl. N13.4 + 15 further codes closed out this
#       session)
#   (b) EVERY operation across the WHOLE ADMISSION is minor
#       (all_opcs_excluded, NEW here) — a wider-scope version of (a):
#       row_fully_minor checks only the one operation this specific
#       episode was attributed to; all_opcs_excluded checks the
#       entire admission. Direct match to the actual test being
#       applied: "would this surgery be deemed to require PCA/NCA" —
#       if every operation present is on the minor list, none of them
#       would be.
#   (c) no operation precedes it at all (proc_category is NA — also
#       covers exclude_missing_proc_category and its OpDate fallback,
#       both of which are strict subsets of this same condition, so
#       not added separately)
#   (d) exclude_pca_before_procedure — a preceding operation WAS found
#       via OpDate (proc_category is NOT NA), but the actual
#       anaesthesia timing shows PCA substantially predates it (>2h
#       buffer, validated median gap ~14.6h across 328 episodes) —
#       the labelled procedure technically precedes the episode by
#       date, but the timing evidence says it isn't genuinely what's
#       driving PCA need.
#   (e) exclude_short_anaesthesia (<20min, NEW here) — originally
#       built as a data-quality sanity check rather than a deliberate
#       severity threshold, but a genuinely brief anaesthetic does
#       plausibly indicate a minor procedure — included per explicit
#       decision.
#
# DELIBERATELY NOT INCLUDED, despite being part of exclude_minor_procedure_any —
# these answer a different question ("is this episode usable for
# clean PSM matching"), not "was the labelled procedure severe enough
# to require PCA":
#   - has_multiple_operations: about WHICH operation drove PCA use,
#     not whether either operation was severe enough. A child with
#     two real operations is still genuinely surgical.
#   - exclude_long_anaesthesia (>16h): flagged as a likely data-entry
#     error, and if taken at face value points the OPPOSITE direction
#     (longer anaesthesia suggests MORE severe surgery, not less).
#   - exclude_late_pca_start (>24h after): a timing anomaly, not a
#     statement about procedure severity.
#   - exclude_short_pca (<2h duration): ambiguous — could reflect a
#     genuinely minor procedure, or an unrelated early switch to a
#     different analgesic. Not a clean severity signal either way.
#
# DELIBERATELY A NEW FIELD, NOT AN OVERWRITE of case_type — every
# population definition (S3/S4, PSM) currently filters on
# case_type == "surgical"; overwriting it would silently change every
# locked population count. case_type_refined is available for use
# where this distinction matters, without disturbing anything already
# built on the original case_type.
# ============================================================

pca_episodes <- pca_episodes |>
  mutate(
    case_type_refined = if_else(
      case_type == "surgical" &
        (coalesce(row_fully_minor, FALSE) |
           coalesce(all_opcs_excluded, FALSE) |
           coalesce(exclude_short_anaesthesia, FALSE) |
           is.na(proc_category) |
           coalesce(exclude_pca_before_procedure, FALSE)),
      "medical",
      case_type
    )
  )

cat("\ncase_type vs case_type_refined:\n")
print(table(pca_episodes$case_type, pca_episodes$case_type_refined, useNA = "ifany"))
cat("\nEpisodes reclassified surgical -> medical:",
    sum(pca_episodes$case_type == "surgical" &
          pca_episodes$case_type_refined == "medical", na.rm = TRUE), "\n")

# ============================================================
# FLAG 7 - Era
# ============================================================

pca_episodes <- pca_episodes |>
  mutate(
    year = year(episode_start),
    era  = case_when(
      year <= 2018 ~ "Early (2014-2018)",
      year >= 2019 ~ "Recent (2019-2026)",
      TRUE         ~ NA_character_
    )
  )

cat("\nEra distribution:\n")
print(table(pca_episodes$era, pca_episodes$drug))

# ============================================================
# FLAG 8 - Admission type
#
# Original admission_type KEPT UNCHANGED below for backward
# compatibility (a separate analysis depends on its exact 4
# categories). New fields added alongside:
#   - is_trauma: standalone flag (trauma can co-occur with any
#     proc_category, not just orthopaedic — no longer swallowed
#     into a single "Trauma" bucket that hides this)
#   - medical_specialty_grouped: collapses ADMITTING_SPECIALTY's
#     60+ raw values into clinically meaningful groups
#   - admission_type_detailed: full granularity — surgical uses
#     proc_category (NOW episode-level, from 2_cleaningscript.R
#     6B), medical uses medical_specialty_grouped
#   - likely_nonopioid_nausea_confound: oncology/haematology/GI
#     admissions flagged directly, since these plausibly have
#     non-opioid causes of nausea (chemo, underlying GI disease)
#     that a burden-percentage threshold cannot reliably separate
#     from true opioid-driven nausea
# ============================================================

specialty_groups <- c(
  "PAEDIATRIC MEDICAL ONCOLOGY"         = "Oncology",
  "MEDICAL ONCOLOGY"                    = "Oncology",
  "PAEDIATRIC CLINICAL HAEMATOLOGY"     = "Haematology/Oncology",
  "CLINICAL HAEMATOLOGY"                = "Haematology/Oncology",
  "PAEDIATRIC GASTROENTEROLOGY"         = "Gastroenterology",
  "GASTROENTEROLOGY"                    = "Gastroenterology",
  "PAEDIATRIC GASTROINTESTINAL SURGERY" = "Gastroenterology",
  "PAEDIATRIC RESPIRATORY MEDICINE"     = "Respiratory medicine",
  "PAEDIATRIC NEUROLOGY"                = "Neurology",
  "PAEDIATRIC INTENSIVE CARE"           = "PICU",
  "PAEDIATRICS"                         = "General paediatrics",
  "NEONATOLOGY"                         = "Neonatology",
  "WELL BABIES"                         = "Neonatology",
  "ACCIDENT & EMERGENCY"                = "A&E (unclear cause)"
)

# ============================================================
# Confirmed PRE-OP SYMPTOMATIC — patient-level evidence of
# pre-existing nausea/vomiting (antiemetic given before treatment
# started), rather than inferring it from specialty/category alone.
# catches genuine cases a category-based
# exclusion misses (e.g. non-neurosurgery patient already
# symptomatic) and avoids wrongly excluding category members who were
# NOT actually pre-op symptomatic (most neurosurgery patients).
# Combined with category flag as OR: the two catch non-overlapping
# cases (43 pre-op-only, 156 category-only, 51 both, out of 1023
# oxycodone episodes checked).
#
# This flag is anchored to episode_start rather than
# anaesthesia_stop ("before surgery"), since anchoring to
# anaesthesia_stop is unreliable for medical episodes: 68.5%
# (998/1457) either have anaesthesia_stop == NA (which would leave
# the flag silently defaulting to FALSE — "unknown", not "confirmed
# not symptomatic") or have it populated but irrelevant to this
# specific episode (proc_category NA — the anaesthetic belongs to an
# unrelated event on the admission). Anchoring to episode_start gives
# 100% coverage and more precisely answers the question this flag is
# meant to answer: was the child already symptomatic before opioid
# treatment started, not specifically before an operation. For
# surgical episodes this makes little practical difference
# (episode_start and anaesthesia_stop are normally close together —
# PCA typically starts shortly after surgery ends); for medical
# episodes it's what makes the flag usable at all.
# ============================================================

preop_antiemetic_events <- non_pca_mar_clean |>
  filter(given == TRUE, tolower(GENERIC_NAME) %in% tolower(antiemetic_drugs)) |>
  select(PAT_ENC_CSN_ID, ADMIN_TIME)

preop_symptomatic_lookup <- pca_episodes |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start) |>
  left_join(preop_antiemetic_events, by = "PAT_ENC_CSN_ID", relationship = "many-to-many") |>
  mutate(
    treatment_start_date = as.Date(episode_start),
    admin_date            = as.Date(ADMIN_TIME)
  ) |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  summarise(
    confirmed_preop_symptomatic = any(admin_date < treatment_start_date, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nconfirmed_preop_symptomatic — coverage check (should be ~100%\n")
cat("computable now, vs the previous anaesthesia_stop-anchored version):\n")
cat("Total episodes:", nrow(pca_episodes), "\n")
cat("Episodes where this flag is now computable (episode_start always\n")
cat("present, so effectively all):", sum(!is.na(pca_episodes$episode_start)), "\n\n")

pca_episodes <- pca_episodes |>
  left_join(preop_symptomatic_lookup, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(
    is_trauma = str_detect(coalesce(PFMD_PROC_OPCS_CODES, ""),
                           "^W1[0-9]|^W2[0-9]|^W30"),

    medical_specialty_grouped = recode(
      ADMITTING_SPECIALTY, !!!specialty_groups,
      .default = "Other medical/unclassified"
    ),

    admission_type_detailed = case_when(
      case_type == "surgical" & !is.na(proc_category) ~ proc_category,
      case_type == "surgical"                          ~ "Surgical (unclassified OPCS)",
      case_type == "medical"                            ~ paste0("Medical: ", medical_specialty_grouped),
      TRUE                                               ~ "Unclassified"
    ),

    confirmed_preop_symptomatic = replace_na(confirmed_preop_symptomatic, FALSE),

    likely_nonopioid_nausea_confound = medical_specialty_grouped %in%
      c("Oncology", "Haematology/Oncology", "Gastroenterology") |
      admission_type_detailed %in% c("Neurosurgery", "Vascular/lines"),

    # FINAL combined flag — use THIS for antiemetic burden analysis
    nonopioid_nausea_confound_final = likely_nonopioid_nausea_confound | confirmed_preop_symptomatic,

    # ORIGINAL — UNCHANGED, kept for backward compatibility (separate analysis)
    admission_type = case_when(
      str_detect(coalesce(PFMD_PROC_OPCS_CODES, ""),
                 "^W1[0-9]|^W2[0-9]|^W30") ~ "Trauma",
      proc_category == "Orthopaedic/spine" &
        !str_detect(coalesce(PFMD_PROC_OPCS_CODES, ""),
                    "^W1[0-9]|^W2[0-9]|^W30") ~ "Elective orthopaedic",
      case_type == "surgical" ~ "Other surgical",
      TRUE ~ "Medical"
    )
  )

cat("\nAdmission type distribution (original, unchanged):\n")
print(table(pca_episodes$admission_type, pca_episodes$drug))

cat("\nAdmission type distribution (detailed, NEW):\n")
print(table(pca_episodes$admission_type_detailed, pca_episodes$drug))

cat("\nTrauma flag distribution (NEW, standalone):\n")
print(table(pca_episodes$is_trauma))

cat("\nLikely non-opioid nausea confound (oncology/haem/GI) (NEW):\n")
print(table(pca_episodes$likely_nonopioid_nausea_confound))

# ============================================================
# ============================================================
# FLAG 9 - Final side-effect classifications
#
# Single source of truth for antiemetic/antipruritic/naloxone-
# reversal flags — read from here by both 10_comparative_PSM.R and
# 12_safety.R.
# ============================================================

cat("\nAdding final side-effect classifications...\n")

# naloxone_reclassifications.rds is not produced by any script in
# this pipeline — it's a manually curated file from chart review of
# naloxone administration events (see methodology notes on
# distinguishing pruritus-dose from respiratory-reversal-dose
# naloxone). It must exist on disk before running this script.
naloxone_reclassifications <- readRDS("1 - data/3 - processed_data/naloxone_reclassifications.rds")
missing_cols <- setdiff(c("PAT_ENC_CSN_ID", "ADMIN_TIME", "reclassified_to"),
                        names(naloxone_reclassifications))
if (length(missing_cols) > 0) {
  stop("naloxone_reclassifications.rds is missing column(s): ",
       paste(missing_cols, collapse = ", "))
}

non_pca_mar_reclassified <- non_pca_mar_clean |>
  left_join(
    naloxone_reclassifications |> select(PAT_ENC_CSN_ID, ADMIN_TIME, reclassified_to),
    by = c("PAT_ENC_CSN_ID", "ADMIN_TIME")
  ) |>
  mutate(naloxone_type = if_else(!is.na(reclassified_to), reclassified_to, naloxone_type)) |>
  select(-reclassified_to)

# ============================================================
# ANTIEMETIC DEFINITION
#
# Confound exclusion is confirmed_preop_symptomatic only. A
# specialty/category-level exclusion (e.g. flagging all neurosurgery
# patients) is deliberately not used, since it's too blunt — a
# neurosurgery patient can have genuine opioid-driven nausea
# unrelated to their neurological condition, and a blanket category
# exclusion would wrongly discard them.
#
# Episodes are split into three groups by anaesthesia-timing
# reliability (antiemetic_group, below): Group A has no anaesthesia
# data at all; Group B has anaesthesia data that reliably corresponds
# to this episode; Group C has anaesthesia data that exists but has
# been shown not to reliably correspond to this specific episode (see
# has_reliable_anaesthesia_timing) and is excluded from the reactive
# antiemetic denominator entirely, since post-anaesthesia antiemetic
# timing can't be interpreted without a trustworthy anchor.
#
# No burden-percentage ceiling is applied. A single fixed threshold
# was tested against rotation-away rate as an external validation
# check and showed no natural breakpoint at any threshold tried
# (25/50/75%/median) — rotation-away rate climbs smoothly and
# monotonically with burden rather than showing a genuine cliff, so no
# cutoff is used.
#
# The unadjusted rate by drug should be interpreted alongside the
# PSM-matched comparison in 10_comparative_PSM.R rather than on its own.
# ============================================================
#
# GROUP B — episodes with a real anaesthetic (anaesthesia_stop AND
# proc_category both populated): reactive if
#   (a) antiemetic given >12h after anaesthesia_stop, AND
#   (b) no antiemetic given in the 7 days before anaesthesia_start
#       (not episode_start, not anaesthesia_stop — anchoring to
#       anaesthesia_start specifically avoids catching routine
#       intraoperative antiemetic prophylaxis: anchoring to
#       episode_start/anaesthesia_stop instead gives a 63-68%
#       "pre-existing" rate, implausibly high and dominated by
#       routine dosing rather than genuine symptoms, whereas
#       anchoring to anaesthesia_start gives a clinically plausible
#       9.1%).
#
# GROUP A — episodes with no anaesthesia data (anaesthesia_stop or
# proc_category is NA): reactive if
#   (a) any antiemetic given during the episode, AND
#   (b) no antiemetic given in the 7 days before episode_start
# ============================================================

antiemetic_events <- non_pca_mar_reclassified |>
  filter(given == TRUE, tolower(GENERIC_NAME) %in% tolower(antiemetic_drugs)) |>
  select(PAT_ENC_CSN_ID, ADMIN_TIME)

stopifnot(exists("anaesthesia_events_clean"))

anaesthesia_timing_for_antiemetic <- anaesthesia_events_clean |>
  select(PAT_ENC_CSN_ID, anaesthesia_stop, anaesthesia_duration_mins) |>
  distinct(PAT_ENC_CSN_ID, .keep_all = TRUE) |>
  mutate(anaesthesia_start = anaesthesia_stop - (anaesthesia_duration_mins * 60))

pca_episodes <- pca_episodes |>
  left_join(anaesthesia_timing_for_antiemetic |> select(PAT_ENC_CSN_ID, anaesthesia_start),
            by = "PAT_ENC_CSN_ID") |>
  mutate(
    # ============================================================
    # antiemetic_group
    #
    # Split by anaesthesia-timing reliability rather than simply by
    # whether a GA occurred, decoupled from proc_category/
    # case_type_refined (a different question — procedure severity,
    # not GA occurrence):
    #   Group A — no anaesthesia_stop at all
    #   Group B — anaesthesia_stop present and reliable
    #   Group C — anaesthesia_stop present but unreliable for this
    #     episode (has_reliable_anaesthesia_timing == FALSE)
    #
    # Group C episodes are excluded from the reactive antiemetic
    # calculation entirely, rather than folded into either rule.
    # Neither rule's core assumption holds for them: Group A's "any
    # antiemetic during episode, no 12h filter" rule would wrongly
    # count ordinary early post-op PONV as reactive for episodes that
    # did have a real anaesthetic, and Group B's 12h-after-
    # anaesthesia-stop filter can't be trusted when anaesthesia_stop
    # itself doesn't reliably belong to this episode. Post-anaesthesia
    # antiemetic use can't be classified as reactive without a
    # trustworthy timing anchor.
    antiemetic_group = case_when(
      is.na(anaesthesia_stop) ~ "A",
      !has_reliable_anaesthesia_timing ~ "C",
      TRUE ~ "B"
    )
  )

cat("\nantiemetic_group split:\n")
print(table(pca_episodes$antiemetic_group))
cat("Group C — real anaesthesia but unreliable timing, EXCLUDED from\n")
cat("reactive antiemetic denominator entirely:",
    sum(pca_episodes$antiemetic_group == "C"), "\n\n")

# ---- GROUP A ----
group_a_episodes <- pca_episodes |> filter(antiemetic_group == "A")

group_a_excluded <- group_a_episodes |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start) |>
  left_join(antiemetic_events, by = "PAT_ENC_CSN_ID", relationship = "many-to-many") |>
  mutate(hrs_before = as.numeric(difftime(episode_start, ADMIN_TIME, units = "hours")),
         excl = !is.na(ADMIN_TIME) & hrs_before > 0 & hrs_before <= 168) |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  summarise(excluded_preexisting = any(excl, na.rm = TRUE), .groups = "drop")

group_a_any_during <- group_a_episodes |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end) |>
  left_join(antiemetic_events, by = "PAT_ENC_CSN_ID", relationship = "many-to-many") |>
  filter(!is.na(ADMIN_TIME), ADMIN_TIME >= episode_start, ADMIN_TIME <= episode_end) |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  summarise(any_antiemetic_during = TRUE, .groups = "drop")

group_a_result <- group_a_episodes |>
  select(PAT_ENC_CSN_ID, episode_id) |>
  left_join(group_a_excluded, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  left_join(group_a_any_during, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(
    excluded_preexisting = replace_na(excluded_preexisting, FALSE),
    any_antiemetic_during = replace_na(any_antiemetic_during, FALSE),
    reactive_antiemetic_final = any_antiemetic_during & !excluded_preexisting
  ) |>
  select(PAT_ENC_CSN_ID, episode_id, reactive_antiemetic_final)

# ---- GROUP B ----
group_b_episodes <- pca_episodes |> filter(antiemetic_group == "B")

group_b_excluded <- group_b_episodes |>
  select(PAT_ENC_CSN_ID, episode_id, anaesthesia_start) |>
  left_join(antiemetic_events, by = "PAT_ENC_CSN_ID", relationship = "many-to-many") |>
  mutate(hrs_before_start = as.numeric(difftime(anaesthesia_start, ADMIN_TIME, units = "hours")),
         excl = !is.na(ADMIN_TIME) & hrs_before_start > 0 & hrs_before_start <= 168) |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  summarise(excluded_preexisting = any(excl, na.rm = TRUE), .groups = "drop")

group_b_qualifying <- group_b_episodes |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end, anaesthesia_stop) |>
  left_join(antiemetic_events, by = "PAT_ENC_CSN_ID", relationship = "many-to-many") |>
  filter(!is.na(ADMIN_TIME), ADMIN_TIME >= episode_start, ADMIN_TIME <= episode_end) |>
  mutate(hrs_post_stop = as.numeric(difftime(ADMIN_TIME, anaesthesia_stop, units = "hours")),
         qualifies = hrs_post_stop > 12) |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  summarise(any_qualifying = any(qualifies, na.rm = TRUE), .groups = "drop")

group_b_result <- group_b_episodes |>
  select(PAT_ENC_CSN_ID, episode_id) |>
  left_join(group_b_excluded, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  left_join(group_b_qualifying, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(
    excluded_preexisting = replace_na(excluded_preexisting, FALSE),
    any_qualifying = replace_na(any_qualifying, FALSE),
    reactive_antiemetic_final = any_qualifying & !excluded_preexisting
  ) |>
  select(PAT_ENC_CSN_ID, episode_id, reactive_antiemetic_final)

antiemetic_final_result <- bind_rows(group_a_result, group_b_result)

cat("\n--- Antiemetic (FINAL, two-group design) ---\n")
cat("Group A:", sum(group_a_episodes$PAT_ENC_CSN_ID %in% group_a_result$PAT_ENC_CSN_ID), "episodes,",
    sum(group_a_result$reactive_antiemetic_final), "reactive (",
    round(100*mean(group_a_result$reactive_antiemetic_final), 1), "%)\n")
cat("Group B:", nrow(group_b_result), "episodes,",
    sum(group_b_result$reactive_antiemetic_final), "reactive (",
    round(100*mean(group_b_result$reactive_antiemetic_final), 1), "%)\n\n")

# ---- Antipruritic (no time filter, unchanged) ----
antipruritic_events <- non_pca_mar_reclassified |>
  filter(given == TRUE,
         tolower(GENERIC_NAME) == "chlorphenamine" |
           naloxone_type %in% c("low_dose_pruritus", "low_dose_pruritus_reclassified")) |>
  select(PAT_ENC_CSN_ID, ADMIN_TIME)

antipruritic_qualifying <- pca_episodes |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end) |>
  left_join(antipruritic_events, by = "PAT_ENC_CSN_ID", relationship = "many-to-many") |>
  filter(!is.na(ADMIN_TIME), ADMIN_TIME >= episode_start, ADMIN_TIME <= episode_end) |>
  distinct(PAT_ENC_CSN_ID, episode_id) |>
  mutate(any_antipruritic_final = TRUE)

# ---- Naloxone reversal (minimal/no cutoff, unchanged) ----
reversal_events <- non_pca_mar_reclassified |>
  filter(given == TRUE, naloxone_type == "reversal_respiratory") |>
  select(PAT_ENC_CSN_ID, ADMIN_TIME)

reversal_qualifying <- pca_episodes |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end) |>
  left_join(reversal_events, by = "PAT_ENC_CSN_ID", relationship = "many-to-many") |>
  filter(!is.na(ADMIN_TIME), ADMIN_TIME >= episode_start, ADMIN_TIME <= episode_end) |>
  distinct(PAT_ENC_CSN_ID, episode_id) |>
  mutate(any_naloxone_reversal_final = TRUE)

# ---- Join everything onto pca_episodes ----
pca_episodes <- pca_episodes |>
  left_join(antiemetic_final_result, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  left_join(antipruritic_qualifying, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  left_join(reversal_qualifying, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  mutate(
    # ============================================================
    # reactive_antiemetic_final — Group C preserved as NA
    #
    # antiemetic_final_result only contains rows for Group A/B
    # (Group C deliberately excluded — see antiemetic_group
    # derivation above). A naive replace_na(..., FALSE) here would
    # silently convert Group C's genuine "excluded, unknown" NA into
    # a false "no reactive antiemetic" — the SAME bug pattern already
    # caught once during testing for confirmed_preop_symptomatic. Group C
    # must stay NA: excluded from the denominator entirely, not
    # counted as a negative.
    # ============================================================
    reactive_antiemetic_final = case_when(
      antiemetic_group == "C" ~ NA,
      TRUE ~ replace_na(reactive_antiemetic_final, FALSE)
    ),
    any_antipruritic_final      = replace_na(any_antipruritic_final, FALSE),
    any_naloxone_reversal_final = replace_na(any_naloxone_reversal_final, FALSE),

    # Composite — Group C's NA correctly propagates through | here
    # too (NA | TRUE = TRUE, NA | FALSE = NA) rather than being
    # silently coerced
    any_side_effect_final = reactive_antiemetic_final | any_antipruritic_final | any_naloxone_reversal_final
  )

cat("Reactive antiemetic (FINAL, two-group design, no ceiling):",
    sum(pca_episodes$reactive_antiemetic_final, na.rm = TRUE), "/",
    sum(!is.na(pca_episodes$reactive_antiemetic_final)),
    sprintf("(%.1f%%)\n", 100*mean(pca_episodes$reactive_antiemetic_final, na.rm = TRUE)))
cat("Excluded from antiemetic denominator (Group C, unreliable timing):",
    sum(is.na(pca_episodes$reactive_antiemetic_final)), "\n")
cat("Antipruritic (final):", sum(pca_episodes$any_antipruritic_final), "\n")
cat("Naloxone reversal (final):", sum(pca_episodes$any_naloxone_reversal_final), "\n")
cat("Composite side effect (final):", sum(pca_episodes$any_side_effect_final, na.rm = TRUE), "\n")

# ============================================================
# POPULATION SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("ANALYSIS POPULATION SUMMARY\n")
cat("============================================================\n")
cat("\nSection 1 - All episodes:", nrow(pca_episodes), "\n")

cat("\nSection 2 - Oxycodone descriptive:\n")
cat("  Layer 1:", pca_episodes |>
      filter(drug == "oxycodone", is_first_drug_episode,
             !concurrent_epidural) |> nrow(), "\n")
cat("  Layer 2:", pca_episodes |>
      filter(drug == "oxycodone", is_first_drug_episode,
             !concurrent_epidural, case_type == "surgical",
             !had_chemotherapy) |> nrow(), "\n")
cat("  Layer 3:", pca_episodes |>
      filter(drug == "oxycodone", is_first_drug_episode,
             !concurrent_epidural, case_type == "surgical",
             !had_chemotherapy, !exclude_minor_procedure_any) |> nrow(), "\n")

cat("\nSection 3 - PSM (morphine vs oxycodone):\n")
cat("  Morphine:", pca_episodes |>
      filter(drug == "morphine", is_first_line, is_first_episode,
             !concurrent_epidural, case_type == "surgical",
             !exclude_minor_procedure_any) |> nrow(), "\n")
cat("  Oxycodone:", pca_episodes |>
      filter(drug == "oxycodone", is_first_line, is_first_episode,
             !concurrent_epidural, case_type == "surgical",
             !exclude_minor_procedure_any) |> nrow(), "\n")

cat("\nSensitivity - include caudals (epidural MAR only group):\n")
cat("  Morphine:", pca_episodes |>
      filter(drug == "morphine", is_first_line, is_first_episode,
             !concurrent_epidural | neuraxial_epidural_mar_only,
             case_type == "surgical",
             !exclude_minor_procedure_any_strict) |> nrow(), "\n")
cat("  Oxycodone:", pca_episodes |>
      filter(drug == "oxycodone", is_first_line, is_first_episode,
             !concurrent_epidural | neuraxial_epidural_mar_only,
             case_type == "surgical",
             !exclude_minor_procedure_any_strict) |> nrow(), "\n")

cat("\nSection 4 - Safety:", nrow(pca_episodes), "\n")
cat("============================================================\n")

# ============================================================
# PAIN SCORE OUTCOMES — PEAK AND TIME-WEIGHTED MEAN
#
# Peak pain score and time-weighted mean pain are computed once here,
# for every episode in the full cohort, and saved as permanent
# columns on pca_episodes_flagged. Any downstream script that needs
# pain data (10_comparative_PSM.R, 14_sensitivity.R, and others)
# reads these columns directly via ordinary left_join/semi_join/
# filter, rather than each recomputing the same metric independently
# for whatever subset it happens to need — a single source of truth
# for pain-score methodology across the project.
#
# time_weighted_mean_pain weights each reading by how long it applied
# (time to the next reading, or to episode_end for the last reading)
# rather than treating every reading as an equal, independent point
# regardless of duration. This avoids a floor-effect problem that a
# plain per-episode median suffers from, since most individual
# readings sit at 0 and would swamp a simple median regardless of
# real differences in the tail. Time before the first reading in an
# episode is excluded — no score exists to assign it, and
# backfilling would be an assumption, not a measurement.
# ============================================================

cat("\nComputing pain score outcomes (peak, time-weighted mean)...\n")

stopifnot(exists("pain_scores_clean"))

pain_episode_windows <- pca_episodes |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end) |>
  distinct()

pain_readings <- pain_episode_windows |>
  left_join(
    pain_scores_clean |> select(PAT_ENC_CSN_ID, score_time, score_value),
    by = "PAT_ENC_CSN_ID", relationship = "many-to-many"
  ) |>
  filter(!is.na(score_time), !is.na(score_value),
         score_time >= episode_start, score_time <= episode_end) |>
  arrange(PAT_ENC_CSN_ID, episode_id, score_time)

peak_pain_all <- pain_readings |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  summarise(peak_pain_score = max(score_value, na.rm = TRUE), .groups = "drop") |>
  mutate(peak_pain_score = if_else(is.infinite(peak_pain_score), NA_real_, peak_pain_score))

time_weighted_pain_all <- pain_readings |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  mutate(
    next_time = lead(score_time, default = NA),
    duration_hrs = if_else(
      is.na(next_time),
      as.numeric(difftime(first(episode_end), score_time, units = "hours")),
      as.numeric(difftime(next_time, score_time, units = "hours"))
    )
  ) |>
  ungroup() |>
  filter(duration_hrs >= 0) |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  summarise(
    time_weighted_mean_pain = sum(score_value * duration_hrs) / sum(duration_hrs),
    .groups = "drop"
  ) |>
  filter(!is.na(time_weighted_mean_pain))

pca_episodes <- pca_episodes |>
  left_join(peak_pain_all, by = c("PAT_ENC_CSN_ID", "episode_id")) |>
  left_join(time_weighted_pain_all, by = c("PAT_ENC_CSN_ID", "episode_id"))

cat("Peak pain score computed for:", sum(!is.na(pca_episodes$peak_pain_score)),
    "/", nrow(pca_episodes), "episodes\n")
cat("Time-weighted mean pain score computed for:",
    sum(!is.na(pca_episodes$time_weighted_mean_pain)),
    "/", nrow(pca_episodes), "episodes\n")

# ============================================================
# SAVE
# ============================================================

pca_episodes_flagged <- pca_episodes

saveRDS(pca_episodes,
        "1 - data/3 - processed_data/pca_episodes_flagged.rds")

cat("\nFlagged episodes saved\n")
cat("Total episodes:", nrow(pca_episodes), "\n")
cat("Columns:", ncol(pca_episodes), "\n")
cat("\nEPIDURAL DETECTION SUMMARY:\n")
cat("  Primary (concurrent_epidural):\n")
cat("    had_lda_epidural + neuraxial_intrathecal_mar\n")
cat("  Sensitivity (neuraxial_epidural_mar_only):\n")
cat("    epidural MAR without LDA — likely caudals\n")
cat("\nFLAGS PRODUCED BY THIS SCRIPT:\n")
cat("  had_lda_epidural, neuraxial_intrathecal_mar\n")
cat("  neuraxial_epidural_mar, neuraxial_epidural_mar_only\n")
cat("  had_intrathecal_opioid, neuraxial_mar_only\n")
cat("  exclude_minor_procedure_any_strict\n")
cat("  is_trauma, medical_specialty_grouped, admission_type_detailed\n")
cat("  likely_nonopioid_nausea_confound\n")
cat("  exclude_pca_before_procedure\n")
cat("\nproc_category / PFMD_PROC_OPCS_CODES now arrive from\n")
cat("2_cleaningscript.R (episode-level attribution) — NOT joined in\n")
cat("this script anymore. Confirm 2_cleaningscript.R's 6A/6B are in\n")
cat("place before running this script, or these columns will be\n")
cat("missing/NA throughout.\n")
