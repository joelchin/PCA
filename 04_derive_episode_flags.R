# ============================================================
# 04 - Derive episode flags and outcome definitions
# ============================================================
# Builds the episode-level flags used by the rest of the pipeline:
# episode sequence, concurrent neuraxial and peripheral blocks,
# intraoperative drug flags, exclusion criteria, admission type, and the
# final antiemetic, antipruritic and naloxone-reversal outcomes, plus
# peak and time-weighted mean pain. Output: pca_episodes_flagged.rds.
#
# Concurrent epidural analgesia (concurrent_epidural) is defined from the
# LDA epidural table only, which records continuous catheters. Spinal
# anaesthesia (intrathecal route in the intraoperative MAR) is not
# counted as concurrent epidural; intrathecal opioid is a propensity
# score covariate. Epidural-route MAR entries with no LDA record are
# probably single-shot caudal blocks and are kept as the separate flag
# neuraxial_epidural_mar_only (sensitivity analysis).
#
# proc_category and PFMD_PROC_OPCS_CODES arrive on pca_episodes already
# attached at episode level from 02_clean_data.R (section 6D), so
# operations_clean is not joined again here.
# ============================================================

library(tidyverse)
library(lubridate)
source(file.path(getOption("pca_code_dir", "."), "drug_lists.R"))

# ============================================================
# LOAD DATA
# ============================================================

pca_episodes             <- readRDS("1 - data/3 - processed_data/pca_episodes.rds")
admissions_clean         <- readRDS("1 - data/3 - processed_data/admissions_clean.rds")
operations_clean         <- readRDS("1 - data/3 - processed_data/operations_clean.rds")
nerve_blocks_clean       <- readRDS("1 - data/3 - processed_data/nerve_blocks_clean.rds")
anaesthesia_events_clean <- readRDS("1 - data/3 - processed_data/anaesthesia_events_clean.rds")
non_pca_mar_clean        <- readRDS("1 - data/3 - processed_data/non_pca_mar_clean.rds")
pain_scores_clean        <- readRDS("1 - data/3 - processed_data/pain_scores_clean.rds")
intraop_mar_clean        <- readRDS("1 - data/3 - processed_data/intraop_mar_clean.rds")
lda_epidural_clean       <- readRDS("1 - data/3 - processed_data/lda_epidural_clean.rds")

cat("All data loaded\n")
cat("Episodes to flag:", nrow(pca_episodes), "\n")
cat("proc_category attached at episode level:",
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
  arrange(PAT_ENC_CSN_ID, episode_start) |>
  group_by(PAT_ENC_CSN_ID) |>
  mutate(
    prior_drug    = lag(drug),
    is_first_line = is.na(prior_drug),
    is_rotation   = !is_first_line & drug != prior_drug
  ) |>
  ungroup()

cat("First episodes (any drug):",    sum(pca_episodes$is_first_episode), "\n")
cat("First line episodes:",          sum(pca_episodes$is_first_line), "\n")
cat("Rotation episodes:",            sum(pca_episodes$is_rotation), "\n")

# ============================================================
# FLAG 2 - Concurrent nerve block and neuraxial flags
# ============================================================
# Sources:
#   1. nerve_blocks table (peripheral blocks; epidurals are documented in
#      the LDA system, so it holds no neuraxial entries)
#   2. intraoperative MAR route: intrathecal (spinal) and epidural
#   3. LDA epidural table (continuous catheter)
#
# concurrent_epidural = had_lda_epidural (source 3 only)
# neuraxial_epidural_mar_only = epidural MAR route, no LDA record and no
#   intrathecal route (probable caudal; sensitivity flag)

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

# Anaesthesia encounter ID, needed to link intraoperative MAR
if (!"anaesthesia_csn_id" %in% names(pca_episodes)) {
  pca_episodes <- pca_episodes |>
    left_join(
      anaesthesia_events_clean |>
        select(PAT_ENC_CSN_ID, anaesthesia_csn_id),
      by = "PAT_ENC_CSN_ID"
    )
  cat("anaesthesia_csn_id added to pca_episodes\n")
}

# Source 2: intraoperative MAR route (intrathecal vs epidural)
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

# Source 3: LDA epidural table (continuous catheter)
lda_epidural_csns <- lda_epidural_clean |>
  distinct(PAT_ENC_CSN_ID) |>
  mutate(had_lda_epidural = TRUE)

pca_episodes <- pca_episodes |>
  left_join(lda_epidural_csns, by = "PAT_ENC_CSN_ID") |>
  mutate(
    had_lda_epidural = replace_na(had_lda_epidural, FALSE),
    
    concurrent_epidural = had_lda_epidural,
    
    neuraxial_epidural_mar_only = neuraxial_epidural_mar &
      !had_lda_epidural &
      !neuraxial_intrathecal_mar
  )

cat("Episodes with any peripheral block:",
    sum(pca_episodes$had_any_block), "\n")
cat("Episodes with neuraxial (nerve blocks table):",
    sum(pca_episodes$had_neuraxial), "\n")
cat("Episodes with intrathecal MAR (spinal):",
    sum(pca_episodes$neuraxial_intrathecal_mar), "\n")
cat("Episodes with epidural MAR route (includes caudals):",
    sum(pca_episodes$neuraxial_epidural_mar), "\n")
cat("Episodes with LDA epidural (continuous catheter):",
    sum(pca_episodes$had_lda_epidural), "\n")
cat("Episodes with concurrent epidural (LDA catheter):",
    sum(pca_episodes$concurrent_epidural), "\n")
cat("Episodes with epidural MAR route only, no LDA (probable caudals):",
    sum(pca_episodes$neuraxial_epidural_mar_only), "\n")
cat("Episodes with peripheral block only:",
    sum(pca_episodes$had_peripheral &
          !pca_episodes$concurrent_epidural), "\n")
cat("\nConcurrent epidural by drug:\n")
print(table(pca_episodes$drug, pca_episodes$concurrent_epidural))

# ============================================================
# FLAG 3 - Intraoperative drug flags
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
# FLAG 4 - Predefined exclusion flags
# ============================================================

cat("\nAdding the predefined exclusion flags...\n")

# all_opcs_excluded: TRUE only if every code on every operation row of
# the admission is on opcs_exclude (drug_lists.R). All codes on a row are
# checked, not just the first, so an admission with a substantive
# companion code is not wrongly excluded.

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
# ============================================================
# Flags episodes that start more than 2 hours BEFORE the admission's
# anaesthesia stop (implausible PCA-to-procedure timing); the opposite
# direction to exclude_late_pca_start. The 2-hour buffer allows for PCA
# being programmed in theatre before the end of the anaesthetic.
#
# Applied to single-operation admissions only: the extract holds one
# anaesthetic record per admission, so with several operations the
# recorded anaesthesia_stop cannot be attributed to a specific
# operation. Multi-operation admissions are excluded anyway by
# has_multiple_operations.
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

cat("\nPredefined exclusion flag summary:\n")
cat("Non-substantive OPCS only:",
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
cat("(exclude_minor_procedure_any is finalised after the demographics join, below)\n")

# ============================================================
# has_reliable_anaesthesia_timing
# ============================================================
# The extract holds one anaesthetic record per admission, so in
# multi-operation admissions anaesthesia_stop may belong to a different
# operation from the one an episode followed. The flag is TRUE for
# single-operation admissions, and for multi-operation admissions only
# when the episode's attributed operation (attributed_op_order, from
# 02_clean_data.R section 6D) is the operation whose anaesthesia record
# was matched by ID. It is used to place episodes in antiemetic groups
# A/B/C below.
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

cat("\nhas_reliable_anaesthesia_timing by drug:\n")
print(table(pca_episodes$drug, pca_episodes$has_reliable_anaesthesia_timing))

# ============================================================
# FLAG 5 - Demographics
# ============================================================
# Only the admissions demographics are joined; proc_category already
# comes from 02_clean_data.R at episode level.
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

cat("Demographics joined\n")
cat("proc_category present:", "proc_category" %in% names(pca_episodes), "\n")
cat("PFMD_PROC_OPCS_CODES present:", "PFMD_PROC_OPCS_CODES" %in% names(pca_episodes), "\n")

# ============================================================
# Exclusion flags that depend on case_type
# ============================================================
# These need case_type, which exists on pca_episodes only after the
# demographics join above.

# Episodes with no anaesthesia_stop (so exclude_pca_before_procedure
# cannot fire) whose only operation has no preceding date
pca_episodes <- pca_episodes |>
  mutate(
    exclude_pca_before_procedure_opdate_fallback =
      !has_multiple_operations &
      is.na(anaesthesia_stop) &
      case_type == "surgical" &
      is.na(proc_category)
  )

cat("\nPCA before procedure, OpDate fallback (anaesthesia_stop NA):",
    sum(pca_episodes$exclude_pca_before_procedure_opdate_fallback, na.rm = TRUE), "\n")

# Surgical episodes with no procedure category cannot enter a
# procedure-matched comparison
pca_episodes <- pca_episodes |>
  mutate(
    exclude_missing_proc_category = case_type == "surgical" & is.na(proc_category)
  )

cat("Missing proc_category:",
    sum(pca_episodes$exclude_missing_proc_category, na.rm = TRUE), "\n")

# exclude_minor_procedure_any and exclude_minor_procedure_any_strict are
# finalised here because every contributing flag now exists.
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

cat("\nAny predefined exclusion:",
    sum(pca_episodes$exclude_minor_procedure_any), "\n")
cat("Any predefined exclusion (strict version):",
    sum(pca_episodes$exclude_minor_procedure_any_strict), "\n")

cat("\nSurgical episodes with NA proc_category not caught by\n")
cat("exclude_minor_procedure_any (should be 0):\n")
cat(sum(pca_episodes$case_type == "surgical" &
          is.na(pca_episodes$proc_category) &
          !pca_episodes$exclude_minor_procedure_any, na.rm = TRUE), "\n")

# ============================================================
# case_type_refined
# ============================================================
# case_type (from 02_clean_data.R) is surgical if any operation is
# recorded for the admission, whatever its size. case_type_refined
# reclassifies a surgical episode as medical when any of the following
# holds:
#   - its attributed operation is entirely minor (row_fully_minor)
#   - every operation in the admission is minor (all_opcs_excluded)
#   - it is shorter than 20 minutes of anaesthesia
#   - no operation precedes it (proc_category is NA)
#   - it starts more than 2 hours before anaesthesia stop
#     (exclude_pca_before_procedure)
# It is a separate field so that case_type, which the population
# definitions use, is unchanged.
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
# FLAG 6 - Era
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
# FLAG 7 - Admission type
# ============================================================

pca_episodes <- pca_episodes |>
  mutate(
    admission_type = case_when(
      str_detect(coalesce(PFMD_PROC_OPCS_CODES, ""),
                 "^W1[0-9]|^W2[0-9]|^W30") ~ "Trauma",
      proc_category == "Orthopaedic/spine" ~ "Elective orthopaedic",
      case_type == "surgical" ~ "Other surgical",
      TRUE ~ "Medical"
    )
  )

cat("\nAdmission type by drug:\n")
print(table(pca_episodes$admission_type, pca_episodes$drug))

# ============================================================
# ============================================================
# FLAG 8 - Final side-effect outcomes
# ============================================================
# Single definition of the antiemetic, antipruritic and naloxone-reversal
# outcomes, read by the downstream analysis scripts.
# ============================================================

cat("\nAdding final side-effect classifications...\n")

# naloxone_reclassifications.rds is not produced by this pipeline. It
# holds the manual chart-review reclassification of naloxone
# administrations (pruritus dose vs respiratory reversal) and must exist
# in the processed-data folder before this script is run.
naloxone_reclassifications <- readRDS("1 - data/3 - processed_data/naloxone_reclassifications.rds")

non_pca_mar_reclassified <- non_pca_mar_clean |>
  left_join(
    naloxone_reclassifications |> select(PAT_ENC_CSN_ID, ADMIN_TIME, reclassified_to),
    by = c("PAT_ENC_CSN_ID", "ADMIN_TIME")
  ) |>
  mutate(naloxone_type = if_else(!is.na(reclassified_to) & tolower(GENERIC_NAME) %in% "naloxone",
                                 reclassified_to, naloxone_type)) |>
  select(-reclassified_to)

# ============================================================
# Reactive antiemetic outcome
# ============================================================
# Episodes are split into three groups by the reliability of the
# anaesthesia timing (antiemetic_group):
#
#   A - no anaesthesia_stop recorded. Reactive if an antiemetic was
#       given during the episode and none in the 7 days before
#       episode_start.
#   B - anaesthesia_stop present and reliable for this episode. Reactive
#       if an antiemetic was given more than 12 hours after
#       anaesthesia_stop (and within the episode) and none in the 7 days
#       before anaesthesia_start. Anchoring to anaesthesia_start avoids
#       counting routine intraoperative prophylaxis as pre-existing
#       symptoms.
#   C - anaesthesia_stop present but not reliably linked to this episode
#       (has_reliable_anaesthesia_timing is FALSE). Left as NA and
#       excluded from the denominator, because post-anaesthesia
#       antiemetic timing cannot be interpreted without a trustworthy
#       anchor.
#
# No cap on the proportion of shifts with antiemetic is applied.
# Antiemetic drugs: antiemetic_drugs in drug_lists.R (dexamethasone
# excluded).
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
    antiemetic_group = case_when(
      is.na(anaesthesia_stop) ~ "A",
      !has_reliable_anaesthesia_timing ~ "C",
      TRUE ~ "B"
    )
  )

cat("\nantiemetic_group split:\n")
print(table(pca_episodes$antiemetic_group))
cat("Group C (unreliable anaesthesia timing, excluded from the antiemetic denominator):",
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

cat("\n--- Reactive antiemetic, groups A and B ---\n")
cat("Group A:", sum(group_a_episodes$PAT_ENC_CSN_ID %in% group_a_result$PAT_ENC_CSN_ID), "episodes,",
    sum(group_a_result$reactive_antiemetic_final), "reactive (",
    round(100*mean(group_a_result$reactive_antiemetic_final), 1), "%)\n")
cat("Group B:", nrow(group_b_result), "episodes,",
    sum(group_b_result$reactive_antiemetic_final), "reactive (",
    round(100*mean(group_b_result$reactive_antiemetic_final), 1), "%)\n\n")

# ---- Antipruritic (chlorphenamine or low-dose naloxone; no time filter) ----
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

# ---- Naloxone reversal (400 microgram formulation) ----
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
    # Group C stays NA (excluded from the denominator), not FALSE
    reactive_antiemetic_final = case_when(
      antiemetic_group == "C" ~ NA,
      TRUE ~ replace_na(reactive_antiemetic_final, FALSE)
    ),
    any_antipruritic_final      = replace_na(any_antipruritic_final, FALSE),
    any_naloxone_reversal_final = replace_na(any_naloxone_reversal_final, FALSE),
    
    # Composite of the three outcomes (NA for group C unless another
    # component is TRUE); the manuscript reports the components separately
    any_side_effect_final = reactive_antiemetic_final | any_antipruritic_final | any_naloxone_reversal_final
  )

cat("Reactive antiemetic:",
    sum(pca_episodes$reactive_antiemetic_final, na.rm = TRUE), "/",
    sum(!is.na(pca_episodes$reactive_antiemetic_final)),
    sprintf("(%.1f%%)\n", 100*mean(pca_episodes$reactive_antiemetic_final, na.rm = TRUE)))
cat("Excluded from antiemetic denominator (group C):",
    sum(is.na(pca_episodes$reactive_antiemetic_final)), "\n")
cat("Antipruritic:", sum(pca_episodes$any_antipruritic_final), "\n")
cat("Naloxone reversal:", sum(pca_episodes$any_naloxone_reversal_final), "\n")
cat("Composite side effect:", sum(pca_episodes$any_side_effect_final, na.rm = TRUE), "\n")

# ============================================================
# Pain score outcomes: peak and time-weighted mean
# ============================================================
# Computed once for every episode and saved on pca_episodes_flagged, so
# that downstream scripts read the same values. Readings within
# [episode_start, episode_end] are used. The time-weighted mean weights
# each reading by the time until the next reading (or episode_end for
# the last reading) and excludes time before the first reading, for
# which no score exists. This avoids the floor effect of a plain median,
# since most individual readings are zero.
# ============================================================

cat("\nComputing pain score outcomes (peak, time-weighted mean)...\n")

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

saveRDS(pca_episodes,
        "1 - data/3 - processed_data/pca_episodes_flagged.rds")

cat("\nFlagged episodes saved\n")
cat("Total episodes:", nrow(pca_episodes), "\n")
cat("Columns:", ncol(pca_episodes), "\n")
