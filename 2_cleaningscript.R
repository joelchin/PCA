# ============================================================
# 2 — Clean raw tables and build PCA/NCA episodes
# Paediatric PCA/NCA study — Addenbrooke's Hospital (CUH NHS FT)
# Author: J Chin
# ============================================================
#
# Purpose
#   Cleans each raw EPIC table loaded by 1_load_raw_data.R, builds the episode
#   table (pca_episodes: one row per PCA/NCA episode), attributes a procedure
#   category to every episode, and saves all cleaned tables.
#
# Inputs
#   In session (from script 1): admissions, pca_mar, non_pca_mar, nerve_blocks,
#     operations, anaesthesia_events, anaesthesia_ax_sdes, lda_* tables,
#     pain_scores, intraop_mar.
#   Files: drug_lists.R (drug name lists; edit them there, not here) and
#     1 - data/3 - processed_data/1_shifts.rds (12-hour shift grid; created if
#     missing, see README).
#
# Outputs (1 - data/3 - processed_data/, .rds)
#   admissions_clean, pca_mar_clean, pca_episodes, non_pca_mar_clean,
#   nerve_blocks_clean, operations_clean, anaesthesia_events_clean,
#   anaesthesia_ax_clean, lda_epidural_clean, lda_urethral_clean,
#   lda_epidural_assess_clean, pain_scores_clean, intraop_mar_clean.
#
# Key definitions
#   - Identifiers: MRN and name are dropped from admissions on load. Date of
#     birth is used only to derive age and is then dropped, so exact dates of
#     birth are not carried into any saved file.
#   - Shift grid: 12-hour blocks anchored at 08:00 from 2014-01-01. If any
#     timestamp falls after the last block the grid is extended to the end of
#     the current calendar year, so shift matching cannot silently return NA.
#     This is the only step that writes back to a reference file.
#   - Procedure category (Sections 6A-6D): derived from the first *substantive*
#     OPCS code on each operation record (not just the first code listed) and
#     attributed to each PCA episode from the operation nearest in time before
#     the episode began, rather than the admission's first operation. Category
#     names are restricted to the original 16 categories.
#   - Runtime: about 30 minutes (shift matching runs in parallel via furrr).
#
# Run order: 1 -> 2 -> 3 -> 4 -> ... (see RUN_ALL_1_to_16.R)
# ============================================================

library(tidyverse)
library(lubridate)
library(furrr)
library(future)
library(tictoc)

# Load drug lists
source("drug_lists.R")

# ============================================================
# PARALLEL PROCESSING SETUP
# ============================================================

available_cores <- future::availableCores()
options(parallelly.fork.enable = TRUE)
plan(multisession, workers = available_cores - 1)

cat("Using", available_cores - 1, "cores for parallel processing\n")

# ============================================================
# SHIFT MATCHING FUNCTION
# ============================================================

match_date_to_shift_number <- function(date_time, shifts) {
  suppressPackageStartupMessages(library(lubridate))
  date_time <- as.POSIXct(date_time, format = "%Y-%m-%d %H:%M:%S")
  matched_shift <- shifts |>
    filter(date_time >= start_time & date_time <= end_time)
  if (nrow(matched_shift) > 0) {
    return(matched_shift$shift_number[1])
  } else {
    return(NA)
  }
}

# ============================================================
# LOAD SHIFTS TABLE — AUTO-EXTEND IF DATA CEILING HAS MOVED PAST IT
#
# 1_shifts.rds is a static, pre-generated 12h-block grid anchored at
# 2014-01-01 08:00 that does not auto-extend on its own. Any
# timestamp falling after its last block would otherwise be silently
# mapped to NA by match_date_to_shift_number() rather than raising an
# error, and any PCA episode with an NA shift would then silently
# drop out of episode_shifts/episode_metrics downstream while
# remaining present in pca_episodes — a population mismatch between
# the two objects that's easy to miss.
#
# To avoid this, the table is checked against the current data's
# ceiling on every run and extended forward to the end of the current
# calendar year whenever needed, so a new data refresh can never
# silently outrun the table's coverage.
# ============================================================

shifts_file <- "1 - data/3 - processed_data/1_shifts.rds"
if (file.exists(shifts_file)) {
  shifts <- readRDS(shifts_file)
} else {
  # No shift grid on disk: start a new one. Blocks are 12 hours long, anchored at
  # 2014-01-01 08:00 (clock time as returned by read_excel, labelled UTC) and
  # numbered from 1. The extension step below fills it forward to the end of the
  # current year. An existing grid is never replaced by this branch.
  dir.create("1 - data/3 - processed_data", recursive = TRUE, showWarnings = FALSE)
  shifts <- tibble(
    shift_number = 1L,
    start_time   = as.POSIXct("2014-01-01 08:00:00", tz = "UTC"),
    end_time     = as.POSIXct("2014-01-01 08:00:00", tz = "UTC") + hours(12) - seconds(1)
  )
  cat("1_shifts.rds not found — starting a new shift grid at 2014-01-01 08:00\n")
}
cat("Shifts loaded:", nrow(shifts), "shifts\n")

target_ceiling <- as.POSIXct(
  paste0(year(Sys.Date()), "-12-31 23:59:59"),
  tz = tz(shifts$end_time)
)

if (max(shifts$end_time) < target_ceiling) {

  last_shift_number <- max(shifts$shift_number)
  last_end_time     <- max(shifts$end_time)

  new_start_times <- seq(
    from = last_end_time + 1,
    to   = target_ceiling,
    by   = "12 hours"
  )

  new_shifts <- tibble(
    shift_number = seq(last_shift_number + 1,
                       last_shift_number + length(new_start_times)),
    start_time = new_start_times,
    end_time   = new_start_times + hours(12) - seconds(1)
  )

  shifts <- bind_rows(shifts, new_shifts)

  saveRDS(shifts, "1 - data/3 - processed_data/1_shifts.rds")

  cat("Shifts table extended:", nrow(new_shifts), "new blocks added",
      "(previous coverage capped at", format(last_end_time), ")\n")
  cat("   New coverage ceiling:", format(max(shifts$end_time)), "\n")

} else {
  cat("Shifts table coverage already extends to",
      format(max(shifts$end_time)), "— no extension needed\n")
}

tic()

# ============================================================
# 1. ADMISSIONS
# ============================================================

admissions_clean <- admissions |>
  select(-PAT_MRN_ID, -PAT_NAME, -INPATIENT_DATA_ID) |>
  mutate(
    ADM_DATE       = as.POSIXct(ADM_DATE),
    DISCHARGE_DATE = as.POSIXct(DISCHARGE_DATE),
    DOB            = as.POSIXct(DOB),
    WEIGHT_KG      = as.numeric(WEIGHT_KG),
    age_days       = as.numeric(
      difftime(ADM_DATE, DOB, units = "days")
    ),
    age_group = case_when(
      age_days < 28       ~ "Neonate (<28 days)",
      AGE_IN_MONTHS < 12  ~ "Infant (28 days-1yr)",
      AGE < 6             ~ "Young child (1-5yr)",
      AGE < 12            ~ "Child (6-11yr)",
      AGE < 16            ~ "Adolescent (12-15yr)",
      TRUE                ~ NA_character_
    ),
    LOS_days = as.numeric(
      difftime(DISCHARGE_DATE, ADM_DATE, units = "days")
    )
  ) |>
  select(-DOB) |>   # DOB is only needed to derive age_days above; not carried forward
  distinct(PAT_ENC_CSN_ID, .keep_all = TRUE) |>
  mutate(
    admission_shift = future_map_dbl(
      ADM_DATE,
      ~ match_date_to_shift_number(.x, shifts),
      .progress = TRUE
    ),
    discharge_shift = future_map_dbl(
      DISCHARGE_DATE,
      ~ match_date_to_shift_number(.x, shifts),
      .progress = TRUE
    )
  )

cat("Admissions clean:", nrow(admissions_clean), "rows\n")
cat("\nAge group distribution:\n")
print(table(admissions_clean$age_group, useNA = "always"))

# ============================================================
# 2. PCA MAR
# ============================================================
pca_mar_clean <- pca_mar |>
  mutate(
    TAKEN_TIME       = as.POSIXct(TAKEN_TIME),
    Discontinue_Time = as.POSIXct(Discontinue_Time),

    drug = case_when(
      str_detect(toupper(MEDICATION), "MORPHINE/KETAMINE") ~ "morphine+ketamine",
      str_detect(toupper(MEDICATION), "OXYCODONE")         ~ "oxycodone",
      str_detect(toupper(MEDICATION), "FENTANYL")          ~ "fentanyl",
      str_detect(toupper(MEDICATION), "MORPHINE")          ~ "morphine",
      TRUE                                                  ~ "other"
    ),

    morphine_mg = as.numeric(
      str_extract(MED_ORDER,
                  "(?i)morphine(?: sulfate)? ([0-9.]+) mg", group = 1)
    ),
    ketamine_mg = as.numeric(
      str_extract(MED_ORDER,
                  "(?i)ketamine ([0-9.]+) mg", group = 1)
    ),
    oxycodone_mg = as.numeric(
      str_extract(MED_ORDER,
                  "(?i)oxycodone(?: hydrochloride)? ([0-9.]+) mg",
                  group = 1)
    ),
    diamorphine_mg = as.numeric(
      str_extract(MED_ORDER,
                  "(?i)diamorphine(?: hydrochloride)? ([0-9.]+) mg",
                  group = 1)
    ),
    fentanyl_mcg = as.numeric(
      str_remove_all(
        str_extract(MED_ORDER,
                    "(?i)fentanyl ([0-9,]+) mcg", group = 1),
        ","
      )
    ),

    weight_based = str_detect(MED_ORDER, "(?i)mg/kg|mcg/kg"),

    lockout_mins = as.numeric(
      str_extract(PROGRAMME, "(?i)(?<=lockout interval: )\\d+")
    ),

    mode = case_when(
      lockout_mins <= 7   ~ "PCA",
      lockout_mins >= 10  ~ "NCA",
      TRUE                ~ "unknown"
    ),

    bolus_dose_ml = as.numeric(
      str_extract(PROGRAMME, "(?i)(?<=bolus dose: )\\d+\\.?\\d*")
    ),

    # Continuous background infusion rate (mL/hr), parsed from the
    # PROGRAMME text. Covers 98.5% of rows; the raw INFUSION_RATE
    # column only covers 13.9% and is too sparse to use directly.
    background_rate_ml = as.numeric(
      str_extract(PROGRAMME, "(?i)(?<=continuous infusion rate:\\s)[0-9.]+")
    ),

    # Loading dose (mL), parsed from PROGRAMME. "None" (most orders)
    # correctly yields NA rather than 0 — any downstream
    # had_loading_dose flag should be based on !is.na(), not > 0.
    loading_dose_ml = as.numeric(
      str_extract(PROGRAMME, "(?i)(?<=loading dose:\\s)[0-9.]+")
    ),

    entry_duration_hrs = as.numeric(
      difftime(Discontinue_Time, TAKEN_TIME, units = "hours")
    ),

    duration_flag = case_when(
      entry_duration_hrs < 0  ~ "negative - check",
      TRUE                    ~ "ok"
    )
  )

# ============================================================
# CONCENTRATION EXTRACTION (mg per mL)
#
# Needed by 7_dosing.R. Two more direct sources were considered and
# rejected: PROGRAMME's "Each Xml bolus should deliver Y" statement
# covers only 2/20,470 rows, and DOSE/DOSE_UNIT is ~99.7% missing
# despite DOSE_UNIT occasionally containing "mg/mL". Concentration is
# instead derived from MED_ORDER: oxycodone_mg/morphine_mg already
# hold the pharmacy-calculated total mg dissolved in that syringe
# (weight-adjusted at preparation time), so dividing by diluent
# volume gives concentration directly. weight_based orders ("X mg/kg"
# text) are excluded, since there's no fixed total mg to divide.
# ============================================================

cat("\nExtracting drug concentration (mg/mL) from MED_ORDER...\n")

pca_mar_clean <- pca_mar_clean |>
  mutate(
    diluent_volume_ml = as.numeric(
      str_extract(MED_ORDER, "(?i)(?<=in )[\\s\\S]*?([0-9.]+)\\s*mL", group = 1)
    ),
    morphine_conc_mg_per_ml = if_else(
      !weight_based & drug == "morphine" & !is.na(morphine_mg) & !is.na(diluent_volume_ml) & diluent_volume_ml > 0,
      morphine_mg / diluent_volume_ml,
      NA_real_
    ),
    oxycodone_conc_mg_per_ml = if_else(
      !weight_based & drug == "oxycodone" & !is.na(oxycodone_mg) & !is.na(diluent_volume_ml) & diluent_volume_ml > 0,
      oxycodone_mg / diluent_volume_ml,
      NA_real_
    )
  )

cat("Morphine concentration extracted:",
    sum(!is.na(pca_mar_clean$morphine_conc_mg_per_ml)), "/",
    sum(pca_mar_clean$drug == "morphine"), "\n")
cat("Oxycodone concentration extracted:",
    sum(!is.na(pca_mar_clean$oxycodone_conc_mg_per_ml)), "/",
    sum(pca_mar_clean$drug == "oxycodone"), "\n")
cat("Weight-based orders excluded (no fixed concentration):",
    sum(pca_mar_clean$weight_based, na.rm = TRUE), "\n\n")

cat("PCA MAR clean:", nrow(pca_mar_clean), "rows\n")
cat("\nDrug distribution:\n")
print(table(pca_mar_clean$drug))
cat("\nMode distribution:\n")
print(table(pca_mar_clean$mode))
cat("\nMAR Action distribution:\n")
print(table(pca_mar_clean$MAR_ACTION))
cat("\nDuration flags:\n")
print(table(pca_mar_clean$duration_flag))
cat("\nBackground rate extracted:",
    sum(!is.na(pca_mar_clean$background_rate_ml)), "/",
    nrow(pca_mar_clean),
    sprintf("(%.1f%%)\n", 100*mean(!is.na(pca_mar_clean$background_rate_ml))))
cat("Loading dose specified:",
    sum(!is.na(pca_mar_clean$loading_dose_ml)), "/",
    nrow(pca_mar_clean),
    sprintf("(%.1f%%)\n", 100*mean(!is.na(pca_mar_clean$loading_dose_ml))))
# ============================================================
# 3. PCA EPISODES
# ============================================================

pca_episodes <- pca_mar_clean |>

  group_by(PAT_ENC_CSN_ID, ORDER_MED_ID) |>
  summarise(
    drug             = first(drug),
    mode             = first(mode),
    lockout_mins     = first(lockout_mins),
    bolus_dose_ml    = first(bolus_dose_ml),
    order_start      = min(TAKEN_TIME, na.rm = TRUE),
    order_end        = max(Discontinue_Time, na.rm = TRUE),
    n_syringes       = n(),
    has_missing_stop = any(is.na(Discontinue_Time)),
    .groups = "drop"
  ) |>

  arrange(PAT_ENC_CSN_ID, order_start) |>
  group_by(PAT_ENC_CSN_ID) |>
  mutate(
    gap_hrs     = as.numeric(
      difftime(order_start, lag(order_end), units = "hours")
    ),
    drug_change = drug != lag(drug, default = first(drug)),
    new_episode = is.na(gap_hrs) | gap_hrs > 12 | drug_change,
    episode_id  = cumsum(new_episode)
  ) |>

  group_by(PAT_ENC_CSN_ID, episode_id) |>
  summarise(
    drug                 = first(drug),
    mode                 = first(mode),
    lockout_mins         = first(lockout_mins),
    bolus_dose_ml        = first(bolus_dose_ml),
    episode_start        = min(order_start, na.rm = TRUE),
    episode_end          = max(order_end, na.rm = TRUE),
    episode_duration_hrs = as.numeric(
      difftime(
        max(order_end, na.rm = TRUE),
        min(order_start, na.rm = TRUE),
        units = "hours"
      )
    ),
    n_orders             = n(),
    n_syringes           = sum(n_syringes),
    has_missing_stop     = any(has_missing_stop),
    drug_switch          = n_distinct(drug) > 1,
    .groups = "drop"
  ) |>

  mutate(
    episode_end = if_else(
      is.infinite(episode_end),
      as.POSIXct(NA),
      episode_end
    ),
    episode_duration_hrs = if_else(
      is.infinite(episode_duration_hrs) | episode_duration_hrs < 0,
      NA_real_,
      episode_duration_hrs
    )
  ) |>

  mutate(
    episode_start_shift = future_map_dbl(
      episode_start,
      ~ match_date_to_shift_number(.x, shifts),
      .progress = TRUE
    ),
    episode_end_shift = future_map_dbl(
      episode_end,
      ~ match_date_to_shift_number(.x, shifts),
      .progress = TRUE
    )
  )

cat("PCA episodes:", nrow(pca_episodes), "\n")
cat("Episodes with missing stop:", sum(pca_episodes$has_missing_stop), "\n")
cat("Median syringes per episode:",
    median(pca_episodes$n_syringes), "\n")
cat("Median orders per episode:",
    median(pca_episodes$n_orders), "\n")
cat("\nEpisode duration summary (hrs):\n")
print(summary(pca_episodes$episode_duration_hrs))
cat("\nDrug distribution:\n")
print(table(pca_episodes$drug))
cat("\nMode distribution:\n")
print(table(pca_episodes$mode))

# ============================================================
# 4. NON-PCA MAR
# ============================================================

non_pca_mar_clean <- non_pca_mar |>
  mutate(
    ADMIN_TIME = as.POSIXct(ADMIN_TIME),

    naloxone_type = case_when(
      tolower(GENERIC_NAME) == "naloxone" &
        str_detect(MEDICATION, "40MICR")  ~ "low_dose_pruritus",
      tolower(GENERIC_NAME) == "naloxone" &
        str_detect(MEDICATION, "400MIC")  ~ "reversal_respiratory",
      tolower(GENERIC_NAME) == "naloxone" ~ "naloxone_unknown",
      TRUE                                ~ NA_character_
    ),

    drug_category = case_when(
      tolower(GENERIC_NAME) %in% tolower(analgesic_drugs)   ~ "analgesic",
      tolower(GENERIC_NAME) %in% tolower(side_effect_drugs) ~ "side_effect_marker",
      tolower(GENERIC_NAME) %in% tolower(laxative_drugs)    ~ "laxative",
      tolower(GENERIC_NAME) %in% tolower(chemo_drugs)       ~ "chemotherapy",
      TRUE                                                   ~ "other"
    ),

    analgesic_class = case_when(
      tolower(GENERIC_NAME) == "paracetamol"                      ~ "paracetamol",
      tolower(GENERIC_NAME) %in% tolower(nsaid_drugs)             ~ "NSAID",
      tolower(GENERIC_NAME) %in% tolower(opioid_drugs)            ~ "opioid",
      tolower(GENERIC_NAME) %in% tolower(neuropathic_drugs)       ~ "neuropathic",
      tolower(GENERIC_NAME) %in% tolower(adjuvant_drugs)          ~ "adjuvant",
      TRUE                                                         ~ NA_character_
    ),

    given = case_when(
      str_detect(tolower(MAR_ACTION), "given")                ~ TRUE,
      str_detect(tolower(MAR_ACTION), "patient-administered") ~ TRUE,
      str_detect(tolower(MAR_ACTION),
                 "not given|held|missed")                     ~ FALSE,
      TRUE                                                     ~ NA
    )
  ) |>
  mutate(
    med_admin_shift = future_map_dbl(
      ADMIN_TIME,
      ~ match_date_to_shift_number(.x, shifts),
      .progress = TRUE
    )
  )

cat("Non-PCA MAR clean:", nrow(non_pca_mar_clean), "rows\n")
cat("\nDrug category distribution:\n")
print(table(non_pca_mar_clean$drug_category))
cat("\nAnalgesic class distribution:\n")
print(table(non_pca_mar_clean$analgesic_class))
cat("\nNaloxone type:\n")
print(table(non_pca_mar_clean$naloxone_type, useNA = "always"))
cat("\nGiven distribution:\n")
print(table(non_pca_mar_clean$given, useNA = "always"))

# ============================================================
# 5. NERVE BLOCKS
# ============================================================

nerve_blocks_clean <- nerve_blocks |>
  mutate(
    block_start_time = as.POSIXct(block_start_time),
    block_end_time   = as.POSIXct(block_end_time),
    block_category = case_when(
      is.na(block_type) | block_type == ""              ~ "unspecified",
      str_detect(tolower(block_type),
                 paste(neuraxial_blocks,
                       collapse = "|"))                 ~ "neuraxial",
      str_detect(tolower(block_type),
                 paste(truncal_blocks,
                       collapse = "|"))                 ~ "truncal",
      str_detect(tolower(block_type),
                 paste(lower_limb_blocks,
                       collapse = "|"))                 ~ "lower_limb",
      str_detect(tolower(block_type),
                 paste(upper_limb_blocks,
                       collapse = "|"))                 ~ "upper_limb",
      TRUE                                               ~ "other"
    ),
    technique = case_when(
      str_detect(tolower(injection_technique),
                 "catheter|continuous")                 ~ "catheter",
      str_detect(tolower(injection_technique),
                 "single")                              ~ "single_shot",
      TRUE                                               ~ "unknown"
    )
  )

cat("Nerve blocks clean:", nrow(nerve_blocks_clean), "rows\n")
cat("\nBlock category:\n")
print(table(nerve_blocks_clean$block_category))
cat("\nTechnique:\n")
print(table(nerve_blocks_clean$technique))

# ============================================================
# 6. OPERATIONS
# ============================================================

operations_clean <- operations |>
  mutate(
    primary_procedure = OP_ORDER == 1,
    primary_opcs  = str_extract(
      PFMD_PROC_OPCS_CODES, "^[A-Z]\\d+\\.?\\d*"
    ),
    opcs_chapter  = str_extract(primary_opcs, "^[A-Z]")
  )

cat("Operations clean (base extraction):", nrow(operations_clean), "rows\n")

# ============================================================
# 6A. primary_opcs_substantive — RESOLVING BUNDLED OPCS CODES
#
# Operation records can bundle more than one OPCS code under a single
# row (25.7% of rows, 1,200/4,663) when multiple procedures are
# performed together. A simple first-code extraction (primary_opcs
# above) is a positional choice, not a clinical one, and can select
# an incidental code ahead of the genuinely substantive procedure —
# for example a row coded "N13.4 \r\n T30.9", where T30.9 (exploratory
# laparotomy) is the substantive operation but sits second.
#
# primary_opcs_substantive addresses this directly: it takes the
# first code in the row that is not on opcs_exclude (the
# minor-procedure list, sourced from drug_lists.R so this script and
# 4_flags.R share one definition), falling back to primary_opcs only
# when every code in the row is minor — in that case the row is
# non-substantive regardless of which code is picked.
#
# primary_opcs is retained unchanged above for reference;
# primary_opcs_substantive is what proc_category and the 6D
# episode-attribution step are built from.
# ============================================================

# opcs_treat_as_minor_for_category — codes to skip when searching for
# a substantive companion code, used only for proc_category
# selection. Deliberately separate from opcs_exclude (which drives
# PSM eligibility) — these codes are judged incidental on clinical
# grounds (standalone diagnostic/minor procedures unlikely to
# independently explain multi-day PCA use), using the same standard
# already applied when opcs_exclude itself was built.
#
# Defined here as a standalone vector rather than inline inside the
# mutate() call below: opcs_exclude has ~31 elements, and dplyr
# recycles a vector of that length down the rows of operations_clean
# rather than passing each row the full list, so an inline %in% check
# inside map2_chr would silently compare against one arbitrary code
# per row instead of the full minor-code list.
opcs_treat_as_minor_for_category <- c(
  opcs_exclude,
  "N13.4",   # biopsy testis — oncology-adjacent, often bundled with a genuinely substantive companion code

  # ---- additional minor/incidental codes ----
  "A76.5",   # splanchnic sympathetic nerve block — an analgesic
             # technique, not a driving procedure
  "H02.4",   # laparoscopic incidental appendicectomy — self-labelled
  "H44.6",   # EUA rectum, manual evacuation — minor, brief
  "M13.1",   # needle biopsy kidney, US-guided — standalone diagnostic
  "T87.2",   # excision biopsy cervical lymph node — oncology-diagnostic,
             # same pattern as N13.4 (often bundled with Hickman/BMA)
  "S15.1",   # biopsy head lesion — standalone diagnostic
  "S15.2",   # skin biopsy — standalone diagnostic
  "P09.2",   # drainage vulval lesion — minor
  "P31.2",   # drainage pouch of Douglas — minor
  "P27.3",   # cystovaginoscopy (+biopsy abdominal mass) — diagnostic
  "F36.2",   # biopsy tonsil — standalone diagnostic

  # ---- tube/stoma-maintenance cluster (Abdominal/GI) ----
  # device servicing, not new operative trauma — closure/exchange of
  # an existing gastrostomy/jejunostomy or reduction of a stoma
  # prolapse, not a fresh source of pain
  "G34.4",   # closure gastrostomy
  "G34.5",   # removal GT tube + GT insertion OGD
  "G60.3",   # closure jejunostomy
  "G75.5",   # reduction ileostomy prolapse

  # ---- identified via systematic review of oncology-patient episodes ----
  "X55.8",   # examination under anaesthetic alone — correctly
             # non-substantive but originally missing from this list,
             # so row_fully_minor/case_type_refined never picked it up
  "A52.5",   # thoracic epidural — an anaesthetic technique, not an
             # operation; miscoded as if it were a procedure
  "M45.5",   # rigid cystoscopy + catheter insertion for micturating
             # cystogram — diagnostic imaging procedure, not surgery

  # X55.1 was initially grouped with the amputation/mass-excision
  # cluster below under the assumption it meant mass excision, like
  # X53.1/X53.2. Review of further instances shows it's consistently
  # coded for "biopsy abdominal mass" alone — a diagnostic sampling
  # procedure, the same oncology-diagnostic pattern as N13.4, not
  # genuine major surgery — hence its place here instead.
  "X55.1",

  # X27.3 — excision supernumerary great toe, bundled in the same row
  # with a literal "OPERATION NOT FOUND" placeholder. Isolated
  # supernumerary toe excision is not typically a PCA-level
  # procedure clinically, and the co-occurring placeholder suggests
  # this admission's real driving procedure may simply be missing
  # from operations_clean, rather than the toe excision being the
  # genuine substantive surgery.
  "X27.3"
)

operations_clean <- operations_clean |>
  mutate(
    all_opcs_in_row = str_split(coalesce(PFMD_PROC_OPCS_CODES, ""), "\r\n"),
    primary_opcs_substantive = map2_chr(
      all_opcs_in_row, primary_opcs,
      function(codes, fallback) {
        # Each candidate is checked against a letter+digits pattern
        # before being considered, since placeholder text like "N/A"
        # (used when a coder had no code to enter) would otherwise be
        # accepted as a literal code — its chapter-letter extraction
        # matches the leading "N" purely by coincidence, which lands
        # genuinely uncoded admissions in "Male genitalia" (chapter N)
        # for no clinically meaningful reason. One case in the PSM
        # population illustrates the problem directly: an "N09"-
        # bucketed episode whose real procedure was "failed
        # laparoscopy converted to open / emergency appendicectomy",
        # with PFMD_PROC_OPCS_CODES literally containing "N/A" rather
        # than a real code.
        valid_codes <- codes[str_detect(codes, "^[A-Z]\\d+\\.?\\d*") & codes != ""]
        non_minor <- valid_codes[!(valid_codes %in% opcs_treat_as_minor_for_category)]
        if (length(non_minor) > 0) non_minor[1] else fallback
      }
    ),
    substantive_opcs_chapter = str_extract(primary_opcs_substantive, "^[A-Z]"),
    # row_fully_minor — EVERY code on this specific operation row is
    # on the minor list, i.e. nothing in this operation record rises
    # above incidental. Used downstream (Flag 6, 4_flags.R) to build
    # case_type_refined — episodes whose attributed operation is
    # entirely minor get treated as effectively non-surgical, since
    # the labelled procedure can't plausibly explain PCA need.
    row_fully_minor = map_lgl(
      all_opcs_in_row,
      ~ all(.x[.x != ""] %in% opcs_treat_as_minor_for_category)
    )
  )

cat("\nprimary_opcs_substantive computed\n")
cat("Rows where substantive code differs from plain first code:",
    sum(operations_clean$primary_opcs_substantive != operations_clean$primary_opcs,
        na.rm = TRUE), "\n\n")

# ============================================================
# 6B. proc_category — NOW DERIVED FROM primary_opcs_substantive
#
# Base chapter mapping (unchanged categories/logic), applied to the
# CORRECTED substantive code rather than the positional first code.
# ============================================================

operations_clean <- operations_clean |>
  mutate(
    proc_category = case_when(
      substantive_opcs_chapter == "A"                ~ "Neurosurgery",
      substantive_opcs_chapter == "B"                ~ "Endocrine surgery",
      substantive_opcs_chapter %in% c("C", "D")     ~ "Ophthalmology/ENT",
      substantive_opcs_chapter == "E"                ~ "Respiratory/ENT",
      substantive_opcs_chapter == "F"                ~ "Dental/oral",
      substantive_opcs_chapter == "G"                ~ "Abdominal/GI",
      substantive_opcs_chapter == "H"                ~ "Colorectal",
      substantive_opcs_chapter %in% c("J", "K")     ~ "Hepatobiliary",
      substantive_opcs_chapter == "L"                ~ "Vascular/lines",
      substantive_opcs_chapter == "M"                ~ "Urology",
      substantive_opcs_chapter == "N"                ~ "Male genitalia",
      substantive_opcs_chapter %in% c("P", "Q")     ~ "Gynaecology",
      substantive_opcs_chapter %in% c("S", "T")     ~ "Skin/trauma/general",
      substantive_opcs_chapter %in% c("V", "W")     ~ "Orthopaedic/spine",
      substantive_opcs_chapter == "X"                ~ "Unclassified",
      TRUE                                            ~ "Other"
    )
  )

cat("Procedure category (from substantive code):\n")
print(table(operations_clean$proc_category))

# ============================================================
# 6C. proc_category — REMAINING MANUAL OVERRIDES
#
# Once primary_opcs_substantive correctly surfaces the real driving
# code in bundled rows (6A/6B above), most categories that might seem
# to need protecting (orchidopexy under "Male genitalia", major
# vascular codes, maxillofacial codes) turn out to already be correct
# under the base chapter mapping. Only a handful of genuine
# exceptions remain, handled below. Every category used here is one
# of the original 16 — no new categories are introduced.
# ============================================================

operations_clean <- operations_clean |>
  mutate(
    proc_category = case_when(

      # N13.4 (biopsy testis) — genuinely incidental, oncology-
      # adjacent (60% pre-dates admission's anaesthesia, chemo-
      # positive). Only fires when N13.4 IS the substantive code
      # (i.e. no non-minor companion exists in that row) — the
      # T30.9 case above no longer reaches this line at all, since
      # primary_opcs_substantive already resolved to T30.9 for that
      # row, correctly landing it in Skin/trauma/general via 6B.
      primary_opcs_substantive == "N13.4" ~ "Unclassified",

      # Unclassified (chapter X) split: the DDH/foot deformity
      # correction cluster is genuine orthopaedic surgery
      primary_opcs_substantive %in% c("X22.1", "X22.2", "X22.3", "X22.5",
                                      "X24.2", "X24.3", "X25.1", "X25.4",
                                      "X48.1") ~ "Orthopaedic/spine",

      # Amputation/mass-excision cluster — genuine major surgery,
      # rehomed to an existing category rather than an invented one;
      # Skin/trauma/general already holds comparable content (lymph
      # node excisions, laparotomies, mass-related procedures). X55.1
      # is deliberately excluded here — see opcs_treat_as_minor_for_
      # category above, where it's classified as a biopsy rather than
      # an excision, the same pattern as N13.4.
      primary_opcs_substantive %in% c("X09.3", "X09.5", "X10.1", "X10.9",
                                      "X12.6", "X14.8", "X53.1",
                                      "X53.2") ~ "Skin/trauma/general",

      # The base chapter mapping above has no explicit rule for
      # chapter O, so O-codes would otherwise fall through to the
      # generic "Other" catch-all. Only two O-codes exist in this
      # data:
      # O27.1 — genuine orthopaedic surgery ("reconstruction joint
      # stabilisation" / "reconstruction knee joint stabilisation")
      primary_opcs_substantive == "O27.1" ~ "Orthopaedic/spine",
      # O56.2 — "portacath removal", functionally identical to L94.8
      # (same chapter-L code, already on opcs_exclude) just coded
      # under chapter O in this instance — categorised to match
      # where L94.3/L94.8 already live
      primary_opcs_substantive == "O56.2" ~ "Vascular/lines",

      # ---- Clinical recategorisation ----
      # These two overrides move genuinely substantive procedures to
      # the category a clinical reader would actually expect, rather
      # than the strict OPCS chapter letter — same
      # principle in both cases: OPCS classifies by anatomical
      # system/structure operated on, not by which specialty performs
      # or clinically "owns" the procedure.
      #
      # N09.2 (orchidopexy one stage) / N09.4 (laparoscopic
      # orchidopexy second stage) — OPCS chapter N ("Male genitalia")
      # by anatomical structure, but orchidopexy is a standard
      # paediatric UROLOGY procedure in clinical practice.
      primary_opcs_substantive %in% c("N09.2", "N09.4") ~ "Urology",

      # L75.1 (craniotomy resection/ligation congenital
      # arteriovenous malformation, posterior fossa) — OPCS chapter L
      # ("Vascular/lines") because an AVM is technically a vascular
      # malformation, but the actual operation is a posterior fossa
      # CRANIOTOMY performed by neurosurgeons — clinically
      # NEUROSURGERY, not vascular surgery.
      primary_opcs_substantive == "L75.1" ~ "Neurosurgery",

      # X27.3 (excision supernumerary great toe) deliberately gets no
      # override here and falls through to the base chapter-X mapping
      # ("Unclassified"). Two reasons: it's not typically a PCA-level
      # procedure clinically — brief, low-morbidity, day-case in
      # nature — and the same row bundles a literal "OPERATION NOT
      # FOUND" placeholder, suggesting this admission's operation
      # record is already incomplete and the real driving procedure
      # may be missing entirely. It's also included in
      # opcs_treat_as_minor_for_category above so it's treated as
      # non-explanatory for case_type_refined as well.

      TRUE ~ proc_category
    )
  )

cat("\nproc_category final overrides applied\n")
cat("Final procedure category distribution:\n")
print(table(operations_clean$proc_category))

# ============================================================
# 6D. EPISODE-LEVEL PROCEDURE ATTRIBUTION
#
# Each PCA episode is attributed to whichever operation is nearest in
# time to it (by OpDate), rather than to the first operation recorded
# for the admission as a whole. This matters for multi-operation
# admissions (15.9% of S4): a simple admission-level attribution
# would apply the same procedure category to every episode in the
# admission regardless of which operation the episode actually
# followed, whereas episode-level attribution lets each episode carry
# the category of the operation it's clinically linked to — affecting
# 234/867 multi-op episodes (29.4% of multi-op, 4.3% of all S4).
#
# Reads proc_category as computed in 6A-6C above (from
# primary_opcs_substantive, which already accounts for bundled codes
# within a single operation row).
#
# S3/PSM is unaffected by this — those admissions are already
# restricted to a single operation via has_multiple_operations.
#
# One consequence worth noting: proc_category can genuinely differ
# across episodes within the same admission. That's the intended
# behaviour, not an inconsistency — it reflects each episode's actual
# clinical context rather than treating the admission as a single
# undifferentiated unit.
#
# Limitation: OpDate is date-only (0% of operations_clean rows have a
# non-midnight time component), so same-day multiple operations can't
# be perfectly ordered within the day. Ties are broken by lowest
# OP_ORDER, consistent with how "primary" is defined elsewhere in
# this pipeline.
#
# 4_flags.R's demographics section relies on proc_category and
# PFMD_PROC_OPCS_CODES arriving already attached to pca_episodes from
# this section, and does not re-join operations_clean itself.
# ============================================================

cat("\nAttributing proc_category at episode level (nearest preceding\n")
cat("operation by OpDate, not admission-level first-row)...\n\n")

all_ops_dated <- operations_clean |>
  select(PAT_ENC_CSN_ID, OpDate, OP_ORDER, proc_category, PFMD_PROC_OPCS_CODES,
         row_fully_minor) |>
  mutate(OpDate = as.Date(OpDate))

episode_nearest_op <- pca_episodes |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start) |>
  mutate(episode_start_date = as.Date(episode_start)) |>
  left_join(all_ops_dated, by = "PAT_ENC_CSN_ID", relationship = "many-to-many") |>
  filter(OpDate <= episode_start_date | is.na(OpDate)) |>
  group_by(PAT_ENC_CSN_ID, episode_id) |>
  slice_max(OpDate, n = 1, with_ties = TRUE) |>
  slice_min(OP_ORDER, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(PAT_ENC_CSN_ID, episode_id, proc_category, PFMD_PROC_OPCS_CODES,
         row_fully_minor, attributed_op_order = OP_ORDER)

pca_episodes <- pca_episodes |>
  left_join(episode_nearest_op, by = c("PAT_ENC_CSN_ID", "episode_id"))

cat("Episode-level proc_category attached to pca_episodes\n")
cat("  Episodes with a category (operation found):",
    sum(!is.na(pca_episodes$proc_category)), "\n")
cat("  NA (genuinely medical, no operation this admission):",
    sum(is.na(pca_episodes$proc_category)), "\n\n")
cat("Episode-level procedure category distribution:\n")
print(table(pca_episodes$proc_category, useNA = "ifany"))

# ============================================================
# 7. ANAESTHESIA EVENTS
# ============================================================

anaesthesia_events_clean <- anaesthesia_events |>
  mutate(across(
    c(anaesthesia_start, induction_of_anaesthesia, intubation,
      spinal_in, invasive_lines_inserted, anaesthesia_ready,
      procedure_start, procedure_finish, extubation,
      line_safety_complete, anaesthesia_stop),
    as.POSIXct
  )) |>
  mutate(
    anaesthesia_duration_mins = as.numeric(
      difftime(anaesthesia_stop, anaesthesia_start, units = "mins")
    ),
    procedure_duration_mins = as.numeric(
      difftime(procedure_finish, procedure_start, units = "mins")
    ),
    was_intubated = !is.na(intubation),
    had_spinal    = !is.na(spinal_in)
  )

cat("Anaesthesia events clean:",
    nrow(anaesthesia_events_clean), "rows\n")
cat("\nAnaesthesia duration summary (mins):\n")
print(summary(anaesthesia_events_clean$anaesthesia_duration_mins))

# ============================================================
# 7B. ANAESTHESIA ASSESSMENT SDEs — ASA STATUS & PROCEDURE URGENCY
#
# Joins directly onto operations via anae_CSN_id (SDE file) =
# ANAESTHESIA_CSN_ID (operations, raw) — a clean 100% match
# (4,417/4,417), so no bridge via anaesthesia_events_clean is needed
# here, unlike the LDA epidural CSN mismatch handled below.
# ============================================================

anaesthesia_ax_clean <- anaesthesia_ax_sdes |>
  left_join(
    operations |> select(ANAESTHESIA_CSN_ID, PAT_ENC_CSN_ID, OP_ORDER),
    by = c("anae_CSN_id" = "ANAESTHESIA_CSN_ID")
  ) |>
  mutate(
    asa_grade = str_extract(ASA_status, "\\d"),
    asa_grade = factor(asa_grade, levels = as.character(1:5), ordered = TRUE),
    proc_urgency = factor(
      str_to_lower(procedure_urgency),
      levels = c("elective", "expedited", "urgent", "immediate"),
      ordered = TRUE
    ),
    # NEW: binary collapse — anything not elective is emergency
    urgency_binary = case_when(
      proc_urgency == "elective"                              ~ "elective",
      proc_urgency %in% c("expedited", "urgent", "immediate")  ~ "emergency",
      TRUE                                                      ~ NA_character_
    )
  )

cat("\nUrgency (binary) distribution:\n")
print(table(anaesthesia_ax_clean$urgency_binary, useNA = "always"))

# ============================================================
# 7C. LDA CATHETERS (EPIDURAL & URETHRAL)
#
# Both sheets use a far-future sentinel (~2157-11-19) in
# removal_instant to mean "removal not documented" rather than NA.
# Converted to an explicit flag + NA below — otherwise duration
# calculations downstream would silently produce ~131-year durations.
#
# Separately, a small number of rows (19 epidural, unchecked for
# urethral) have placement_instant missing while removal_instant is
# genuine — flagged as placement_not_documented rather than dropped,
# since these rows still confirm the catheter existed.
#
# Note: the LDA epidural and epidural assessment tables cap at
# ~2026-05-15/18 in this data refresh, unlike every other table,
# which extends to mid-July 2026. This is a genuine gap in the
# source extract rather than a processing artefact here.
# ============================================================

lda_epidural_clean <- lda_epidural |>
  rename(PAT_ENC_CSN_ID = pat_enc_CSN_id) |>   # standardise naming — raw table uses inconsistent casing
  mutate(
    placement_not_documented = is.na(placement_instant),
    removal_not_documented   = year(removal_instant) > 2100,
    removal_instant          = if_else(removal_not_documented, as.POSIXct(NA), removal_instant)
  )

lda_urethral_clean <- lda_urethral |>
  rename(PAT_ENC_CSN_ID = pat_enc_CSN_id) |>   # same fix — confirm this table has the identical column name first
  mutate(
    placement_not_documented = is.na(placement_instant),
    removal_not_documented   = year(removal_instant) > 2100,
    removal_instant          = if_else(removal_not_documented, as.POSIXct(NA), removal_instant)
  )

cat("\nLDA Epidural clean:", nrow(lda_epidural_clean), "rows\n")
lda_epidural_clean |>
  summarise(
    total = n(),
    fully_complete         = sum(!placement_not_documented & !removal_not_documented),
    placement_missing_only = sum(placement_not_documented & !removal_not_documented),
    removal_missing_only   = sum(!placement_not_documented & removal_not_documented),
    both_missing           = sum(placement_not_documented & removal_not_documented)
  ) |>
  print()

cat("\nLDA Urethral clean:", nrow(lda_urethral_clean), "rows\n")
lda_urethral_clean |>
  summarise(
    total = n(),
    fully_complete         = sum(!placement_not_documented & !removal_not_documented),
    placement_missing_only = sum(placement_not_documented & !removal_not_documented),
    removal_missing_only   = sum(!placement_not_documented & removal_not_documented),
    both_missing           = sum(placement_not_documented & removal_not_documented)
  ) |>
  print()

# LDA epidural assessments — no sentinel/duration fields to fix,
# just carried through as-is for now
lda_epidural_assess_clean <- lda_epidural_assess

cat("\nLDA Epidural Assess clean:", nrow(lda_epidural_assess_clean), "rows\n")

# ============================================================
# 8. PAIN SCORES
# ============================================================

pain_scores_clean <- pain_scores |>
  mutate(
    dow     = wday(score_time, label = TRUE, week_start = 1),
    weekend = ifelse(dow %in% c("Sat", "Sun"), "Weekend", "Weekday"),
    pain_category = case_when(
      score_value == 0  ~ "No pain",
      score_value <= 3  ~ "Mild (1-3)",
      score_value <= 6  ~ "Moderate (4-6)",
      score_value >= 7  ~ "Severe (7-10)"
    )
  ) |>
  mutate(
    pain_score_shift = future_map_dbl(
      score_time,
      ~ match_date_to_shift_number(.x, shifts),
      .progress = TRUE
    )
  )

cat("Pain scores clean:", nrow(pain_scores_clean), "rows\n")
cat("\nPain category distribution:\n")
print(table(pain_scores_clean$pain_category))
cat("\nWeekend distribution:\n")
print(table(pain_scores_clean$weekend))

# ============================================================
# 9. ADD CASE TYPE TO ADMISSIONS
# ============================================================

admissions_clean <- admissions_clean |>
  mutate(
    case_type = case_when(
      PAT_ENC_CSN_ID %in% operations_clean$PAT_ENC_CSN_ID ~ "surgical",
      TRUE                                                  ~ "medical"
    )
  )

cat("\nCase type split:\n")
print(table(admissions_clean$case_type))

# ============================================================
# 9B. INTRAOPERATIVE MAR
# ============================================================

intraop_mar_clean <- intraop_mar |>
  mutate(
    PAT_ENC_CSN_ID = as.numeric(PAT_ENC_CSN_ID),
    generic_name   = tolower(str_squish(generic_name)),
    mar_route      = tolower(str_squish(coalesce(mar_route, "")))
  ) |>
  filter(!is.na(generic_name))

cat("Intraoperative MAR clean:", nrow(intraop_mar_clean), "rows\n")
cat("\nOp phase distribution:\n")
print(table(intraop_mar_clean$op_phase_desc))
cat("\nTop intraop drugs:\n")
intraop_mar_clean |>
  count(generic_name, sort = TRUE) |>
  head(20) |>
  print()
cat("\nRoute distribution:\n")
intraop_mar_clean |>
  count(mar_route, sort = TRUE) |>
  head(10) |>
  print()


# ============================================================
# 10. SAVE ALL CLEANED TABLES
# ============================================================

dir.create("1 - data/3 - processed_data",
           showWarnings = FALSE, recursive = TRUE)

saveRDS(admissions_clean,
        "1 - data/3 - processed_data/admissions_clean.rds")
saveRDS(pca_mar_clean,
        "1 - data/3 - processed_data/pca_mar_clean.rds")
saveRDS(pca_episodes,
        "1 - data/3 - processed_data/pca_episodes.rds")
saveRDS(non_pca_mar_clean,
        "1 - data/3 - processed_data/non_pca_mar_clean.rds")
saveRDS(nerve_blocks_clean,
        "1 - data/3 - processed_data/nerve_blocks_clean.rds")
saveRDS(operations_clean,
        "1 - data/3 - processed_data/operations_clean.rds")
saveRDS(anaesthesia_events_clean,
        "1 - data/3 - processed_data/anaesthesia_events_clean.rds")
saveRDS(anaesthesia_ax_clean,
        "1 - data/3 - processed_data/anaesthesia_ax_clean.rds")
saveRDS(lda_epidural_clean,
        "1 - data/3 - processed_data/lda_epidural_clean.rds")
saveRDS(lda_urethral_clean,
        "1 - data/3 - processed_data/lda_urethral_clean.rds")
saveRDS(lda_epidural_assess_clean,
        "1 - data/3 - processed_data/lda_epidural_assess_clean.rds")
saveRDS(pain_scores_clean,
        "1 - data/3 - processed_data/pain_scores_clean.rds")
saveRDS(intraop_mar_clean,
        "1 - data/3 - processed_data/intraop_mar_clean.rds")

toc()
cat("\nAll cleaned tables saved to 1 - data/3 - processed_data/\n")
cat("\nFiles saved:\n")
list.files("1 - data/3 - processed_data/")

# Reset to sequential processing
plan(sequential)
