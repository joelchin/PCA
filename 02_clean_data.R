# ============================================================
# 02 - Clean data
# ============================================================
# Cleans each raw table loaded by 01_load_raw_data.R, builds the PCA/NCA
# episodes and saves the cleaned tables to 1 - data/3 - processed_data/.
# Slow (about 30 minutes, dominated by the shift-matching steps), so
# downstream scripts reuse the saved outputs.
#
# proc_category and PFMD_PROC_OPCS_CODES are attached to pca_episodes at
# episode level in section 6D (nearest preceding operation), not by an
# admission-level join; 04_derive_episode_flags.R relies on this.
# ============================================================

library(tidyverse)
library(lubridate)
library(furrr)
library(future)
library(tictoc)

source(file.path(getOption("pca_code_dir", "."), "drug_lists.R"))

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
# SHIFT TABLE
# ============================================================
# 12-hour blocks anchored at 2014-01-01 08:00. If the table does not
# exist it is created; it is then extended forward to the end of the
# calendar year of the latest timestamp in the raw data, so that no
# timestamp falls outside the grid (unmatched timestamps would be NA and
# the episode would silently drop out of the shift-level tables).

shifts_path <- "1 - data/3 - processed_data/1_shifts.rds"
dir.create(dirname(shifts_path), showWarnings = FALSE, recursive = TRUE)

if (file.exists(shifts_path)) {
  shifts <- readRDS(shifts_path)
} else {
  grid_start <- as.POSIXct("2014-01-01 08:00:00", tz = "UTC")
  shifts <- tibble(
    shift_number = 1L,
    start_time   = grid_start,
    end_time     = grid_start + hours(12) - seconds(1)
  )
}
cat("Shifts loaded:", nrow(shifts), "shifts\n")

latest_timestamp <- max(
  as.POSIXct(c(admissions$ADM_DATE, admissions$DISCHARGE_DATE)),
  as.POSIXct(c(pca_mar$TAKEN_TIME, pca_mar$Discontinue_Time)),
  as.POSIXct(non_pca_mar$ADMIN_TIME),
  as.POSIXct(pain_scores$score_time),
  na.rm = TRUE
)

target_ceiling <- as.POSIXct(
  paste0(year(latest_timestamp), "-12-31 23:59:59"),
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
  saveRDS(shifts, shifts_path)
  
  cat("Shifts table extended:", nrow(new_shifts), "new blocks added\n")
  cat("   New coverage ceiling:", format(max(shifts$end_time)), "\n")
  
} else {
  cat("Shifts table already extends to", format(max(shifts$end_time)), "\n")
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
    
    # Background infusion rate (mL/h) parsed from the PROGRAMME text;
    # the raw INFUSION_RATE column is too sparse to use.
    background_rate_ml = as.numeric(
      str_extract(PROGRAMME, "(?i)(?<=continuous infusion rate:\\s)[0-9.]+")
    ),
    
    # Prescribed loading dose (mL); "None" gives NA, not 0
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
# CONCENTRATION (mg per mL), used by 07_calculate_doses.R
# ============================================================
# Taken from MED_ORDER: oxycodone_mg / morphine_mg hold the total mg in
# the syringe, so dividing by the diluent volume gives the concentration.
# DOSE/DOSE_UNIT and the PROGRAMME text are too sparsely populated to
# use. Weight-based orders ("X mg/kg") have no fixed total mg and are
# excluded.

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
# 6A. primary_opcs_substantive: resolving bundled OPCS codes
# ============================================================
# A single operation row can carry several OPCS codes. Taking the first
# code is positional rather than clinical (e.g. "N13.4 \r\n T30.9" would
# pick the incidental testis biopsy over the laparotomy), so
# primary_opcs_substantive takes the first code on the row that is not
# on the minor list, falling back to primary_opcs only when every code
# is minor. proc_category and the episode-level attribution in 6D are
# built from it; primary_opcs itself is kept unchanged.
#
# opcs_treat_as_minor_for_category is used only to choose proc_category
# and for case_type_refined. It extends opcs_exclude (which defines
# eligibility for the comparative population) with codes judged
# incidental on clinical grounds. It is defined as a standalone vector
# because an inline vector of this length would be recycled down the
# rows inside mutate().

opcs_treat_as_minor_for_category <- c(
  opcs_exclude,
  # Diagnostic or sampling procedures
  "N13.4",   # biopsy testis
  "M13.1",   # needle biopsy kidney
  "T87.2",   # excision biopsy cervical lymph node
  "S15.1",   # biopsy head lesion
  "S15.2",   # skin biopsy
  "F36.2",   # biopsy tonsil
  "P27.3",   # cystovaginoscopy
  "X55.1",   # biopsy abdominal mass
  "X55.8",   # examination under anaesthetic
  "M45.5",   # rigid cystoscopy for micturating cystogram
  # Minor drainage or evacuation
  "H44.6",   # EUA rectum, manual evacuation
  "P09.2",   # drainage vulval lesion
  "P31.2",   # drainage pouch of Douglas
  # Analgesic or anaesthetic techniques coded as procedures
  "A76.5",   # splanchnic sympathetic nerve block
  "A52.5",   # thoracic epidural
  # Incidental procedures
  "H02.4",   # incidental appendicectomy at laparoscopy
  "X27.3",   # supernumerary toe excision, bundled with an
             # "operation not found" placeholder
  # Tube and stoma maintenance
  "G34.4",   # closure gastrostomy
  "G34.5",   # removal and reinsertion of gastrostomy tube
  "G60.3",   # closure jejunostomy
  "G75.5"    # reduction ileostomy prolapse
)

operations_clean <- operations_clean |>
  mutate(
    all_opcs_in_row = str_split(coalesce(PFMD_PROC_OPCS_CODES, ""), "\r\n"),
    primary_opcs_substantive = map2_chr(
      all_opcs_in_row, primary_opcs,
      function(codes, fallback) {
        # Only strings that look like an OPCS code are candidates;
        # placeholders such as "N/A" would otherwise be read as chapter N.
        valid_codes <- codes[str_detect(codes, "^[A-Z]\\d+\\.?\\d*") & codes != ""]
        non_minor <- valid_codes[!(valid_codes %in% opcs_treat_as_minor_for_category)]
        if (length(non_minor) > 0) non_minor[1] else fallback
      }
    ),
    substantive_opcs_chapter = str_extract(primary_opcs_substantive, "^[A-Z]"),
    # TRUE when every code on this operation row is on the minor list;
    # used in 04_derive_episode_flags.R for case_type_refined
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
# 6B. proc_category from primary_opcs_substantive (OPCS chapter)
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
# 6C. proc_category: manual overrides
# ============================================================
# Exceptions to the chapter mapping. No new categories are introduced.

operations_clean <- operations_clean |>
  mutate(
    proc_category = case_when(
      
      # Testis biopsy is incidental when it is the only substantive code
      primary_opcs_substantive == "N13.4" ~ "Unclassified",
      
      # Deformity-correction codes in chapter X are orthopaedic surgery
      primary_opcs_substantive %in% c("X22.1", "X22.2", "X22.3", "X22.5",
                                      "X24.2", "X24.3", "X25.1", "X25.4",
                                      "X48.1") ~ "Orthopaedic/spine",
      
      # Amputation and mass-excision codes in chapter X
      primary_opcs_substantive %in% c("X09.3", "X09.5", "X10.1", "X10.9",
                                      "X12.6", "X14.8", "X53.1",
                                      "X53.2") ~ "Skin/trauma/general",
      
      # Chapter O has no rule in the chapter mapping
      primary_opcs_substantive == "O27.1" ~ "Orthopaedic/spine",
      primary_opcs_substantive == "O56.2" ~ "Vascular/lines",
      
      # Clinical rather than anatomical grouping: orchidopexy is
      # paediatric urology, and craniotomy for AVM is neurosurgery
      primary_opcs_substantive %in% c("N09.2", "N09.4") ~ "Urology",
      primary_opcs_substantive == "L75.1" ~ "Neurosurgery",
      
      TRUE ~ proc_category
    )
  )

cat("\nproc_category final overrides applied\n")
cat("Final procedure category distribution:\n")
print(table(operations_clean$proc_category))

# ============================================================
# 6D. Episode-level procedure attribution
# ============================================================
# Each episode is attributed to the operation nearest before it by
# OpDate, rather than to the first operation of the admission, so that
# episodes in multi-operation admissions carry the category of the
# operation they followed. OpDate is date-only, so same-day operations
# are ordered by lowest OP_ORDER. proc_category can therefore differ
# between episodes of one admission. The primary comparative population
# is restricted to single-operation admissions, so is unaffected.

cat("\nAttributing proc_category at episode level...\n\n")

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
# 7B. Anaesthesia assessment: ASA grade
# ============================================================
# Joins to operations via anae_CSN_id = ANAESTHESIA_CSN_ID.

anaesthesia_ax_clean <- anaesthesia_ax_sdes |>
  left_join(
    operations |> select(ANAESTHESIA_CSN_ID, PAT_ENC_CSN_ID, OP_ORDER),
    by = c("anae_CSN_id" = "ANAESTHESIA_CSN_ID")
  ) |>
  mutate(
    asa_grade = str_extract(ASA_status, "\\d"),
    asa_grade = factor(asa_grade, levels = as.character(1:5), ordered = TRUE)
  )

cat("\nASA grade distribution:\n")
print(table(anaesthesia_ax_clean$asa_grade, useNA = "always"))

# ============================================================
# 7C. LDA catheters (epidural and urethral)
# ============================================================
# Both sheets use a far-future sentinel date (year > 2100) in
# removal_instant for "removal not documented". It is converted to NA
# with an explicit flag, otherwise durations would be about 130 years.
# Rows with no placement time but a genuine removal time are flagged
# (placement_not_documented) rather than dropped, since they still
# confirm that a catheter existed.

lda_epidural_clean <- lda_epidural |>
  rename(PAT_ENC_CSN_ID = pat_enc_CSN_id) |>   # raw column name has inconsistent casing
  mutate(
    placement_not_documented = is.na(placement_instant),
    removal_not_documented   = year(removal_instant) > 2100,
    removal_instant          = if_else(removal_not_documented, as.POSIXct(NA), removal_instant)
  )

lda_urethral_clean <- lda_urethral |>
  rename(PAT_ENC_CSN_ID = pat_enc_CSN_id) |>
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

# Epidural assessments need no cleaning
lda_epidural_assess_clean <- lda_epidural_assess

cat("\nLDA Epidural Assess clean:", nrow(lda_epidural_assess_clean), "rows\n")

# ============================================================
# 8. PAIN SCORES
# ============================================================

pain_scores_clean <- pain_scores |>
  mutate(
    pain_score_shift = future_map_dbl(
      score_time,
      ~ match_date_to_shift_number(.x, shifts),
      .progress = TRUE
    )
  )

cat("Pain scores clean:", nrow(pain_scores_clean), "rows\n")

# ============================================================
# 9. CASE TYPE (surgical if any operation is recorded for the admission)
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
print(list.files("1 - data/3 - processed_data/"))

plan(sequential)
plan(sequential)