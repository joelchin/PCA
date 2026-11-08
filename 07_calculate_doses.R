# ============================================================
# 07 - Opioid dose per episode (background infusion)
# ============================================================
# Builds the episode-level background dose used by later scripts
# (episode_dose.rds), from the syringe records in pca_mar_clean:
#
#   background_rate_mg_hr   = background_rate_ml x concentration (mg/mL)
#   syringe_background_mg   = background_rate_mg_hr x syringe duration
#   total_background_mg     = sum over the syringes of the episode
#
# The consumption metric is this delivered background dose per kg per
# hour. It is not the drug content of the syringe label (oxycodone_mg or
# morphine_mg in pca_mar_clean), which is the nominal fill strength and
# overstates what was infused. Bolus doses are added in
# 08_calculate_true_total_dose.R. The prescribed loading dose is a
# contingency, not a delivered dose, and is not added to any total.
# ============================================================

library(tidyverse)
library(lubridate)

stopifnot(
  exists("pca_mar_clean"),
  exists("pca_episodes_flagged"),
  exists("episode_metrics")
)

# ============================================================
# STEP 1 - Syringe-level dose dataset
# ============================================================

cat("--- Building syringe-level dose dataset ---\n\n")

# Morphine orders in the ">50 kg" formulation ("morphine 50 mg in 50 mL
# paediatric PCA/NCA (> 50 kg)") have no concentration from the
# extraction in 02_clean_data.R; the concentration is 1 mg/mL.
pca_mar_clean <- pca_mar_clean |>
  mutate(
    morphine_conc_mg_per_ml = case_when(
      !is.na(morphine_conc_mg_per_ml) ~ morphine_conc_mg_per_ml,
      drug == "morphine" &
        MED_ORDER == "morphine 50 mg in 50 mL paediatric PCA/NCA (> 50 kg)" ~ 1.0,
      TRUE ~ morphine_conc_mg_per_ml
    )
  )

n_missing_after <- sum(is.na(pca_mar_clean$morphine_conc_mg_per_ml[
  pca_mar_clean$drug == "morphine"]))
cat(sprintf("Morphine syringes with missing concentration: %d\n\n", n_missing_after))

# Assign each syringe record to its episode: TAKEN_TIME must fall within
# episode_start to episode_end and the drug must match. If a syringe
# matches more than one episode, the closest start is used.

episode_ends <- pca_episodes_flagged |>
  select(PAT_ENC_CSN_ID, episode_id, episode_start, episode_end,
         WEIGHT_KG, AGE, age_group, drug)

cat("Linking syringe records to episodes via time overlap...\n")

syringe_episode_link <- pca_mar_clean |>
  filter(drug %in% c("oxycodone", "morphine")) |>
  select(PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME, drug) |>
  left_join(
    episode_ends,
    by = c("PAT_ENC_CSN_ID", "drug"),
    relationship = "many-to-many"
  ) |>
  filter(
    TAKEN_TIME >= episode_start,
    is.na(episode_end) | TAKEN_TIME <= episode_end
  ) |>
  group_by(PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME, drug) |>
  slice_min(abs(as.numeric(TAKEN_TIME - episode_start)), n = 1) |>
  ungroup() |>
  select(PAT_ENC_CSN_ID, ORDER_MED_ID, TAKEN_TIME, drug, episode_id,
         episode_end, WEIGHT_KG, AGE, age_group)

cat(sprintf("Syringe records linked: %d / %d\n",
            nrow(syringe_episode_link),
            sum(pca_mar_clean$drug %in% c("oxycodone", "morphine"))))
cat("\n")

dose_data <- pca_mar_clean |>
  filter(drug %in% c("oxycodone", "morphine")) |>
  left_join(syringe_episode_link,
            by = c("PAT_ENC_CSN_ID", "ORDER_MED_ID", "TAKEN_TIME", "drug"),
            relationship = "many-to-many") |>
  filter(!is.na(episode_id)) |>
  mutate(
    concentration_mg_per_ml = case_when(
      drug == "oxycodone" ~ oxycodone_conc_mg_per_ml,
      drug == "morphine" ~ morphine_conc_mg_per_ml,
      TRUE ~ NA_real_
    ),
    # Missing Discontinue_Time is imputed with episode_end
    syringe_end = if_else(
      is.na(Discontinue_Time),
      episode_end,
      Discontinue_Time
    ),
    syringe_duration_hrs = as.numeric(
      difftime(syringe_end, TAKEN_TIME, units = "hours")
    ),
    # Durations over 7 days or negative are treated as data errors
    syringe_duration_hrs = if_else(
      syringe_duration_hrs > 168 | syringe_duration_hrs < 0,
      NA_real_,
      syringe_duration_hrs
    ),
    background_rate_mg_hr = background_rate_ml * concentration_mg_per_ml,
    background_rate_mg_kg_hr = background_rate_mg_hr / WEIGHT_KG,
    bolus_dose_mg = bolus_dose_ml * concentration_mg_per_ml,
    bolus_dose_mg_kg = bolus_dose_mg / WEIGHT_KG,
    loading_dose_mg = loading_dose_ml * concentration_mg_per_ml,
    loading_dose_mg_kg = loading_dose_mg / WEIGHT_KG,
    syringe_background_mg = background_rate_mg_hr * syringe_duration_hrs,
    syringe_background_mg_kg = syringe_background_mg / WEIGHT_KG
  )

cat("Syringe records built:\n")
cat(sprintf(" Oxycodone: %d\n", sum(dose_data$drug == "oxycodone")))
cat(sprintf(" Morphine: %d\n", sum(dose_data$drug == "morphine")))
cat(sprintf(" Missing concentration: %d\n",
            sum(is.na(dose_data$concentration_mg_per_ml))))
cat(sprintf(" Missing syringe duration after imputation: %d\n",
            sum(is.na(dose_data$syringe_duration_hrs))))
cat("\n")

# ============================================================
# STEP 2 - Episode-level dose summary
# ============================================================

cat("--- Aggregating to episode level ---\n\n")

episode_dose <- dose_data |>
  group_by(PAT_ENC_CSN_ID, episode_id, drug, AGE, age_group, WEIGHT_KG) |>
  summarise(
    n_syringes = n(),
    # Median of the per-syringe prescribed rates (does not weight by how
    # long each syringe ran)
    background_rate_mg_kg_hr_median = median(background_rate_mg_kg_hr,
                                             na.rm = TRUE),
    # Prescribed bolus dose (first syringe)
    bolus_dose_mg_kg_first = first(bolus_dose_mg_kg),
    # Prescribed loading dose (first syringe); not part of delivered dose
    had_loading_dose = any(!is.na(loading_dose_mg_kg)),
    loading_dose_mg_kg = first(loading_dose_mg_kg[!is.na(loading_dose_mg_kg)]),
    # Background dose delivered: sum over syringes of rate x duration
    total_background_mg = sum(syringe_background_mg, na.rm = TRUE),
    total_background_mg_kg = sum(syringe_background_mg_kg, na.rm = TRUE),
    # Episode duration from dose data (sum of syringe durations)
    episode_duration_dose_hrs = sum(syringe_duration_hrs, na.rm = TRUE),
    daily_background_mg_kg = if_else(
      episode_duration_dose_hrs > 0,
      total_background_mg_kg / episode_duration_dose_hrs * 24,
      NA_real_
    ),
    concentration_mg_per_ml = first(concentration_mg_per_ml),
    .groups = "drop"
  ) |>
  mutate(
    # Duration-weighted background rate (mg/kg/h)
    background_mg_kg_hr_true = if_else(
      episode_duration_dose_hrs > 0,
      total_background_mg_kg / episode_duration_dose_hrs,
      NA_real_
    ),
    # Implausible daily background doses are excluded below
    flag_implausible_dose = daily_background_mg_kg > 2 |
      daily_background_mg_kg < 0
  )

cat(sprintf("Episode dose summaries: %d\n", nrow(episode_dose)))
cat(sprintf(" Implausible doses flagged: %d\n",
            sum(episode_dose$flag_implausible_dose, na.rm = TRUE)))
cat("\n")

# Remove implausible daily doses (above 2 mg/kg/day or negative)
episode_dose_clean <- episode_dose |>
  filter(!flag_implausible_dose | is.na(flag_implausible_dose))

# ============================================================
# SAVE
# ============================================================

saveRDS(episode_dose_clean,
        "1 - data/3 - processed_data/episode_dose.rds")

# Later scripts use the name episode_dose
episode_dose <- episode_dose_clean

cat("episode_dose.rds saved\n")
cat("Background dose is prescribed rate x syringe duration; bolus doses\n")
cat("are added in 08_calculate_true_total_dose.R.\n")
