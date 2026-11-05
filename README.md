# Paediatric PCA/NCA opioid analysis pipeline (oxycodone vs morphine)

R code for a retrospective, single-centre, propensity score-matched comparison
of oxycodone and morphine patient- and nurse-controlled analgesia (PCA/NCA) in
children, using an electronic health record (EPIC) extract.

> **No data are included in this repository.** The analysis uses NHS
> patient-level data that cannot be shared. This repository contains analysis
> code only. Manuscript citation: _to be added on publication_.

## Requirements

- R >= 4.4
- Packages: `tidyverse`, `lubridate`, `readxl`, `MatchIt`, `logistf`,
  `sandwich`, `gtsummary`, `flextable`, `officer`, `furrr`, `future`, `tictoc`, `zip`
  (script 01 installs missing CRAN binaries; `MatchIt` and `logistf` may need a
  manual install on locked-down machines).

## Configuration

Paths and extract file names are read from environment variables so that no
site-specific paths live in the code. Copy `.Renviron.example` into your
user-level `.Renviron` and edit:

| Variable | Meaning |
|---|---|
| `PCA_DATA_DIR` | Folder containing the raw extract and all derived data (outside this repo) |
| `PCA_FILE_PREFIX` / `PCA_FILE_SUFFIX` | Extract files are expected as `<prefix><table name><suffix>` |

Optional local-only files (never committed) in `PCA_DATA_DIR/manual_corrections/`:

- `gender_fix_pat_ids.txt`: patient IDs, one per line, for a documented manual sex correction (script 06)
- `corrupted_counter_csns.txt`: encounter IDs, one per line, with a confirmed corrupted pump counter (script 08)

Also required, produced outside this pipeline: `naloxone_reclassifications.rds`
(manually curated from chart review; used by script 04).

## Pipeline

Run from a fresh R session, in order. Intermediate outputs are written as
`.rds` files under `PCA_DATA_DIR`.

| Script | Purpose |
|---|---|
| `01_load_raw_data.R` | Read raw extract, check inputs, set up packages |
| `02_clean_data.R` | Clean each table; build PCA/NCA episodes (slow; cached) |
| `03_build_shift_table.R` | Build 12-hour shift-level analysis table and episode metrics |
| `04_derive_episode_flags.R` | Exclusion flags, concurrent-epidural detection, side-effect outcome flags |
| `05_load_processed_data.R` | Load processed data at the start of a session; verification checks |
| `06_classify_cohorts.R` | Surgical/non-surgical classification; analysis populations |
| `07_calculate_doses.R` | Opioid dose and consumption from syringe records |
| `08_calculate_true_total_dose.R` | Total dose (background + bolus) from pump flowsheets |
| `09_comparative_psm.R` | Propensity score matching and primary comparison |
| `10_opioid_consumption.R` | Consumption in the matched cohort |
| `11_oxycodone_descriptive.R` | Descriptive analysis of all oxycodone episodes |
| `12_safety_outcomes.R` | Safety outcomes |
| `13_subgroup_analyses.R` | Subgroup analyses |
| `14_sensitivity_analyses.R` | Sensitivity analyses |
| `15_rotation_analysis.R` | Opioid rotation analysis |
| `16_publication_tables.R` | Manuscript tables |
| `17_pain_score_analysis.R` | Pain score completeness, bias check and sensitivity analysis |
| `18_final_outcomes.R` | Final matched-cohort outcome summary with effect sizes |

`drug_lists.R` (shared drug, nerve-block and OPCS classification lists) is sourced by scripts 02, 04 and 05 and must sit in the same folder as them (`R/`).

## Data governance

Analysis was performed within the NHS Trust's secure environment. Scripts are
written to print aggregate output only. Do not commit data, extracts, logs or
outputs containing patient-level information.

## Licence

_To be added._
