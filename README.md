# Paediatric PCA/NCA opioid analysis pipeline (oxycodone vs morphine)

R code for a retrospective, single-centre, propensity score-matched comparison
of oxycodone and morphine patient- and nurse-controlled analgesia (PCA/NCA) in
children, using an electronic health record (EPIC) extract.

> **No data are included in this repository.** The analysis uses NHS
> patient-level data that cannot be shared. This repository contains analysis
> code only. Manuscript citation: _to be added on publication_.

## Requirements

- R >= 4.4 (the manuscript analyses were run in R 4.5.3)
- Packages: `tidyverse`, `lubridate`, `readxl`, `MatchIt`, `logistf`,
  `survival`, `survminer`, `gtsummary`, `flextable`, `officer`, `furrr`,
  `future`, `tictoc`, `zip` (optionally `ragg` for the figure export).
  Script 01 installs missing CRAN binaries; `MatchIt`, `logistf` and
  `survminer` may need a manual install on locked-down machines.

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
(manual chart-review reclassification of naloxone administrations; used by script 04).

## Pipeline

Run from a fresh R session, with the working directory set to this folder, and
run the scripts in order in the same session. Intermediate outputs are written
as `.rds` files under `PCA_DATA_DIR`. Scripts 09 to 12 use objects left in the
session by the previous scripts (`psm_data`, `matched_data`, `match_out`, and
the Table 4 objects), so run them together.

| Script | Purpose |
|---|---|
| `01_load_raw_data.R` | Read raw extract, check inputs, set up packages |
| `02_clean_data.R` | Clean each table; build PCA/NCA episodes (slow; cached) |
| `03_build_shift_table.R` | Build 12-hour shift-level table and episode metrics |
| `04_derive_episode_flags.R` | Exclusion flags, concurrent-epidural detection, outcome definitions, pain scores |
| `05_load_processed_data.R` | Load the processed data at the start of an analysis session |
| `06_classify_cohorts.R` | Analysis populations and population flow counts |
| `07_calculate_doses.R` | Background infusion dose per episode |
| `08_calculate_true_total_dose.R` | Total dose (background + bolus) from pump flowsheets |
| `09_comparative_psm.R` | Propensity score matching |
| `10_time_to_effect.R` | Time to a sustained pain score of 3 or below (Figure 2) |
| `11_rotation_analysis.R` | Opioid rotation (Table 4) |
| `12_publication_tables.R` | Manuscript tables; comparison with Howard et al. (2010) |

`drug_lists.R` (shared drug, nerve-block and OPCS classification lists) is
sourced by scripts 02 to 05 and must sit in the same folder as the scripts.

## Data governance

Analysis was performed within the NHS Trust's secure environment. Scripts are
written to print aggregate output only. Do not commit data, extracts, logs or
outputs containing patient-level information.

## Licence

_To be added._
