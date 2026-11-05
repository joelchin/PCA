# ============================================================
# DRUG CLASSIFICATION LISTS — Paediatric PCA Study
# Addenbrooke's Hospital
# Author: J Chin
# ============================================================
#
# Shared drug, block and OPCS classification lists. Sourced by:
#   - 02_clean_data.R
#   - 04_derive_episode_flags.R
#   - 05_load_processed_data.R
#
# Edit lists HERE: changes propagate to all scripts.
# After editing, re-run 02_clean_data.R and everything downstream,
# UNLESS the change is only to an analysis-stage list (antiemetics,
# antipruritics etc.): then re-run from 03_build_shift_table.R.
# ============================================================

# ============================================================
# PCA DRUGS
# Classified from MEDICATION field in pca_mar
# Order matters in the 02_clean_data.R case_when — morphine+ketamine
# must be detected before morphine
# ============================================================

pca_drugs <- c(
  "morphine",
  "morphine+ketamine",
  "oxycodone",
  "fentanyl",
  "diamorphine"       # zero episodes — newer formulary
)

# ============================================================
# NON-PCA MAR — ANALGESIC CLASSIFICATION
# ============================================================

# Simple analgesics
simple_analgesics <- c(
  "paracetamol"
)

# NSAIDs — high confidence analgesic use only
# Deliberately excluded:
#   aspirin        — antiplatelet, not analgesic in paeds
#   indometacin    — PDA closure indication
#   meloxicam      — rheumatology indication
#   mefenamic acid — gynaecology indication
nsaid_drugs <- c(
  "ibuprofen",
  "diclofenac sodium",
  "diclofenac",
  "naproxen",
  "parecoxib"
)

# Opioids (ward/non-PCA)
opioid_drugs <- c(
  "morphine sulfate",
  "oxycodone",
  "oxycodone hydrochloride",
  "fentanyl",
  "codeine phosphate",
  "codeine",
  "dihydrocodeine",
  "tramadol",
  "alfentanil",
  "diamorphine hydrochloride",
  "methadone",
  "pethidine",
  "remifentanil"
)

# Neuropathic agents
neuropathic_drugs <- c(
  "gabapentin",
  "pregabalin",
  "amitriptyline hydrochloride"
)

# Adjuvants
adjuvant_drugs <- c(
  "ketamine",
  "esketamine",
  "clonidine",
  "dexmedetomidine",
  "magnesium sulfate",
  "baclofen"
)

# Combined analgesic list — used in 02_clean_data.R drug_category
analgesic_drugs <- c(
  simple_analgesics,
  nsaid_drugs,
  opioid_drugs,
  neuropathic_drugs,
  adjuvant_drugs
)

# ============================================================
# SIDE EFFECT MARKERS
# Confirmed from EPIC standard paediatric postop order set
#
# Dexamethasone EXCLUDED:
#   - Examination showed repeated scheduled dosing every 6-8hrs
#   - Predominantly intraop prophylaxis or anti-inflammatory
#   - 9,063 total administrations; pattern inconsistent with PONV
#   - Cannot reliably attribute to opioid-induced PONV
#
# All restricted to active episode shifts in 03_build_shift_table.R
# via inner_join — preop/intraop/post-episode excluded
#
# Primary analysis uses shift 2 onwards (Window B):
#   Shift 1 excluded — reflects prophylactic prescribing
#   at induction not established opioid side effects
#   Shift 1 antiemetic rate 55-71% vs shift 2+ rate 38-49%
#   confirming massive prophylaxis inflation in shift 1
# ============================================================

# Antiemetics — unambiguous postop use
antiemetic_drugs <- c(
  "ondansetron",
  "cyclizine",
  "metoclopramide",
  "droperidol",
  "haloperidol"
)

# Antipruritic
# Note: low-dose naloxone (40mcg) also used as pruritus proxy
# Distinguished by MEDICATION formulation in 02_clean_data.R
# naloxone_type = "low_dose_pruritus" for 40mcg formulation
antipruritic_drugs <- c(
  "chlorphenamine",
  "naloxone"          # 40mcg formulation only — see naloxone_type
)

# Opioid reversal — distinguished by formulation
# NALOXONE 40MICR  = low_dose_pruritus  (n=130)
# NALOXONE 400MIC  = reversal_respiratory (n=30) ← key safety outcome
reversal_drugs <- c(
  "naloxone"          # 400mcg formulation only — see naloxone_type
)

# Combined side effect marker list
# Used in 02_clean_data.R for drug_category = "side_effect_marker"
side_effect_drugs <- c(
  antiemetic_drugs,
  "chlorphenamine",
  "naloxone"
)

# ============================================================
# LAXATIVES
# Opioid-induced constipation proxy
# Expected temporal pattern: increases with shift_in_episode
# Early = pre-existing or non-opioid cause
# Late = more likely opioid-induced constipation
# ============================================================

laxative_drugs <- c(
  "macrogol 3350 + electrolytes",
  "lactulose",
  "senna",
  "docusate sodium",
  "bisacodyl",
  "sodium picosulfate"
)

# ============================================================
# CHEMOTHERAPY MARKERS
# Used to flag admissions with concurrent chemotherapy
# Chemotherapy-induced nausea (CINV) confounds antiemetic
# attribution — cannot distinguish CINV from opioid PONV
# Patients flagged for exclusion from primary side effect
# analysis (had_chemotherapy flag in 04_derive_episode_flags.R)
# ============================================================

chemo_drugs <- c(
  "methotrexate",
  "mercaptopurine",
  "cyclophosphamide",
  "vincristine",
  "doxorubicin",
  "cisplatin",
  "carboplatin",
  "temozolomide",
  "imatinib",
  "dasatinib",
  "crizotinib",
  "filgrastim",       # G-CSF — marker of chemotherapy
  "lenograstim"       # G-CSF — marker of chemotherapy
)

# ============================================================
# NERVE BLOCK CATEGORIES
# Used in 02_clean_data.R nerve_blocks cleaning section
# Classified by block_type field (free text — str_detect)
# ============================================================

neuraxial_blocks <- c(
  "epidural", "spinal", "intrathecal"
)

truncal_blocks <- c(
  "paravertebral", "intercostal", "serratus",
  "tap", "rectus", "ilioinguinal", "ring",
  "penile", "wound"
)

lower_limb_blocks <- c(
  "femoral", "sciatic", "popliteal", "saphenous",
  "tibial", "lumbar plexus", "fascia iliaca",
  "fascial iliaca", "3 in 1", "adductor"
)

upper_limb_blocks <- c(
  "axillary", "brachial", "median", "ulnar",
  "radial", "interscalene", "supraclavicular"
)

# ============================================================
# DAVID'S OPCS EXCLUSION LIST
# Non-substantive procedures — excluded from PSM population
# (Section 3 only — NOT from Sections 1, 2, or 4)
#
# Rule: exclude if ALL operations in admission are on this list
# If patient also had a real operation — keep them
# ============================================================

# ============================================================
# opcs_exclude — David's non-substantive procedure list
#
# Lives here (drug_lists.R) as the single source of truth, sourced
# by both 2_cleaningscript.R (primary_opcs_substantive extraction,
# Section 6A) and 4_flags.R (all_opcs_excluded, Flag 5) — previously
# duplicated only in 4_flags.R, moved here this session to avoid the
# two definitions drifting apart.
# ============================================================

# ============================================================
# opcs_exclude — David's non-substantive procedure list
#
# Lives here (drug_lists.R) as the single source of truth, sourced
# by both 2_cleaningscript.R (primary_opcs_substantive extraction,
# Section 6A) and 4_flags.R (all_opcs_excluded, Flag 5) — previously
# duplicated only in 4_flags.R, moved here this session to avoid the
# two definitions drifting apart.
# ============================================================

opcs_exclude <- c(
  "A54.2",              # Intrathecal chemotherapy
  "A55.9",              # Lumbar puncture
  "E42.1",              # Permanent tracheostomy
  "E42.3",              # Temporary tracheostomy
  "E49.2",              # Flexible bronchoscopy with lavage
  "E49.9",              # Flexible bronchoscopy
  "E85.1",              # Intubation and tracheostomy
  "G45.1",              # OGD with biopsy
  "G45.9",              # OGD
  "G47.8",              # Nasogastric tube insertion
  "G64.3",              # Nasojejunal tube insertion
  "H22.9",              # Colonoscopy
  "L91.1",              # Hickman line insertion
  "L91.2",              # Central line insertion
  "L91.3",              # Central venous catheter removal
  "L91.4",              # Hickman line removal
  "L91.5",              # Long line insertion
  "L91.8",              # Tunnelled venous catheter removal
  "L94.3",              # Portacath insertion
  "L94.8",              # Portacath removal
  "L99.7",              # PICC line insertion
  "M47.9",              # Urinary catheter insertion
  "O56.2",              # Portacath removal — chapter-O duplicate of
  # L94.8 (same clinical procedure, coded under
  # a different chapter in this instance; found
  # during the O-chapter completeness check,
  # [session date])
  "S43.3",              # Suture removal
  "S43.4",              # Suture removal
  "T34.9",              # Abdominal drain insertion
  "U05.2",              # MRI head
  "U05.4",              # CT lumbar spine
  "U05.5",              # MRI spine
  "U21.1",              # MRI
  "W36.5",              # Bone marrow aspirate/trephine
  "X56.8",              # Trial of extubation
  "X59.9"                # Operation not found
)

# ============================================================
# INTRAOPERATIVE DRUG DETECTION PATTERNS
# Used in 04_derive_episode_flags.R for intraop drug flags
# These are PSM COVARIATES not exclusion criteria
# Consistent with POPCORN trial pragmatic design
# ============================================================

# Detect by str_detect on generic_name (lowercase)
remifentanil_pattern   <- "remifentanil"
intraop_opioid_pattern <- "morphine|fentanyl|alfentanil|diamorphine|oxycodone"
ketamine_pattern       <- "ketamine"
clonidine_pattern      <- "clonidine"
neuraxial_route_pattern <- "intrathecal|epidural"

cat("✅ Drug lists loaded\n")