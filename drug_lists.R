# ============================================================
# Drug, nerve-block and OPCS classification lists
# ============================================================
# Sourced by 02_clean_data.R, 03_build_shift_table.R,
# 04_derive_episode_flags.R and 05_load_processed_data.R.
# Matching is on lower-case generic names. After editing a list,
# re-run from 02_clean_data.R.
# ============================================================

# ---- Analgesic classes (non-PCA MAR) ----

simple_analgesics <- c("paracetamol")

# NSAIDs with a clear analgesic indication only. Aspirin (antiplatelet),
# indometacin (PDA closure), meloxicam (rheumatology) and mefenamic acid
# (gynaecology) are deliberately omitted.
nsaid_drugs <- c(
  "ibuprofen",
  "diclofenac sodium",
  "diclofenac",
  "naproxen",
  "parecoxib"
)

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

neuropathic_drugs <- c(
  "gabapentin",
  "pregabalin",
  "amitriptyline hydrochloride"
)

adjuvant_drugs <- c(
  "ketamine",
  "esketamine",
  "clonidine",
  "dexmedetomidine",
  "magnesium sulfate",
  "baclofen"
)

analgesic_drugs <- c(
  simple_analgesics,
  nsaid_drugs,
  opioid_drugs,
  neuropathic_drugs,
  adjuvant_drugs
)

# ---- Side-effect markers ----

# Antiemetics used as the proxy for nausea/vomiting. Dexamethasone is
# excluded: it is given on the order set but is also used for
# inflammation and chemotherapy-induced nausea, so its administration
# cannot be attributed to opioid-related nausea.
antiemetic_drugs <- c(
  "ondansetron",
  "cyclizine",
  "metoclopramide",
  "droperidol",
  "haloperidol"
)

# Antipruritic proxy: chlorphenamine, or low-dose (40 microgram)
# naloxone. Naloxone formulations are separated in 02_clean_data.R
# (naloxone_type): 40 microgram = pruritus, 400 microgram = respiratory
# reversal.
side_effect_drugs <- c(
  antiemetic_drugs,
  "chlorphenamine",
  "naloxone"
)

# Laxatives (descriptive only; not an analysed outcome)
laxative_drugs <- c(
  "macrogol 3350 + electrolytes",
  "lactulose",
  "senna",
  "docusate sodium",
  "bisacodyl",
  "sodium picosulfate"
)

# ---- Nerve-block categories (str_detect on free-text block_type) ----

neuraxial_blocks <- c("epidural", "spinal", "intrathecal")

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

# ---- OPCS codes for non-substantive procedures ----
# An admission is excluded from the primary comparative population
# only if every code on every operation record is on this list.
opcs_exclude <- c(
  "A54.2",  # Intrathecal chemotherapy
  "A55.9",  # Lumbar puncture
  "E42.1",  # Permanent tracheostomy
  "E42.3",  # Temporary tracheostomy
  "E49.2",  # Flexible bronchoscopy with lavage
  "E49.9",  # Flexible bronchoscopy
  "E85.1",  # Intubation and tracheostomy
  "G45.1",  # OGD with biopsy
  "G45.9",  # OGD
  "G47.8",  # Nasogastric tube insertion
  "G64.3",  # Nasojejunal tube insertion
  "H22.9",  # Colonoscopy
  "L91.1",  # Hickman line insertion
  "L91.2",  # Central line insertion
  "L91.3",  # Central venous catheter removal
  "L91.4",  # Hickman line removal
  "L91.5",  # Long line insertion
  "L91.8",  # Tunnelled venous catheter removal
  "L94.3",  # Portacath insertion
  "L94.8",  # Portacath removal
  "L99.7",  # PICC line insertion
  "M47.9",  # Urinary catheter insertion
  "O56.2",  # Portacath removal (coded under chapter O)
  "S43.3",  # Suture removal
  "S43.4",  # Suture removal
  "T34.9",  # Abdominal drain insertion
  "U05.2",  # MRI head
  "U05.4",  # CT lumbar spine
  "U05.5",  # MRI spine
  "U21.1",  # MRI
  "W36.5",  # Bone marrow aspirate/trephine
  "X56.8",  # Trial of extubation
  "X59.9"   # Operation not found
)
