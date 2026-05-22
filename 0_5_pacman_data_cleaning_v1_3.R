# =============================================================================
# Script 00.5 v1.3: FoxInsight Data Cleaning & Metadata Assembly
# =============================================================================
# Changes from v1.2:
#   - amantadine added as derived binary covariate from PDMBCurrMedSym (999 → NA)
#   - PDMBCurrMedSym added to source variable traceability columns
#   - amantadine added to missingness report
# =============================================================================
# Purpose:
#   Assemble a single clean metadata file from three FoxInsight source tables
#   (PDMB, General, PD Length) and compute all derived variables needed
#   downstream. Output is used by Script 01 in place of the Stagman TSV.
#
# Inputs:
#   foxinsight_pdmb.csv      — PDMB survey (one row per participant; anchor file)
#   foxinsight_general.csv   — General survey (longitudinal; multiple rows per participant)
#   foxinsight_pd_length.csv — PD duration survey (longitudinal; multiple rows per participant)
#
# Outputs:
#   foxinsight_metadata_clean.csv   — For manual inspection
#   foxinsight_metadata_clean.RData — For Script 01 pipeline (contains: meta_clean)
#
# Derived variables computed here:
#   BMI              — from HeightCm + WeightKgs (general, closest visit to PDMB age)
#   disease_duration — current_age - InitPDDiagAge (from pd_length)
#   Bristol_Avg      — weighted mean stool type from PDMBStool1–7 frequency columns
#   azilect          — binary 0/1 from PDMBCurrMedAzl (999 → NA)
#   amantadine       — binary 0/1 from PDMBCurrMedSym (999 → NA)
#   levodopa         — binary 0/1 from PDMBCurrMedSinam + PDMBCurrMedParc (either = 1)
#   parkinsons       — "PD case" / "Control" from CurrPDDiag (general)
#   sex              — "Male" / "Female" from PDMBSex
#   batch            — "Main" / "Pilot" from PDMBStudy
#
# Antibiotic exclusion:
#   PDMBAntiBio == 1 (confirmed use) → excluded from output
#   PDMBAntiBio == 2 (unsure)        → retained (conservative inclusion)
#   PDMBAntiBio == 999               → treated as missing, retained
#
# Baseline selection for longitudinal files:
#   Priority 1: Row where schedule_of_activities == "REG" (confirmed baseline)
#   Priority 2: Row with lowest days_elapsed (assumed baseline; schedule_interpolated = 1)
#   This mirrors the logic in Script 01 v1.3 Section 8.3.4.
# =============================================================================

# --- 1. SETUP ---
.libPaths("/arc/project/st-silkec-1/rsandboxlib")
setwd("/arc/project/st-silkec-1/pacman")

library(dplyr)
library(tidyr)

META_DIR    <- "/arc/project/st-silkec-1/pacman/microbiome_metadata"

# --- 2. LOAD SOURCE FILES ---
pdmb      <- read.csv(file.path(META_DIR, "foxinsight_pdmb.csv"),      stringsAsFactors = FALSE)
general   <- read.csv(file.path(META_DIR, "foxinsight_general.csv"),    stringsAsFactors = FALSE)
pd_length <- read.csv(file.path(META_DIR, "foxinsight_pd_length.csv"),  stringsAsFactors = FALSE)

cat("--- SOURCE FILE ROW COUNTS ---\n")
cat("PDMB:      ", nrow(pdmb),      "rows,", n_distinct(pdmb$fox_insight_id),      "unique IDs\n")
cat("General:   ", nrow(general),   "rows,", n_distinct(general$fox_insight_id),   "unique IDs\n")
cat("PD length: ", nrow(pd_length), "rows,", n_distinct(pd_length$fox_insight_id), "unique IDs\n")

# Confirm ID column is present in all three
stopifnot("fox_insight_id missing from pdmb"      = "fox_insight_id" %in% names(pdmb))
stopifnot("fox_insight_id missing from general"   = "fox_insight_id" %in% names(general))
stopifnot("fox_insight_id missing from pd_length" = "fox_insight_id" %in% names(pd_length))

# =============================================================================
# SECTION 3: ANTIBIOTIC EXCLUSION
# Applied first — excluded participants never enter the join.
# Only PDMBAntiBio == 1 (confirmed use in prior 6 months) is excluded.
# PDMBAntiBio == 2 (unsure) is retained.
# =============================================================================

n_before_abx <- nrow(pdmb)
pdmb <- pdmb %>% filter(PDMBAntiBio != 1 | is.na(PDMBAntiBio))
n_after_abx  <- nrow(pdmb)

cat("\n--- ANTIBIOTIC EXCLUSION ---\n")
cat("Excluded (confirmed antibiotic use, PDMBAntiBio == 1):",
    n_before_abx - n_after_abx, "\n")
cat("Retained after exclusion:", n_after_abx, "\n")
cat("  of which unsure (PDMBAntiBio == 2):",
    sum(pdmb$PDMBAntiBio == 2, na.rm = TRUE), "\n")

# Create antibiotic certainty flag for sensitivity analysis:
#   abx_certain_no = 1 → confirmed no antibiotic use (PDMBAntiBio == 0)
#   abx_certain_no = 0 → unsure (PDMBAntiBio == 2) or missing/skipped
# Sensitivity analysis in downstream scripts can re-run restricted to
# abx_certain_no == 1 to confirm results are robust to the unsure group.
pdmb <- pdmb %>%
  mutate(abx_certain_no = case_when(
    PDMBAntiBio == 0   ~ 1L,   # confirmed no use
    PDMBAntiBio == 2   ~ 0L,   # unsure — retained but flagged
    PDMBAntiBio == 999 ~ NA_integer_,
    TRUE               ~ NA_integer_
  ))

cat("\nabx_certain_no flag (1=confirmed no use, 0=unsure, NA=missing):\n")
print(table(pdmb$abx_certain_no, useNA = "ifany"))

# =============================================================================
# SECTION 4: BASELINE SELECTION FOR LONGITUDINAL FILES
# Select one row per participant from general and pd_length using:
#   Priority 1: schedule_of_activities == "REG" (confirmed baseline)
#   Priority 2: earliest days_elapsed (assumed baseline)
# =============================================================================

select_baseline <- function(df, label) {
  stopifnot(
    "schedule_of_activities column missing" = "schedule_of_activities" %in% names(df),
    "days_elapsed column missing"           = "days_elapsed" %in% names(df)
  )
  
  reg_rows <- df %>%
    filter(trimws(as.character(schedule_of_activities)) == "REG") %>%
    distinct(fox_insight_id, .keep_all = TRUE) %>%
    mutate(schedule_interpolated = 0L)
  
  ids_with_reg <- unique(reg_rows$fox_insight_id)
  
  fallback_rows <- df %>%
    filter(!fox_insight_id %in% ids_with_reg) %>%
    group_by(fox_insight_id) %>%
    arrange(days_elapsed) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(schedule_interpolated = 1L)
  
  out <- bind_rows(reg_rows, fallback_rows)
  
  cat("\n---", label, "baseline selection ---\n")
  cat("  Confirmed baseline (REG):            ", sum(out$schedule_interpolated == 0), "\n")
  cat("  Assumed baseline (days_elapsed):     ", sum(out$schedule_interpolated == 1), "\n")
  cat("  Total rows carried forward:          ", nrow(out), "\n")
  
  out
}

# Filter to PDMB participants only before baseline selection (speeds up processing)
pdmb_ids <- unique(pdmb$fox_insight_id)

general_baseline <- general %>% filter(fox_insight_id %in% pdmb_ids) %>%
  select_baseline("General")

# InitPDDiagAge is sparsely recorded — participants fill it in at one visit only,
# not necessarily the first. Bidirectional fill propagates the value to all rows
# before age-matching, ensuring the closest-age row has the value available.
pd_length_baseline <- pd_length %>%
  filter(fox_insight_id %in% pdmb_ids) %>%
  group_by(fox_insight_id) %>%
  arrange(age) %>%
  fill(InitPDDiagAge, .direction = "downup") %>%
  left_join(pdmb %>% select(fox_insight_id, PDMBAgeYrs), by = "fox_insight_id") %>%
  mutate(age_diff = abs(age - PDMBAgeYrs)) %>%
  slice_min(age_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-age_diff, -age)

cat("\n--- PD length baseline selection ---\n")
cat("  Matched by closest age to PDMB collection\n")
cat("  Total rows carried forward:", nrow(pd_length_baseline), "\n")

# =============================================================================
# SECTION 5: JOIN
# PDMB is the anchor. Left join preserves all post-exclusion PDMB participants.
# schedule_interpolated columns renamed to avoid collision.
# =============================================================================

# Load CurrPDDiag — separate file, one row per participant, no deduplication needed
pd_diag <- read.csv(file.path(META_DIR, "CurrPDDiag.csv"), stringsAsFactors = FALSE)
cat("CurrPDDiag loaded:", nrow(pd_diag), "participants\n")

# Columns to carry from general (drop redundant schedule/days cols after use)
general_keep <- general_baseline %>%
  left_join(pd_diag, by = "fox_insight_id") %>%
  select(fox_insight_id, schedule_interpolated,
         CurrPDDiag,
         HeightCm, WeightKgs, WeightLbs,
         Sex, SexAtBirth)

# Columns to carry from pd_length
pd_length_keep <- pd_length_baseline %>%
  select(fox_insight_id,
         InitPDDiag, InitPDDiagAge, AgeAtEnrollment, YearsWithPD,
         WeightKgs_pdl = WeightKgs, WeightLbs_pdl = WeightLbs)

combined <- pdmb %>%
  left_join(general_keep,   by = "fox_insight_id") %>%
  left_join(pd_length_keep, by = "fox_insight_id")

cat("\n--- POST-JOIN ---\n")
cat("Rows:", nrow(combined), "(should equal", n_after_abx, ")\n")
stopifnot("Join inflated rows — duplicate IDs in general or pd_length" =
            nrow(combined) == n_after_abx)

# =============================================================================
# SECTION 6: DERIVED VARIABLES
# =============================================================================

# --- 6.1 Recode sex ---
# PDMBSex: 1 = Male, 2 = Female
combined <- combined %>%
  mutate(sex = case_when(
    PDMBSex == 1 ~ "Male",
    PDMBSex == 2 ~ "Female",
    TRUE         ~ NA_character_
  ))

# --- 6.2 Recode PD status ---
# CurrPDDiag: 1 = PD case, 0 = Control
combined <- combined %>%
  mutate(parkinsons = case_when(
    CurrPDDiag == 1 ~ "PD case",
    CurrPDDiag == 0 ~ "Control",
    TRUE            ~ NA_character_
  ))

cat("\n--- PD STATUS ---\n")
print(table(combined$parkinsons, useNA = "ifany"))

# --- 6.3 Recode batch ---
# PDMBStudy: 1 = Pilot, 2 = Main
combined <- combined %>%
  mutate(batch = case_when(
    PDMBStudy == 1 ~ "Pilot",
    PDMBStudy == 2 ~ "Main",
    TRUE           ~ NA_character_
  ))

# --- 6.4 BMI ---
# Use HeightCm from general and WeightKgs from general (collected at same visit).
# WeightKgs_pdl (from pd_length) used as fallback if general weight is missing.
# BMI = weight(kg) / (height(m))^2
combined <- combined %>%
  mutate(
    weight_kg_final = case_when(
      !is.na(WeightKgs) & WeightKgs > 0              ~ WeightKgs,
      !is.na(WeightKgs_pdl) & WeightKgs_pdl > 0      ~ WeightKgs_pdl,
      TRUE                                            ~ NA_real_
    ),
    height_m = HeightCm / 100,
    BMI = ifelse(
      !is.na(weight_kg_final) & !is.na(height_m) & height_m > 0,
      round(weight_kg_final / height_m^2, 2),
      NA_real_
    )
  )

# Plausibility check — flag implausible BMI values
n_implausible_bmi <- sum(combined$BMI < 12 | combined$BMI > 70, na.rm = TRUE)
if (n_implausible_bmi > 0) {
  warning(n_implausible_bmi, " participant(s) have implausible BMI (<12 or >70) — set to NA")
  combined$BMI[!is.na(combined$BMI) & (combined$BMI < 12 | combined$BMI > 70)] <- NA
}

cat("\n--- BMI SUMMARY ---\n")
print(summary(combined$BMI))
cat("Missing BMI:", sum(is.na(combined$BMI)), "\n")

# --- 6.5 Disease duration ---
# current_age comes from PDMBAgeYrs (age at sample collection)
# InitPDDiagAge comes from pd_length
# Negative values = data entry error → NA
combined <- combined %>%
  mutate(
    current_age      = PDMBAgeYrs,
    disease_duration = PDMBAgeYrs - InitPDDiagAge,
    disease_duration = ifelse(disease_duration < 0, NA_real_, disease_duration)
  )

n_neg_dur <- sum((combined$PDMBAgeYrs - combined$InitPDDiagAge) < 0, na.rm = TRUE)
if (n_neg_dur > 0) {
  cat("WARNING:", n_neg_dur,
      "participant(s) had negative disease duration (age < InitPDDiagAge) — set to NA\n")
}

cat("\n--- DISEASE DURATION SUMMARY (years) ---\n")
print(summary(combined$disease_duration))

# --- 6.6 Medication binaries ---
# 999 = skipped/not answered → NA
# 1 = yes, 0 = no

recode_med <- function(x) {
  case_when(x == 1 ~ 1L, x == 0 ~ 0L, TRUE ~ NA_integer_)
}

combined <- combined %>%
  mutate(
    azilect    = recode_med(PDMBCurrMedAzl),
    amantadine = recode_med(PDMBCurrMedSym),
    # Levodopa: any formulation (Sinemet or Parcopa)
    levodopa = as.integer(recode_med(PDMBCurrMedSinam) == 1 |
                            recode_med(PDMBCurrMedParc)  == 1)
  )

# levodopa NA if both source columns are NA
combined <- combined %>%
  mutate(levodopa = ifelse(
    is.na(PDMBCurrMedSinam) & is.na(PDMBCurrMedParc), NA_integer_, levodopa
  ))

cat("\n--- MEDICATION SUMMARY ---\n")
cat("Azilect (1=yes, 0=no, NA=missing):\n")
print(table(combined$azilect, useNA = "ifany"))
cat("Azilect prevalence among non-missing:",
    round(mean(combined$azilect == 1, na.rm = TRUE) * 100, 1), "%\n")
cat("\nAmantadine (1=yes, 0=no, NA=missing):\n")
print(table(combined$amantadine, useNA = "ifany"))
cat("Amantadine prevalence among non-missing:",
    round(mean(combined$amantadine == 1, na.rm = TRUE) * 100, 1), "%\n")
cat("\nLevodopa (1=yes, 0=no, NA=missing):\n")
print(table(combined$levodopa, useNA = "ifany"))

# --- 6.7 Bristol Average Score ---
# PDMBStool1–7 encode frequency of each Bristol stool type.
# Frequency levels match the Stagman source mapping used in Scripts 01.5 and 02.
# Bristol types 1–7 weight the mean stool consistency.
freq_map    <- c("never" = 0, "< 3x/week" = 1.5, "3-6x/week" = 4.5,
                 "daily" = 7, "2-3x/day" = 17.5, "> 3x/day" = 25)
stool_cols  <- paste0("PDMBStool", 1:7)

# Check for unexpected values
unexpected_bristol <- combined %>%
  select(all_of(stool_cols)) %>%
  summarise(across(everything(),
                   ~sum(!as.character(.) %in% c(names(freq_map), NA, ""), na.rm = TRUE)))
if (any(unexpected_bristol > 0)) {
  warning("Unexpected values in Bristol stool columns — check coding:\n")
  print(unexpected_bristol[, unexpected_bristol > 0, drop = FALSE])
}

bristol_mat <- combined %>%
  mutate(across(all_of(stool_cols),
                ~freq_map[as.character(.)],
                .names = "n_{.col}")) %>%
  select(starts_with("n_PDMBStool")) %>%
  as.matrix()

w_sum   <- as.numeric(bristol_mat %*% 1:7)
total_f <- rowSums(bristol_mat, na.rm = TRUE)

combined <- combined %>%
  mutate(Bristol_Avg = ifelse(total_f > 0, w_sum / total_f, NA_real_))

cat("\n--- BRISTOL AVERAGE SCORE ---\n")
print(summary(combined$Bristol_Avg))
cat("Missing Bristol_Avg:", sum(is.na(combined$Bristol_Avg)), "\n")

# =============================================================================
# SECTION 7: SELECT AND RENAME FINAL COLUMNS
# Keep only columns needed downstream + key source columns for traceability.
# =============================================================================

meta_clean <- combined %>%
  mutate(FoxDEN_ID = fox_insight_id) %>%
  select(
    # Identity
    FoxDEN_ID, fox_insight_id,
    # Key grouping variables
    parkinsons, sex, batch,
    # Demographics
    current_age, BMI,
    # Disease
    disease_duration, InitPDDiagAge, YearsWithPD,
    # Medications
    azilect, amantadine, levodopa,
    # Gut
    Bristol_Avg,
    # Antibiotic flags
    # PDMBAntiBio: raw source (0=no, 1=yes [excluded], 2=unsure, 999=skipped)
    # abx_certain_no: 1=confirmed no use, 0=unsure — use for sensitivity analysis
    PDMBAntiBio, abx_certain_no,
    # Baseline provenance
    schedule_interpolated,
    # Source variables retained for traceability
    PDMBAgeYrs, HeightCm, WeightKgs, WeightKgs_pdl,
    CurrPDDiag, PDMBSex, PDMBStudy,
    PDMBCurrMedAzl, PDMBCurrMedSym, PDMBCurrMedSinam, PDMBCurrMedParc,
    # All Bristol source columns
    all_of(stool_cols)
  )

# =============================================================================
# SECTION 8: FINAL REPORT & COMPLETENESS CHECK
# =============================================================================

cat("\n===== FINAL METADATA SUMMARY =====\n")
cat("Total participants (post antibiotic exclusion):", nrow(meta_clean), "\n")
cat("\nBy PD status:\n")
print(table(meta_clean$parkinsons, useNA = "ifany"))
cat("\nBy sex:\n")
print(table(meta_clean$sex, useNA = "ifany"))
cat("\nBy batch:\n")
print(table(meta_clean$batch, useNA = "ifany"))
cat("\nBaseline assignment (general):\n")
print(table(meta_clean$schedule_interpolated,
            dnn = "schedule_interpolated (0=REG, 1=fallback)"))

cat("\n--- MISSINGNESS FOR KEY ANALYSIS VARIABLES ---\n")
key_vars <- c("parkinsons", "sex", "BMI", "current_age", "disease_duration",
              "azilect", "amantadine", "levodopa", "Bristol_Avg", "abx_certain_no")
miss_summary <- sapply(key_vars, function(v) sum(is.na(meta_clean[[v]])))
print(data.frame(variable = names(miss_summary), n_missing = miss_summary,
                 pct_missing = round(miss_summary / nrow(meta_clean) * 100, 1),
                 row.names = NULL))

# =============================================================================
# SECTION 9: EXPORT
# =============================================================================

write.csv(meta_clean,
          file.path(META_DIR, "foxinsight_metadata_clean.csv"),
          row.names = FALSE, na = "")

save(meta_clean, file = file.path(META_DIR, "foxinsight_metadata_clean.RData"))

cat("\n--- SCRIPT 00.5 COMPLETE ---\n")
cat("Outputs saved:\n")
cat("  CSV:   microbiome_metadata/foxinsight_metadata_clean.csv\n")
cat("  RData: microbiome_metadata/foxinsight_metadata_clean.RData\n")
cat("Proceed to Script 01 v1.9 which loads foxinsight_metadata_clean.RData\n")
cat("in place of the Stagman TSV.\n")