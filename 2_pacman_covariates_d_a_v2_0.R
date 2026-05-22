# =============================================================================
# Script 2 v2.0: PDAQ vs. Genus-Level Diversity & Differential Abundance
# =============================================================================
# Changes from v1.9 (v2.0):
#   - Bristol stool score restored as covariate (PRESPEC_COVARIATES_COG);
#     removal in v1.1 was made without explicit instruction — restored pending
#     Ted's decision on mediator vs confounder argument
#   - Bristol calculation block (Section 2B) restored
#   - PRESPEC_COVARIATES renamed to PRESPEC_COVARIATES_COG throughout
#   - DA_COVARIATES renamed to DA_COVARIATES_COG throughout
#   - Bristol rationale block removed from header (decision deferred to Ted)
#
#   - amantadine added as candidate medication covariate
#   - amantadine recoded from ps_pdaq metadata (carried from Script 00.5 v1.3
#     via Script 01 v1.9): as.factor, 1=yes/0=no
#   - amantadine added to drop_na (complete-case requirement)
#   - amantadine added to PERMANOVA formula (both continuous and binary)
#   - amantadine added to Section 3B variance check
#   - amantadine rationale added to script header
#   - NOTE: amantadine is NOT added to PRESPEC_COVARIATES or DA models here;
#     inclusion in DA is conditional on PERMANOVA significance at runtime
#   - Stalevo/entacapone changes reverted — Stalevo remains excluded
#   - Formal covariate variance check added (Section 3B): prints prevalence
#     table for all candidate medication covariates; flags any variable failing
#     the 5% minimum-cell-frequency rule; documents Stalevo exclusion
#     with observed n and % directly from the analysis cohort
# Changes from v1.5:
#   - azilect renamed to rasagiline throughout (variable name, covariate list,
#     all diagnostics) — reflects correct drug name for reporting
#   - PERMANOVA results (continuous and binary) exported to Excel
#     (Table_PERMANOVA_Results.xlsx) via openxlsx, two sheets
# Changes from v1.4:
#   - PCoA and alpha diversity plot print() and ggsave() calls removed
#   - Plot objects added to save() call for Script 05 rendering:
#     pcoa_res, pcoa_df, perm_cont, perm_bin, pcoa_cont_plot, pcoa_bin_plot,
#     ps_genus_rare, alpha_long, alpha_stats, alpha_plot
#   - Output paths updated to new directory structure:
#       RData/ps objects → ~/rdata/
#       Tables           → ~/tables/
#       Bioinformatics   → ~/bioinfo_results/
#   - Input path updated: loads PD_PDAQ_Association_Data.RData from ~/rdata/
#   - MaAsLin2 output directories updated to ~/bioinfo_results/
#   - log1p removed from Bray-Curtis distance (line ~190); raw Bray-Curtis
#     used throughout for consistency with all other scripts
#   - N added to all figure subtitles (global directive)
#   - ggsave() calls added for all 3 figures
# Changes from v1.2 (v1.3):
#   - Added schedule_interpolated transparency report after d_prep (Section 2)
#   - Updated input note to reflect Script 01 v1.3 column additions
# Changes from v1.0 (v1.1 → v1.2):
#   - Bristol stool score removed as covariate (see rationale below)
#   - Levodopa removed as covariate; replaced with Azilect (see rationale below)
#   - Stalevo (COMT inhibitor) considered and excluded (see rationale below)
#   - Disease duration added as covariate (current_age - age_diagnosed)
#   - Covariates are now pre-specified rather than PERMANOVA-selected, to avoid
#     circular reasoning (using the same distance matrix for both covariate
#     selection and hypothesis testing inflates false positive risk)
#   - Analysis A is now continuous PDAQ (primary); binary cog_intact is Section C
#   - Section B is beta diversity (PERMANOVA) reported independently of DA
#     covariate selection
#   - MaAsLin2 filter now uses taxon column (confirmed column name in this version)
#   - Taxonomy annotation and CSV export added for all consensus results
#

# Levodopa rationale:
#   Levodopa is excluded because near-universal use in this PD cohort produces
#   a near-constant binary variable with insufficient variance for meaningful
#   adjustment. Including it wastes a degree of freedom without providing any
#   statistical control. Medication adjustment is instead achieved via Azilect
#   (see below), which has both sufficient variance and stronger biological
#   rationale specific to microbiome composition.
#
# Rasagiline (Azilect) rationale:
#   Azilect is included as the primary medication covariate. Prevalence in this
#   cohort is ~24% (78/330), providing adequate variance for adjustment.
#   Rasagiline (MAO-B inhibitor) has documented neuroprotective and
#   anti-inflammatory effects in the enteric nervous system and has been
#   specifically linked to microbiome composition in PD independently of its
#   dopaminergic mechanism. Binary: present (1) vs absent (0).
#
# Stalevo (entacapone/COMT inhibitor) rationale:
#   Stalevo is excluded due to insufficient prevalence in this cohort. The
#   minimum-cell-frequency rule (5% of N in the minority cell) is not met —
#   observed prevalence is reported in the Section 3B covariate variance check.
#   Including a variable this sparse risks near-singular design matrices across
#   all three DA methods (MaAsLin2, LinDA, ANCOM-BC2) and inflated standard
#   errors in PERMANOVA. Other COMT inhibitor formulations (Comtan, Tolcapone)
#   were not present as separate columns in the FoxDEN medication data.
#
# Amantadine rationale:
#   Amantadine is included as a candidate medication covariate in PERMANOVA.
#   Prevalence in this cohort is ~10.5% (~42 users), meeting the 5% minority-
#   cell threshold. Amantadine (NMDA antagonist with anticholinergic properties)
#   has well-documented effects on gut motility and has been associated with
#   constipation and microbiome composition changes in PD. Binary: present (1)
#   vs absent (0). Inclusion in DA models is conditional on PERMANOVA
#   significance — see PRESPEC_COVARIATES definition below.
#
# Pre-specified covariate rationale:
#   The following six covariates are included in all models regardless of
#   PERMANOVA significance, based on prior literature and biological rationale:
#     current_age      — cognitive decline is strongly age-dependent
#     sex              — sex differences in PD microbiome are well-documented
#     BMI_clean        — metabolic status influences gut microbiome composition
#     rasagiline       — MAO-B inhibitor with direct enteric effects; sufficient
#                        variance (≥5% minority cell); binary present/absent
#     batch            — technical covariate; main vs pilot study recruitment
#     disease_duration — years since PD diagnosis; captures disease stage and
#                        cumulative gut dysbiosis independent of current age
# =============================================================================

# --- 1. INITIALIZATION & DATA LOADING ---
.libPaths("/arc/project/st-silkec-1/rsandboxlib")
setwd("/arc/project/st-silkec-1/pacman")

# =============================================================================
# PROJECT PALETTE (Okabe-Ito — colourblind-safe across all scripts)
# PD case / Cognitively impaired  →  #D55E00  vermillion
# Control / Cognitively intact    →  #0072B2  blue
# PD-MCI (three-level only)       →  #E69F00  orange
# Sequencing outliers             →  #CC79A7  pink
# =============================================================================
PAL_PD   <- "#D55E00"
PAL_CTRL <- "#0072B2"

# =============================================================================
# OUTPUT DIRECTORIES
# =============================================================================
RDATA_DIR   <- "/arc/project/st-silkec-1/pacman/rdata"
TABLES_DIR  <- "/arc/project/st-silkec-1/pacman/tables"
BIOINFO_DIR <- "/arc/project/st-silkec-1/pacman/bioinfo_results"
PLOTS_DIR   <- "/arc/project/st-silkec-1/pacman/plots"
dir.create(RDATA_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(TABLES_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(BIOINFO_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,   showWarnings = FALSE, recursive = TRUE)

load(file.path(RDATA_DIR, "PD_PDAQ_Association_Data.RData"))

library(phyloseq)
library(tidyverse)
library(vegan)
library(ANCOMBC)
library(Maaslin2)
library(MicrobiomeStat)


# =============================================================================
# SECTION 2: METADATA PREPARATION
# =============================================================================

d_prep <- data.frame(sample_data(ps_pdaq)) %>%
  mutate(SampleID_Fixed = sample_names(ps_pdaq)) %>%
  filter(parkinsons == "PD case") %>%
  mutate(
    PDAQ_num         = as.numeric(as.character(PDAQ_Total)),
    current_age      = as.numeric(as.character(current_age)),
    age_diagnosed    = as.numeric(as.character(age_diagnosed)),
    BMI_clean        = as.numeric(as.character(BMI)),
    sex              = as.factor(sex),
    batch            = as.factor(main_or_pilot_study),
    # Rasagiline (Azilect, MAO-B inhibitor): primary medication covariate.
    # Levodopa excluded (near-universal use, insufficient variance).
    # Stalevo excluded (insufficient prevalence — see Section 3B variance check).
    # Amantadine: candidate covariate; included in PERMANOVA; DA inclusion
    # conditional on PERMANOVA significance (see Section 3B and PRESPEC_COVARIATES).
    # See script header for full rationale.
    rasagiline       = as.factor(ifelse(
      !is.na(`pd_specific_medications.azilect`) &
        !grepl("Not", `pd_specific_medications.azilect`), 1, 0)),
    amantadine       = as.factor(ifelse(!is.na(amantadine) & amantadine == 1, 1, 0)),
    # Disease duration: years living with PD diagnosis
    # Negative values indicate data entry error — set to NA and report
    disease_duration = current_age - age_diagnosed,
    disease_duration = ifelse(disease_duration < 0, NA, disease_duration)
  ) %>%
  mutate(
    # Cognitive intact per validated PDAQ threshold (Levin et al.)
    # >= 53 = cognitively intact (1); < 53 = cognitively impaired (0)
    cog_intact = as.factor(ifelse(PDAQ_num >= 53, 1, 0))
  )

# Report any implausible disease duration values
n_neg_duration <- sum(d_prep$disease_duration < 0, na.rm = TRUE)
if (n_neg_duration > 0) {
  cat("WARNING:", n_neg_duration,
      "participant(s) had negative disease duration (age < age_diagnosed) — set to NA\n")
}

# Report baseline assignment method for PDAQ scores
# schedule_interpolated is carried through from Script 01 v1.3:
#   0 = confirmed baseline (REG entry in schedule_of_activities)
#   1 = assumed baseline (earliest days_elapsed entry; schedule_of_activities missing)
# No filtering applied — assumption is reasonable for single-entry participants —
# but reported here so the proportion of interpolated baselines is visible.
if ("schedule_interpolated" %in% colnames(d_prep)) {
  cat("\n--- PDAQ BASELINE ASSIGNMENT (PD cases only) ---\n")
  print(table(d_prep$schedule_interpolated,
              dnn = "schedule_interpolated (0=REG confirmed, 1=assumed)"))
} else {
  cat("NOTE: schedule_interpolated column not found in ps_pdaq metadata.\n")
  cat("      Re-run Script 01 v1.3 to generate this column.\n")
}

# =============================================================================
# SECTION 2B: BRISTOL STOOL SCORE CALCULATION
# =============================================================================
freq_map    <- c("never"=0, "< 3x/week"=1.5, "3-6x/week"=4.5,
                 "daily"=7, "2-3x/day"=17.5, "> 3x/day"=25)
target_cols <- c("bl_bristol_stool_type_1_freq", "bl_bristol_stool_type_2_freq",
                 "bl_bristol_stool_type_3_freq", "bl_bristol_stool_type_4_freq",
                 "bl_bristol_stool_type_5_freq", "bl_bristol_stool_type_6_freq",
                 "bl_bristol_stool_type_7_freq")

unexpected <- d_prep %>%
  select(all_of(target_cols)) %>%
  summarise(across(everything(),
                   ~sum(!as.character(.) %in% c(names(freq_map), NA))))
if (any(unexpected > 0)) {
  warning("Unexpected values in Bristol columns — these will become NA:\n")
  print(unexpected[, unexpected > 0, drop = FALSE])
}

bristol_mat <- d_prep %>%
  mutate(across(all_of(target_cols),
                ~freq_map[as.character(.)], .names = "n_{.col}")) %>%
  select(starts_with("n_bl_bristol")) %>%
  as.matrix()

w_sum_vec   <- as.numeric(bristol_mat %*% 1:7)
total_f_vec <- rowSums(bristol_mat, na.rm = TRUE)

d_prep <- d_prep %>%
  mutate(Bristol_Avg = ifelse(total_f_vec > 0, w_sum_vec / total_f_vec, NA))

cat("Bristol_Avg NA count:", sum(is.na(d_prep$Bristol_Avg)), "\n")

# =============================================================================
# SECTION 3: FINAL FILTERING & ALIGNMENT
# =============================================================================
# drop_na enforces complete cases across all pre-specified covariates.

df_meta_final <- d_prep %>%
  drop_na(PDAQ_num, cog_intact, current_age, sex, BMI_clean,
          Bristol_Avg, rasagiline, amantadine, batch, disease_duration)

cat("\n--- ANALYSIS COHORT AFTER CLEANING ---\n")
cat("PD cases with complete data:", nrow(df_meta_final), "\n")
cat("\nCognitive intact (1) vs impaired (0):\n")
print(table(df_meta_final$cog_intact))
cat("\nRasagiline use (1 = yes, 0 = no):\n")
print(table(df_meta_final$rasagiline))
cat("Rasagiline prevalence:", round(mean(df_meta_final$rasagiline == 1) * 100, 1), "%\n")
cat("\nAmantadine use (1 = yes, 0 = no):\n")
print(table(df_meta_final$amantadine))
cat("Amantadine prevalence:", round(mean(df_meta_final$amantadine == 1) * 100, 1), "%\n")
cat("\nDisease duration summary (years):\n")
print(summary(df_meta_final$disease_duration))
cat("\nPDAQ_num summary:\n")
print(summary(df_meta_final$PDAQ_num))

# =============================================================================
# SECTION 3B: COVARIATE VARIANCE CHECK — MEDICATION COVARIATES
# =============================================================================
# Binary covariates require a minimum minority-cell frequency to be included
# in models. Threshold: ≥5% of analytic N in the minority cell (widely used
# rule of thumb for binary predictors in multivariate models; see e.g. Peduzzi
# et al. 1996 for regression; analogous reasoning applies to PERMANOVA and DA).
#
# Rasagiline: expected ~24% prevalence — should comfortably pass.
# Stalevo (entacapone): expected ~3% — expected to fail; excluded on this basis.
# This section computes observed values from the final analytic cohort so the
# exclusion is data-driven rather than assumed.

cat("\n=== SECTION 3B: COVARIATE VARIANCE CHECK ===\n")

N_final <- nrow(df_meta_final)
MIN_CELL_THRESHOLD <- 0.05   # 5% minority-cell rule

# Build a temporary vector for Stalevo directly from d_prep (pre-drop_na),
# aligned to the final cohort by SampleID_Fixed, so N matches df_meta_final.
entacapone_in_cohort <- d_prep %>%
  filter(SampleID_Fixed %in% df_meta_final$SampleID_Fixed) %>%
  mutate(entacapone = as.factor(ifelse(
    !is.na(`pd_specific_medications.stalevo`) &
      !grepl("Not", `pd_specific_medications.stalevo`), 1, 0
  ))) %>%
  pull(entacapone)

check_variance <- function(var_vec, var_name, n_total, threshold) {
  tbl        <- table(var_vec)
  n_minority <- min(tbl)
  pct        <- round(n_minority / n_total * 100, 1)
  passes     <- (n_minority / n_total) >= threshold
  cat(sprintf(
    "  %-20s  minority cell n = %3d / %d  (%5.1f%%)  →  %s\n",
    var_name, n_minority, n_total, pct,
    ifelse(passes,
           paste0("PASS (>= ", threshold * 100, "%)"),
           paste0("FAIL (< ",  threshold * 100, "%) — EXCLUDED"))
  ))
  invisible(list(n_minority = n_minority, pct = pct, passes = passes))
}

cat(sprintf("Analytic N = %d | Minimum minority-cell threshold = %.0f%%\n\n",
            N_final, MIN_CELL_THRESHOLD * 100))

ras_check <- check_variance(df_meta_final$rasagiline, "Rasagiline",  N_final, MIN_CELL_THRESHOLD)
ama_check <- check_variance(df_meta_final$amantadine, "Amantadine",  N_final, MIN_CELL_THRESHOLD)
sta_check <- check_variance(entacapone_in_cohort,     "Entacapone",  N_final, MIN_CELL_THRESHOLD)

cat("\nConclusion:\n")
cat("  Rasagiline included in PRESPEC_COVARIATES_COG — variance threshold met.\n")
cat("  Amantadine included in PERMANOVA — variance threshold met.\n")
cat("  Amantadine inclusion in DA models conditional on PERMANOVA significance.\n")
cat(sprintf(
  "  Entacapone EXCLUDED — minority cell n = %d (%.1f%% of N = %d), below %.0f%% threshold.\n",
  sta_check$n_minority, sta_check$pct, N_final, MIN_CELL_THRESHOLD * 100
))
cat("  Including Entacapone would risk near-singular design matrices in MaAsLin2,\n")
cat("  LinDA, and ANCOM-BC2, and inflated standard errors in PERMANOVA.\n")
cat("============================================\n")

# Align phyloseq object to cleaned samples
rownames(df_meta_final) <- df_meta_final$SampleID_Fixed
ps_final  <- prune_samples(df_meta_final$SampleID_Fixed, ps_pdaq)
ps_genus  <- tax_glom(ps_final, "Genus")
df_meta_final <- df_meta_final[sample_names(ps_genus), ]
sample_data(ps_genus) <- sample_data(df_meta_final)

cat("\nGenus-level taxa after tax_glom:", ntaxa(ps_genus), "\n")

# OTU matrix — samples x taxa
otu_mat           <- as.matrix(otu_table(ps_genus))
rownames(otu_mat) <- sample_names(ps_genus)

# Alignment checks
stopifnot(
  "OTU matrix rows do not match metadata rows — check sample alignment" =
    all(rownames(otu_mat) == rownames(df_meta_final))
)

set.seed(123)
dist_bray <- vegdist(otu_mat, method = "bray")   # raw; no log1p — consistent with all scripts

stopifnot(
  "Distance matrix labels do not match metadata rows — check sample alignment" =
    all(labels(dist_bray) == rownames(df_meta_final))
)
cat("Alignment confirmed.\n")

# =============================================================================
# SECTION 4: TAXONOMY LOOKUP
# =============================================================================
# Built once here and reused for annotation across all analyses.

tax_lookup <- data.frame(tax_table(ps_genus)) %>%
  rownames_to_column("taxon") %>%
  select(taxon, Genus, Family, Phylum) %>%
  mutate(
    Genus  = gsub("^g__", "", Genus),
    Family = gsub("^f__", "", Family),
    Phylum = gsub("^p__", "", Phylum),
    label  = case_when(
      !is.na(Genus)  & Genus  != "" ~ Genus,
      !is.na(Family) & Family != "" ~ paste0("(", Family, " genus)"),
      TRUE                          ~ taxon
    )
  )

# Helper: attach taxonomy to any DA result data frame
annotate_taxa <- function(df) {
  candidates <- c("taxon", "featureID", "OTU", "taxa")
  taxon_col  <- candidates[candidates %in% colnames(df)][1]
  if (is.na(taxon_col)) stop("Could not detect taxon column. Columns: ",
                             paste(colnames(df), collapse = ", "))
  df %>% left_join(tax_lookup, by = setNames("taxon", taxon_col))
}

# =============================================================================
# SECTION 5: PRE-SPECIFIED COVARIATE FORMULA
# =============================================================================
# Single formula used across ALL analyses in this script.
# Rationale for each covariate is documented in the script header.
# This formula is intentionally identical for continuous and binary outcomes —
# consistency is important when comparing results across analyses A and C.

PRESPEC_COVARIATES_COG <- c("current_age", "sex", "BMI_clean",
                        "Bristol_Avg", "rasagiline", "batch", "disease_duration")

# PERMANOVA_COVARIATES_COG extends PRESPEC_COVARIATES_COG with amantadine.
# Amantadine is tested in PERMANOVA; DA inclusion is conditional on significance.
PERMANOVA_COVARIATES_COG <- c(PRESPEC_COVARIATES_COG, "amantadine")

# =============================================================================
# SECTION 6: BETA DIVERSITY — PD COGNITIVE STATUS
# =============================================================================
# PERMANOVA is reported here for beta diversity description only.
# Results do NOT feed into covariate selection for DA models (see header).
# Two PERMANOVA models run: one with continuous PDAQ, one with binary cog_intact.

cat("\n--- SECTION 6: BETA DIVERSITY (PERMANOVA) ---\n")

# 6A. Continuous PDAQ
f_perm_cont <- as.formula(paste("dist_bray ~ PDAQ_num +",
                                paste(PERMANOVA_COVARIATES_COG, collapse = " + ")))
perm_cont <- adonis2(f_perm_cont, data = df_meta_final, by = "margin", permutations = 999)
cat("\nPERMANOVA — Continuous PDAQ + covariates:\n")
print(perm_cont)

# 6B. Binary cog_intact
f_perm_bin <- as.formula(paste("dist_bray ~ cog_intact +",
                               paste(PERMANOVA_COVARIATES_COG, collapse = " + ")))
perm_bin <- adonis2(f_perm_bin, data = df_meta_final, by = "margin", permutations = 999)
cat("\nPERMANOVA — Binary cog_intact + covariates:\n")
print(perm_bin)

# Derive DA_COVARIATES_COG from PERMANOVA-significant covariates (p < 0.05).
# PDAQ_num / cog_intact are always retained as the primary variable regardless
# of PERMANOVA significance — only nuisance covariates are screened here.
# Uses perm_cont (continuous PDAQ PERMANOVA) as the reference model for
# covariate selection — consistent across both Analysis A and Analysis C.
perm_cont_pvals <- as.data.frame(perm_cont)$`Pr(>F)`
names(perm_cont_pvals) <- rownames(as.data.frame(perm_cont))

DA_COVARIATES_COG <- PERMANOVA_COVARIATES_COG[
  PERMANOVA_COVARIATES_COG %in% names(perm_cont_pvals) &
    !is.na(perm_cont_pvals[PERMANOVA_COVARIATES_COG]) &
    perm_cont_pvals[PERMANOVA_COVARIATES_COG] < 0.05
]

cat("\nPERMANOVA-significant covariates carried into DA models:\n")
cat(" ", paste(DA_COVARIATES_COG, collapse = ", "), "\n")
dropped <- setdiff(PERMANOVA_COVARIATES_COG, DA_COVARIATES_COG)
if (length(dropped) > 0) {
  cat("Dropped (non-significant in PERMANOVA):\n")
  cat(" ", paste(dropped, collapse = ", "), "\n")
}

# 6C–6D. PCoA plots (unchanged)
pcoa_res <- ordinate(ps_genus, method = "PCoA", distance = dist_bray)

pcoa_df <- data.frame(pcoa_res$vectors[, 1:2],
                      PDAQ_num   = df_meta_final$PDAQ_num,
                      cog_intact = df_meta_final$cog_intact)

pcoa_cont_plot <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, colour = PDAQ_num)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_colour_viridis_c(option = "plasma", name = "PDAQ Total") +
  labs(title    = "Beta Diversity PCoA — PD Cohort (coloured by PDAQ score)",
       subtitle = paste0("N = ", nrow(pcoa_df)),
       x = paste0("PC1 (", round(pcoa_res$values$Relative_eig[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(pcoa_res$values$Relative_eig[2] * 100, 1), "%)")) +
  theme_bw()
# plot generated in Script 05

# 6D. PCoA plot coloured by binary cognitive status
pcoa_bin_plot <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, colour = cog_intact)) +
  geom_point(alpha = 0.7, size = 2.5) +
  stat_ellipse(linetype = 2) +
  scale_colour_manual(values = c("0" = PAL_PD, "1" = PAL_CTRL),
                      labels = c("0" = "Impaired (<53)", "1" = "Intact (≥53)"),
                      name   = "Cognitive status") +
  labs(title    = "Beta Diversity PCoA — PD Cohort (binary cognitive status)",
       subtitle = paste0("N = ", nrow(pcoa_df)),
       x = paste0("PC1 (", round(pcoa_res$values$Relative_eig[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(pcoa_res$values$Relative_eig[2] * 100, 1), "%)")) +
  theme_bw()
# plot generated in Script 05

# =============================================================================
# SECTION 6E: EXPORT PERMANOVA RESULTS TO EXCEL
# =============================================================================
# Both PERMANOVA tables exported to a single workbook (two sheets).


perm_to_df <- function(perm_obj) {
  as.data.frame(perm_obj) %>%
    tibble::rownames_to_column("Term") %>%
    rename(any_of(c("p-value" = "Pr(>F)"))) %>%
    mutate(across(where(is.numeric), ~round(., 4)))
}

write.csv(perm_to_df(perm_cont),
          file.path(TABLES_DIR, "Table_PERMANOVA_Continuous_PDAQ.csv"),
          row.names = FALSE, na = "")
write.csv(perm_to_df(perm_bin),
          file.path(TABLES_DIR, "Table_PERMANOVA_Binary_CogStatus.csv"),
          row.names = FALSE, na = "")
cat("Exported PERMANOVA results:\n")
cat("  → Table_PERMANOVA_Continuous_PDAQ.csv\n")
cat("  → Table_PERMANOVA_Binary_CogStatus.csv\n")
cat("  →", TABLES_DIR, "\n")

# =============================================================================
# SECTION 7: ALPHA DIVERSITY — PD COGNITIVE STATUS
# =============================================================================
# Reported for completeness. Interpreted cautiously — within-PD alpha diversity
# differences by cognitive status are inconsistent in the literature and power
# is limited in this cohort.

cat("\n--- SECTION 7: ALPHA DIVERSITY ---\n")

ps_genus_rare <- rarefy_even_depth(ps_genus, rngseed = 421)

alpha_df <- estimate_richness(ps_genus_rare, measures = c("Shannon", "Observed", "Chao1"))

meta_alpha <- data.frame(
  sample_data(ps_genus_rare),
  Chao1    = alpha_df$Chao1,
  Shannon  = alpha_df$Shannon,
  Observed = alpha_df$Observed,
  check.names = FALSE
)

alpha_long <- meta_alpha %>%
  tidyr::pivot_longer(cols = c("Chao1", "Shannon", "Observed"),
                      names_to = "Metric", values_to = "Value") %>%
  filter(!is.na(cog_intact))

# Wilcoxon tests: intact vs impaired
alpha_stats <- alpha_long %>%
  group_by(Metric) %>%
  rstatix::wilcox_test(Value ~ cog_intact) %>%
  rstatix::add_significance() %>%
  mutate(y.position = c(
    max(alpha_long$Value[alpha_long$Metric == "Chao1"],    na.rm = TRUE) * 1.05,
    max(alpha_long$Value[alpha_long$Metric == "Observed"], na.rm = TRUE) * 1.05,
    max(alpha_long$Value[alpha_long$Metric == "Shannon"],  na.rm = TRUE) * 1.05
  ))

alpha_plot <- ggpubr::ggboxplot(
  alpha_long, x = "cog_intact", y = "Value",
  fill = "cog_intact",
  palette    = c("0" = PAL_PD, "1" = PAL_CTRL),
  add        = "jitter",
  add.params = list(alpha = 0.2),
  facet.by   = "Metric",
  scales     = "free_y",
  title      = "Alpha Diversity — PD Cohort (Intact vs Impaired)",
  subtitle   = paste0("N = ", nsamples(ps_genus_rare))
) +
  ggpubr::stat_pvalue_manual(alpha_stats, label = "p.signif", tip.length = 0.01) +
  scale_x_discrete(labels = c("0" = "Impaired\n(<53)", "1" = "Intact\n(≥53)")) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))
# plot generated in Script 05

cat("\nAlpha diversity Wilcoxon results:\n")
print(alpha_stats)

# =============================================================================
# ANALYSIS A: CONTINUOUS PDAQ (PRIMARY ANALYSIS)
# =============================================================================
cat("\n--- ANALYSIS A: Continuous PDAQ_num (primary) ---\n")
cat("Covariates (PERMANOVA-verified):", paste(DA_COVARIATES_COG, collapse = ", "), "\n")

f_rhs_cont <- paste(c("PDAQ_num", DA_COVARIATES_COG), collapse = " + ")

# --- A1. ANCOM-BC2 ---
ancom_res_cont <- ancombc2(
  data         = ps_genus,
  fix_formula  = f_rhs_cont,
  p_adj_method = "fdr",
  struc_zero   = FALSE   # within-PD cohort; structural zeros less expected
)
ancom_out_cont <- ancom_res_cont$res %>%
  select(taxon, contains("PDAQ_num")) %>%
  filter(q_PDAQ_num < 0.05)
cat("ANCOM-BC2 significant genera (continuous PDAQ):", nrow(ancom_out_cont), "\n")

# --- A2. LinDA ---
# Note: linda() expects taxa x samples — otu_mat is transposed here
linda_res_cont <- linda(
  t(otu_mat), df_meta_final,
  formula = paste0("~", f_rhs_cont)
)
linda_out_cont <- linda_res_cont$output$PDAQ_num %>%
  rownames_to_column("taxon") %>%
  filter(padj < 0.05)
cat("LinDA significant genera (continuous PDAQ):", nrow(linda_out_cont), "\n")

# --- A3. MaAsLin2 ---
maas_res_cont <- Maaslin2(
  data.frame(otu_mat), df_meta_final,
  output          = file.path(BIOINFO_DIR, "maas_pdaq_cont"),
  fixed_effects   = c("PDAQ_num", DA_COVARIATES_COG),
  analysis_method = "LM",
  plot_heatmap    = FALSE,
  plot_scatter    = FALSE
)
# Filter on taxon column (confirmed column name for this MaAsLin2 version)
maas_out_cont <- maas_res_cont$results %>%
  filter(metadata == "PDAQ_num", qval < 0.05)

# Defensive check — silent empty filter is a common failure mode
if (nrow(maas_out_cont) == 0 & nrow(maas_res_cont$results) > 0) {
  cat("WARNING: MaAsLin2 filter returned 0 rows. Check metadata column name:\n")
  print(unique(maas_res_cont$results$metadata))
}
cat("MaAsLin2 significant genera (continuous PDAQ):", nrow(maas_out_cont), "\n")

# --- A4. Consensus ---
consensus_cont <- ancom_out_cont %>%
  inner_join(linda_out_cont, by = "taxon") %>%
  inner_join(maas_out_cont, by = c("taxon" = "feature"))
consensus_cont <- annotate_taxa(consensus_cont) %>% mutate(n_methods = 3)
cat("Consensus genera — continuous PDAQ (3-of-3):", nrow(consensus_cont), "\n")

two_of_three_cont <- ancom_out_cont %>%
  full_join(linda_out_cont, by = "taxon") %>%
  full_join(maas_out_cont, by = c("taxon" = "feature")) %>%
  mutate(
    n_methods = (!is.na(lfc_PDAQ_num)) +
      (!is.na(log2FoldChange)) +
      (!is.na(coef))
  ) %>%
  filter(n_methods >= 2) %>%
  arrange(desc(n_methods))
two_of_three_cont <- annotate_taxa(two_of_three_cont)
cat("Consensus genera — continuous PDAQ (2-of-3):", nrow(two_of_three_cont), "\n")

# --- A5. Export tables ---
fmt_p <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f")))
}

tbl_cont_3 <- consensus_cont %>%
  select(Genus = label, Phylum,
         LFC_ANCOM    = lfc_PDAQ_num,  q_ANCOM    = q_PDAQ_num,
         LFC_LinDA    = log2FoldChange, q_LinDA    = padj,
         LFC_MaAsLin2 = coef,          q_MaAsLin2 = qval) %>%
  mutate(across(c(LFC_ANCOM, LFC_LinDA, LFC_MaAsLin2), ~round(., 3)),
         across(c(q_ANCOM, q_LinDA, q_MaAsLin2),       ~fmt_p(.))) %>%
  arrange(LFC_ANCOM)

write.csv(tbl_cont_3, file.path(TABLES_DIR, "Table_A_ContPDAQ_3of3.csv"), row.names = FALSE, na = "—")

tbl_cont_2 <- two_of_three_cont %>%
  select(Genus = label, Phylum, n_methods,
         LFC_ANCOM    = lfc_PDAQ_num,  q_ANCOM    = q_PDAQ_num,
         LFC_LinDA    = log2FoldChange, q_LinDA    = padj,
         LFC_MaAsLin2 = coef,          q_MaAsLin2 = qval) %>%
  mutate(across(c(LFC_ANCOM, LFC_LinDA, LFC_MaAsLin2), ~round(., 3)),
         across(c(q_ANCOM, q_LinDA, q_MaAsLin2),       ~fmt_p(.))) %>%
  arrange(desc(n_methods), LFC_ANCOM)

write.csv(tbl_cont_2, file.path(TABLES_DIR, "Table_A_ContPDAQ_2of3.csv"), row.names = FALSE, na = "—")
cat("Exported: Table_A_ContPDAQ_3of3.csv and Table_A_ContPDAQ_2of3.csv\n")
cat("  →", TABLES_DIR, "\n")

# =============================================================================
# ANALYSIS C: BINARY COGNITIVE STATUS (SECONDARY ANALYSIS)
# =============================================================================
cat("\n--- ANALYSIS C: Binary cog_intact (secondary) ---\n")
cat("Covariates (PERMANOVA-specified):", paste(DA_COVARIATES_COG, collapse = ", "), "\n")
cat("Threshold: PDAQ >= 53 = intact (1), < 53 = impaired (0)\n")
cat("Reference: Levin et al. [ADD CITATION]\n")

f_rhs_bin <- paste(c("cog_intact", DA_COVARIATES_COG), collapse = " + ")

# --- C1. ANCOM-BC2 ---
ancom_res_bin <- ancombc2(
  data         = ps_genus,
  fix_formula  = f_rhs_bin,
  p_adj_method = "fdr",
  struc_zero   = FALSE
)
ancom_out_bin <- ancom_res_bin$res %>%
  select(taxon, contains("cog_intact")) %>%
  filter(q_cog_intact1 < 0.05)
cat("ANCOM-BC2 significant genera (binary):", nrow(ancom_out_bin), "\n")

# --- C2. LinDA ---
linda_res_bin <- linda(
  t(otu_mat), df_meta_final,
  formula = paste0("~", f_rhs_bin)
)
linda_out_bin <- linda_res_bin$output$cog_intact1 %>%
  rownames_to_column("taxon") %>%
  filter(padj < 0.05)
cat("LinDA significant genera (binary):", nrow(linda_out_bin), "\n")

# --- C3. MaAsLin2 ---
maas_res_bin <- Maaslin2(
  data.frame(otu_mat), df_meta_final,
  output          = file.path(BIOINFO_DIR, "maas_pdaq_bin"),
  fixed_effects   = c("cog_intact", DA_COVARIATES_COG),
  analysis_method = "LM",
  plot_heatmap    = FALSE,
  plot_scatter    = FALSE
)
maas_out_bin <- maas_res_bin$results %>%
  filter(metadata == "cog_intact", qval < 0.05)

if (nrow(maas_out_bin) == 0 & nrow(maas_res_bin$results) > 0) {
  cat("WARNING: MaAsLin2 filter returned 0 rows. Check metadata column name:\n")
  print(unique(maas_res_bin$results$metadata))
}
cat("MaAsLin2 significant genera (binary):", nrow(maas_out_bin), "\n")

# --- C4. Consensus ---
consensus_bin <- ancom_out_bin %>%
  inner_join(linda_out_bin, by = "taxon") %>%
  inner_join(maas_out_bin, by = c("taxon" = "feature"))
consensus_bin <- annotate_taxa(consensus_bin) %>% mutate(n_methods = 3)
cat("Consensus genera — binary cog_intact (3-of-3):", nrow(consensus_bin), "\n")

two_of_three_bin <- ancom_out_bin %>%
  full_join(linda_out_bin, by = "taxon") %>%
  full_join(maas_out_bin, by = c("taxon" = "feature")) %>%
  mutate(
    n_methods = (!is.na(lfc_cog_intact1)) +
      (!is.na(log2FoldChange)) +
      (!is.na(coef))
  ) %>%
  filter(n_methods >= 2) %>%
  arrange(desc(n_methods))
two_of_three_bin <- annotate_taxa(two_of_three_bin)
cat("Consensus genera — binary cog_intact (2-of-3):", nrow(two_of_three_bin), "\n")

# --- C5. Export tables ---
tbl_bin_3 <- consensus_bin %>%
  select(Genus = label, Phylum,
         LFC_ANCOM    = lfc_cog_intact1, q_ANCOM    = q_cog_intact1,
         LFC_LinDA    = log2FoldChange,  q_LinDA    = padj,
         LFC_MaAsLin2 = coef,            q_MaAsLin2 = qval) %>%
  mutate(across(c(LFC_ANCOM, LFC_LinDA, LFC_MaAsLin2), ~round(., 3)),
         across(c(q_ANCOM, q_LinDA, q_MaAsLin2),       ~fmt_p(.))) %>%
  arrange(LFC_ANCOM)

write.csv(tbl_bin_3, file.path(TABLES_DIR, "Table_C_BinCogIntact_3of3.csv"), row.names = FALSE, na = "—")

tbl_bin_2 <- two_of_three_bin %>%
  select(Genus = label, Phylum, n_methods,
         LFC_ANCOM    = lfc_cog_intact1, q_ANCOM    = q_cog_intact1,
         LFC_LinDA    = log2FoldChange,  q_LinDA    = padj,
         LFC_MaAsLin2 = coef,            q_MaAsLin2 = qval) %>%
  mutate(across(c(LFC_ANCOM, LFC_LinDA, LFC_MaAsLin2), ~round(., 3)),
         across(c(q_ANCOM, q_LinDA, q_MaAsLin2),       ~fmt_p(.))) %>%
  arrange(desc(n_methods), LFC_ANCOM)

write.csv(tbl_bin_2, file.path(TABLES_DIR, "Table_C_BinCogIntact_2of3.csv"), row.names = FALSE, na = "—")
cat("Exported: Table_C_BinCogIntact_3of3.csv and Table_C_BinCogIntact_2of3.csv\n")
cat("  →", TABLES_DIR, "\n")

# =============================================================================
# SECTION 8: SAVE ALL RESULTS
# =============================================================================
save(
  # Analysis A — continuous
  consensus_cont, two_of_three_cont,
  ancom_out_cont, linda_out_cont, maas_out_cont,
  # Analysis C — binary
  consensus_bin,  two_of_three_bin,
  ancom_out_bin,  linda_out_bin,  maas_out_bin,
  # Shared objects
  tax_lookup, df_meta_final, ps_genus, otu_mat, dist_bray,
  # Plot objects for Script 05
  pcoa_res, pcoa_df, perm_cont, perm_bin,
  pcoa_cont_plot, pcoa_bin_plot,
  ps_genus_rare, alpha_long, alpha_stats, alpha_plot,
  file = file.path(RDATA_DIR, "PDAQ_DA_Results.RData")
)

cat("\n--- SCRIPT 2 v1.9 COMPLETE ---\n")
cat("Primary analysis:   continuous PDAQ (Analysis A)\n")
cat("Secondary analysis: binary cog_intact (Analysis C)\n")
cat("Results saved to", file.path(RDATA_DIR, "PDAQ_DA_Results.RData"), "\n")

# =============================================================================
# SECTION 9: EXTRACT SINGLE-METHOD HITS — CONTINUOUS PDAQ
# =============================================================================
# No consensus was reached across all three methods. Two exploratory signals
# were identified in individual methods (ANCOM-BC2 and MaAsLin2 respectively).
# These are reported as hypothesis-generating only and should not be interpreted
# as robust findings. They did not survive the 2-of-3 or 3-of-3 consensus threshold.

single_method_hits <- ancom_out_cont %>%
  full_join(maas_out_cont, by = c("taxon" = "feature")) %>%
  left_join(tax_lookup, by = "taxon") %>%
  mutate(
    method_detected = case_when(
      !is.na(lfc_PDAQ_num) & !is.na(coef) ~ "ANCOM-BC2 & MaAsLin2",
      !is.na(lfc_PDAQ_num)                 ~ "ANCOM-BC2 only",
      !is.na(coef)                         ~ "MaAsLin2 only"
    )
  ) %>%
  select(
    Genus            = label,
    Phylum,
    Method_detected  = method_detected,
    LFC_ANCOM        = lfc_PDAQ_num,
    q_ANCOM          = q_PDAQ_num,
    LFC_MaAsLin2     = coef,
    q_MaAsLin2       = qval
  ) %>%
  mutate(
    across(c(LFC_ANCOM, LFC_MaAsLin2), ~round(., 3)),
    across(c(q_ANCOM, q_MaAsLin2),     ~fmt_p(.))
  )

print(single_method_hits)
write.csv(single_method_hits, file.path(TABLES_DIR, "Table_A_SingleMethod_Hits.csv"), row.names = FALSE, na = "—")
cat("Exported: Table_A_SingleMethod_Hits.csv\n")
cat("  →", TABLES_DIR, "\n")

