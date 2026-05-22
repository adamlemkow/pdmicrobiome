# =============================================================================
# Script 3 v2.3: Longitudinal PDAQ Prediction from Baseline Microbiome
# Two-stage slope approach
# =============================================================================
# Changes from v2.2:
#   - Sellimonas ggsave() removed — plot generated in Script 05
#   - p_sellimonas, plot_df, plot_df_sens, sellimonas_taxon added to save()
#   - ids_4plus, slopes_4plus, stage2_sens added to save()
#   - Sensitivity scatter data (plot_df_sens) built at end of script using
#     maas_slope_sensitivity norm file; NULL if file not yet generated
#   - RData output updated to v2.3
#   - Output paths updated to new directory structure:
#       RData/ps objects → ~/rdata/
#       Tables           → ~/tables/
#       Bioinformatics   → ~/bioinfo_results/
#       Plots            → ~/plots/
#   - Input path updated: loads PD_PDAQ_Association_Data.RData from ~/rdata/
#   - MaAsLin2 output directories updated (maas_slope, maas_slope_sensitivity)
#   - Sellimonas norm_file path updated to resolve from ~/bioinfo_results/
#   - RData output filename updated to v2.2
#   - N is present in Sellimonas scatter caption (only figure); directive met
# Changes from v2.0 (v2.1):
#   - Single-visit participants explicitly excluded from Stage 2
#   - Section 7: Sellimonas scatter plot added
#   - Section 9: Sensitivity analysis (4+ visits) added
#   - LinDA outlier.pct argument corrected
#   - Output RData updated to v2.1
#
# Research question addressed:
#   Q3 — Does baseline gut microbiome composition predict the rate of
#         cognitive decline (PDAQ slope over time) in PD patients?
#
# Inputs:  ~/rdata/PD_PDAQ_Association_Data.RData (ps_pdaq, from Script 1)
#          foxden_pdaq15_raw.csv (all PDAQ timepoints)
#
# Outputs: ~/rdata/03_Longitudinal_PDAQ_Results_v2.2.RData
#          ~/tables/Table_03_Slopes_Summary.csv     (slope distribution)
#          ~/tables/Table_03_Consensus_2of2.csv     (MaAsLin2 + LinDA hits)
#          ~/tables/Table_03_SingleMethod.csv       (marginal single-method hits)
#          ~/plots/Figure_03_Sellimonas_slope.pdf   (scatter: abundance vs slope)
#
# Stage 1 model:
#   PDAQ_clean ~ time_years + (1 + time_years | clean_id)
#   Fitted via lme4::lmer on all valid PDAQ timepoints per participant.
#   Participants with <2 timepoints excluded from Stage 2.
#   time_years = time_months / 12; slope unit = PDAQ points per year.
#
# Stage 2 model:
#   taxon_abundance ~ pdaq_slope + current_age + sex + BMI_clean + batch
#   MaAsLin2: LM method, TSS + log normalisation, BH FDR correction.
#   LinDA: adaptive winsorisation (outlier.pct = 0.03), BH FDR correction.
#   Consensus: taxa significant in both methods (q < 0.05).
#
# Method note:
#   ANCOM-BC2 excluded — no support for continuous numeric outcomes in
#   the version available in this library path.
# =============================================================================

# --- SECTION 1: INITIALIZATION ---
.libPaths("/arc/project/st-silkec-1/rsandboxlib")
setwd("/arc/project/st-silkec-1/pacman")

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
library(lme4)
library(Maaslin2)
library(MicrobiomeStat)   # LinDA
library(ggplot2)
library(ggpubr)

# =============================================================================
# SECTION 2: LOAD ALL PDAQ TIMEPOINTS
# =============================================================================
# Load the full longitudinal PDAQ file — all timepoints per participant.
# Script 01 retains only the baseline row; this script uses all rows.

pdaq_raw <- read.csv(
  "/arc/project/st-silkec-1/pacman/other_foxden/foxden_pdaq15_raw.csv",
  stringsAsFactors = FALSE
)

cat("--- RAW PDAQ FILE ---\n")
cat("Total rows:", nrow(pdaq_raw), "\n")
cat("Unique participants:", n_distinct(pdaq_raw$fox_insight_id), "\n")

daily_cols <- c("DailyRead", "DailyClock", "DailyMoney", "DailyInstruct",
                "DailyProblem", "DailyExplain", "DailyErrand", "DailyMap",
                "DailyNumber", "DailyMany", "DailyLearn", "DailyFinance",
                "DailyThought", "DailyDiscuss", "DailyRememb")

# Recode value 5 (prefer not to answer) to NA before summing.
# Prevents impossible scores above 60 (max valid = 15 items x 4 = 60).
pdaq_items <- pdaq_raw[, daily_cols]
pdaq_items[pdaq_items == 5] <- NA
pdaq_raw$PDAQ_clean <- rowSums(pdaq_items, na.rm = FALSE)

pdaq_raw <- pdaq_raw %>%
  mutate(
    clean_id = as.character(as.numeric(gsub("[^0-9]", "", fox_insight_id))),
    time_months = case_when(
      trimws(as.character(schedule_of_activities)) == "REG" ~ 0,
      !is.na(suppressWarnings(as.numeric(as.character(schedule_of_activities)))) ~
        suppressWarnings(as.numeric(as.character(schedule_of_activities))),
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(PDAQ_clean), !is.na(time_months))

cat("\nAfter filtering (complete PDAQ_clean + valid timepoint):\n")
cat("Total rows:", nrow(pdaq_raw), "\n")
cat("Unique participants:", n_distinct(pdaq_raw$clean_id), "\n")

# =============================================================================
# SECTION 3: BASELINE MICROBIOME PREPARATION
# =============================================================================

meta_baseline <- data.frame(sample_data(ps_pdaq)) %>%
  filter(parkinsons == "PD case") %>%
  mutate(
    clean_id         = rownames(.),
    current_age      = as.numeric(as.character(current_age)),
    BMI_clean        = as.numeric(as.character(BMI)),
    sex              = as.factor(sex),
    batch            = as.factor(main_or_pilot_study),
  ) %>%
  select(clean_id, current_age, BMI_clean, sex, batch) %>%
  drop_na()

cat("\n--- BASELINE MICROBIOME COHORT ---\n")
cat("PD cases with complete baseline covariates:", nrow(meta_baseline), "\n")

ps_base  <- prune_samples(meta_baseline$clean_id, ps_pdaq)
ps_genus <- tax_glom(ps_base, "Genus")
otu_mat  <- as.matrix(otu_table(ps_genus))
rownames(otu_mat) <- sample_names(ps_genus)

cat("Genus-level taxa after tax_glom:", ntaxa(ps_genus), "\n")

# =============================================================================
# SECTION 4: TAXONOMY LOOKUP & HELPERS
# =============================================================================

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

annotate_taxa <- function(df) {
  candidates <- c("taxon", "feature", "featureID", "OTU", "taxa")
  taxon_col  <- candidates[candidates %in% colnames(df)][1]
  if (is.na(taxon_col)) stop("Could not detect taxon column. Columns: ",
                             paste(colnames(df), collapse = ", "))
  df %>% left_join(tax_lookup, by = setNames("taxon", taxon_col))
}

fmt_p <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f")))
}

# =============================================================================
# SECTION 5: STAGE 1 — ESTIMATE PER-PERSON PDAQ SLOPES
# =============================================================================
# Fit a linear mixed model to all PDAQ timepoints within the analytical cohort.
# Random intercept + random slope by participant allows each person to have
# their own baseline PDAQ level AND their own rate of change over time.
# The fixed slope (time_years coefficient) captures the population mean;
# the BLUP (best linear unbiased predictor) for each participant's slope
# is extracted as their individual cognitive trajectory estimate.
#
# Participants with only one timepoint cannot contribute meaningful slope
# information — their BLUP would shrink entirely to the population mean,
# which is not a real individual estimate. These are excluded from Stage 2.

pdaq_cohort <- pdaq_raw %>%
  filter(clean_id %in% meta_baseline$clean_id)

cat("\n--- STAGE 1: PER-PERSON PDAQ SLOPE ESTIMATION ---\n")
cat("Total PDAQ observations in cohort:", nrow(pdaq_cohort), "\n")
cat("Unique participants:", n_distinct(pdaq_cohort$clean_id), "\n")

visits_n <- pdaq_cohort %>%
  group_by(clean_id) %>%
  summarise(n_visits = n(), .groups = "drop")

cat("Participants with only 1 timepoint (excluded from Stage 2):",
    sum(visits_n$n_visits == 1), "\n")
cat("Participants with 2+ timepoints (slope-estimable):",
    sum(visits_n$n_visits >= 2), "\n")
cat("Participants with 4+ timepoints (well-characterised slopes):",
    sum(visits_n$n_visits >= 4), "\n")

# time_months scaled to years so slope = PDAQ points per year
pdaq_cohort <- pdaq_cohort %>%
  mutate(time_years = time_months / 12)

lmer_fit <- lmer(
  PDAQ_clean ~ time_years + (1 + time_years | clean_id),
  data    = pdaq_cohort,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 2e5))
)

cat("\n--- Population-level PDAQ trajectory (lmer) ---\n")
print(summary(lmer_fit)$coefficients)

# Extract per-person slopes: fixed slope + individual random effect
re          <- ranef(lmer_fit)$clean_id
fixed_slope <- fixef(lmer_fit)["time_years"]

slopes_df <- re %>%
  rownames_to_column("clean_id") %>%
  rename(re_intercept = `(Intercept)`, re_slope = time_years) %>%
  mutate(pdaq_slope = fixed_slope + re_slope) %>%
  select(clean_id, pdaq_slope)

# Exclude single-visit participants — their slope is not a real estimate
slopes_df <- slopes_df %>%
  filter(clean_id %in% visits_n$clean_id[visits_n$n_visits >= 2])

cat("\n--- Individual PDAQ slope distribution (points per year) ---\n")
cat("Participants with estimable slopes (2+ visits):", nrow(slopes_df), "\n")
print(summary(slopes_df$pdaq_slope))
cat("Participants with negative slope (declining):",
    sum(slopes_df$pdaq_slope < 0), "\n")
cat("Participants with positive slope (improving):",
    sum(slopes_df$pdaq_slope > 0), "\n")

# Export slope summary table
tbl_slopes <- slopes_df %>%
  left_join(visits_n, by = "clean_id") %>%
  arrange(pdaq_slope)

write.csv(tbl_slopes, file.path(TABLES_DIR, "Table_03_Slopes_Summary.csv"),
          row.names = FALSE, na = "—")
cat("Exported: Table_03_Slopes_Summary.csv\n")
cat("  →", TABLES_DIR, "\n")

# =============================================================================
# SECTION 6: STAGE 2 — DIFFERENTIAL ABUNDANCE VS PDAQ SLOPE
# =============================================================================
# One row per participant. Outcome = individual PDAQ slope (points/year).
# Predictors = baseline genus-level OTU abundances + pre-specified covariates.
# Single-visit participants are absent from slopes_df and excluded automatically.

stage2_meta <- slopes_df %>%
  inner_join(meta_baseline, by = "clean_id") %>%
  column_to_rownames("clean_id")

cat("\n--- STAGE 2: DIFFERENTIAL ABUNDANCE vs PDAQ SLOPE ---\n")
cat("Participants entering Stage 2:", nrow(stage2_meta), "\n")

# Align OTU matrix to Stage 2 participants
otu_stage2 <- otu_mat[rownames(stage2_meta), , drop = FALSE]

# Prevalence filter at 10% — appropriate for cross-sectional N ~290
prev_threshold <- 0.10
prevalence     <- colSums(otu_stage2 > 0) / nrow(otu_stage2)
otu_filt       <- otu_stage2[, prevalence >= prev_threshold, drop = FALSE]
cat("Prevalence filter (>=10%):", ncol(otu_filt), "taxa retained,",
    ncol(otu_stage2) - ncol(otu_filt), "removed\n")

DA_COVARIATES <- c("current_age", "sex", "BMI_clean", "batch")

# --- 6A. MaAsLin2 ---
cat("\n--- Running MaAsLin2 ---\n")
maas_res <- Maaslin2(
  input_data      = data.frame(otu_filt),
  input_metadata  = stage2_meta,
  output          = file.path(BIOINFO_DIR, "maas_slope"),
  fixed_effects   = c("pdaq_slope", DA_COVARIATES),
  analysis_method = "LM",
  normalization   = "TSS",
  transform       = "LOG",
  plot_heatmap    = FALSE,
  plot_scatter    = FALSE
)

maas_out <- maas_res$results %>%
  filter(metadata == "pdaq_slope", qval < 0.05) %>%
  rename(taxon = feature)

cat("MaAsLin2 significant genera (q < 0.05):", nrow(maas_out), "\n")

# --- 6B. LinDA ---
cat("\n--- Running LinDA ---\n")
linda_res <- linda(
  feature.dat      = data.frame(t(otu_filt)),
  meta.dat         = stage2_meta,
  formula          = paste("~", paste(c("pdaq_slope", DA_COVARIATES),
                                      collapse = " + ")),
  feature.dat.type = "count",
  p.adj.method     = "BH",
  alpha            = 0.05,
  outlier.pct      = 0.03
)

linda_out <- linda_res$output$pdaq_slope %>%
  rownames_to_column("taxon") %>%
  filter(padj < 0.05)

cat("LinDA significant genera (padj < 0.05):", nrow(linda_out), "\n")

# =============================================================================
# SECTION 7: CONSENSUS & RESULTS
# =============================================================================

# --- 2-of-2 consensus ---
consensus_2of2 <- inner_join(
  maas_out  %>% select(taxon, coef_maas = coef,           q_maas  = qval),
  linda_out %>% select(taxon, lfc_linda = log2FoldChange,  q_linda = padj),
  by = "taxon"
)

cat("\n--- CONSENSUS (both MaAsLin2 + LinDA, q < 0.05) ---\n")
cat("Consensus hits:", nrow(consensus_2of2), "\n")

# --- Single-method hits ---
all_taxa      <- union(maas_out$taxon, linda_out$taxon)
single_method <- all_taxa[!all_taxa %in% consensus_2of2$taxon]

single_out <- bind_rows(
  maas_out  %>% filter(taxon %in% single_method) %>%
    mutate(method = "MaAsLin2") %>%
    select(taxon, method, coef, qval),
  linda_out %>% filter(taxon %in% single_method) %>%
    mutate(method = "LinDA", coef = log2FoldChange, qval = padj) %>%
    select(taxon, method, coef, qval)
)

cat("Single-method marginal hits:", length(single_method), "\n")

# --- Annotate ---
consensus_annotated <- consensus_2of2 %>%
  annotate_taxa() %>%
  select(Genus = label, Phylum,
         LFC_MaAsLin2 = coef_maas, q_MaAsLin2 = q_maas,
         LFC_LinDA    = lfc_linda,  q_LinDA    = q_linda) %>%
  mutate(across(starts_with("LFC"), ~ round(.x, 3)),
         across(starts_with("q_"),  fmt_p)) %>%
  arrange(LFC_MaAsLin2)

single_annotated <- single_out %>%
  annotate_taxa() %>%
  select(Genus = label, Phylum, Method = method, coef, qval) %>%
  mutate(coef = round(coef, 3), qval = fmt_p(qval)) %>%
  arrange(Method, coef)

cat("\nConsensus genera (MaAsLin2 + LinDA):\n")
print(consensus_annotated)

cat("\nMarginal single-method hits:\n")
print(single_annotated)

# --- Export ---
write.csv(consensus_annotated,
          file.path(TABLES_DIR, "Table_03_Consensus_2of2.csv"), row.names = FALSE, na = "—")
cat("Exported: Table_03_Consensus_2of2.csv\n")

write.csv(single_annotated,
          file.path(TABLES_DIR, "Table_03_SingleMethod.csv"), row.names = FALSE, na = "—")
cat("Exported: Table_03_SingleMethod.csv\n")
cat("  →", TABLES_DIR, "\n") 


# =============================================================================
# SECTION 8: SELLIMONAS VISUALISATION
# =============================================================================
# Builds scatter plot of baseline Sellimonas abundance vs individual PDAQ slope.
# Plot object passed to Script 05 for final rendering.

# Initialise as NULL — overwritten below if Sellimonas found and norm file exists
sellimonas_taxon <- NULL
plot_df          <- NULL
p_sellimonas     <- NULL

sellimonas_taxon <- tax_lookup %>%
  filter(grepl("Sellimonas", label, ignore.case = TRUE)) %>%
  pull(taxon)

if (length(sellimonas_taxon) == 0) {
  cat("NOTE: Sellimonas not found in tax_lookup — scatter plot skipped\n")
} else {
  cat("Sellimonas taxon ID:", sellimonas_taxon, "\n")
  maas_norm_path <- file.path(BIOINFO_DIR,
                              "maas_slope/features/filtered_data_norm_transformed.tsv")
  cat("Norm file exists:", file.exists(maas_norm_path), "\n")
  
  if (file.exists(maas_norm_path)) {
    maas_norm <- read.delim(maas_norm_path, row.names = 1, check.names = FALSE)
    
    sellimonas_abund <- data.frame(
      clean_id       = rownames(maas_norm),
      sellimonas_log = maas_norm[, sellimonas_taxon]
    )
    plot_df <- sellimonas_abund %>%
      inner_join(slopes_df, by = "clean_id")
    
    p_sellimonas <- ggplot(plot_df, aes(x = sellimonas_log, y = pdaq_slope)) +
      geom_point(alpha = 0.5, size = 2, colour = "#D55E00") +
      geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.8) +
      labs(
        title    = "Baseline Sellimonas Abundance vs PDAQ Cognitive Slope",
        subtitle = paste0("N = ", nrow(plot_df), " PD cases (2+ visits)"),
        x        = "Sellimonas (log-normalised relative abundance)",
        y        = "Individual PDAQ slope (points/year)",
        caption  = paste0(
          "Linear fit with 95% CI. Negative slope = cognitive decline. ",
          "N = ", nrow(plot_df), ".")
      ) +
      theme_bw()
    
    print(p_sellimonas)
    cat("Sellimonas scatter plot built (N =", nrow(plot_df), ")\n")
  } else {
    cat("NOTE: MaAsLin2 norm file not found — Sellimonas scatter skipped\n")
  }
}

# =============================================================================
# SECTION 9: SENSITIVITY — RESTRICT TO PARTICIPANTS WITH 4+ VISITS
# =============================================================================
# Tests whether the consensus hit(s) are robust to excluding participants
# with sparse longitudinal data (2–3 visits), whose individual slopes are
# less reliably estimated. A consistent direction with nominal p < 0.05 is
# the threshold here — FDR correction will be conservative given the reduced N.

cat("\n--- SECTION 9: SENSITIVITY ANALYSIS (4+ visits only) ---\n")

ids_4plus    <- visits_n$clean_id[visits_n$n_visits >= 4]
stage2_sens  <- stage2_meta[rownames(stage2_meta) %in% ids_4plus, ]
otu_sens     <- otu_filt[rownames(otu_filt) %in% ids_4plus, ]

cat("Participants in sensitivity analysis (4+ visits):", nrow(stage2_sens), "\n")

maas_sens <- Maaslin2(
  input_data      = data.frame(otu_sens),
  input_metadata  = stage2_sens,
  output          = file.path(BIOINFO_DIR, "maas_slope_sensitivity"),
  fixed_effects   = c("pdaq_slope", DA_COVARIATES),
  analysis_method = "LM",
  normalization   = "TSS",
  transform       = "LOG",
  plot_heatmap    = FALSE,
  plot_scatter    = FALSE
)

# Check each consensus hit in the sensitivity model
cat("\nSensitivity results for consensus genera:\n")
sens_consensus <- maas_sens$results %>%
  filter(metadata == "pdaq_slope",
         feature  %in% consensus_2of2$taxon) %>%
  left_join(tax_lookup, by = c("feature" = "taxon")) %>%
  select(Genus = label, coef, pval, qval) %>%
  mutate(coef = round(coef, 3),
         pval = fmt_p(pval),
         qval = fmt_p(qval),
         direction_consistent = ifelse(
           sign(coef) == sign(consensus_annotated$LFC_MaAsLin2[
             match(Genus, consensus_annotated$Genus)]),
           "yes", "NO — check"))

print(sens_consensus)
cat("(Direction consistent = same sign as primary analysis)\n")

# =============================================================================
# SECTION 10: SAVE RESULTS
# =============================================================================
# Build sensitivity scatter data for Script 05 (4+ visit Sellimonas plot)
if (exists("sellimonas_taxon") && length(sellimonas_taxon) > 0 &&
    file.exists(file.path(BIOINFO_DIR, "maas_slope_sensitivity/features/filtered_data_norm_transformed.tsv"))) {
  maas_norm_sens <- read.delim(
    file.path(BIOINFO_DIR, "maas_slope_sensitivity/features/filtered_data_norm_transformed.tsv"),
    row.names = 1, check.names = FALSE)
  slopes_4plus <- slopes_df %>%
    filter(clean_id %in% ids_4plus)
  sellimonas_abund_sens <- data.frame(
    clean_id       = rownames(maas_norm_sens),
    sellimonas_log = maas_norm_sens[, sellimonas_taxon]
  )
  plot_df_sens <- sellimonas_abund_sens %>%
    inner_join(slopes_4plus, by = "clean_id")
} else {
  plot_df_sens <- NULL
  slopes_4plus <- NULL
  cat("NOTE: Sensitivity scatter data not built — norm file missing or Sellimonas not found\n")
}

save(
  lmer_fit, slopes_df, visits_n,
  maas_res, maas_out,
  linda_res, linda_out,
  consensus_2of2, consensus_annotated,
  single_out, single_annotated,
  maas_sens, sens_consensus,
  tax_lookup, meta_baseline, otu_mat, otu_filt,
  stage2_meta, stage2_sens, pdaq_cohort,
  ids_4plus, slopes_4plus,
  sellimonas_taxon,
  plot_df, plot_df_sens,
  p_sellimonas,
  file = file.path(RDATA_DIR, "03_Longitudinal_PDAQ_Results_v2.3.RData")
)

cat("\n--- SCRIPT 3 v2.3 COMPLETE ---\n")
cat("Stage 1: per-person PDAQ slopes via lmer (2+ visit participants only)\n")
cat("Stage 2: MaAsLin2 + LinDA vs pdaq_slope (N =", nrow(stage2_meta), ")\n")
cat("Sensitivity: 4+ visit participants (N =", nrow(stage2_sens), ")\n")
cat("Results saved to", file.path(RDATA_DIR, "03_Longitudinal_PDAQ_Results_v2.3.RData"), "\n")

