# =============================================================================
# Script 1.5 v1.5: PD vs Control — Genus-Level Differential Abundance
# =============================================================================
# Changes from v1.4 (v1.5):
#   - all_covariates renamed to DA_COVARIATES_POP throughout
#   - sig_covariates renamed to DA_COVARIATES_POP_SIG throughout
#   - Naming convention now mirrors Script 02 (_POP vs _COG suffix)
# Changes from v1.3:
#   - 10% prevalence filter added to otu_mat and ps_genus before DA modelling
#     (consistent with Scripts 02 and 03); ps_genus_full retained unfiltered
#     for PERMANOVA so Bray-Curtis reflects total community composition
#   - annotate_taxa() candidates updated to include "feature" (MaAsLin2 col name)
#   - Completion message version number corrected to v1.4
# Changes from v1.2:
#   - Volcano plot print() and ggsave() removed — plot generated in Script 05
#   - volcano_df and p_volcano added to save() call
#   - Script 05 loads PD_vs_Control_DA_Results.RData to render volcano
#   - Output paths updated to new directory structure:
#       RData/ps objects → ~/rdata/
#       Tables           → ~/tables/
#       Bioinformatics   → ~/bioinfo_results/
#   - Input path updated: loads PD_vs_Control_Cleaned.RData from ~/rdata/
#   - MaAsLin2 output directory updated to ~/bioinfo_results/maas_pd_vs_control
#   - No figures in this script; N directive not applicable here
# Changes from v1.0 (v1.1):
#   - Fixed input section reference: PD_vs_Control_Cleaned.RData is saved in
#     Script 01 Section 8.2, not Section 9.2 as previously documented
#   - Added schedule_interpolated awareness note: this script loads ps_clean
#     which carries no PDAQ or schedule metadata, so baseline selection from
#     Script 01 v1.3 does not affect this script directly. Noted for clarity.
#
# Research question addressed:
#   Q2 — Which genera are differentially abundant between PD patients and
#         healthy controls?
#
# Inputs:  ~/rdata/PD_vs_Control_Cleaned.RData (saved by Script 1 Section 8.2)
# Outputs: ~/rdata/PD_vs_Control_DA_Results.RData
#          ~/tables/Table_PDvsControl_3of3.csv
#          ~/tables/Table_PDvsControl_2of3.csv
#          ~/bioinfo_results/maas_pd_vs_control/
#
# Analysis approach:
#   - Covariates selected via PERMANOVA (data-driven, avoids assumption)
#   - Three-method consensus DA: ANCOM-BC2, LinDA, MaAsLin2
#   - Results reported at 3-of-3 and 2-of-3 consensus thresholds
#   - Note: levodopa excluded as covariate — controls are not on PD medication
#   - Note: Bristol_Avg retained here (PD vs Control); mediation argument is
#     less clear-cut than within-PD analyses. Flag for biostatistician review.
# =============================================================================

# --- 1. INITIALIZATION ---
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

load(file.path(RDATA_DIR, "PD_vs_Control_Cleaned.RData"))

library(phyloseq)
library(tidyverse)
library(vegan)
library(ANCOMBC)
library(Maaslin2)
library(MicrobiomeStat)

# --- 2. METADATA PREPARATION ---
# Both PD and Control participants included.
# No levodopa covariate — not applicable to controls.
freq_map    <- c("never"=0, "< 3x/week"=1.5, "3-6x/week"=4.5,
                 "daily"=7, "2-3x/day"=17.5, "> 3x/day"=25)
target_cols <- c("bl_bristol_stool_type_1_freq", "bl_bristol_stool_type_2_freq",
                 "bl_bristol_stool_type_3_freq", "bl_bristol_stool_type_4_freq",
                 "bl_bristol_stool_type_5_freq", "bl_bristol_stool_type_6_freq",
                 "bl_bristol_stool_type_7_freq")

d_prep <- data.frame(sample_data(ps_clean)) %>%
  mutate(SampleID_Fixed = sample_names(ps_clean)) %>%
  mutate(
    current_age = as.numeric(as.character(current_age)),
    BMI_clean   = as.numeric(as.character(BMI)),
    sex         = as.factor(sex),
    batch       = as.factor(main_or_pilot_study),
    parkinsons  = as.factor(parkinsons)
  )

# Bristol calculation
# Verify columns exist — Stagman uses bl_bristol_stool_type_X_freq
missing_bristol <- setdiff(target_cols, colnames(d_prep))
if (length(missing_bristol) > 0) {
  stop("Bristol frequency columns not found in metadata:\n",
       paste(missing_bristol, collapse = "\n"),
       "\nCheck column names with: colnames(data.frame(sample_data(ps_clean)))")
}

unexpected <- d_prep %>%
  select(all_of(target_cols)) %>%
  summarise(across(everything(), ~sum(!as.character(.) %in% c(names(freq_map), NA))))
if (any(unexpected > 0)) {
  warning("Unexpected values in Bristol frequency columns — these will become NA:\n")
  print(unexpected[, unexpected > 0, drop = FALSE])
}

# Diagnostic: warn if >50% NA in any Bristol column
bristol_na_check <- d_prep %>%
  select(all_of(target_cols)) %>%
  summarise(across(everything(), ~mean(is.na(.))))
if (any(bristol_na_check > 0.5)) {
  warning(">50% NA in at least one Bristol column — check Stagman data quality:\n")
  print(bristol_na_check)
}

bristol_mat <- d_prep %>%
  mutate(across(all_of(target_cols), ~freq_map[as.character(.)], .names = "n_{.col}")) %>%
  select(starts_with("n_bl_bristol")) %>%
  as.matrix()

w_sum_vec   <- as.numeric(bristol_mat %*% 1:7)
total_f_vec <- rowSums(bristol_mat, na.rm = TRUE)

# --- 3. FILTERING & ALIGNMENT ---
df_meta_final <- d_prep %>%
  mutate(Bristol_Avg = ifelse(total_f_vec > 0, w_sum_vec / total_f_vec, NA)) %>%
  drop_na(current_age, sex, BMI_clean, Bristol_Avg, batch, parkinsons)

cat("\n--- FULL COHORT AFTER CLEANING ---\n")
cat("Total samples:", nrow(df_meta_final), "\n")
cat("\nBreakdown by PD status:\n")
print(table(df_meta_final$parkinsons))
cat("\nBreakdown by sex:\n")
print(table(df_meta_final$sex))

# Align phyloseq object
rownames(df_meta_final) <- df_meta_final$SampleID_Fixed
ps_final  <- prune_samples(df_meta_final$SampleID_Fixed, ps_clean)
ps_genus  <- tax_glom(ps_final, "Genus")
df_meta_final <- df_meta_final[sample_names(ps_genus), ]
sample_data(ps_genus) <- sample_data(df_meta_final)

cat("\nGenus-level taxa after tax_glom:", ntaxa(ps_genus), "\n")

# OTU matrix — samples x taxa (pre-filter, for PERMANOVA only)
otu_mat_full           <- as.matrix(otu_table(ps_genus))
rownames(otu_mat_full) <- sample_names(ps_genus)

# =============================================================================
# PREVALENCE FILTER — 10% threshold (consistent with Scripts 02 and 03)
# Applied to both ps_genus and otu_mat so all three DA methods operate on
# the same taxa. ps_genus_full retained unfiltered for PERMANOVA so
# Bray-Curtis reflects total community composition.
# =============================================================================
prev_threshold <- 0.10
prevalence     <- colSums(otu_mat_full > 0) / nrow(otu_mat_full)
keep_taxa      <- names(prevalence[prevalence >= prev_threshold])

ps_genus_full  <- ps_genus
ps_genus       <- prune_taxa(keep_taxa, ps_genus)
otu_mat        <- as.matrix(otu_table(ps_genus))
rownames(otu_mat) <- sample_names(ps_genus)

cat(sprintf("Prevalence filter (>=10%%): %d taxa retained, %d removed\n",
            length(keep_taxa),
            ntaxa(ps_genus_full) - length(keep_taxa)))

# Alignment checks
stopifnot(
  "OTU matrix rows do not match metadata rows" =
    all(rownames(otu_mat) == rownames(df_meta_final))
)

set.seed(123)
dist_bray <- vegdist(otu_mat_full, method = "bray")   # unfiltered for PERMANOVA

stopifnot(
  "Distance matrix labels do not match metadata rows" =
    all(labels(dist_bray) == rownames(df_meta_final))
)
cat("Alignment confirmed.\n")

# =============================================================================
# SECTION 4: PERMANOVA — COVARIATE SELECTION
# =============================================================================
# Run full PERMANOVA with all candidate covariates plus parkinsons as predictor.
# PERMANOVA-significant covariates (p < 0.05) are carried into DA models.
# This ensures covariate selection is data-driven and not assumed.
# Note: Bristol is included despite potential mediator status — flag for
# biostatistician review.

cat("\n--- PERMANOVA: Full model with all candidate covariates ---\n")
perm_full <- adonis2(dist_bray ~ parkinsons + current_age + sex + BMI_clean +
                       Bristol_Avg + batch,
                     data = df_meta_final, by = "margin")
print(perm_full)

# Extract significant covariates (p < 0.05), always keeping parkinsons
DA_COVARIATES_POP <- c("current_age", "sex", "BMI_clean", "Bristol_Avg", "batch")
DA_COVARIATES_POP_SIG <- DA_COVARIATES_POP[
  !is.na(perm_full[DA_COVARIATES_POP, "Pr(>F)"]) &
    perm_full[DA_COVARIATES_POP, "Pr(>F)"] < 0.05
]

f_rhs <- paste(c("parkinsons", DA_COVARIATES_POP_SIG), collapse = " + ")
cat("\nFormula for DA models (PERMANOVA-filtered covariates):", f_rhs, "\n")

# =============================================================================
# SECTION 5: CONSENSUS DIFFERENTIAL ABUNDANCE — PD vs Control
# =============================================================================

# 5A. ANCOM-BC2
# struc_zero = TRUE is valid here — we have a meaningful biological grouping
# (PD vs Control) where a taxon may be entirely absent in one group.
cat("\n--- Running ANCOM-BC2 ---\n")
ancom_res <- ancombc2(data = ps_genus, fix_formula = f_rhs,
                      p_adj_method = "fdr", struc_zero = TRUE,
                      group = "parkinsons")
ancom_out <- ancom_res$res %>%
  select(taxon, contains("parkinsonsPD case")) %>%
  filter(`q_parkinsonsPD case` < 0.05)
cat("ANCOM-BC2 significant genera:", nrow(ancom_out), "\n")

# 5B. LinDA — expects taxa x samples
cat("\n--- Running LinDA ---\n")
linda_res <- linda(t(otu_mat), df_meta_final,
                   formula = paste0("~", f_rhs))
linda_out <- linda_res$output$`parkinsonsPD case` %>%
  rownames_to_column("taxon") %>%
  filter(padj < 0.05)
cat("LinDA significant genera:", nrow(linda_out), "\n")

# 5C. MaAsLin2 — expects samples x taxa
# analysis_method = "LM" appropriate for pre-normalised/rarefied data.
# Switch to "NEGBIN" or "ZINB" if working with raw counts.
cat("\n--- Running MaAsLin2 ---\n")
maas_res <- Maaslin2(
  data.frame(otu_mat), df_meta_final,
  file.path(BIOINFO_DIR, "maas_pd_vs_control"),
  fixed_effects = c("parkinsons", DA_COVARIATES_POP_SIG),
  analysis_method = "LM",
  plot_heatmap    = FALSE,
  plot_scatter    = FALSE
)
maas_out_full <- maas_res$results
maas_out <- maas_res$results %>%
  filter(metadata == "parkinsons", qval < 0.05)

# Defensive check — silent empty filter is a common failure mode
if (nrow(maas_out) == 0 & nrow(maas_res$results) > 0) {
  cat("WARNING: MaAsLin2 filter returned 0 rows. Check metadata column name:\n")
  print(unique(maas_res$results$metadata))
}
cat("MaAsLin2 significant genera:", nrow(maas_out), "\n")

# =============================================================================
# SECTION 6: TAXONOMY LOOKUP
# =============================================================================
# Built once here and reused for annotation of all consensus result tables.

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
  candidates <- c("taxon", "feature", "featureID", "OTU", "taxa")
  taxon_col  <- candidates[candidates %in% colnames(df)][1]
  if (is.na(taxon_col)) stop("Could not detect taxon column. Columns: ",
                             paste(colnames(df), collapse = ", "))
  df %>% left_join(tax_lookup, by = setNames("taxon", taxon_col))
}

# Helper: format p/q values for publication tables
fmt_p <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f")))
}

# =============================================================================
# SECTION 7: CONSENSUS RESULTS
# =============================================================================

# 7.1 Strict 3-of-3 consensus
consensus_3 <- ancom_out %>%
  inner_join(linda_out, by = "taxon") %>%
  inner_join(maas_out, by = c("taxon" = "feature"))
consensus_3 <- annotate_taxa(consensus_3) %>% mutate(n_methods = 3)
cat("\nConsensus genera — PD vs Control (all 3 methods):", nrow(consensus_3), "\n")
print(consensus_3)

# 7.2 Permissive 2-of-3 consensus
two_of_three <- ancom_out %>%
  full_join(linda_out, by = "taxon") %>%
  full_join(maas_out, by = c("taxon" = "feature")) %>%
  mutate(
    n_methods = (!is.na(`lfc_parkinsonsPD case`)) +
      (!is.na(log2FoldChange)) +
      (!is.na(coef))
  ) %>%
  filter(n_methods >= 2) %>%
  arrange(desc(n_methods))
two_of_three <- annotate_taxa(two_of_three)
cat("Consensus genera — PD vs Control (2+ methods):", nrow(two_of_three), "\n")
print(two_of_three)

# =============================================================================
# SECTION 8: EXPORT TABLES
# =============================================================================

tbl_3of3 <- consensus_3 %>%
  select(Genus = label, Phylum,
         LFC_ANCOM    = `lfc_parkinsonsPD case`,
         q_ANCOM      = `q_parkinsonsPD case`,
         LFC_LinDA    = log2FoldChange,
         q_LinDA      = padj,
         LFC_MaAsLin2 = coef,
         q_MaAsLin2   = qval) %>%
  mutate(across(c(LFC_ANCOM, LFC_LinDA, LFC_MaAsLin2), ~round(., 3)),
         across(c(q_ANCOM,   q_LinDA,   q_MaAsLin2),   ~fmt_p(.))) %>%
  arrange(LFC_ANCOM)

write.csv(tbl_3of3, file.path(TABLES_DIR, "Table_PDvsControl_3of3.csv"), row.names = FALSE, na = "—")

tbl_2of3 <- two_of_three %>%
  select(Genus = label, Phylum, n_methods,
         LFC_ANCOM    = `lfc_parkinsonsPD case`,
         q_ANCOM      = `q_parkinsonsPD case`,
         LFC_LinDA    = log2FoldChange,
         q_LinDA      = padj,
         LFC_MaAsLin2 = coef,
         q_MaAsLin2   = qval) %>%
  mutate(across(c(LFC_ANCOM, LFC_LinDA, LFC_MaAsLin2), ~round(., 3)),
         across(c(q_ANCOM,   q_LinDA,   q_MaAsLin2),   ~fmt_p(.))) %>%
  arrange(desc(n_methods), LFC_ANCOM)

write.csv(tbl_2of3, file.path(TABLES_DIR, "Table_PDvsControl_2of3.csv"), row.names = FALSE, na = "—")
cat("Exported: Table_PDvsControl_3of3.csv and Table_PDvsControl_2of3.csv\n")
cat("  →", file.path(TABLES_DIR), "\n")

# =============================================================================
# SECTION 8.5: VOLCANO PLOT — PD vs Control DA Overview
# =============================================================================

library(ggrepel)

# Full MaAsLin2 results (add this line right after your Maaslin2() call if not
# already present, before the qval filter):
# maas_out_full <- maas_res$results

volcano_df <- maas_out_full %>%
  filter(metadata == "parkinsons") %>%
  rename(taxon = feature) %>%
  left_join(tax_lookup, by = "taxon") %>%
  mutate(
    neg_log10_q = -log10(qval + 1e-10),
    consensus = case_when(
      taxon %in% consensus_3$taxon    ~ "3-of-3 consensus",
      taxon %in% two_of_three$taxon   ~ "2-of-3 consensus",
      qval < 0.05                     ~ "Significant (single method)",
      TRUE                            ~ "Not significant"
    ),
    consensus = factor(consensus, levels = c(
      "3-of-3 consensus",
      "2-of-3 consensus",
      "Significant (single method)",
      "Not significant"
    )),
    # Label all consensus hits; pad to top 15 with next-best if needed
    is_consensus = taxon %in% c(consensus_3$taxon, two_of_three$taxon),
    lab = case_when(
      is_consensus                         ~ label,
      rank(-neg_log10_q) <= 15 & !is_consensus ~ label,
      TRUE                                 ~ NA_character_
    )
  )

pal_volcano <- c(
  "3-of-3 consensus"           = "#D55E00",
  "2-of-3 consensus"           = "#E69F00",
  "Significant (single method)"= "#0072B2",
  "Not significant"            = "grey75"
)

size_volcano <- c(
  "3-of-3 consensus"           = 4,
  "2-of-3 consensus"           = 3.5,
  "Significant (single method)"= 2.5,
  "Not significant"            = 1.2
)

p_volcano <- ggplot(volcano_df,
                    aes(x = coef, y = neg_log10_q,
                        colour = consensus, size = consensus)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text",
           x     = max(volcano_df$coef, na.rm = TRUE),
           y     = -log10(0.05) + 0.15,
           label = "FDR = 0.05", hjust = 1,
           size  = 3.2, colour = "grey40") +
  geom_point(alpha = 0.75) +
  geom_text_repel(
    aes(label = lab),
    size          = 3.2,
    fontface      = "italic",
    max.overlaps  = 20,
    box.padding   = 0.5,
    point.padding = 0.3,
    show.legend   = FALSE,
    segment.colour = "grey50",
    segment.size   = 0.3
  ) +
  scale_colour_manual(values = pal_volcano, name = NULL) +
  scale_size_manual(values = size_volcano, guide = "none") +
  scale_x_continuous(
    limits = c(min(volcano_df$coef, na.rm = TRUE) - 0.1,
               max(volcano_df$coef, na.rm = TRUE) + 0.1)
  ) +
  labs(
    title    = "Differential Abundance — PD vs Healthy Control",
    subtitle = paste0("N = ", nrow(df_meta_final),
                      "  |  ", ntaxa(ps_genus), " genera tested",
                      "  |  1 genus at 3-of-3 consensus, 10 at 2-of-3"),
    x        = "Log fold change (PD vs Control)",
    y        = expression(-log[10](q~value))
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(size = 10),
    legend.position = "bottom",
    legend.text     = element_text(size = 10)
  )

# plot generated in Script 05

# =============================================================================
# SECTION 9: SAVE RESULTS
# =============================================================================
save(consensus_3, two_of_three, ancom_out, linda_out, maas_out, maas_out_full,
     tax_lookup, df_meta_final, ps_genus, otu_mat,
     volcano_df, p_volcano,
     file = file.path(RDATA_DIR, "PD_vs_Control_DA_Results.RData"))

cat("\n--- SCRIPT 1.5 v1.4 COMPLETE ---\n")
cat("Results saved to", file.path(RDATA_DIR, "PD_vs_Control_DA_Results.RData"), "\n")
cat("Proceed to Script 2 for PDAQ cognitive association analyses.\n")