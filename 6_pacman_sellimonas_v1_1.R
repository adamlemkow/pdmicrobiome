# =============================================================================
# Script 06: Sellimonas Descriptive Characterization v1.1
# =============================================================================
# Changes from v1.0:
#   - All three ggsave() calls updated to save both PNG (300 dpi) and PDF (vector)
#   - bg = "white" added to all ggsave() calls for consistency
#   - Completion message updated to reflect PNG + PDF outputs
#
# Context:
#   Sellimonas (Lachnospiraceae, Firmicutes) was identified as the sole
#   2-of-2 consensus hit (MaAsLin2 + LinDA) in Script 03 v2.1, associating
#   higher baseline abundance with faster PDAQ cognitive decline (avg LFC
#   ~-0.67; sensitivity N=248: coef=-0.753, p<0.001, q=0.011).
#
#   This script characterizes that finding descriptively across four dimensions:
#     A. Prevalence and relative abundance — full cohort vs PD cases only
#     B. Prevalence and relative abundance — cognitively intact vs impaired PD
#     C. Correlation with continuous baseline PDAQ score
#     D. Publication-ready figures for A, B, and C
#
# Inputs:
#   rdata/PD_vs_Control_Cleaned.RData      (ps_clean — full cohort, from Script 01)
#   rdata/PD_PDAQ_Association_Data.RData   (ps_pdaq  — PDAQ subset, from Script 01)
#   rdata/03_Longitudinal_PDAQ_Results_v2_1.RData  (tax_lookup, slopes_df, from Script 03)
#
# Outputs:
#   tables/Table_06_Sellimonas_Prevalence.csv   (Sections A & B summary table)
#   tables/Table_06_Sellimonas_Correlation.csv  (Section C correlation results)
#   plots/Figure_06_Sellimonas_Cohort.png       (violin/jitter: full cohort vs PD)
#   plots/Figure_06_Sellimonas_CogStatus.png    (violin/jitter: intact vs impaired)
#   plots/Figure_06_Sellimonas_PDAQ.png         (scatter: abundance vs PDAQ score)
#
# Notes:
#   - Relative abundance is computed by TSS (proportional) from raw counts.
#   - log1p transformation applied only for correlation and visualisation;
#     raw TSS retained in summary tables for biological interpretability.
#   - Cognitive status threshold: PDAQ_Total >= 53 = intact (Levin et al.).
#   - Colourblind-safe Okabe-Ito palette consistent with pipeline convention.
# =============================================================================

# =============================================================================
# SECTION 1: INITIALIZATION
# =============================================================================

.libPaths("/arc/project/st-silkec-1/rsandboxlib")
setwd("/arc/project/st-silkec-1/pacman")

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# Project palette (Okabe-Ito — consistent with Scripts 01, 02, 05)
PAL_PD      <- "#D55E00"   # vermillion  — PD / impaired
PAL_CTRL    <- "#0072B2"   # blue        — control / intact
PAL_PDMCI   <- "#E69F00"   # orange      — not used here; retained for reference
PAL_OUTLIER <- "#CC79A7"   # pink        — not used here; retained for reference

PLOTS_DIR  <- "/arc/project/st-silkec-1/pacman/plots"
TABLES_DIR <- "/arc/project/st-silkec-1/pacman/tables"
RDATA_DIR  <- "/arc/project/st-silkec-1/pacman/rdata"

for (d in c(PLOTS_DIR, TABLES_DIR, RDATA_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# =============================================================================
# SECTION 2: LOAD DATA
# =============================================================================

cat("--- Loading RData objects ---\n")

load(file.path(RDATA_DIR, "PD_vs_Control_Cleaned.RData"))            # ps_clean
load(file.path(RDATA_DIR, "PD_PDAQ_Association_Data.RData"))          # ps_pdaq
load(file.path(RDATA_DIR, "03_Longitudinal_PDAQ_Results_v2.3.RData")) # tax_lookup, slopes_df, meta_baseline

cat("ps_clean samples (full cohort):", nsamples(ps_clean), "\n")
cat("ps_pdaq samples  (PDAQ subset):", nsamples(ps_pdaq),  "\n")

# =============================================================================
# SECTION 3: IDENTIFY SELLIMONAS TAXON ID
# =============================================================================

sellimonas_id <- tax_lookup %>%
  filter(grepl("Sellimonas", label, ignore.case = TRUE)) %>%
  pull(taxon)

if (length(sellimonas_id) == 0) stop("Sellimonas not found in tax_lookup — check taxonomy table.")
if (length(sellimonas_id) > 1)  {
  cat("WARNING: multiple Sellimonas entries found — using first.\n")
  print(tax_lookup %>% filter(grepl("Sellimonas", label, ignore.case = TRUE)))
  sellimonas_id <- sellimonas_id[1]
}

cat("Sellimonas taxon ID:", sellimonas_id, "\n")

# Taxonomy annotation for plot labels
sellimonas_label <- tax_lookup %>%
  filter(taxon == sellimonas_id) %>%
  mutate(full_label = paste0("Sellimonas (", Family, ", ", Phylum, ")")) %>%
  pull(full_label)

cat("Sellimonas annotation:", sellimonas_label, "\n")

# =============================================================================
# SECTION 4: HELPER FUNCTIONS
# =============================================================================

# Extract Sellimonas relative abundance from a phyloseq object.
# Returns a data frame with: sample_id, raw_counts, rel_abund, present.
# Gloms to genus level first to match the DA analysis.
extract_sellimonas <- function(ps) {
  ps_g  <- tax_glom(ps, "Genus")
  otu   <- as.data.frame(otu_table(ps_g))
  if (taxa_are_rows(ps_g)) otu <- t(otu)
  otu   <- as.data.frame(otu)
  
  if (!sellimonas_id %in% colnames(otu)) {
    stop("Sellimonas taxon ID not found in OTU table after tax_glom.")
  }
  
  rel   <- otu / rowSums(otu)
  
  data.frame(
    sample_id  = rownames(otu),
    raw_counts = otu[[sellimonas_id]],
    rel_abund  = rel[[sellimonas_id]],
    stringsAsFactors = FALSE
  ) %>%
    mutate(present = raw_counts > 0)
}

# Summary table for a grouping variable
summarise_sellimonas <- function(df, group_var) {
  df %>%
    filter(!is.na(.data[[group_var]])) %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      N                     = n(),
      n_present             = sum(present),
      pct_present           = round(mean(present) * 100, 1),
      median_rel_abund_pct  = round(median(rel_abund) * 100, 4),
      mean_rel_abund_pct    = round(mean(rel_abund)   * 100, 4),
      median_rel_pos_only   = round(median(rel_abund[present]) * 100, 4),
      .groups = "drop"
    ) %>%
    rename(Group = 1)
}

# =============================================================================
# SECTION 5A: FULL COHORT vs PD CASES
# =============================================================================

cat("\n===== SECTION A: Full cohort vs PD cases =====\n")

sell_full <- extract_sellimonas(ps_clean) %>%
  left_join(
    data.frame(sample_data(ps_clean)) %>%
      rownames_to_column("sample_id") %>%
      select(sample_id, parkinsons),
    by = "sample_id"
  )

cat("Total samples:", nrow(sell_full), "\n")
cat("PD cases:", sum(sell_full$parkinsons == "PD case", na.rm = TRUE), "\n")
cat("Controls:", sum(sell_full$parkinsons == "Control",  na.rm = TRUE), "\n")

summ_A <- summarise_sellimonas(sell_full, "parkinsons")
cat("\nPrevalence & relative abundance by PD status:\n")
print(summ_A)

wt_A <- wilcox.test(rel_abund ~ parkinsons, data = sell_full)
cat("Wilcoxon p (rel abund, PD vs Control):", round(wt_A$p.value, 4), "\n")

# =============================================================================
# SECTION 5B: INTACT vs IMPAIRED PD (cross-sectional PDAQ)
# =============================================================================

cat("\n===== SECTION B: Cognitively intact vs impaired PD =====\n")

meta_pdaq_df <- data.frame(sample_data(ps_pdaq)) %>%
  rownames_to_column("sample_id") %>%
  mutate(
    PDAQ_num   = as.numeric(as.character(PDAQ_Total)),
    cog_status = case_when(
      PDAQ_num >= 53 ~ "Intact (>=53)",
      PDAQ_num <  53 ~ "Impaired (<53)",
      TRUE           ~ NA_character_
    )
  ) %>%
  select(sample_id, parkinsons, PDAQ_num, cog_status)

sell_pdaq <- extract_sellimonas(ps_pdaq) %>%
  left_join(meta_pdaq_df, by = "sample_id")

# Restrict to PD cases for cog status comparison
sell_pd <- sell_pdaq %>% filter(parkinsons == "PD case")

cat("PD cases with PDAQ data:", nrow(sell_pd), "\n")
cat("Cognitively intact  (>=53):", sum(sell_pd$cog_status == "Intact (>=53)",  na.rm = TRUE), "\n")
cat("Cognitively impaired (<53):", sum(sell_pd$cog_status == "Impaired (<53)", na.rm = TRUE), "\n")
cat("Missing cog_status:", sum(is.na(sell_pd$cog_status)), "\n")

summ_B <- summarise_sellimonas(sell_pd, "cog_status")
cat("\nPrevalence & relative abundance by cognitive status:\n")
print(summ_B)

wt_B <- wilcox.test(rel_abund ~ cog_status,
                    data = filter(sell_pd, !is.na(cog_status)))
cat("Wilcoxon p (rel abund, intact vs impaired):", round(wt_B$p.value, 4), "\n")

# =============================================================================
# SECTION 5C: CORRELATION WITH CONTINUOUS PDAQ SCORE
# =============================================================================

cat("\n===== SECTION C: Correlation with continuous PDAQ score =====\n")

sell_corr <- sell_pd %>%
  filter(!is.na(PDAQ_num)) %>%
  mutate(log_rel = log1p(rel_abund))

cat("N with both Sellimonas and PDAQ score:", nrow(sell_corr), "\n")
cat("N with Sellimonas = 0 (zero-inflated):",
    sum(sell_corr$raw_counts == 0), "(", round(mean(sell_corr$raw_counts == 0)*100,1), "%)\n")

# Spearman — primary (non-parametric, appropriate given zero inflation)
sp <- cor.test(sell_corr$log_rel, sell_corr$PDAQ_num, method = "spearman")
cat("\nSpearman rho:", round(sp$estimate, 4), "  p =", round(sp$p.value, 4), "\n")

# Pearson — reported for completeness
pe <- cor.test(sell_corr$log_rel, sell_corr$PDAQ_num, method = "pearson")
cat("Pearson r:   ", round(pe$estimate, 4), "  p =", round(pe$p.value, 4), "\n")

# =============================================================================
# SECTION 6: EXPORT SUMMARY TABLES
# =============================================================================

# Combine A and B with a section label
tbl_prev <- bind_rows(
  summ_A %>% mutate(Comparison = "Full cohort by PD status",  .before = 1),
  summ_B %>% mutate(Comparison = "PD cases by cognitive status", .before = 1)
)

write.csv(tbl_prev,
          file.path(TABLES_DIR, "Table_06_Sellimonas_Prevalence.csv"),
          row.names = FALSE, na = "—")
cat("\nExported: tables/Table_06_Sellimonas_Prevalence.csv\n")

corr_tbl <- data.frame(
  Method      = c("Spearman", "Pearson"),
  Statistic   = c(round(sp$estimate, 4), round(pe$estimate, 4)),
  p_value     = c(round(sp$p.value,   4), round(pe$p.value,   4)),
  N           = nrow(sell_corr),
  Outcome     = "Baseline PDAQ_Total (continuous)",
  Predictor   = "log1p(Sellimonas relative abundance)"
)

write.csv(corr_tbl,
          file.path(TABLES_DIR, "Table_06_Sellimonas_Correlation.csv"),
          row.names = FALSE, na = "—")
cat("Exported: tables/Table_06_Sellimonas_Correlation.csv\n")

# =============================================================================
# SECTION 7: FIGURES
# =============================================================================

# --- 7A: Full cohort vs PD cases (violin + jitter) ---
fig_A <- sell_full %>%
  filter(!is.na(parkinsons)) %>%
  mutate(
    log_rel    = log1p(rel_abund),
    parkinsons = factor(parkinsons, levels = c("Control", "PD case"))
  ) %>%
  ggplot(aes(x = parkinsons, y = log_rel, fill = parkinsons)) +
  geom_violin(trim = TRUE, alpha = 0.7, colour = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4,
              aes(colour = parkinsons)) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     label.x.npc = "center",
                     label.y.npc = "top") +
  scale_fill_manual(values  = c("Control" = PAL_CTRL, "PD case" = PAL_PD)) +
  scale_colour_manual(values = c("Control" = PAL_CTRL, "PD case" = PAL_PD)) +
  labs(
    x       = NULL,
    y       = "Sellimonas abundance (log1p relative abundance)",
    title   = "Sellimonas prevalence: Controls vs PD cases",
    caption = paste0("Points = individual participants. Wilcoxon rank-sum test.\n",
                     "N = ", sum(sell_full$parkinsons == "Control", na.rm=TRUE),
                     " controls, ",
                     sum(sell_full$parkinsons == "PD case", na.rm=TRUE), " PD cases.")
  ) +
  theme_pubr() +
  theme(legend.position = "none")

ggsave(file.path(PLOTS_DIR, "Figure_06_Sellimonas_Cohort.png"),
       fig_A, width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave(file.path(PLOTS_DIR, "Figure_06_Sellimonas_Cohort.pdf"),
       fig_A, width = 5, height = 5, units = "in", bg = "white")
cat("Saved: Figure_06_Sellimonas_Cohort (PNG + PDF)\n")

# --- 7B: Intact vs Impaired PD (violin + jitter) ---
fig_B <- sell_pd %>%
  filter(!is.na(cog_status)) %>%
  mutate(
    log_rel    = log1p(rel_abund),
    cog_status = factor(cog_status, levels = c("Intact (>=53)", "Impaired (<53)"))
  ) %>%
  ggplot(aes(x = cog_status, y = log_rel, fill = cog_status)) +
  geom_violin(trim = TRUE, alpha = 0.7, colour = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4,
              aes(colour = cog_status)) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     label.x.npc = "center",
                     label.y.npc = "top") +
  scale_fill_manual(values  = c("Intact (>=53)" = PAL_CTRL,
                                "Impaired (<53)" = PAL_PD)) +
  scale_colour_manual(values = c("Intact (>=53)" = PAL_CTRL,
                                 "Impaired (<53)" = PAL_PD)) +
  labs(
    x       = NULL,
    y       = "Sellimonas abundance (log1p relative abundance)",
    title   = "Sellimonas: cognitively intact vs impaired PD",
    caption = paste0("PD cases only. PDAQ ≥53 = intact (Levin et al.).\n",
                     "N = ",
                     sum(sell_pd$cog_status == "Intact (>=53)",  na.rm=TRUE),
                     " intact, ",
                     sum(sell_pd$cog_status == "Impaired (<53)", na.rm=TRUE),
                     " impaired.")
  ) +
  theme_pubr() +
  theme(legend.position = "none")

ggsave(file.path(PLOTS_DIR, "Figure_06_Sellimonas_CogStatus.png"),
       fig_B, width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave(file.path(PLOTS_DIR, "Figure_06_Sellimonas_CogStatus.pdf"),
       fig_B, width = 5, height = 5, units = "in", bg = "white")
cat("Saved: Figure_06_Sellimonas_CogStatus (PNG + PDF)\n")

# --- 7C: Continuous PDAQ scatter ---
# Annotate points by cognitive status for visual context
fig_C <- sell_corr %>%
  mutate(
    cog_status = factor(
      ifelse(is.na(cog_status), "Unknown", cog_status),
      levels = c("Intact (>=53)", "Impaired (<53)", "Unknown")
    )
  ) %>%
  ggplot(aes(x = log_rel, y = PDAQ_num, colour = cog_status)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "grey30", linewidth = 0.8,
              inherit.aes = FALSE,
              aes(x = log_rel, y = PDAQ_num)) +
  stat_cor(method = "spearman",
           inherit.aes = FALSE,
           aes(x = log_rel, y = PDAQ_num),
           label.x.npc = "left",
           label.y.npc = "top") +
  scale_colour_manual(
    values = c("Intact (>=53)"  = PAL_CTRL,
               "Impaired (<53)" = PAL_PD,
               "Unknown"        = "grey60"),
    name = "Cognitive status"
  ) +
  labs(
    x       = "Sellimonas abundance (log1p relative abundance)",
    y       = "Baseline PDAQ score",
    title   = "Sellimonas abundance vs baseline PDAQ score",
    caption = paste0("PD cases only. N = ", nrow(sell_corr), ". ",
                     "Spearman rho = ", round(sp$estimate, 3),
                     ", p = ", round(sp$p.value, 4), ".\n",
                     "PDAQ ≥53 = cognitively intact (Levin et al.). ",
                     "Horizontal dashed line at threshold."),
    colour  = "Cognitive status"
  ) +
  geom_hline(yintercept = 53, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  theme_pubr() +
  theme(legend.position = "right")

ggsave(file.path(PLOTS_DIR, "Figure_06_Sellimonas_PDAQ.png"),
       fig_C, width = 7, height = 5, units = "in", dpi = 300, bg = "white")
ggsave(file.path(PLOTS_DIR, "Figure_06_Sellimonas_PDAQ.pdf"),
       fig_C, width = 7, height = 5, units = "in", bg = "white")
cat("Saved: Figure_06_Sellimonas_PDAQ (PNG + PDF)\n")

# =============================================================================
# SECTION 8: CONSOLE SUMMARY
# =============================================================================

cat("\n")
cat("=============================================================\n")
cat("  SCRIPT 06 v1.1 COMPLETE — SELLIMONAS CHARACTERIZATION     \n")
cat("=============================================================\n\n")

cat("--- A: Full cohort by PD status ---\n")
print(summ_A)
cat("Wilcoxon p:", round(wt_A$p.value, 4), "\n\n")

cat("--- B: PD cases by cognitive status ---\n")
print(summ_B)
cat("Wilcoxon p:", round(wt_B$p.value, 4), "\n\n")

cat("--- C: Correlation with continuous PDAQ (PD cases) ---\n")
cat("N =", nrow(sell_corr), "| Zero-inflated:",
    round(mean(sell_corr$raw_counts == 0)*100, 1), "%\n")
cat("Spearman rho =", round(sp$estimate, 4), "  p =", round(sp$p.value, 4), "\n")
cat("Pearson   r  =", round(pe$estimate, 4), "  p =", round(pe$p.value, 4), "\n\n")

cat("Outputs saved to:\n")
cat("  tables/Table_06_Sellimonas_Prevalence.csv\n")
cat("  tables/Table_06_Sellimonas_Correlation.csv\n")
cat("  plots/Figure_06_Sellimonas_Cohort.png + .pdf\n")
cat("  plots/Figure_06_Sellimonas_CogStatus.png + .pdf\n")
cat("  plots/Figure_06_Sellimonas_PDAQ.png + .pdf\n")