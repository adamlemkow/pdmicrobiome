# Script 01: PACMAN diversity analysis v1.9
# Changes from v1.8:
#   - amantadine added to fi_covars select — carries new column from Script 00.5
#     v1.3 through into ps_pdaq metadata for Script 02

# =============================================================================

.libPaths("/arc/project/st-silkec-1/rsandboxlib")

library(vegan)
library(rstatix)
library(dplyr)
library(tidyr)
library(phyloseq)
library(reshape2)
library(ggplot2)
library(ggpubr)

setwd("/arc/project/st-silkec-1/pacman")

# =============================================================================
# PROJECT PALETTE (Okabe-Ito — colourblind-safe across all scripts)
# PD case / Cognitively impaired  →  #D55E00  vermillion
# Control / Cognitively intact    →  #0072B2  blue
# PD-MCI (three-level only)       →  #E69F00  orange
# Sequencing outliers             →  #CC79A7  pink
# =============================================================================
PAL_PD      <- "#D55E00"
PAL_CTRL    <- "#0072B2"
PAL_OUTLIER <- "#CC79A7"

# =============================================================================
# OUTPUT DIRECTORIES
# =============================================================================
RDATA_DIR  <- "/arc/project/st-silkec-1/pacman/rdata"
PLOTS_DIR  <- "/arc/project/st-silkec-1/pacman/plots"
dir.create(RDATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# --- SECTION 1: Load and Filter Metadata ---
# 1.1 Load raw survey data
# Stagman is retained as the metadata source for Scripts 01 and 01.5 (PD vs Control).
# FoxInsight metadata (Script 00.5 output) is used from Script 02 onwards.
meta_raw <- read.delim("/arc/project/st-silkec-1/pacman/other_foxden/Stagman_2024_Survey_Data_2024-11-15.tsv", sep="\t")

# 1.2 Filter for stool samples and valid IDs
# Antibiotic exclusion applied here to match Script 00.5 logic:
#   PDMBAntiBio == 1 (confirmed use) → excluded
#   PDMBAntiBio == 2 (unsure)        → retained (conservative inclusion)
#   PDMBAntiBio == 999 / NA          → retained
n_pre_abx  <- nrow(meta_raw %>% filter(sample_type == "Stool"))
meta_stool <- meta_raw %>%
  filter(sample_type == "Stool") %>%
  filter(!is.na(FoxDEN_ID) & FoxDEN_ID != "" & FoxDEN_ID != " ") %>%
  filter(!is.na(parkinsons) & parkinsons != "") %>%
  filter(is.na(antibiotic_six_months) | antibiotic_six_months != "yes") %>%
  distinct(FoxDEN_ID, .keep_all = TRUE)

cat(sprintf("Antibiotic exclusion (Script 01): %d removed, %d retained\n",
            n_pre_abx - nrow(meta_stool), nrow(meta_stool)))

rownames(meta_stool) <- meta_stool$FoxDEN_ID

# --- SECTION 2: Load OTU Data ---
# 2.1 Load OTU table and handle NAs
otu_wide <- read.table("/arc/project/st-silkec-1/pacman/other_foxden/Pdmb_Otu_Data_2024-11-15.tsv",
                       header = TRUE, check.names = FALSE)
otu_wide[is.na(otu_wide)] = 0

# --- SECTION 3: Align OTU Data with Metadata ---
# 3.1 Subset OTU data to stool samples
otu_stool_only <- otu_wide %>%
  filter(FoxDEN_ID %in% meta_stool$FoxDEN_ID) %>%
  distinct(FoxDEN_ID, .keep_all = TRUE)

rownames(otu_stool_only) <- otu_stool_only$FoxDEN_ID
otu_numeric <- as.matrix(otu_stool_only[, grepl("GCF_", colnames(otu_stool_only))])

# 3.2 Find common samples between OTU and Metadata
common_samples <- intersect(rownames(otu_numeric), rownames(meta_stool))
otu_final_samples <- otu_numeric[common_samples, ]

# --- SECTION 4: Load and Align Taxonomy ---
# 4.1 Load assignments and match with OTU features
tax_raw <- read.delim("/arc/project/st-silkec-1/pacman/other_foxden/Otu_Assignments_2024-11-15.tsv", sep="\t")
rownames(tax_raw) <- tax_raw$OTU
common_taxa <- intersect(colnames(otu_final_samples), rownames(tax_raw))

otu_final_aligned <- otu_final_samples[, common_taxa]
tax_final_aligned <- as.matrix(tax_raw[common_taxa, -c(1,2)])

# --- SECTION 5: Assemble Phyloseq Object & Report Full Cohort ---
# 5.1 Create main phyloseq object
ps <- phyloseq(
  otu_table(otu_final_aligned, taxa_are_rows = FALSE),
  sample_data(meta_stool[common_samples, ]),
  tax_table(tax_final_aligned)
)

# 5.2 Report full cohort sample sizes
# These numbers describe everyone with a stool sample and valid metadata before
# any outlier removal or PDAQ subsetting. Use these as your "full cohort" N in-text.
cat("\n--- FULL COHORT SAMPLE SIZES (stool samples with metadata) ---\n")
cat("Total samples:", nsamples(ps), "\n")
cat("Total taxa (OTUs):", ntaxa(ps), "\n")
cat("\nBreakdown by PD status:\n")
print(table(sample_data(ps)$parkinsons))
cat("\nBreakdown by sex:\n")
print(table(sample_data(ps)$sex))
cat("\nBreakdown by batch (main vs pilot):\n")
print(table(sample_data(ps)$main_or_pilot_study))

# 5.3 Rarefy for diversity analyses
# Rarefaction normalises sequencing depth across samples. Samples below the
# minimum depth threshold are dropped automatically — we report how many here.
ps_rare <- rarefy_even_depth(ps, rngseed = 421)
cat("\nSamples retained after rarefaction:", nsamples(ps_rare), "\n")
cat("(Samples dropped due to low depth:", nsamples(ps) - nsamples(ps_rare), ")\n")

# =============================================================================
# SECTION 6: BETA DIVERSITY & OUTLIER DETECTION
# Beta diversity runs first so outliers can be identified and removed before
# alpha diversity is calculated. This ensures all analyses use the same
# cleaned sample set.
# =============================================================================

# 6.1 Initial PCoA on full rarefied dataset — identify outliers visually
otu_math_matrix <- as.matrix(otu_table(ps))
bray_dist <- vegan::vegdist(otu_math_matrix, method = "bray")
pcoa_res  <- ordinate(ps, method = "PCoA", distance = bray_dist)

pcoa_df <- data.frame(pcoa_res$vectors[, 1:2],
                      parkinsons  = sample_data(ps)$parkinsons,
                      LibrarySize = sample_sums(ps))

line_sample_names <- rownames(pcoa_df)[pcoa_df$Axis.1 > 0.4]
main_sample_names <- setdiff(sample_names(ps), line_sample_names)

cat("\nOutliers identified (Axis.1 > 0.4):", length(line_sample_names), "\n")
print(line_sample_names)

# 6.2 Plot initial PCoA with outlier threshold marked
pcoa_initial_plot <- plot_ordination(ps, pcoa_res, color = "parkinsons") +
  geom_point(alpha = 0.5, size = 2) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = PAL_OUTLIER) +
  annotate("text", x = 0.4, y = 0, label = "Outlier threshold",
           color = PAL_OUTLIER, hjust = 1, size = 3.5) +
  scale_color_manual(values = c("Control" = PAL_CTRL, "PD case" = PAL_PD)) +
  labs(title    = "Initial PCoA — Full Cohort (Pre-outlier Removal)",
       subtitle = paste0("N = ", nsamples(ps))) +
  theme_pubr()
# plot generated in Script 05

# 6.3 Justify outlier removal via sequencing depth
# Outlier samples are expected to have anomalously high sequencing depth,
# which distorts Bray-Curtis distances and pulls them away from the main cluster.
depth_df <- data.frame(
  SampleID    = sample_names(ps),
  LibrarySize = sample_sums(ps),   # pre-rarefaction raw counts
  Group       = ifelse(sample_names(ps) %in% line_sample_names, "Outlier", "Main Cluster")
)

depth_wilcox <- wilcox.test(LibrarySize ~ Group, data = depth_df)
depth_p      <- signif(depth_wilcox$p.value, 3)
depth_label  <- paste0("Wilcoxon p = ", depth_p)

depth_plot <- ggboxplot(depth_df, x = "Group", y = "LibrarySize",
                        fill = "Group", palette = c(PAL_CTRL, PAL_OUTLIER),
                        add = "jitter", add.params = list(alpha = 0.3),
                        title    = "Sequencing Depth: Outliers vs. Main Cluster (Pre-rarefaction)",
                        subtitle = paste0("N = ", nrow(depth_df),
                                          "  (outliers: ", sum(depth_df$Group == "Outlier"), ")")) +
  annotate("text", x = 1.5,
           y = max(depth_df$LibrarySize, na.rm = TRUE) * 1.05,
           label = depth_label, size = 4)
# plot generated in Script 05

# 6.4 Final beta diversity on cleaned dataset (outliers removed)
# Raw Bray-Curtis used (no log1p) for consistency across all scripts.
ps_main        <- prune_samples(main_sample_names, ps_rare)
otu_main_clean <- as.matrix(otu_table(ps_main))
bray_final     <- vegan::vegdist(otu_main_clean, method = "bray")
pcoa_final     <- ordinate(ps_main, method = "PCoA", distance = bray_final)
metadata_final <- as(sample_data(ps_main), "data.frame")
permanova_final <- vegan::adonis2(bray_final ~ parkinsons, data = metadata_final, permutations = 999)

cat("\nPERMANOVA — Full Cohort PD vs Control (post-outlier removal):\n")
print(permanova_final)

# Extract PERMANOVA stats for plot annotation
perm_r2      <- round(permanova_final["Model", "R2"], 3)
perm_p       <- permanova_final["Model", "Pr(>F)"]
perm_p_label <- ifelse(perm_p < 0.001, "p < 0.001", paste0("p = ", signif(perm_p, 3)))
perm_label   <- paste0("PERMANOVA: R² = ", perm_r2, ", ", perm_p_label)

pcoa_final_plot <- plot_ordination(ps_main, pcoa_final, color = "parkinsons") +
  geom_point(alpha = 0.5, size = 2) +
  stat_ellipse(aes(group = parkinsons), linetype = 2) +
  scale_color_manual(values = c("Control" = PAL_CTRL, "PD case" = PAL_PD)) +
  labs(title    = "Beta Diversity PCoA — Full Cohort (PD vs Control, Cleaned)",
       subtitle = paste0("N = ", nsamples(ps_main)),
       caption  = perm_label) +
  theme_pubr() +
  theme(plot.caption = element_text(hjust = 0.5, size = 10, face = "italic"))
# plot generated in Script 05

# =============================================================================
# SECTION 7: ALPHA DIVERSITY
# Runs on ps_main (outlier-cleaned) to match the beta diversity sample set.
# =============================================================================

# 7.1 Estimate richness metrics
alpha_df <- estimate_richness(ps_main, measures = c("Shannon", "Observed", "Chao1"))

meta_alpha <- data.frame(
  sample_data(ps_main),
  Chao1    = alpha_df$Chao1,
  Shannon  = alpha_df$Shannon,
  Observed = alpha_df$Observed,
  check.names = FALSE
)

# 7.2 Statistical testing (Wilcoxon rank-sum, PD vs Control)
alpha_long <- meta_alpha %>%
  tidyr::pivot_longer(cols = c("Chao1", "Shannon", "Observed"),
                      names_to = "Metric", values_to = "Value") %>%
  filter(!is.na(parkinsons))

alpha_stats <- alpha_long %>%
  group_by(Metric) %>%
  wilcox_test(Value ~ parkinsons) %>%
  add_significance() %>%
  mutate(y.position = c(
    max(alpha_long$Value[alpha_long$Metric == "Chao1"],    na.rm = TRUE) * 1.05,
    max(alpha_long$Value[alpha_long$Metric == "Observed"], na.rm = TRUE) * 1.05,
    max(alpha_long$Value[alpha_long$Metric == "Shannon"],  na.rm = TRUE) * 1.05
  ))

# 7.3 Plot alpha diversity
alpha_plot <- ggboxplot(alpha_long, x = "parkinsons", y = "Value", fill = "parkinsons",
                        palette = c("Control" = PAL_CTRL, "PD case" = PAL_PD),
                        add = "jitter",
                        add.params = list(alpha = 0.2), facet.by = "Metric",
                        scales = "free_y",
                        title    = "Alpha Diversity — Full Cohort (PD vs Control)",
                        subtitle = paste0("N = ", nsamples(ps_main))) +
  stat_pvalue_manual(alpha_stats, label = "p.signif", tip.length = 0.01) +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12))
# plot generated in Script 05

# =============================================================================
# SECTION 8: DATA EXPORT AND PDAQ INTEGRATION
# =============================================================================

# 8.1 Clean Sample IDs (standardize to numeric) using un-rarefied ps for Script 02
ps_clean  <- prune_samples(main_sample_names, ps)
clean_ids <- gsub("[^0-9]", "", sample_names(ps_clean))
clean_ids <- as.character(as.numeric(clean_ids))
sample_names(ps_clean) <- clean_ids

# 8.2 Save Track A: PD vs Control (full cleaned cohort, un-rarefied)
# Also saves plot objects and supporting data for Script 05 figure generation
save(ps_clean, ps_rare,
     pcoa_res, pcoa_df, line_sample_names,
     pcoa_initial_plot,
     depth_df, depth_wilcox, depth_p, depth_label,
     depth_plot,
     pcoa_final, permanova_final, perm_r2, perm_p, perm_label,
     pcoa_final_plot,
     alpha_long, alpha_stats, alpha_plot,
     file = file.path(RDATA_DIR, "PD_vs_Control_Cleaned.RData"))

# 8.3 Load PDAQ, fix IDs, and calculate total score
pdaq_data <- read.csv("/arc/project/st-silkec-1/pacman/other_foxden/foxden_pdaq15_raw.csv")

# 8.3.1 ID repair: if ID column is NA, pull from row names
if (all(is.na(pdaq_data$fox_insight_id))) {
  pdaq_data <- read.csv("/arc/project/st-silkec-1/pacman/other_foxden/foxden_pdaq15_raw.csv", row.names = 1)
  pdaq_data$fox_insight_id <- rownames(pdaq_data)
}

# 8.3.2 Calculate PDAQ_Total from the 15 daily items
# Value 5 (prefer not to answer) is recoded to NA before summing — it is
# not a valid difficulty rating and would inflate scores above the 60-point
# maximum (15 items x 4 = 60). na.rm = FALSE then enforces strict
# completeness: any participant with a missing or PNTA item receives NA
# and is excluded. This avoids inflated or partial scores.
daily_cols <- c("DailyRead", "DailyClock", "DailyMoney", "DailyInstruct",
                "DailyProblem", "DailyExplain", "DailyErrand", "DailyMap",
                "DailyNumber", "DailyMany", "DailyLearn", "DailyFinance",
                "DailyThought", "DailyDiscuss", "DailyRememb")

pdaq_items <- pdaq_data[, daily_cols]
pdaq_items[pdaq_items == 5] <- NA
pdaq_data$PDAQ_Total <- rowSums(pdaq_items, na.rm = FALSE)

# 8.3.3 Clean IDs
pdaq_data <- pdaq_data[!is.na(pdaq_data$PDAQ_Total), ]
pdaq_data$clean_id <- as.character(as.numeric(gsub("[^0-9]", "", pdaq_data$fox_insight_id)))

pdaq_data <- pdaq_data %>%
  mutate(
    months_since_baseline = case_when(
      trimws(as.character(schedule_of_activities)) == "REG" ~ 0,
      !is.na(suppressWarnings(as.numeric(as.character(schedule_of_activities)))) ~
        suppressWarnings(as.numeric(as.character(schedule_of_activities))),
      TRUE ~ NA_real_
    )
  )

# 8.3.4 Filter to stool cohort FIRST, then do baseline selection
# This prevents match() from grabbing wrong visits from unrelated participants
pdaq_data <- pdaq_data %>%
  filter(clean_id %in% sample_names(ps_clean))

pdaq_reg <- pdaq_data %>%
  filter(trimws(as.character(schedule_of_activities)) == "REG") %>%
  distinct(clean_id, .keep_all = TRUE) %>%
  mutate(schedule_interpolated = 0L)

ids_with_reg <- unique(pdaq_reg$clean_id)

pdaq_fallback <- pdaq_data %>%
  filter(!clean_id %in% ids_with_reg) %>%
  group_by(clean_id) %>%
  arrange(days_elapsed) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    months_since_baseline = 0,
    schedule_interpolated = 1L
  )

pdaq_data <- bind_rows(pdaq_reg, pdaq_fallback)

cat("\n--- PDAQ BASELINE SELECTION ---\n")
cat("Confirmed baseline (REG):", sum(pdaq_data$schedule_interpolated == 0), "\n")
cat("Assumed baseline (days_elapsed fallback):", sum(pdaq_data$schedule_interpolated == 1), "\n")
cat("Total PDAQ entries carried forward:", nrow(pdaq_data), "\n")

# 8.3.5 Intersect with stool sample IDs and prune
common_pdaq_ids <- intersect(sample_names(ps_clean), pdaq_data$clean_id)
cat("\nFinal validated participants (complete PDAQ, non-duplicate):", length(common_pdaq_ids), "\n")

if (length(common_pdaq_ids) == 0) {
  stop("CRITICAL ERROR: No matching IDs remain after filtering.")
}

ps_pdaq <- prune_samples(common_pdaq_ids, ps_clean)

# 8.4 Map PDAQ_Total, months_since_baseline, and schedule_interpolated into phyloseq metadata
meta_pdaq   <- data.frame(sample_data(ps_pdaq))
score_index <- match(sample_names(ps_pdaq), pdaq_data$clean_id)
meta_pdaq$PDAQ_Total             <- pdaq_data$PDAQ_Total[score_index]
meta_pdaq$months_since_baseline  <- pdaq_data$months_since_baseline[score_index]
meta_pdaq$schedule_interpolated  <- pdaq_data$schedule_interpolated[score_index]
sample_data(ps_pdaq) <- sample_data(meta_pdaq)

# 8.4.5 Merge FoxInsight covariates into ps_pdaq
# Stagman is used for PD vs Control analyses (Scripts 01, 01.5) but lacks
# antibiotic data and has less complete covariate coverage than FoxInsight.
# For PDAQ analyses (Script 02, PD-only), we replace Stagman-derived covariates
# with FoxInsight versions from Script 00.5: BMI, disease_duration, azilect,
# levodopa, Bristol_Avg, abx_certain_no.
# Participants not in meta_clean (controls, unmatched IDs) retain Stagman values.

load("/arc/project/st-silkec-1/pacman/microbiome_metadata/foxinsight_metadata_clean.RData")

fi_covars <- meta_clean %>%
  mutate(clean_id = as.character(as.numeric(gsub("[^0-9]", "", FoxDEN_ID)))) %>%
  select(clean_id, BMI, disease_duration, azilect, amantadine, levodopa,
         Bristol_Avg, abx_certain_no, current_age)

meta_pdaq_fi <- data.frame(sample_data(ps_pdaq)) %>%
  tibble::rownames_to_column("clean_id") %>%
  left_join(fi_covars, by = "clean_id", suffix = c("_stagman", "_fi")) %>%
  mutate(
    BMI         = coalesce(BMI_fi,         BMI_stagman),
    current_age = coalesce(current_age_fi, current_age_stagman)
  ) %>%
  select(-BMI_stagman, -BMI_fi, -current_age_stagman, -current_age_fi)

rownames(meta_pdaq_fi) <- meta_pdaq_fi$clean_id
meta_pdaq_fi$clean_id  <- NULL

sample_data(ps_pdaq) <- sample_data(meta_pdaq_fi)

# 8.5 Report PDAQ subset sample sizes
# These are participants with both a valid stool sample AND a complete PDAQ score.
# This is the cohort used in Script 02.
cat("\n--- PDAQ SUBSET SAMPLE SIZES ---\n")
cat("Total samples with PDAQ scores:", nsamples(ps_pdaq), "\n")
cat("\nBaseline assignment method:\n")
print(table(sample_data(ps_pdaq)$schedule_interpolated,
            dnn = "schedule_interpolated (0=REG confirmed, 1=days_elapsed fallback)"))
cat("\nBreakdown by PD status:\n")
print(table(sample_data(ps_pdaq)$parkinsons))
cat("\nBreakdown by sex:\n")
print(table(sample_data(ps_pdaq)$sex))
cat("\nBreakdown by batch:\n")
print(table(sample_data(ps_pdaq)$main_or_pilot_study))
cat("\nPDAQ Total score summary (all participants):\n")
print(summary(as.numeric(sample_data(ps_pdaq)$PDAQ_Total)))
cat("\nPDAQ Total score by PD status:\n")
pdaq_by_group <- data.frame(sample_data(ps_pdaq)) %>%
  mutate(PDAQ_num = as.numeric(as.character(PDAQ_Total))) %>%
  group_by(parkinsons) %>%
  summarise(N      = n(),
            Mean   = round(mean(PDAQ_num,   na.rm = TRUE), 1),
            SD     = round(sd(PDAQ_num,     na.rm = TRUE), 1),
            Median = median(PDAQ_num,        na.rm = TRUE),
            Min    = min(PDAQ_num,           na.rm = TRUE),
            Max    = max(PDAQ_num,           na.rm = TRUE))
print(pdaq_by_group)

# Wilcoxon test: do PDAQ scores differ between PD and Control?
# Expected: PD cases score lower (more functional impairment). Confirms construct validity.
pdaq_meta_df <- data.frame(sample_data(ps_pdaq)) %>%
  mutate(PDAQ_num = as.numeric(as.character(PDAQ_Total)))
cat("\nWilcoxon test — PDAQ scores PD vs Control:\n")
print(wilcox.test(PDAQ_num ~ parkinsons, data = pdaq_meta_df))

# 8.6 Save deferred to after Section 9 — ps_pdaq_pd, ps_pdaq_rare, and all
# PDAQ diversity objects (alpha_pdaq_*, pcoa_pdaq*, permanova_pdaq*) are
# created in Section 9 and must exist before the save() call.

# =============================================================================
# SECTION 9: PDAQ SUBSET — REPRESENTATIVENESS CHECK & DIVERSITY REANALYSIS
# =============================================================================
# Narrative: We confirm the PDAQ subset is representative of the full cohort
# before presenting it as the basis for cognitive association analyses in Script 02.
# No outlier removal is applied — the full-cohort PCoA threshold (Axis.1 > 0.38)
# is not valid on the subset as PCoA axes will shift with a different sample set.

cat("\n--- SECTION 9: PDAQ SUBSET REANALYSIS ---\n")
cat("Full cohort N:", nsamples(ps_clean), "\n")
cat("PDAQ subset N:", nsamples(ps_pdaq), "\n")

# --- 9.1 REPRESENTATIVENESS CHECK ---
# NOTE: This check is not applicable in the current pipeline.
# Because Script 00.5 uses foxinsight_pdmb.csv as the metadata anchor, every
# participant in the pipeline has a stool sample collected at the same visit as
# their PDMB survey. PDAQ (DailyActivity survey) was collected as part of the
# same FoxInsight enrollment process. As a result, all 366 participants in
# ps_clean also have complete PDAQ data — ps_clean and ps_pdaq are identical.
# A formal representativeness comparison requires a larger "full cohort" from
# which the PDAQ subset was drawn; that distinction does not exist here.
# If a future version of this pipeline includes participants without PDAQ data
# (e.g. by relaxing the PDMB anchor), this section should be reinstated.

# cat_test helper retained for potential future use
cat_test <- function(tbl, label) {
  cat("\n--- Representativeness:", label, "---\n")
  expected <- chisq.test(tbl)$expected
  if (any(expected < 5)) {
    cat("(Expected cell count < 5 detected — using Fisher's exact test)\n")
    print(fisher.test(tbl))
  } else {
    print(chisq.test(tbl))
  }
}

# --- 9.2 Rarefy PDAQ subset — PD cases only ---
# Filter to PD cases before rarefaction — cognitive status comparison is PD-only
ps_pdaq_pd <- prune_samples(
  sample_data(ps_pdaq)$parkinsons == "PD case", ps_pdaq
)

# Derive cog_intact before rarefaction so it's available in sample_data
meta_pd <- data.frame(sample_data(ps_pdaq_pd)) %>%
  mutate(
    PDAQ_num   = as.numeric(as.character(PDAQ_Total)),
    cog_intact = factor(ifelse(PDAQ_num >= 53, "Intact (≥53)", "Impaired (<53)"),
                        levels = c("Impaired (<53)", "Intact (≥53)"))
  )
rownames(meta_pd) <- rownames(data.frame(sample_data(ps_pdaq_pd)))
sample_data(ps_pdaq_pd) <- sample_data(meta_pd)

ps_pdaq_rare <- rarefy_even_depth(ps_pdaq_pd, rngseed = 421)

cat("\nPD cases carried into Section 9:", nsamples(ps_pdaq_pd), "\n")
cat("Cognitive intact (≥53):", sum(meta_pd$cog_intact == "Intact (≥53)", na.rm = TRUE), "\n")
cat("Cognitive impaired (<53):", sum(meta_pd$cog_intact == "Impaired (<53)", na.rm = TRUE), "\n")

# --- 9.3 ALPHA DIVERSITY: Intact vs Impaired in PD Cohort ---
alpha_pdaq_df <- estimate_richness(ps_pdaq_rare, measures = c("Shannon", "Observed", "Chao1"))

meta_alpha_pdaq <- data.frame(
  sample_data(ps_pdaq_rare),
  Chao1    = alpha_pdaq_df$Chao1,
  Shannon  = alpha_pdaq_df$Shannon,
  Observed = alpha_pdaq_df$Observed,
  check.names = FALSE
)

alpha_pdaq_long <- meta_alpha_pdaq %>%
  tidyr::pivot_longer(cols = c("Chao1", "Shannon", "Observed"),
                      names_to = "Metric", values_to = "Value") %>%
  filter(!is.na(cog_intact))

alpha_pdaq_stats <- alpha_pdaq_long %>%
  group_by(Metric) %>%
  wilcox_test(Value ~ cog_intact) %>%
  add_significance() %>%
  mutate(y.position = c(
    max(alpha_pdaq_long$Value[alpha_pdaq_long$Metric == "Chao1"],    na.rm = TRUE) * 1.05,
    max(alpha_pdaq_long$Value[alpha_pdaq_long$Metric == "Observed"], na.rm = TRUE) * 1.05,
    max(alpha_pdaq_long$Value[alpha_pdaq_long$Metric == "Shannon"],  na.rm = TRUE) * 1.05
  ))

alpha_pdaq_plot <- ggboxplot(alpha_pdaq_long, x = "cog_intact", y = "Value",
                             fill = "cog_intact",
                             palette = c("Impaired (<53)" = PAL_PD, "Intact (≥53)" = PAL_CTRL),
                             add = "jitter", add.params = list(alpha = 0.2),
                             facet.by = "Metric", scales = "free_y",
                             title    = "Alpha Diversity — PD Cohort (Intact vs Impaired Cognition)",
                             subtitle = paste0("N = ", nsamples(ps_pdaq_rare))) +
  stat_pvalue_manual(alpha_pdaq_stats, label = "p.signif", tip.length = 0.01) +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12))
# plot generated in Script 05

# --- 9.4 BETA DIVERSITY: Intact vs Impaired in PD Cohort ---
# Raw Bray-Curtis used (no log1p) for consistency across all scripts.
otu_pdaq_mat <- as.matrix(otu_table(ps_pdaq_rare))
bray_pdaq    <- vegan::vegdist(otu_pdaq_mat, method = "bray")
pcoa_pdaq    <- ordinate(ps_pdaq_rare, method = "PCoA", distance = bray_pdaq)
meta_pdaq_df <- as(sample_data(ps_pdaq_rare), "data.frame")

permanova_pdaq <- vegan::adonis2(bray_pdaq ~ cog_intact, data = meta_pdaq_df, permutations = 999)
cat("\n--- PERMANOVA: Intact vs Impaired Cognition (PD Cohort) ---\n")
print(permanova_pdaq)

perm_pdaq_r2      <- round(permanova_pdaq["Model", "R2"], 3)
perm_pdaq_p       <- permanova_pdaq["Model", "Pr(>F)"]
perm_pdaq_p_label <- ifelse(perm_pdaq_p < 0.001, "p < 0.001", paste0("p = ", signif(perm_pdaq_p, 3)))
perm_pdaq_label   <- paste0("PERMANOVA: R² = ", perm_pdaq_r2, ", ", perm_pdaq_p_label)

pcoa_pdaq_plot <- plot_ordination(ps_pdaq_rare, pcoa_pdaq, color = "cog_intact") +
  geom_point(alpha = 0.5, size = 2) +
  stat_ellipse(aes(group = cog_intact), linetype = 2) +
  scale_color_manual(values = c("Impaired (<53)" = PAL_PD, "Intact (≥53)" = PAL_CTRL),
                     name = "Cognitive status") +
  labs(title    = "Beta Diversity PCoA — PD Cohort (Intact vs Impaired Cognition)",
       subtitle = paste0("N = ", nsamples(ps_pdaq_rare)),
       caption  = perm_pdaq_label) +
  theme_pubr() +
  theme(plot.caption = element_text(hjust = 0.5, size = 10, face = "italic"))
# plot generated in Script 05

cat("\n--- SECTION 9 COMPLETE ---\n")
cat("Proceed to Script 02 for PDAQ cognitive association analyses.\n")
# =============================================================================
# 8.6 Save Track B: PDAQ association data
# Placed here (after Section 9) because ps_pdaq_pd, ps_pdaq_rare, and all
# PDAQ diversity objects are created in Section 9.
# Also saves cognition cohort diversity plot objects for Script 05.
# =============================================================================
save(ps_pdaq,
     ps_pdaq_pd, ps_pdaq_rare,
     alpha_pdaq_long, alpha_pdaq_stats, alpha_pdaq_plot,
     pcoa_pdaq, permanova_pdaq, perm_pdaq_r2, perm_pdaq_p, perm_pdaq_label,
     pcoa_pdaq_plot,
     file = file.path(RDATA_DIR, "PD_PDAQ_Association_Data.RData"))

cat("\n--- SCRIPT 01 COMPLETE: DATA READY FOR SCRIPT 02 ---\n")
print(ps_pdaq)
