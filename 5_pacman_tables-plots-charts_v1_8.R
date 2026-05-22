# =============================================================================
# Script 05: PACMAN Tables, Plots & Charts v1.8
# =============================================================================
# Changes from v1.7 (v1.8):
#   - Bristol_num added to meta_all mutate block
#   - Bristol stool score row added to Table 1 (before batch row)
# Changes from v1.6:
#   - Section C8 corrected to use alpha_pdaq_plot (was alpha_plot — name
#     collision with Script 01 PD vs Control alpha diversity object)
#   - PDF output removed from save_plot() — PNG only per global directive
#   - alpha_pdaq_plot loaded from PDAQ_DA_Results.RData (Script 02 renamed
#     alpha objects to alpha_pdaq_* to avoid collision)
#   - All analysis plots consolidated here from Scripts 01, 01.5, 02, 03
#   - New RData loads: PD_vs_Control_DA_Results.RData,
#     03_Longitudinal_PDAQ_Results_v2.3.RData
#   - New figures added:
#       Figure C1: Pre-outlier removal PCoA (Script 01)
#       Figure C2: Sequencing depth outlier justification (Script 01)
#       Figure C3: Post-outlier PCoA PD vs Control (Script 01)
#       Figure C4: Full cohort alpha diversity PD vs Control (Script 01)
#       Figure C5: Volcano plot PD vs Control DA (Script 01.5)
#       Figure C6: Beta diversity PCoA — continuous PDAQ (Script 02)
#       Figure C7: Beta diversity PCoA — binary cog_intact (Script 02)
#       Figure C8: Alpha diversity — cog intact vs impaired (Script 02)
#       Figure C9: Sellimonas vs PDAQ slope — 2+ visits (Script 03)
#       Figure C10: Sellimonas vs PDAQ slope — 4+ visit sensitivity (Script 03)
#   - All figures saved as .png (300 dpi)
#   - ggrepel added to library calls for volcano plot
#
# Sections:
#   1.   Setup & data loading
#   2.   PDAQ cognitive bins summary table
#   3.   PDAQ score distribution histogram
#   A.   Participant flowchart
#   B.   Demographic table (Table 1)
#   C1.  Pre-outlier PCoA
#   C2.  Sequencing depth boxplot
#   C3.  Post-outlier PCoA — PD vs Control
#   C4.  Full cohort alpha diversity
#   C5.  Volcano plot — PD vs Control DA
#   C6.  Cognition cohort PCoA — continuous PDAQ
#   C7.  Cognition cohort PCoA — binary cog_intact
#   C8.  Cognition cohort alpha diversity
#   C9.  Sellimonas vs cognitive trajectory (2+ visits)
#   C10. Sellimonas vs cognitive trajectory (4+ visit sensitivity)
# =============================================================================

# =============================================================================
# SECTION 1: SETUP & DATA LOADING
# =============================================================================

.libPaths("/arc/project/st-silkec-1/rsandboxlib")

library(phyloseq)
library(tidyverse)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(grid)
library(gridExtra)

setwd("/arc/project/st-silkec-1/pacman")

RDATA_DIR <- "/arc/project/st-silkec-1/pacman/rdata"
PLOTS_DIR <- "/arc/project/st-silkec-1/pacman/plots"
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = TRUE)
cat("Output directory:", PLOTS_DIR, "\n")

load(file.path(RDATA_DIR, "PD_vs_Control_Cleaned.RData"))
load(file.path(RDATA_DIR, "PD_vs_Control_DA_Results.RData"))
load(file.path(RDATA_DIR, "PD_PDAQ_Association_Data.RData"))
load(file.path(RDATA_DIR, "PDAQ_DA_Results.RData"))
load(file.path(RDATA_DIR, "03_Longitudinal_PDAQ_Results_v2.3.RData"))

# =============================================================================
# PROJECT PALETTE (Okabe-Ito — colourblind-safe across all scripts)
# PD case / Cognitively impaired  →  #D55E00  vermillion
# Control / Cognitively intact    →  #0072B2  blue
# PD-MCI (three-level only)       →  #E69F00  orange
# Sequencing outliers             →  #CC79A7  pink
# =============================================================================
PAL_PD      <- "#D55E00"
PAL_CTRL    <- "#0072B2"
PAL_PDMCI   <- "#E69F00"
PAL_OUTLIER <- "#CC79A7"

# Helper: save both PNG and PDF
save_plot <- function(plot_obj, filename, width, height) {
  ggsave(file.path(PLOTS_DIR, paste0(filename, ".png")),
         plot = plot_obj, width = width, height = height, dpi = 300, bg = "white")
  cat("Saved:", filename, "\n")
}

# =============================================================================
# SECTION 2: PDAQ COGNITIVE BINS SUMMARY TABLE
# =============================================================================

# Ensure PDAQ categories are synced (53 and 31 cutoffs)
meta_plot <- data.frame(sample_data(ps_pdaq)) %>%
  filter(parkinsons == "PD case") %>%
  mutate(
    PDAQ_num = as.numeric(as.character(PDAQ_Total)),
    Cognitive_Status = case_when(
      PDAQ_num >= 53                   ~ "Intact",
      PDAQ_num >= 31 & PDAQ_num < 53   ~ "PD-MCI",
      PDAQ_num < 31                    ~ "PDD",
      TRUE                             ~ NA_character_
    )
  ) %>%
  mutate(Cognitive_Status = factor(Cognitive_Status,
                                   levels = c("Intact", "PD-MCI", "PDD")))

cat("\n--- PDAQ COGNITIVE BINS SUMMARY TABLE ---\n")
pdaq_summary <- meta_plot %>%
  filter(!is.na(Cognitive_Status)) %>%
  group_by(Cognitive_Status) %>%
  summarise(
    Count      = n(),
    Percentage = round((n() / nrow(meta_plot %>%
                                     filter(!is.na(Cognitive_Status)))) * 100, 1),
    Mean_Score = round(mean(PDAQ_num, na.rm = TRUE), 2),
    Min        = min(PDAQ_num, na.rm = TRUE),
    Max        = max(PDAQ_num, na.rm = TRUE)
  )
print(as.data.frame(pdaq_summary))
cat("------------------------------------------\n")

# =============================================================================
# SECTION 3: PDAQ SCORE DISTRIBUTION HISTOGRAM
# =============================================================================

p_hist <- ggplot(meta_plot %>% filter(!is.na(Cognitive_Status)),
                 aes(x = PDAQ_num, fill = Cognitive_Status)) +
  geom_histogram(binwidth = 2, color = "white", alpha = 0.9) +
  scale_fill_manual(values = c("Intact"  = PAL_CTRL,
                               "PD-MCI" = PAL_PDMCI,
                               "PDD"    = PAL_PD)) +
  theme_minimal(base_size = 14) +
  labs(
    title    = "Distribution of PDAQ Scores in PD Cohort",
    subtitle = paste0("Total N = ", sum(pdaq_summary$Count)),
    x        = "PDAQ Total Score",
    y        = "Number of Participants",
    fill     = "Clinical Status"
  ) +
  geom_vline(xintercept = c(31, 53), linetype = "dashed", color = "red")

print(p_hist)
save_plot(p_hist, "PDAQ_Score_Distribution_Histogram", width = 8, height = 5)
cat("Histogram saved.\n")

# =============================================================================
# SECTION A: PARTICIPANT FLOWCHART
# =============================================================================
# Two-track cascade:
#   Track A — Full FoxDen stool cohort (Script 01): 649 → outlier removal → 607
#   Track B — PDMB/FoxInsight filter (Script 00.5): 418 → antibiotic exclusion → 392
#   Convergence: stool cohort ∩ PDMB + complete PDAQ → 362 (355 PD + 7 Control)
#   PD-only branch:
#     355 PD → complete covariates → 351 (Script 02: cross-sectional DA)
#     351    → 2+ PDAQ visits      → 292 (Script 03: longitudinal prediction)
#
# All N values derived from live objects. Fixed upstream values (n_pdmb_raw,
# n_abx_excluded, n_outliers, n_foxden_raw) are hard-coded from confirmed
# console output and annotated with their source.
# =============================================================================

cat("\n--- SECTION A: PARTICIPANT FLOWCHART ---\n")

# --- A1. Derive N values ---

meta_pdaq_full <- data.frame(sample_data(ps_pdaq)) %>%
  mutate(PDAQ_num = as.numeric(as.character(PDAQ_Total)))

# From ps_clean (full FoxDen cohort, post-outlier removal, Script 01)
n_foxden_clean <- nsamples(ps_clean)                                  # 607
meta_clean_all <- data.frame(sample_data(ps_clean))
n_foxden_pd    <- sum(meta_clean_all$parkinsons == "PD case", na.rm = TRUE)  # 403
n_foxden_ctrl  <- sum(meta_clean_all$parkinsons == "Control", na.rm = TRUE)  # 204

# From ps_pdaq (PDAQ cohort, Script 01)
n_pdaq_total   <- nsamples(ps_pdaq)                                   # 362
n_pdaq_pd      <- sum(meta_pdaq_full$parkinsons == "PD case", na.rm = TRUE)  # 355
n_pdaq_ctrl    <- sum(meta_pdaq_full$parkinsons == "Control", na.rm = TRUE)  # 7

# From df_meta_final (Script 02 analysis cohort)
n_analysis_pd  <- nrow(df_meta_final)                                 # 351

# From Script 03 RData (longitudinal cohort, 2+ visits)
n_longitudinal <- nrow(stage2_meta)                                   # 292

# Hard-coded upstream values — confirmed from console output
n_foxden_raw   <- 649L   # Script 01: total stool samples with metadata
n_foxden_pd_raw  <- 435L   # Script 01: PD cases pre-outlier removal
n_foxden_ctrl_raw <- 214L  # Script 01: controls pre-outlier removal
n_outliers     <- 42L    # Script 01: outliers identified (Axis.1 > 0.4)
n_pdmb_raw     <- 418L   # Script 00.5: PDMB participants
n_abx_excluded <- 26L    # Script 00.5: antibiotic exclusions
n_pdmb_post_abx <- n_pdmb_raw - n_abx_excluded                       # 392

# Derived
n_missing_covars <- n_pdaq_pd - n_analysis_pd                        # 4

cat("N values for flowchart:\n")
cat("  FoxDen stool raw:                    ", n_foxden_raw,
    "(", n_foxden_pd_raw, "PD +", n_foxden_ctrl_raw, "Control )\n")
cat("  Outliers removed:                    ", n_outliers, "\n")
cat("  FoxDen cleaned:                      ", n_foxden_clean,
    "(", n_foxden_pd, "PD +", n_foxden_ctrl, "Control )\n")
cat("  PDMB raw:                            ", n_pdmb_raw, "\n")
cat("  Antibiotic exclusions:               ", n_abx_excluded, "\n")
cat("  PDMB post-abx:                       ", n_pdmb_post_abx, "\n")
cat("  PDAQ cohort (intersection):          ", n_pdaq_total,
    "(", n_pdaq_pd, "PD +", n_pdaq_ctrl, "Control )\n")
cat("  PD cases, missing covariates:        ", n_missing_covars, "\n")
cat("  Script 02 analysis cohort (PD):      ", n_analysis_pd, "\n")
cat("  Script 03 longitudinal cohort (PD):  ", n_longitudinal, "\n")

# --- A2. Build flowchart ---
# Layout: two-column approach.
#   Left column  (x = 0.50): main cascade boxes
#   Right column (x = 1.60): exclusion side-boxes
#   Analysis annotation boxes sit to the right of the PD-branch boxes
#
# Colour scheme (Okabe-Ito):
#   Teal   (#C3E6D8 / #009E73) — raw data entry boxes
#   Blue   (#CCE5F6 / #0072B2) — retained cohort steps
#   Yellow (#FDF6C3 / #E69F00) — exclusion side-boxes
#   Orange (#FAD7C3 / #D55E00) — final analysis cohorts
#   Purple (#E8D5F5 / #7B2D8B) — analysis annotation banners

make_box <- function(id, x, y, w, h, label,
                     fill = "#CCE5F6", colour = "#0072B2") {
  data.frame(id = id, x = x, y = y, w = w, h = h,
             label = label, fill = fill, colour = colour,
             stringsAsFactors = FALSE)
}

fmt_n <- function(n) ifelse(is.na(n), "FILL IN", as.character(n))

# Y positions (top to bottom): 9.2, 8.0, 6.8, 5.6, 4.4, 3.2, 2.0
# Box height: 0.70 throughout for uniformity
BH <- 0.70   # box height
XL <- 1.42   # main column x centre
XR <- 2.52   # exclusion column x centre
WL <- 0.90   # main box width
WR <- 0.62   # exclusion box width
WA <- 0.72   # analysis annotation box width
XA <- 0.42   # analysis annotation x centre (left side)

boxes <- bind_rows(
  
  # ── Track A: FoxDen stool cohort ──────────────────────────────────────────
  make_box(1, XL, 9.2, WL, BH,
           paste0("FoxDen stool cohort\nN = ", n_foxden_raw,
                  "  (", n_foxden_pd_raw, " PD + ", n_foxden_ctrl_raw, " Control)"),
           fill = "#C3E6D8", colour = "#009E73"),
  
  make_box(2, XL, 8.0, WL, BH,
           paste0("After outlier removal\n(Bray-Curtis PCoA, Axis.1 > 0.4)\nN = ",
                  n_foxden_clean,
                  "  (", n_foxden_pd, " PD + ", n_foxden_ctrl, " Control)"),
           fill = "#CCE5F6", colour = "#0072B2"),
  
  # ── Convergence: PDAQ cohort ───────────────────────────────────────────────
  make_box(3, XL, 6.4, WL, BH,
           paste0("PDAQ cohort\n(stool \u2229 PDMB, complete PDAQ-15)\nN = ",
                  n_pdaq_total,
                  "  (", n_pdaq_pd, " PD + ", n_pdaq_ctrl, " Control)"),
           fill = "#CCE5F6", colour = "#0072B2"),
  
  # ── PD-only branch ────────────────────────────────────────────────────────
  make_box(4, XL, 5.0, WL, BH,
           paste0("PD cases\n(complete covariates)\nN = ", n_analysis_pd),
           fill = "#FAD7C3", colour = "#D55E00"),
  
  make_box(5, XL, 3.6, WL, BH,
           paste0("PD cases, 2+ PDAQ visits\n(estimable cognitive slope)\nN = ",
                  n_longitudinal),
           fill = "#FAD7C3", colour = "#D55E00"),
  
  # ── Track B: PDMB filter (feeds into convergence box) ────────────────────
  make_box(6, XR, 9.2, WR, BH,
           paste0("PDMB participants\nN = ", n_pdmb_raw),
           fill = "#C3E6D8", colour = "#009E73"),
  
  make_box(7, XR, 8.0, WR, BH,
           paste0("After abx exclusion\nN = ", n_pdmb_post_abx),
           fill = "#CCE5F6", colour = "#0072B2"),
  
  # ── Exclusion side-boxes ──────────────────────────────────────────────────
  make_box(8, XR, 6.4, WR, BH,
           paste0("7 controls\nnot used analytically"),
           fill = "#FDF6C3", colour = "#E69F00"),
  
  make_box(9, XR, 5.0, WR, BH,
           paste0("Excluded:\nmissing covariates\nN = ", n_missing_covars),
           fill = "#FDF6C3", colour = "#E69F00"),
  
  make_box(10, XR, 3.6, WR, BH,
           paste0("Excluded:\nsingle-visit only\nN = ", n_analysis_pd - n_longitudinal),
           fill = "#FDF6C3", colour = "#E69F00"),
  
  # ── Analysis annotation boxes ─────────────────────────────────────────────
  make_box(11, XA, 8.0, WA, BH,
           "Scripts 01 & 01.5\nFull cohort diversity\n& PD vs Control DA",
           fill = "#E8D5F5", colour = "#7B2D8B"),
  
  make_box(12, XA, 5.0, WA, BH,
           "Script 02\nCross-sectional\ncognitive DA",
           fill = "#E8D5F5", colour = "#7B2D8B"),
  
  make_box(13, XA, 3.6, WA, BH,
           "Script 03\nLongitudinal\ncognitive prediction",
           fill = "#E8D5F5", colour = "#7B2D8B")
)

# Downward arrows along main column
arrows_main <- data.frame(
  x    = XL,
  xend = XL,
  y    = c(9.2 - BH/2, 8.0 - BH/2, 6.4 - BH/2, 5.0 - BH/2),
  yend = c(8.0 + BH/2, 6.4 + BH/2, 5.0 + BH/2, 3.6 + BH/2)
)

# Downward arrow along Track B
arrows_trackB <- data.frame(
  x    = XR,
  xend = XR,
  y    = 9.2 - BH/2,
  yend = 8.0 + BH/2
)

# Converging arrow: Track B (post-abx) feeds into PDAQ cohort box
# L-shaped: goes down from box 7 then left to box 3
conv_path <- data.frame(
  x    = c(XR,  XR,  XL + WL/2),
  y    = c(8.0 - BH/2, 6.4, 6.4),
  group = 1
)

# Dashed exclusion arrows (horizontal, right from main column)
excl_df <- data.frame(
  x    = XL + WL/2,
  xend = XR - WR/2,
  y    = c(6.4, 5.0, 3.6),
  yend = c(6.4, 5.0, 3.6)
)

# Analysis annotation arrows (horizontal, right from analysis cohort boxes)
annot_df <- data.frame(
  x    = XL - WL/2,
  xend = XA + WA/2,
  y    = c(8.0, 5.0, 3.6),
  yend = c(8.0, 5.0, 3.6)
)

p_flowchart <- ggplot() +
  # Boxes
  geom_rect(data = boxes,
            aes(xmin = x - w/2, xmax = x + w/2,
                ymin = y - h/2, ymax = y + h/2,
                fill = fill, colour = colour),
            linewidth = 0.8) +
  scale_fill_identity() +
  scale_colour_identity() +
  # Labels
  geom_text(data = boxes,
            aes(x = x, y = y, label = label),
            size = 3.8, lineheight = 1.2, fontface = "plain") +
  # Main column downward arrows
  geom_segment(data = arrows_main,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               colour = "#0072B2", linewidth = 0.7) +
  # Track B downward arrow
  geom_segment(data = arrows_trackB,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               colour = "#009E73", linewidth = 0.7) +
  # Converging L-path from Track B into PDAQ cohort box
  geom_path(data = conv_path,
            aes(x = x, y = y, group = group),
            arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
            colour = "#009E73", linewidth = 0.7) +
  # Dashed exclusion arrows
  geom_segment(data = excl_df,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.15, "cm"), type = "open"),
               colour = "#E69F00", linewidth = 0.55, linetype = "dashed") +
  # Analysis annotation arrows
  geom_segment(data = annot_df,
               aes(x = x, xend = xend, y = y, yend = yend),
               colour = "#7B2D8B", linewidth = 0.55, linetype = "dotted") +
  scale_x_continuous(limits = c(0.0, 3.20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(3.1, 9.7), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(12, 12, 12, 12))

print(p_flowchart)
save_plot(p_flowchart, "Figure_Participant_Flowchart", width = 10, height = 8)
cat("Flowchart saved.\n")

# =============================================================================
# SECTION B: DEMOGRAPHIC TABLE (TABLE 1)
# =============================================================================
# Three columns: Full PDAQ Cohort | PD Cases | Controls
# Continuous variables: median [IQR]
# Categorical variables: n (%)
# Between-group p-values (PD vs Control):
#   Wilcoxon rank-sum for continuous; chi-squared or Fisher's exact for categorical
# Disease duration and azilect are PD-only (Controls show —)
# =============================================================================

cat("\n--- SECTION B: DEMOGRAPHIC TABLE (TABLE 1) ---\n")

# --- B1. Prepare data ---

freq_map_b    <- c("never"=0, "< 3x/week"=1.5, "3-6x/week"=4.5,
                   "daily"=7, "2-3x/day"=17.5, "> 3x/day"=25)
target_cols_b <- c("bl_bristol_stool_type_1_freq", "bl_bristol_stool_type_2_freq",
                   "bl_bristol_stool_type_3_freq", "bl_bristol_stool_type_4_freq",
                   "bl_bristol_stool_type_5_freq", "bl_bristol_stool_type_6_freq",
                   "bl_bristol_stool_type_7_freq")

meta_raw <- data.frame(sample_data(ps_pdaq))
bmat_b   <- sapply(target_cols_b, function(col) freq_map_b[as.character(meta_raw[[col]])])
w_b      <- as.numeric(bmat_b %*% 1:7)
tot_b    <- rowSums(bmat_b, na.rm = TRUE)

meta_all <- meta_raw %>%
  mutate(
    PDAQ_num         = as.numeric(as.character(PDAQ_Total)),
    current_age      = as.numeric(as.character(current_age)),
    BMI_num          = as.numeric(as.character(BMI)),
    sex_clean        = as.character(sex),
    batch_clean      = as.character(main_or_pilot_study),
    parkinsons_clean = as.character(parkinsons),
    disease_duration = suppressWarnings(
      as.numeric(as.character(current_age)) -
        as.numeric(as.character(age_diagnosed))
    ),
    disease_duration = ifelse(disease_duration < 0, NA, disease_duration),
    azilect_bin      = ifelse(
      !is.na(`pd_specific_medications.azilect`) &
        !grepl("Not", `pd_specific_medications.azilect`), 1L, 0L),
    cog_intact       = case_when(
      PDAQ_num >= 53 ~ "Intact (\u226553)",
      PDAQ_num <  53 ~ "Impaired (<53)",
      TRUE           ~ NA_character_
    ),
    Bristol_num      = ifelse(tot_b > 0, w_b / tot_b, NA_real_)
  )

meta_pd   <- meta_all %>% filter(parkinsons_clean == "PD case")
meta_ctrl <- meta_all %>% filter(parkinsons_clean == "Control")

# --- B2. Helper functions ---

med_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("\u2014")
  sprintf("%.1f [%.1f\u2013%.1f]",
          median(x), quantile(x, 0.25), quantile(x, 0.75))
}

pct_n <- function(x, level) {
  n_total <- sum(!is.na(x))
  n_level <- sum(x == level, na.rm = TRUE)
  if (n_total == 0) return("\u2014")
  sprintf("%d (%.1f%%)", n_level, 100 * n_level / n_total)
}

wilcox_p <- function(x, g) {
  df <- data.frame(x = x, g = g) %>% filter(!is.na(x), !is.na(g))
  if (length(unique(df$g)) < 2) return(NA_real_)
  tryCatch(wilcox.test(x ~ g, data = df, exact = FALSE)$p.value,
           error = function(e) NA_real_)
}

chisq_p <- function(x, g) {
  tbl <- table(x, g)
  if (any(dim(tbl) < 2)) return(NA_real_)
  expected <- suppressWarnings(chisq.test(tbl)$expected)
  tryCatch(
    if (any(expected < 5)) fisher.test(tbl, simulate.p.value = TRUE)$p.value
    else chisq.test(tbl)$p.value,
    error = function(e) NA_real_
  )
}

fmt_p <- function(p) {
  if (is.na(p))    return("\u2014")
  if (p < 0.001)   return("<0.001")
  sprintf("%.3f", p)
}

# --- B3. Assemble rows ---

build_row <- function(variable, all_val, pd_val, ctrl_val, p_val = NA) {
  data.frame(Variable      = variable,
             `Full Cohort` = all_val,
             `PD Cases`    = pd_val,
             `Controls`    = ctrl_val,
             `p-value`     = fmt_p(p_val),
             check.names = FALSE, stringsAsFactors = FALSE)
}

table1 <- bind_rows(
  
  build_row("N",
            as.character(nrow(meta_all)),
            as.character(nrow(meta_pd)),
            as.character(nrow(meta_ctrl))),
  
  build_row("Age, years \u2014 median [IQR]",
            med_iqr(meta_all$current_age),
            med_iqr(meta_pd$current_age),
            med_iqr(meta_ctrl$current_age),
            wilcox_p(meta_all$current_age, meta_all$parkinsons_clean)),
  
  build_row("Female sex \u2014 n (%)",
            pct_n(meta_all$sex_clean,  "Female"),
            pct_n(meta_pd$sex_clean,   "Female"),
            pct_n(meta_ctrl$sex_clean, "Female"),
            chisq_p(meta_all$sex_clean, meta_all$parkinsons_clean)),
  
  build_row("BMI, kg/m\u00b2 \u2014 median [IQR]",
            med_iqr(meta_all$BMI_num),
            med_iqr(meta_pd$BMI_num),
            med_iqr(meta_ctrl$BMI_num),
            wilcox_p(meta_all$BMI_num, meta_all$parkinsons_clean)),
  
  build_row("Disease duration, years \u2014 median [IQR] \u00b9",
            med_iqr(meta_all$disease_duration),
            med_iqr(meta_pd$disease_duration),
            "\u2014",
            NA_real_),
  
  build_row("PDAQ Total score \u2014 median [IQR]",
            med_iqr(meta_all$PDAQ_num),
            med_iqr(meta_pd$PDAQ_num),
            med_iqr(meta_ctrl$PDAQ_num),
            wilcox_p(meta_all$PDAQ_num, meta_all$parkinsons_clean)),
  
  build_row("Cognitively intact (PDAQ \u226553) \u2014 n (%)",
            pct_n(meta_all$cog_intact,  "Intact (\u226553)"),
            pct_n(meta_pd$cog_intact,   "Intact (\u226553)"),
            pct_n(meta_ctrl$cog_intact, "Intact (\u226553)"),
            chisq_p(meta_all$cog_intact, meta_all$parkinsons_clean)),
  
  build_row("Azilect (rasagiline) use \u2014 n (%) \u00b9",
            pct_n(meta_all$azilect_bin, 1),
            pct_n(meta_pd$azilect_bin,  1),
            "\u2014",
            NA_real_),
  
  build_row("Bristol stool score \u2014 median [IQR]",
            med_iqr(meta_all$Bristol_num),
            med_iqr(meta_pd$Bristol_num),
            med_iqr(meta_ctrl$Bristol_num),
            wilcox_p(meta_all$Bristol_num, meta_all$parkinsons_clean)),
  
  build_row("Main study cohort \u2014 n (%)",
            pct_n(meta_all$batch_clean,  "Main"),
            pct_n(meta_pd$batch_clean,   "Main"),
            pct_n(meta_ctrl$batch_clean, "Main"),
            chisq_p(meta_all$batch_clean, meta_all$parkinsons_clean))
)

# --- B4. Print to console ---

cat("\nTable 1. Demographic and Clinical Characteristics\n")
cat(strrep("=", 80), "\n")
print(format(table1, justify = "left"), row.names = FALSE)
cat(strrep("-", 80), "\n")
cat("\u00b9 PD cases only. Disease duration = age at collection \u2212 age at diagnosis.\n")
cat("  Azilect prevalence denominator excludes missing values.\n")
cat("  p-values: Wilcoxon rank-sum (continuous); \u03c7\u00b2 or Fisher's exact (categorical).\n")
cat("  \u2014 indicates variable not applicable for that group.\n\n")

# --- B5. Export CSV ---

write.csv(table1,
          file.path(PLOTS_DIR, "Table1_Demographics.csv"),
          row.names = FALSE, na = "\u2014")
cat("Table 1 exported to", file.path(PLOTS_DIR, "Table1_Demographics.csv"), "\n")

# --- B6. Optional: export as PNG via tableGrob ---
# Uncomment to activate. Requires gridExtra (loaded above).

# tbl_grob <- tableGrob(
#   table1, rows = NULL,
#   theme = ttheme_minimal(
#     base_size = 9,
#     core    = list(fg_params = list(hjust = 0, x = 0.05)),
#     colhead = list(fg_params = list(fontface = "bold"))
#   )
# )
# ggsave(file.path(PLOTS_DIR, "Table1_Demographics.png"),
#        plot = grid.arrange(tbl_grob),
#        width = 12, height = 5, dpi = 300, bg = "white")
# cat("Table 1 PNG saved to", file.path(PLOTS_DIR, "Table1_Demographics.png"), "\n")


# =============================================================================
# Phylum-Level Composition Bar Chart: PD vs Healthy Control
# =============================================================================
# Prerequisite: ps_clean must be in the environment
#   load("PD_vs_Control_Cleaned.RData")
#
# Output: Figure_Phylum_Composition_PD_vs_Control.png
# =============================================================================

# --- 1. Agglomerate to phylum and compute per-sample relative abundance ---

ps_phylum <- tax_glom(ps_clean, taxrank = "Phylum")
ps_rel    <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

phylum_long <- psmelt(ps_rel) %>%
  mutate(
    Phylum = gsub("^p__", "", Phylum),
    Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unclassified", Phylum),
    Group  = case_when(
      parkinsons == "PD case" ~ "Parkinson's Disease",
      parkinsons == "Control" ~ "Healthy Control",
      TRUE                    ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

# --- 2. Identify phyla >= 1% mean abundance (across both groups combined) ---

phyla_keep <- phylum_long %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(mean_abund >= 0.01) %>%
  pull(Phylum)

cat("Phyla retained (>= 1%):", paste(phyla_keep, collapse = ", "), "\n")

# --- 3. Collapse, summarise, and order ---

phylum_plot_df <- phylum_long %>%
  mutate(Phylum_grouped = ifelse(Phylum %in% phyla_keep, Phylum, "Other (<1%)")) %>%
  group_by(Sample, Group, Phylum_grouped) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, Phylum_grouped) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# Order: largest-to-smallest overall, Other always last
phylum_order <- phylum_plot_df %>%
  group_by(Phylum_grouped) %>%
  summarise(total = sum(Mean_Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  filter(Phylum_grouped != "Other (<1%)") %>%
  pull(Phylum_grouped) %>%
  c("Other (<1%)")

phylum_plot_df <- phylum_plot_df %>%
  mutate(
    Phylum_grouped = factor(Phylum_grouped, levels = rev(phylum_order)),
    Group          = factor(Group, levels = c("Parkinson's Disease", "Healthy Control"))
  )

# --- 4. Colour palette (Paul Tol muted, colourblind-friendly; grey for Other) ---

n_phyla   <- length(phylum_order) - 1
tol_muted <- c("#332288", "#88CCEE", "#44AA99", "#117733",
               "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
pal        <- c(tol_muted[seq_len(n_phyla)], "#AAAAAA")
names(pal) <- phylum_order

# --- 5. Plot ---

p_phylum <- ggplot(phylum_plot_df,
                   aes(x = Group, y = Mean_Abundance, fill = Phylum_grouped)) +
  geom_bar(stat = "identity", position = "stack",
           colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = pal, breaks = phylum_order, name = "Phylum") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Phylum-Level Gut Microbiome Composition",
    x     = NULL,
    y     = "Mean Relative Abundance"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x        = element_text(size = 12, face = "bold"),
    legend.position    = "right",
    legend.key.size    = unit(0.5, "cm"),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold")
  )

print(p_phylum)
save_plot(p_phylum, "Figure_Phylum_Composition_PD_vs_Control", width = 7, height = 6)
cat("Saved: Figure_Phylum_Composition_PD_vs_Control\n")

# --- 6. Summary table ---

phylum_long %>%
  mutate(Phylum_grouped = ifelse(Phylum %in% phyla_keep, Phylum, "Other (<1%)")) %>%
  group_by(Sample, Group, Phylum_grouped) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, Phylum_grouped) %>%
  summarise(Mean_Pct = round(mean(Abundance, na.rm = TRUE) * 100, 2), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Group, values_from = Mean_Pct) %>%
  arrange(match(Phylum_grouped, phylum_order)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# =============================================================================
# Genus-Level Composition Bar Chart: PD vs Healthy Control (Top 15 Genera)
# =============================================================================
# Prerequisite: ps_clean must be in the environment
#   load("PD_vs_Control_Cleaned.RData")
#
# Approach: top 15 genera by mean relative abundance across both groups
# combined; all remaining genera collapsed into "Other".
#
# Output: Figure_Genus_Composition_PD_vs_Control.png
# =============================================================================

# --- 1. Agglomerate to genus and compute per-sample relative abundance ---

ps_genus_comp <- tax_glom(ps_clean, taxrank = "Genus")
ps_rel        <- transform_sample_counts(ps_genus_comp, function(x) x / sum(x))

genus_long <- psmelt(ps_rel) %>%
  mutate(
    Genus = gsub("^g__", "", Genus),
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus),
    Group = case_when(
      parkinsons == "PD case" ~ "Parkinson's Disease",
      parkinsons == "Control" ~ "Healthy Control",
      TRUE                    ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

# --- 2. Identify top 15 genera by mean abundance across both groups ---

top15_genera <- genus_long %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

cat("Top 15 genera retained:\n")
print(top15_genera)

# --- 3. Collapse, summarise, and order ---

genus_plot_df <- genus_long %>%
  mutate(Genus_grouped = ifelse(Genus %in% top15_genera, Genus, "Other")) %>%
  group_by(Sample, Group, Genus_grouped) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, Genus_grouped) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# Order: largest-to-smallest overall, Other always last
genus_order <- genus_plot_df %>%
  group_by(Genus_grouped) %>%
  summarise(total = sum(Mean_Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  filter(Genus_grouped != "Other") %>%
  pull(Genus_grouped) %>%
  c("Other")

genus_plot_df <- genus_plot_df %>%
  mutate(
    Genus_grouped = factor(Genus_grouped, levels = rev(genus_order)),
    Group         = factor(Group, levels = c("Parkinson's Disease", "Healthy Control"))
  )

# --- 4. Colour palette ---
# 15 named genera + Other = 16 colours.
# Uses Paul Tol's high-contrast + muted schemes combined; grey for Other.

pal_16 <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
  "#DDCC77", "#CC6677", "#882255", "#AA4499", "#6699CC",
  "#661100", "#D9A755", "#4477AA", "#66CCEE", "#228833",
  "#AAAAAA"  # Other
)
names(pal_16) <- genus_order

# --- 5. Plot ---

p_genus <- ggplot(genus_plot_df,
                  aes(x = Group, y = Mean_Abundance, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack",
           colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = pal_16, breaks = genus_order,
                    name = "Genus") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Genus-Level Gut Microbiome Composition (Top 15)",
    x     = NULL,
    y     = "Mean Relative Abundance"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x        = element_text(size = 12, face = "bold"),
    legend.position    = "right",
    legend.key.size    = unit(0.45, "cm"),
    legend.text        = element_text(face = "italic", size = 9),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold")
  )

print(p_genus)
save_plot(p_genus, "Figure_Genus_Composition_PD_vs_Control", width = 8, height = 6)
cat("Saved: Figure_Genus_Composition_PD_vs_Control\n")

# --- 6. Summary table ---

genus_long %>%
  mutate(Genus_grouped = ifelse(Genus %in% top15_genera, Genus, "Other")) %>%
  group_by(Sample, Group, Genus_grouped) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, Genus_grouped) %>%
  summarise(Mean_Pct = round(mean(Abundance, na.rm = TRUE) * 100, 2),
            .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Group, values_from = Mean_Pct) %>%
  arrange(match(Genus_grouped, genus_order)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# =============================================================================
# SECTION C1: PRE-OUTLIER REMOVAL PCoA
# =============================================================================
cat("\n--- SECTION C1: Pre-outlier PCoA ---\n")
print(pcoa_initial_plot)
save_plot(pcoa_initial_plot, "Figure_C1_PCoA_PreOutlierRemoval", width = 7, height = 5)

# =============================================================================
# SECTION C2: SEQUENCING DEPTH OUTLIER JUSTIFICATION
# =============================================================================
cat("\n--- SECTION C2: Sequencing depth ---\n")
print(depth_plot)
save_plot(depth_plot, "Figure_C2_SequencingDepth_Outliers", width = 6, height = 5)

# =============================================================================
# SECTION C3: POST-OUTLIER PCoA — PD vs CONTROL
# =============================================================================
cat("\n--- SECTION C3: Post-outlier PCoA PD vs Control ---\n")
print(pcoa_final_plot)
save_plot(pcoa_final_plot, "Figure_C3_PCoA_PDvsControl_Cleaned", width = 7, height = 5)

# =============================================================================
# SECTION C4: FULL COHORT ALPHA DIVERSITY — PD vs CONTROL
# =============================================================================
cat("\n--- SECTION C4: Full cohort alpha diversity ---\n")
print(alpha_plot)
save_plot(alpha_plot, "Figure_C4_Alpha_PDvsControl", width = 9, height = 5)

# =============================================================================
# SECTION C5: VOLCANO PLOT — PD vs CONTROL DA
# =============================================================================
cat("\n--- SECTION C5: Volcano plot PD vs Control ---\n")

# Rebuild volcano with updated labelling: 2-of-3 and 3-of-3 consensus only
# (no single-method labels per thesis directive)
volcano_df_plot <- volcano_df %>%
  mutate(
    # Reclassify single-method hits as "Not significant" for display —
    # only 2-of-3 and 3-of-3 consensus hits are highlighted per thesis directive
    consensus = case_when(
      taxon %in% consensus_3$taxon  ~ "3-of-3 consensus",
      taxon %in% two_of_three$taxon ~ "2-of-3 consensus",
      TRUE                          ~ "Not significant"
    ),
    consensus = factor(consensus, levels = c(
      "3-of-3 consensus",
      "2-of-3 consensus",
      "Not significant"
    )),
    lab = case_when(
      taxon %in% consensus_3$taxon  ~ label,
      taxon %in% two_of_three$taxon ~ label,
      TRUE                          ~ NA_character_
    )
  )

p_volcano_05 <- ggplot(volcano_df_plot,
                       aes(x = coef, y = neg_log10_q,
                           colour = consensus, size = consensus)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text",
           x     = max(volcano_df_plot$coef, na.rm = TRUE),
           y     = -log10(0.05) + 0.15,
           label = "FDR = 0.05", hjust = 1,
           size  = 3.2, colour = "grey40") +
  geom_point(alpha = 0.75) +
  geom_text_repel(
    aes(label = lab),
    size           = 3.2,
    fontface       = "italic",
    max.overlaps   = 20,
    box.padding    = 0.5,
    point.padding  = 0.3,
    show.legend    = FALSE,
    segment.colour = "grey50",
    segment.size   = 0.3
  ) +
  scale_colour_manual(
    values = c("3-of-3 consensus" = PAL_PD,
               "2-of-3 consensus" = PAL_PDMCI,
               "Not significant"  = "grey75"),
    name = NULL
  ) +
  scale_size_manual(
    values = c("3-of-3 consensus" = 4,
               "2-of-3 consensus" = 3.5,
               "Not significant"  = 1.2),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(min(volcano_df_plot$coef, na.rm = TRUE) - 0.1,
               max(volcano_df_plot$coef, na.rm = TRUE) + 0.1)
  ) +
  labs(
    title    = "Differential Abundance — PD vs Healthy Control",
    subtitle = paste0("N = ", nrow(consensus_3), " at 3-of-3, ",
                      nrow(two_of_three), " at 2-of-3 consensus"),
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

print(p_volcano_05)
save_plot(p_volcano_05, "Figure_C5_Volcano_PDvsControl", width = 8, height = 6)

# =============================================================================
# SECTION C6: COGNITION COHORT PCoA — CONTINUOUS PDAQ
# =============================================================================
cat("\n--- SECTION C6: Cognition cohort PCoA (continuous PDAQ) ---\n")
print(pcoa_cont_plot)
save_plot(pcoa_cont_plot, "Figure_C6_PCoA_ContinuousPDAQ", width = 7, height = 5)

# =============================================================================
# SECTION C7: COGNITION COHORT PCoA — BINARY COG_INTACT
# =============================================================================
cat("\n--- SECTION C7: Cognition cohort PCoA (binary cog status) ---\n")
print(pcoa_bin_plot)
save_plot(pcoa_bin_plot, "Figure_C7_PCoA_BinaryCogStatus", width = 7, height = 5)

# =============================================================================
# SECTION C8: COGNITION COHORT ALPHA DIVERSITY
# =============================================================================
cat("\n--- SECTION C8: Cognition cohort alpha diversity ---\n")
print(alpha_pdaq_plot)
save_plot(alpha_pdaq_plot, "Figure_C8_Alpha_CogIntactVsImpaired", width = 9, height = 5)

# =============================================================================
# SECTION C9: SELLIMONAS VS COGNITIVE TRAJECTORY — 2+ VISITS (N=292)
# =============================================================================
cat("\n--- SECTION C9: Sellimonas vs cognitive trajectory (2+ visits) ---\n")

if (!is.null(plot_df) && nrow(plot_df) > 0) {
  p_sell_primary <- ggplot(plot_df, aes(x = sellimonas_log, y = pdaq_slope)) +
    geom_point(alpha = 0.5, size = 2, colour = PAL_CTRL) +
    geom_smooth(method = "lm", se = TRUE,
                colour = "#005A91", linewidth = 0.8) +
    stat_cor(method = "pearson",
             label.x.npc = "left", label.y.npc = "top") +
    labs(
      x        = "Sellimonas abundance (TSS + log normalised)",
      y        = "Individual PDAQ slope (points / year)",
      title    = "Baseline Sellimonas vs Cognitive Trajectory",
      subtitle = paste0("Primary analysis  |  N = ", nrow(plot_df),
                        "  (2+ PDAQ visits)"),
      caption  = "Negative slope = faster cognitive decline. Shaded band = 95% CI."
    ) +
    theme_pubr() +
    theme(plot.caption = element_text(hjust = 0, size = 9, face = "italic"))
  print(p_sell_primary)
  save_plot(p_sell_primary, "Figure_C9_Sellimonas_Slope_Primary", width = 6, height = 5)
} else {
  cat("WARNING: plot_df not found — re-run Script 03 v2.3 to generate\n")
}

# =============================================================================
# SECTION C10: SELLIMONAS VS COGNITIVE TRAJECTORY — 4+ VISIT SENSITIVITY
# =============================================================================
cat("\n--- SECTION C10: Sellimonas vs cognitive trajectory (4+ visit sensitivity) ---\n")

if (!is.null(plot_df_sens) && nrow(plot_df_sens) > 0) {
  p_sell_sens <- ggplot(plot_df_sens, aes(x = sellimonas_log, y = pdaq_slope)) +
    geom_point(alpha = 0.5, size = 2, colour = PAL_PD) +
    geom_smooth(method = "lm", se = TRUE,
                colour = "#993C1D", linewidth = 0.8) +
    stat_cor(method = "pearson",
             label.x.npc = "left", label.y.npc = "top") +
    labs(
      x        = "Sellimonas abundance (TSS + log normalised)",
      y        = "Individual PDAQ slope (points / year)",
      title    = "Baseline Sellimonas vs Cognitive Trajectory",
      subtitle = paste0("Sensitivity analysis  |  N = ", nrow(plot_df_sens),
                        "  (4+ PDAQ visits)"),
      caption  = "Negative slope = faster cognitive decline. Shaded band = 95% CI."
    ) +
    theme_pubr() +
    theme(plot.caption = element_text(hjust = 0, size = 9, face = "italic"))
  print(p_sell_sens)
  save_plot(p_sell_sens, "Figure_C10_Sellimonas_Slope_Sensitivity", width = 6, height = 5)
} else {
  cat("WARNING: plot_df_sens not found — re-run Script 03 v2.3 to generate\n")
}

# =============================================================================
# SCRIPT 05 COMPLETE
# =============================================================================

cat("\n--- SCRIPT 05 v1.8 COMPLETE ---\n")
cat("All outputs saved to:", PLOTS_DIR, "\n")
cat("\nFigures generated:\n")
cat("  Section 3:  PDAQ_Score_Distribution_Histogram\n")
cat("  Section A:  Figure_Participant_Flowchart\n")
cat("  Section B:  Table1_Demographics.csv\n")
cat("  Section C1: Figure_C1_PCoA_PreOutlierRemoval\n")
cat("  Section C2: Figure_C2_SequencingDepth_Outliers\n")
cat("  Section C3: Figure_C3_PCoA_PDvsControl_Cleaned\n")
cat("  Section C4: Figure_C4_Alpha_PDvsControl\n")
cat("  Section C5: Figure_C5_Volcano_PDvsControl\n")
cat("  Section C6: Figure_C6_PCoA_ContinuousPDAQ\n")
cat("  Section C7: Figure_C7_PCoA_BinaryCogStatus\n")
cat("  Section C8: Figure_C8_Alpha_CogIntactVsImpaired\n")
cat("  Section C9: Figure_C9_Sellimonas_Slope_Primary\n")
cat("  Section C10: Figure_C10_Sellimonas_Slope_Sensitivity\n")
cat("  Composition: Figure_Phylum_Composition_PD_vs_Control\n")
cat("               Figure_Genus_Composition_PD_vs_Control\n")
cat("\nAll figures saved as .png (300 dpi) and .pdf (vector).\n")