# ============================================================
# BACKGROUND DETERMINATION: POWER ANALYSIS
# Site: [INSERT SITE NAME]
# Analyst: [INSERT NAME]
# Date: [INSERT DATE]
# R version: 4.3+
# ============================================================

# ── 1. INSTALL AND LOAD PACKAGES ────────────────────────────

# Run install lines once, then comment them out
# install.packages("pwr")
# install.packages("effsize")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("writexl")
# install.packages("scales")

library(pwr)
library(effsize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(writexl)
library(scales)

# ── 2. LOAD YOUR DATA ───────────────────────────────────────
# Expected CSV structure — one row per sample:
#
# Sample_ID | Population | Arsenic | Barium | Cadmium | Chromium |
# Copper | Lead | Nickel | Selenium | Silver | Zinc
#
# Population column must contain exactly "Background" or "Site"
# Non-detects must be entered as NA (not zero, not half-DL)
# Detected values must be positive numeric concentrations in mg/kg

data <- read.csv("data/data.csv", stringsAsFactors = FALSE)

# Quick check — confirm population labels and sample counts
cat("Sample counts by population:\n")
print(table(data$Population))

cat("\nColumn names in dataset:\n")
print(names(data))

# ── 3. DEFINE METALS ────────────────────────────────────────
# These must match your column names exactly

metals <- c("Arsenic", "Barium", "Cadmium", "Chromium",
            "Copper",  "Lead",   "Nickel",  "Selenium",
            "Silver",  "Zinc")

# Minimum detections required in each group for a meaningful
# power calculation. Metals below this in either group return NA.
min_detections <- 3

# ── 4. SPLIT INTO BACKGROUND AND SITE ───────────────────────

bg_data   <- data %>% filter(Population == "Background")
site_data <- data %>% filter(Population == "Site")

cat("\nBackground samples:", nrow(bg_data))
cat("\nSite samples:", nrow(site_data), "\n")

# ── 5. CALCULATE DETECTED N DIRECTLY FROM RAW DATA ──────────
# sapply loops over each metal and counts non-NA values per group.
# NA = non-detect. Positive numeric = detected.

n_bg_detected   <- sapply(metals, function(m) sum(!is.na(bg_data[[m]])))
n_site_detected <- sapply(metals, function(m) sum(!is.na(site_data[[m]])))
n_bg_total      <- nrow(bg_data)
n_site_total    <- nrow(site_data)

detected_counts <- data.frame(
  Metal             = metals,
  N_BG_Total        = n_bg_total,
  N_Site_Total      = n_site_total,
  N_BG_Detected     = n_bg_detected,
  N_Site_Detected   = n_site_detected,
  Pct_BG_Detected   = round(n_bg_detected   / n_bg_total   * 100, 1),
  Pct_Site_Detected = round(n_site_detected / n_site_total * 100, 1),
  Pct_Overall       = round((n_bg_detected + n_site_detected) /
                              (n_bg_total + n_site_total) * 100, 1),
  row.names         = NULL,
  stringsAsFactors  = FALSE
)

cat("\n── Detectability summary ──────────────────────────────\n")
print(detected_counts)

# Ensure outputs directory exists
if (!dir.exists("outputs")) dir.create("outputs")

write_xlsx(detected_counts, "outputs/detectability_summary.xlsx")
cat("Detectability summary exported to outputs/detectability_summary.xlsx\n")

# ── 6. COMPUTE COHEN'S D FROM RAW DETECTED VALUES ───────────
# Computed on log-transformed detected values only,
# consistent with the normality testing approach in this SOP.
# Returns NA if either group has fewer than min_detections.

compute_d <- function(metal) {

  bg_vals   <- bg_data[[metal]]
  site_vals <- site_data[[metal]]

  bg_clean   <- bg_vals[!is.na(bg_vals) & bg_vals > 0]
  site_clean <- site_vals[!is.na(site_vals) & site_vals > 0]

  if (length(bg_clean) < min_detections |
      length(site_clean) < min_detections) {
    return(list(d = NA, note = "Insufficient detections"))
  }

  bg_log   <- log(bg_clean)
  site_log <- log(site_clean)

  result <- tryCatch(
    cohen.d(site_log, bg_log),
    error = function(e) NULL
  )

  if (is.null(result)) {
    return(list(d = NA, note = "Cohen's d computation failed"))
  }

  return(list(
    d    = abs(result$estimate),
    note = as.character(result$magnitude)
  ))
}

d_results   <- lapply(metals, compute_d)
d_values    <- sapply(d_results, function(x) x$d)
d_magnitude <- sapply(d_results, function(x) x$note)

# ── 7. RUN POWER ANALYSIS FOR EACH METAL ────────────────────

run_power <- function(n1, n2, d, alpha = 0.05) {
  if (is.na(d) | is.na(n1) | is.na(n2)) return(NA)
  if (n1 < min_detections | n2 < min_detections) return(NA)
  if (d == 0) return(NA)

  tryCatch(
    pwr.t2n.test(
      n1          = n1,
      n2          = n2,
      d           = d,
      sig.level   = alpha,
      alternative = "two.sided"
    )$power,
    error = function(e) NA
  )
}

power_vals <- mapply(
  run_power,
  n1 = n_bg_detected,
  n2 = n_site_detected,
  d  = d_values
)

# ── 8. MINIMUM N REQUIRED TO REACH 0.80 POWER ───────────────
# Solves for equal-group n needed for power = 0.80.
# Use this to recommend sample counts for future investigations.

min_n_required <- sapply(d_values, function(d) {
  if (is.na(d) | d == 0) return(NA)
  tryCatch(
    ceiling(pwr.t.test(
      d           = d,
      sig.level   = 0.05,
      power       = 0.80,
      type        = "two.sample",
      alternative = "two.sided"
    )$n),
    error = function(e) NA
  )
})

# ── 9. ASSEMBLE FULL RESULTS TABLE ──────────────────────────

results <- data.frame(
  Metal                = metals,
  N_BG_Total           = n_bg_total,
  N_Site_Total         = n_site_total,
  N_BG_Detected        = n_bg_detected,
  N_Site_Detected      = n_site_detected,
  Pct_BG_Detected      = detected_counts$Pct_BG_Detected,
  Pct_Site_Detected    = detected_counts$Pct_Site_Detected,
  Cohens_d             = round(d_values, 3),
  Effect_Magnitude     = d_magnitude,
  Achieved_Power       = round(power_vals, 3),
  Beta                 = round(1 - power_vals, 3),
  Min_N_Per_Group_0.80 = min_n_required,
  row.names            = NULL,
  stringsAsFactors     = FALSE
) %>%
  mutate(
    Adequacy = case_when(
      is.na(Achieved_Power)             ~ "Insufficient data",
      Achieved_Power >= 0.80            ~ "Adequate",
      Achieved_Power >= 0.50            ~ "Marginal",
      TRUE                              ~ "Underpowered"
    ),
    Flag = case_when(
      N_BG_Detected < min_detections   ~ paste0("BG detections below minimum (n=", N_BG_Detected, ")"),
      N_Site_Detected < min_detections ~ paste0("Site detections below minimum (n=", N_Site_Detected, ")"),
      is.na(Achieved_Power)            ~ "Power could not be computed",
      TRUE                             ~ ""
    )
  )

cat("\n── Power analysis results ─────────────────────────────\n")
print(results %>% select(Metal, N_BG_Detected, N_Site_Detected,
                         Cohens_d, Achieved_Power, Beta,
                         Adequacy, Flag))

write_xlsx(results, "outputs/power_analysis_results.xlsx")
cat("\nFull results exported to outputs/power_analysis_results.xlsx\n")

# ── 10. BAR CHART — ACHIEVED POWER BY METAL ─────────────────

plot_data <- results %>%
  mutate(
    Metal = factor(Metal, levels = Metal[order(
      ifelse(is.na(Achieved_Power), -1, Achieved_Power)
    )])
  )

ggplot(plot_data, aes(x = Achieved_Power, y = Metal, fill = Adequacy)) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = 0.80, linetype = "dashed",
             color = "#C00000", linewidth = 0.9) +
  annotate("text", x = 0.82, y = 0.7,
           label = "Target power = 0.80\n(beta = 0.20)",
           color = "#C00000", size = 3.2, hjust = 0) +
  geom_text(
    aes(label = ifelse(
      is.na(Achieved_Power),
      "Insufficient\ndetections",
      paste0(round(Achieved_Power * 100, 0), "%")
    )),
    hjust = -0.1, size = 3.2, color = "#333333"
  ) +
  scale_fill_manual(values = c(
    "Adequate"           = "#375623",
    "Marginal"           = "#854F0B",
    "Underpowered"       = "#C00000",
    "Insufficient data"  = "#CCCCCC"
  )) +
  scale_x_continuous(
    limits = c(0, 1.2),
    breaks = seq(0, 1, 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title    = "Achieved Statistical Power by Metal",
    subtitle = paste0(
      "Two-sample t-test approximation for Wilcoxon (ARE = 0.955)\n",
      "alpha = 0.05, two-tailed | Detected values only | ",
      "Log-transformed concentrations"
    ),
    x     = "Achieved power (1 - beta)",
    y     = NULL,
    fill  = "Adequacy"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", color = "#1F4E79", size = 13),
    plot.subtitle      = element_text(color = "#595959", size = 9.5),
    legend.position    = "bottom",
    panel.grid.major.y = element_blank()
  )

ggsave("outputs/power_analysis_barplot.png", width = 9, height = 6,
       dpi = 300, bg = "white")
cat("Bar chart saved to outputs/power_analysis_barplot.png\n")

# ── 11. POWER CURVES BY SAMPLE SIZE ─────────────────────────
# Shows how power changes with n for metals with sufficient data.
# Dotted vertical lines show actual average detected n per metal.

metals_with_d <- metals[!is.na(d_values)]
d_with_d      <- d_values[!is.na(d_values)]
metals_excl   <- metals[is.na(d_values)]

n_seq <- seq(3, 150, by = 2)

curve_data <- lapply(seq_along(metals_with_d), function(i) {
  d <- d_with_d[i]
  power_seq <- sapply(n_seq, function(n) {
    tryCatch(
      pwr.t.test(
        n           = n,
        d           = d,
        sig.level   = 0.05,
        type        = "two.sample",
        alternative = "two.sided"
      )$power,
      error = function(e) NA
    )
  })
  data.frame(Metal = metals_with_d[i], n = n_seq, Power = power_seq)
}) %>% bind_rows()

actual_n <- results %>%
  filter(Metal %in% metals_with_d) %>%
  mutate(n_avg = (N_BG_Detected + N_Site_Detected) / 2) %>%
  select(Metal, n_avg)

ggplot(curve_data, aes(x = n, y = Power, color = Metal)) +
  geom_line(linewidth = 0.85) +
  geom_hline(yintercept = 0.80, linetype = "dashed",
             color = "#C00000", linewidth = 0.8) +
  geom_vline(data = actual_n,
             aes(xintercept = n_avg, color = Metal),
             linetype = "dotted", linewidth = 0.6, alpha = 0.7) +
  annotate("text", x = 147, y = 0.83,
           label = "Power = 0.80", color = "#C00000",
           size = 3, hjust = 1) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  labs(
    title    = "Power Curves by Metal and Sample Size",
    subtitle = paste0(
      "Solid lines: power as n increases | Dotted verticals: actual average detected n\n",
      "alpha = 0.05, two-tailed, equal group sizes assumed\n",
      ifelse(length(metals_excl) > 0,
             paste("Excluded (insufficient detections):", paste(metals_excl, collapse = ", ")),
             "All metals included")
    ),
    x     = "Detected sample size per group (n)",
    y     = "Statistical power (1 - beta)",
    color = "Metal"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", color = "#1F4E79", size = 13),
    plot.subtitle = element_text(color = "#595959", size = 9),
    legend.position = "right"
  )

ggsave("outputs/power_curves_by_metal.png", width = 11, height = 6,
       dpi = 300, bg = "white")
cat("Power curves saved to outputs/power_curves_by_metal.png\n")

# ── 12. SESSION INFO ─────────────────────────────────────────

cat("\n── Session info ───────────────────────────────────────\n")
print(sessionInfo())
