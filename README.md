# Power Analysis for Background Determination of Metals in Environmental Forensics

A reproducible R workflow for evaluating the **statistical power** of two-sample comparisons (site vs. background) on an analyte-by-analyte basis. The workflow is designed to support defensible background determinations of metals in soil under regulatory and forensic scrutiny, where a non-significant Wilcoxon test result must be backed by a demonstration that the test actually had a fair chance of detecting a real difference.

---

## What this project does

Given a tidy CSV of soil sample results — one row per sample, one column per metal, with a `Population` column labelled `Background` or `Site` — this workflow:

1. **Counts detections** per metal in each population (non-detects entered as `NA`).
2. **Calculates Cohen's *d*** on log-transformed detected values to express the magnitude of difference between site and background relative to natural variability.
3. **Runs a two-sample power calculation** for each metal using the unequal-n *t*-test as a defensible approximation for the Wilcoxon rank-sum test (asymptotic relative efficiency, ARE = 0.955).
4. **Reports achieved power** (1 − β) per metal and classifies adequacy as **Adequate** (≥ 0.80), **Marginal** (0.50 – 0.80), **Underpowered** (< 0.50), or **Insufficient data** (< 3 detections in either group).
5. **Solves for the minimum sample size per group** needed to reach the conventional power target of 0.80, which directly supports recommendations for follow-up sampling.
6. **Exports** an Excel results table, a detectability summary, a bar chart of achieved power, and a set of power curves showing how power scales with sample size for each metal.

---

## Why this matters

In environmental forensics, the question is rarely *"is there a statistically significant difference between site and background?"* — it is **"can I defend a non-significant result as evidence that the populations are the same?"**

That second question depends entirely on statistical power. Without a power analysis:

- A non-significant Wilcoxon result for a sparsely-detected metal (e.g., silver with only 2 background detections) is **uninterpretable**. It is just as consistent with "no difference exists" as with "a difference exists but the test could not see it through the noise of small *n*."
- Regulators, opposing experts, and courts can readily challenge any claim of equivalence that is not backed by a quantified probability of Type II error (β).
- Recommendations for additional sampling become arbitrary instead of grounded in a calculated minimum-detectable effect.

This workflow converts a vague concern about sample size into a **specific, defensible statement** that can be placed in a report, an expert affidavit, or litigation support documentation: *"Given my detected sample sizes and the observed variability, my test could reliably detect a difference of magnitude X. Anything smaller would have been missed."*

---

## The statistical concepts in plain language

**Alpha (α)** is the probability of declaring two groups different when they are actually the same — the false alarm rate. Conventionally set to 0.05.

**Beta (β)** is the probability of declaring two groups the same when they are actually different — the missed-difference rate. This is the error that matters when defending a non-significant result.

**Power (1 − β)** is the probability of correctly detecting a real difference when one exists. The conventional minimum is 0.80, meaning an 80% chance of catching a real effect.

A useful analogy: imagine trying to hear a faint sound in a noisy room. The "sound" is the difference between site and background. The "noise" is the variability in the data. β is how likely you are to miss the sound. More samples means a better chance of hearing a faint signal; fewer samples means the signal can be present and still go undetected.

---

## Repository contents

```
power-analysis-background-determination/
├── README.md                   # This file
├── LICENSE                     # MIT License
├── .gitignore                  # R / RStudio ignores
├── R/
│   └── Power_Analysis.R        # Main analysis script
├── data/
│   └── example_data_template.csv   # Empty CSV template with required columns
├── shiny/
│   └── README.md               # Plans for the upcoming Shiny app
└── outputs/                    # Created at runtime — generated tables and figures
```

---

## Quick start

### Prerequisites

- R version 4.3 or higher
- The following packages: `pwr`, `effsize`, `ggplot2`, `dplyr`, `tidyr`, `writexl`, `scales`

Install once:

```r
install.packages(c("pwr", "effsize", "ggplot2", "dplyr",
                   "tidyr", "writexl", "scales"))
```

### Input data format

Your CSV must have:

- One row per sample.
- A `Sample_ID` column.
- A `Population` column containing exactly `Background` or `Site` (case sensitive).
- One column per metal, named to match the `metals` vector in the script (`Arsenic`, `Barium`, `Cadmium`, `Chromium`, `Copper`, `Lead`, `Nickel`, `Selenium`, `Silver`, `Zinc`).
- **Non-detects entered as blank or `NA`** — not zero, not the detection limit, and not half the detection limit. The script counts non-`NA` entries as detections.
- Detected values entered as positive numeric concentrations (mg/kg).

A blank template is provided in `data/example_data_template.csv`.

### Running the analysis

From the project root:

```r
setwd("path/to/power-analysis-background-determination")
source("R/Power_Analysis.R")
```

Or open `R/Power_Analysis.R` in RStudio, set the working directory to the project root, and source the file.

### Outputs

The script writes the following to the working directory:

| File | Description |
|------|-------------|
| `detectability_summary.xlsx` | Detection counts and percentages per metal |
| `power_analysis_results.xlsx` | Full results: Cohen's *d*, achieved power, β, adequacy classification, minimum *n* for 0.80 power |
| `power_analysis_barplot.png` | Bar chart of achieved power per metal with 0.80 reference line |
| `power_curves_by_metal.png` | Power-vs-sample-size curves with dotted verticals at the current detected *n* |

---

## Reading the results

The bar chart is the headline. Each bar shows the achieved power for one metal. The red dashed line at 0.80 is the conventional adequacy threshold. Metals to the right of the line carry interpretive weight; metals to the left do not.

The power curves are the supporting evidence. For each metal, follow the dotted vertical line (your current detected *n*) up to where it crosses the solid curve. The *y*-value at that intersection is your achieved power. The curves also let you read off the *n* you would need to reach the 0.80 target — answering the practical question "how many more samples do I need to collect?"

---

## Methodological notes

**Why the *t*-test approximation for a Wilcoxon power calculation.** The asymptotic relative efficiency of the Wilcoxon rank-sum test relative to the *t*-test is 0.955 under normality and is at least 0.864 for any continuous distribution (Hollander, Wolfe, & Chicken, 2014). Using the *t*-test power calculation therefore gives a slightly conservative — and defensible — approximation of Wilcoxon power. For high-scrutiny deliverables where an exact Wilcoxon power is preferred, the `wmwpow` package implements a direct simulation-based calculation.

**Why log-transformed concentrations.** Environmental concentration data is typically right-skewed and approximately log-normal. Cohen's *d* on log-transformed values produces a stable, scale-invariant effect size that is appropriate for the underlying distribution and consistent with the Shapiro-Wilk normality testing performed in the broader background-determination SOP.

**Why "detected *n*" rather than total *n*.** A non-detect contributes no information about the magnitude of difference between populations. The effective sample size for an effect-size comparison is the number of actual measurements, not the number of laboratory tubes that came back. Metals with fewer than three detections in either group are flagged as having insufficient data for a meaningful power calculation — itself a finding worth reporting.

---

## Roadmap — Shiny app (in progress)

The next phase wraps this workflow in an interactive Shiny app that lets a user:

- Upload a CSV directly.
- Map their column names to the expected metals.
- Choose alpha and the target power threshold.
- View the bar chart, power curves, and results table interactively.
- Export everything as a single ZIP.

The app will be deployable to `shinyapps.io` under the project owner's credentials. See `shiny/README.md` for the development plan.

---

## Citation

If you use this workflow, please cite the methodology source for the Wilcoxon ARE argument:

> Hollander, M., Wolfe, D. A., & Chicken, E. (2014). *Nonparametric Statistical Methods* (3rd ed.). Wiley.

And the `pwr` package for the underlying power calculations:

> Champely, S. (2020). *pwr: Basic Functions for Power Analysis*. R package version 1.3-0.

---

## License

This project is released under the MIT License. See `LICENSE` for details.

---

## Author

Mike Dereviankin — environmental forensics workflow development.
