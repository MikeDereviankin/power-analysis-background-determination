# Shiny App — Power Analysis for Background Determination

This folder is the home of the interactive Shiny version of the power analysis workflow. The app wraps the R script in `R/Power_Analysis.R` so that a user with no R experience can run the same analysis from a browser.

## Planned features

- **Upload a CSV** of soil sample data directly through the browser.
- **Column mapping interface** that lets the user point the app at their `Population` column and their metal columns even when names do not match the default vector.
- **Configurable thresholds** — alpha, target power, and minimum detection count are all exposed as sliders / inputs rather than hard-coded.
- **Interactive results table** with sortable columns and adequacy colour coding.
- **Interactive plots** — achieved power bar chart and power-vs-*n* curves rendered with `plotly` so the user can hover for exact values.
- **Single-click export** of all results (Excel + PNGs) packaged as a ZIP.

## Planned file layout

```
shiny/
├── README.md           # This file
├── app.R               # Single-file Shiny app (UI + server)
├── R/
│   └── power_funcs.R   # Helper functions refactored out of Power_Analysis.R
├── www/                # Static assets (CSS, images, favicon)
└── rsconnect/          # Deployment metadata for shinyapps.io (gitignored)
```

## Deployment

The app will be deployed to `shinyapps.io` under the project owner's account.

```r
# One-time setup
install.packages("rsconnect")
rsconnect::setAccountInfo(
  name   = "<your-shinyapps-account>",
  token  = "<your-token>",
  secret = "<your-secret>"
)

# Deploy from this folder
rsconnect::deployApp(appDir = "shiny")
```

## Status

**Not yet built.** Workflow has been validated in the standalone R script (`R/Power_Analysis.R`). Next step is refactoring the script into reusable functions and wiring up the Shiny UI.
