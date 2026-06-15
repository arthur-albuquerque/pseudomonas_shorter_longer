# Shorter Antibiotic Duration for Pseudomonas Bacteremia
Arthur M Albuquerque, Carolina Santolia, Fernando G Zampieri, James M Brophy and Todd C Lee

A two-stage Bayesian fixed-effects meta-analysis testing whether *Pseudomonas aeruginosa* modifies the effect of a 7-day versus 14-day antibiotic course on all-cause mortality in gram-negative bacteremia.

## Files

- [`01_main.R`](01_main.R) — Main analysis script. Loads the data, computes Stage-1 per-trial effect sizes and within-trial interactions, runs the Stage-2 Bayesian fixed-effects meta-analyses, and produces the publication figures.
- [`report.qmd`](report.qmd) — Quarto source for the full analysis report.
- [`report.html`](https://htmlpreview.github.io/?https://github.com/arthur-albuquerque/pseudomonas_shorter_longer/blob/main/report.html) — Full analysis report (click to view in browser).
- [`pseudomonas_data.xlsx`](pseudomonas_data.xlsx) — Trial-level event data 
- [`figures/`](figures) — Output figures: Figure 1 (data-only forest plot) and Figure 2 (pooled posteriors), each as PNG and TIFF
