# Inference with Multi-Outcome Random Forest

This repository provides replication code for our analysis of statistical inference in multi-outcome random forests, focusing on confidence regions and error rate estimation.

## Requirements

- R/RStudio (latest version)
- Required packages: See [requirements.txt](requirements.txt)
- Tested on: MacBook Pro (M3 Max, 64GB RAM)

## Replication Steps

### 1. Clone Repository
```bash
git clone git@github.com:marnare/inferenceMGRF.git
```

### 2. Confidence Regions
📊 [Confidence Regions Notebook](confidence_regions/confidence_regions.ipynb)
- Reproduces Figure 1
- Demonstrates confidence ellipse construction

### 3. Error Rate Analysis
> Note: Set `run_chunk = TRUE` in each file. Expected runtime: ~12 hours per file.

1. **Constant Covariance**
   - [Analysis](confidence_regions/confidenceEllipses_sims_samplesizes.Rmd)
   - [Results](confidence_regions/confidenceEllipses_sims_samplesizes.html)
   - 500 simulated experiments

2. **Heterogeneous Covariance**
   - [Analysis](confidence_regions/confidenceEllipses_sims_samplesizes_personalized.Rmd)
   - [Results](confidence_regions/confidenceEllipses_sims_samplesizes_personalized.html)

3. **Results Tables**
   - [Table Generation](confidence_regions/table_risks.Rmd)
   - Produces Tables 2 & 3

### 4. Visualizations
📈 [Error Rate Plots](confidence_regions/confidenceEllipses_outcomes.Rmd)
- Generates Figure 2
- [View Results](confidence_regions/confidenceEllipses_outcomes.html)

### 5. Application
🔬 [Empirical Analysis](application/application.Rmd)
- Reproduces Figure 3
- [View Results](application/application.html)

## Questions?
Open an issue or contact the authors.

---
For detailed methodology, please see our accompanying paper.
