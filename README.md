Replication code for Inference with Multi-Outcome Random Forest.

Run each script in the sequence given below, assuming you have installed R/Rstudio. Ensure you have installed all packages in [requirements.txt](https://github.com/marnare/inferenceMGRF/blob/main/requirements.txt) Scripts are run on MacBook Pro with Apple M3 Max chip, 64 GB Memory. 

1.  Clone the github repository in the local directory: git clone git@github.com:marnare/inferenceMGRF.git

2.  Confidence Ellipses (Figure 1): [Confidence Regions Jupyer Notebook](https://github.com/marnare/inferenceMGRF/blob/main/confidence_regions/confidence_regions.ipynb)

3.  Type 1 and Type 2 errors (Table 2 and Table 3):

  - Press Command + F and set run_chunk = TRUE wherever indicated. Run [this rmarkdown](https://github.com/marnare/inferenceMGRF/blob/main/confidence_regions/confidenceEllipses_sims_samplesizes.Rmd) that calculates average Type 1 and Type 2 error rates across different sample sizes and correlation of treatment effects under constant covariance assumption (with 500 simulated experiments); [Replicated html version](https://github.com/marnare/inferenceMGRF/blob/main/confidence_regions/confidenceEllipses_sims_samplesizes.html) (Completion time: > 7 hours).
  - Press Command + F and set run_chunk = TRUE wherever indicated. Run [this markdown](https://github.com/marnare/inferenceMGRF/blob/main/confidence_regions/confidenceEllipses_sims_samplesizes_personalized.Rmd) that calculates average Type 1 and Type 2 error rates across different sample sizes and correlation of treatment effects under heterogeneous covariance assumption (with 500 simulated experiments) [Replicated html version](https://github.com/marnare/inferenceMGRF/blob/main/confidence_regions/confidenceEllipses_sims_samplesizes_personalized.html) (Completion time: > 7 h).  
  - Press Command + F and set run_chunk = TRUE wherever indicated. Run [this markdown](https://github.com/marnare/inferenceMGRF/blob/main/confidence_regions/table_risks.Rmd) that formats Tables 2 and 3 and saves them as a tex output. 

4. Figure 2:
  - Press Command + F and set run_chunk = TRUE wherever indicated. Run [this markdown](https://github.com/marnare/inferenceMGRF/blob/main/confidence_regions/confidenceEllipses_outcomes.Rmd) that produces a plot of Type 1 and Type 2 error rates. [Replicated html version](https://github.com/marnare/inferenceMGRF/blob/main/confidence_regions/confidenceEllipses_outcomes.html). 

