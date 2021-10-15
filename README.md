# Sun-et-al-2021-protein-coding-variants-in-human-disease
R package to generate results as per Sun et al,. 2021, Genetic associations of protein-coding variants in human disease


## Functions
* allelic_het_simulator - simulates allelic heterogenity between two population cohorts and computes the resuting IVW-uplift.
* surface_plots - outputs surface plots for the observed (i.e., simulated) and expected (i.e., theoretically predicted) values for IVW-uplift across a range of MAFs and allelic enrichment values.

## Installation
1. install.packages("devtools")
2. library(devtools)
3. install_github("cnfoley/Sun-et-al-2021-protein-coding-variants-in-human-disease", build_vignettes = FALSE)
4. install_github("cnfoley/mrclust", build_vignettes = TRUE)
5. library(Sun-et-al-2021-protein-coding-variants-in-human-disease)
6. browseVignettes("Sun-et-al-2021-protein-coding-variants-in-human-disease")

