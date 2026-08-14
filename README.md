# Replication package: *Male monopolization and reproductive skew in a tolerant multilevel society*

[![DOI](https://zenodo.org/badge/991280096.svg)](https://doi.org/10.5281/zenodo.15582801)

Published in: . . . 

DOI: . . .

Date: 2026-08-14



Analyzing long-term data on wild Guinea baboons to investigate paternity success, reproductive skew, and the effect of age and rank on reproductive success.

The major requirement for our code to run is a working `cmdstanr`/`Stan` setup.
Check out https://mc-stan.org/cmdstanr/ for help on that.

In order to replicate some of the supplementary analyses with respect to reproductive skew, we need also the `SkewCalc` package (https://github.com/ctross/SkewCalc).
The primary analyses, however, do not require that package.

Also, some of our plotting code requires an additional package `viridisLite`.

`install.packages("viridisLite")`

We recommend using RStudio and opening the unzipped folder as an *RStudio project*, which handles working directories.
Alternatively, we unpack the replication package and set our working directory manually to the folder like so:

`setwd("~/Documents/path/to/unzipped")` # adapt according to your local settings

Regardless of how we deal with the working directory, if we then run the following code, we should obtain a `TRUE`:

`file.exists("stan_models/rs_model.stan")`

If that works, we can proceed to the actual scripts.

## Reproductive success models

Script 02 contains the model code to reproduce the results in the main text.
Scripts 00 and 01 contain code for data preparation and prior simulations, and both these scripts are not required to be run in order to regenerate our model results.

## Skew (models)

Script 03 contains code for data preparation and fitting of the reproductive skew models.

# File descriptions

## Scripts

  - `00_dataprep.R`: reading and pre-processing of the data for reproductive success model
  
  - `01_prior_simulation.R`: prior simulations for our reproductive success model
  
  - `02_models.R`: actual reproductive success model (incl figures)
  
  - `03_skew.R`: data processing and models for skew analysis (incl figures)

## Data files

  - `data/domdata.csv` (714 rows):
  
    Dyadic dominance interactions with decided outcome.
    
    * `$date`: date of interaction
    
    * `$party`: party/group in which interaction occured
    
    * `$winner`, `$loser`: male id codes
    

  
  - `data/presence_matrix.csv` (714 rows x 72 columns):
    
    * This table maps male party/group membership to dominance interactions.
    
    * Each row corresponds to one dominance interaction (i.e. we have as many row as in `domdata.csv`).
    
    * A value of 1 in this matrix indicates that the male in the column was present on the date of the interaction (row).
  
    * Each column corresponds to one male-party-stint. For example `M04@six_1` is the first stint of male 04 and it represent his presence in party 5. The next column, `M04@sixi_2` is his second overall stint, but this time in party 6I (sixi).
    
  
  - `data/rs_data.csv` (190 rows):
  
    * `$id`: male id
    
    * `$year`: year
    
    * `$party`: party
    
    * `$nfem`: number of females in the male's unit
    
    * `$age_cat`: prime or non-prime
    
    * `$widename`: maps the row in this table to names which include male party association
    
    * `$ratindex`: auxiliary index to map ratings in the Stan model code to outcome of the female count model
  
  - `data/skewdata.csv` (45 rows):
  
    * `$party`: party
    
    * `$duration`: number of alive days during the study
    
    * `$offspring`: number of offspring sired

## Stan model files
  
  - `stan_models/rs_model.stan`: model code to fit Elo/reproductive success model
  
  - `stan_models/mindex_original.stan`: model code to calculate $M$ index
  
  - `stan_models/mindex_unconstrained.stan`: our adaptation of the $M$ index with unconstrained *alpha* parameter vector
  
## Note on reproducibility

Exact reproducibility is hard to achieve using Stan (see [here](https://mc-stan.org/docs/reference-manual/reproducibility.html) for details).
Nevertheless, we used fixed seeds for R's and Stan's random number generators throughout in an attempt to minimize discrepancies over different computer setups.
For example, the numeric values in table 1 were within a range of $\pm0.09$ over several runs on different machines and for table 2 they were within a range of $\pm0.03$ over several runs.

## Session info

The scripts were successfully run last on 2026-08-14 with:

```
R version 4.6.1 (2026-06-24)
Platform: aarch64-apple-darwin23
Running under: macOS Tahoe 26.5.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] cmdstanr_0.9.0

loaded via a namespace (and not attached):
 [1] vctrs_0.7.3          cli_3.6.6            knitr_1.51           rlang_1.2.0          xfun_0.59            otel_0.2.0          
 [7] processx_3.9.0       generics_0.1.4       tensorA_0.36.2.1     data.table_1.18.4    jsonlite_2.0.0       glue_1.8.1          
[13] backports_1.5.1      distributional_0.8.1 ps_1.9.3             grid_4.6.1           evaluate_1.0.5       tibble_3.3.1        
[19] abind_1.4-8          lifecycle_1.0.5      compiler_4.6.1       posterior_1.7.1      coda_0.19-4.1        pkgconfig_2.0.3     
[25] rstudioapi_0.19.0    lattice_0.22-9       viridisLite_0.4.3    R6_2.6.1             pillar_1.11.1        magrittr_2.0.5      
[31] checkmate_2.3.4      tools_4.6.1          withr_3.0.3          SkewCalc_1.0         matrixStats_1.5.0   
```

