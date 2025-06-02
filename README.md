# Replication package: *Male monopolization and reproductive skew in a tolerant multilevel society*

Analyzing long-term data on wild Guinea baboons to investigate paternity success, reproductive skew, and the effect of age and rank on reproductive success.


The major requirement for our code to run is a working `cmdstanr`/`Stan` setup.
Check out https://mc-stan.org/cmdstanr/ for help on that.

In order to replicate some of the supplementary analyses with respect to reproductive skew, we need also the SkewCalc package (https://github.com/ctross/SkewCalc).
The primary analyses, however, do not require that package.

Also, some of our plotting code requires an additional package `viridisLite`.

`install.packages("viridisLite")`

In order to get the code to work, we need to unpack the replication package into our working directory.
If we then run the following code, we should obtain a `TRUE`:

`file.exists("stan_models/rs_model.stan")`

If that works, we can proceed to the actual scripts.

## Note on reproducibility

Exact reproducibility is hard to achieve using Stan (see [here](https://mc-stan.org/docs/reference-manual/reproducibility.html) for details).
Nevertheless, we used fixed seeds for R's and Stan's random number generators throughout in an attempt to minimize discrepancies over different computer setups.
For example, the numeric values in table 1 were within a range of $\pm0.09$ over several runs on different machines and for table 2 they were within a range of $\pm0.03$ over several runs.

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
  
