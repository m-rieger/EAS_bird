# EAS

*add link to manuscript (DOI)*

This repository contains data and code to reproduce results from *add Title*

## Abstract


## Content
- **OEFS_Paper_Models.R**: script for the complete modelling process.
- /**R**: contains code for helper functions in R. 
- /**data**: contains data ready for analysis.
    - OEFS_Paper_LandscapeNRW.csv: landscape shares of NRW.
    - OEFS_Paper_RawData.csv: raw data from EAS surveys for selected species
    - OEFS_Paper_SiteInfo.csv: raw data of site specific parameters like weather data and landscapes
    - OEFS_Paper_SpeciesInfo.csv: species-specific information for modelling process (family, binomial coefficients, ...)
- /**01_models**: This is where RDS-files of models will be saved.

## Required packages outside of CRAN
- `cmdstanr` package. Install with `devtools::install_github("stan-dev/cmdstanr")` and, for the first time you install it, `cmdstanr::install_cmdstan()`. You can find more information [here](https://mc-stan.org/cmdstanr/index.html).
