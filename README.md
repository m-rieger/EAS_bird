# EAS_bird
![](https://img.shields.io/github/license/m-rieger/EAS_bird) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12518872.svg)](https://doi.org/10.5281/zenodo.12518872)


This repository contains exemplified data and code to reproduce results from Rieger, M.R., Grüneberg, C., Oberhaus, M., Trautmann, S., Parepa, M., Anthes, N. (2024). Bird population trend analyses for a monitoring scheme with a highly structured sampling design. Preprint []().

## Abstract

Population trends derived from systematic monitoring programmes are essential to identify species of conservation concern and to evaluate conservation measures. However, monitoring data pose several challenges for statistical analysis, including spatial bias due to an unbalanced sampling of landscapes or habitats, variation in observer expertise, frequent observer changes, and overdispersion or zero-inflation in the raw data. An additional challenge arises from so-called ‘rolling’ survey designs, where each site is only visited once within each multi-year rotation cycle. We developed a GAMM-based workflow that addresses these challenges and exemplify its application with the highly structured data from the Ecological Area Sampling (EAS) in the German federal state North-Rhine Westphalia (NRW). First, we derive a routine that allows informed decisions about the most appropriate combination of distribution family (Poisson or negative binomial), model covariates (e.g., habitat characteristics), and zero-inflation formulations to reflect species-specific data distributions. Second, we develop a correction factor that buffers population trend estimates for variation in observer expertise as reflected in variation in total bird abundance. Third, we integrate landscape-specific trends that adjust for between-year variation in the representation of habitat or landscape types within the yearly subset of sampled sites. In a consistency check, we found good match between our GAMM-based EAS trends and TRIM-based trends from the standard German common Bird monitoring scheme. The study provides a template script for R statistical software so the workflow can be adapted to other monitoring programmes with comparable survey designs and data structures.
*Keywords: breeding bird monitoring, trend analysis, rolling surveys, spatial bias*


## Content
### raw data
The data subset consists of four species (Common Blackbird, Common Chiffchaff, Common Kestrel, White Wagtail) that were modelled for both biogeographical regions and used as example species in the article.     
- `/data`: contains data ready for analysis.  
    - `EAS_bird_LandscapeNRW.csv`: landscape shares of NRW. Columns: *landscape* (natural region and metropolitan area combination), *region* (biogeographical region), *area_NRW* (area in km² per landscape)  
    - `EAS_bird_RawData.csv`: raw data from EAS surveys for selected species. Columns: *abundance* (number of territories per km²), *exclude* (exclusion of data due to an internal validation process of the LANUV (*yes*)), *ID* (site ID), *year* (survey year), *region* (biogeographical region, atl = atlantic, kon = continental), *landcape* (natural region and metropolitan area combination), *species*, *species_sci* (English and scientific species names. The species 'other' contains the summed up abundances of all species that are not listed in the data set to compute observer effects based on total abundances), *obs.ID* (observer ID)  
    - `EAS_bird_SiteInfo.csv`: raw data of site specific parameters like weather data and landscapes. Columns: *Tmean_*, *Tmin_*, *Tmax_*, *rain_*, *sun_* (long-term averages of mean, minimum and maximum temperatures (°C), precipitation and sunshine duration (h) for spring and winter), *altitude_asl* (altitude above sea level), *arable*, *forest*, *settlement* (portion of site (%) covered by arable land, forests, or settlement areas)   
    - `EAS_bird_SpeciesInfo.csv`: species-specific information for modelling process. Columns: *species*, *species_sci* (English and scientific species names), *modR* (biogeographical region(s) used for modelling), *modZI* (binomial model coefficients), *modFAM* (model family, pois = Poisson, nb = negative binomial, zip = zero-inflated Poisson, zinb = zero-inflated negative binomial), *modFINAL* (indicates whether model is used as final model (TRUE) or to exemplify model selection (FALSE))    
    
### data that will be created
- `/01_models`: This is where RDS-files of models as well as 
    - convergence parameters (`Convergence_*.csv`) 
    - k fold-cross validation results (`Kfold_*.csv`) (optional) and
    - model statistics (e.g. comparison of observed and predicted distribution parameters) (`Statistics_*.csv`) 
    will be saved.
- `/02_output`: This is where 
    - posterior predictions of trends (`PosteriorPredictions_*.csv`) 
    - longterm trends (`LongtermTrend_*.csv`) 
    - coefficient estimates (`CoefficientEstimates_*.csv`)
    - pairwise comparisons (of neighbouring years) (`PairWise_*.csv`)
    - PCA scores (`PCA_scores_*.csv`)
    - rotated PCA loadings (`PCA_rotated_loadings_*.csv`)
    - data used for modelling (`Data_modelling.csv`) and
    - posterior prediction plots per species (`PostPred_*.pdf`) (optional)
    will be saved.

### code
- `EAS_bird.R`: script for data preparation, modelling process, posterior predictive checks and posterior predictions.
- `EAS_bird_plots.R`: script with helper functions for plotting

## Required packages outside of CRAN
- `EAS` package (helper functions to analyse population trends from monitoring data with highly structured sampling designs). Installed by default via `devtools::install_github("m-rieger/EAS")`  
- `cmdstanr` package (optional, for parallelized evaluation, otherwise rstan will be used). Install with `devtools::install_github("stan-dev/cmdstanr")` and, for the first time you install it, `cmdstanr::install_cmdstan()`. You can find more information [here](https://mc-stan.org/cmdstanr/index.html).

