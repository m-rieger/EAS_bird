## Functions used in the OEFS Paper ##
######################################-

## author: Mirjam R. Rieger
## latest update: 2023-06-03 (MR)
  # 2023-06-03 (MR): added add.zeros-function
  # 2023-06-16 (MR): added option (FALSE/TRUE) to obs.eff and weight to get additional dataframe or not
  #                  and updated weight-function to calculate weights per region, too
  # 2023-07-03 (MR): added some details to mod.stat()
  # 2023-07-28 (MR): minor changes to mod.stat and mod.conv


######################################################-
#### 1) weights based on habitat representativity ####
######################################################-

## This function calculates weights based on annual habitat representativity. 
## If a habitat was underrepresented in a given year, the annual weight will be >1
## If a habitat was overrepresented  in a given year, the annual weight will be <1
## If a habitat was represented according to its share in a given year, the annual weight will be 1

## input:
#########-

## df.habitat (default = NULL):       dataframe including present "habitats" and the total "area" of the habitat of the studied region (this is the target share)
## habitat    (default = "habitat"):  columname of habitat in df.habitat and df.data
## area       (default = "area"):     columname of area (share) in df.habitat
## df.data    (default = NULL):       dataframe including site "ID" for each "year" the site was surveyed and the assigned "habitat".
##                                      you can use your raw data per species since the function condenses the df to its unique site-year combination
## ID         (default = "ID"):       columname of site ID in df.data
## year       (default = "year"):     columname of year in df.data
## by_spec    (default = FALSE):      logical indicating whether weights should be calculated per species 
##                                      (e.g., if observations were excluded species-specific resulting in different annual site-subsets per species)
## species    (default = "species"):  columname of species name in df.data (only needed if by_spec = TRUE)
## add.df     (default = TRUE):       if TRUE, returns a dataframe with weights per habitat and year (and species)

## output:
##########-

## the function weight returns a dataframe including the implemented raw data "df.data" as well as
   ## "weight" per habitat and year (and species if by_spec = TRUE) and the default columns "habitat", "ID", "year" (and "species"). 
   ## If by_spec = FALSE, the weights per habitat and year are the same for all species. 
## The function optionally returns the dataframe "df.w" containing weights per habitat and year (and species if by_spec = TRUE). This dataframe has the following columns:
   ## year:     year of survey
   ## habitat:  habitat
   ## N.sites:  number of surveyed sites per habitat and year (and species)
   ## actual:   actual annual share of the habitat: number of surveyed sites per habitat and year (and species)/total number of surveyed sites per year (and species)
   ## target:   target share of the habitat
   ## weight:   weight for each habitat-year (or habitat-year-species) combination
   ## (species: species name (only if by_spec = TRUE))


weight <- function(df.habitat = NULL, habitat = "habitat", area = "area",
                   df.data = NULL, ID = "ID", year = "year",
                   by_spec = FALSE, species = "species",
                   by_reg = FALSE, region = "region",
                   add.df = TRUE) {
  
  require(dplyr)
  
  ## check if data structure is correct
  #####################################-
  
  if(is.null(df.habitat)) stop("You need to define df.habitat.")
  if(is.null(df.data))    stop("You need to define df.data.")

  if(!(habitat %in% colnames(df.habitat))) stop(paste0("The column '", habitat, "' is missing in df.habitat."))
  if(!(area %in% colnames(df.habitat)))    stop(paste0("The column '", area, "' is missing in df.habitat."))
  if(!(habitat %in% colnames(df.data)))    stop(paste0("The column '", habitat, "' is missing in df.data."))
  if(!(ID %in% colnames(df.data)))         stop(paste0("The column '", ID, "' is missing in df.data."))
  if(!(year %in% colnames(df.data)))       stop(paste0("The column '", year, "' is missing in df.data."))
  
  if(by_spec) if(!(species %in% colnames(df.data)))  stop(paste0("The column '", species, "' is missing in df.data."))
  if(by_reg)  {
    if(!(region  %in% colnames(df.data)))    stop(paste0("The column '", region, "' is missing in df.data."))
    if(!(region  %in% colnames(df.habitat))) stop(paste0("The column '", region, "' is missing in df.habitat."))
  }
  
  ## calculate target share
  #########################-
  df.habitat$habitat <- df.habitat[, habitat]
  df.habitat$area    <- df.habitat[, area]
  
  ## target habitat per region (if by_reg = TRUE)
  if(by_reg) {
    df.habitat$region    <- df.habitat[, region]
    df.habitat <- df.habitat %>% group_by(habitat, region) %>%
      summarise(area = sum(area), .groups = "drop") %>%
      ungroup() %>%
      group_by(region) %>%
      mutate(target = area/sum(area))
  }
  
  ## target habitat (if by_reg = FALSE)
  if(by_reg == FALSE) {
    df.habitat <- df.habitat %>% group_by(habitat) %>%
    summarise(area = sum(area), .groups = "drop")
    df.habitat$target  <- df.habitat$area/sum(df.habitat$area)
  }
      
  ## calculate actual share
  #########################-
  df.data$habitat <- df.data[, habitat]
  df.data$ID      <- df.data[, ID]
  df.data$year    <- df.data[, year]
  
  ## define columns based on by_reg and by_spec
  if(by_reg == TRUE) {
    df.data$region    <- df.data[, region]
    
    col.list <- col.list1 <- c("habitat", "region")
    col.list2 <- c("year", "region")
  }
  
  if(by_reg == FALSE) {
    col.list <- col.list1 <- c("habitat")
    col.list2 <- c("year")
  }
  
  if(by_spec == TRUE) {
    df.data$species    <- df.data[, species]
    
    col.list <- c(col.list, "species")
    col.list2 <- c(col.list2, "species")
  }
  
  ## number of sites and actual share per habitat and year (and species) (and region)
  df.data2 <- unique(df.data[, c("ID", "year", col.list)])
  N.sites <- df.data2 %>% 
    group_by_at(c("year", col.list)) %>% 
    summarize(N.sites = n(), .groups = "drop") %>%
    ungroup() %>%
    group_by_at(col.list2) %>%
    mutate(actual = N.sites/sum(N.sites))
  
  ## check whether habitats are missing
  #####################################-
  
  if(by_reg == FALSE) {
    hab1 <- unique(df.habitat$habitat)
    hab2 <- unique(df.data2$habitat)
    
    for (i in 1:length(hab1)) {
      if(!(hab1[i] %in% hab2)) warning(paste0("Habitat ", hab1[i], " is not present in df.data. The dataset does not cover all present habitats."))
    }
    
    for (i in 1:length(hab2)) {
      if(!(hab2[i] %in% hab1)) stop(paste0("Habitat ", hab2[i], " is missing in df.habitat. Unable to calculate weights."))
    }
  }

  
  if(by_reg == TRUE) {
    
    for(r in unique(df.habitat$region)) {
      
      hab1 <- unique(df.habitat$habitat[df.habitat$region == r])
      hab2 <- unique(df.data2$habitat[df.data2$region == r])
      
      for (i in 1:length(hab1)) {
        if(!(hab1[i] %in% hab2)) warning(paste0("Region '", r, "': Habitat ", hab1[i], " is not present in df.data. The dataset does not cover all present habitats."))
      }
      
      for (i in 1:length(hab2)) {
        if(!(hab2[i] %in% hab1)) stop(paste0("Region '", r, "': Habitat ", hab2[i], " is missing in df.habitat. Unable to calculate weights."))
      }
    }
  }
  
  
  ## calculate weight
  ###################-
  df.w        <- left_join(N.sites, df.habitat[, c(col.list1, "target")], by = c(col.list1))
  df.w$weight <- df.w$target/df.w$actual
  
  ## add to dataframe
  ###################-
  df.weight <- left_join(df.data, df.w[, c("year", col.list, "weight")], by = c("year", col.list))

  ## return data
  if(add.df) df.w <<- df.w
  return(df.weight)

}

# df.W <- weight(df.habitat = df.habitat, habitat = "habitat", area = "Fläche",
#        df.data = df.data, ID = "FS_ID", year = "Jahr", by_spec = FALSE, species = "Name")
#   

################################################################################-
################################################################################-

## raw data for obs.eff
#######################-

# df.obs <- dat[, c("Name", "Jahr", "FS_ID", "Abundanz", "KartQ_prop", "Kart.Q25")]
# 
# ID <- "FS_ID"
# abundance <- "Abundanz"
# year = "Jahr"
# df.data <- df.obs
# OE <- 0.25



######################################################-
#### 2) observer effects based on total abundance ####
######################################################-

## This function calculates observer effects based on the proportional total abundance 
## (= total abundance per site and year / mean total abundance per site). 
## For each site and year, observer effects are classified in three classes using a defined threshold (default = 0.25) 
## If proportional abundance < (1-threshold), the effect is classified as "negative".
## If proportional abundance > (1+threshold), the effect is classified as "positive".
## If (1-threshold) <= proportional abundance <= (1+threshold), the effect is classified as "none".

## input:
#########-

## data      (default = NULL):        dataframe including "abundance" per site "ID" and survey "year" for all species (= raw data)
## ID        (default = "ID"):        columname of site ID in data
## abundance (default = "abundance"): columname of abundance in data
## year      (default = "year"):      columname of year in data
## OE        (default = 0.25):        threshold for classifying observer effects, should be a value >0 and <1.
## add.df    (default = TRUE):        if TRUE, returns a dataframe with observer effect per site and year

## output:
##########-

## the function obs.eff returns a dataframe including the implemented raw data "data" 
  ## as well as the additionally column for observer effect ("OE_xx") per site and year and the default column "abundance".

## The function optionally returns the dataframe "df.oe" containing observer effects per site and year if df.oe = TRUE. 
## This dataframe has the following columns:
  ## ID:            site ID
  ## year:          year of survey
  ## Ab.total:      total abundance per site and year (summed up across all species)
  ## mean.Ab.total: mean total abundance (Ab.total) per site 
  ## Ab.prop:       proportional abundance per site and year (= Ab.total / mean.Ab.total per site)
  ## OE_xx:         observer effect levels "negative", "none", and "positive" per site and year. xx is the used threshold in %.



obs.eff <- function(data = NULL, ID = "ID", abundance = "abundance", year = "year",
                    OE = 0.25,
                    add.df = TRUE) {
  
  require(dplyr)
  
  if(is.null(data))   stop("You need to define data")

  if(!(ID %in% colnames(data)))        stop(paste0("The column '", ID, "' is missing in data"))
  if(!(abundance %in% colnames(data))) stop(paste0("The column '", abundance, "' is missing in data"))
  if(!(year %in% colnames(data)))      stop(paste0("The column '", year, "' is missing in data"))
  if(OE <= 0 | OE >= 1)                stop("The threshold OE must be between 0 and 1.")
  
  ## define columns
  # data$ID        <- data[, ID]
  data$abundance <- data[, abundance]
  # data$year      <- data[, year]
  
  ## total and proportional abundance per site
  ############################################-
  df.oe <- data %>% group_by_at(.vars = c(ID, year)) %>% 
    summarize(Ab.total = sum(abundance), .groups = "drop") %>%
    group_by_at(.vars = ID) %>%
    mutate(mean.Ab.total = mean(Ab.total)) %>%
    ungroup() %>%
    mutate(Ab.prop = Ab.total/mean.Ab.total)
    
  ## define categories based on threshold OE
  ##########################################-
  df.oe$OE                         <- "none"
  df.oe$OE[df.oe$Ab.prop < (1-OE)] <- "negative"
  df.oe$OE[df.oe$Ab.prop > (1+OE)] <- "positive"
  
  df.oe$OE <- as.factor(df.oe$OE)
  
  ## change column names
  colnames(df.oe)[colnames(df.oe) == "OE"] <- paste0("OE_", 100*OE)
  
  ## return data
  if(add.df) df.oe <<- df.oe
  
  df.obs.eff <- left_join(data, df.oe[, c(ID, year, paste0("OE_", 100*OE))], by = c(ID, year))
  return(df.obs.eff)
  
}


# df.OE <- obs.eff(data = df.obs, ID = "FS_ID", abundance = "Abundanz", year = "Jahr")


######################################################-
#### 3) View model statistics (propZ, mean, ...)  ####
######################################################-


## input:
#########-

## model.list   (default = NULL):              list with different models (of one species) for comparisons
## model.name   (default = names(model.list)): names of models (e.g. family and zero-inflated coefficients)
## response     (default = NULL):              name of response variable (character)
## plot.stats   (default = FALSE):             if TRUE, graphical output is given
## spec         (default = NULL):              name of species (for title of plot)

## output:
##########-

## the function mod.stat returns a dataframe including model statistics (different quantiles of simulated values, observed value, BayesP-value) 
## for the proportion of zeros, minimum, maximum, mean, median, and SD of the response variable for different models

## The function optionally returns a graphical output inlcuding the observed value of the raw data 
## as well as median with 50% and 95% CrI of the simulated values if plot.stats = TRUE. 

mod.stat <- function(model.list = NULL, model.name = NULL, response = NULL,
                     plot.stats = FALSE, spec = NULL) {
  
  require(brms)
  ## define propZ-function
  prop_zero   <- function(z) sum(z == 0) / length(z)  
  
  if(is.null(model.list))  stop("You need to define the model-list.")
  if(!is.list(model.list)) stop("Model-list must be a list.")
  if(is.null(model.name)) model.name <- names(model.list)
  if(is.null(model.name)) model.name <- as.character(1:length(model.list)) # if the list is not named and still NULL
  if(is.null(response))    stop("You need to define the response.")

  ## define df.
  df.modS  <- data.frame(Stats = c("prop_zero", "min", "max", "mean", "median", "sd"))
  fun.list <-                  list(prop_zero,   min,   max,   mean,   median,   sd)
  
  df.modS.full <- NULL
  
  ## loop for different models (for comparison)
  for (m in 1:length(model.list)) {
    
    #mod <- model.list[m]
    ## define observed (yobs) and simulated (yrep) response
    yrep <- posterior_predict(model.list[[m]])
    yobs <- model.list[[m]]$data[, response]
  
    ## loop for different stats
    for (i in 1:length(fun.list)) {
      
      fun <- fun.list[[i]]
      
      stat_yrep <- apply(yrep, 1, fun)
      
      df.modS$lwr2.5[i]   <- quantile(stat_yrep, probs = 0.025)
      df.modS$lwr5[i]     <- quantile(stat_yrep, probs = 0.05)
      df.modS$lwr25[i]    <- quantile(stat_yrep, probs = 0.25) 
      df.modS$fit[i]      <- quantile(stat_yrep, probs = 0.5)
      df.modS$upr75[i]    <- quantile(stat_yrep, probs = 0.75) 
      df.modS$upr95[i]    <- quantile(stat_yrep, probs = 0.95)
      df.modS$upr97.5[i]  <- quantile(stat_yrep, probs = 0.975)
      df.modS$obs[i]      <- fun(yobs)
      df.modS$BayesP[i]   <- (length(stat_yrep[stat_yrep > fun(yobs)])/length(stat_yrep)) + (length(stat_yrep[stat_yrep == fun(yobs)])/length(stat_yrep)/2)
      df.modS$model       <- as.character(model.name[m])
      df.modS$species     <- spec
    }

    df.modS.full <- rbind(df.modS.full, df.modS)
  }
  
  
  if(plot.stats == TRUE) {
    require(ggplot2) 
    df.modS.full$Stats <- factor(df.modS.full$Stats, levels = c("prop_zero", "min", "max", "mean", "median", "sd"))
    
    print(ggplot(df.modS.full) +
      geom_hline(aes(yintercept = obs), color = "blue", linewidth = 2) +
      geom_pointrange(aes(x = model, ymin = lwr2.5, y = fit, ymax = upr97.5), 
                      color = "grey30", linewidth = 1) +
      geom_pointrange(aes(x = model, ymin = lwr25, y = fit, ymax = upr75), 
                      color = "grey30", linewidth = 2, fatten = 8) +
      ylab(paste0("function(", response, ")")) +
      xlab("model") +
      ggtitle(paste0(sp, " - model statistics of response variable '", response, "'"), 
              subtitle = "observed value/raw data (blue line) and 
simulated values (median with 50% (thick) and  95% (thin) CrI)") +
      facet_wrap(~ Stats, scales = "free_y") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  }
  
  return(df.modS.full)

}

#Effs <- mod.stat(model.list = list(mod2, mod2), response = "Ab.c", plot.stats = TRUE)


## df.weight (df.obs.eff) which is saved will overwrite any object called 'df.weight (df.obs.eff)
# e.g. add save.dat == TRUE/FALSE

###########################################-
#### 4) add zero abundances if missing ####
###########################################-

## This function adds zero abundances for each species-ID.year combination which is not present

## input:
#########-

## data         (default = NULL):        dataframe including "abundance" per "ID.year" for all "species" (= raw data)
## ID.year      (default = "ID.year"):   columname of ID.year in data
## species      (default = "species"):   columname of species names in data
## abundance    (default = "abundance"): columname of abundance in data
## ID.year.dep  (default = NULL): columns in data which values depend on ID.year (values will be duplicated for added zero abundances)
## species.dep  (default = NULL): columns in data which values depend on species (values will be duplicated for added zero abundances)

## output:
##########-

## The function returns the original dataframe extended by the newly added zero abundances.

add.zeros <- function(data = NULL, ID.year = "ID.year", species = "species", abundance = "abundance",
                      ID.year.dep = NULL,
                      species.dep = NULL) {
  
  require(dplyr)
  
  if(is.null(data))  stop("You need to define 'data'.")
  if(!ID.year   %in% colnames(data)) stop(paste0("Column '", ID.year, "' is missing in data."))
  if(!species   %in% colnames(data)) stop(paste0("Column '", species, "' is missing in data."))
  if(!abundance %in% colnames(data)) stop(paste0("Column '", abundance, "' is missing in data."))
  
  if(!is.null(ID.year.dep)) for(i in ID.year.dep) if(!i %in% colnames(data)) stop(paste0("Column '", i, "' is missing in data."))
  if(!is.null(species.dep)) for(i in species.dep) if(!i %in% colnames(data)) stop(paste0("Column '", i, "' is missing in data."))
  
  if(unique(is.na(data[,abundance])) != FALSE) warning(paste0("Column '", abundance, "' in data has NAs - these are replaced by zero abundances."))
  
  # create all possible ID.year-species combinations
  df.null <- as.data.frame(expand.grid("ID.year" = levels(data$ID.year),
                                       "species" = levels(data$species)))
  
  # add ID.year and species dependencies
  df.null <- left_join(df.null, unique(df.raw[, c(ID.year, ID.year.dep)]), by = ID.year)
  df.null <- left_join(df.null, unique(df.raw[, c(species, species.dep)]), by = species)
  
  # join with data
  df.zero   <- left_join(df.null, data, by = c(ID.year, ID.year.dep, species, species.dep))
  
  # replace NA abundances with 0
  df.zero[is.na(df.zero[abundance]),abundance] <- 0
  
  # get number of added zeros
  print(paste0(nrow(df.zero)-nrow(data), " zeros were added to data."))
  
  # get columns with NAs (all columns without dependencies)
  cols <- NULL
  cols <- colnames(data)[!colnames(data) %in% c(ID.year, ID.year.dep, species, species.dep, abundance)]
  
  if(!is.null(cols)) warning(paste0("Column '", cols, "' in data has NAs for all added zeros. Check whether it needs to be added to 'ID.year.dep' or 'species.dep' or add values manually.
                                    "))
  
  
  return(df.zero)

}


######################################################-
#### 5) View model convergence (Rhat, NEff, ...)  ####
######################################################-


## input:
#########-

## model.list   (default = NULL):              list with different models (of one species) for comparisons
## model.name   (default = names(model.list)): names of models (e.g. family and zero-inflated coefficients)
## response     (default = NULL):              name of response variable (character)
## plot.conv    (default = FALSE):             if TRUE, graphical output is given
## spec         (default = NULL):              name of species (for title of plot)

## output:
##########-

## the function mod.conv returns a dataframe including model convergence parameters, e.g.,
## maximum Rhat, minimum Neff/Ntot, maximum MCSE/SD, the number of divergent transitions and maximum treedepth

## The function optionally returns a graphical output if plot.stats = TRUE. 

mod.conv <- function(model.list = NULL, model.name = NULL, td = NULL,
                     plot.conv = FALSE, spec = NULL) {
  
  require(brms)
  require(tidyverse)
  require(bayestestR)

  if(is.null(model.list))  stop("You need to define the model-list.")
  if(!is.list(model.list)) stop("Model-list must be a list.")
  if(is.null(model.name)) model.name <- names(model.list)
  if(is.null(model.name)) model.name <- as.character(1:length(model.list)) # if the list is not named and still NULL
  if(is.null(td))          stop("You need to define the maximum treedepth 'td'.")
  
  df.modS.full <- NULL
  
  for (m in 1:length(model.list)) {
    
    df.modS  <- data.frame(conv = c("max. Rhat", "min. Neff/Ntot", "max. MCSE/SD", "div.trans", "max. td"), value = NA, Nthres = NA,
                           base = c(1, 1, 0, 0, 0),
                           threshold = c(1.01, 0.1, 0.1, 1, td),
                           color = "grey30")
    
    # rhat, neff, mcse
    df.modS$value[df.modS$conv == "max. Rhat"]  <-  max(rhat(model.list[[m]]))
    df.modS$Nthres[df.modS$conv == "max. Rhat"] <-  sum(rhat(model.list[[m]]) >= 1.01)
    
    df.modS$value[df.modS$conv == "min. Neff/Ntot"]   <- min(neff_ratio(model.list[[m]]))
    df.modS$Nthres[df.modS$conv == "min. Neff/Ntot"]  <- sum(neff_ratio(model.list[[m]]) <= 0.1)
    
    ps.df <- as.data.frame(apply(X = as.data.frame(model.list[[m]]), MARGIN = 2, FUN = sd))
    colnames(ps.df)[1] <- "sd"; ps.df$Parameter <- row.names(ps.df)
    mcse <- left_join(mcse(model.list[[m]], effects = "all", component = "all"), ps.df, by = "Parameter")
    
    df.modS$value[df.modS$conv == "max. MCSE/SD"]  <- max(mcse$MCSE / mcse$sd)
    df.modS$Nthres[df.modS$conv == "max. MCSE/SD"] <- sum((mcse$MCSE / mcse$sd) >= 0.1)
    
    # div. trans, treedepth
    nuts <- nuts_params(model.list[[m]])
    df.modS$value[df.modS$conv == "div.trans"]   <- max(nuts$Value[nuts$Parameter == "divergent__"])
    df.modS$value[df.modS$conv == "max. td"]    <- max(nuts$Value[nuts$Parameter == "treedepth__"])
    df.modS$Nthres[df.modS$conv == "div.trans"]   <- length(nuts$Value[nuts$Parameter == "divergent__" & nuts$Value > 1])
    df.modS$Nthres[df.modS$conv == "max. td"]    <- length(nuts$Value[nuts$Parameter == "treedepth__" & nuts$Value > td])
    df.modS$model          <- as.character(model.name[m])
    df.modS$species        <- spec
    
    df.modS.full <- rbind(df.modS.full, df.modS)
    
  }  
  
  df.modS.full$color[df.modS.full$Nthres >= 1] <- "orange" 
  
  if(plot.conv == TRUE) {
    require(ggplot2) 
    df.modS.full$conv <- factor(df.modS.full$conv, levels = c("max. Rhat", "min. Neff/Ntot", "max. MCSE/SD", "div.trans", "max. td"))
    
    print(ggplot(df.modS.full) +
            geom_hline(aes(yintercept = base), color = "grey20", lty = "solid") +
            geom_hline(aes(yintercept = threshold), color = "grey20", lty = "dashed") +
            geom_point(aes(x = model, y = value,), color = df.modS.full$color, size = 2) +
            ylab(paste0("value")) +
            xlab("model") +
            ggtitle(paste0(spec, " - model convergence"), 
                    subtitle = "As long as the values are between the solid and dashed line (grey dots), 
model convergence is fine. The closer to the solid line, the better.") +
            facet_wrap(~ conv, scales = "free_y") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  }
  
  return(select(df.modS.full, select = - color))
  
}
