###########################-
#### OEFS Paper Script ####
###########################-

## author: Mirjam R. Rieger
## latest update: 2023-06-16 (MR)
  # 2023-06-03 (MR): created script based on MA script
  # 2023-06-16 (MR): initial commit GitHub
  # 2023-07-03 (MR): add posterior predictions


#### 0) about this script ####
##############################-

  # siehe "Anleitung zur Brutvogel-Trendanalyse im Rahmen der  OEFS.csv"


#### 1) preparations ####
#########################-

rm(list = ls())

## used packages
library(tidyverse) # for dataset handling and graphs
library(brms)  # for models (brm())
#library(bayestestR) # for mcse()
library(pracma) # for factor loading rotation (PCA)
library(parallel) # detectCores()

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
csr <- require(cmdstanr) # for core/chain parallelisation, if not installed, chains will not be parallized
set_cmdstan_path(path = "~/cmdstan")

# define backend and cores based on whether cmdstanr is installed or not
cores <- detectCores() # get the number of cores for parallelization

if(csr == TRUE) {
  backend <- "cmdstanr"
  thread  <- floor((cores-1)/2)
  ncores  <- 2
}

if(csr == FALSE) {
  backend = "rstan"
  ncores <- min(cores-1, nc) 
  thread <- NULL
}


## load functions
prop_zero   <- function(z) sum(z == 0) / length(z)
source("./R/OEFS_Paper_Functions.R")

## define first (year1) and last year (yearx)
year1 <- 2002
yearx <- 2020

## define some model parameters
## --> do you want to do a (faster) test run, e.g. with fewer chains and iterations (TRUE) or the original run with (FALSE)?
test <- TRUE

if(test == TRUE) {
  nc <- 4           # no. of chains
  ni <- 2000        # no. of iterations (incl. warm-up)
  nw <- 500         # no. of iterations (warm-up only)
  a.delta <- 0.99   # adapt_delta value
  max.td  <- 12     # max_treedepth value
}

if(test == FALSE) {
  nc <- 8           # no. of chains
  ni <- 5500        # no. of iterations (incl. warm-up)
  nw <- 1500        # no. of iterations (warm-up only)
  a.delta <- 0.9999 # adapt_delta value
  max.td  <- 15     # max_treedepth value
}

## define species subset 
# in case you want to model only some species, e.g. c("Common Blackbird", "White Wagtail")
# if NULL, all available species will be used
Spec.sub <- NULL 

## read raw data files
df.raw  <- read.csv("./data/OEFS_Paper_RawData.csv",      sep = ",", encoding = "latin1")

## read additional files
df.PCA  <- read.csv("./data/OEFS_Paper_SiteInfo.csv",     sep = ",", encoding = "latin1")
df.land <- read.csv("./data/OEFS_Paper_LandscapeNRW.csv", sep = ",", encoding = "latin1")
df.spec <- read.csv("./data/OEFS_Paper_SpeciesInfo.csv",  sep = ",", encoding = "latin1")


#### 2) prepare raw data ####
#############################-

## adjust structure of dataframes
df.raw  <- df.raw  %>% mutate_at(c("ID", "region", "landscape", "species", "species_sci"), factor)
df.PCA  <- df.PCA  %>% mutate_at(c("ID"), factor)
df.land <- df.land %>% mutate_at(c("landscape", "region"), factor)
df.spec <- df.spec %>% mutate_at(c("species", "species_sci", "modR", "modZI", "modFAM"), factor)

## exclude years outside the range of year1-yearx
df.raw <- droplevels(subset(df.raw, year >= year1 & year <= yearx))

## add unique variable for each site and year (ID.year)
df.raw$ID.year <- as.factor(paste0(df.raw$ID, ".", df.raw$year))

## check number of observations per species (should be the same, otherwise zeros are missing)
table(df.raw$species)
# add zeros
dat <- add.zeros(data = df.raw,
                 ID.year.dep = c("ID", "year", "region", "landscape", "obs.ID"),
                 species.dep = c("species_sci")) # for details, see function-script
# check columns with NAs: replace NAs in 'exclude' with 'no'
dat$exclude[is.na(dat$exclude)] <- "no"

## exclude implausible values (exclude = "yes")
dat <- droplevels(subset(dat, exclude != "yes"))

## remove objects to free memory
rm(df.raw)

#### 2.1) observer effect ####
#############################-

## add total abundance per ID and year and mean total abundance per ID to delineate observer effects
## based on the default deviation of 25%
dat <- obs.eff(data = dat,
               ID = "ID",               # default column name
               abundance = "abundance", # default column name
               year = "year",           # default column name
               OE = 0.25,               # default threshold of 25% to delineate observer effect
               add.df = TRUE)           # default: returns an additional dataframe with observer effects per site and year

#### 2.2) further adaptions ####
###############################-

## round up .5 abundances
dat$ab.c <- ceiling(dat$abundance)

## scale year to year1 = 0 (for faster runtime)
dat$year.s <- dat$year-year1

## exclude species not present in df.spec 
# (here: 'others' which were only needed for total abundances to delineate observer effects)
for (sp in levels(dat$species)) {
  if (!(sp %in% df.spec$species)) {
    warning(paste0("Species '", sp, "' not present in df.spec. Species was excluded from modelling process."))
     
    ## remove missing species from dat
    dat <- droplevels(subset(dat, species != sp))
  }
}


#### 2.3) Artenliste erstellen ####
##################################-

spec.list <- sort(as.character(unique(dat$species)))

#### 3) PCA ####
################-
# this is done for all IDs which are included within the chosen time span to guarantee comparability between species

df.PCA <- left_join(summarize(dat %>% group_by(ID, year, region), .groups = "drop"), df.PCA, by = "ID")
df.PCA <- as.data.frame(df.PCA)

## run PCA
set.seed(3171); mod.PCA   <- prcomp(select(df.PCA, -c(ID, year, region)), scale = T)

## check summary
summary(mod.PCA)$importance

## plot variance explained by PCs
names(mod.PCA$sdev) <- as.character(1:length(mod.PCA$sdev)) # assures PCs are numbered in following graph
screeplot(mod.PCA, xlab = "Principal Component", main = "Variance of PCs"); abline(1,0)

## check scatter of PCs against each other (you might want to check more/other PCs)
par(mfrow = c(2,2))
plot(mod.PCA$x[,1], mod.PCA$x[,2],  xlab = "PC1", ylab = "PC2", pch = 16); abline(h = 0, v = 0, lty = 2)
plot(mod.PCA$x[,1], mod.PCA$x[,3],  xlab = "PC1", ylab = "PC3", pch = 16); abline(h = 0, v = 0, lty = 2)
plot(mod.PCA$x[,2], mod.PCA$x[,3],  xlab = "PC2", ylab = "PC3", pch = 16); abline(h = 0, v = 0, lty = 2)

## check correlation between original variables and PC scores
par(mfrow = c(1,1))
biplot(mod.PCA, choices = c(1, 2), col = c("grey90","black")); abline(h = 0, v = 0, lty = 2)
biplot(mod.PCA, choices = c(1, 3), col = c("grey90","black")); abline(h = 0, v = 0, lty = 2)
biplot(mod.PCA, choices = c(2, 3), col = c("grey90","black")); abline(h = 0, v = 0, lty = 2)

## factor rotation for differentiated PC meaning (for the first 3 PCs)
pca.raw <- mod.PCA$rotation[, 1:3] %*% diag(mod.PCA$sdev, 3, 3)
pca.rot <- varimax(pca.raw)$loadings
rot.df  <- as.data.frame(pca.rot[, 1:3])
rot.df$Covariate  <- rownames(rot.df)
pca.inverse       <- t(pinv(pca.raw))
pca.RotatedScores <- scale(select(df.PCA, -c(ID, year, region))) %*% pca.inverse

## Store scores in raw data frame:
df.PCA$PC1  <- mod.PCA$x[,1]          # unrotated Scores for PC1, PC2, PC3
df.PCA$PC2  <- mod.PCA$x[,2]
df.PCA$PC3  <- mod.PCA$x[,3]
df.PCA$PC1r <- pca.RotatedScores[,1]  # rotated Scores for PC1, PC2, PC3
df.PCA$PC2r <- pca.RotatedScores[,2]
df.PCA$PC3r <- pca.RotatedScores[,3]

## combine dat with PCA scores
dat2 <- left_join(dat, df.PCA, by = c("ID", "year", "region"))

## remove objects to free memory
rm(pca.inverse)
rm(pca.RotatedScores)
rm(pca.raw)
rm(mod.PCA)
rm(rot.df)

#### 4) weights ####
####################-
# this is done for each species separately since the exclusion of OK = 0 values is species-specific

#### 4.1) weights for both biogeographical regions combined ####
################################################################-
dat3 <- weight(df.habitat = df.land,  # landscape data
               habitat = "landscape", # column name for habitat
               area = "area_NRW",     # column name for area
               df.data = dat2,        # raw data
               ID = "ID",             # default column name
               year = "year",         # default column name
               by_spec = TRUE,        # calculate weight by species (e.g., if number of sites differs between species due to exclusion)
               species = "species",   # default column name
               by_reg = FALSE,        # default
               add.df = TRUE)         # default: returns an additional dataframe with weights per habitat and year (and species)

# rename "weight" to prevent overwriting in the next step
colnames(dat3)[colnames(dat3) == "weight"] <- "weight.l"

#### 4.2) weights for each biogeographical region separately ####
#################################################################-
## this is needed for species which are only modelled for one biogeographical region to guarantee a
## representative weighting of landscapes per region

dat4 <- weight(df.habitat = df.land,  # landscape data
               habitat = "landscape", # column name for habitat
               area = "area_NRW",     # column name for area
               df.data = dat3,        # raw data
               ID = "ID",             # default column name
               year = "year",         # default column name
               by_spec = TRUE,        # calculate weight by species (e.g., if number of sites differs between species due to exclusion)
               species = "species",   # default column name
               by_reg = TRUE,         # calculate weights separately per region
               region = "region",     # default column name
               add.df = TRUE)         # default: returns an additional dataframe with weights per habitat and year (and species)

# rename "weight" 
colnames(dat4)[colnames(dat4) == "weight"] <- "weight.l.r"

#### 5) Model ####
##################-

## ECXLUDE for testing only
###########################################-
df.spec[5, 1] <- "Common Chiffchaff"
df.spec[5, 3] <- "both"
df.spec[5, 4] <- "R"
df.spec[5, 5] <- "zip"
##########################################-

## create list of species which you want to model, either

# ... automatically based on whether it can be modelled or not
if (is.null(Spec.sub)) spec.list <- sort(as.character(unique(df.spec$species[!is.na(df.spec$modFAM)])))

# ... or based on the species subset which was defined in the beginning
if (!is.null(Spec.sub)) {
  spec.list <- Spec.sub
  no.mod <- NULL
  
  # check whether the selected species can be modelled
  for (sp in spec.list) {
    
    # is species present in dataset?
    if (!(sp %in% df.spec$species)) {
      warning(paste0("No modelling approach available for ", sp, "."))
      no.mod <- c(no.mod, sp)
    }
  }
  
  # exclude species which can't be modelled
  spec.list <- spec.list[which(c(!(spec.list %in% no.mod)))]

}



## for-loop for all species in spec.list

for (sp in spec.list) {
  
  ## get information from df.spec
  mod.bgr     <- unique(df.spec$modR[df.spec$species == sp])
  mod.fam.l   <- df.spec$modFAM[df.spec$species == sp]
  zi.type.l   <- df.spec$modZI[df.spec$species == sp] 
  
  ## subset dataframe and specify weight based on modelled biogeographical region(s) (mod.bgr)
  if(mod.bgr %in% c("both")) {
    S.dat <- droplevels(dat4[dat4$species == sp,])
    S.dat$weight <- S.dat$weight.l
  }
  
  if(mod.bgr %in% c("atl", "kon")) {
    S.dat <- droplevels(dat4[dat4$species == sp & dat4$region == as.character(mod.bgr),])
    S.dat$weight <- S.dat$weight.l.r
  }
  
  
  ## for-loop for different model per species (for comparison)
  for(m in 1:length(mod.fam.l)) {
    
    mod.fam   <- mod.fam.l[m]
    zi.type   <- zi.type.l[m]
    
    #### 5.1) define model formula
    ##############################-
    
    ## 5.1a) for poisson and negative binomial (without zero inflation)
    if(mod.fam %in% c("pois", "nb")) {
      
      if(zi.type == "none") {
        mod.form <- bf(ab.c|weights(weight, scale = FALSE) ~ 
                         OE_25 + 
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                         (1|ID))
      }
    }
    
    ## 5.1b) for zero-inflated poisson and zero-inflated negative binomial
    if(mod.fam %in% c("zip", "zinb")) {
      
      if(zi.type == "none") {
        mod.form <- bf(ab.c|weights(weight, scale = FALSE) ~ 
                         OE_25 + 
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                         (1|ID),
                       zi ~ 1)
      }
      
      if(zi.type == "R") {
        mod.form <- bf(ab.c|weights(weight, scale = FALSE) ~ 
                         OE_25 + 
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                         (1|ID),
                       zi ~ region)
      }
      
      if(zi.type == "PCA") {
        mod.form <- bf(ab.c|weights(weight, scale = FALSE) ~ 
                         OE_25 + 
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                         (1|ID),
                       zi ~ PC1r + PC2r + PC3r)
      }
      
      if(zi.type == "F") {
        mod.form <- bf(ab.c|weights(weight, scale = FALSE) ~ 
                         OE_25 + 
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                         (1|ID),
                       zi ~ (1|ID))
      }
      
      if(zi.type == "RF") {
        mod.form <- bf(ab.c|weights(weight, scale = FALSE) ~ 
                         OE_25 + 
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                         (1|ID),
                       zi ~ region + (1|ID))
      }
      
      if(zi.type == "PCAF") {
        mod.form <- bf(ab.c|weights(weight, scale = FALSE) ~ 
                         OE_25 + 
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                         (1|ID),
                       zi ~ PC1r + PC2r + PC3r + (1|ID))
      }
      
      if(zi.type == "PCARF") {
        mod.form <- bf(ab.c|weights(weight, scale = FALSE) ~ 
                         OE_25 + 
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                         (1|ID),
                       zi ~ region + PC1r + PC2r + PC3r + (1|ID))
      }
      
    }
    
    #### 5.2) update model formula (based on modelled biogeographical region(s))
    ############################################################################-
    
    if(mod.bgr %in% c("both")) {
      S.dat <- droplevels(dat4[dat4$species == sp,])
      S.dat$weight <- S.dat$weight.l
      # update model formula
      mod.form <- update(mod.form, ~ . + s(year.s, by = region) + region)
      
    }
    
    if(mod.bgr %in% c("atl", "kon")) {
      S.dat <- droplevels(dat4[dat4$species == sp & dat4$region == as.character(mod.bgr),])
      S.dat$weight <- S.dat$weight.l.r
      # update model formula
      mod.form <- update(mod.form, ~ . + s(year.s))
    }
    
    #### 5.3) define model family
    #############################-
    
    if(mod.fam == "pois") fam <- poisson(link = "log")
    if(mod.fam == "nb")   fam <- negbinomial(link = "log")
    if(mod.fam == "zip")  fam <- zero_inflated_poisson(link = "log")
    if(mod.fam == "zinb") fam <- zero_inflated_negbinomial(link = "log")
    
    
    #### 5.4) set priors
    ####################-
    
    priors <- get_prior(mod.form, 
                        data = S.dat, family = fam)
    priors$prior[1] <- "normal(0, 2.5)"
    
    
    #### 5.5) run model
    ###################-
    
    mod <- brm(mod.form, 
               data = S.dat, family = fam,
               cores = ncores, iter = ni, warmup = nw, chains = nc, 
               backend = backend, threads = threading(thread),
               thin = 1,
               control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
    
   
    ## save model
    saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type, mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
    
  } # end of for-loop (model)
} # end of for-loop (species)

  
## in fact it should work but I need to check the models on the 16 core compute. Afterwards add some model evaluations and graphs


for(sp in spec.list){
  
  mod.list <- list.files(path = "./01_models", pattern = sp)
  mods <- list()
  for(i in 1:length(mod.list))  mods[[i]]  <- readRDS(paste0("./01_models/", mod.list[i]))

  print(sp)
  print(mod.stat(model.list = mods, response = "ab.c", plot.stats = T))
  print(mod.conv(model.list = mods, td = max.td, plot.conv = T))
  
} 

