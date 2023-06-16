###########################-
#### OEFS Paper Script ####
###########################-

## author: Mirjam R. Rieger
## latest update: 2023-06-03 (MR)
  # 2023-06-03 (MR): created script based on MA script


#### 0) about this script ####
##############################-

  # siehe "Anleitung zur Brutvogel-Trendanalyse im Rahmen der  OEFS.csv"


#### 1) preparations ####
#########################-

rm(list = ls())

## set working directory
setwd("D:/Uni/Masterarbeit/Paper/Paper_Skript") ## exclude this later

## used packages
library(tidyverse) # for dataset handling and graphs
library(brms)  # for models (brm())
library(bayestestR) # for mcse()
library(pracma) # for factor loading rotation (PCA)

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
csr <- require(cmdstanr) # for core/chain parallelisation, if not installed, chains will not be parallized
#set_cmdstan_path(path = "~/cmdstan")

## load functions
prop_zero   <- function(z) sum(z == 0) / length(z)
source("./OEFS_Paper_Functions.R")

## define first (year1) and last year (yearx)
year1 <- 2002
yearx <- 2020

## define parameters for models
cores <- 16 # Anzahl Kerne/Threads des verwendeten Computers
nc <- 8     # Anzahl Ketten - nicht aendern
ni <- 5500  # Anzahl Wiederholungen (inkl. warm-up) - nicht aendern
nw <- 1500  # Anzahl Wiederholungen im warm-up - nicht aendern
a.delta <- 0.9999 # adapt_delta   - nicht aendern, außer R empfiehlt eine Erhoehung
max.td  <- 15   # max_treedepth - nicht aendern, außer R empfiehlt eine Erhoehung

## define species subset 
# in case you want to model only some species, e.g. c("Common Blackbird", "White Wagtail")
# if NULL, all available species will be used
Spec.sub <- NULL 

## read raw data files
df.raw       <- read.csv("../Paper_RawData/OEFS_Paper_RawData.csv", sep = ",", encoding = "latin1")

## read additional files
df.PCA  <- read.csv("../Paper_RawData/OEFS_Paper_SiteInfo.csv",     sep = ",", encoding = "latin1")
df.land <- read.csv("../Paper_RawData/OEFS_Paper_LandscapeNRW.csv", sep = ",", encoding = "latin1")
df.spec <- read.csv("../Paper_RawData/OEFS_Paper_SpeciesInfo.csv",  sep = ",", encoding = "latin1")


#### 2) prepare raw data ####
#############################-

## adjust df structure
df.raw  <- df.raw %>% mutate_at(c("ID", "region", "landscape", "species", "species_sci"), factor)
df.PCA  <- df.PCA %>% mutate_at(c("ID"), factor)
df.land <- df.land %>% mutate_at(c("landscape", "region"), factor)
df.spec <- df.spec %>% mutate_at(c("species", "species_sci", "modR", "modZI", "modFAM"), factor)

## exclude previous years
df.raw <- droplevels(subset(df.raw, year >= year1 & year <= yearx))

## add ID.year variable
df.raw$ID.year <- as.factor(paste0(df.raw$ID, ".", df.raw$year))

## check number of observations per species (should be the same, otherwise zeros are missing)
table(df.raw$species)
# add zeros
dat <- add.zeros(data = df.raw,
                 ID.year.dep = c("ID", "year", "region", "landscape", "obs.ID"),
                 species.dep = c("species_sci")) # for details, see function-script
# check columns with NAs: replace NAs in exclude with "no"
dat$exclude[is.na(dat$exclude)] <- "no"

## exclude implausible values (exclude = "yes")
dat <- droplevels(subset(dat, exclude != "yes"))

## remove objects to free memory
rm(df.raw)

#### 2a) observer effect ####
#############################-

## add total abundance per ID and year and mean total abundance per ID to delineate observer effects
## based on the default deviation of 25%
dat <- obs.eff(data = dat,
               OE = 0.25) # default threshold of 25% to delineate observer effect

#### 2b) further adaptions ####
###############################-

## round up .5 abundances
dat$ab.c <- ceiling(dat$abundance)

## scale year to year1 = 0 (for decreased runtime)
dat$year.s <- dat$year-year1

## exclude species not present in df.spec (here: 'others' which were only needed for total abundances to delineate observer effects)
for (sp in levels(dat$species)) {
  if (!(sp %in% df.spec$species)) {
    warning(paste0("Species '", sp, "' not present in df.spec. Species was excluded from modelling process."))
     
    ## remove missing species from dat
    dat <- droplevels(subset(dat, species != sp))
  }
}


#### 2c) Artenliste erstellen ####
##################################-

spec.list <- sort(as.character(unique(dat$species)))

#### 3) PCA ####
################################################################################-
# this is done for all IDs which are included within the chosen time span to guarantee comparability between species

df.PCA <- left_join(summarize(dat %>% group_by(ID, year, region), .groups = "drop"), df.PCA, by = "ID")
df.PCA <- as.data.frame(df.PCA)

## run PCA
set.seed(3171); mod.PCA   <- prcomp(select(df.PCA, -c(ID, year, region)), scale = T)
#saveRDS(mod.PCA, paste0("./02_Modell-Info/PCA_Modell_", year1, "-", yearx, ".RDS")) # save to plot model output in Startpage

## save summary
write.csv(as.data.frame(summary(mod.PCA)$importance), paste0("../Paper_output/PCA_summary_", year1, "-", yearx, ".csv"), 
          row.names = T, fileEncoding = "latin1")

## plot variance explained by PCs
names(mod.PCA$sdev) <- as.character(1:length(mod.PCA$sdev)) # assures PCs are numbered in following graph
screeplot(mod.PCA, xlab = "Principal Component", main = "Variance of PCs"); abline(1,0)

## check scatter of PCs against each other
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

hist(df.PCA$PC1r)
hist(df.PCA$PC2r)
hist(df.PCA$PC3r)

write.csv(df.PCA, paste0("../Paper_output/PCA_scores_"          , year1, "-", yearx, ".csv"), 
          row.names = F, fileEncoding = "latin1")
write.csv(rot.df, paste0("../Paper_output/PCA_rotated_loadings_", year1, "-", yearx, ".csv"), 
          row.names = F, fileEncoding = "latin1")

## combine dat with PCA scores
dat2 <- left_join(dat, df.PCA, by = c("ID", "year", "region"))

## remove objects to free memory
rm(pca.inverse)
rm(pca.RotatedScores)
rm(pca.raw)
rm(mod.PCA)
rm(rot.df)

#### 4) weights ####
################################################################################-
# this is done for each species separately since the exclusion of OK = 0 values is species-specific

## 4.1) weights for both biogeographical regions combined
#########################################################
dat3 <- weight(df.habitat = df.land, habitat = "landscape", area = "area_NRW",
                   df.data = dat2, ID = "ID", year = "year",
                   by_spec = TRUE, species = "species")

dat3$weight1 <- dat3$weight
dat3 <- select(dat3, -c(weight, habitat))

## 4.2) weights for each biogeographical region separately
###########################################################

######
###### CHECK whether this is really true? might be that this should be done separateley
######

dat3$landscape2    <- paste0(dat3$landscape, ".", dat3$region)
df.land$landscape2 <- paste0(df.land$landscape, ".", df.land$region)

dat4 <- weight(df.habitat = df.land, habitat = "landscape2", area = "area_NRW",
               df.data = dat3, ID = "ID", year = "year",
               by_spec = TRUE, species = "species")

dat <- dat4
#### 5) Model ####
################################################################################-


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
  
  # define model family & model
  mod.bgr   <- df.spec$modR[df.spec$species == sp]
  mod.fam   <- df.spec$modFAM[df.spec$species == sp]
  
  zi.type   <- df.spec$modZI[df.spec$species == sp] 
  
  
  # subset dataframe
  S.dat <- droplevels(dat[dat$species == sp,])
  S.dat.R <- NULL
  if(mod.bgr %in% c("atl", "kon")) {
    S.dat.R <- droplevels(S.dat[S.dat$region == as.character(mod.bgr),])
    
  }
  
  # define backend and cores based on whether cmdstanr is installed or not
  if(csr == TRUE) {
    backend <- "cmdstanr"
    thread  <- cores/2-1
    ncores  <- 2
  }
  
  if(csr == FALSE) {
    backend <- "rstan"
    ncores <- min(cores-1, nc) 
  }
  

  mod.form0 <- bf(ab.c|weights(weight1, scale = FALSE) ~ 
                   s(year.s, by = region) + 
                   Kart.Q25 + 
                   poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                   (1|ID) +
                   region,
                 zi ~ PC1r)
  
  mod.form <- bf(ab.c|weights(weight1, scale = FALSE) ~ 
                   s(year.s, by = region) + 
                   Kart.Q25 + 
                   poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                   (1|ID) +
                   region,
                 zi ~ PC1r,
                 family = poisson(link = "log"))
  
  mod.form <- update(mod.form, formula. = ~., family = zero_inflated_poisson(link = "log"))
  identical(mod.form0, mod.form)
  
  ############ update first formula works
  mod.form <- bf(ab.c|weights(weight1, scale = FALSE) ~ 
                    s(year.s, by = region) + 
                    Kart.Q25 + 
                    poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                    (1|ID),
                  zi ~ 1)
  
  mod.form <- update(mod.form, ~ . + region)
  identical(mod.form0, mod.form)
  
  ########### add family works
  mod.form <- bf(ab.c|weights(weight1, scale = FALSE) ~ 
                   s(year.s, by = region) + region + 
                   Kart.Q25 + 
                   poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                   (1|ID),
                 zi ~ 1) + zero_inflated_poisson(link = "log")
  
  mod.form1 <- mod.form0 + zero_inflated_poisson(link = "log")
  
  ## mod.form and mod.form1 do look the same but identical() says FALSE
  
  fam <- zero_inflated_poisson(link = "log")
  mod.form0 <- bf(ab.c|weights(weight1, scale = FALSE) ~ 
                    s(year.s, by = region) + 
                    Kart.Q25 + 
                    poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                    (1|ID) +
                    region,
                  zi ~ 1)
  
  
  priors <- get_prior(mod.form0, 
                      data = S.dat, family = zero_inflated_poisson(link = "log"))
  priors$prior[1] <- "normal(0, 2.5)"
  
  set.seed(3171)
  
  mod <- brm(mod.form0, 
             data = S.dat, family = fam,
             cores = ncores, iter = ni, warmup = nw, chains = nc, 
             backend = backend, threads = threading(thread),
             thin = 1,
             control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
  


  ## 5.1.1 modelling both regions (atl & kon)
  ######################################################################-
  
  if (mod.bgr == "both") {
    
    ## zero-inflated Poisson, zi ~ 1  
    ####################################################################-
    if (mod.fam == "zip" & zi.type == "none") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ 1), 
                          data = S.dat, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ 1), 
                 data = S.dat, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## zero-inflated Poisson, zi ~ region  
    ####################################################################-
    if (mod.fam == "zip" & zi.type == "R") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ region), 
                          data = S.dat, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ region), 
                 data = S.dat, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## zero-inflated Poisson, zi ~ PCs  
    ####################################################################- 
    if (mod.fam == "zip" & zi.type == "PCA") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ PC1r + PC2r + PC3r), 
                          data = S.dat, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ PC1r + PC2r + PC3r), 
                 data = S.dat, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)     
      
    }
    
    
    ## zero-inflated Poisson, zi ~ (1|ID)  
    ####################################################################-
    if (mod.fam == "zip" & zi.type == "F") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ (1|ID)), 
                          data = S.dat, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ (1|ID)), 
                 data = S.dat, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }      
  
    ## zero-inflated Poisson, zi ~ region + (1|ID)  
    ####################################################################-
    if (mod.fam == "zip" & zi.type == "RF") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ region + (1|ID)), 
                          data = S.dat, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ region + (1|ID)), 
                 data = S.dat, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    
    ## zero-inflated Poisson, zi ~ PCs + (1|ID)  
    ####################################################################- 
    if (mod.fam == "zip" & zi.type == "PCAF") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ PC1r + PC2r + PC3r + (1|ID)), 
                          data = S.dat, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ PC1r + PC2r + PC3r + (1|ID)), 
                 data = S.dat, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)     
      
    }
    
    ## zero-inflated Poisson, zi ~ region + PCs + (1|ID)  
    ####################################################################-    
    if (mod.fam == "zip" & zi.type == "PCARF") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ region + PC1r + PC2r + PC3r + (1|ID)), 
                          data = S.dat, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ region + PC1r + PC2r + PC3r + (1|ID)), 
                 data = S.dat, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## zero-inflated negative binomial, zi ~ 1  
    ####################################################################-
    if (mod.fam == "zinb" & zi.type == "none") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ 1), 
                          data = S.dat, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ 1), 
                 data = S.dat, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## zero-inflated negative binomial, zi ~ region  
    ####################################################################-
    if (mod.fam == "zinb" & zi.type == "R") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ region), 
                          data = S.dat, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ region), 
                 data = S.dat, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## zero-inflated negative binomial, zi ~ PCs  
    ####################################################################- 
    if (mod.fam == "zinb" & zi.type == "PCA") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ PC1r + PC2r + PC3r), 
                          data = S.dat, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ PC1r + PC2r + PC3r), 
                 data = S.dat, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)     
      
    }
    
    
    ## zero-inflated negative binomial, zi ~ (1|ID)  
    ####################################################################-
    if (mod.fam == "zinb" & zi.type == "F") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ (1|ID)), 
                          data = S.dat, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ (1|ID)), 
                 data = S.dat, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }      
    
    
    ## zero-inflated negative binomial, zi ~ region + (1|ID)  
    ####################################################################-
    if (mod.fam == "zinb" & zi.type == "RF") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ region + (1|ID)), 
                          data = S.dat, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ region + (1|ID)), 
                 data = S.dat, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## zero-inflated negative binomial, zi ~ PCs + (1|ID)  
    ####################################################################- 
    if (mod.fam == "zinb" & zi.type == "PCAF") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ PC1r + PC2r + PC3r + (1|ID)), 
                          data = S.dat, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ PC1r + PC2r + PC3r + (1|ID)), 
                 data = S.dat, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)     
      
    }
    
    ## zero-inflated negative binomial, zi ~ region + PCs + (1|ID)  
    ####################################################################-    
    if (mod.fam == "zinb" & zi.type == "PCARF") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ region + PC1r + PC2r + PC3r + (1|ID)), 
                          data = S.dat, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ region + PC1r + PC2r + PC3r + (1|ID)), 
                 data = S.dat, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## Poisson without zero-inflation 
    ####################################################################-
    if (mod.fam =="pois"  & zi.type == "none") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID)), 
                          data = S.dat, family = poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID)), 
                 data = S.dat, family = poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)    
      
    }
    
    ## negative binomial without zero-inflation 
    ####################################################################-
    if (mod.fam =="nb"  & zi.type == "none") {
      
      priors <- get_prior(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                               s(year.s, by = region) + region + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID)), 
                          data = S.dat, family = negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight1, scale = FALSE) ~ 
                      s(year.s, by = region) + region +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID)), 
                 data = S.dat, family = negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)    
      
    }  
  }
  
  ## 5.1.2 modelling only one region (atl OR kon)
  ######################################################################-
  
  if (mod.bgr %in% c("atl", "kon")) {
    
    ## zero-inflated Poisson, zi ~ 1  
    ####################################################################-
    if (mod.fam == "zip" & zi.type == "none") {
      
      priors <- get_prior(bf(ab.c|weights(weight, scale = FALSE) ~ 
                               s(year.s) +
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ 1), 
                          data = S.dat.R, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight, scale = FALSE) ~ 
                      s(year.s) +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ 1),
                 data = S.dat.R, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## zero-inflated Poisson, zi ~ PCs  
    ##################################################################- 
    if (mod.fam == "zip" & zi.type == "PCA") {
      
      priors <- get_prior(bf(ab.c|weights(weight, scale = FALSE) ~ 
                               s(year.s) +
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ PC1r + PC2r + PC3r), 
                          data = S.dat.R, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight, scale = FALSE) ~ 
                      s(year.s) +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ PC1r + PC2r + PC3r), 
                 data = S.dat.R, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)     
      
    }
    
    ## zero-inflated Poisson, zi ~ (1|ID)  
    ####################################################################-
    if (mod.fam == "zip" & zi.type == "F") {
      
      priors <- get_prior(bf(ab.c|weights(weight, scale = FALSE) ~ 
                               s(year.s) +
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ (1|ID)), 
                          data = S.dat.R, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight, scale = FALSE) ~ 
                      s(year.s) +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ (1|ID)),
                 data = S.dat.R, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## zero-inflated Poisson, zi ~ PCs + (1|ID)  
    ##################################################################- 
    if (mod.fam == "zip" & zi.type == "PCAF") {
      
      priors <- get_prior(bf(ab.c|weights(weight, scale = FALSE) ~ 
                               s(year.s) +
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ PC1r + PC2r + PC3r + (1|ID)), 
                          data = S.dat.R, family = zero_inflated_poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight, scale = FALSE) ~ 
                      s(year.s) +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ PC1r + PC2r + PC3r + (1|ID)), 
                 data = S.dat.R, family = zero_inflated_poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)     
      
    }
    
    ## zero-inflated negative binomial, zi ~ 1  
    #####################################################################-
    if (mod.fam == "zinb" & zi.type == "none") {
      
      priors <- get_prior(bf(ab.c|weights(weight, scale = FALSE) ~ 
                               s(year.s) + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ 1), 
                          data = S.dat.R, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight, scale = FALSE) ~ 
                      s(year.s) +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ 1), 
                 data = S.dat.R, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## zero-inflated negative binomial, zi ~ PCs  
    ##################################################################- 
    if (mod.fam == "zinb" & zi.type == "PCA") {
      
      priors <- get_prior(bf(ab.c|weights(weight, scale = FALSE) ~ 
                               s(year.s) + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ PC1r + PC2r + PC3r), 
                          data = S.dat.R, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight, scale = FALSE) ~ 
                      s(year.s) +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ PC1r + PC2r + PC3r), 
                 data = S.dat.R, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)     
      
    }
    
    
    ## zero-inflated negative binomial, zi ~ (1|ID)  
    #####################################################################-
    if (mod.fam == "zinb" & zi.type == "F") {
      
      priors <- get_prior(bf(ab.c|weights(weight, scale = FALSE) ~ 
                               s(year.s) + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ (1|ID)), 
                          data = S.dat.R, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight, scale = FALSE) ~ 
                      s(year.s) +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ (1|ID)), 
                 data = S.dat.R, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors) 
      
    }
    
    ## zero-inflated negative binomial, zi ~ PCs + (1|ID)  
    ##################################################################- 
    if (mod.fam == "zinb" & zi.type == "PCAF") {
      
      priors <- get_prior(bf(ab.c|weights(weight, scale = FALSE) ~ 
                               s(year.s) + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID),
                             zi ~ PC1r + PC2r + PC3r + (1|ID)), 
                          data = S.dat.R, family = zero_inflated_negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight, scale = FALSE) ~ 
                      s(year.s) +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID),
                    zi ~ PC1r + PC2r + PC3r + (1|ID)), 
                 data = S.dat.R, family = zero_inflated_negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)     
      
    }
    
    ## Poisson without zero-inflation 
    #####################################################-
    if (mod.fam =="pois"  & zi.type == "none") {
      
      priors <- get_prior(bf(ab.c|weights(weight, scale = FALSE) ~ 
                               s(year.s) +
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID)), 
                          data = S.dat.R, family = poisson(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight, scale = FALSE) ~ 
                      s(year.s) +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID)), 
                 data = S.dat.R, family = poisson(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)    
      
    }
    
    ## negative binomial without zero-inflation 
    #####################################################-
    if (mod.fam =="nb"  & zi.type == "none") {
      
      priors <- get_prior(bf(ab.c|weights(weight, scale = FALSE) ~ 
                               s(year.s) + 
                               Kart.Q25 + 
                               poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                               (1|ID)), 
                          data = S.dat.R, family = negbinomial(link = "log"))
      priors$prior[1] <- "normal(0, 2.5)"
      
      set.seed(3171)
      mod <- brm(bf(ab.c|weights(weight, scale = FALSE) ~ 
                      s(year.s) +
                      Kart.Q25 + 
                      poly( PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) + 
                      (1|ID)), 
                 data = S.dat.R, family = negbinomial(link = "log"),
                 cores = ncores, iter = ni, warmup = nw, chains = nc, 
                 backend = backend, threads = threading(thread),
                 thin = 1,
                 control = list(adapt_delta = a.delta, max_treedepth = max.td), prior = priors)    
      
    }  
    
  }
  
  ## save model
  #############
  saveRDS(mod, paste0("./01_Modelle/OEFS_Modell_", sp, "_", zi.type, mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))

}



