##################-
#### EAS bird ####
##################-


#### 1) preparations ####
#########################-

rm(list = ls())

## used packages
pckgs <- c(
    "tidyverse", # for dataset handling and graphs
    "brms", # for models (brm())
    "bayestestR", # for mcse()
    "pracma", # for factor loading rotation (PCA)
    "parallel", # detectCores()
    "ggbeeswarm", # for plots
    "svDialogs", # for pop-up window
    "future" # for parallelization of kfold
)
for (i in pckgs) {
    if (!i %in% installed.packages()) {
        install.packages(i, dependencies = TRUE)
    }
}
sapply(pckgs, require, character.only = TRUE)


csr <- require(cmdstanr) # for core/chain parallelisation, if not installed, chains will not be parallelized and it will take more time
# check out readme for instructions on how to use cmdstanr and cmdstan if you want to use it
if(csr) cs <- try(cmdstan_version(), silent = TRUE) # cmdstan installed?
if("try-error" %in% class(cs)) {
  csr <- FALSE
}
if(csr) check_cmdstan_toolchain(fix = TRUE) # toolchain set up properly?

if(csr) {
  var <- dlgInput("Are you sure that you want to use cmdstanr for parallelized evaluation? 
                  Otherwise, rstan will be used instead [Yes/No]", 
                  default = "")$res
  if(!(var %in% c("Yes", "yes", "Y", "y"))) csr <- FALSE 
}


## load functions and EAS package
prop_zero   <- function(z) sum(z == 0) / length(z)
if (!"EAS" %in% installed.packages()) {
  devtools::install_github("m-rieger/EAS")
}
library(EAS)

## define first (year1) and last year (yearx)
year1 <- 2002
yearx <- 2020

## define baseline year (for index)
yearI <- 2006

## define time range for longterm trend and last year or the range
period <- 12
yearP  <- yearx-1 # here: 2019

## define some model parameters
## --> do you want to do a (faster) test run, e.g. with fewer chains and iterations (TRUE) or the original run with (FALSE)?
test <- TRUE

if (test == TRUE) {
  nc <- 4           # no. of chains
  ni <- 3000        # no. of iterations (incl. warm-up)
  nw <- 1000        # no. of iterations (warm-up only)
  a.delta <- 0.95   # adapt_delta value
  max.td  <- 10     # max_treedepth value
}

if (test == FALSE) {
  nc <- 8           # no. of chains
  ni <- 5500        # no. of iterations (incl. warm-up)
  nw <- 1500        # no. of iterations (warm-up only)
  a.delta <- 0.9999 # adapt_delta value
  max.td  <- 15     # max_treedepth value
}


# define backend and cores based on whether cmdstanr is installed or not
cores <- detectCores() # get the number of cores for parallelization

if (csr == TRUE) {
  backend <- "cmdstanr"
  ncores  <- 4
  thread  <- floor(cores/ncores)
}

if (csr == FALSE) {
  backend <- "rstan"
  ncores  <- min(cores, nc)
  thread  <- NULL
}

## define species subset
# in case you want to model only some species, e.g. c("Common Blackbird", "White Wagtail")
# if NULL, all available species will be used
Spec.sub <- c("Common Blackbird", # nb
             "Common Kestrel", # pois
             "Eurasian Jay", # pois
             #"Heckenbraunelle", # zinb F
             "Eurasian Blackcap", # zinb none
             "Eurasian Magpie", # zinb PCA
             "Carrion Crow", # zinb PCAF
             "Marsh Tit", # zinb PCARF
             "Common Chiffchaff",  # zinb R
             # "Zaunkönig", # zinb R
             "Great Tit", # zinb RF
             #"Nachtigall", # zip F
             "Common Chaffinch", # zip none
             #"Gartenrotschwanz", # zip PCA
             "White Wagtail", # zip PCAF
             "Great Spotted Woodpecker", # zip PCARF
             "Willow Tit", # zip R
             "Common House Martin") # zip RF

## read raw data files
df.raw  <- read.csv("./data/EAS_bird_RawData.csv",      sep = ",", encoding = "latin1")

## read additional files
# landuse and climate variables
df.PCA  <- read.csv("./data/EAS_bird_SiteInfo.csv",     sep = ",", encoding = "latin1")
# landscape shares
df.land <- read.csv("./data/EAS_bird_LandscapeNRW.csv", sep = ",", encoding = "latin1")
# species information (model type, family)
df.spec <- read.csv("./data/EAS_bird_SpeciesInfo.csv",  sep = ",", encoding = "latin1")

## create folders for models and output
if(!dir.exists("01_models")) dir.create("01_models")
if(!dir.exists("02_output")) dir.create("02_output")

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
                 ID.year.dep = c("ID", "year", "region", "landscape", "obs.ID", "habitat", "metro"),
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

## use 2*abundance with offset 2
dat$ab2 <- 2*dat$abundance
dat$off <- 2

## add natural region
dat$nat.region <- dat$habitat
dat$nat.region[dat$metro == "yes"] <- "metro"
dat$nat.region <- as.factor(dat$nat.region)

## scale year to year1 = 0 (for faster runtime)
dat$year.s <- dat$year - year1

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
set.seed(3171); mod.PCA <- prcomp(select(df.PCA, -c(ID, year, region)), scale = TRUE)

## check summary
summary(mod.PCA)$importance

## plot variance explained by PCs
names(mod.PCA$sdev) <- as.character(1:length(mod.PCA$sdev)) # assures PCs are numbered in following graph
screeplot(mod.PCA, xlab = "Principal Component", main = "Variance of PCs"); abline(1,0)

## check scatter of PCs against each other (you might want to check more/other PCs)
plot(mod.PCA$x[,1], mod.PCA$x[,2],  xlab = "PC1", ylab = "PC2", pch = 16); abline(h = 0, v = 0, lty = 2)
plot(mod.PCA$x[,1], mod.PCA$x[,3],  xlab = "PC1", ylab = "PC3", pch = 16); abline(h = 0, v = 0, lty = 2)
plot(mod.PCA$x[,2], mod.PCA$x[,3],  xlab = "PC2", ylab = "PC3", pch = 16); abline(h = 0, v = 0, lty = 2)

## check correlation between original variables and PC scores
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
dat <- left_join(dat, df.PCA, by = c("ID", "year", "region"))

## remove objects to free memory
rm(pca.inverse)
rm(pca.RotatedScores)
rm(pca.raw)
rm(mod.PCA)
rm(rot.df)


###### calculation of nat.region composition
df.land$nat.region <- df.land$habitat
df.land$nat.region[df.land$metro == "yes"] <- "metro"

df.land.ov <- df.land %>% group_by(nat.region) %>%
  summarize(area_ov = sum(area_NRW), .groups = "drop")
df.land.ov$prop.area <- round(df.land.ov$area_ov/sum(df.land.ov$area_ov), 4)

df.land.bgr <- df.land %>% group_by(nat.region, region) %>%
  summarize(area_ov = sum(area_NRW), .groups = "drop") %>%
  group_by(region) %>%
  mutate(area_bgr = sum(area_ov)) %>%
  ungroup()

df.land.bgr$prop.area <- round(df.land.bgr$area_ov/df.land.bgr$area_bgr, 4)


#### 4) weights ####
####################-
# this is done for each species separately since the exclusion of OK = 0 values is species-specific

# #### 4.1) weights for both biogeographical regions combined ####
# ################################################################-
# dat <- weight(df.habitat = df.land,   # landscape data
#                habitat = "landscape", # column name for habitat
#                area = "area_NRW",     # column name for area
#                df.data = dat,         # raw data
#                ID = "ID",             # default column name
#                year = "year",         # default column name
#                by_spec = TRUE,        # calculate weight by species (e.g., if number of sites differs between species due to exclusion)
#                species = "species",   # default column name
#                by_reg = FALSE,        # default
#                add.df = TRUE)         # default: returns an additional dataframe with weights per habitat and year (and species)
# # remark: here, KB.yes and SB.yes are not present in the dataset
# 
# # rename "weight" to prevent overwriting in the next step
# colnames(dat)[colnames(dat) == "weight"] <- "weight.l"
# 
# #### 4.2) weights for each biogeographical region separately ####
# #################################################################-
# ## this is needed for species which are only modelled for one biogeographical region to guarantee a
# ## representative weighting of landscapes per region
# 
# dat <- weight(df.habitat = df.land,  # landscape data
#                habitat = "landscape", # column name for habitat
#                area = "area_NRW",     # column name for area
#                df.data = dat,        # raw data
#                ID = "ID",             # default column name
#                year = "year",         # default column name
#                by_spec = TRUE,        # calculate weight by species (e.g., if number of sites differs between species due to exclusion)
#                species = "species",   # default column name
#                by_reg = TRUE,         # calculate weights separately per region
#                region = "region",     # default column name
#                add.df = TRUE)         # default: returns an additional dataframe with weights per habitat and year (and species)
# # remark: here, several landscape combinations are not present in 'atl' and 'kon'
# 
# # rename "weight"
# colnames(dat)[colnames(dat) == "weight"] <- "weight.l.r"

#### 5) Model ####
##################-

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


## prior fun
prior.fun <- function(mod.form, data, family) {
  priors <- get_prior(mod.form, data = S.dat, family = fam)
  priors$prior[priors$dpar == "" & priors$class == "b" & priors$coef != ""] <- "normal(0, 2.5)"
  priors$prior[priors$coef %in% c("polyPC1r21", "polyPC1r22", 
                                  "polyPC2r21", "polyPC2r22", 
                                  "polyPC3r21", "polyPC3r22")] <- "normal(0, 10)"
  priors$prior[priors$prior == "" & priors$dpar == "zi" & priors$class == "b" & priors$coef != ""] <- "normal(0, 2.5)"
  priors$prior[priors$coef %in% c("PC1r", "PC2r", "PC3r")] <- "normal(0, 10)"
  
  return(priors)
}

#### 5.1) define model formula
##############################-
mod.form <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                 OE_25 +
                 offset(log(off)) +
                 poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                 (1|ID))

mod.form.n <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                   OE_25 +
                   offset(log(off)) +
                   poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                   (1|ID),
                 zi ~ 1)

mod.form.R <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                   OE_25 +
                   offset(log(off)) +
                   poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                   (1|ID),
                 zi ~ nat.region)

mod.form.PC <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                    OE_25 +
                    offset(log(off)) +
                    poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                    (1|ID),
                  zi ~ PC1r + PC2r + PC3r)

mod.form.F <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                   OE_25 +
                   offset(log(off)) +
                   poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                   (1|ID),
                 zi ~ (1|ID))

mod.form.RF <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                    OE_25 +
                    offset(log(off)) +
                    poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                    (1|ID),
                  zi ~ nat.region + (1|ID))

mod.form.PCF <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                     OE_25 +
                     offset(log(off)) +
                     poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                     (1|ID),
                   zi ~ PC1r + PC2r + PC3r + (1|ID))

mod.form.PCRF <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                      OE_25 +
                      offset(log(off)) +
                      poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                      (1|ID),
                    zi ~ nat.region + PC1r + PC2r + PC3r + (1|ID))


## set global model options
options(
  brms.backend = backend, 
  brms.threads = threading(thread),
  control = list(adapt_delta = a.delta, max_treedepth = max.td),
  cmdstanr_write_stan_file_dir = "./01_models/stan_dir"
)

## for-loop for all species in spec.list
# spec.list <- spec.list[50:61]
## base priors
p <- pF <- c(set_prior("normal(0, 2.5)", class = "b", 
                       coef = c("nat.regionB", "nat.regionKB", "nat.regionKM", 
                                "nat.regionmetro", "nat.regionSB", "nat.regionST", 
                                "OE_25none", "OE_25positive", 
                                "syear.s:nat.regionA_1", "syear.s:nat.regionB_1", 
                                "syear.s:nat.regionKB_1", "syear.s:nat.regionKM_1", 
                                "syear.s:nat.regionmetro_1", "syear.s:nat.regionSB_1", 
                                "syear.s:nat.regionST_1")),
             set_prior("normal(0, 10)", class = "b", 
                       coef = c("polyPC1r21", "polyPC2r21", "polyPC3r21", 
                                "polyPC1r22", "polyPC2r22", "polyPC3r22")))
pR <- pRF <- c(p, set_prior("normal(0, 2.5)", class = "b", dpar = "zi",
                            coef = c("nat.regionB", "nat.regionKB", "nat.regionKM", 
                                     "nat.regionmetro", "nat.regionSB", "nat.regionST")))
pPC <- pPCF <- c(p, set_prior("normal(0, 10)", class = "b", dpar = "zi",
                              coef = c("PC1r", "PC2r", "PC3r")))
pPCR <- pPCRF <- c(pR, set_prior("normal(0, 10)", class = "b", dpar = "zi",
                                 coef = c("PC1r", "PC2r", "PC3r")))

spec.listB <- spec.list[spec.list %in% df.spec$species[df.spec$modR == "both"]]
spec.listA <- spec.list[spec.list %in% df.spec$species[df.spec$modR == "atl"]]
spec.listC <- spec.list[spec.list %in% df.spec$species[df.spec$modR == "con"]]

## models for region both ------------------------------------------------------
##------------------------------------------------------------------------------
mod.bgr <- "both"

## pois none -------------------------------------------------------------------
mod.fam <- "pois"; zi.type <- "none"; fam <- poisson(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = p) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## nb none ---------------------------------------------------------------------
mod.fam <- "nb"; zi.type <- "none"; fam <- negbinomial(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = p) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zip none --------------------------------------------------------------------
mod.fam <- "zip"; zi.type <- "none"; fam <- zero_inflated_poisson(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.n,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = p) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zip R -----------------------------------------------------------------------
mod.fam <- "zip"; zi.type <- "R"; fam <- zero_inflated_poisson(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.R,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pR) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zip PCA ---------------------------------------------------------------------
mod.fam <- "zip"; zi.type <- "PCA"; fam <- zero_inflated_poisson(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.PC,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pPC) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zip F -----------------------------------------------------------------------
mod.fam <- "zip"; zi.type <- "F"; fam <- zero_inflated_poisson(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.F,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pF) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zip RF ----------------------------------------------------------------------
mod.fam <- "zip"; zi.type <- "RF"; fam <- zero_inflated_poisson(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.RF,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pRF) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zip PCAF --------------------------------------------------------------------
mod.fam <- "zip"; zi.type <- "PCAF"; fam <- zero_inflated_poisson(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.PCF,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pPCF) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zip PCARF -------------------------------------------------------------------
mod.fam <- "zip"; zi.type <- "PCARF"; fam <- zero_inflated_poisson(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.PCRF,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pPCRF) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zinb none -------------------------------------------------------------------
mod.fam <- "zinb"; zi.type <- "none"; fam <- zero_inflated_negbinomial(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.n,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = p) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zinb R ----------------------------------------------------------------------
mod.fam <- "zinb"; zi.type <- "R"; fam <- zero_inflated_negbinomial(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.R,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pR) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zinb PCA --------------------------------------------------------------------
mod.fam <- "zinb"; zi.type <- "PCA"; fam <- zero_inflated_negbinomial(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.PC,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pPC) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zinb F ----------------------------------------------------------------------
mod.fam <- "zinb"; zi.type <- "F"; fam <- zero_inflated_negbinomial(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.F,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pF) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zinb RF ---------------------------------------------------------------------
mod.fam <- "zinb"; zi.type <- "RF"; fam <- zero_inflated_negbinomial(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.RF,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pRF) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zinb PCAF -------------------------------------------------------------------
mod.fam <- "zinb"; zi.type <- "PCAF"; fam <- zero_inflated_negbinomial(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.PCF,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pPCF) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

## zinb PCARF ------------------------------------------------------------------
mod.fam <- "zinb"; zi.type <- "PCARF"; fam <- zero_inflated_negbinomial(link = "log")

for (sp in spec.listB) {
  
  mod <- brm(mod.form.PCRF,
             cores = ncores, 
             iter = ni,
             warmup = nw, 
             chains = nc,
             data = droplevels(dat[dat$species == sp,]), 
             family = fam, prior = pPCRF) 
  stanfit <- mod$fit
  
  plan(multisession); mc.cores = parallel::detectCores()
  mod <- add_criterion(mod, "kfold", K = 16, chains = 1, threads = threading(1))
  plan(sequential)
  
  saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, 
                      "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  
}

# for (sp in spec.list) {
# 
#   ## get information from df.spec
#   mod.bgr   <- unique(df.spec$modR[df.spec$species == sp])
# 
#   ## subset dataframe and specify weight based on modelled biogeographical region(s) (mod.bgr)
#   if (mod.bgr %in% c("both")) {
#     S.dat <- droplevels(dat[dat$species == sp,])
#     S.dat$weight <- S.dat$weight.l
#   }
# 
#   if (mod.bgr %in% c("atl", "kon")) {
#     S.dat <- droplevels(dat[dat$species == sp & dat$region == as.character(mod.bgr),])
#     S.dat$weight <- S.dat$weight.l.r
#   }
# 
#   
#   ## pois none
#   mod.fam <- "pois"; zi.type <- "none"; fam <- poisson(link = "log")
#   priors <- prior.fun(mod.form, S.dat, fam)
#   
#   mod <- brm(mod.form,
#              cores = ncores, 
#              iter = ni,
#              warmup = nw, 
#              chains = nc,
#              data = S.dat, family = fam, prior = priors,
#              cmdstanr_write_stan_file_dir = "./01_models/stan__dir") 
#   
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
# 
#   ## nb none  
#   mod.fam <- "nb"; fam <- negbinomial(link = "log")
#   priors <- prior.fun(mod.form, S.dat, fam)
#   mod <- update(mod, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zip none
#   mod.fam <- "zip"; fam <- zero_inflated_poisson(link = "log")
#   priors <- prior.fun(mod.form.n, S.dat, fam)
#   mod <- update(mod, mod.form.n, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zip R
#   zi.type <- "R"
#   priors <- prior.fun(mod.form.R, S.dat, fam)
#   mod <- update(mod, mod.form.R, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zip PC
#   zi.type <- "PCA"
#   priors <- prior.fun(mod.form.PC, S.dat, fam)
#   mod <- update(mod, mod.form.PC, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zip F
#   zi.type <- "F"
#   priors <- prior.fun(mod.form.F, S.dat, fam)
#   mod <- update(mod, mod.form.F, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zip RF
#   zi.type <- "RF"
#   priors <- prior.fun(mod.form.RF, S.dat, fam)
#   mod <- update(mod, mod.form.RF, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zip PCF
#   zi.type <- "PCAF"
#   priors <- prior.fun(mod.form.PCF, S.dat, fam)
#   mod <- update(mod, mod.form.PCF, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zip PCRF
#   zi.type <- "PCARF"
#   priors <- prior.fun(mod.form.PCRF, S.dat, fam)
#   mod <- update(mod, mod.form.PCRF, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zinb none
#   zi.type <- "none"
#   mod.fam <- "zinb"; fam <- zero_inflated_negbinomial(link = "log")
#   priors <- prior.fun(mod.form.n, S.dat, fam)
#   mod <- update(mod, mod.form.n, family = fam, prior = priors, cores = ncores,
#                 file = paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   # saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zinb R
#   zi.type <- "R"
#   priors <- prior.fun(mod.form.R, S.dat, fam)
#   mod <- update(mod, mod.form.R, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zinb PC
#   zi.type <- "PCA"
#   priors <- prior.fun(mod.form.PC, S.dat, fam)
#   mod <- update(mod, mod.form.PC, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zinb F
#   zi.type <- "F"
#   priors <- prior.fun(mod.form.F, S.dat, fam)
#   mod <- update(mod, mod.form.F, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zinb RF
#   zi.type <- "RF"
#   priors <- prior.fun(mod.form.RF, S.dat, fam)
#   mod <- update(mod, mod.form.RF, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zinb PCF
#   zi.type <- "PCAF"
#   priors <- prior.fun(mod.form.PCF, S.dat, fam)
#   mod <- update(mod, mod.form.PCF, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
#   ## zinb PCRF
#   zi.type <- "PCARF"
#   priors <- prior.fun(mod.form.PCRF, S.dat, fam)
#   mod <- update(mod, mod.form.PCRF, family = fam, prior = priors, cores = ncores)
#   saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
#   
# } # end of for-loop (species)
# 

#### 6) posterior predictive checks and model convergence ####
##############################################################-

## check model convergence and posterior predictions of response variable
## to chose the best fitting model structure

df.stat  <- data.frame()
df.conv  <- data.frame()
df.waic  <- data.frame()
df.loo   <- data.frame()
df.kfold <- data.frame()

for (sp in spec.list){

  mod.list <- list.files(path = "./01_models", pattern = sp)
  mods <- list()
  for (i in 1:length(mod.list))  {
    m <- unlist(strsplit(mod.list[i], split = "_"))[3]
    mods[[m]]  <- readRDS(paste0("./01_models/", mod.list[i]))
    # print(mod.list[i])
    # print(mods[[m]]$formula)
    # print(mods[[m]]$family)
  }

  # ## posterior predictions of response variable
  # jpeg(paste0("./01_models/mod_selection/", sp, "_stats.jpeg"), width = 700, height = 400)
  # tmp.stat <- mod.stat(model.list = mods, model.name = names(mods),
  #                response = "ab2",
  #                plot.stats = TRUE,
  #                spec = sp)
  # dev.off()
  # 
  # df.stat <- rbind(df.stat, tmp.stat)
  # 
  # ## model convergence
  # jpeg(paste0("./01_models/mod_selection/", sp, "_conv.jpeg"), width = 700, height = 400)
  # tmp.conv <- mod.conv(model.list = mods, model.name = names(mods),
  #                td = max.td,
  #                plot.conv = TRUE,
  #                spec = sp)
  # dev.off()
  # 
  # df.conv <- rbind(df.conv, tmp.conv)
  # 
  # tmp.waic <- as.data.frame(waic(mods[["nonepois"]], mods[["nonenb"]], 
  #                  mods[["nonezip"]], mods[["nonezinb"]],
  #                  mods[["Rzip"]], mods[["Rzinb"]],
  #                  mods[["PCAzip"]], mods[["PCAzinb"]],
  #                  mods[["Fzip"]], mods[["Fzinb"]],
  #                  mods[["RFzip"]], mods[["RFzinb"]],
  #                  mods[["PCAFzip"]], mods[["PCAFzinb"]],
  #                  mods[["PCARFzip"]], mods[["PCARFzinb"]])[["diffs"]])
  # 
  # tmp.waic$model <- rownames(tmp.waic)
  # tmp.waic$model <- gsub('.*\\\"([^\"]+)\\\".*', '\\1', tmp.waic$model)
  # tmp.waic$species <- sp
  # 
  # df.waic <- rbind(df.waic, tmp.waic)
  
  tmp.loo <- as.data.frame(loo(mods[["nonepois"]], mods[["nonenb"]], 
                                 mods[["nonezip"]], mods[["nonezinb"]],
                                 mods[["Rzip"]], mods[["Rzinb"]],
                                 mods[["PCAzip"]], mods[["PCAzinb"]],
                                 mods[["Fzip"]], mods[["Fzinb"]],
                                 mods[["RFzip"]], mods[["RFzinb"]],
                                 mods[["PCAFzip"]], mods[["PCAFzinb"]],
                                 mods[["PCARFzip"]], mods[["PCARFzinb"]])[["diffs"]])
  
  tmp.loo$model <- rownames(tmp.loo)
  tmp.loo$model <- gsub('.*\\\"([^\"]+)\\\".*', '\\1', tmp.loo$model)
  tmp.loo$species <- sp
  
  df.loo <- rbind(df.loo, tmp.loo)
  
  
  # tmp.kfold <- as.data.frame(loo_compare(mods[["nonepois"]], mods[["nonenb"]], 
  #                                        mods[["nonezip"]], mods[["nonezinb"]],
  #                                        mods[["Rzip"]], mods[["Rzinb"]],
  #                                        mods[["PCAzip"]], mods[["PCAzinb"]],
  #                                        mods[["Fzip"]], mods[["Fzinb"]],
  #                                        mods[["RFzip"]], mods[["RFzinb"]],
  #                                        mods[["PCAFzip"]], mods[["PCAFzinb"]],
  #                                        mods[["PCARFzip"]], mods[["PCARFzinb"]],
  #                                        criterion = "kfold"))
  # 
  # tmp.kfold$model <- rownames(tmp.kfold)
  # tmp.kfold$model <- gsub('.*\\\"([^\"]+)\\\".*', '\\1', tmp.kfold$model)
  # tmp.kfold$species <- sp
  # 
  # df.kfold <- rbind(df.kfold, tmp.kfold)

}

write.csv(df.conv, "./01_models/mod_selection/00_ModSelection_Convergence.csv", 
          row.names = F, fileEncoding = "latin1")
write.csv(df.stat, "./01_models/mod_selection/00_ModSelection_Statistics.csv", 
          row.names = F, fileEncoding = "latin1")
write.csv(df.waic, "./01_models/mod_selection/00_ModSelection_WAIC.csv", 
          row.names = F, fileEncoding = "latin1")
write.csv(df.loo, "./01_models/mod_selection/00_ModSelection_LOO.csv", 
          row.names = F, fileEncoding = "latin1")
write.csv(df.kfold, "./01_models/mod_selection/00_ModSelection_kfold16.csv", 
          row.names = F, fileEncoding = "latin1")

# plan(multisession)
# 
# for (sp in spec.list){
#   
#   mod.list <- list.files(path = "./01_models", pattern = sp)
#   mods <- list()
#   
#   for (i in 1:length(mod.list))  {
#     m <- unlist(strsplit(mod.list[i], split = "_"))[3]
#     mods[[m]]  <- readRDS(paste0("./01_models/", mod.list[i]))
#     # mods[[m]] <- brms::add_criterion(mods[[m]], "kfold", K = 8, 
#     #                                  chains = 1, threads = threading(2), overwrite = T)
#     # saveRDS(mods[[m]], paste0("./01_models/", mod.list[i]))
#   }
#   
#   tmp.kfold <- as.data.frame(loo_compare(mods[["nonepois"]], mods[["nonenb"]], 
#               mods[["nonezip"]], mods[["nonezinb"]],
#               mods[["Rzip"]], mods[["Rzinb"]],
#               mods[["PCAzip"]], mods[["PCAzinb"]],
#               mods[["Fzip"]], mods[["Fzinb"]],
#               mods[["RFzip"]], mods[["RFzinb"]],
#               mods[["PCAFzip"]], mods[["PCAFzinb"]],
#               mods[["PCARFzip"]], mods[["PCARFzinb"]],
#               criterion = "kfold"))
#   
#   tmp.kfold$model <- rownames(tmp.kfold)
#   tmp.kfold$model <- gsub('.*\\\"([^\"]+)\\\".*', '\\1', tmp.kfold$model)
#   tmp.kfold$species <- sp
#   
#   df.kfold <- rbind(df.kfold, tmp.kfold)
#   
# }
#   
# plan(sequential)

## try: instead of sd: median of absolute deviation (MAD):
  # abs. difference of value to mean, median of these differences

#### 7) posterior predictions ####
##################################-

## get posterior predictions (trends, indices, longterm trends)
newdat.total <- NULL
lt.total     <- NULL

for (sp in spec.list) {

  ## check final model
  if (nrow(df.spec[df.spec$species == sp & df.spec$modFINAL == "yes",]) > 1) stop(paste0("There is more than one final model for '", sp, "'."))
  if (nrow(df.spec[df.spec$species == sp & df.spec$modFINAL == "yes",]) < 1) stop(paste0("There is no final model for '", sp, "'."))

  mod.bgr <- df.spec$modR[df.spec$species == sp   & df.spec$modFINAL == "yes"]
  mod.fam <- df.spec$modFAM[df.spec$species == sp & df.spec$modFINAL == "yes"]
  zi.type <- df.spec$modZI[df.spec$species == sp  & df.spec$modFINAL == "yes"] 

  ## get final model and raw data
  mod <- readRDS(paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", mod.bgr, "_", year1, "-", yearx, ".RDS"))
  S.dat <- mod$data
  S.dat$year <- S.dat$year.s + year1

  ## add region (for plotting colors per region)
  if (mod.bgr %in% c("atl", "kon")) S.dat$region <- mod.bgr

  ## create newdat for model predictions
  if (mod.bgr == "both") {
    newdat <- expand.grid(year.s = seq(year1-year1, yearx-year1, length = 4*(yearx-year1)+1),
                          nat.region = levels(as.factor(S.dat$nat.region)),
                          #region = levels(as.factor(S.dat$region)),
                          OE_25  = "none",
                          off = 1, # abundance per 1 km²
                          PC1r   = mean(S.dat$PC1r),  
                          PC2r   = mean(S.dat$PC2r),  
                          PC3r   = mean(S.dat$PC3r))
  }

  if (mod.bgr %in% c("atl", "kon")) {
    newdat <- expand.grid(year.s = seq(year1-year1, yearx-year1, length = 4*(yearx-year1)+1),
                          nat.region = levels(as.factor(S.dat$nat.region)),
                          OE_25  = "none",
                          off = 1, # abundance per 1 km²
                          PC1r   = mean(S.dat$PC1r),
                          PC2r   = mean(S.dat$PC2r),
                          PC3r   = mean(S.dat$PC3r))
  }


  newdat$year <- newdat$year.s + year1

  ## create newdat2 file for slope significance
  epsilon <- 0.0001
  newdat2  <- newdat  %>% mutate(year.s = year.s + epsilon)

  ## extract predictions and newdat in lists, add overall predictions weighted by proportion of atl and kon in NRW
  pred.list   <- NULL
  pred2.list  <- NULL
  newdat.list <- NULL

  # if (mod.bgr %in% c("both", "atl")) {
  #   pred.list[["atl"]]       <- as.data.frame(t(fitted(mod, newdata = newdat,  summary = FALSE, re_formula = NA))[newdat$region == "atl",])
  #   pred2.list[["atl"]]      <- as.data.frame(t(fitted(mod, newdata = newdat2,  summary = FALSE, re_formula = NA))[newdat$region == "atl",])
  #   newdat.list[["atl"]]     <- newdat[newdat$region == "atl",]
  # }
  # 
  # if (mod.bgr %in% c("both", "kon")) {
  #   pred.list[["kon"]]       <- as.data.frame(t(fitted(mod, newdata = newdat,  summary = FALSE, re_formula = NA))[newdat$region == "kon",])
  #   pred2.list[["kon"]]      <- as.data.frame(t(fitted(mod, newdata = newdat2,  summary = FALSE, re_formula = NA))[newdat$region == "kon",])
  #   newdat.list[["kon"]]     <- newdat[newdat$region == "kon",]
  # }

  for(nr in unique(newdat$nat.region)) {
    pred.list[[nr]]       <- as.data.frame(t(fitted(mod, newdata = newdat,  summary = FALSE, re_formula = NA))[newdat$nat.region == nr,])
    pred2.list[[nr]]      <- as.data.frame(t(fitted(mod, newdata = newdat2,  summary = FALSE, re_formula = NA))[newdat$nat.region == nr,])
    newdat.list[[nr]]     <- newdat[newdat$nat.region == nr,]    
    
  }
  
  if (mod.bgr == "both") {

    pred.list[["atl"]]   <- pred.list[["A"]]*0.117 + 
      pred.list[["B"]]*0.209 + 
      pred.list[["KB"]]*0.051 + 
      pred.list[["KM"]]*0.151 + 
      #pred.list[["SB"]]*0.0001 + 
      pred.list[["ST"]]*0.351 + 
      pred.list[["metro"]]*0.121 
    
    pred2.list[["atl"]]   <- pred2.list[["A"]]*0.117 + 
      pred2.list[["B"]]*0.209 + 
      pred2.list[["KB"]]*0.051 + 
      pred2.list[["KM"]]*0.151 + 
      #pred2.list[["SB"]]*0.0001 + 
      pred2.list[["ST"]]*0.351 + 
      pred2.list[["metro"]]*0.121 
    
    pred.list[["kon"]]   <- pred.list[["A"]]*0.021 + 
      pred.list[["B"]]*0.092 + 
      pred.list[["KB"]]*0.201 + 
      pred.list[["KM"]]*0.004 + 
      pred.list[["SB"]]*0.647 + 
      pred.list[["ST"]]*0.012 + 
      pred.list[["metro"]]*0.023 
    
    pred2.list[["kon"]]   <- pred2.list[["A"]]*0.021 + 
      pred2.list[["B"]]*0.092 + 
      pred2.list[["KB"]]*0.201 + 
      pred2.list[["KM"]]*0.004 + 
      pred2.list[["SB"]]*0.647 + 
      pred2.list[["ST"]]*0.012 + 
      pred2.list[["metro"]]*0.023 
      
    pred.list[["overall"]]  <- pred.list[["atl"]]*0.555 + pred.list[["kon"]]*0.445
    pred2.list[["overall"]]  <- pred2.list[["atl"]]*0.555 + pred2.list[["kon"]]*0.445
    newdat.list[["overall"]] <- newdat.list[["kon"]] <- newdat.list[["atl"]] <- newdat.list[["A"]]
    newdat.list[["overall"]]$nat.region <- "overall"
    newdat.list[["atl"]]$nat.region <- "atl"
    newdat.list[["kon"]]$nat.region <- "kon"
  }

  ## loop through regions
  for (r in names(pred.list)) {
    pred.r   <- pred.list[[r]]
    pred2.r  <- pred2.list[[r]]
    newdat.r <- newdat.list[[r]]

    ## trend & 95% CrI
    ##################-
    newdat.r$lwr <- apply(X = pred.r, MARGIN = 1, FUN = quantile, prob = 0.025)
    newdat.r$fit <- apply(X = pred.r, MARGIN = 1, FUN = quantile, prob = 0.5)
    newdat.r$upr <- apply(X = pred.r, MARGIN = 1, FUN = quantile, prob = 0.975)

    ## Calculate the difference and divide by epsilon - this is analogous to dx/dt
    diff <- (pred2.r - pred.r) /epsilon

    ## using the differences, we calculate the mean and lower and upper quantiles
    newdat.r$lwr.diff  <- apply(diff, MARGIN = 1, FUN = quantile, prob = 0.025)
    newdat.r$fit.diff  <- apply(diff, MARGIN = 1, FUN = quantile, prob = 0.5)
    newdat.r$upr.diff  <- apply(diff, MARGIN = 1, FUN = quantile, prob = 0.975)

    newdat.r$incr     <- NA; newdat.r$decr     <- NA
    newdat.r$fit.incr <- NA; newdat.r$fit.decr <- NA

    ## select slope parts with sign. increase/decrease (95% CrI of slope does not include 0)
    newdat.r$incr[newdat.r$fit.diff > 0 & newdat.r$lwr.diff > 0] <- newdat.r$fit.diff[newdat.r$fit.diff > 0 & newdat.r$lwr.diff > 0]
    newdat.r$decr[newdat.r$fit.diff < 0 & newdat.r$upr.diff < 0] <- newdat.r$fit.diff[newdat.r$fit.diff < 0 & newdat.r$upr.diff < 0]
    newdat.r$fit.incr[!is.na(newdat.r$incr)] <- newdat.r$fit[!is.na(newdat.r$incr)] # this is needed for plotting only
    newdat.r$fit.decr[!is.na(newdat.r$decr)] <- newdat.r$fit[!is.na(newdat.r$decr)] # this is needed for plotting only

    ## index & 95% CrI
    ##################-
    pred.i <- matrix(nrow = nrow(pred.r), ncol = ncol(pred.r))

    for (i in 1:ncol(pred.r)) {pred.i[,i] <- pred.r[,i]/pred.r[newdat.r$year == yearI,i]}

    pred.i <- pred.r/apply(X = pred.r[newdat.r$year == yearI,], MARGIN = 1, FUN = quantile, probs = 0.5)

    newdat.r$lwr.i <- apply(X = pred.i, MARGIN = 1, FUN = quantile, probs = 0.025)
    newdat.r$fit.i <- apply(X = pred.i, MARGIN = 1, FUN = quantile, probs = 0.5)
    newdat.r$upr.i <- apply(X = pred.i, MARGIN = 1, FUN = quantile, probs = 0.975)

    ## longterm change & 95% CrI
    ############################-
    diff.lt <- pred.r[newdat.r$year == yearP,] - pred.r[newdat.r$year == yearP-period,]

    lt.tbl <- data.frame(nat.region = r)
    lt.tbl$lwr    <- apply(X = diff.lt, MARGIN = 1, FUN = quantile, prob = 0.025)
    lt.tbl$lwr25  <- apply(X = diff.lt, MARGIN = 1, FUN = quantile, prob = 0.25)
    lt.tbl$fit    <- apply(X = diff.lt, MARGIN = 1, FUN = quantile, prob = 0.5)
    lt.tbl$upr75  <- apply(X = diff.lt, MARGIN = 1, FUN = quantile, prob = 0.75)
    lt.tbl$upr    <- apply(X = diff.lt, MARGIN = 1, FUN = quantile, prob = 0.975)
    lt.tbl$BayesP <- length(diff.lt[diff.lt > 0]) / length(diff.lt)


    ## merge data
    newdat.r$species <- sp
    lt.tbl$species   <- sp
    lt.tbl$longterm  <- paste0(yearP - period, "-", yearP)
    newdat.total     <- rbind(newdat.total, newdat.r)
    lt.total         <- rbind(lt.total, lt.tbl)

  } ## end of region loop
  
  write.csv(newdat.total, paste0("./02_output/PosteriorPreditions_", year1, "-", yearx ,".csv"),
            row.names = F, fileEncoding = "latin1")
  write.csv(lt.total,     paste0("./02_output/LongtermTrend_", year1, "-", yearx , "_", yearP - period, "-", yearP, ".csv"),
            row.names = F, fileEncoding = "latin1")
  
  ## plot posterior predictions
  #############################-

  ## plot abundance trend
  g <- ggplot() +

    # add horizontal line at y = 0
    geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.5) +

    # add raw data (colored dots)
    geom_beeswarm(data = S.dat, aes(x = year, y = ab2/2, color = nat.region),
                  size = 0.5, cex = 0.3, alpha = 0.5) +

    # add trend
    geom_ribbon(data = newdat.total[newdat.total$species == sp,],
                aes(x = year, ymin = lwr, ymax = upr, color = nat.region, fill = nat.region),
                alpha = 0.3) +
    geom_line(data = newdat.total[newdat.total$species == sp,],
              aes(x = year, y = fit, color = nat.region),
              lwd = 1) +

    # add sign. slope
    geom_line(aes(year, fit.incr, group = nat.region),
              data = newdat.total[newdat.total$species == sp,],
              col = "black", lty = 1, linewidth = 1, na.rm = TRUE) +
    geom_line(aes(year, fit.decr, group = nat.region), 
              data = newdat.total[newdat.total$species == sp,],
              col = "black", lty = 1, linewidth = 1, na.rm = TRUE) +


    # add labs
    ylim(0, max(c(S.dat$ab2/2, newdat.total$upr[newdat.total$species == sp]))) +
    labs(x = "year", y = "abundance per km²",
         title = paste0(sp, " (", mod.fam, ", ", zi.type, ")"),
         subtitle = "black = period of robust increase/decrease") +
    facet_wrap(~nat.region) +
    theme_classic()

  # add region specific colors
  if (mod.bgr == "both") {
    g <- g +
      scale_color_manual(values = c("atl" = "blue", "kon" = "orange", "overall" = "grey40")) +
      scale_fill_manual( values = c("atl" = "blue", "kon" = "orange", "overall" = "grey40"))
  }

  if (mod.bgr == "atl") {
    g <- g +
      scale_color_manual(values = c("atl" = "blue")) +
      scale_fill_manual( values = c("atl" = "blue"))
  }

  if (mod.bgr == "kon") {
    g <- g +
      scale_color_manual(values = c("kon" = "orange")) +
      scale_fill_manual( values = c("kon" = "orange"))
  }

  print(g)


  ## plot index
  g <- ggplot(newdat.total[newdat.total$species == sp,],
              aes(x = year, y = fit.i, ymin = lwr.i, ymax = upr.i)) +

    # add horizontal lines at y = 0 and y = 1
    geom_hline(yintercept = 0, show.legend = FALSE, linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = 1, show.legend = FALSE, linetype = "dashed", linewidth = 0.5) +

    # add trend
    geom_ribbon(aes(color = region, fill = region), alpha = 0.2) +
    geom_line(aes(color = region), lwd = 1) +

    # add labs
    ylim(0, max(newdat.total$upr.i[newdat.total$species == sp])) +
    labs(y = "Index", x = "survey year",
         title = paste0(sp, " (", mod.fam, ", ", zi.type, ")")) + 
    theme_classic()

  # add region specific colors
  if (mod.bgr == "both") {
    g <- g +
      scale_color_manual(values = c("atl" = "blue", "kon" = "orange", "overall" = "grey40")) +
      scale_fill_manual( values = c("atl" = "blue", "kon" = "orange", "overall" = "grey40"))
  }

  if (mod.bgr == "atl") {
    g <- g +
      scale_color_manual(values = c("atl" = "blue")) +
      scale_fill_manual( values = c("atl" = "blue"))
  }

  if (mod.bgr == "kon") {
    g <- g +
      scale_color_manual(values = c("kon" = "orange")) +
      scale_fill_manual( values = c("kon" = "orange"))
  }

  print(g)

  ## plot longterm trend
  g <- ggplot() +

    # add horizontal line at 0
    geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.5) +

    # add longterm trend
    geom_pointrange(aes(x = region, y = fit, ymin = lwr, ymax = upr, color = region),
                    data = lt.total[lt.total$species == sp,],
                    linewidth = 1, fatten = 8) +
    geom_pointrange(aes(x = region, y = fit, ymin = lwr25, ymax = upr75, color = region),
                    data = lt.total[lt.total$species == sp,],
                    linewidth = 2, fatten = 8) +

    # add labs
    labs(x = "biogeographical region",
         y = paste0(period, "-year abundance difference \n (", yearP - period, "-", yearP, ")"),
         title = paste0(sp, " (", mod.fam, ", ", zi.type, ")"),
         subtitle = "Longterm trend with 50% (thick) and 95% (thin) CrI") +
      theme_classic() +
      theme(legend.position = "none")

  # add region specific colors
  if(mod.bgr == "both") {
    g <- g +
      scale_color_manual(values = c("atl" = "blue", "kon" = "orange", "overall" = "grey40"))
  }

  if (mod.bgr == "atl") {
    g <- g +
      scale_color_manual(values = c("atl" = "blue"))
  }

  if (mod.bgr == "kon") {
    g <- g +
      scale_color_manual(values = c("kon" = "orange"))
  }

  print(g)


} ## end of species loop
