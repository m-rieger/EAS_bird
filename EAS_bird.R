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
    "svDialogs" # for pop-up window
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
source("./EAS_bird_plots.R")

## define first (year1) and last year (yearx)
year1 <- 2002
yearx <- 2020

## define baseline year (for index)
yearI <- 2006

## define time range for longterm trend and last year or the range
period <- 12
yearP  <- yearx-1 # here: 2019

## do you want to save plots per species as pdf (TRUE) or not (FALSE)?
PDF <- TRUE

## define some model parameters
## --> do you want to do a model selection (TRUE) or only model the 
##     final model (FALSE) 
mod.sel <- FALSE

## --> do you want to do a (faster) test run, e.g. with fewer iterations (TRUE) 
##     or the original run with (FALSE)?
test <- TRUE

if (test == TRUE) {
  nc <- 4           # no. of chains
  ni <- 2000        # no. of iterations (incl. warm-up)
  nw <- 1000        # no. of iterations (warm-up only)
  a.delta <- 0.95   # adapt_delta value
  max.td  <- 10     # max_treedepth value
}

if (test == FALSE) {
  nc <- 4           # no. of chains
  ni <- 3000        # no. of iterations (incl. warm-up)
  nw <- 1000        # no. of iterations (warm-up only)
  a.delta <- 0.99 # adapt_delta value
  max.td  <- 10     # max_treedepth value
}


# define backend and cores based on whether cmdstanr is installed or not
cores <- detectCores() # get the number of cores for parallelization

if (csr == TRUE) {
  backend <- "cmdstanr"
  ncores  <- 4
  thread  <- cores/ncores
}

if (csr == FALSE) {
  backend <- "rstan"
  ncores  <- min(cores-1, nc)
  thread  <- NULL
}

## define species subset
# in case you want to model only some species, e.g. c("Common Blackbird", "White Wagtail")
# if NULL, all available species will be used
Spec.sub <- NULL

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

## based on the default deviation of 20% (for comparison)
dat <- obs.eff(data = dat,
               ID = "ID",               # default column name
               abundance = "abundance", # default column name
               year = "year",           # default column name
               OE = 0.20,               # default threshold of 25% to delineate observer effect
               add.df = TRUE)           # default: returns an additional dataframe with observer effects per site and year

## based on the default deviation of 30% (for comparison)
dat <- obs.eff(data = dat,
               ID = "ID",               # default column name
               abundance = "abundance", # default column name
               year = "year",           # default column name
               OE = 0.30,               # default threshold of 25% to delineate observer effect
               add.df = TRUE)           # default: returns an additional dataframe with observer effects per site and year

#### 2.2) further adaptions ####
###############################-

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

## save rotated loadings
write.csv(df.PCA, paste0("./02_output/PCA_scores_", year1, "-", yearx, ".csv"), 
          row.names = F, fileEncoding = "latin1")
write.csv(rot.df, paste0("./02_output/PCA_rotated_loadings_", year1, "-", yearx, ".csv"), 
          row.names = F, fileEncoding = "latin1")

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
# NOTE: weights were removed from the analysis but they might be handy if you 
# want to check whether your habitat/landscape/natural regions are balanced 
# (weight ~ 1) or not (overrepresented: weight < 1, underrepresented: weight > 1)

#### 4.1) weights for both biogeographical regions combined ####
################################################################-
dat <- weight(df.habitat = df.land,   # landscape data
               habitat = "landscape", # column name for habitat
               area = "area_NRW",     # column name for area
               df.data = dat,         # raw data
               ID = "ID",             # default column name
               year = "year",         # default column name
               by_spec = TRUE,        # calculate weight by species (e.g., if number of sites differs between species due to exclusion)
               species = "species",   # default column name
               by_reg = FALSE,        # default
               add.df = TRUE)         # default: returns an additional dataframe with weights per habitat and year (and species)
# remark: here, KB.yes and SB.yes are not present in the dataset

# rename "weight" to prevent overwriting in the next step
colnames(dat)[colnames(dat) == "weight"] <- "weight.l"

#### 4.2) weights for each biogeographical region separately ####
#################################################################-
## this is needed for species which are only modelled for one biogeographical region to guarantee a
## representative weighting of landscapes per region

dat <- weight(df.habitat = df.land,  # landscape data
               habitat = "landscape", # column name for habitat
               area = "area_NRW",     # column name for area
               df.data = dat,        # raw data
               ID = "ID",             # default column name
               year = "year",         # default column name
               by_spec = TRUE,        # calculate weight by species (e.g., if number of sites differs between species due to exclusion)
               species = "species",   # default column name
               by_reg = TRUE,         # calculate weights separately per region
               region = "region",     # default column name
               add.df = TRUE)         # default: returns an additional dataframe with weights per habitat and year (and species)
# remark: here, several landscape combinations are not present in 'atl' and 'kon'

# rename "weight"
colnames(dat)[colnames(dat) == "weight"] <- "weight.l.r"

#### 5) Model ####
##################-

## create list of species which you want to model, either

# ... automatically based on whether it can be modelled or not
if (is.null(Spec.sub)) {
  spec.list <- sort(as.character(unique(df.spec$species[!is.na(df.spec$modFAM)])))
}

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

## save raw data
write.csv(dat[dat$species %in% spec.list,], "./02_output/Data_modelling.csv",
          fileEncoding = "latin1")

## for-loop for all species in spec.list

for (sp in spec.list) {

  ## get information from df.spec
  mod.bgr   <- unique(df.spec$modR[df.spec$species == sp])
  
  ## list of model family and zi-type for model selection
  if(mod.sel == TRUE) {
    mod.fam.l <- df.spec$modFAM[df.spec$species == sp]
    zi.type.l <- df.spec$modZI[df.spec$species == sp]    
  }
  
  ## only one model family and zi-type for final models
  if(mod.sel == FALSE) {
    mod.fam.l <- df.spec$modFAM[df.spec$species == sp & df.spec$modFINAL == TRUE]
    zi.type.l <- df.spec$modZI[df.spec$species == sp & df.spec$modFINAL == TRUE]    
  }

  ## subset dataframe 
  S.dat <- droplevels(dat[dat$species == sp,])

  ## for-loop for different model per species (for comparison)
  for (m in 1:length(mod.fam.l)) {

    mod.fam <- mod.fam.l[m]
    zi.type <- zi.type.l[m]

    #### 5.1) define model formula
    ##############################-

    ## 5.1a) for poisson and negative binomial (without zero inflation)
    if (mod.fam %in% c("pois", "nb")) {

      if (zi.type == "none") {
        mod.form <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                         OE_25 +
                         offset(log(off)) +
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                         (1|ID))
      }
    }

    ## 5.1b) for zero-inflated poisson and zero-inflated negative binomial
    if (mod.fam %in% c("zip", "zinb")) {

      if (zi.type == "none") {
        mod.form <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                         OE_25 +
                         offset(log(off)) +
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                         (1|ID),
                       zi ~ 1)
      }

      if (zi.type == "R") {
        mod.form <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                         OE_25 +
                         offset(log(off)) +
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                         (1|ID),
                       zi ~ nat.region)
      }

      if (zi.type == "PCA") {
        mod.form <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                         OE_25 +
                         offset(log(off)) +
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                         (1|ID),
                       zi ~ PC1r + PC2r + PC3r)
      }

      if (zi.type == "F") {
        mod.form <- bf(ab2 ~ s(year.s, by = nat.region) + nat.region +
                         OE_25 +
                         offset(log(off)) +
                         poly(PC1r, 2) + poly(PC2r, 2) + poly(PC3r, 2) +
                         (1|ID),
                       zi ~ (1|ID))
      }
    }

    #### 5.2) define model family
    #############################-

    if (mod.fam == "pois") fam <- poisson(link = "log")
    if (mod.fam == "nb")   fam <- negbinomial(link = "log")
    if (mod.fam == "zip")  fam <- zero_inflated_poisson(link = "log")
    if (mod.fam == "zinb") fam <- zero_inflated_negbinomial(link = "log")


    #### 5.3) set priors
    ####################-
    ## set more specific priors (less narrow for PCs)

    priors <- get_prior(mod.form,
                        data = S.dat, family = fam)

    priors$prior[priors$dpar == "" & priors$class == "b" & 
                   priors$coef != ""] <- "normal(0, 2.5)"
    priors$prior[priors$coef %in% c("polyPC1r21", "polyPC1r22", 
                                    "polyPC2r21", "polyPC2r22", 
                                    "polyPC3r21", "polyPC3r22")] <- "normal(0, 10)"
    priors$prior[priors$prior == "" & priors$dpar == "zi" & 
                   priors$class == "b" & priors$coef != ""] <- "normal(0, 2.5)"
    priors$prior[priors$coef %in% c("PC1r", "PC2r", "PC3r")] <- "normal(0, 10)"

    #### 5.4) run model
    ###################-

    mod <- brm(mod.form,
               data = S.dat, family = fam,
               cores = ncores, iter = ni, warmup = nw, chains = nc,
               backend = backend, threads = threading(thread),
               thin = 1,
               control = list(adapt_delta = a.delta, max_treedepth = max.td), 
               prior = priors) 
    
    ## add kfold in case of model selection
    if(mod.sel == TRUE) {
    plan(multisession); mc.cores = parallel::detectCores()
    mod <- add_criterion(mod, "kfold", K = 16, chains = 4, threads = threading(1))
    plan(sequential)      
    }


    ## save model
    saveRDS(mod, paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", 
                        mod.bgr, "_", year1, "-", yearx, ".RDS"))

  } # end of for-loop (model)
} # end of for-loop (species)


#### 6) posterior predictive checks and model convergence ####
##############################################################-

## check model convergence and posterior predictions of response variable
## to chose the best fitting model structure
df.conv <- data.frame()
df.stat <- data.frame()
df.kfold <- data.frame()

## loop through species
for (sp in spec.list){

  ## load models
  mod.list <- list.files(path = "./01_models", pattern = sp)
  mods <- list()
  for (i in 1:length(mod.list))  {
    m <- unlist(strsplit(mod.list[i], split = "_"))[3]
    mods[[m]]  <- readRDS(paste0("./01_models/", mod.list[i]))
    print(mod.list[i])
    print(mods[[m]]$formula)
    print(mods[[m]]$family)
  }

  ## posterior predictions of response variable
  tmp.stat <- print(mod.stat(model.list = mods, model.name = names(mods),
                 response = "ab2",
                 plot.stats = TRUE,
                 spec = sp))

  ## model convergence
  tmp.conv <- print(mod.conv(model.list = mods, model.name = names(mods),
                 td = max.td,
                 plot.conv = TRUE,
                 spec = sp))
  
  ## kfold selection criterion
  if(mod.sel == TRUE) {
    ## get kfolds
    kfoldL <- list()
    for(m in names(mods)) kfoldL[[m]] <- mods[[m]]$criteria$kfold
    
    tmp.kfold <- as.data.frame(loo_compare(kfoldL))
    tmp.kfold$model <- row.names(tmp.kfold)
    
    ## add upr and lwr (1.96*SE, 1*SE)
    tmp.kfold <- tmp.kfold %>% rowwise() %>%
      mutate(lwr = elpd_diff - 1.96*se_diff,
             upr = elpd_diff + 1.96*se_diff,
             lwr5 = elpd_diff - se_diff,
             upr5 = elpd_diff + se_diff) %>% ungroup()
    
    ## add npar based on model matrix (to get most parsimonious model)
    for(m in names(mods)) {
      tmp.kfold$npar[tmp.kfold$model == m] <- length(as_draws_list(mods[[m]])[[1]])
    }
    
    tmp.kfold$species <- sp
    
    ## define best and most parsimonious model
    tmp.kfold$top <- ifelse(tmp.kfold$upr5 >= 0, "yes SE", 
                     ifelse(tmp.kfold$upr >= 0, "yes 1.96*SE", "no"))
    
    ## IMPORTANT: check issues with model convergence of final models and if any, 
    # select the next best and most parsimonious model

    tmp.kfold <- tmp.kfold %>% 
      mutate("final SE" = ifelse(npar == min(npar[top == "yes SE"]), TRUE, FALSE),
             "final 1.96*SE" = ifelse(npar == min(npar[top %in% c("yes 1.96*SE", "yes SE")]), TRUE, FALSE)) %>%
      ungroup()
    
    ## plot model selection
    print(ggplot(tmp.kfold) +
      geom_hline(aes(yintercept = 0), lty = "dashed") +
      geom_pointrange(aes(x = model, y = elpd_diff, 
                          ymin = elpd_diff -1.96*se_diff, ymax = elpd_diff +1.96*se_diff, 
                          color = `final SE`, pch = top)) +
      geom_pointrange(aes(x = model, y = elpd_diff, 
                          ymin = elpd_diff -se_diff, ymax = elpd_diff +se_diff, 
                          color = `final SE`, pch = top), lwd = 1) +
      scale_color_viridis_d("final mod", end = 0.8, option = "inferno") +
      scale_shape_manual("within 0 diff.", values = c("yes SE" = 19, "yes 1.96*SE" = 13, "no" = 1)) +
      ggtitle(sp) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)))
    
    df.kfold <- rbind(df.kfold, tmp.kfold)
  }
  df.stat <- rbind(df.stat, tmp.stat)
  df.conv <- rbind(df.conv, tmp.conv)
  
}

write.csv(df.stat, paste0("./01_models/Statistics_", year1, "-", yearx ,".csv"), row.names = F)
write.csv(df.conv, paste0("./01_models/Convergence_", year1, "-", yearx ,".csv"), row.names = F)
if(mod.sel) write.csv(df.kfold, paste0("./01_models/Kfold_", year1, "-", yearx ,".csv"), row.names = F)

#### 7) posterior predictions ####
##################################-

#### 7.1) get posterior predictions (trends, indices, longterm trends) ####
###########################################################################-
Coef.total <- NULL
newdat.total <- NULL
lt.total     <- NULL
pw.total     <- NULL

## original names of coefficients
cL <- c("b_Intercept", 
        "b_nat.regionB", "b_nat.regionKB", "b_nat.regionKM", 
        "b_nat.regionmetro", "b_nat.regionSB", "b_nat.regionST",
        "b_OE_25none", "b_OE_25positive", 
        "b_polyPC1r21", "b_polyPC2r21", "b_polyPC3r21", 
        "b_polyPC1r22", "b_polyPC2r22", "b_polyPC3r22", 
        "bs_syear.s:nat.regionA_1", "bs_syear.s:nat.regionB_1", 
        "bs_syear.s:nat.regionKB_1", "bs_syear.s:nat.regionKM_1", 
        "bs_syear.s:nat.regionmetro_1", "bs_syear.s:nat.regionSB_1", 
        "bs_syear.s:nat.regionST_1",
        "sds_syear.snat.regionA_1", "sds_syear.snat.regionB_1", 
        "sds_syear.snat.regionKB_1", "sds_syear.snat.regionKM_1", 
        "sds_syear.snat.regionmetro_1", "sds_syear.snat.regionSB_1", 
        "sds_syear.snat.regionST_1",
        "sd_ID__Intercept", 
        "b_zi_Intercept", 
        "b_zi_nat.regionB", "b_zi_nat.regionKB", "b_zi_nat.regionKM", 
        "b_zi_nat.regionmetro", "b_zi_nat.regionSB", "b_zi_nat.regionST",
        "b_zi_PC1r", "b_zi_PC2r", "b_zi_PC3r",
        "sd_ID__zi_Intercept")

## understandable names of coefficients (with region in zi-model)
cRL   <- c("intercept [negative, A]", 
           "region [B]", "region [KB]", "region [KM]", 
           "region [metro]", "region [SB]", "region [ST]", 
           "observer effect [none]", "observer effect [positive]", 
           "PC1 - linear", "PC2 - linear", "PC3 - linear", 
           "PC1 - quadratic", "PC2 - quadratic", "PC3 - quadratic", 
           "survey year:region [A]", "survey year:region [B]", 
           "survey year:region [KB]", "survey year:region [KM]", 
           "survey year:region [metro]", "survey year:region [SB]", 
           "survey year:region [ST]", 
           "survey year:region [A] - smoother", "survey year:region [B] - smoother",                                     
           "survey year:region [KB] - smoother", "survey year:region [KM] - smoother",                                     
           "survey year:region [metro] - smoother", "survey year:region [SB] - smoother",                                     
           "survey year:region [ST] - smoother",                                   
           "ID intercept", 
           "zi - intercept [A]", 
           "zi - region [B]", "zi - region [KB]", "zi - region [KM]", 
           "zi - region [metro]", "zi - region [SB]", "zi - region [ST]", 
           "zi - PC1 - linear", "zi - PC2 - linear", "zi - PC3 - linear", 
           "zi - ID intercept")

## understandable names of coefficients (without region in zi-model)
cPL  <- cRL
cPL[cPL == "zi - intercept [A]"] <- "zi - intercept"

for (sp in spec.list) {

  ## check final model
  if (nrow(df.spec[df.spec$species == sp & df.spec$modFINAL == "yes",]) > 1) stop(paste0("There is more than one final model for '", sp, "'."))
  if (nrow(df.spec[df.spec$species == sp & df.spec$modFINAL == "yes",]) < 1) stop(paste0("There is no final model for '", sp, "'."))

  mod.bgr <- df.spec$modR[df.spec$species == sp   & df.spec$modFINAL == "yes"]
  mod.fam <- df.spec$modFAM[df.spec$species == sp & df.spec$modFINAL == "yes"]
  zi.type <- df.spec$modZI[df.spec$species == sp  & df.spec$modFINAL == "yes"] 

  ## get final model and raw data
  mod <- readRDS(paste0("./01_models/mod_", sp, "_", zi.type,  mod.fam, "_R", 
                        mod.bgr, "_", year1, "-", yearx, ".RDS"))
  S.dat <- mod$data
  S.dat$year <- S.dat$year.s + year1
  
  ## coefficient estimates
  mod.coef <- cL[which(c(cL %in% colnames(as.data.frame(mod))))]

  Coef.Est <- as.data.frame(mod)[,mod.coef] # extract from here
  
  Coef.df <- data.frame(Coef = mod.coef)
  
  
  for (i in mod.coef) {
    
    Coef.df$lwr[Coef.df$Coef == i]    <- quantile(Coef.Est[, c(i)], probs = 0.025)
    Coef.df$fit[Coef.df$Coef == i]    <- quantile(Coef.Est[, c(i)], probs = 0.5)
    Coef.df$upr[Coef.df$Coef == i]    <- quantile(Coef.Est[, c(i)], probs = 0.975)
    Coef.df$BayesP[Coef.df$Coef == i] <- mean(Coef.Est[, c(i)] > 0)
    
  }
  
  Coef.df$species    <- sp
  Coef.df$modR    <- mod.bgr
  Coef.df$modZI   <- zi.type
  Coef.df$modFAM  <- mod.fam
  Coef.df$period  <- as.character(paste0(year1, " bis ", yearx))
  Coef.df$coefficient <- NA
  
  for (c in Coef.df$Coef) {
    
    ifelse(zi.type %in% c("R"), 
           Coef.df$coefficient[Coef.df$Coef == c] <- cRL[cL == c],
           Coef.df$coefficient[Coef.df$Coef == c] <- cPL[cL == c])
    
  }

  Coef.total <- rbind(Coef.total, Coef.df)

  ## create newdat for model predictions
  newdat <- expand.grid(year.s = seq(year1-year1, yearx-year1, length = 4*(yearx-year1)+1),
                        nat.region = levels(as.factor(S.dat$nat.region)),
                        #region = levels(as.factor(S.dat$region)),
                        OE_25  = "none",
                        off = 1) # abundance per 1 km²

  ## extract mean PC values (per nat. region)
  newdat <- newdat %>% group_by(nat.region) %>%
    mutate(PC1r = mean(S.dat$PC1r[S.dat$nat.region == first(nat.region)]),
           PC2r = mean(S.dat$PC2r[S.dat$nat.region == first(nat.region)]),
           PC3r = mean(S.dat$PC3r[S.dat$nat.region == first(nat.region)]))

  newdat$year <- newdat$year.s + year1

  ## extract predictions and newdat per nat. region in lists, add overall
  # predictions and per biogeographical region weighted by proportion of
  # nat. regions in NRW (or biog. regions)
  predL   <- list()
  newdatL <- list()

  for(nr in unique(newdat$nat.region)) {
    predL[[nr]]  <- as.data.frame(t(fitted(mod, newdata = newdat, summary = FALSE,
                                           re_formula = NA))[newdat$nat.region == nr,])
    newdatL[[nr]] <- newdat[newdat$nat.region == nr,]

  }

  predL[["atl"]]   <- predL[["A"]]*0.117 +
    predL[["B"]]*0.209 +
    predL[["KB"]]*0.051 +
    predL[["KM"]]*0.151 +
    #predL[["SB"]]*0.0001 +
    predL[["ST"]]*0.351 +
    predL[["metro"]]*0.121


  predL[["kon"]]   <- predL[["A"]]*0.021 +
    predL[["B"]]*0.092 +
    predL[["KB"]]*0.201 +
    predL[["KM"]]*0.004 +
    predL[["SB"]]*0.647 +
    predL[["ST"]]*0.012 +
    predL[["metro"]]*0.023

  predL[["overall"]]  <- predL[["atl"]]*0.555 + predL[["kon"]]*0.445
  newdatL[["overall"]] <- newdatL[["kon"]] <- newdatL[["atl"]] <- newdatL[["A"]]
  newdatL[["overall"]]$nat.region <- "overall"
  newdatL[["atl"]]$nat.region <- "atl"
  newdatL[["kon"]]$nat.region <- "kon"


  ## loop through different predictions
  for (r in names(predL)) {
    pred.r   <- predL[[r]]
    newdat.r <- newdatL[[r]]

    ## trend & 95% CrI
    ##################-
    newdat.r$lwr <- apply(X = pred.r, MARGIN = 1, FUN = quantile, prob = 0.025)
    newdat.r$fit <- apply(X = pred.r, MARGIN = 1, FUN = quantile, prob = 0.5)
    newdat.r$upr <- apply(X = pred.r, MARGIN = 1, FUN = quantile, prob = 0.975)

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

    ## pairwise comparison & 95% CrI
    ################################-
    pw.tbl <- data.frame(nat.region = r, year = year1:(yearx-1))
    newdat.r$pw.slope <- NA # needed for plotting
    newdat.r$pw.slope.i <- NA # needed for plotting

    for (y in year1:(yearx-1)) {

      diff.pw <- pred.r[newdat.r$year == y+1,] - pred.r[newdat.r$year == y,]
      diff.pw.i <- pred.i[newdat.r$year == y+1,] - pred.i[newdat.r$year == y,]

      pw.tbl$lwr[pw.tbl$year == y]    <- apply(X = diff.pw, MARGIN = 1, FUN = quantile, prob = 0.025)
      pw.tbl$fit[pw.tbl$year == y]    <- apply(X = diff.pw, MARGIN = 1, FUN = quantile, prob = 0.5)
      pw.tbl$upr[pw.tbl$year == y]    <- apply(X = diff.pw, MARGIN = 1, FUN = quantile, prob = 0.975)
      pw.tbl$BayesP[pw.tbl$year == y] <- length(diff.pw[diff.pw > 0])/length(diff.pw)
      pw.tbl$lwr.i[pw.tbl$year == y]    <- apply(X = diff.pw.i, MARGIN = 1, FUN = quantile, prob = 0.025)
      pw.tbl$fit.i[pw.tbl$year == y]    <- apply(X = diff.pw.i, MARGIN = 1, FUN = quantile, prob = 0.5)
      pw.tbl$upr.i[pw.tbl$year == y]    <- apply(X = diff.pw.i, MARGIN = 1, FUN = quantile, prob = 0.975)
      pw.tbl$BayesP.i[pw.tbl$year == y] <- length(diff.pw.i[diff.pw.i > 0])/length(diff.pw.i)

      if(pw.tbl$BayesP[pw.tbl$year == y] >= 0.975) {
        newdat.r$pw.slope[floor(newdat.r$year) == y] <- "incr"
      }

      if(pw.tbl$BayesP[pw.tbl$year == y] <= 0.025) {
        newdat.r$pw.slope[floor(newdat.r$year) == y] <- "decr"
      }

      if(pw.tbl$BayesP.i[pw.tbl$year == y] >= 0.975) {
        newdat.r$pw.slope.i[floor(newdat.r$year) == y] <- "incr"
      }

      if(pw.tbl$BayesP.i[pw.tbl$year == y] <= 0.025) {
        newdat.r$pw.slope.i[floor(newdat.r$year) == y] <- "decr"
      }

    }

    ## merge data
    newdat.r$species <- sp
    lt.tbl$species   <- sp
    pw.tbl$species   <- sp
    lt.tbl$longterm  <- paste0(yearP - period, "-", yearP)
    newdat.total     <- rbind(newdat.total, newdat.r)
    lt.total         <- rbind(lt.total, lt.tbl)
    pw.total         <- rbind(pw.total, pw.tbl)

  } ## end of region loop

  write.csv(newdat.total, paste0("./02_output/PosteriorPredictions_", year1, "-", yearx ,".csv"),
            row.names = F, fileEncoding = "latin1")
  write.csv(lt.total,     paste0("./02_output/LongtermTrend_", year1, "-", yearx , "_", yearP - period, "-", yearP, ".csv"),
            row.names = F, fileEncoding = "latin1")
  write.csv(pw.total,     paste0("./02_output/PairWise_", year1, "-", yearx , "_", yearP - period, "-", yearP, ".csv"),
            row.names = F, fileEncoding = "latin1")
  write.csv(Coef.total, paste0("./02_output/CoefficientEstimates_", year1, "-", yearx ,".csv"), 
            row.names = F, fileEncoding = "latin1")

  #### 7.2) plot posterior predictions (trends, indices, longterm trends) ####
  ############################################################################-
  
  ## S.dat for atl, kon and overall (for plotting reasons)
  S.dat2 <- left_join(S.dat, unique(df.PCA[, c("ID", "region")]), by = "ID")
  S.dat2$nat.region <- S.dat2$region
  tmp <- S.dat2 
  tmp$nat.region <- "overall"
  S.dat2 <- rbind(S.dat2, tmp); rm(tmp)
  
  ## fit.slope (for plotting reasons of incr or decr)
  newdat.total$fit.slope[!is.na(newdat.total$pw.slope)] <- newdat.total$fit[!is.na(newdat.total$pw.slope)]
  newdat.total$fit.slope.i[!is.na(newdat.total$pw.slope.i)] <- newdat.total$fit.i[!is.na(newdat.total$pw.slope.i)]
  
  if(PDF) pdf(file = paste0("./02_output/PostPred_", sp, "_",  year1, "-", yearx ,".pdf"))
  
  ## nat. regions
  ###############-
  tmp <- newdat.total[newdat.total$species == sp & 
                        !newdat.total$nat.region %in% c("atl", "kon", "overall"),]
  
  ## plot abundance trend
  plot.ab(df.raw = S.dat, newdat = tmp, group = "nat.region")
  
  ## plot index
  plot.in(newdat = tmp, group = "nat.region")
  
  ## plot longterm
  tmp <- lt.total[lt.total$species == sp & 
                    !lt.total$nat.region %in% c("atl", "kon", "overall"),]
  
  plot.lt(newdat = tmp, group = "nat.region")
  
  ## atl, kon, overall
  ####################-
  tmp <- newdat.total[newdat.total$species == sp & 
                        newdat.total$nat.region %in% c("atl", "kon", "overall"),]
   
  ## plot abundance trend
  plot.ab(df.raw = S.dat2, newdat = tmp, group = "region")

  ## plot index
  plot.in(newdat = tmp, group = "region")
  
  ## plot longterm
  tmp <- lt.total[lt.total$species == sp & 
                    lt.total$nat.region %in% c("atl", "kon", "overall"),]
  
  plot.lt(newdat = tmp, group = "region")

  if(PDF) dev.off()

} ## end of species loop
