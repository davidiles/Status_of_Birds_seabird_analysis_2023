# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------
my_packs = c('tidyverse','readxl','RColorBrewer','viridis','jagsUI','mgcv','ggrepel','scales')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())
theme_set(theme_bw())
setwd("D:/iles_ECCC/Seabirds/Status_of_Birds_seabird_analysis_2022/")
`%!in%` <- Negate(`%in%`)

# ************************************************
# Read in data
# ************************************************

# --------------------
# Complete counts
# --------------------

ccounts = read_xlsx("data/PACIFIC_DATA_2022.xlsx", sheet = 1) %>% subset(Region == "Pacific")
# Convert pairs to individuals
ccounts$Count[ccounts$Count_Type == "Pairs"] = ccounts$Count[ccounts$Count_Type == "Pairs"]*2
ccounts$Count_Type <- "Individuals"

spp_vec_ccounts = unique(ccounts$Species)

# --------------------
# Plot-based counts
# --------------------

pcounts = read_xlsx("data/PACIFIC_DATA_2022.xlsx", sheet = 2) %>% subset(Region == "Pacific") %>%
  subset(Species != "unidentified")
pcounts$Count <- as.numeric(pcounts$Count)
pcounts <- na.omit(pcounts)

# Temporarily omit STPE (not sure if plots are appropriate for this species)
pcounts <- subset(pcounts, Species != "STPE")

pcounts <- subset(pcounts,
                   (Species == "ANMU" & Plot_Type == "spplots: anmuplots") |
                     (Species == "CAAU" & Plot_Type == "spplots: caauplots") |
                     (Species == "RHAU" & Plot_Type == "spplots: rhauplots") |
                     (Species == "TUPU" & Plot_Type == "spplots: tupuplots") )

spp_vec_pcounts = unique(pcounts$Species)

# ------------------------------------------------
# Code for Bayesian analysis
# ------------------------------------------------

# The jags script to fit the model
sink("scripts/jags_files/seabird_gam.jags")
cat("
    model {
  
  # --------------------------------
  # Colony-level intercepts (note that GAM structure includes intercept for colony 1)
  # --------------------------------
  
  C[1] <- 0
  for(k in 2:ncolony){
    C[k] ~ dnorm(0,0.01) 
  }
  
  # --------------------------------
  # GAM smooth
  # --------------------------------
  
  taubeta <- pow(sdbeta,-2) # prior on precision of gam coefficients
  sdbeta ~ dunif(0,5)
  
  nk1 <- nknots-1
  nk2 <- ((nknots*2)-2)
  B.X[1] ~ dnorm(0,0.01)
  
  ## prior for s(year)... 
  K1 <- S1[1:nk1,1:nk1] * lambda[1]  + S1[1:nk1,(nknots:nk2)] * lambda[2]
  B.X[(2:nknots)] ~ dmnorm(zero[(2:nknots)],K1) 
  
  #K1 is the prior on the precisions of the mnorm B.X values (mean GAM parameters for a species)
  ## smoothing parameter
  
  for(i in 1:2) {
    lambda[i] ~ dgamma(0.05,0.005)
    rho[i] <- log(lambda[i])
  } # i
  
  for(j in 1:nknots){ # Computation of GAM components
    
    for(k in 1:ncolony){
      beta.X[k,j] ~ dnorm(B.X[j],taubeta)
      for ( i in 1:nyearspred ){
        X.part[i,j,k] <- beta.X[k,j]*(X[i,j])
      } # i
      
    } # k
  } # j
  
  # --------------------------------
  # Combine pieces
  # --------------------------------
  
  # Residual overdispersion
  sdnoise ~ dunif(0,2)
  taunoise <- pow(sdnoise,-2)
  
  for (i in 1:nyearspred){
    for(k in 1:ncolony){
      yeareffect[i,k] <- sum(X.part[i,1:nknots,k])
      eta[i,k] <- C[k] + yeareffect[i,k]
      log_mu[i,k] ~ dnorm(eta[i,k],taunoise)
    } # k
  } # i
  
  # --------------------------------
  # Likelihood
  # --------------------------------
  
  for (i in 1:ncounts) { 
    count[i] ~ dpois(exp(log_mu[year[i],colony[i]]))                              
  }
  
}
    ",fill = TRUE)
sink()

# ------------------------------------------------
# Loop through species and conduct analysis
# ------------------------------------------------

spp_vec <- union(spp_vec_ccounts, spp_vec_pcounts)
spp_vec <- spp_vec[spp_vec %!in% c("BLOY","COMU","BRCO","HOPU")]

for (spp in spp_vec){
  
  spdat_ccount = spdat_pcount = NULL
  
  # Extract relevant data
  spdat_ccount = subset(ccounts, Species == spp & Year >= 1970)
  spdat_pcount = subset(pcounts, Species == spp & Year >= 1970)
  
  # If there are plot-based counts for an island, use those.  Otherwise, use "complete" colony counts
  if (nrow(spdat_pcount)>0) spdat_ccount <- subset(spdat_ccount, Colony_Name %!in% spdat_pcount$Colony_Name)
  
  spdat_ccount$Survey_Type <- "Colony Counts"
  spdat_pcount$Survey_Type <- "Plot Counts"
  spdat <- dplyr::bind_rows(spdat_ccount,spdat_pcount) %>%
    
    # Remove data prior to 1970
    subset(Year >= 1970)
  
  colonies_to_include = spdat %>%
    group_by(Colony_Name) %>%
    summarize(mean_count = mean(Count),
              min_count = min(Count),
              max_count = max(Count),
              first_survey = min(Year),
              last_survey = max(Year),
              n_surveys = length(unique(Year)),
              Survey_Type = Survey_Type[1]) %>%
    
    # Limit analysis to colonies that were surveyed at least twice, 
    # and within the first and last 10 years of surveys
    subset(n_surveys >= 3 & first_survey <= 1994 & last_survey >= 2010)
  
  if (nrow(colonies_to_include) == 0) next
  
  spdat <- subset(spdat, Colony_Name %in% colonies_to_include$Colony_Name)

  first_year = min(spdat$Year)
  final_year = max(spdat$Year)
  
  # ------------------------------------------------
  # Prepare data for analysis
  # ------------------------------------------------
  
  # Tables that link colony numbers to colony names
  spdat$colony_numeric <- as.integer(factor(spdat$Colony_Name))
  colony_name_table = unique(spdat[,c("Colony_Name","colony_numeric")])
  
  # Tables that link year index to actual year
  year_table = data.frame(Year = min(spdat$Year):2021)
  year_table$year_numeric = 1:nrow(year_table)
  
  # Data for import into jags
  nknots = 6
  year <- spdat$Year - min(year_table$Year) + 1
  ymax <- max(year_table$year_numeric)
  nyears = length(1:ymax)
  colony = spdat$colony_numeric
  ncolony <- max(colony)
  count <- spdat$Count
  ncounts = length(count)
  
  # Use jagam to prepare basis functions
  nyearspred = length(1:ymax)
  preddat = data.frame(yrs = 1:ymax,count = 1)
  form = as.formula(paste("count ~ s(yrs,k =",nknots,")"))
  gamprep = jagam(formula = form,
                  data = preddat,
                  file = "scripts/jags_files/tempgam.txt",
                  centred = T)
  
  # Package data into a list for JAGS
  jags_data = list(X = gamprep$jags.data$X,
                   S1 = gamprep$jags.data$S1,
                   zero = gamprep$jags.data$zero,
                   colony = colony,
                   ncounts = ncounts,
                   ncolony = ncolony,
                   count = count,
                   nknots = nknots,
                   nyearspred = nyearspred,
                   year = year)
  
  # Fit model using JAGS
  parameters.to.save = c("sdnoise","sdbeta","C","yeareffect","log_mu","eta")
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter = 20000,
              n.burnin = 10000,
              n.thin = 10,
              model.file = "scripts/jags_files/seabird_gam.jags",
              n.chains = 3,
              parallel = TRUE)
  
  # Save results for each species
  save.image(paste0("output/Pacific/model_results/",spp,".RData"))
  
}

