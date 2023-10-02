# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------
my_packs = c('tidyverse','readxl','RColorBrewer','viridis','jagsUI','mgcv','ggrepel','scales')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())
theme_set(theme_bw())
setwd("~/iles_ECCC/Seabirds/Status_of_Birds_seabird_analysis_2022/")
`%!in%` <- Negate(`%in%`)

# ------------------------------------------------
# Read in data
# ------------------------------------------------

atl = read_xlsx("data/ATLANTIC_QC_AR_DATA_2022.xlsx", sheet = 1)

# Convert pairs to individuals
atl$Count[atl$Count_Type == "Pairs"] = atl$Count[atl$Count_Type == "Pairs"]*2
atl$Count_Type[atl$Count_Type == "Pairs"] = "Individuals"
atl = subset(atl, Count_Type == "Individuals") # Removes several AOS rows

atl_summary = atl %>%
  group_by(Species,Colony_Name,Count_Type) %>%
  summarize(first_year = min(Year),
            last_year = max(Year),
            n_counts = length(unique(Year)),
            min_count = min(Count),
            max_count = max(Count),
            mean_count = mean(Count))

write.csv(atl_summary, "output/Atlantic/dataviz/0_Atlantic_summary.csv",row.names = FALSE)

atl_species = unique(atl$Species) %>% sort()

# ------------------------------------------------
# Code for Bayesian analysis
# ------------------------------------------------

# The jags script to fit the model
sink("scripts/jags_files/seabird_gam.jags")
cat("
    model {
  
  # --------------------------------
  # Fixed effects and observation
  # --------------------------------
  
  # Colony-level intercepts
  C[1] <- 0
  for(k in 2:ncolony){
    C[k] ~ dnorm(0,0.01) 
  }
  
  sdnoise ~ dunif(0,2)
  taunoise <- pow(sdnoise,-2)
  
  for (i in 1:ncounts) { 
    count[i] ~ dpois(mu[i])                                 # response  
    mu[i] <-  exp(eta1[i])                                  # expected response
    eta1[i] ~ dnorm(eta[i],taunoise)             # overdispersion
    eta[i] <- C[colony[i]] + yeareffect[year[i],colony[i]]  # process model - colony intercept and year-smoother
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
  
  for (i in 1:nyearspred){
    for(k in 1:ncolony){
      yeareffect[i,k] <- sum(X.part[i,1:nknots,k])
    } # k
  } # i
  
  # --------------------------------
  # Derived parameters: predictions of annual expected counts
  # --------------------------------
  
  for (i in 1:nyearspred){
    for(j in 1:nknots){
      X.partpred[i,j] <- B.X[j]*(X[i,j])
    } #j
    
    x.gampred[i] <- sum(X.partpred[i,1:nknots])
    
    for(k in 1:ncolony){
      etapred[i,k] <- exp(C[k] + yeareffect[i,k])
    }
    
  } #i
}
    ",fill = TRUE)
sink()

# ------------------------------------------------
# Loop through species and conduct analysis
# ------------------------------------------------

for (spp in atl_species){

  # Extract relevant data
  spdat = subset(atl, Species == spp)

  # Limit analysis to colonies that were surveyed at least twice, and within the first and last 10 years of surveys
  first_year = min(spdat$Year)
  final_year = max(spdat$Year)

  colonies_to_include = spdat %>%
    group_by(Colony_Name) %>%
    summarize(mean_count = mean(Count),
              first_survey = min(Year),
              last_survey = max(Year),
              n_surveys = length(unique(Year)))
  
  # ------------------------------------------------
  # Some custom data cleaning for certain species
  # ------------------------------------------------
  
  if (spp != "ROST" & spp != "NOFU" & spp != "CATE"){
    # Remove the smallest colonies, which are strongly inflating uncertainty
    colonies_to_include = colonies_to_include %>% subset(mean_count > 100)
  }
  
  #if (spp == "RBGU"){
  #  # Omit the smallest RBGU colonies that are strongly influencing uncertainty
  #  colonies_to_include = subset(colonies_to_include, mean_count > 1000)
  #}
  
  if (spp == "LESP"){
     # Omit the smallest LESP colonies that are strongly influencing uncertainty
     colonies_to_include = subset(colonies_to_include, mean_count > 1000)
  }

  spdat = subset(spdat, Colony_Name %in% colonies_to_include$Colony_Name)

  # Skip species if there is not enough data
  if (nrow(spdat) == 0 | length(unique(spdat$Colony_Name)) < 2) next

  # Tables that link colony numbers to colony names
  spdat$colony_numeric <- as.integer(factor(spdat$Colony_Name))
  colony_name_table = unique(spdat[,c("Colony_Name","colony_numeric")])

  # Tables that link year index to actual year
  year_table = data.frame(Year = 1950:2021)
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
  parameters.to.save = c("testnoise","sdnoise","sdbeta","C","beta.X","etapred")
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter = 200000,
              n.burnin = 100000,
              n.thin = 100,
              model.file = "scripts/jags_files/seabird_gam.jags",
              n.chains = 3,
              parallel = TRUE)

  # Save results for each species
  save.image(paste0("~/iles_ECCC/Seabirds/Status_of_Birds_seabird_analysis_2022/output/Atlantic/model_results/",spp,".RData"))

  # -----------------------------------------------------
  # Plot colony trajectories
  # -----------------------------------------------------
  
  # Extract predictions in dataframe format
  N_samples = reshape2::melt(out$sims.list$etapred) %>%
    rename(samp = Var1, year_number = Var2, colony_number = Var3, N_pred = value) %>%
    full_join(colony_name_table, by = c("colony_number" = "colony_numeric")) %>%
    full_join(year_table, by = c("year_number" = "year_numeric")) %>%
    add_column(Species = spp)
  
  N_summary_colony = N_samples %>% group_by(Colony_Name, Year) %>%
    summarize(q05 = quantile(N_pred,0.05),
              q500 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q95 = quantile(N_pred,0.95))
  
  colors = viridis(10)[c(3,6,9)]
  
  # Dynamics from 1970 onwards
  colony_plot_1970 = ggplot() +
    geom_ribbon(data = subset(N_summary_colony, Year >= 1970), aes(x = Year, ymin = q05, ymax = q95), fill = "gray85", col = "transparent")+
    geom_line(data = subset(N_summary_colony, Year >= 1970), aes(x = Year, y = q500), col = "black")+
    geom_vline(data = subset(spdat, Count == 0 & Year >= 1970),aes(xintercept = Year), col = "red", size = 1.5, alpha = 0.3)+
    geom_point(data = subset(atl, Species == spp & Year >= 1970),aes(x = Year, y = Count, col = Count_Type))+
    xlab("Year")+
    ylab("Population index")+
    scale_color_manual(values = colors,name = "Count type", drop = FALSE)+
    facet_wrap(Colony_Name~., scales = "free_y")+
    ggtitle(paste0(spp," colony-level trajectories"))+
    scale_y_continuous(labels = comma)
  
  print(colony_plot_1970)
  
  N_samples_regional = N_samples %>% group_by(Species,samp,Year) %>%
    summarize(N_pred = sum(N_pred))
  
  N_summary_regional = N_samples_regional %>%
    group_by(Species,Year) %>%
    summarize(q05 = quantile(N_pred,0.05),
              q500 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q95 = quantile(N_pred,0.95))
  
  # Dynamics from 1970 onwards
  regional_plot_1970 = ggplot() +
    geom_ribbon(data = subset(N_summary_regional, Year >= 1970), aes(x = Year, ymin = q05, ymax = q95), fill = "gray85", col = "transparent")+
    geom_line(data = subset(N_summary_regional, Year >= 1970), aes(x = Year, y = q500), col = "black")+
    xlab("Year")+
    ylab("Population index")+
    scale_y_continuous(labels = comma)+
    ggtitle(paste0(spp," regional trajectory"))
  
  regional_plot_1970 
}

