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
    eta1[i] ~ dnorm(eta[i],taunoise)                        # overdispersion
    eta[i] <- C[colony[i]] + yeareffect[year[i],colony[i]]  # process model - colony intercept and year-smoother
  }
  
  # --------------------------------
  # GAM smooth
  # --------------------------------
  
  taubeta <- pow(sdbeta,-2) # prior on precision of gam coefficients
  sdbeta ~ dunif(0,1)
  
  nk1 <- nknots-1
  nk2 <- ((nknots*2)-2)
  # B.X[1] ~ dnorm(0,0.01)
  B.X[1] <- 0
  
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
    
    for(k in 1:ncolony){
      etapred[i,k] <- exp(C[k] + yeareffect[i,k])
    }
    
    etasum[i] <- sum(etapred[i,1:ncolony])
    
  } #i
  
  # Overall trend
  
}
    ",fill = TRUE)
sink()

# ------------------------------------------------
# Simulate different trajectories at a series of colonies
# ------------------------------------------------

sim_results <- colony_sim_results <- data.frame()

for (sim_run in 1:250){
  
  ncolony <- 10
  nyears <- 50
  
  spdat <- expand.grid(year = 1:nyears,colony = 1:ncolony, log_eta = NA, mu=NA, y = NA)
  
  sigma <- 0.15
  colony_intercepts <- runif(ncolony,1000,50000)
  
  colony_betas <- matrix(NA,nrow = ncolony,ncol = 3)
  colony_betas[,1] <- rnorm(ncolony,0,0.01)
  colony_betas[,2] <- rnorm(ncolony,-0.0002,0.0001)
  colony_betas[,3] <- runif(ncolony,-0.000025,0.00001)
  
  for (i in 1:nrow(spdat)){
    j = spdat$colony[i]
    spdat$log_eta[i] <- log(colony_intercepts[j]) + colony_betas[j,1]*spdat$year[i] + colony_betas[j,2]*spdat$year[i]^2 + colony_betas[j,3]*spdat$year[i]^3
    spdat$mu[i] <- exp(spdat$log_eta[i] + rnorm(1,0,sigma))
    spdat$y[i] <- rpois(1,spdat$mu[i])
  }
  
  # Remove 75% of data
  #to_remove <- sample(1:nrow(spdat),round(nrow(spdat)*0.75))
  #spdat$y[to_remove] <- NA
  
  ggplot(spdat)+
    geom_line(aes(x = year, y = exp(log_eta)))+
    geom_point(aes(x = year, y = y))+
    facet_wrap(colony~., scales = "free")
  
  # Total index is assumed to be the sum of annual expected (median) counts
  regional_total <- spdat %>%
    group_by(year) %>%
    summarize(Nsum = sum(exp(log_eta)))
  
  # Data for import into jags
  nknots = 6
  year <- spdat$year
  ymax <- nyears
  nyears = nyears
  colony = spdat$colony
  ncolony <- ncolony
  count <- as.numeric(spdat$y)
  ncounts = nrow(spdat)
  
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
  parameters.to.save = c("sdnoise","sdbeta","C","beta.X","etapred")
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter = 10000,
              n.burnin = 5000,
              n.thin = 5,
              model.file = "scripts/jags_files/seabird_gam.jags",
              n.chains = 3,
              parallel = TRUE)
  
  # -----------------------------------------------------
  # Plot colony trajectories
  # -----------------------------------------------------
  
  # Extract predictions in dataframe format
  N_samples = reshape2::melt(out$sims.list$etapred) %>%
    rename(samp = Var1, year = Var2, colony = Var3, N_pred = value) 
  
  N_summary_colony = N_samples %>% group_by(colony,year) %>%
    summarize(q05 = quantile(N_pred,0.05),
              q500 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q95 = quantile(N_pred,0.95))
  
  colors = viridis(10)[c(3,6,9)]
  
  # Dynamics from 1970 onwards
  colony_plot = ggplot() +
    geom_ribbon(data = N_summary_colony, aes(x = year, ymin = q05, ymax = q95), fill = "gray85", col = "transparent")+
    geom_line(data = spdat,aes(x = year, y = exp(log_eta)))+
    geom_point(data = spdat,aes(x = year, y = y))+
    xlab("Year")+
    ylab("Population index")+
    facet_wrap(colony~., scales = "free_y")+
    scale_y_continuous(labels = comma)
  
  N_samples_regional = N_samples %>% 
    group_by(samp,year) %>%
    summarize(N_pred = sum(N_pred))
  
  N_summary_regional = N_samples_regional %>%
    group_by(year) %>%
    summarize(q05 = quantile(N_pred,0.05),
              q50 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q95 = quantile(N_pred,0.95))
  
  # Regional dynamics
  regional_plot = ggplot() +
    geom_ribbon(data = N_summary_regional, aes(x = year, ymin = q05, ymax = q95), fill = "gray85", col = "transparent")+
    geom_line(data = N_summary_regional, aes(x = year, y = q50), col = "black")+
    geom_line(data = regional_total, aes(x = year, y = Nsum), col = "blue")+
    xlab("Year")+
    ylab("Population index")+
    scale_y_continuous(labels = comma)
  
  # Trend comparison
  trend_true <- log(regional_total$Nsum[nyears]/regional_total$Nsum[1])
  trend_samples <- log(subset(N_samples_regional,year == nyears)$N_pred/subset(N_samples_regional,year == 1)$N_pred)
  trend_est_mean <- mean(trend_samples)
  trend_est_median <- median(trend_samples)
  trend_est_lcl <- quantile(trend_samples,0.025) %>% as.numeric()
  trend_est_ucl <- quantile(trend_samples,0.975) %>% as.numeric()
  
  run_result <- data.frame(true = trend_true,
                           est_mean = trend_est_mean,
                           est_median = trend_est_median,
                           lcl = trend_est_lcl,
                           ucl = trend_est_ucl,
                           bias = trend_est_median - trend_true,
                           cov = trend_est_lcl < trend_true & trend_est_ucl > trend_true)
  
  sim_results <- rbind(sim_results,run_result)
  hist(sim_results$bias, main = paste0("Mean bias = ",round(mean(sim_results$bias),2),"\nCoverage = ",round(mean(sim_results$cov),2)))
  
  print(sim_run)
  
  # Colony-level trends
  colony_trend_samples <- log(out$sims.list$etapred[,nyears,]/out$sims.list$etapred[,1,]) %>%
    reshape2::melt() %>%
    rename(samp = Var1, colony = Var2, trend = value)
  
  colony_trend_estimates <- colony_trend_samples %>%
    group_by(colony) %>%
    summarize(trend_mean = mean(trend),
              trend_med = median(trend))
  
  colony_trend_true <- subset(spdat, year == nyears)$log_eta - subset(spdat, year == 1)$log_eta
  
  colony_trend_biases <- data.frame(sim_run = sim_run,
                                    colony = 1:ncolony,
                                    bias = colony_trend_estimates$trend_mean - colony_trend_true)
  
  colony_sim_results <- rbind(colony_sim_results,colony_trend_biases)
  
}

# Is bias induced by summation, or is it present at each colony too?

mean(colony_sim_results$bias)
