
my_packs = c('tidyverse','readxl','RColorBrewer','viridis','jagsUI','mgcv','ggrepel','scales','ggthemes')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

stub <- function() {}
thisPath <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  if (length(grep("^-f$", cmdArgs)) > 0) {
    # R console option
    normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  } else if (length(grep("^--file=", cmdArgs)) > 0) {
    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
  } else if (Sys.getenv("RSTUDIO") == "1") {
    # RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {
    # 'source'd via R console
    dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {
    stop("Cannot find file path")
  }
}

dirname <- thisPath()
setwd(dirname)

`%!in%` <- Negate(`%in%`)

# ------------------------------------------------
# ggplot theme
# ------------------------------------------------

CustomTheme <- theme_set(theme_bw())
CustomTheme <- theme_update(legend.key = element_rect(colour = NA), 
                            legend.key.height = unit(1.2, "line"),
                            panel.grid.major = element_line(colour = 'transparent'),
                            panel.grid.minor = element_line(colour = 'transparent'),
                            panel.border = element_rect(linetype = "solid",
                                                        colour = "black",
                                                        size = 1, fill = NA),
                            axis.line = element_line(colour = "black"),
                            strip.text = element_text(size = 12, colour = "black"),
                            strip.background = element_rect(colour = "black",
                                                            fill = "lightblue2",
                                                            linetype = "solid"),
                            axis.title.y = element_text(margin = margin(0,10,0,0)),
                            axis.title.x = element_text(margin = margin(10,0,0,0)),
                            panel.background = element_rect(fill = "white"))

# ----------------------------------------------------------
# Part 1: Simulate an entire 50-year time series for each of 9 colonies
# ----------------------------------------------------------

trend_results <- data.frame()
for (run in 1:100){

  set.seed(run)
  ncolony <- 9
  nyears <- 50
  
  N_matrix <- matrix(NA,nrow=ncolony,ncol = nyears)
  
  # Initial abundance
  N_matrix[,1] <- rlnorm(ncolony, meanlog = log(10000), sdlog = 2) %>% sort()
  
  # Additive annual variance different among colonies
  sd_additive <- runif(ncolony,0,0.1) #%>% sort() %>% rev()
  
  # Set a colony-specific variance in annual growth rates
  process_sd = runif(ncolony,0,0.3) #%>% sort() %>% rev()
  
  # Magnitude of observation error
  obs_sd = 0.2
  
  # Simulate trajectories at each colony
  for (i in 1:ncolony){

    # ----------------------------------------------------------
    # Code to generate log-linear trajectories
    # ----------------------------------------------------------
    
    trend <- rnorm(1,0,0.03)
    
    # Generate a random trajectory at the colony (random walk)
    for (t in 2:nyears) N_matrix[i,t] <- exp(log(N_matrix[i,t-1]) + trend)
    
    # ----------------------------------------------------------
    # Code to generate complex random walks, where some colonies 'switch' their trend part way through the simulation
    # ----------------------------------------------------------
    
    # # switch dynamics for some colonies
    # switch = sample(0:1,1)
    # 
    # r_mean = runif(1,-0.02,0.01)
    # r_switch <- runif(1,-0.1,0.05)
    # 
    # # Generate a random trajectory at the colony (random walk)
    # for (t in 2:nyears){
    #   N_matrix[i,t] <- exp(log(N_matrix[i,t-1]) + rnorm(1,mean=r_mean, sd = process_sd[i]))
    #   if (switch == 1 & (t > nyears/2))N_matrix[i,t] <- exp(log(N_matrix[i,t-1]) + rnorm(1,mean=r_mean+r_switch, sd = process_sd[i]))
    # }
    # 
    
    # ----------------------------------------------------------
    # Additive log-normal error layered on top
    # ----------------------------------------------------------
    
    N_matrix[i,] <- exp(rnorm(nyears,log(N_matrix[i,]),sd_additive[i]))
  }
  
  # Convert to dataframe (to plot with ggplot)
  N_df <- reshape2::melt(N_matrix) %>% 
    rename(Colony = Var1,Year=Var2,N=value)
  
  # Plot dynamics (on log10 scale)
  ggplot()+
    geom_line(aes(x = 1:nyears, y = colSums(N_matrix), col = "Regional Sum"),linewidth = 1)+
    geom_line(data = N_df, aes(x = Year, y = N, col = factor(Colony)))+
    theme_few()+
    scale_y_continuous(labels = comma, trans = "log10")+
    scale_color_manual(values = c(brewer.pal(ncolony,"Spectral"),"black"), name = "Colony")+
    ggtitle("Simulated trajectories at each of 9 colonies")+
    ylab("Abundance")+
    xlab("Year")
  
  ggplot()+
    geom_line(aes(x = 1:nyears, y = colSums(N_matrix), col = "Regional Sum"),linewidth = 1)+
    #geom_line(data = N_df, aes(x = Year, y = N, col = factor(Colony)))+
    theme_few()+
    scale_y_continuous(labels = comma, trans = "log10")+
    #scale_color_manual(values = c(brewer.pal(ncolony,"Spectral"),"black"), name = "Colony")+
    ggtitle("Simulated trajectories at each of 9 colonies")+
    ylab("Abundance")+
    xlab("Year")
  
  # ----------------------------------------------------------
  # Part 2: Simulate intermittent surveys (3-5 surveys at each colony)
  # ----------------------------------------------------------
  
  N_df$lambda_obs <- rlnorm(nrow(N_df),meanlog = log(N_df$N) - 0.5*obs_sd^2, sdlog = obs_sd)
  N_df$SurveyCount <- NA
  
  for (i in 1:ncolony){
    
    survey_years <- c()
    
    # Simulate one count in first 10 years of surveys
    survey_years <- c(survey_years, sample(1:10,1))
    
    # Simulate one count in final 10 years of surveys
    survey_years <- c(survey_years, sample(nyears:(nyears-10),1))
    
    # Simulate 0-4 additional surveys
    survey_years <- c(survey_years, sample(11:(nyears-11),sample(0:4,1)))
  
    # Poisson observations
    N_df$SurveyCount[(N_df$Colony == i) & 
                       (N_df$Year %in% survey_years)] <- rpois(n = length(survey_years), 
                                                               lambda = N_df$lambda_obs[(N_df$Colony == i) & (N_df$Year %in% survey_years)])
    
  }
  
  # Plot survey counts at each of 5 colonies
  ggplot()+
    geom_line(data = N_df, aes(x = Year, y = N), col = "gray80")+
    geom_point(data = N_df, aes(x = Year, y = SurveyCount), col = "black")+
    theme_few()+
    facet_wrap(Colony~., scales = "free_y")+
    scale_y_continuous(trans="log10", labels = comma)+
    ggtitle("Simulated surveys at each of 9 colonies")
  
  # ----------------------------------------------------------
  # Part 3: Fit model to simulated survey data
  # ----------------------------------------------------------
  
  spdat <- na.omit(N_df) %>% dplyr::select(Colony,Year,SurveyCount)
  
  # Data for import into jags
  nknots = 4
  year <- spdat$Year
  ymax <- nyears
  colony = spdat$Colony
  count <- spdat$SurveyCount
  ncounts = length(count)
  
  # Use jagam to prepare basis functions
  nyearspred = length(1:ymax)
  preddat = data.frame(yrs = 1:ymax,count = 1)
  form = as.formula(paste("count ~ s(yrs,k =",nknots,")"))
  gamprep = jagam(formula = form,
                  data = preddat,
                  file = "tempgam.txt",
                  centred = T)
  
  # Basis functions (just to visualize)
  #gam_example <- gam(count~s(yrs,k = nknots,bs = "tp"),data = preddat,centred=T)
  #model_matrix <- predict(gam_example, type = "lpmatrix")
  #matplot(preddat$yrs, model_matrix[,-1], type = "l", lty = 2)
  
  
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
  parameters.to.save = c("sdnoise","sdbeta","C","beta.X","population_index")
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter = 110000,
              n.burnin = 10000,
              n.thin = 10,
              model.file = "seabird_GAMM.jags",
              n.chains = 3,
              parallel = TRUE)
  
  out$mcmc.info$elapsed.mins # 0.25 mins
  
  # ----------------------------------------------------------
  # Part 4: Summarize predictions and compare to true (i.e., simulated) trajectories
  # ----------------------------------------------------------
  
  # Extract predictions in dataframe format
  fit_samples = reshape2::melt(out$sims.list$population_index) %>%
    rename(samp = Var1, Year = Var2, Colony = Var3, N_pred = value)
  
  N_summary_colony = fit_samples %>% 
    group_by(Colony, Year) %>%
    summarize(q025 = quantile(N_pred,0.025),
              q50 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q975 = quantile(N_pred,0.975)) %>%
    
    # Join with true values
    full_join(N_df)
  
  # Plot estimates
  ggplot()+
    geom_ribbon(data = N_summary_colony, aes(x = Year, ymin = q025, ymax = q975), alpha = 0.2, fill = "dodgerblue")+
    geom_line(data = N_summary_colony, aes(x = Year, y = q50, col = "Estimate"))+
    
    geom_line(data = N_df, aes(x = Year, y = N, col = "True Trajectory"))+
    geom_point(data = N_df, aes(x = Year, y = SurveyCount, col = "Observed Count"))+
    theme_few()+
    facet_wrap(Colony~., scales = "free")+
    scale_y_continuous(trans="log10", labels = comma)+
    ggtitle("Simulated trajectories at each of 9 colonies")+
    scale_color_manual(values = c("dodgerblue","black","red"), name = "")+
    ylab("Index of abundance")
  
  # ----------------------------------------------------------
  # Part 5: Calculate regional total
  # ----------------------------------------------------------
  
  # True regional total
  regional_df <- N_df %>%
    group_by(Year) %>%
    summarize(N = sum(N))
  
  # Estimated regional total
  regional_samples = fit_samples %>% 
    group_by(Year,samp) %>%
    summarize(N_pred = sum(N_pred)) 
  
  # Summary (mean and 95% CI)
  N_summary_regional <- regional_samples %>%
    group_by(Year) %>%
    summarize(q025 = quantile(N_pred,0.025),
              q50 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q975 = quantile(N_pred,0.975)) %>%
    
    # Join with true values
    full_join(regional_df)
  
  # Trend estimate
  baseline_year <- 1
  trend_true <- 100*((regional_df$N[regional_df$Year == nyears]/regional_df$N[regional_df$Year == baseline_year])^(1/(nyears-baseline_year))-1)
  
  trend_est <- 100*((regional_samples$N_pred[regional_samples$Year == nyears]/regional_samples$N_pred[regional_samples$Year == baseline_year])^(1/(nyears-baseline_year))-1)
  trend_est <- quantile(trend_est,c(0.025,0.5,0.975))
  
  ggplot()+
    geom_ribbon(data = N_summary_regional, aes(x = Year, ymin = q025, ymax = q975), alpha = 0.2, fill = "dodgerblue")+
    geom_line(data = N_summary_regional, aes(x = Year, y = q50, col = "Estimate"))+
    
    geom_line(data = regional_df, aes(x = Year, y = N, col = "True Trajectory"))+
    theme_few()+
    scale_y_continuous(labels = comma)+
    ggtitle("Regional trajectory")+
    scale_color_manual(values = c("dodgerblue","black","red"), name = "")+
    ylab("Index of abundance")+
    geom_text(aes(x = 0, 
                  y = max(c(N_summary_regional$q975,N_summary_regional$N))), 
              label = paste0("True trend = ",round(trend_true,2),"% per year\nEst trend = ",round(trend_est[2],2),"% (",round(trend_est[1],2)," to ",round(trend_est[3],2),")"), hjust=0)
    
  # ----------------------------------------------------------
  # Append results for this simulation run to dataframe
  # ----------------------------------------------------------
  
  trend_results <- rbind(trend_results,data.frame(run = run,
                                                  trend_true = trend_true,
                                                  trend_est_q025 = trend_est[1],
                                                  trend_est_q500 = trend_est[2],
                                                  trend_est_q975 = trend_est[3],
                                                  cov = trend_true > trend_est[1] & trend_true < trend_est[3],
                                                  max_Rhat = max(out$Rhat$population_index)))
  
  # ----------------------------------------------------------
  # Plot results
  # ----------------------------------------------------------
  
  lim = range(trend_results[,c("trend_true","trend_est_q025","trend_est_q975")])
  lim = c(-6,6)
  trend_plot <- ggplot(data = trend_results, aes(x = trend_true, y = trend_est_q500, ymin = trend_est_q025, ymax = trend_est_q975,col=cov))+
    geom_abline(intercept=0,slope=1,col="gray85")+
    geom_errorbar(width=0)+
    geom_point()+
    coord_cartesian(ylim=lim,xlim=lim)+
    theme_bw()+
    xlab("True (simulated) regional trend")+
    ylab("Estimated regional trend")+
    scale_color_manual(values=c("red","dodgerblue"), name = "Coverage")+
    ggtitle("Fit with GAMM")
  print(trend_plot)
  
}

# ----------------------------------------------------------
# Summarize results across repeated simulations
# ----------------------------------------------------------

mean(trend_results$cov) # coverage
mean(trend_results$trend_est_q500 - trend_results$trend_true) # accuracy
mean(trend_results$trend_est_q975 - trend_results$trend_est_q025) # precision

lim = c(-6,6)
trend_plot <- ggplot(data = trend_results, aes(x = trend_true, y = trend_est_q500, ymin = trend_est_q025, ymax = trend_est_q975,col=cov))+
  geom_abline(intercept=0,slope=1,col="gray85")+
  geom_errorbar(width=0)+
  geom_point()+
  coord_cartesian(ylim=lim,xlim=lim)+
  theme_bw()+
  xlab("True (simulated) regional trend")+
  ylab("Estimated regional trend")+
  scale_color_manual(values=c("red","dodgerblue"), name = "Coverage")+
  ggtitle("Fit with GAMM")
print(trend_plot)