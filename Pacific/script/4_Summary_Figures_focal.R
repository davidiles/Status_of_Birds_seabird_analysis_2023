# Custom analysis summary/figures for Laurie Wilson

# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs = c('tidyverse','readxl','RColorBrewer','viridis','jagsUI','mgcv','ggrepel','scales')

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
setwd("../")

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

# ------------------------------------------------
# Read in data
# ------------------------------------------------

# Complete counts
ccounts = read_xlsx("data/PACIFIC_DATA_2023.xlsx", sheet = 1) %>% subset(Region == "Pacific")
spp_vec = unique(ccounts$Species)

# ------------------------------------------------
# Prepare species codes / labels / generation lengths
# ------------------------------------------------

GenLength <- read_xlsx("../other_files/Bird_et_al_2020_AppendixS4_GenLength.xlsx") %>%
  select('Scientific name','GenLength')
SOBC_species <- read_xlsx("../other_files/SOCB_species.xlsx")

seabird_codes_names <- read_xlsx("../other_files/seabird_names_2023.xlsx") %>%
  left_join(SOBC_species) %>%
  left_join(GenLength, by = c('scientific_name' = 'Scientific name')) %>%
  mutate(GenLength = round(GenLength))

seabird_codes_names %>% as.data.frame()

# Empty dataframes to store results
Trend_Estimates <- data.frame()
Annual_Indices <- data.frame()
Colony_Indices <- data.frame()

spp_vec = c("ANMU","CAAU","RHAU","TUPU")

for (spp in spp_vec){
  
  # Load model results
  file = paste0("output/model_results/",spp,".RData")
  if(!file.exists(file)) next 
  load(file)
  
  # ------------------------------------------------
  # Most recent colony-level abundance estimate for plot-based surveys
  # (needed to appropriately "weight" plot-based estimates at regional scale)
  # ------------------------------------------------
  
  plot_based_colonies <- unique(spdat_pcount$Colony_Name)
  
  # Read in 'total colony abundance' estimates
  pabund = read_xlsx("data/PACIFIC_DATA_2023.xlsx", sheet = 3) %>%
    subset(Species == spp) %>%
    group_by(Colony_Name) %>%
    
    # Select most recent "total abundance" estimate
    arrange(Year_of_Size_Estimate)%>%
    filter(row_number()==n())
  
  pabund$Total_Size_of_Colony[pabund$Comments == "Pairs"] <- pabund$Total_Size_of_Colony[pabund$Comments == "Pairs"]*2
  pabund <- pabund %>% 
    select(Species, Colony_Name, Total_Size_of_Colony,Year_of_Size_Estimate) %>%
    rename(Count = Total_Size_of_Colony)
  
  # ------------------------------------------------
  # Identify reasonable starting point for trend estimates
  # ------------------------------------------------
  
  colony_summary <- spdat %>%
    group_by(Colony_Name,Survey_Type) %>%
    summarize(mean_count = max(Count),
              first_survey = min(Year),
              final_survey = max(Year)) %>%
    arrange(desc(mean_count)) %>%
    ungroup()
  
  # For plot-based surveys, substitute most recent colony-level count
  if (nrow(pabund)>0){
    for (colony in unique(pabund$Colony_Name)){
      colony_summary$mean_count[which(colony_summary$Colony_Name == colony)] <- pabund$Count[which(pabund$Colony_Name == colony)]
    }
  }
  
  colony_summary <- colony_summary %>%
    mutate(include = mean_count >= (sum(mean_count,na.rm = TRUE)*0.2))
  
  # If plot-based counts exist, consider those to be "high quality" surveys as well
  colony_summary$include[colony_summary$Survey_Type == "Plot Counts"] <- TRUE
  
  if (sum(colony_summary$include)==0) colony_summary$include = TRUE
  
  # Omit Frederick Island for ANMU, because survey protocols have changed
  if (spp == "ANMU"){
    colony_summary$include[which(colony_summary$Colony_Name == "Frederick Island (Susk Gwaii)")] <- FALSE
  }
  if (spp == "CAAU"){
    colony_summary$include[which(colony_summary$Colony_Name == "Frederick Island (Susk Gwaii)")] <- FALSE
  }
  
  empirical_start <- floor(mean(subset(colony_summary, include)$first_survey))
  empirical_end <- ceiling(mean(subset(colony_summary, include)$final_survey))
  
  empirical_start <- max(c(empirical_start, min(spdat$Year)))
  empirical_end <- min(c(empirical_end, max(spdat$Year)))
  
  # ------------------------------------------------
  # Decide which colonies to include in regional summary
  # ------------------------------------------------
  
  colonies_to_include_for_regional = colonies_to_include %>%
    subset(n_surveys >= 2 &
             first_survey <= (empirical_start+10) &
             last_survey >= (empirical_end-10))
  colonies_to_include_for_regional
  
  # ------------------------------------------------
  # (Approximate) proportion of population included in regional summary
  # ------------------------------------------------
  
  colony_summary$regional <- FALSE
  colony_summary$regional[colony_summary$Colony_Name %in% colonies_to_include_for_regional$Colony_Name] <- TRUE
  
  proportion <- sum(colony_summary$mean_count[colony_summary$regional],na.rm = TRUE)/sum(colony_summary$mean_count,na.rm = TRUE)
  
  coverage_cat <- "Low"
  if (proportion > 0.25) coverage_cat <- "Medium"
  if (proportion > 0.50) coverage_cat <- "High"
  
  # ------------------------------------------------
  # Generate plots of colony-level trends with observed data overlaid
  # ------------------------------------------------
  eta_samples = reshape2::melt(out$sims.list$population_index) %>%
    rename(samp = Var1, year_number = Var2, colony_number = Var3, N_pred = value) %>%
    full_join(colony_name_table, by = c("colony_number" = "colony_numeric")) %>%
    full_join(year_table, by = c("year_number" = "year_numeric")) %>%
    add_column(Species = spp)
  
  # General summaries
  N_summary_colony = eta_samples %>% 
    group_by(Colony_Name, Year) %>%
    summarize(q025 = quantile(N_pred,0.025),
              q50 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q975 = quantile(N_pred,0.975))
  
  # Colonies monitored by plots
  pcol <- unique(spdat_pcount$Colony_Name)
  
  # pcounts2 <- subset(pcounts,
  #                   (Species == "ANMU" & Plot_Type == "spplots: anmuplots") |
  #                     (Species == "CAAU" & Plot_Type == "spplots: caauplots") |
  #                     (Species == "RHAU" & Plot_Type == "spplots: rhauplots") |
  #                     (Species == "TUPU" & Plot_Type == "spplots: tupuplots") ) %>%
  #   subset(Species == spp)
  
  # Plot 500 samples from the posterior (to visualize spread in trajectories)
  samples_to_plot <- sample(unique(eta_samples$samp),500)
  
  # Colony-level dynamics
  colony_plot_freeaxis = ggplot() +
    
    geom_ribbon(data = subset(N_summary_colony, Colony_Name %in% pcol), aes(x = Year, ymin = q025, ymax = q975),alpha = 0.01, fill = "dodgerblue", col = "black", linetype = 2, linewidth = 0.5)+
    geom_line(data = subset(eta_samples, samp %in% samples_to_plot & Colony_Name %in% pcol), aes(x = Year, y = N_pred, col = factor(samp)),alpha = 0.1)+
    geom_line(data = subset(N_summary_colony, Colony_Name %in% pcol),  aes(x = Year, y = q50), linewidth = 1, col = "black")+
    
    geom_point(data = subset(spdat, Colony_Name %in% pcol),aes(x = Year, y = Count / n_plots, size = n_plots))+
    xlab("Year")+
    ylab("Population index")+
    scale_color_manual(values=rep("dodgerblue",length(unique(eta_samples$samp))), guide = "none")+
    scale_size_continuous(name = "Number of plots surveyed", range = c(1,5))+
    facet_wrap(Colony_Name~., scales = "free_y")+
    scale_y_continuous(labels = comma, trans = "log10")+
    ggtitle(seabird_codes_names$english_name[seabird_codes_names$Species == spp])+
    theme(legend.position="bottom")
  
  #print(colony_plot_freeaxis)
  
  png(paste0("output/model_results/figures/",spp,"_Plot_Surveys.png"), width = 8, height = 6, units = "in", res = 600)
  print(colony_plot_freeaxis)
  dev.off()
  
  # ***********************************************************
  # Re-scale colony-level estimates for plot-based surveys
  # ***********************************************************
  
  # Predictions of expected counts (log scale)
  eta_samples = reshape2::melt(out$sims.list$eta) %>%
    rename(samp = Var1, year_number = Var2, colony_number = Var3, eta = value) %>%
    full_join(colony_name_table, by = c("colony_number" = "colony_numeric")) %>%
    full_join(year_table, by = c("year_number" = "year_numeric")) %>%
    add_column(Species = spp)
  
  # Predictions of expected counts (log scale) + noise
  log_mu_samples = reshape2::melt(out$sims.list$log_mu) %>%
    rename(samp = Var1, year_number = Var2, colony_number = Var3, log_mu = value) %>%
    full_join(colony_name_table, by = c("colony_number" = "colony_numeric")) %>%
    full_join(year_table, by = c("year_number" = "year_numeric")) %>%
    add_column(Species = spp)
  
  # Predictions of noise
  log_sd_samples = reshape2::melt(out$sims.list$sdnoise) %>%
    rename(samp = Var1, colony_number = Var2, log_sd = value) %>%
    full_join(colony_name_table, by = c("colony_number" = "colony_numeric")) %>%
    add_column(Species = spp)
  
  fit_samples <- full_join(eta_samples,log_mu_samples)
  fit_samples$resid <- fit_samples$log_mu - fit_samples$eta
  
  # placeholder for reweighted abundances
  fit_samples$eta_reweighted <- fit_samples$eta
  
  # -------------------------------------------------------
  # Loop through colonies and recalculate estimates of eta
  # -------------------------------------------------------
  
  for (colony in plot_based_colonies){
    
    # Year of most recent abundance estimate ("baseline year")
    byear <- subset(pabund, Colony_Name == colony)$Year_of_Size_Estimate
    
    # If byear is too early, pin it to the earliest year of fit_samples
    if (byear < min(fit_samples$Year)) byear = min(fit_samples$Year)
    
    byear_i <- which(fit_samples$Colony_Name == colony & fit_samples$Year == byear) # indices in N_samples df
    
    # Assume that the abundance estimate represents exp(log_mu), where log_mu = eta + resid
    # We want to focus our inference on eta
    N_est <- pabund$Count[which(pabund$Colony_Name == colony)]
    
    # Let N_est = exp(log_mu). Also let log(N_est) = log_mu = eta + resid. Therefore, eta = log(N_est) - resid
    fit_samples$eta_reweighted[byear_i] <- log(N_est) - fit_samples$resid[byear_i]
    
    # Now recalculate the entire time series
    for (y in unique(fit_samples$Year)){
      year_i <- which(fit_samples$Colony_Name == colony & fit_samples$Year == y) # indices
      
      # eta[i] = eta[byear] + r  ;   r = eta[i] - eta[byear]
      r_estimate <- fit_samples$eta[year_i] - fit_samples$eta[byear_i]
      fit_samples$eta_reweighted[year_i] <- fit_samples$eta_reweighted[byear_i] + r_estimate
    }
    
  }
  
  fit_samples$N_pred <- exp(fit_samples$eta_reweighted)
  
  # ------------------------------------------------
  # Colony-level summary
  # ------------------------------------------------
  
  N_summary_colony = fit_samples %>% 
    mutate(Species = spp) %>%
    group_by(Species,Colony_Name, Year) %>%
    summarize(q025 = quantile(N_pred,0.025),
              q50 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q975 = quantile(N_pred,0.975))
  
  N_summary_colony$regional <- "No"
  N_summary_colony$regional[N_summary_colony$Colony_Name %in% colonies_to_include_for_regional$Colony_Name] <- "Yes"
  
  Colony_Indices <- rbind(Colony_Indices,N_summary_colony)
  N_summary_colony <- subset(N_summary_colony,regional == "Yes")
  
  # -------------------------------------------------------
  # Calculate colony-level trends
  # -------------------------------------------------------
  
  Colony_Trends <- data.frame()
  
  for (colony in unique(fit_samples$Colony_Name)){
    
    first_survey <- subset(colony_summary, Colony_Name == colony)$first_survey
    most_recent_survey <- subset(colony_summary, Colony_Name == colony)$final_survey
    
    for (baseline_year in c(first_survey,most_recent_survey-10)){
      
      if (baseline_year < min(fit_samples$Year)) next
      
      baseline = subset(fit_samples, Year == baseline_year & Colony_Name == colony)$N_pred
      most_recent = subset(fit_samples, Year == most_recent_survey & Colony_Name == colony)$N_pred
      
      tmp <- data.frame(N_pred_baseline = subset(fit_samples, Year == baseline_year & Colony_Name == colony)$N_pred,
                        N_pred_recent = subset(fit_samples, Year == most_recent_survey & Colony_Name == colony)$N_pred) %>%
        mutate(percent_change = 100*(N_pred_recent - N_pred_baseline)/N_pred_baseline,
               trend = 100 * ((N_pred_recent/N_pred_baseline)^(1/(most_recent_survey-baseline_year))-1))
      
      # Summarize percent change in table format
      colony_trend = data.frame(
        
        species = spp,
        colony = colony,
        years = paste0(baseline_year,"-",most_recent_survey),
        year_start = baseline_year,
        year_end = most_recent_survey,
        
        trend = quantile(tmp$trend,0.500,na.rm = TRUE),
        trend_lci = quantile(tmp$trend,0.025,na.rm = TRUE),
        trend_uci = quantile(tmp$trend,0.975,na.rm = TRUE),
        
        percent_change = quantile(tmp$percent_change,0.500,na.rm = TRUE),
        percent_change_lci = quantile(tmp$percent_change,0.025,na.rm = TRUE),
        percent_change_uci = quantile(tmp$percent_change,0.975,na.rm = TRUE),
        
        prob_decrease_0 = mean(tmp$percent_change < 0,na.rm = TRUE),
        prob_decrease_25 = mean(tmp$percent_change<= -25,na.rm = TRUE),
        prob_decrease_30 = mean(tmp$percent_change<= -30,na.rm = TRUE),
        prob_decrease_50 = mean(tmp$percent_change<= -50,na.rm = TRUE),
        
        prob_increase_0 = mean(tmp$percent_change >= 0,na.rm = TRUE),
        prob_increase_33 = mean(tmp$percent_change>= (100/3),na.rm = TRUE),
        prob_increase_100 = mean(tmp$percent_change>= 100,na.rm = TRUE),
        
        prob_LD = mean(tmp$percent_change<= -50,na.rm = TRUE),
        prob_MD = mean(tmp$percent_change>= -50 & tmp$percent_change<= -25,na.rm = TRUE),
        prob_LC = mean(tmp$percent_change>= -25 & tmp$percent_change<= (100/3),na.rm = TRUE),
        prob_MI = mean(tmp$percent_change>= (100/3) & tmp$percent_change<= 100,na.rm = TRUE),
        prob_LI = mean(tmp$percent_change>= 100,na.rm = TRUE))
      
      Colony_Trends <- rbind(Colony_Trends,colony_trend)
      
      
    }
  }
  
  write.csv(Colony_Trends, file = paste0("output/model_results/tables/Colony_Trends_",spp,".csv"), row.names = FALSE)
  
  # ------------------------------------------------
  # Regional summary
  # ------------------------------------------------
  
  fit_samples = subset(fit_samples, Colony_Name %in% colonies_to_include_for_regional$Colony_Name)
  
  N_samples_regional = fit_samples %>% 
    group_by(Species,samp,Year) %>%
    summarize(N_pred = sum(N_pred))
  
  N_summary_regional <- N_samples_regional %>%
    group_by(Year) %>%
    summarize(q025 = quantile(N_pred,0.025),
              q50 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q975 = quantile(N_pred,0.975)) %>%
    mutate(Colony_Name = "Pacific Region", regional = "Yes") %>%
    bind_rows(N_summary_colony) %>%
    subset(regional == "Yes")
  
  # Order colonies by relative abundance in most recent year
  relabund <- N_summary_regional %>%
    subset(Year == max(N_summary_regional$Year)) %>%
    arrange(q50)
  N_summary_regional$Colony_Name <- factor(N_summary_regional$Colony_Name,
                                           levels = relabund$Colony_Name)
  
  
  colors = rev(viridis(nrow(relabund)))
  colors[length(colors)] <- "black"
  regional_plot <- ggplot(data = N_summary_regional, 
                          aes(x = Year, 
                              ymin=q025, 
                              y = q50,
                              ymax = q975,
                              label = Colony_Name,
                              col = Colony_Name,
                              fill = Colony_Name))+
    
    # Entire trajectory
    geom_ribbon(alpha = 0.4, linewidth = 0.5, col = "transparent")+
    geom_line(linewidth = 1)+
    
    scale_y_continuous(labels = comma)+
    scale_color_manual(values=colors, name = "")+
    scale_fill_manual(values=colors, name = "")+
    
    ylab("Index of Abundance")+
    ggtitle(seabird_codes_names$english_name[seabird_codes_names$Species == spp])
  regional_plot
  
  png(paste0("output/model_results/figures/",spp,"_colony_Plots_For_Laurie.png"), width = 6, height = 4, units = "in", res = 600)
  print(regional_plot)
  dev.off()
  
}

# Save annual colony-level indices
write.csv(Colony_Indices, file = "output/model_results/tables/SOCB_2023_Pacific_Indices_colony.csv", row.names = FALSE)

