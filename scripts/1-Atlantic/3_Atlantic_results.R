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
setwd("../../")

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

# ************************************************
# Read in data
# ************************************************

# ------------------------------------------------
# Read in data
# ------------------------------------------------

atl = read_xlsx("data/ATLANTIC_QC_AR_DATA_2023_updated.xlsx", sheet = 1)

# Convert pairs to individuals
atl$Count[atl$Count_Type == "Pairs"] = atl$Count[atl$Count_Type == "Pairs"]*2
atl$Count_Type[atl$Count_Type == "Pairs"] = "Individuals"
atl = subset(atl, Count_Type == "Individuals") # Removes several AOS rows

atl_species = unique(atl$Species)
#atl_species = atl_species[-which(atl_species %in% c("COEI","CATE","RBGU"))]

# ------------------------------------------------
# Prepare species codes / labels / generation lengths
# ------------------------------------------------

GenLength <- read_xlsx("data/Bird_et_al_2020_AppendixS4_GenLength.xlsx") %>%
  select('Scientific name','GenLength')
SOBC_species <- read_xlsx("from_collaborators/SOBC_Template_From_Birds_Canada/SOCB_species.xlsx")

seabird_codes_names <- read_xlsx("data/seabird_names_2023.xlsx") %>%
  left_join(SOBC_species) %>%
  left_join(GenLength, by = c('scientific_name' = 'Scientific name')) %>%
  mutate(GenLength = round(GenLength))

seabird_codes_names %>% as.data.frame()

# ************************************************
# Loop through species and extract/process results
# ************************************************

# Empty dataframes to store results
Trend_Estimates <- data.frame()
Annual_Indices <- data.frame()

for (spp in atl_species){
  
  # Load model results
  file = paste0("output/Atlantic/model_results/",spp,".RData")
  if(!file.exists(file)) next 
  load(file)
  
  # Extract predictions in dataframe format
  fit_samples = reshape2::melt(out$sims.list$population_index) %>%
    rename(samp = Var1, year_number = Var2, colony_number = Var3, N_pred = value) %>%
    full_join(colony_name_table, by = c("colony_number" = "colony_numeric")) %>%
    full_join(year_table, by = c("year_number" = "year_numeric")) %>%
    add_column(Species = spp)
  
  # ------------------------------------------------
  # Identify reasonable start/end points for 'reliable' trend estimates
  # ------------------------------------------------
  
  colony_summary <- spdat %>%
    group_by(Colony_Name) %>%
    summarize(mean_count = mean(Count,na.rm = TRUE),
              first_survey = min(Year),
              final_survey = max(Year)) %>%
    arrange(desc(mean_count)) %>%
    ungroup() %>%
    mutate(include = mean_count >= (sum(mean_count,na.rm = TRUE)*0.2))
  
  if (sum(colony_summary$include,na.rm = TRUE)==0) colony_summary$include <- TRUE
  
  empirical_start <- ceiling(mean(subset(colony_summary, include)$first_survey))
  empirical_end <- floor(mean(subset(colony_summary, include)$final_survey))
  
  empirical_start <- max(c(empirical_start, min(spdat$Year),1970))
  empirical_end <- min(c(empirical_end, max(spdat$Year)))
  
  if (spp == "ROST") empirical_end = 2022
  
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
  # Colony-level summary
  # ------------------------------------------------
  
  N_summary_colony = fit_samples %>% 
    group_by(Colony_Name, Year) %>%
    summarize(q050 = quantile(N_pred,0.025),
              q50 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q950 = quantile(N_pred,0.975))
  
  N_summary_colony$regional <- "No"
  N_summary_colony$regional[N_summary_colony$Colony_Name %in% colonies_to_include_for_regional$Colony_Name] <- "Yes"
  
  # ------------------------------------------------
  # Regional summary
  # ------------------------------------------------
  
  first_year_for_summary = empirical_start
  final_year_for_summary = empirical_end
  
  fit_samples = subset(fit_samples, Colony_Name %in% colonies_to_include_for_regional$Colony_Name)
  
  N_samples_regional = fit_samples %>% 
    group_by(Species,samp,Year) %>%
    summarize(N_pred = sum(N_pred))
  
  # ----------------------------------------------
  # Calculate "LOESS" smooth for Birds Canada based on posterior medians
  # ----------------------------------------------
  
  pop_index_regional <- N_samples_regional %>%
    subset(Year >= first_year_for_summary & Year <= final_year_for_summary) %>%
    group_by(Species,Year) %>%
    summarize(index = quantile(N_pred,0.500))
  
  LOESS = loess(log(index)~Year, data=pop_index_regional, span=0.55, na.action = na.exclude)
  
  #add the predicted LOESS values back to the original dataframe with a new column name
  pop_index_regional$LOESS_index <-exp(predict(LOESS) )
  
  #ggplot()+
  #       geom_point(data = pop_index_regional, aes(x = Year, y = index))+
  #       geom_line(data = pop_index_regional, aes(x = Year, y = LOESS_index))
  
  # ----------------------------------------------
  # Calculate percent change relative to a baseline year for each species
  # ----------------------------------------------
  
  spp_GenLength = subset(seabird_codes_names, Species == spp)$GenLength
  short_term_trend_start <- max(first_year_for_summary,final_year_for_summary - spp_GenLength*3)
  short_term_trend_start <- final_year_for_summary - spp_GenLength*3
  
  # Baseline years for trend summary
  Baseline_Years <- c(max(c(1970,first_year_for_summary)),short_term_trend_start,final_year_for_summary-10)
  
  # This is the year from which comparisons are made
  for (i in 1:length(Baseline_Years)){
    
    # ------------------------------------------
    # Trend "period" labels
    # ------------------------------------------
    
    baseline_year = Baseline_Years[i]
    if (i == 1) trend_type = "all years"
    if (i == 2) trend_type = "3Gen-Recent"
    if (i == 3) trend_type = "10-years"
    
    baseline = subset(N_samples_regional, Year == baseline_year) %>%
      mutate(N_pred_baseline = N_pred)
    
    tmp = full_join(N_samples_regional, baseline[,c("samp","N_pred_baseline")]) %>%
      mutate(percent_change = 100*(N_pred - N_pred_baseline)/N_pred_baseline,
             trend = 100 * ((N_pred/N_pred_baseline)^(1/(Year-baseline_year))-1)) %>%
      dplyr::select(-N_pred_baseline) %>%
      add_column(baseline_year = baseline_year) %>%
      subset(Year == final_year_for_summary)
    
    if (baseline_year < first_year_for_summary) tmp <- tmp %>% mutate(percent_change = NA, trend = NA)
    
    # Summarize percent change in table format
    species_trend = data.frame(
      
      # Use fields as required by NatureCounts
      results_code = "SCMP",
      version = 2023,
      area_code = "Atlantic",
      species_code = seabird_codes_names$species_code[seabird_codes_names$Species == spp],
      species_id = seabird_codes_names$species_id[seabird_codes_names$Species == spp],
      season = "breeding",
      period = trend_type,
      years = paste0(baseline_year,"-",mean(tmp$Year,na.rm = TRUE)),
      year_start = baseline_year,
      year_end = mean(tmp$Year),
      trnd = quantile(tmp$trend,0.500,na.rm = TRUE),
      index_type = NA,
      upper_ci = quantile(tmp$trend,0.975,na.rm = TRUE),
      lower_ci = quantile(tmp$trend,0.025,na.rm = TRUE),
      stderr = NA,
      model_type = "GAMM",
      model_fit = NA,
      percent_change = quantile(tmp$percent_change,0.500,na.rm = TRUE),
      percent_change_low = quantile(tmp$percent_change,0.025,na.rm = TRUE),
      percent_change_high = quantile(tmp$percent_change,0.975,na.rm = TRUE),
      
      prob_decrease_0 = mean(tmp$percent_change < 0,na.rm = TRUE),
      prob_decrease_25 = mean(tmp$percent_change<= -25,na.rm = TRUE),
      prob_decrease_30 = mean(tmp$percent_change<= -30,na.rm = TRUE),
      prob_decrease_50 = mean(tmp$percent_change<= -50,na.rm = TRUE),
      
      prob_increase_0 = mean(tmp$percent_change >= 0,na.rm = TRUE),
      prob_increase_33 = mean(tmp$percent_change>= (100/3),na.rm = TRUE),
      prob_increase_100 = mean(tmp$percent_change>= 100,na.rm = TRUE),
      
      reliability = NA,
      precision_num = NA,
      precision_cat = NA,
      coverage_num = proportion,
      coverage_cat = coverage_cat,
      
      goal = NA,
      goal_lower = NA,
      sample_size = NA,
      sample_total = NA,
      subtitle = NA,
      
      prob_LD = mean(tmp$percent_change<= -50,na.rm = TRUE),
      prob_MD = mean(tmp$percent_change>= -50 & tmp$percent_change<= -25,na.rm = TRUE),
      prob_LC = mean(tmp$percent_change>= -25 & tmp$percent_change<= (100/3),na.rm = TRUE),
      prob_MI = mean(tmp$percent_change>= (100/3) & tmp$percent_change<= 100,na.rm = TRUE),
      prob_LI = mean(tmp$percent_change>= 100,na.rm = TRUE))
    
    Trend_Estimates <- rbind(Trend_Estimates,species_trend)
    
    # ------------------------------------------
    # Population indices
    # ------------------------------------------
    
    if (baseline_year < first_year_for_summary) next
    
    species_indices <- N_samples_regional %>%
      subset(Year >= baseline_year & Year <= final_year_for_summary) %>%
      group_by(Species,Year) %>%
      summarize(results_code = "SCMP",
                area_code = "Atlantic",
                year = mean(Year),
                season = "breeding",
                period = trend_type,
                species_code = seabird_codes_names$species_code[seabird_codes_names$Species == spp],
                species_id = seabird_codes_names$species_id[seabird_codes_names$Species == spp],
                index = quantile(N_pred,0.500),
                stderr = NA,
                stdev = sd(N_pred),
                upper_ci = quantile(N_pred,0.975),
                lower_ci = quantile(N_pred,0.025))
    
    species_indices <- left_join(species_indices, pop_index_regional[,c("Year","LOESS_index")])
    
    Annual_Indices <- rbind(Annual_Indices,species_indices)
    
    
  }
  
  # ------------------------------------------------
  # Plot colony-level trajectories
  # ------------------------------------------------
  
  N_summary_colony$regional <- factor(N_summary_colony$regional, levels = c("Yes","No"))
  
  # Colony-level dynamics
  colony_plot_freeaxis = ggplot() +
    
    # Full time series
    geom_ribbon(data = N_summary_colony, aes(x = Year, ymin = q050, ymax = q950), fill = "gray80",col = "transparent")+
    geom_line(data = N_summary_colony, aes(x = Year, y = q50), col = "gray60")+
    
    geom_ribbon(data = subset(N_summary_colony, Year >= first_year_for_summary & Year <= final_year_for_summary), aes(x = Year, ymin = q050, ymax = q950, fill = regional), col = "transparent")+
    geom_line(data = subset(N_summary_colony, Year >= first_year_for_summary & Year <= final_year_for_summary), aes(x = Year, y = q50, col = regional))+
    
    geom_point(data = spdat,aes(x = Year, y = Count))+
    xlab("Year")+
    ylab("Population index")+
    scale_color_manual(values = c("darkblue","gray60"),name = "Included in\nRegional\nSummary", drop = FALSE)+
    scale_fill_manual(values = c("dodgerblue","gray80"),name = "Included in\nRegional\nSummary", drop = FALSE)+
    
    facet_wrap(Colony_Name~., scales = "free_y")+
    scale_y_continuous(labels = comma, trans = "log10")+
    ggtitle(seabird_codes_names$english_name[seabird_codes_names$Species == spp])
  
  png(paste0("output/Atlantic/model_results/figures/",spp,"_colony_trajectories.png"), width = 12, height = 6, units = "in", res = 600)
  print(colony_plot_freeaxis)
  dev.off()
  
  # ------------------------------------------------
  # Plot regional trajectory
  # ------------------------------------------------
  
  # Plot 500 samples from the posterior (to visualize spread in trajectories)
  samples_to_plot <- sample(unique(N_samples_regional$samp),500)
  
  
  regional_summary <- N_samples_regional %>%
    group_by(Year) %>%
    summarize(N = median(N_pred),
              N_lcl = quantile(N_pred,0.025),
              N_ucl = quantile(N_pred,0.975))
  
  ylim = regional_summary %>%
    subset(Year >= first_year_for_summary & Year <= final_year_for_summary) %>%
    summarize(min = min(N_lcl)*0.5,
              max = max(N_ucl)*2)
  regional_plot <- ggplot()+
    
    # Entire trajectory
    geom_ribbon(data = regional_summary, aes(x = Year, ymin=N_lcl, ymax = N_ucl),alpha = 0.2, fill = "gray80", col = "gray60", linetype = 2, linewidth = 0.5)+
    geom_line(data = regional_summary, aes(x = Year, y=N), linewidth = 1, col = "gray60")+
    
    # Reliable survey window
    geom_ribbon(data = subset(regional_summary, Year >= first_year_for_summary & Year <= final_year_for_summary), aes(x = Year, ymin=N_lcl, ymax = N_ucl),alpha = 0.1, fill = "dodgerblue", col = "black", linetype = 2, linewidth = 0.5)+
    geom_line(data = subset(N_samples_regional, samp %in% samples_to_plot & Year >= first_year_for_summary & Year <= final_year_for_summary), aes(x = Year, y = N_pred, col = factor(samp)),alpha = 0.1)+
    geom_line(data = subset(regional_summary, Year >= first_year_for_summary & Year <= final_year_for_summary), aes(x = Year, y=N), linewidth = 1, col = "black")+
    
    coord_cartesian(ylim = c(ylim$min,ylim$max))+
    scale_y_continuous(labels = comma)+
    scale_color_manual(values=rep("blue",length(unique(N_samples_regional$samp))), guide = "none")+
    ylab("Index of Abundance")+
    ggtitle(seabird_codes_names$english_name[seabird_codes_names$Species == spp])
  
  png(paste0("output/Atlantic/model_results/figures/",spp,"_regional.png"), width = 6, height = 4, units = "in", res = 600)
  print(regional_plot)
  dev.off()
  
}

# Assign precision categories
Trend_Estimates$precision_cat <- "Medium"
Trend_Estimates$precision_cat[which((Trend_Estimates$upper_ci - Trend_Estimates$lower_ci)> 6.7)] <- "Low"
Trend_Estimates$precision_cat[which((Trend_Estimates$upper_ci - Trend_Estimates$lower_ci)< 3.5)] <- "High"

# ------------------------------------------------
# Output trend/change estimates
# ------------------------------------------------

write.csv(Trend_Estimates, file = "output/Atlantic/model_results/tables/SOCB_2023_Atlantic_Trends.csv", row.names = FALSE)

# ------------------------------------------------
# Output annual regional indices
# ------------------------------------------------

write.csv(Annual_Indices, file = "output/Atlantic/model_results/tables/SOCB_2023_Atlantic_Indices.csv", row.names = FALSE)

# ------------------------------------------------
# Summary plots
# ------------------------------------------------

Trend_Estimates <- left_join(Trend_Estimates, SOBC_species)

lim = max(abs(Trend_Estimates[,c("upper_ci","lower_ci")]),na.rm = TRUE)
lim = c(-lim,lim)

allspecies_trends_highquality_plot <- ggplot(data = subset(Trend_Estimates,period == "all years"), aes(y = english_name, x = trnd, xmin = lower_ci, xmax = upper_ci))+
  geom_vline(xintercept = 0, linetype = 2, col = "gray80", linewidth = 1)+
  geom_errorbarh(height=0)+
  geom_point()+
  geom_text(aes(y = english_name, x = trnd, label = years),vjust=-0.5, size = 3)+
  theme_bw()+
  ylab("Species")+
  xlab("Trend (% change per year)")+
  ggtitle("Full time series")+
  coord_cartesian(xlim=lim)

png(paste0("output/Atlantic/model_results/figures/0_allspecies_trends_highquality.png"), width = 6, height = 6, units = "in", res = 600)
print(allspecies_trends_highquality_plot)
dev.off()

allspecies_trends_3Gen_plot <- ggplot(data = subset(Trend_Estimates,period == "3Gen-Recent"), aes(y = english_name, x = trnd, xmin = lower_ci, xmax = upper_ci))+
  geom_vline(xintercept = 0, linetype = 2, col = "gray80", linewidth = 1)+
  geom_errorbarh(height=0)+
  geom_point()+
  geom_text(aes(y = english_name, x = trnd, label = years),vjust=-0.5, size = 3)+
  theme_bw()+
  ylab("Species")+
  xlab("Trend (% change per year)")+
  ggtitle("3Gen-Recent")+
  coord_cartesian(xlim=lim)

png(paste0("output/Atlantic/model_results/figures/0_allspecies_trends_3Gen.png"), width = 6, height = 6, units = "in", res = 600)
print(allspecies_trends_3Gen_plot)
dev.off()

allspecies_trends_10yr_plot <- ggplot(data = subset(Trend_Estimates,period == "10-years"), aes(y = english_name, x = trnd, xmin = lower_ci, xmax = upper_ci))+
  geom_vline(xintercept = 0, linetype = 2, col = "gray80", linewidth = 1)+
  geom_errorbarh(height=0)+
  geom_point()+
  geom_text(aes(y = english_name, x = trnd, label = years),vjust=-0.5, size = 3)+
  theme_bw()+
  ylab("Species")+
  xlab("Trend (% change per year)")+
  ggtitle("10-years")+
  coord_cartesian(xlim=lim)

png(paste0("output/Atlantic/model_results/figures/0_allspecies_trends_10yr.png"), width = 6, height = 6, units = "in", res = 600)
print(allspecies_trends_10yr_plot)
dev.off()