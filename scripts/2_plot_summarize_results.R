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

atl = read_xlsx("data/ATLANTIC_QC_AR_DATA_2022.xlsx", sheet = 1)

# Convert pairs to individuals
atl$Count[atl$Count_Type == "Pairs"] = atl$Count[atl$Count_Type == "Pairs"]*2
atl$Count_Type[atl$Count_Type == "Pairs"] = "Individuals"
atl = subset(atl, Count_Type == "Individuals") # Removes several AOS rows

atl_species = unique(atl$Species)

# ------------------------------------------------
# Loop through species and extract/plot predictions for each, separately
# ------------------------------------------------

# Empty dataframe to store trend/trajectory summaries for all species
allspecies_percent_change_summary = data.frame()
allspecies_N_summary_regional = data.frame()

for (spp in atl_species){
  
  # Load model results
  file = paste0("~/iles_ECCC/Seabirds/Status_of_Birds_seabird_analysis_2022/output/Atlantic/model_results/",spp,".RData")
  if(!file.exists(file)) next 
  load(file)
  
  # Extract predictions in dataframe format
  N_samples = reshape2::melt(out$sims.list$etapred) %>%
    rename(samp = Var1, year_number = Var2, colony_number = Var3, N_pred = value) %>%
    full_join(colony_name_table, by = c("colony_number" = "colony_numeric")) %>%
    full_join(year_table, by = c("year_number" = "year_numeric")) %>%
    add_column(Species = spp)
  
  # ------------------------------------------------
  # Colony-level summaries
  # ------------------------------------------------
  
  N_summary_colony = N_samples %>% group_by(Colony_Name, Year) %>%
    summarize(q05 = quantile(N_pred,0.05),
              q50 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q95 = quantile(N_pred,0.95))
  
  colors = viridis(10)[c(3,6,9)]
  
  # Dynamics from 1970 onwards
  colony_plot_1970 = ggplot() +
    geom_ribbon(data = subset(N_summary_colony, Year >= 1970), aes(x = Year, ymin = q05, ymax = q95), fill = "gray85", col = "transparent")+
    geom_line(data = subset(N_summary_colony, Year >= 1970), aes(x = Year, y = q50), col = "black")+
    geom_vline(data = subset(spdat, Count == 0 & Year >= 1970),aes(xintercept = Year), col = "red", size = 1.5, alpha = 0.3)+
    geom_point(data = subset(atl, Species == spp & Year >= 1970),aes(x = Year, y = Count, col = Count_Type))+
    xlab("Year")+
    ylab("Population index")+
    scale_color_manual(values = colors,name = "Count type", drop = FALSE)+
    facet_wrap(Colony_Name~., scales = "free_y")+
    scale_y_continuous(labels = comma)+
    ggtitle(paste0(spp," colony-level trajectories"))
  
  pdf(paste0("output/Atlantic/model_results/figures/",spp,"_colony_1970.pdf"), width = 20,height = 10)
  print(colony_plot_1970)
  dev.off()
  
  # Dynamics from 1990 onwards
  colony_plot_1990 = ggplot() +
    geom_ribbon(data = subset(N_summary_colony, Year >= 1990), aes(x = Year, ymin = q05, ymax = q95), fill = "gray85", col = "transparent")+
    geom_line(data = subset(N_summary_colony, Year >= 1990), aes(x = Year, y = q50), col = "black")+
    geom_vline(data = subset(spdat, Count == 0 & Year >= 1990),aes(xintercept = Year), col = "red", size = 1.5, alpha = 0.3)+
    geom_point(data = subset(atl, Species == spp & Year >= 1990),aes(x = Year, y = Count, col = Count_Type))+
    xlab("Year")+
    ylab("Population index")+
    scale_color_manual(values = colors,name = "Count type", drop = FALSE)+
    facet_wrap(Colony_Name~., scales = "free_y")+
    scale_y_continuous(labels = comma)+
    ggtitle(paste0(spp," colony-level trajectories"))
  
  pdf(paste0("output/Atlantic/model_results/figures/",spp,"_colony_1990.pdf"), width = 20,height = 10)
  print(colony_plot_1990)
  dev.off()
  
  # ------------------------------------------------
  # Regional summary
  # ------------------------------------------------
  
  # For calculating regional summaries, only include colonies with at least 1 survey prior to 2001 and 1 survey after 2010
  colonies_to_include_for_regional = colonies_to_include %>%
    subset(first_survey < 2001 & last_survey > 2010)
  
  N_samples = subset(N_samples, Colony_Name %in% colonies_to_include_for_regional$Colony_Name)
  
  N_samples_regional = N_samples %>% group_by(Species,samp,Year) %>%
    summarize(N_pred = sum(N_pred))
  
  N_summary_regional = N_samples_regional %>%
    group_by(Species,Year) %>%
    summarize(q05 = quantile(N_pred,0.05),
              q50 = quantile(N_pred,0.500),
              mean = mean(N_pred),
              q95 = quantile(N_pred,0.95))
  
  # Dynamics from 1970 onwards
  regional_plot_1970 = ggplot() +
    geom_ribbon(data = subset(N_summary_regional, Year >= 1970), aes(x = Year, ymin = q05, ymax = q95), fill = "gray85", col = "transparent")+
    geom_line(data = subset(N_summary_regional, Year >= 1970), aes(x = Year, y = q50), col = "black")+
    xlab("Year")+
    ylab("Population index")+
    scale_y_continuous(labels = comma)+
    ggtitle(paste0(spp," regional trajectory since 1970"))
  
  pdf(paste0("output/Atlantic/model_results/figures/",spp,"_regional_1970.pdf"))
  print(regional_plot_1970)
  dev.off()
  
  # Dynamics from 1990 onwards
  regional_plot_1990 = ggplot() +
    geom_ribbon(data = subset(N_summary_regional, Year >= 1990), aes(x = Year, ymin = q05, ymax = q95), fill = "gray85", col = "transparent")+
    geom_line(data = subset(N_summary_regional, Year >= 1990), aes(x = Year, y = q50), col = "black")+
    xlab("Year")+
    scale_y_continuous(labels = comma)+
    ylab("Population index")+
    ggtitle(paste0(spp," regional trajectory since 1990"))
  
  pdf(paste0("output/Atlantic/model_results/figures/",spp,"_regional_1990.pdf"))
  print(regional_plot_1990)
  dev.off()
  
  # ----------------------------------------------
  # Calculate percent change relative to baseline year for each species
  # ----------------------------------------------
  
  # This is the year from which comparisons are made
  for (baseline_year in c(1970,1991)){
    
    baseline = subset(N_samples_regional, Year == baseline_year) %>%
      mutate(N_pred_baseline = N_pred)
    
    tmp = full_join(N_samples_regional, baseline[,c("samp","N_pred_baseline")]) %>%
      mutate(percent_change = 100*(N_pred - N_pred_baseline)/N_pred_baseline) %>%
      dplyr::select(-N_pred_baseline)
    
    # Summarize percent change in table format
    percent_change_summary_baseline = tmp %>%
      add_column(baseline_year = baseline_year) %>%
      group_by(Species,Year,baseline_year) %>%
      summarize(prob_decline = mean(percent_change < 0),
                
                # Change on percent scale
                percent_change_q025 = quantile(percent_change,0.025),
                percent_change_q05 = quantile(percent_change,0.05),
                percent_change_q10 = quantile(percent_change,0.1),
                percent_change_q50 = quantile(percent_change,0.5),
                percent_change_mean = mean(percent_change),
                percent_change_q90 = quantile(percent_change,0.9),
                percent_change_q95 = quantile(percent_change,0.05),
                percent_change_q975 = quantile(percent_change,0.975)) %>%
      subset(Year >= baseline_year) 
    
    # Assign a "current Status" to the species (for plotting)
    current_status = subset(percent_change_summary_baseline, Year == 2021)
    Status = "Uncertain"
    # If there is a > 90% chance the species has declined, colour it red
    if (current_status$percent_change_q90 < 0) Status = "Decrease"
    # If there is a > 90% chance the species has increased, colour it blue
    if (current_status$percent_change_q10 > 0) Status = "Increase"
    percent_change_summary_baseline$Status = Status
    
    
    # Append results to "full" dataframe for all species
    allspecies_percent_change_summary = rbind(allspecies_percent_change_summary,percent_change_summary_baseline)
    
    allspecies_N_summary_regional = rbind(allspecies_N_summary_regional,N_summary_regional)
    
    # ----------------------------------------------
    # Plot percent change relative to baseline for the species
    # ----------------------------------------------
    
    ylim = max(c(-90,900,abs(percent_change_summary_baseline$percent_change_q025),abs(percent_change_summary_baseline$percent_change_q975)))
    ylim = log((ylim)/100 + 1)
    y_labels = data.frame(pchange = c(-90,-50,0,100,900),labels = c("-90%","-50%","","+100%","+900%"))
    y_labels$r = log((y_labels$pchange)/100 + 1)
    
    
    fillcol = "gray70"
    # If there is a > 90% chance the species has declined, colour it red
    if (current_status$percent_change_q90 < 0) fillcol = "#f94b55"
    # If there is a > 90% chance the species has increased, colour it blue
    if (current_status$percent_change_q10 > 0) fillcol = "cornflowerblue"
    
    spp_pchange_plot = ggplot(percent_change_summary_baseline) +
      geom_hline(yintercept = 0, col = "gray95", size = 1.5)+
      geom_ribbon(aes(x = Year, ymin = log((percent_change_q10)/100 + 1), ymax = log((percent_change_q90)/100 + 1), linetype = Species),fill=fillcol, alpha = 0.2)+
      geom_line(aes(x = Year, y = log((percent_change_q50)/100 + 1), linetype = Species), size = 1.5, col = fillcol)+
      scale_linetype_manual(values = rep(1,length(unique(percent_change_summary_baseline$Species))), guide = "none")+
      ylab(paste0("Percent change relative to ",baseline_year))+
      xlab("Year")+
      ggtitle(paste0(spp," percent change since ",baseline_year))+
      scale_y_continuous(breaks = y_labels$r, labels = y_labels$labels)+
      coord_cartesian(ylim=c(-ylim,ylim),xlim = c(baseline_year,max(percent_change_summary_baseline$Year)))
    print(spp_pchange_plot)
    
    pdf(paste0("output/Atlantic/model_results/figures/",spp,"_percentchange_since_",baseline_year,".pdf"))
    print(spp_pchange_plot)
    dev.off()
    
    
  }
  
}

# ------------------------------------------------
# Output trend estimates
# ------------------------------------------------
write.csv(allspecies_percent_change_summary, file = "output/Atlantic/model_results/tables/Atlantic_change_estimates_2022.csv", row.names = FALSE)

# ------------------------------------------------
# Output annual regional indices
# ------------------------------------------------
write.csv(allspecies_N_summary_regional, file = "output/Atlantic/model_results/tables/Atlantic_annual_indices_2022.csv", row.names = FALSE)

# ------------------------------------------------
# Generate a figure of the trajectory of individual species groups
# (in terms of percent change since 1970)
# ------------------------------------------------

# baseline year to plot
for (byear in c(1970,1991)){
  
  # Y axis labels for figures (convert to percent change, but present on logarithmic scale so that a 9x increase looks the same as a 9x decrease)
  ylim = max(c(-90,900,abs(subset(allspecies_percent_change_summary, baseline_year == byear)$percent_change_q50)))
  ylim = log((ylim)/100 + 1)
  y_labels = data.frame(pchange = c(-90,-50,0,100,900),labels = c("-90%","-50%","","+100%","+900%"))
  y_labels$r = log((y_labels$pchange)/100 + 1)
  
  # --------------------------------------
  # Plot trajectory of each species with uncertainty included
  # --------------------------------------
  allspecies_percent_change_summary$Status = factor(allspecies_percent_change_summary$Status, levels = c("Increase","Uncertain","Decrease"))
  
  plot_species_facets = ggplot(subset(allspecies_percent_change_summary, baseline_year == byear)) +
    geom_hline(yintercept = 0, col = "gray95", size = 1.5)+
    geom_ribbon(aes(x = Year, ymin = log((percent_change_q10)/100 + 1), ymax = log((percent_change_q90)/100 + 1), linetype = Species, fill = Status), alpha = 0.2)+
    geom_line(aes(x = Year, y = log((percent_change_q50)/100 + 1), linetype = Species, col = Status), size = 1.5)+
    scale_linetype_manual(values = rep(1,length(unique(allspecies_percent_change_summary$Species))), guide = "none")+
    scale_fill_manual(values = c("cornflowerblue","gray70","#f94b55"), name = "Current abundance\nrelative to baseline year", guide = "none")+
    scale_color_manual(values = c("cornflowerblue","gray70","#f94b55"), name = "Current abundance\nrelative to baseline year", guide = "none")+
    ylab(paste0("Percent change relative to ",byear))+
    xlab("Year")+
    facet_wrap(Species~.)+
    ggtitle(paste0("Percent change relative to ",byear))+
    scale_y_continuous(breaks = y_labels$r, labels = y_labels$labels)+
    coord_cartesian(ylim=c(-ylim,ylim),xlim = c(byear,max(allspecies_percent_change_summary$Year)))
  print(plot_species_facets)
  
  pdf(paste0("output/Atlantic/model_results/figures/0_percent_change_speciesfacets_",byear,".pdf"), width = 7, height = 7)
  print(plot_species_facets)
  dev.off()
  
  # --------------------------------------
  # All species on single figure, without uncertainty shown
  # --------------------------------------
  plot_allspecies = ggplot(subset(allspecies_percent_change_summary, baseline_year == byear)) +
    geom_hline(yintercept = 0, col = "gray95", size = 1.5)+
    geom_line(aes(x = Year, y = log((percent_change_q50)/100 + 1), linetype = Species, col = Status), size = 1.5)+
    geom_label_repel(data = subset(allspecies_percent_change_summary, baseline_year == byear & Year == 2021), aes(x = Year, y = log((percent_change_q50)/100 + 1), col = Status, label = Species), hjust = 0, fontface = "bold", direction = "y", nudge_x = 2, force = 2, size = 3, segment.alpha = 0.3)+
    scale_linetype_manual(values = rep(1,length(unique(allspecies_percent_change_summary$Species))), guide = "none")+
    scale_fill_manual(values = c("cornflowerblue","gray70","#f94b55"), name = "Current abundance\nrelative to baseline year", guide = "none")+
    scale_color_manual(values = c("cornflowerblue","gray70","#f94b55"), name = "Current abundance\nrelative to baseline year", guide = "none")+
    ylab(paste0("Percent change relative to ",byear))+
    xlab("Year")+
    ggtitle(paste0("Percent change relative to ",byear))+
    scale_y_continuous(breaks = y_labels$r, labels = y_labels$labels)+
    coord_cartesian(ylim=c(-ylim,ylim),xlim = c(byear,2025))
  print(plot_allspecies)
  
  pdf(paste0("output/Atlantic/model_results/figures/0_percent_change_allspecies_",byear,".pdf"), width = 7, height = 5)
  print(plot_allspecies)
  dev.off()
}

