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
allspecies_regional_change_samples = data.frame()

# Remove COEI from species list
atl_species = atl_species[-which(atl_species %in% c("COEI","CATE","RBGU"))]

# Extract predictions for each species
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
  
  # ------------------------------------------------
  # Regional summary
  # ------------------------------------------------
  
  # For calculating regional summaries, only include colonies with at least 1 survey prior to 2001 and 1 survey after 2010
  colonies_to_include_for_regional = colonies_to_include %>%
    subset(first_survey <= 1990 & last_survey >= 2000)
  
  N_samples = subset(N_samples, Colony_Name %in% colonies_to_include_for_regional$Colony_Name)
  
  N_samples_regional = N_samples %>% 
    group_by(Species,samp,Year) %>%
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
  
  # ----------------------------------------------
  # Calculate percent change relative to baseline year for each species
  # ----------------------------------------------
  
  # This is the year from which comparisons are made
  for (baseline_year in c(1970)){
    
    baseline = subset(N_samples_regional, Year == baseline_year) %>%
      mutate(N_pred_baseline = N_pred)
    
    tmp = full_join(N_samples_regional, baseline[,c("samp","N_pred_baseline")]) %>%
      mutate(percent_change = 100*(N_pred - N_pred_baseline)/N_pred_baseline) %>%
      dplyr::select(-N_pred_baseline) %>%
      add_column(baseline_year = baseline_year)
    
    # Summarize percent change in table format
    percent_change_summary_baseline = tmp %>% 
      group_by(Species,Year,baseline_year) %>%
      summarize(prob_decline = mean(percent_change < 0),
                
                # Change on percent scale
                percent_change_q025 = quantile(percent_change,0.025),
                percent_change_q05 = quantile(percent_change,0.05),
                percent_change_q10 = quantile(percent_change,0.1),
                percent_change_q50 = quantile(percent_change,0.5),
                percent_change_mean = mean(percent_change),
                percent_change_q90 = quantile(percent_change,0.9),
                percent_change_q95 = quantile(percent_change,0.95),
                percent_change_q975 = quantile(percent_change,0.975)) %>%
      subset(Year >= baseline_year) 

    # Append results to "full" dataframe for all species
    allspecies_percent_change_summary = rbind(allspecies_percent_change_summary,percent_change_summary_baseline)
    allspecies_N_summary_regional = rbind(allspecies_N_summary_regional,N_summary_regional)
    allspecies_regional_change_samples = rbind(allspecies_regional_change_samples,tmp)
    
    # ----------------------------------------------
    # Plot percent change relative to baseline for the species
    # ----------------------------------------------

    ylim = max(c(-90,900,abs(percent_change_summary_baseline$percent_change_q025),abs(percent_change_summary_baseline$percent_change_q975)))
    ylim = log((ylim)/100 + 1)
    y_labels = data.frame(pchange = c(-90,-50,0,100,900),labels = c("-90%","-50%","","+100%","+900%"))
    y_labels$r = log((y_labels$pchange)/100 + 1)


    fillcol = "gray70"
    
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
byear = 1970

# Y axis labels for figures (convert to percent change, but present on logarithmic scale so that a 9x increase looks the same as a 9x decrease)
ylim = max(c(-80,400,abs(subset(allspecies_percent_change_summary, baseline_year == byear)$percent_change_q50)))
ylim = log((ylim)/100 + 1)
y_labels = data.frame(pchange = c(-80,-50,0,100,400),labels = c("-80%","-50%","No change","+100%","+400%"))
y_labels$r = log((y_labels$pchange)/100 + 1)

# --------------------------------------
# Stacked bar chart illustrating species change categories
# --------------------------------------

species_change_status = allspecies_percent_change_summary %>%
  subset(Year == 2021 & baseline_year == byear)

species_change_status$Status = NA
species_change_status$Status[which(species_change_status$percent_change_mean <= -50)] = "Large Decrease"
species_change_status$Status[which(species_change_status$percent_change_mean >= -50 & species_change_status$percent_change_mean <= -25)] = "Moderate Decrease"
species_change_status$Status[which(species_change_status$percent_change_mean >= -25 & species_change_status$percent_change_mean <= 0)] = "Small Decrease"
species_change_status$Status[which(species_change_status$percent_change_mean >= 0 & species_change_status$percent_change_mean <= 33)] = "Small Increase"
species_change_status$Status[which(species_change_status$percent_change_mean >= 33 & species_change_status$percent_change_mean <= 100)] = "Moderate Increase"
species_change_status$Status[which(species_change_status$percent_change_mean >= 100)] = "Large Increase"

species_change_status$Status = factor(species_change_status$Status, 
                                      levels = c("Large Decrease",
                                                 "Moderate Decrease",
                                                 "Small Decrease",       
                                                 "Large Increase",
                                                 "Moderate Increase",
                                                 "Small Increase"))
# Counts of each status category
status_count = species_change_status %>%
  group_by(Status, .drop = FALSE) %>%
  count(Status, .drop = FALSE)

status_count$n[which(status_count$Status %in% c("Large Decrease","Moderate Decrease","Small Decrease"))] = -status_count$n[which(status_count$Status %in% c("Large Decrease","Moderate Decrease","Small Decrease"))]

# Status colors
status_colors = data.frame(Status = c("Large Decrease","Moderate Decrease","Small Decrease",
                                      "Large Increase","Moderate Increase","Small Increase"),
                           color = c("#a10012","#ff873f","#fce57a","#3a7af4","#5eadef","#aff9f2"))

plot_status = ggplot(status_count, aes(x = 1, y = n, fill = Status))+
  geom_bar(position = "stack", stat = "identity", width = 0.5)+
  geom_segment(aes(x = 0,xend = 2.5,y = 0, yend = 0))+
  
  geom_text(aes(x = 1.3, y = c(-4.5,-2.5,-0.5,1,2.5,6.5),
                label = c("Large Decrease","Moderate Decrease","Small Decrease",
                          "Small Increase","Moderate Increase","Large Increase")),
            hjust = 0, col = "gray60", size = 3)+
  
  
  scale_y_continuous(limits = c(-11,11), breaks = seq(-12,12,2), labels = c(seq(12,2,-2),0,seq(2,12,2)))+
  coord_cartesian(xlim = c(0.75,5))+
  scale_fill_manual(values = status_colors$color, guide = "none")+
  ylab("Number of species")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_blank()
  )+
  ggtitle(paste0("Status relative to ",byear))

print(plot_status)

pdf(paste0("output/Atlantic/model_results/figures/0_status_",byear,".pdf"), width = 5, height = 5)
print(plot_status)
dev.off()

# --------------------------------------
# All species on single figure, without uncertainty shown
# --------------------------------------

# Calculate the average (mean) percent change across species, to overlay as a black line
mean_change_df = allspecies_regional_change_samples

# Convert to logarithmic scale
mean_change_df$log_pchange = log(mean_change_df$percent_change/100 + 1)

# Average
mean_change_df = mean_change_df %>%
  subset(Year >= byear & baseline_year == byear) %>%
  group_by(samp,Year) %>%
  summarize(N_pred = sum(N_pred),
            log_pchange = mean(log_pchange))

# Convert back to non-log (percent change) scale
mean_change_df$percent_change = 100 * (exp(mean_change_df$log_pchange)-1)

mean_change_summary = mean_change_df %>%
  group_by(Year) %>%
  summarize(
    # Total abundance of all seabirds
    N_pred_q025 = quantile(N_pred,0.025),
    N_pred_q05 = quantile(N_pred,0.05),
    N_pred_q10 = quantile(N_pred,0.1),
    N_pred_q50 = quantile(N_pred,0.5),
    N_pred_mean = mean(N_pred),
    N_pred_q90 = quantile(N_pred,0.9),
    N_pred_q95 = quantile(N_pred,0.05),
    N_pred_q975 = quantile(N_pred,0.975),
    
    # Change on percent scale
    percent_change_q025 = quantile(percent_change,0.025),
    percent_change_q05 = quantile(percent_change,0.05),
    percent_change_q10 = quantile(percent_change,0.1),
    percent_change_q50 = quantile(percent_change,0.5),
    percent_change_mean = mean(percent_change),
    percent_change_q90 = quantile(percent_change,0.9),
    percent_change_q95 = quantile(percent_change,0.05),
    percent_change_q975 = quantile(percent_change,0.975)) %>%
  add_column(Species = "All species")

# Status categories (for labeling and colouring)
allspecies_percent_change_summary = full_join(allspecies_percent_change_summary,species_change_status[,c("Species","Status")])

plot_allspecies = ggplot(subset(allspecies_percent_change_summary, baseline_year == byear)) +
  geom_hline(yintercept = 0, col = "gray80", size = 1, linetype = 1)+
  geom_line(aes(x = Year, y = log((percent_change_q50)/100 + 1), linetype = Species, col = Status), size = 0.5)+
  
  geom_ribbon(data = mean_change_summary, aes(x = Year, ymin = log((percent_change_q025)/100 + 1),ymax = log((percent_change_q975)/100 + 1)), alpha = 0.1)+
  
  geom_line(data = mean_change_summary, aes(x = Year, y = log((percent_change_q50)/100 + 1)), size = 1, col = "black")+
  geom_text_repel(data = subset(allspecies_percent_change_summary, baseline_year == byear & Year == 2021), aes(x = Year, y = log((percent_change_q50)/100 + 1), col = Status, label = Species), hjust = 0, fontface = "bold", direction = "y", nudge_x = 2, force = 2, size = 2, segment.alpha = 0.3, segment.size = 0.2)+
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

# --------------------------------------
# Total number of birds
# --------------------------------------

plot_total_birds = ggplot(data = mean_change_summary, 
                          aes(x = Year, ymin = N_pred_q025, ymax = N_pred_q975, y = N_pred_q50)) +
  
  geom_line()+
  geom_ribbon(alpha = 0.1)+
  
  ylab("Total seabirds")+
  xlab("Year")+
  ggtitle("Total seabirds in monitored colonies")+
  coord_cartesian(xlim = c(byear,2025))+
  scale_y_continuous(labels = comma)

pdf(paste0("output/Atlantic/model_results/figures/0_total_monitored_birds_",byear,".pdf"), width = 5, height = 3)
print(plot_total_birds)
dev.off()
