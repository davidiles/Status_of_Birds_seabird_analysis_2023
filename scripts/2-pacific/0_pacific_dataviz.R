# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------
my_packs = c('tidyverse','readxl','RColorBrewer','viridis','jagsUI','mgcv','ggrepel')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())
theme_set(theme_bw())
setwd("D:/iles_ECCC/Seabirds/Status_of_Birds_seabird_analysis_2022/")
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
# "Complete" counts
# ------------------------------------------------

pac_cc = read_xlsx("data/PACIFIC_DATA_2022.xlsx", sheet = 1) %>% subset(Region == "Pacific")
# Convert pairs to individuals
pac_cc$Count[pac_cc$Count_Type == "Pairs"] = pac_cc$Count[pac_cc$Count_Type == "Pairs"]*2
pac_cc$Count_Type <- "Individuals"
pac_species = unique(pac_cc$Species)

# For each species, plot time series
for (spp in pac_species){

  # Extract relevant data
  spdat = subset(pac_cc, Species == spp)
  if (nrow(spdat) == 0) next

  # Plot
  colors = viridis(10)[c(3,6,9)]
  plot = ggplot(spdat,aes(x = Year, y = Count, col = Count_Type)) +
    geom_vline(data = subset(spdat, Count == 0),aes(xintercept = Year), col = "red", size = 1.5, alpha = 0.3)+
    geom_point()+
    xlab("Year")+
    ggtitle(spp)+
    facet_wrap(Colony_Name~., scales = "free_y")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  pdf(paste0("output/Pacific/dataviz/CompleteCounts_",spp,".pdf"), width = 20, height = 10)
  print(plot)
  dev.off()

}


# ------------------------------------------------
# "Plot-based" counts
# ------------------------------------------------

pac_plot = read_xlsx("data/PACIFIC_DATA_2022.xlsx", sheet = 2) %>% subset(Region == "Pacific") %>%
  subset(Species != "unidentified")
pac_plot$Count <- as.numeric(pac_plot$Count)
pac_plot <- na.omit(pac_plot)

# Temporarily omit STPE (not sure if plots are appropriate for this species)
pac_plot <- subset(pac_plot, Species != "STPE")

pac_species = unique(pac_plot$Species)

# Summary of counts
pp_sum <- pac_plot %>%
  group_by(Species,Colony_Name,Year,Plot_Type,Plot) %>%
  summarize(n = n(),
            mean = mean(Count))

# For each species, plot time series
for (spp in pac_species){
  
  # Extract relevant data
  spdat = subset(pac_plot, Species == spp)
  
  # Match plot type to species
  if (spp == "ANMU") spdat = subset(spdat, Plot_Type == "spplots: anmuplots")
  if (spp == "CAAU") spdat = subset(spdat, Plot_Type == "spplots: caauplots")
  if (spp == "RHAU") spdat = subset(spdat, Plot_Type == "spplots: rhauplots")
  if (spp == "TUPU") spdat = subset(spdat, Plot_Type == "spplots: tupuplots")
  
  if (nrow(spdat) == 0) next
  
  spdat_mean <- spdat %>%
    group_by(Colony_Name,Year) %>%
    summarize(mean_count = mean(Count),
              n_plots = n())
  
  # Plot
  plot = ggplot(spdat_mean,aes(x = Year, y = mean_count, size = n_plots)) +
    geom_vline(data = subset(spdat_mean, mean_count == 0),aes(xintercept = Year), col = "red", size = 1.5, alpha = 0.3)+
    geom_point()+
    xlab("Year")+
    ylab("Mean count per plot")+
    ggtitle(spp)+
    facet_wrap(Colony_Name~., scales = "free_y")+
    scale_size_continuous(name = "Number of plots", limits = c(1,max(spdat_mean$n_plots)))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  pdf(paste0("output/Pacific/dataviz/PlotSurveys_",spp,".pdf"), width = 20, height = 10)
  print(plot)
  dev.off()
  
}