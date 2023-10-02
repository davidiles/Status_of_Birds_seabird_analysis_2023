# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------
my_packs = c('tidyverse','readxl','RColorBrewer','viridis','jagsUI','mgcv','ggrepel')

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

# For each species, plot time series
for (spp in atl_species){

  # Extract relevant data
  spdat = subset(atl, Species == spp)
  if (nrow(spdat) == 0) next

  # Plot
  colors = viridis(10)[c(3,6,9)]
  plot = ggplot(spdat,aes(x = Year, y = Count, col = Count_Type)) +
    geom_vline(data = subset(spdat, Count == 0),aes(xintercept = Year), col = "red", size = 1.5, alpha = 0.3)+
    geom_point()+
    xlab("Year")+
    ggtitle(spp)+
    scale_color_manual(values = colors,name = "Count type", drop = FALSE)+
    facet_wrap(Colony_Name~., scales = "free_y")

  pdf(paste0("output/Atlantic/dataviz/",spp,".pdf"), width = 20, height = 10)
  print(plot)
  dev.off()

}