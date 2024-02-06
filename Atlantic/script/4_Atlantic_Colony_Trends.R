# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------
my_packs = c('tidyverse','readxl','RColorBrewer','viridis','jagsUI','mgcv','ggrepel','scales',
             'ggpubr')

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

atl = read_xlsx("data/ATLANTIC_QC_AR_DATA_2023_updated.xlsx", sheet = 1)

# Convert pairs to individuals
atl$Count[atl$Count_Type == "Pairs"] = atl$Count[atl$Count_Type == "Pairs"]*2
atl$Count_Type[atl$Count_Type == "Pairs"] = "Individuals"
atl = subset(atl, Count_Type == "Individuals") # Removes several AOS rows

atl_species = unique(atl$Species)
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

# ************************************************************************************************
# Select a species to generate plots for
# ************************************************************************************************

spp <- "ATPU"

# Load model results
file = paste0("output/model_results/",spp,".RData")
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

# ------------------------------------------------
# Decide which colonies to include in regional summary
# ------------------------------------------------

colonies_to_include_for_regional = colonies_to_include %>%
  subset(n_surveys >= 2 &
           first_survey <= (empirical_start+10) &
           last_survey >= (empirical_end-10))
colonies_to_include_for_regional

# ------------------------------------------------
# Plots for each colony
# ------------------------------------------------

# Vector of colonies to plot
colonies_to_plot <- colonies_to_include_for_regional$Colony_Name

colony_samples <- subset(fit_samples, Colony_Name %in% colonies_to_plot)

colony_summary <- colony_samples %>%
  group_by(Colony_Name, Year) %>%
  summarize(q025 = quantile(N_pred,0.025),
            q50 = quantile(N_pred,0.500),
            mean = mean(N_pred),
            q975 = quantile(N_pred,0.975))

colony_plots_freeaxis <- list()

for (colony in colonies_to_plot){
  
  cdat <- subset(spdat, Colony_Name == colony)
  colony_pred <- subset(colony_summary, 
                        Colony_Name == colony & 
                          Year >= min(cdat$Year) & 
                          Year <= max(cdat$Year))
  
  samples_to_plot <- subset(fit_samples,
                            Colony_Name == colony &
                              samp %in% sample(unique(fit_samples$samp),500) &
                              Year >= min(cdat$Year) & 
                              Year <= max(cdat$Year)) 
  
  
  cplot <- ggplot() +
    geom_line(data = samples_to_plot, aes(x = Year, y = N_pred, col = factor(samp)), alpha = 0.05)+
    geom_ribbon(data = colony_pred, aes(x = Year, ymin=q025, ymax = q975),
                alpha = 0.1, fill = "dodgerblue", col = "black", linetype = 2, size = 0.1)+
    geom_line(data = colony_pred, aes(x = Year, y=q50), linewidth = 1, col = "black", alpha = 0.8)+
    geom_point(data = cdat, aes(x = Year, y = Count), col = "black", size = 4)+
    scale_y_continuous(labels = comma, limit = c(0,max(c(colony_pred$q975,cdat$Count, samples_to_plot$N_pred))))+
    scale_x_continuous(limits = range(colony_summary$Year))+
    scale_color_manual(values=rep("blue",length(unique(samples_to_plot$samp))), guide = "none")+
    ggtitle(colony)+
    theme(axis.text = element_text(size=14))
  
  if (colony == colonies_to_plot[1]){
    cplot <- cplot +
      ylab("Estimate of Abundance")+
      xlab("Year")
  } else{
    cplot <- cplot +
      ylab("")+
      xlab("")
  }
  
  colony_plots_freeaxis[[colony]] <- cplot
  
  
}

colony_plots_freeaxis <- ggarrange(plotlist = colony_plots_freeaxis, align = "hv")
png(paste0("output/model_results/figures/",spp,"_colonies_freeaxis.png"), width = 16, height = 16*(5/9), units = "in", res = 600)
print(colony_plots_freeaxis)
dev.off()

colony_plots_fixedaxis <- list()
ylim <- c(0,max(subset(spdat, Colony_Name %in% colonies_to_plot)$Count)*1.5)

for (colony in colonies_to_plot){
  
  cdat <- subset(spdat, Colony_Name == colony)
  colony_pred <- subset(colony_summary, 
                        Colony_Name == colony & 
                          Year >= min(cdat$Year) & 
                          Year <= max(cdat$Year))
  
  samples_to_plot <- subset(fit_samples,
                            Colony_Name == colony &
                              samp %in% sample(unique(fit_samples$samp),500) &
                              Year >= min(cdat$Year) & 
                              Year <= max(cdat$Year)) 
  
  
  cplot <- ggplot() +
    geom_line(data = samples_to_plot, aes(x = Year, y = N_pred, col = factor(samp)), alpha = 0.05)+
    geom_ribbon(data = colony_pred, aes(x = Year, ymin=q025, ymax = q975),
                alpha = 0.1, fill = "dodgerblue", col = "black", linetype = 2, size = 0.1)+
    geom_line(data = colony_pred, aes(x = Year, y=q50), linewidth = 1, col = "black", alpha = 0.8)+
    geom_point(data = cdat, aes(x = Year, y = Count), col = "black", size = 4)+
    scale_y_continuous(labels = comma, limit = ylim)+
    scale_x_continuous(limits = range(colony_summary$Year))+
    scale_color_manual(values=rep("blue",length(unique(samples_to_plot$samp))), guide = "none")+
    ggtitle(colony)+
    theme(axis.text = element_text(size=14))
  
  if (colony == colonies_to_plot[1]){
    cplot <- cplot +
      ylab("Estimate of Abundance")+
      xlab("Year")
  } else{
    cplot <- cplot +
      ylab("")+
      xlab("")
  }
  
  colony_plots_fixedaxis[[colony]] <- cplot
  
  
}

colony_plots_fixedaxis <- ggarrange(plotlist = colony_plots_fixedaxis, align = "hv")

png(paste0("output/model_results/figures/",spp,"_colonies_fixedaxis.png"), width = 16, height = 16*(5/9), units = "in", res = 600)
print(colony_plots_fixedaxis)
dev.off()

