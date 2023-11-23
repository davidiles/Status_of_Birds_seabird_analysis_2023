# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------
my_packs = c('tidyverse','readxl','scales')

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
                                                        linewidth = 1, fill = NA),
                            axis.line = element_line(colour = "black"),
                            strip.text = element_text(size = 12, colour = "black"),
                            strip.background = element_rect(colour = "black",
                                                            fill = "lightblue2",
                                                            linetype = "solid"),
                            axis.title.y = element_text(margin = margin(0,10,0,0)),
                            axis.title.x = element_text(margin = margin(10,0,0,0)),
                            panel.background = element_rect(fill = "white"))

# ------------------------------------------------
# Read in available data
# ------------------------------------------------

# Arctic database, updated by DTI on 2023-11-05
dat1 = read_xlsx("data/Arctic/ArcticSeabirdColonyDatabase_edited.xlsx", sheet = 1)

# Replicate plot data from Sarah Gutowsky and Tony Gaston (at Prince Leopold Island)
dat2 = read_xlsx("data/Arctic/Arctic_PrinceLeopoldIsland_edited.xlsx", sheet = 3) %>%
  
  # Make consistent with Arctic seabird database colony names
  mutate(Colony_Name = "PRINCE LEOPOLD IS") 

# Data compiled by Adam Smith in 2012
dat3 = read_xlsx("data/Arctic/ArcticSeabird_SOCB2012_edited.xlsx", sheet = 1)

# ********************
# NOTES:
#  - DTI edited the seabird colony database to include new surveys of NOFU in Mallory et al. 2020
# ********************

# ------------------------------------------------
# Plot data in dat1
# ------------------------------------------------

spp_vec <- c("NOFU")

spp <- "BLKI"

sp_dat <- subset(dat1, species == spp)
sp_dat$individuals <- as.numeric(sp_dat$individuals)
sp_dat <- subset(sp_dat, year >= 1970)
max_y <- max(sp_dat$individuals)
max_x <- max(sp_dat$year)

complete_count_plot <- ggplot(sp_dat, aes(x = year, y = individuals, label = comma(individuals)))+
  geom_point() +
  geom_text(size=2, col = "gray50", hjust = 0, vjust=-1) +
  
  facet_wrap(location~.)+
  scale_y_continuous(labels = comma, limits = c(0,max_y*1.5))+
  coord_cartesian(xlim=c(1970,2023))+
  ggtitle(spp)
complete_count_plot

# ------------------------------------------------
# Plot "replicate plot" data at PLI
# ------------------------------------------------

sp_dat <- subset(dat2, Species == spp)
sp_dat$count <- as.numeric(sp_dat$Total_Size_of_Colony)
sp_dat$year <- as.numeric(sp_dat$Year_of_Size_Estimate)

replicate_count_plot <- ggplot(sp_dat, aes(x = year, y = count))+
  geom_point() +
  scale_y_continuous(labels = comma)+
  ggtitle(spp)

replicate_count_plot
