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

# ------------------------------------------------
# Read in data
# ------------------------------------------------

# --------------------
# Complete 'island total' counts
# --------------------

ccounts = read_xlsx("data/PACIFIC_DATA_2023.xlsx", sheet = 1) %>% subset(Region == "Pacific")

# Convert pairs to individuals
ccounts$Count[ccounts$Count_Type == "Pairs"] = ccounts$Count[ccounts$Count_Type == "Pairs"]*2
ccounts$Count_Type <- "Individuals"

spp_vec_ccounts = unique(ccounts$Species)

# --------------------
# Plot-based counts
# --------------------

pcounts = read_xlsx("data/PACIFIC_DATA_2023.xlsx", sheet = 2) %>% subset(Region == "Pacific") %>%
  subset(Species != "unidentified")
pcounts$Count <- as.numeric(pcounts$Count)
pcounts <- na.omit(pcounts)

# Temporarily omit STPE (not sure if plots are appropriate for this species)
pcounts <- subset(pcounts, Species != "STPE")

pcounts <- subset(pcounts,
                   (Species == "ANMU" & Plot_Type == "spplots: anmuplots") |
                     (Species == "CAAU" & Plot_Type == "spplots: caauplots") |
                     (Species == "RHAU" & Plot_Type == "spplots: rhauplots") |
                     (Species == "TUPU" & Plot_Type == "spplots: tupuplots") )

spp_vec_pcounts = unique(pcounts$Species)

# ------------------------------------------------
# Loop through species and conduct analysis
# ------------------------------------------------

spp_vec <- union(spp_vec_ccounts, spp_vec_pcounts)
#spp_vec = spp_vec[-which(spp_vec %in% c("BLOY","STPE","GWGU","PECO"))]
spp_vec <- c("ANMU","CAAU","RHAU","TUPU")
for (spp in spp_vec){
  
  #if (file.exists(paste0("output/Pacific/model_results/",spp,".RData"))) next
  print(spp)
  spdat_ccount = spdat_pcount = NULL
  
  # ------------------------------------------------
  # Extract "plot based" survey data, if it exists for this species
  # ------------------------------------------------
  
  spdat_pcount <- subset(pcounts,
                    (Species == "ANMU" & Plot_Type == "spplots: anmuplots") |
                      (Species == "CAAU" & Plot_Type == "spplots: caauplots") |
                      (Species == "RHAU" & Plot_Type == "spplots: rhauplots") |
                      (Species == "TUPU" & Plot_Type == "spplots: tupuplots") ) %>%
    subset(Species == spp) %>%
    group_by(Species, Colony_Name, Year) %>%
    summarize(Count = sum(Count),
              n_plots = length(Plot))

  # ------------------------------------------------
  # Extract "total colony abundance" estimates, if it exists
  # ------------------------------------------------
  
  spdat_ccount = subset(ccounts, Species == spp)
  spdat_ccount$n_plots <- 1
  
  # ------------------------------------------------
  # Extract "total colony abundance" estimates, if it exists
  # ------------------------------------------------
  
  # If there are plot-based counts for an island, use those in model.  Otherwise, use "complete" colony counts
  if (nrow(spdat_pcount)>0) spdat_ccount <- subset(spdat_ccount, Colony_Name %!in% spdat_pcount$Colony_Name)
  
  # ------------------------------------------------
  # Join plot-based and 'colony-level' survey data
  #     - select which colonies to include in analysis
  # ------------------------------------------------
  
  spdat_ccount$Survey_Type <- "Colony Counts"
  spdat_pcount$Survey_Type <- "Plot Counts"
  spdat <- dplyr::bind_rows(spdat_ccount,spdat_pcount) %>%
    dplyr::select(-'Source (PI) / point of contact',-'Comments',-'Region')
  
  colonies_to_include = spdat %>%
    group_by(Colony_Name) %>%
    summarize(mean_count = mean(Count),
              first_survey = min(Year),
              last_survey = max(Year),
              n_surveys = length(unique(Year)),
              Survey_Type = unique(Survey_Type)) %>%
    subset(n_surveys > 1)
    
  if (nrow(colonies_to_include) == 0) next
  
  spdat <- subset(spdat, Colony_Name %in% colonies_to_include$Colony_Name)
  
  # ------------------------------------------------
  # Range of years for analysis
  # ------------------------------------------------
  
  first_year = min(spdat$Year)
  final_year = max(spdat$Year)
  
  # ------------------------------------------------
  # Data viz
  # ------------------------------------------------
  
  ggplot(spdat)+
    geom_point(aes(x = Year, y = Count))+
    scale_y_continuous(labels = comma)+
    facet_wrap(Colony_Name~., scales = "free_y")
  
  # ------------------------------------------------
  # Tables linking colony names to colony numbers, and year to "numeric years"
  # ------------------------------------------------
  
  # Tables that link colony numbers to colony names
  spdat$colony_numeric <- as.integer(factor(spdat$Colony_Name))
  colony_name_table = unique(spdat[,c("Colony_Name","colony_numeric")]) %>% arrange(colony_numeric)
  
  # Tables that link year index to actual year
  year_table = data.frame(Year = min(spdat$Year):2023)
  year_table$year_numeric = 1:nrow(year_table)
  
  # ------------------------------------------------
  # Fit model in JAGS and save
  # ------------------------------------------------
  
  # Data for import into jags
  nknots = 5
  year <- spdat$Year - min(year_table$Year) + 1
  ymax <- max(year_table$year_numeric)
  nyears = length(1:ymax)
  colony = spdat$colony_numeric
  ncolony <- max(colony)
  count <- spdat$Count
  ncounts = length(count)
  log_offset = log(spdat$n_plots)
  
  # Use jagam to prepare basis functions
  nyearspred = length(1:ymax)
  preddat = data.frame(yrs = 1:ymax,count = 1)
  form = as.formula(paste("count ~ s(yrs,k =",nknots,")"))
  gamprep = jagam(formula = form,
                  data = preddat,
                  file = "scripts/2-Pacific/tempgam.txt",
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
                   year = year,
                   log_offset = log_offset)
  
  # Fit model using JAGS
  parameters.to.save = c("sdnoise","sdbeta","C","beta.X","log_mu","eta","population_index")
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter = 1100000,
              n.burnin = 100000,
              n.thin = 1000,
              model.file = "scripts/2-Pacific/Pacific_GAMM.jags",
              n.chains = 3,
              parallel = TRUE)
  
  # Save results for each species
  save.image(paste0("output/Pacific/model_results/",spp,".RData"))
  
}

