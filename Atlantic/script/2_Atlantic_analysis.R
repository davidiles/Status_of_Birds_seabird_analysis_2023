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
# Read in data
# ------------------------------------------------

atl = read_xlsx("data/ATLANTIC_QC_AR_DATA_2023_updated.xlsx", sheet = 1)

# Convert pairs to individuals
atl$Count[atl$Count_Type == "Pairs"] = atl$Count[atl$Count_Type == "Pairs"]*2
atl$Count_Type[atl$Count_Type == "Pairs"] = "Individuals"
atl = subset(atl, Count_Type == "Individuals") # Removes several AOS rows

atl_summary = atl %>%
  group_by(Species,Colony_Name,Count_Type) %>%
  summarize(first_year = min(Year),
            last_year = max(Year),
            n_counts = length(unique(Year)),
            min_count = min(Count),
            max_count = max(Count),
            mean_count = mean(Count))

write.csv(atl_summary, "output/Atlantic_summary.csv",row.names = FALSE)

atl_species = unique(atl$Species) %>% sort()

# ------------------------------------------------
# Loop through species and conduct analysis
# ------------------------------------------------

for (spp in atl_species){
  
  print(spp)
  #if (file.exists(paste0("output/model_results/",spp,".RData"))) next
  
  # Extract relevant data
  spdat = subset(atl, Species == spp)
  
  first_year = min(spdat$Year)
  final_year = max(spdat$Year)
  
  colonies_to_include = spdat %>%
    group_by(Colony_Name) %>%
    summarize(mean_count = mean(Count),
              first_survey = min(Year),
              last_survey = max(Year),
              n_surveys = length(unique(Year)))
  
  # ------------------------------------------------
  # Some custom data cleaning for certain species
  # ------------------------------------------------
  
  if (spp != "ROST" & spp != "NOFU" & spp != "CATE"){
    # Remove the smallest colonies, which are strongly inflating uncertainty
    colonies_to_include = colonies_to_include %>% subset(mean_count > 100)
  }
  
  if (spp == "LESP"){
    # Omit the smallest LESP colonies that are strongly influencing uncertainty
    colonies_to_include = subset(colonies_to_include, mean_count > 1000)
    spdat <- subset(spdat, !(Year == 1978 & Colony_Name == "Bon Portage, NS"))
  }
  
  spdat = subset(spdat, Colony_Name %in% colonies_to_include$Colony_Name)
  
  # Skip species if there is not enough data
  if (nrow(spdat) == 0 | length(unique(spdat$Colony_Name)) < 2) next
  
  # Tables that link colony numbers to colony names
  spdat$colony_numeric <- as.integer(factor(spdat$Colony_Name))
  colony_name_table = unique(spdat[,c("Colony_Name","colony_numeric")])
  
  # Tables that link year index to actual year
  year_table = data.frame(Year = min(spdat$Year):max(spdat$Year))
  year_table$year_numeric = 1:nrow(year_table)
  
  # ------------------------------------------------
  # Fit model in JAGS and save
  # ------------------------------------------------
  
  # Data for import into jags
  nknots = 4
  year <- spdat$Year - min(year_table$Year) + 1
  ymax <- max(year_table$year_numeric)
  nyears = length(1:ymax)
  colony = spdat$colony_numeric
  ncolony <- max(colony)
  count <- spdat$Count
  ncounts = length(count)
  
  # Use jagam to prepare basis functions
  nyearspred = length(1:ymax)
  preddat = data.frame(yrs = 1:ymax,count = 1)
  form = as.formula(paste("count ~ s(yrs,k =",nknots,")"))
  gamprep = jagam(formula = form,
                  data = preddat,
                  file = "script/tempgam.txt",
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
                   year = year)
  
  # Fit model using JAGS
  parameters.to.save = c("sdnoise","sdbeta","C","beta.X","population_index")
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter = 1100000,
              n.burnin = 100000,
              n.thin = 1000,
              model.file = "script/Atlantic_GAMM.jags",
              n.chains = 3,
              parallel = TRUE) # 12 min
  
  # Save results for each species
  save.image(paste0("output/model_results/",spp,".RData"))
  
}

