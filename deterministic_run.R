source("R/read.report.R")
source("R/gettables.R")
library(tidyverse)
library(readxl)
library(hydradata)
source("R/read.report.R")  # make sure read.report2.R and gettables.R are sourced
source("R/gettables.R")

#rm(list = ls())

# 1. Load your operating model inputs and outputs
hydraDataList <- readRDS("Sarah_files/hydra_sim_GBself_5bin.rds")
hydraDataList$observedBiomass<- hydraDataList[["observedBiomass"]]%>%  
  filter(survey==1)
hydraDataList$observedSurvSize<-hydraDataList[["observedSurvSize"]]%>%  
  filter(survey==1)
hydraDataList$observedSurvDiet<-hydraDataList[["observedSurvDiet"]]%>%  
  filter(species %in% c(1, 4))
hydraDataList$Nsurveys<-(Nsurveys=1)
hydraDataList$Nsurvey_obs<-(Nsurvey_obs=166)
hydraDataList$Nsurvey_size_obs<-(Nsurvey_size_obs=166)
hydraDataList$Ndietprop_obs<-(Ndietprop_obs=89)
hydraDataList[["observedSurvDiet"]][["inpN"]]<-(inpN=25)
hydraDataList[["observedBiomass"]][["cv"]]<-(cv=0.2)
hydraDataList[["observedCatch"]][["cv"]]<-(cv=0.05)
hydraDataList[["observedCatchSize"]][["inpN"]]<-(inpN=25)
hydraDataList[["observedSurvSize"]][["inpN"]]<-(inpN=25)

repfile <- "OM_scenarios/OM_case1/hydra_sim.rep"
output <- reptoRlist(repfile)

obs_surveyB <- hydraDataList$observedBiomass %>% 
  filter(survey==1)
obs_catchB <- hydraDataList$observedCatch

biorows <- dim(obs_surveyB)[1]
catchrows <- dim(obs_catchB)[1]

#create a table with estimated and observed values
indexfits <- gettables(repfile, biorows, catchrows)

# select survey biomass
biomass<-indexfits[[1]] %>%
  mutate(species = hydraDataList$speciesList[species])
# select catch
catch<-indexfits[[2]]%>%
  mutate(species = hydraDataList$speciesList[species])

# create simulation object
set.seed(23)
sim_data <- NULL
isim <- 1

for (isim in 1:1) {  
  
  # replace index with simulated data
  
  ##### REPLACE CATCH DATA FROM OM estimates ######
  hydraDataList[["observedCatch"]][["catch"]]<- catch$pred_catch
   sim_data[[isim]] <- hydraDataList
  
 ##### REPLACE SURVEY BIOMASS DATA FROM OM estimates######
  # if you have one survey
  hydraDataList[["observedBiomass"]][["biomass"]]<- biomass$pred_bio 
   sim_data[[isim]] <- hydraDataList
  
  #### REPLACE CATCH SIZE COMPOSITION DATA FROM OM estimates####
  
  obs_catch <- hydraDataList$observedCatchSize %>% tibble()
  obs_catch<-obs_catch %>% pivot_longer(cols=7:ncol(.), names_to = "lenbin") %>%
    mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")),
           label = rep("catch",nrow(.)))# %>% filter(value != -999)
  obs_catch$value[which(obs_catch$value == -999)] = 0.00001
  
  pred_catch<-output$pred_catch_size
  obs_catch$pred_catch<-pred_catch
  
  obs_catch <-obs_catch %>% 
    mutate(value=pred_catch) %>%
    select(-pred_catch)
 
  # Convertir de nuevo a formato ancho 
  hydraDataList$observedCatchSize <- obs_catch %>%
    pivot_wider(names_from = lenbin, values_from = value,
                names_prefix = "sizebin") %>%
    arrange(fishery, species, year) %>%
    select(-label)
  
  
 #### SIMULATE SURVEY SIZE COMPOSITION DATA ####
  # remove  %>% filter(survey==1) if you have 2 surveys
  
  sample_size<-25
  
  obsdiet_comp <- hydraDataList$observedSurvDiet %>% tibble()
  diet_combined <- obsdiet_comp %>%
    group_by(year, species, sizebin) %>%
    mutate(w = inpN) %>%
    summarise(
      InpN = sum(w),
      across(c(Atlantic_cod, Atlantic_herring, Atlantic_mackerel, Spiny_dogfish, allotherprey),
             ~ weighted.mean(.x, w, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(survey = 1) %>%
    relocate(survey, .before = year)   # put survey back in the first column
  
  colnames(diet_combined) <- c("survey", "year", "species", "sizebin",
                               "inpN", "Atlantic_cod", "Atlantic_herring", 
                               "Atlantic_mackerel", "Spiny_dogfish", "allotherprey")
  
  obsdiet_comp <- diet_combined %>%
    pivot_longer(
      cols = c(Atlantic_cod, Atlantic_herring, Atlantic_mackerel, Spiny_dogfish, allotherprey),
      names_to = "prey",
      values_to = "value"
    ) %>%
    mutate(
      value = ifelse(value == -999, 0.000001, value)
    )
  
  simulated_diet <- obsdiet_comp %>%
    group_by(species, year, sizebin) %>%
    mutate(
      norm_value = value / sum(value)  # in case values don’t sum exactly to 1
    ) %>%
    summarise(
      prey = prey,
      sim_matrix = rmultinom(1, size = sample_size, prob = norm_value),
      .groups = "drop"
    ) %>%
    mutate(
      sim_count = as.vector(sim_matrix),
      simulated_prop = sim_count / sample_size
    ) %>%
    select(species, year, sizebin, prey, simulated_prop)
  
  
  obsdiet_comp$value = simulated_diet$simulated_prop
  hydraDataList$observedSurvDiet <- obsdiet_comp %>%
    mutate(inpN = 25) %>%
    pivot_wider(names_from = "prey")
  
  #write.csv(hydraDataList$observedSurvDiet, file = "survvvvv.csv", row.names = T)
  sim_data[[isim]] <- hydraDataList
}

# Save the deterministic dataset
saveRDS(hydraDataList, "sims/deterministic/deterministic_sim_data.rds")

#### WRITE tsDat FUNCTION ####
source("R/write_tsDatFile.R")
source("R/read.report.R")
nsim<-1 

hydraDataList2 <- readRDS("sims/deterministic/deterministic_sim_data.rds")

listOfParameters<-list()
listOfParameters$outDir<-paste0(getwd(),"/","sims","/","deterministic","/")
listOfParameters$outputFilename<-"hydra_sim_deterministic"
listOfParameters$fillLength <- 2000

write_tsDatFile(hydraDataList2,listOfParameters)

######### deterministic estimability fits

library(dplyr)
library(tidyverse)
library(ggplot2)
source("R/write_tsDatFile.R")
source("R/read.report.R")

###############################
###### PLOTS #################
#############################

## FISHING MORTALITY OM estimate

repfile_OM <- "OM_scenarios/OM/hydra_sim.rep"
output_OM<-reptoRlist(repfile_OM)


repfile_DET <- "sims/deterministic/hydra_sim.rep"
output_DET<-reptoRlist(repfile_DET)

stepperyr <- output_OM$Nstepsyr
if (length(stepperyr)==0) stepperyr <- nrow(output_OM$EstBsize)/hydraDataList$Nyrs/length(hydraDataList$speciesList)

nlen<-ncol(output_OM$EstBsize)
est_SSB_OM <- output_OM$EstBsize %>% 
  as.data.frame() %>% 
  pivot_longer(cols=1:ncol(.), names_to = "ilen", names_prefix = "V") %>% 
  mutate(species = rep(hydraDataList$speciesList, each = hydraDataList$Nyrs*nlen*stepperyr),
         year  = rep(rep(1:(hydraDataList$Nyrs*stepperyr),each=nlen), length(hydraDataList$speciesList)),
         year = (1-1/stepperyr) + year / stepperyr) %>%
  filter(ilen == 5) 

est_SSB_DET<- output_DET$EstBsize %>% 
  as.data.frame() %>% 
  pivot_longer(cols=1:ncol(.), names_to = "ilen", names_prefix = "V") %>% 
  mutate(species = rep(hydraDataList$speciesList, each = hydraDataList$Nyrs*nlen*stepperyr),
         year  = rep(rep(1:(hydraDataList$Nyrs*stepperyr),each=nlen), length(hydraDataList$speciesList)),
         year = (1-1/stepperyr) + year / stepperyr) %>%
  filter(ilen == 5)

# Add source label
est_SSB_OM <- est_SSB_OM %>% mutate(source = "OM")
est_SSB_DET <- est_SSB_DET %>% mutate(source = "Deterministic")

# Combine
est_SSB_all <- bind_rows(est_SSB_OM, est_SSB_DET)

# Plot
ggplot(est_SSB_all, aes(x = year, y = value, color = source, linetype = source)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    title = "Spawning Stock Biomass (SSB): OM vs Deterministic",
    x = "Year",
    y = "SSB",
    color = "Run",
    linetype = "Run"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")




est_F_OM <- output_OM$Fyr %>% 
  as.data.frame() %>% 
  pivot_longer(cols=3:ncol(.), names_to = "year", names_prefix = "V") %>% 
  rename(species = "V1",
         fleet = "V2") %>% 
  mutate(year = as.numeric(year)-2,
         species = hydraDataList$speciesList[species])

est_F_DET<- output_DET$Fyr %>% 
  as.data.frame() %>% 
  pivot_longer(cols=3:ncol(.), names_to = "year", names_prefix = "V") %>% 
  rename(species = "V1",
         fleet = "V2") %>% 
  mutate(year = as.numeric(year)-2,
         species = hydraDataList$speciesList[species])


est_F_OM <- est_F_OM %>% mutate(source = "OM")
est_F_DET <- est_F_DET %>% mutate(source = "Deterministic")

# Combine into one dataframe
F_combined <- bind_rows(est_F_OM, est_F_DET)
F_combined <- F_combined %>%
  filter(value > 1e-10)

# Plot
F_OM_DET<-ggplot(F_combined, aes(x = year, y = value, color = source, linetype = source)) +
  geom_line(size = 1) +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    title = "Fishing Mortality: Operating Model vs Deterministic Run",
    x = "Year",
    y = "Fishing Mortality (F)",
    color = "Scenario",
    linetype = "Scenario"
  ) +
  theme_minimal()

# Reshape to wide format: one column for each source
F_scatter <- F_combined %>%
  select(species, year, source, value) %>%
  pivot_wider(names_from = source, values_from = value)

# Plot OM vs Deterministic F
ggplot(F_scatter, aes(x = OM, y = Deterministic)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
  facet_wrap(~ species, scales = "free") +
  labs(
    title = "OM vs Deterministic Estimates of Fishing Mortality",
    x = "OM Fishing Mortality (F)",
    y = "Deterministic Fishing Mortality (F)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

#ggsave(filename = "F_OM_DET.jpeg",
#      plot = F_OM_DET,
#     width = 10, height = 8, units = "in", dpi = 300)


##########################
#### RECRUITMENT

est_rec_OM <- output_OM$EstRec %>% 
  as.data.frame() %>% 
  mutate(species = hydraDataList$speciesList) %>% 
  select(species, everything()) %>% 
  pivot_longer(cols = -species, names_to = "year") %>% 
  mutate(year = as.integer(str_remove(year, "V")),
         log_rec = ifelse(value > 0,log(value),NA))

est_rec_DET <- output_DET$EstRec %>% 
  as.data.frame() %>% 
  mutate(species = hydraDataList$speciesList) %>% 
  select(species, everything()) %>% 
  pivot_longer(cols = -species, names_to = "year") %>% 
  mutate(year = as.integer(str_remove(year, "V")),
         log_rec = ifelse(value > 0,log(value),NA))

est_rec_OM <- est_rec_OM %>% mutate(source = "OM")
est_rec_DET <- est_rec_DET %>% mutate(source = "Deterministic")

# Combine into one dataframe
rec_combined <- bind_rows(est_rec_OM, est_rec_DET)

# Plot
REC_OM_DET<-ggplot(rec_combined, aes(x = year, y = value, color = source, linetype = source)) +
  geom_line(size = 1) +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    title = "Estimated Recruitment: Operating Model vs Deterministic Run",
    x = "Year",
    y = "Recruitment (R)",
    color = "Scenario",
    linetype = "Scenario"
  ) +
  theme_minimal()

#ggsave(filename = "REC_OM_DET.jpeg",
 #     plot = REC_OM_DET,
  #   width = 10, height = 8, units = "in", dpi = 300)

# Reshape to wide format: one column for each source
rec_scatter <- rec_combined %>%
  select(species, year, source, value) %>%
  pivot_wider(names_from = source, values_from = value)

# Plot OM vs Deterministic Recruitment
ggplot(rec_scatter, aes(x = OM, y = Deterministic)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
  facet_wrap(~ species, scales = "free") +
  labs(
    title = "OM vs Deterministic Estimates of Recruitment",
    x = "OM Recruitment",
    y = "Deterministic Recruitment"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

##########################
#### SPAWNING BIOMASS

# Read the full file
lines_OM <- readLines(repfile_OM)
lines_DET <- readLines(repfile_DET)

# Find the header line
start_line_OM <- grep("Year Species SSB", lines_OM)
start_line_DET <- grep("Year Species SSB", lines_DET)

n_rows<-(hydraDataList[["Nyrs"]]*hydraDataList[["Nspecies"]])

# Read the data starting just after the header
ssb_data_OM <- read.table(text = lines_OM[(start_line_OM + 1):(start_line_OM + n_rows)],
                          col.names = c("Year", "Species", "SSB"))                          

ssb_data_DET <- read.table(text = lines_DET[(start_line_DET + 1):(start_line_DET + n_rows)],
                          col.names = c("Year", "Species", "SSB"))                          

ssb_data_OM <- ssb_data_OM %>% mutate(Source = "OM")
ssb_data_DET <- ssb_data_DET %>% mutate(Source = "Deterministic")

# Combine both datasets
ssb_combined <- bind_rows(ssb_data_OM, ssb_data_DET)
ssb_combined <- ssb_combined %>%
  mutate(Species = hydraDataList$speciesList[Species])
# Optional: add species names
# ssb_combined$Species <- hydraDataList$speciesList[ssb_combined$Species]

# Plot
SSB_OM_DET<-ggplot(ssb_combined, aes(x = Year, y = SSB, color = Source, linetype = Source)) +
  geom_line(size = 1) +
  facet_wrap(~ Species, scales = "free_y") +
  labs(
    title = "Spawning Stock Biomass (SSB) by Species: OM vs Deterministic",
    x = "Year",
    y = "SSB (units)",
    color = "Model",
    linetype = "Model"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

#ggsave(filename = "SSB_OM_DET.jpeg",
#     plot = SSB_OM_DET,
#   width = 10, height = 8, units = "in", dpi = 300)


# Reshape to wide format: one column for each source
ssb_scatter <- ssb_combined %>%
  select(Species, Year, Source, SSB) %>%
  pivot_wider(names_from = Source, values_from = SSB)

# Plot OM vs Deterministic SSB
ggplot(ssb_scatter, aes(x = OM, y = Deterministic)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
  facet_wrap(~ Species, scales = "free") +
  labs(
    title = "OM vs Deterministic Estimates of Spawning Stock Biomass (SSB)",
    x = "OM SSB",
    y = "Deterministic SSB"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

#########################
####### depletion
#########################

ssb_status <- ssb_combined %>%
  group_by(Species, Source) %>%
  summarize(
    SSB_first = round(SSB[Year == min(Year)], 2),
    SSB_last = round(SSB[Year == max(Year)], 2),
    Status = round(SSB_last / SSB_first, 3),
    .groups = "drop"
  )


ssb_summary <- ssb_combined %>%
  pivot_wider(names_from = Source, values_from = SSB) %>%
  group_by(Species) %>%
  summarize(
    MeanDiff = mean(Deterministic - OM, na.rm = TRUE),
    RelBias  = mean((Deterministic - OM) / OM, na.rm = TRUE) * 100,
    RMSE     = sqrt(mean((Deterministic - OM)^2, na.rm = TRUE)),
    Correlation = cor(OM, Deterministic, use = "complete.obs")
  )


#MeanDiff: raw difference across years
#RelBias: average % over- or under-estimation
#RMSE: error magnitude
#Correlation: shape similarity over time


# Function to read a .par file with vector-style parameter blocks
source("R/read_par.R")

# Read both .par files
par_OM <- read_par_blocked("OM_scenarios/OM/hydra_sim.par")
par_DET <- read_par_blocked("sims/deterministic/hydra_sim.par")

# Compare shared parameters
#common_params <- intersect(names(par_OM), names(par_DET))
# Choose only specific parameters to compare
selected_params <- c("ln_yr1N", "ln_avg_recruitment", "recruitment_devs", "avg_F", "fishsel_pars", "ln_fishery_q", "survey_selpars")  # add any others here
common_params <- intersect(selected_params, intersect(names(par_OM), names(par_DET)))

param_comparison <- purrr::map_dfr(common_params, function(pname) {
  x <- par_OM[[pname]]
  y <- par_DET[[pname]]
  
  # Ensure matching length
  len <- min(length(x), length(y))
  x <- x[1:len]
  y <- y[1:len]
  
  tibble(
    Parameter = pname,
    Index = seq_len(len),
    OM = x,
    DET = y,
    Diff = y - x,
    RelError = (y - x) / x * 100
  )
})

# Summary by parameter
param_summary <- param_comparison %>%
  group_by(Parameter) %>%
  summarize(
    N = n(),
    MeanDiff = mean(Diff),
    MeanRelError = mean(RelError, na.rm = TRUE),
    RMSE = sqrt(mean((Diff)^2)),
    .groups = "drop"
  )

# Output example
knitr::kable(head(param_summary, 10), digits = 4, caption = "Parameter Error Summary (first 10)")







































#Read deterministic outputs

# Replace with deterministic .rep file path
repfile_det <- "sims/deterministic/hydra_sim.rep"  # Example: first simulation
OMrepfile <- "OM_scenarios/OM/hydra_sim.rep" # Original Operating Model .rep

# Source helper functions
source("R/read.report2.R")
source("R/gettables.R")

# Read Operating Model outputs
biorows <- nrow(hydraDataList$observedBiomass %>% filter(survey == 1))
catchrows <- nrow(hydraDataList$observedCatch)

OM_indexfits <- gettables(OMrepfile, biorows, catchrows)
DET_indexfits <- gettables(repfile_det, biorows, catchrows)

# Extract survey and catch
OM_survey <- OM_indexfits[[1]]
OM_catch <- OM_indexfits[[2]]

DET_survey <- DET_indexfits[[1]]
DET_catch <- DET_indexfits[[2]]

#Compare Catch Fits
catch_comp <- OM_catch %>%
  mutate(
    OM_catch = pred_catch,
    est_catch = DET_catch$pred_catch,
    residual = est_catch - OM_catch
  )

# Compare Survey Fits
survey_comp <- OM_survey %>%
  mutate(
    OM_biomass = pred_bio,
    est_biomass = DET_survey$pred_bio,
    residual = est_biomass - OM_biomass
  )

species_names <- hydraDataList$speciesList
catch_comp <- catch_comp %>%
  mutate(Species = species_names[species])  # ¡asignas el nombre usando el índice!

survey_comp <- survey_comp %>%
  mutate(Species = species_names[species])


catch_plot <- catch_comp %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = log(OM_catch), color = "OM Catch"), linetype = "dashed", size = 1) +
  geom_line(aes(y = log(est_catch), color = "EM Catch"), size = 1) +
  facet_wrap(~Species, scales = "free_y", ncol = 2) +
  labs(
    title = "Time Series of OM vs EM estimated catch (Deterministic data)",
    x = "Year",
    y = "log(Catch)",
    color = NULL
  ) +
  scale_color_manual(values = c("OM Catch" = "black", "EM Catch" = "blue")) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )

species_names <- hydraDataList$speciesList

survey_comp <- survey_comp %>%
  mutate(Species = species_names[species])


survey_plot <- survey_comp %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = log(OM_biomass), color = "OM Biomass"), linetype = "dashed", size = 1) +
  geom_line(aes(y = log(est_biomass), color = "EM Biomass"), size = 1) +
  facet_wrap(~Species, scales = "free_y", ncol = 2) +
  labs(
    title = "Time Series of OM vs EM estiated survey biomass (Deterministic data)",
    x = "Year",
    y = "log(Biomass)",
    color = NULL  # Leyenda sin título
  ) +
  scale_color_manual(values = c("OM Biomass" = "black", "EM Biomass" = "green4")) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Nombres de especies bonitos
    legend.position = "bottom",                           # Leyenda abajo
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )
# --- Save Catch Plot ---
#ggsave(filename = "Deterministic_Catch.jpeg",
#      plot = catch_plot,
#     width = 10, height = 8, units = "in", dpi = 300)

# --- Save Survey Biomass Plot ---
#ggsave(filename = "Deterministic_Biomass.jpeg",
#      plot = survey_plot,
#     width = 10, height = 8, units = "in", dpi = 300)


# Quick Metrics

catch_metrics <- catch_comp %>%
  summarise(
    RMSE = sqrt(mean(residual^2)),
    MAE = mean(abs(residual)),
    Bias = mean(residual)
  )

survey_metrics <- survey_comp %>%
  summarise(
    RMSE = sqrt(mean(residual^2)),
    MAE = mean(abs(residual)),
    Bias = mean(residual)
  )

print("Catch Metrics:")
print(catch_metrics)

print("Survey Biomass Metrics:")
print(survey_metrics)

#Quick Residual Histograms

ggplot(catch_comp, aes(x = residual)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
  labs(title = "Catch Residuals (log scale)", x = "Residual (log scale)", y = "Count") +
  theme_minimal()

ggplot(survey_comp, aes(x = residual)) +
  geom_histogram(bins = 30, fill = "green", alpha = 0.5) +
  labs(title = "Survey Biomass Residuals (log scale)", x = "Residual (log scale)", y = "Count") +
  theme_minimal()

# Optional: Residual Time Series

catch_comp %>%
  ggplot(aes(x = year, y = residual, group = interaction(species, fishery))) +
  geom_line() +
  facet_wrap(~species, scales = "free_y") +
  labs(title = "Catch Residuals Over Time", x = "Year", y = "Residual (log scale)") +
  theme_minimal()

survey_comp %>%
  ggplot(aes(x = year, y = residual, group = interaction(species, survey))) +
  geom_line() +
  facet_wrap(~species, scales = "free_y") +
  labs(title = "Survey Biomass Residuals Over Time", x = "Year", y = "Residual (log scale)") +
  theme_minimal()


# --- Read Objective Function Value from .par file ---
parfile_det <- "sims/deterministic/hydra_sim.par"  # Adjust path if needed

# Read the first line only
first_line <- readLines(parfile_det, n = 1)

# Extract the Objective Function Value
objective_value <- as.numeric(sub(".*Objective function value = ([0-9\\.eE\\+-]+).*", "\\1", first_line))

# Extract the Maximum Gradient Component
max_gradient <- as.numeric(sub(".*Maximum gradient component = ([0-9\\.eE\\+-]+).*", "\\1", first_line))

# --- Print results ---
cat("Final Objective Function Value:", objective_value, "\n")
cat("Final Maximum Gradient Component:", max_gradient, "\n\n")

# --- Set thresholds ---
objective_threshold <- 1e-2  # Objective function should be very close to zero
gradient_threshold <- 1e-4   # Gradient should be very small for good convergence

# --- Check Objective Function ---
if (objective_value < objective_threshold) {
  cat("PASS: Objective function value is very close to zero. Good fit.\n")
} else {
  cat("WARNING: Objective function value is higher than expected. Check your fit.\n")
}

# --- Check Maximum Gradient ---
if (max_gradient < gradient_threshold) {
  cat("PASS: Maximum gradient is small. Optimization properly converged.\n")
} else {
  cat("WARNING: Maximum gradient is too large! Possible convergence problems.\n")
}

