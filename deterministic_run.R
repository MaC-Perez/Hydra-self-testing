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
hydraDataList$Ndietprop_obs<-(Ndietprop_obs=112)
hydraDataList[["observedSurvDiet"]][["inpN"]]<-(inpN=25)
hydraDataList[["observedBiomass"]][["cv"]]<-(cv=0.2)
hydraDataList[["observedCatchSize"]][["inpN"]]<-(inpN=25)
hydraDataList[["observedSurvSize"]][["inpN"]]<-(inpN=25)

repfile <- "OM_scenarios/OM/hydra_sim.rep"
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
  
  obs_survey <- hydraDataList$observedSurvSize  %>% filter(survey == 1)  %>% tibble()
  obs_survey <- obs_survey %>% pivot_longer(cols=6:ncol(.), names_to = "lenbin") %>% #filter(value != -999)%>%
    
    mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")),
           label = rep("survey",nrow(.)))
  obs_survey$value[which(obs_survey$value == -999)] = 0.00001
  
  pred_surv<-output$pred_survey_size
  obs_survey$pred_surv<-pred_surv
  
  obs_survey <-obs_survey %>% 
    mutate(value=pred_surv) %>%
    select(-pred_surv)
  
# Convertir de vuelta a formato ancho
  hydraDataList$observedSurvSize <- obs_survey %>%
    pivot_wider(names_from = lenbin, values_from = value,
                names_prefix = "sizebin") %>%
    arrange(survey, species, year) %>%
    select(-label)
  
  #### REPLACE DIET COMPOSITION DATA from mu om estimates####
  
  obsdiet_comp <- hydraDataList$observedSurvDiet %>% tibble()
  obsdiet_comp<-obsdiet_comp %>% pivot_longer(cols=6:ncol(.), names_to = "prey") %>%
    mutate(#lenbin = as.integer(str_remove(lenbin, "V")),
      #species = hydraDataList$speciesList[species],
      label = rep("diet",nrow(.)))
  obsdiet_comp$value[which(obsdiet_comp$value == -999)] = 0.000001
  obsdiet_comp <- select(obsdiet_comp, -label)
  
  pred_diet<-output$pred_dietprop
  if (length(pred_diet)!=nrow(obsdiet_comp)) obsdiet_comp <- obsdiet_comp %>% filter(value != 0)
  nll_diet<-output$nll_dietprop
  obsdiet_comp$value<-pred_diet
  
  hydraDataList$observedSurvDiet<-obsdiet_comp
  hydraDataList$observedSurvDiet<-hydraDataList$observedSurvDiet %>% pivot_wider(names_from = "prey")
  
}

# Save the deterministic dataset
saveRDS(hydraDataList, "deterministic_sim_data.rds")

#### WRITE tsDat FUNCTION ####
source("R/write_tsDatFile.R")
source("R/read.report.R")

hydraDataList2 <- readRDS("deterministic_sim_data.rds")

listOfParameters<-list()
listOfParameters$outDir<-paste0(getwd(),"/","sims","/","deterministic","/")
listOfParameters$outputFilename<-"hydra_sim_deterministic"
listOfParameters$fillLength <- 2000

write_tsDatFile(hydraDataList2,listOfParameters)

######### deterministic estimability fits

library(dplyr)
library(ggplot2)

#Read your deterministic outputs

# Replace with your deterministic .rep file path
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
ggsave(filename = "Deterministic_Catch.jpeg",
       plot = catch_plot,
       width = 10, height = 8, units = "in", dpi = 300)

# --- Save Survey Biomass Plot ---
ggsave(filename = "Deterministic_Biomass.jpeg",
       plot = survey_plot,
       width = 10, height = 8, units = "in", dpi = 300)


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

