source("R/read.report.R")
source("R/gettables.R")
library(tidyverse)
#library(magrittr)
#library(FSA)
library(readxl)
library(hydradata)
# rm(list = ls())

#surv_dietcomps<-read_xlsx("2surv_dietcomps.xlsx", col_names = TRUE)

### Read observed and estimated values, Hydra data list from Sarahs 4 species scenario
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
output<-reptoRlist(repfile)

#repfile<-read.table("inputs/initial_run/pmse_predvals.out", header = FALSE, skip=2, nrow=336)
#colnames(repfile)<-c("survey", "year", "spp", "area", "pred_survey")

#### READ CATCH AND SURVEY OBSERVED BIOMASS ####
# 2 surveys use this line 
#obs_surveyB <- hydraDataList$observedBiomass
# 1 survey use this line
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

for (isim in 1:100) {  
  
  # replace index with simulated data
  
  ##### SIMULATE CATCH DATA ######
  hydraDataList$observedCatch <- obs_catchB %>%
    mutate(catch = rnorm(nrow(.), log(indexfits[[2]]$pred_catch),indexfits[[2]]$cv)) # sd=0.000001
  
  sim_data[[isim]] <- hydraDataList
  
  
  ##### SIMULATE SURVEY BIOMASS DATA ######
  # if you have one survey
  hydraDataList$observedBiomass <- obs_surveyB %>%
    mutate(biomass = rnorm(nrow(.), log(indexfits[[1]]$pred_bio),indexfits[[1]]$cv)) # sd=0.000001
  # store simulated object
  sim_data[[isim]] <- hydraDataList
  
  
  #### SIMULATE CATCH SIZE COMPOSITION DATA ####
  
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
  
  
  # Agrupar por combinaciones únicas de fishery, species y year
  group_keys <- obs_catch %>%
    select(fishery, species, year) %>%
    distinct()
  
  temporal1 = numeric()
  # Iterar por grupo
  for (i in 1:nrow(group_keys)) {
    this_group <- group_keys[i, ]
    
    # Subgrupo para esa combinación
    subgroup <- obs_catch %>%
      filter(
        fishery == this_group$fishery,
        species == this_group$species,
        year == this_group$year
      ) %>%
      arrange(lenbin)
    
    probs <- subgroup$value
    probs <- probs / sum(probs)  # asegurar que sumen 1
    sim_counts <- rmultinom(1, size = 100, prob = probs)
    simulated_props <- as.numeric(sim_counts) / 100
    
    temporal1 <- c(temporal1, simulated_props)
  }
  
  # Reemplazar valores simulados
  obs_catch$value <- temporal1
  
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
  
  # Agrupar combinaciones únicas
  group_keys <- obs_survey %>%
    select(survey, species, year) %>%
    distinct()
  
  temporal1 = numeric()
 
  # Simulación por grupo
  for (i in 1:nrow(group_keys)) {
    this_group <- group_keys[i, ]
    
    subgroup <- obs_survey %>%
      filter(
        survey == this_group$survey,
        species == this_group$species,
        year == this_group$year
      ) %>%
      arrange(lenbin)
     
    probs <- subgroup$value
    probs <- probs / sum(probs)
    sim_counts <- rmultinom(1, size = 100, prob = probs)
    simulated_props <- as.numeric(sim_counts) / 100
    
    temporal1 <- c(temporal1, simulated_props)
  }
  
  # Reasignar valores simulados
  obs_survey$value <- temporal1
  
  # Convertir de vuelta a formato ancho
  hydraDataList$observedSurvSize <- obs_survey %>%
    pivot_wider(names_from = lenbin, values_from = value,
    names_prefix = "sizebin") %>%
    arrange(survey, species, year) %>%
                select(-label)
  
  #to check if I am getting the correct values
  #write.csv(surv_size, file = "surv_size.csv", row.names = T)
  
  
  #### SIMULATE DIET COMPOSITION DATA ####
  
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
  
  hydraDataList$observedSurvDiet
  
  temporal1 = numeric()
  # Lista única de combinaciones de grupos
  group_keys <- obsdiet_comp %>%
    select(survey, species, year, sizebin) %>%
    distinct()
  
    # Iterar por combinación única
  for (i in 1:nrow(group_keys)) {
    this_group <- group_keys[i, ]
    
    # Filtrar solo las filas de esa combinación
    subgroup <- obsdiet_comp %>%
      filter(
        survey == this_group$survey,
        species == this_group$species,
        year == this_group$year,
        sizebin == this_group$sizebin
      )
    
    # Obtener proporciones de presa y simular
    probs <- subgroup$value
    probs <- probs / sum(probs)  # normalizar
    simulated_counts <- rmultinom(1, size = 100, prob = probs)
    simulated_props <- as.numeric(simulated_counts) / 100
    
    temporal1 <- c(temporal1, simulated_props)
  }
  
  # replace data with simulated data
  
  obsdiet_comp$value = temporal1
  hydraDataList$observedSurvDiet<-obsdiet_comp
  hydraDataList$observedSurvDiet<-hydraDataList$observedSurvDiet %>% pivot_wider(names_from = "prey")
  
}

# change simulated catch and survey biomass data from log scale to the original scale 
for (isim in 1:100) {  
  
  sim_data[[isim]][["observedBiomass"]][["biomass"]]<-exp(sim_data[[isim]][["observedBiomass"]][["biomass"]])
  sim_data[[isim]][["observedCatch"]][["catch"]]<-exp(sim_data[[isim]][["observedCatch"]][["catch"]])
  
  }

#write.csv(indexfits[[1]], file = "original.csv", row.names = T)

# save the simulated data object 
#write_rds(sim_data, "sim_data.rds")
#write_rds(sim_data, "sim_data_1survey.rds")

#### WRITE tsDat FUNCTION ####
source("R/write_tsDatFile.R")
source("R/read.report.R")

#read original observations (hydraDataList) and simulated data sets (hydraDataList2)
#hydraDataList <- readRDS("inputs/hydra_sim_GBself_5bin.rds")
#hydraDataList2 <- readRDS("sim_data.rds")
hydraDataList2 <- readRDS("sim_data_1survey.rds")


listOfParameters<-list()
listOfParameters$outDir<-paste0(getwd(),"/","sims","/")
listOfParameters$outputFilename<-"hydra_sim"
listOfParameters$fillLength <- 2000

for (nsim in 1:100){ 
  write_tsDatFile(hydraDataList2[[nsim]],listOfParameters)
}


#### RUN THE MODEL WITH 100 SIMS ####
library(here)
dir<-here()
#dir<-paste0(dir,"/","sims","/","initial")
dir<-paste0(dir,"/","sims","/","OM")

setwd(dir)

#multiple calls to 'system()' given different folders/filenames.

#isim<-1
#system("cp sims/hydra_sim1-ts.dat sims/hydra_sim_GBself_5bin-ts.dat")
#file.copy(from="sims/hydra_sim1-ts", to="sims/hydra_sim_GBself_5bin-ts")
nsim<-1
for (nsim in 1:100) ### DO NOT FORGET TO CHANGE THE NAME OF THE FILE hydra_sim_GBself_5bin-ts.dat FOR hydra_GBself_5bin_simdata-ts.dat
{
  file.copy(from=paste0("hydra_sim",nsim,"-ts.dat"), to= "hydra_GBself_5bin_simdata-ts.dat", overwrite = TRUE)
  system("./hydra_sim -ind hydra_sim_GBself_5bin.dat -ainp hydra_sim_GBself_5bin.pin")
  file.copy(from = "hydra_sim.rep", to = paste0("rep/hydra_sim",nsim,".rep"))
  file.copy(from = "hydra_sim.par", to = paste0("par/hydra_sim",nsim,".par"))
  file.copy(from = "HYDRA_~1.cor", to = paste0("cor/hydra_sim",nsim,".cor"))
  file.copy(from = "HYDRA_~1.std", to = paste0("std/hydra_sim",nsim,".std"))
  file.copy(from = "pmse_predvals.out", to = paste0("out/pmse_predvals",nsim,".out"))
}



#######
## PLOT FITTED DATA
####
source("R/read.report.R")
source("R/gettables.R")
library(tidyverse)

hydraDataList <- readRDS("Sarah_files/hydra_sim_GBself_5bin.rds")
obs_surveyB <- hydraDataList$observedBiomass %>% 
  filter(survey == 1)

obs_catchB <- hydraDataList$observedCatch

biorows <- dim(obs_surveyB)[1]
catchrows <- dim(obs_catchB)[1]

#READ FITS FROM OM
OMrepfile <- "OM_scenarios/OM/hydra_sim.rep"
OMoutput<-reptoRlist(OMrepfile)

#pred_catch and pred_bio are my estimated values from the OM 
indexfits <- gettables(OMrepfile, biorows, catchrows)
OM_catch <- indexfits[[2]]
OM_survey <- indexfits[[1]]

#READ FITS FROM RUNS WITH 100 DATA SETS
rep_files_sim <- paste0("sims/OM/rep/hydra_sim", 1:100, ".rep")
#repfile <- "inputs/initial_run/hydra_sim.rep"

# Read all rep files into a named list
rep_simoutputs <- lapply(rep_files_sim, function(file) {
  cat("Reading:", file, "\n")
  reptoRlist(file)
})


#############
### CATCH
#############

# Estimated OM catch vs estimate EM 100 fits (from rep files)

#Extract EM catch for all runs
all_catch_fits <- purrr::map_dfr(seq_along(rep_simoutputs), function(i) {
  cat("Extracting from:", rep_files_sim[i], "\n")
  
  # Use gettables() on each .rep file path
  tablist <- gettables(rep_files_sim[i], biorows, catchrows)
  
  catchdf <- tablist[[2]] %>%
    mutate(
      obs = catch + 1e-07,
      pred_EM = pred_catch + 1e-07,
      log_obs = log(obs),
      log_pred_EM = log(pred_EM),
      log_lo = log_obs - 1.96 * cv,
      log_hi = log_obs + 1.96 * cv,
      species = hydraDataList$speciesList[species],
      isim = i  # simulation index
    )
  
  return(catchdf)
})

# extract OM catch estimates
OM_catch_plot <- OM_catch %>%
       mutate(
             obs = catch + 1e-07,
             pred_OM = pred_catch + 1e-07,
             log_pred_OM = log(pred_OM),
             species = hydraDataList$speciesList[species],
             isim = 0  # identifier for OM
         )

OM_EMest_catch<-ggplot() +
       # EM simulation lines
       geom_line(data = all_catch_fits,
                                 aes(x = year, y = log_pred_EM, group = isim, color = "EM estimates"),
                                 alpha = 0.3, linewidth = 0.5) +
       # OM line
       geom_line(data = OM_catch_plot,
                                 aes(x = year, y = log_pred_OM, color = "OM estimate"),
                                 linewidth = 1.2) +
       facet_wrap(~species, scales = "free_y") +
       labs(x = "Year", y = "Log(Catch)", color = "Model") +
       theme_minimal()

print(OM_EMest_catch)

#ggsave("OM_EMest_catch.jpeg",
 #            plot = OM_EMest_catch,
  #          width = 10, height = 8, units = "in", dpi = 300)


#############
### SURVEY
#############

# Estimated OM survey vs estimate EM 100 fits (from rep files)

#Extract EM survey biomass for all runs
all_survey_fits <- purrr::map_dfr(seq_along(rep_simoutputs), function(i) {
  cat("Extracting survey biomass from:", rep_files_sim[i], "\n")
  
  tablist <- gettables(rep_files_sim[i], biorows, catchrows)
  
  survey_df <- tablist[[1]] %>%
    mutate(
      pred_EM = pred_bio + 1e-07,
      log_pred_EM = log(pred_EM),
      species = hydraDataList$speciesList[species],
      isim = i
    )
  
  return(survey_df)
})

## extract OM survey biomass estimates
OM_survey_plot <- OM_survey %>%
  mutate(
    pred_OM = pred_bio + 1e-07,
    log_pred_OM = log(pred_OM),
    species = hydraDataList$speciesList[species],
    isim = 0  # identifier for OM
  )

OM_EMest_survey<-ggplot() +
  # EM simulation lines
  geom_line(data = all_survey_fits,
            aes(x = year, y = log_pred_EM, group = isim, color = "EM estimates"),
            alpha = 0.3, linewidth = 0.5) +
  # OM line
  geom_line(data = OM_survey_plot,
            aes(x = year, y = log_pred_OM, color = "OM estimate"),
            linewidth = 1.2) +
  facet_wrap(~species, scales = "free_y") +
  labs(x = "Year", y = "Log(Survey Biomass)", color = "Model") +
  theme_minimal()

print(OM_EMest_survey)

#ggsave("OM_EMest_survey.jpeg",
 #                  plot = OM_EMest_survey,
  #               width = 10, height = 8, units = "in", dpi = 300)
       
       
#boxplot catch prediction
# median of the 100 preds and the 25 and 75 percentiles (IQR), whiskers +- 1.5x IQR (IQR being 50%) and outliers shown as indivisul points
# If for a species in a given year:
#Q1 = 7.0
#Q3 = 7.8
#IQR = 0.8
#Then:
#Whiskers would extend up to:
#Lower: Q1 − 1.5 × 0.8 = 5.8
#Upper: Q3 + 1.5 × 0.8 = 9.0
#Any points outside [5.8, 9.0] would appear as outliers.

catch_boxplot<- ggplot(all_catch_fits, aes(x = factor(year), y = log_pred_EM)) +
  geom_boxplot(fill = "lightblue", outlier.size = 0.5, outlier.alpha = 0.3) +
  facet_wrap(~species, scales = "free_y") +
  labs(x = "Year", y = "Log(Catch)", title = "Catch Predictions: EM Boxplots by Year") +
  theme_minimal()

print(catch_boxplot)

#ggsave("catch_boxplot.jpeg",
#                  plot = catch_boxplot,
#               width = 10, height = 8, units = "in", dpi = 300)

#boxplot survey prediction
survey_boxplot<-ggplot(all_survey_fits, aes(x = factor(year), y = log_pred_EM)) +
  geom_boxplot(fill = "lightgreen", outlier.size = 0.5, outlier.alpha = 0.3) +
  facet_wrap(~species, scales = "free_y") +
  labs(x = "Year", y = "Log(Survey Biomass)", title = "Survey Biomass: EM Boxplots by Year") +
  theme_minimal()

print(survey_boxplot)

#ggsave("survey_boxplot.jpeg",
#       plot = survey_boxplot,
#       width = 10, height = 8, units = "in", dpi = 300)


# estimated catch from 100 EM vs. observed true catch from the operating model (OM)

est_summary_catch <- all_catch_fits %>%
  group_by(year, species, fishery) %>%
  summarise(
    mean_log_est = mean(log_pred_EM, na.rm = TRUE),
    lo = quantile(log_pred_EM, 0.025, na.rm = TRUE),
    hi = quantile(log_pred_EM, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

true_summary_catch <- OM_catch %>%
  mutate(
    log_true = log(catch + 1e-07),
    species = hydraDataList$speciesList[species]  # optional: convert to names
  ) %>%
  group_by(year, species, fishery) %>%
  summarise(mean_log_true = mean(log_true, na.rm = TRUE), .groups = "drop")


plot_data_catch <- left_join(est_summary_catch, true_summary_catch,
                             by = c("year", "species", "fishery"))

fleet3_plot <-ggplot(plot_data_catch, aes(x = year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = mean_log_est), color = "blue", size = 1) +
  geom_line(aes(y = mean_log_true), color = "black", linetype = "dashed", size = 1) +
  facet_wrap(~species, scales = "free_y") + 
  labs(
    title = "Mean Estimated Catch vs. Operating Model Catch",
    x = "Year",
    y = "log(Catch)",
    caption = "Blue: Estimated mean ± 95% CI | Black dashed: OM catch"
  ) +
  theme_minimal()

print(fleet3_plot)

#ggsave("fleet3_fitcatch_plot.jpeg",
#      plot = fleet3_plot,
#     width = 10, height = 8, units = "in", dpi = 300)



#############
### SURVEY BIOMASS
#############
# estimated catch from 100 EM vs. observed true catch from the operating model (OM)

# Estimated survey biomass summary (from 100 rep files)
est_summary_survey <- all_survey_fits %>%
  group_by(year, species, survey) %>%
  summarise(
    mean_log_est = mean(log_pred_EM, na.rm = TRUE),
    lo = quantile(log_pred_EM, 0.025, na.rm = TRUE),
    hi = quantile(log_pred_EM, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# "True" survey biomass from OM .rep file

true_summary_survey <- OM_survey %>%
  mutate(
    log_true = log(biomass + 1e-07),
    species = hydraDataList$speciesList[species]
  ) %>%
  group_by(year, species, survey) %>%
  summarise(
    mean_log_true = mean(log_true, na.rm = TRUE),
    .groups = "drop"
  )

# Join estimated and true summaries
plot_data_survey <- left_join(est_summary_survey, true_summary_survey,
                              by = c("year", "species", "survey"))

# Plot
survey_plot <- ggplot(plot_data_survey, aes(x = year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "green", alpha = 0.3) +
  geom_line(aes(y = mean_log_est), color = "darkgreen", size = 1) +
  geom_line(aes(y = mean_log_true), color = "black", linetype = "dashed", size = 1) +
  facet_wrap(~species, scales = "free_y") +
  labs(
    title = "Mean Estimated Survey Biomass vs. Operating Model Survey Biomass",
    x = "Year",
    y = "log(Biomass)",
    caption = "Green: Estimated mean ± 95% CI | Black dashed: OM survey biomass"
  ) +
  theme_minimal()

print(survey_plot)

#ggsave("survey_fit_plot.jpeg", plot = survey_plot, width = 10, height = 8, units = "in", dpi = 300)

#################
#### PERFORMANCE METRICS
##################

# Per species/fleet/year or just overall
#RMSE Error between the model estimate and OM value. Shows the average size of errors — larger values mean greater deviation.

OM_vs_EM <- all_catch_fits %>%
  left_join(OM_catch_plot %>% 
              select(year, species, fishery, pred_OM, log_pred_OM),  # only bring what you need
            by = c("year", "species", "fishery")) %>%
  mutate(
    error = log_pred_EM - log_pred_OM,
    abs_error = abs(error),
    squared_error = error^2
  )

OM_vs_EM <- OM_vs_EM %>%
  select(fishery, area, year, species, catch, cv, residual, nll, pred_EM, log_pred_EM, pred_OM, log_pred_OM, isim, error,  abs_error, squared_error)

#Summary by simulation
error_summary <- OM_vs_EM %>%
  group_by(isim) %>%
  summarise(
    RMSE = sqrt(mean(squared_error, na.rm = TRUE)),
    MAE  = mean(abs_error, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(error_summary, aes(x = RMSE)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of RMSE across EM simulations")


# RMSE over time by species
error_by_year_species <- OM_vs_EM %>%
  group_by(year, species) %>%
  summarise(
    RMSE = sqrt(mean(squared_error, na.rm = TRUE)),
    MAE = mean(abs_error, na.rm = TRUE),
    .groups = "drop"
  )

# Plot RMSE over time
ggplot(error_by_year_species, aes(x = year, y = RMSE)) +
  geom_line(color = "darkred") +
  facet_wrap(~species, scales = "free_y") +
  labs(title = "RMSE Over Time by Species", y = "RMSE", x = "Year") +
  theme_minimal()

#RMSE per species across simulations
error_by_sim_species <- OM_vs_EM %>%
  group_by(isim, species) %>%
  summarise(
    RMSE = sqrt(mean(squared_error, na.rm = TRUE)),
    .groups = "drop"
  )

ggplot(error_by_sim_species, aes(x = species, y = RMSE)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "RMSE per Species Across Simulations", x = "Species", y = "RMSE") +
  theme_minimal()


###################
#### SURVEY
###################


OM_vs_EM_survey <- all_survey_fits %>%
    left_join(OM_survey_plot %>%
              select(year, species, survey, pred_OM, log_pred_OM),
            by = c("year", "species", "survey")) %>%
  mutate(
    error = log_pred_EM - log_pred_OM,
    abs_error = abs(error),
    squared_error = error^2
  )

OM_vs_EM_survey <- OM_vs_EM_survey %>%
  select(survey, year, species, biomass, cv, pred_bio, pred_EM, log_pred_EM, pred_OM, log_pred_OM, isim, error, abs_error, squared_error)

#Summary by simulation
error_summary_survey <- OM_vs_EM_survey %>%
  group_by(isim) %>%
  summarise(
    RMSE = sqrt(mean(squared_error, na.rm = TRUE)),
    MAE = mean(abs_error, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(error_summary_survey, aes(x = RMSE)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of RMSE for Survey Biomass (EM)", x = "RMSE")

#RMSE over time by species
error_by_year_species_survey <- OM_vs_EM_survey %>%
  group_by(year, species) %>%
  summarise(
    RMSE = sqrt(mean(squared_error, na.rm = TRUE)),
    MAE = mean(abs_error, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(error_by_year_species_survey, aes(x = year, y = RMSE)) +
  geom_line(color = "darkgreen") +
  facet_wrap(~species, scales = "free_y") +
  labs(title = "RMSE Over Time by Species (Survey Biomass)", y = "RMSE", x = "Year") +
  theme_minimal()

#RMSE per species across simulations
error_by_sim_species_survey <- OM_vs_EM_survey %>%
  group_by(isim, species) %>%
  summarise(
    RMSE = sqrt(mean(squared_error, na.rm = TRUE)),
    .groups = "drop"
  )

ggplot(error_by_sim_species_survey, aes(x = species, y = RMSE)) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "RMSE per Species Across Simulations (Survey Biomass)",
       x = "Species", y = "RMSE") +
  theme_minimal()




# The mean difference between estimated and true values. Positive = overestimation; negative = underestimation.

mean_diff_catch <- OM_vs_EM %>%
  group_by(species) %>%
  summarise(
    mean_diff = mean(log_pred_EM - log_pred_OM, na.rm = TRUE),  # log-scale difference
    .groups = "drop"
  )

print(mean_diff_catch)

#Positive values = the EM tends to overestimate, Negative values = the EM tends to underestimate

ggplot(mean_diff_catch, aes(x = species, y = mean_diff)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Mean Log Difference (Estimated - True) for Catch",
    x = "Species",
    y = "Mean Log Difference"
  ) +
  theme_minimal()

ggplot(OM_vs_EM, aes(x = error)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~species, scales = "free_y") +
  labs(
    title = "Residual Distribution: Estimated - OM",
    x = "Residual",
    y = "Frequency"
  ) +
  theme_minimal()


mean_diff_survey <- OM_vs_EM_survey %>%
  group_by(species) %>%
  summarise(
    mean_diff = mean(log_pred_EM - log_pred_OM, na.rm = TRUE),
    .groups = "drop"
  )

print(mean_diff_survey)

ggplot(mean_diff_survey, aes(x = species, y = mean_diff)) +
  geom_col(fill = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Mean Log Difference (Estimated - True) for Survey Biomass",
    x = "Species",
    y = "Mean Log Difference"
  ) +
  theme_minimal()


ggplot(OM_vs_EM_survey, aes(x = error)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
  facet_wrap(~species, scales = "free_y") +
  labs(
    title = "Residual Distribution (Survey): Estimated - OM",
    x = "Residual",
    y = "Frequency"
  ) +
  theme_minimal()



































#### SIMULATED DATA PLOTS ####
#browseVignettes("hydradata")

source("R/read.report.R")
source("R/gettables.R")

library(ggforce)
library(tidyverse)
hydraDataList <- readRDS("Sarah_files/hydra_sim_GBself_5bin.rds")
hydraDataList2 <- readRDS("sim_data_1survey.rds")

#### PLOT SIMULATIONS CATCH ####
# Combine and label simulation data
sim_obs_catch<-purrr::map_dfr(hydraDataList2,"observedCatch",.id = "isim") %>%
   mutate(species = hydraDataList$speciesList[species])

#Summarize across simulations (mean and CI)
sim_summary_fleet1 <- sim_obs_catch %>%
  filter(fishery == 1) %>%
  group_by(year, species) %>%
  summarise(
    mean_logcatch = mean(log(catch), na.rm = TRUE),
    lower_logcatch = quantile(log(catch), 0.025, na.rm = TRUE),
    upper_logcatch = quantile(log(catch), 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "Simulated Mean")

# Get the observed catch
obs_catch1 <- hydraDataList$observedCatch %>%
  filter(fishery == 1) %>%
  mutate(
    species = hydraDataList$speciesList[species],
    year = year,
    mean_logcatch = log(catch),
    source = "Observed"
  ) %>%
  select(year, species, mean_logcatch, source)

# Combine summary + observed into one
plot_data_catch <- bind_rows(
  sim_summary_fleet1 %>% select(year, species, mean_logcatch, source),
  obs_catch1
)

# plot sim data 
fleet1plot <- ggplot() +
  # CI ribbon (for simulated only)
  geom_ribbon(data = sim_summary_fleet1,
              aes(x = year, ymin = lower_logcatch, ymax = upper_logcatch),
              fill = "lightblue", alpha = 0.3) +
  # Lines for mean and observed
  geom_line(data = plot_data_catch,
            aes(x = year, y = mean_logcatch, color = source, linetype = source),
            linewidth = 1) +
  # NEW: Observed points
  geom_point(data = plot_data_catch %>% filter(source == "Observed"),
             aes(x = year, y = mean_logcatch),
             color = "black", shape = 16, size = 1.8) +
  facet_wrap(~species, scales = "free", dir = "v") +
  labs(x = "Year",
       y = "Catch (log t)",
       title = "Fleet 1: Simulated catch (mean, CI) vs. Observed",
       color = "Data Source",
       linetype = "Data Source") +
  scale_color_manual(values = c("Simulated Mean" = "blue", "Observed" = "black")) +
  scale_linetype_manual(values = c("Simulated Mean" = "solid", "Observed" = "dashed")) +
  theme_minimal()

print(fleet1plot)

#ggsave("fleet1_simcatch_plot.jpeg",
#       plot = fleet1plot,
#       width = 10, height = 8, units = "in", dpi = 300)

###
## fleet 2
###

# filter and summarize for Fleet 2
sim_summary_fleet2 <- sim_obs_catch %>%
  filter(fishery == 2) %>%
  group_by(year, species, fishery) %>%
  summarise(
    mean_logcatch = mean(log(catch), na.rm = TRUE),
    lower_logcatch = quantile(log(catch), 0.025, na.rm = TRUE),
    upper_logcatch = quantile(log(catch), 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "Simulated Mean")

# Get observed catch for Fleet 2
obs_catch_fleet2 <- hydraDataList$observedCatch %>%
  filter(fishery == 2) %>%
  mutate(
    species = hydraDataList$speciesList[species],
    mean_logcatch = log(catch),
    source = "Observed"
  ) %>%
  select(year, species, fishery, mean_logcatch, source)

# Combine summary and observed data
plot_data_catch2 <- bind_rows(
  sim_summary_fleet2 %>% select(year, species, fishery, mean_logcatch, source),
  obs_catch_fleet2
)

# Plot just Fleet 2
fleet2_plot <- ggplot() +
  # Confidence interval ribbon
  geom_ribbon(data = sim_summary_fleet2,
              aes(x = year, ymin = lower_logcatch, ymax = upper_logcatch, group = species),
              fill = "lightblue", alpha = 0.3) +
  # Lines for simulated mean and observed
  geom_line(data = plot_data_catch2,
            aes(x = year, y = mean_logcatch, color = source, linetype = source),
            linewidth = 1) +
  facet_wrap(~species, scales = "free_y", dir = "v") +
  labs(x = "Year",
       y = "Catch (log t)",
       title = "Fleet 2: Simulated catch (mean, CI) vs. Observed",
       color = "Data Source",
       linetype = "Data Source") +
  scale_color_manual(values = c("Simulated Mean" = "blue", "Observed" = "black")) +
  scale_linetype_manual(values = c("Simulated Mean" = "solid", "Observed" = "dashed")) +
  theme_minimal()

print(fleet2_plot)

#ggsave("fleet2_simcatch_plot.jpeg",
 #      plot = fleet2_plot,
  #     width = 10, height = 8, units = "in", dpi = 300)

# plot sim data 

##############
### biomass
##############

#### PLOT SIMULATIONS BIOMASS ####
# Combine and label simulation data
sim_surv_catch<-purrr::map_dfr(hydraDataList2,"observedBiomass",.id = "isim") %>%
  mutate(species = hydraDataList$speciesList[species])

#Summarize across simulations (mean and CI)
sim_summary_surv <- sim_surv_catch %>%
  filter(survey == 1) %>%
  group_by(year, species) %>%
  summarise(
    mean_logbio = mean(log(biomass), na.rm = TRUE),
    lower_logbio = quantile(log(biomass), 0.025, na.rm = TRUE),
    upper_logbio = quantile(log(biomass), 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "Simulated Mean")

# Get the observed catch
obs_survey <- hydraDataList$observedBiomass %>%
  filter(survey == 1) %>%
  mutate(
    species = hydraDataList$speciesList[species],
    year = year,
    mean_logbio = log(biomass),
    source = "Observed"
  ) %>%
  select(year, species, mean_logbio, source)

# Combine summary + observed into one
plot_data_survey <- bind_rows(
  sim_summary_surv %>% select(year, species, mean_logbio, source),
  obs_survey
)

# plot sim data 
surveyplot <- ggplot() +
  # CI ribbon (for simulated only)
  geom_ribbon(data = sim_summary_surv,
              aes(x = year, ymin = lower_logbio, ymax = upper_logbio),
              fill = "lightblue", alpha = 0.3) +
  # Lines for mean and observed
  geom_line(data = plot_data_survey,
            aes(x = year, y = mean_logbio, color = source, linetype = source),
            linewidth = 1) +
  # NEW: Observed points
  geom_point(data = plot_data_survey %>% filter(source == "Observed"),
             aes(x = year, y = mean_logbio),
             color = "black", shape = 16, size = 1.8) +
  facet_wrap(~species, scales = "free", dir = "v") +
  labs(x = "Year",
       y = "Biomass (log t)",
       title = "Survey 1: Simulated biomass (mean, CI) vs. Observed",
       color = "Data Source",
       linetype = "Data Source") +
  scale_color_manual(values = c("Simulated Mean" = "blue", "Observed" = "black")) +
  scale_linetype_manual(values = c("Simulated Mean" = "solid", "Observed" = "dashed")) +
  theme_minimal()

print(surveyplot)


#ggsave("survey_simbio_plot.jpeg",
#      plot = surveyplot,
#     width = 10, height = 8, units = "in", dpi = 300)











fleet2plot<-sim_obs_catch %>% filter(fishery==1)%>%
  ggplot() +
  aes(x = year, y = log(catch), col = isim) +
  geom_line() +
  facet_wrap(~species, scales = "free",dir="v") +
  #theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "Catch (t)",
       title = "Simulated catch")

print(fleet2plot)

fleet2plot<-sim_obs_catch %>% filter(fishery==2)%>%
  ggplot() +
  aes(x = year, y = (catch), col = isim) +
  geom_line() +
  facet_wrap(~species, scales = "free",dir="v") +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "Catch (t)",
       title = "Simulated catch")

print(fleet2plot)


#### PLOT SIM SURVEY BIOM ####

sim_obs_bio<-purrr::map_dfr(hydraDataList2,"observedBiomass",.id = "isim")
sim_obs_bio<- sim_obs_bio %>%
  mutate(species = hydraDataList$speciesList[species])

surv1plot<-sim_obs_bio %>% filter(survey==1)%>%
  ggplot() +
  aes(x = year, y = (biomass), col = isim) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "Biomass (t)",
       title = "Simulated biomass")

print(surv1plot)


#surv2plot<-sim_obs_bio %>% filter(survey==2)%>%
#  ggplot() +
#  aes(x = year, y = (biomass), col = isim) +
#  geom_line() +
#  facet_wrap(~species, scales = "free") +
#  theme(legend.position = "none") +
#  labs(x = "Year",
#       y = "Biomass (t)",
#       title = "Time series of estimated LN(biomass)")
#
#print(surv2plot)

#### PLOT SIM LENGHT SURVEY ####

sim_surv_lenght<-purrr::map_dfr(hydraDataList2,"observedSurvSize",.id = "isim")
sim_surv_lenght<- sim_surv_lenght %>%
  mutate(species = hydraDataList$speciesList[species])

sim_surv_lenght <- sim_surv_lenght %>% pivot_longer(cols=7:ncol(.), names_to = "lenbin") %>%
  # filter(value != -999)%>%
  mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")))

sp<-1
plot_surv <- list()
especies<-unique(sim_surv_lenght$species)
for (sp in especies) {
  
  temp_size<-sim_surv_lenght %>% filter(species == sp & survey==1) %>%
    group_by(year) %>%
    summarize(mu_ss=mean(inpN))
  
  plot_surv[[sp]] <- sim_surv_lenght %>% filter (species==sp & survey==1) %>%
    ggplot() +
    aes(x=lenbin, y = value) +
    geom_line(aes(col = isim)) +
    facet_wrap(~year, dir="v") +
    geom_text(data=temp_size, aes(x = 4.5, y = 0.5, label = mu_ss), size=3) +
    theme(legend.position = "bottom") +
    labs(col="") +
    guides(col = guide_legend(nrow = 1))
}  

plot_surv$Atlantic_cod
plot_surv$Atlantic_herring
plot_surv$Atlantic_mackerel
plot_surv$Spiny_dogfish

#sp<-1
#plot_surv <- list()
#especies<-unique(sim_surv_lenght$species)
#for (sp in especies) {
#  
#  temp_size<-sim_surv_lenght %>% filter(species == sp & survey==2) %>%
#    group_by(year) %>%
#    summarize(mu_ss=mean(inpN))
#  
#  plot_surv[[sp]] <- sim_surv_lenght %>% filter (species==sp & survey==2) %>%
#    ggplot() +
#    aes(x=lenbin, y = value) +
#    geom_line(aes(col = isim)) +
#    facet_wrap(~year, dir="v") +
#    geom_text(data=temp_size, aes(x = 4.5, y = 0.5, label = mu_ss), size=3) +
#    theme(legend.position = "bottom") +
#    labs(col="") +
#    guides(col = guide_legend(nrow = 1))
#}  

#plot_surv$Atlantic_cod
#plot_surv$Atlantic_herring
#plot_surv$Atlantic_mackerel
#plot_surv$Spiny_dogfish

#### PLOT SIM LENGHT CATCH ####

sim_catch_lenght<-purrr::map_dfr(hydraDataList2,"observedCatchSize",.id = "isim")
sim_catch_lenght<- sim_catch_lenght %>%
  mutate(species = hydraDataList$speciesList[species])

sim_catch_lenght <- sim_catch_lenght %>% pivot_longer(cols=8:ncol(.), names_to = "lenbin") %>%
  # filter(value != -999)%>%
  mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")))

sp<-1
plot_catch <- list()
especies<-unique(sim_catch_lenght$species)
for (sp in especies) {
  
  temp_size<-sim_catch_lenght %>% filter(species == sp & fishery  ==1) %>%
    group_by(year) %>%
    summarize(mu_ss=mean(inpN))
  
  plot_catch[[sp]] <- sim_catch_lenght %>% filter (species==sp & fishery  ==1) %>%
    ggplot() +
    aes(x=lenbin, y = value) +
    geom_line(aes(col = isim)) +
    facet_wrap(~year, dir="v") +
    geom_text(data=temp_size, aes(x = 4.5, y = 0.5, label = mu_ss), size=3) +
    theme(legend.position = "none") + #bottom
    labs(col="") #+
  #guides(col = guide_legend(nrow = 1))
}  

plot_catch$Atlantic_cod
plot_catch$Spiny_dogfish

sp<-1
plot_catch <- list()
especies<-unique(sim_catch_lenght$species)
for (sp in especies) {
  
  temp_size<-sim_catch_lenght %>% filter(species == sp & fishery  ==2) %>%
    group_by(year) %>%
    summarize(mu_ss=mean(inpN))
  
  plot_catch[[sp]] <- sim_catch_lenght %>% filter (species==sp & fishery  ==2) %>%
    ggplot() +
    aes(x=lenbin, y = value) +
    geom_line(aes(col = isim)) +
    facet_wrap(~year, dir="v") +
    geom_text(data=temp_size, aes(x = 4.5, y = 0.5, label = mu_ss), size=3) +
    theme(legend.position = "bottom") +
    labs(col="") +
    guides(col = guide_legend(nrow = 1))
}  

plot_catch$Atlantic_herring
plot_catch$Atlantic_mackerel


#### PLOT SIM DIET COMP ####
##### WE DONT NEED THIS #######
sim_diet<-purrr::map_dfr(hydraDataList2,"observedSurvDiet",.id = "isim")
sim_diet <- sim_diet %>% pivot_longer(cols=7:ncol(.), names_to = "prey") %>%
  rename(prop = `value`)    
sim_diet$type2<-"o" #rep(("o"),each=22810)
sim_diet$sizefit<- paste0(sim_diet$sizebin,".",sim_diet$type2)

sp<-1
nsize <- 5 #hydraDataList2$Nsizebins### problem for each sim
stringbit <- paste0(rep(1:nsize, each=2),c(".o",".e"))
limits_use <- rep("",3*nsize)
breaks_use <- rep(NA,3*nsize)
for (i in 1:nsize) {
  lo <- i*3-2
  hi <- i*3-1
  limits_use[lo:hi] <- paste0(rep(i, 2),c(".o",".e"))
  breaks_use[lo:hi] <- limits_use[lo:hi]
}

especies<-unique(sim_diet$species)
plot_diet <- NULL
for (sp in especies) {
  
  plot_diet[[sp]] <-  sim_diet %>%
    tibble() %>% 
    filter(species == sp & survey == 1) %>%
    mutate(prop = as.numeric(prop)) %>% 
    ggplot(aes(x = sizefit, y = prop, group = type2, fill = prey)) +
    geom_col(position = "fill") +
    scale_x_discrete(limits = limits_use,
                     breaks = breaks_use,
                     labels = limits_use) + 
    coord_flip() +
    facet_wrap(~as.numeric(year)) +
    theme_bw() +
    labs(x = "size & source (o=observed, e=expected)",
         fill = "prey",
         y = "proportion in diet") +
    scale_fill_brewer(type = "qual", palette = 3)
  
}
