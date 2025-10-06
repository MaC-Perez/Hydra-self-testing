source("R/read.report.R")
source("R/gettables.R")
library(ggplot2)
library(dplyr)
library(tidyverse)
library(hydradata)
library(scales)

# READ SARAHS DATA SETS
hydraDataList <- readRDS("Sarah_files/hydra_sim_GBself_5bin.rds")

# READ FITS FROM SARAHS DATA SETS
repfile <- "OM_scenarios/OM_case1/hydra_sim.rep"
output <- reptoRlist(repfile)

obs_surveyB <- hydraDataList$observedBiomass %>% 
  filter(survey == 1)

obs_catchB <- hydraDataList$observedCatch

biorows <- dim(obs_surveyB)[1]
catchrows <- dim(obs_catchB)[1]

#create a table with estimated and observed values
indexfits <- gettables(repfile, biorows, catchrows)

stepperyr <- output$Nstepsyr
if (length(stepperyr)==0) stepperyr <- nrow(output$EstBsize)/hydraDataList$Nyrs/length(hydraDataList$speciesList)

est_bio <- output$EstBsize %>%
  rowSums() %>% 
  as.data.frame() %>% 
  rename(bio = ".") %>% 
  mutate(species = rep(hydraDataList$speciesList, each = hydraDataList$Nyrs*stepperyr),
         year  = rep(1:(hydraDataList$Nyrs*stepperyr),hydraDataList$Nspecies),
         year = (1-1/stepperyr) + year / stepperyr,  #5 time steps per year
         log_bio = ifelse(bio>0,log(bio),NA))

      source  = "Simulation")


library(purrr)
library(dplyr)

# Function to build file paths 
file_paths <- paste0("sims/OM_case1/rep/hydra_sim", 1:30, ".rep") 
# Read each repfile with reptoRlist()
sim_outputs <- lapply(file_paths, reptoRlist)
#save
saveRDS(sim_outputs, "sims/OM_case1/sim_outputs_case1.rds")
#read
sim_outputs <- readRDS("sims/OM_case1/sim_outputs_case1.rds")

# Process function (now works, since EstBsize exists)
process_fit <- function(rep_obj, sim_id) {
  rep_obj$EstBsize %>%
    rowSums() %>% 
    as.data.frame() %>% 
    rename(bio = ".") %>% 
    mutate(
      species = rep(hydraDataList$speciesList, each = hydraDataList$Nyrs * stepperyr),
      year    = rep(1:(hydraDataList$Nyrs * stepperyr), hydraDataList$Nspecies),
      year    = (1 - 1/stepperyr) + year / stepperyr,
      log_bio = ifelse(bio > 0, log(bio), NA),
      sim     = paste0("sim", sim_id),
      source  = "Simulation"
    )
}

# Apply to all 100 sims
fits_all <- purrr::map_dfr(1:30, function(i) {
  process_fit(sim_outputs[[i]], i)
})

# Combine OM + Sims
est_bio_long <- est_bio %>%
  transmute(
    year,
    species,
    pred_bio = bio,
    sim = "OM",
    source = "OM"
  )

fits_all_long <- fits_all %>%
  transmute(
    year,
    species,
    pred_bio = bio,
    sim,
    source
  )

bio_all <- bind_rows(est_bio_long, fits_all_long)


sim_OM_survey <- ggplot(bio_all, aes(x = year, y = pred_bio)) +
  geom_line(data = filter(bio_all, source == "Simulation"),
            aes(group = interaction(sim, species)),
            color = "lightblue", alpha = 0.2) +
  geom_line(data = filter(bio_all, source == "OM"),
            color = "black", size = 1) +
  facet_wrap(~species, scales = "free") +
  theme_minimal() +
  labs(x = "Year",
       y = "Biomass (t)",
       title = "Time series of OM biomass (black) and simulation fits (blue)") +
  scale_y_continuous(labels = label_number(accuracy = 1))

sim_OM_survey

#ggsave("plots/sim_OM_survey.png", plot = sim_OM_survey, width = 10, height = 6, dpi = 300)

# ============================================
# 1. Separate OM and Simulation values
# ============================================
om_df <- bio_all %>%
  filter(source == "OM") %>%
  select(year, species, true_bio = pred_bio)

sim_df <- bio_all %>%
  filter(source == "Simulation") %>%
  select(year, species, sim, est_bio = pred_bio)

# Join OM truth to each simulation estimate
eval_df_bio <- sim_df %>%
  left_join(om_df, by = c("year", "species"))

# ============================================
# 2. Compute metrics
# ============================================

# Bias
bias_df <- eval_df_bio %>%
  group_by(year, species) %>%
  summarise(bias = mean(est_bio - true_bio, na.rm = TRUE), .groups = "drop")

# Relative Bias
relbias_df <- eval_df_bio %>%
  group_by(year, species) %>%
  summarise(rel_bias = mean((est_bio - true_bio) / true_bio, na.rm = TRUE), .groups = "drop")

# Precision (SD across sims)
precision_df <- eval_df_bio %>%
  group_by(year, species) %>%
  summarise(precision = sd(est_bio, na.rm = TRUE), .groups = "drop")

# RMSE
rmse_df <- eval_df_bio %>%
  group_by(year, species) %>%
  summarise(rmse = sqrt(mean((est_bio - true_bio)^2, na.rm = TRUE)), .groups = "drop")

# ============================================
# 3. Combine into one summary table
# ============================================
summary_df_bio <- bias_df %>%
  left_join(relbias_df,  by = c("year", "species")) %>%
  left_join(precision_df, by = c("year", "species")) %>%
  left_join(rmse_df,     by = c("year", "species"))

# ============================================
# 4. Example plots
# ============================================

# Bias over time
bias_plot <- ggplot(bias_df, aes(x = year, y = bias, color = species)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Bias in Estimated Biomass", y = "Bias (est - true)")

# Histogram of relative errors
ggplot(summary_df_bio, aes(x = rel_bias)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Relative Errors by Species (Biomass)",
       x = "Relative Error", y = "Count") +
  facet_wrap(~species, scales = "free_y")

# Density plot
dp_bio <- ggplot(summary_df_bio, aes(x = rel_bias)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Density of Relative Errors by Species (Biomass)",
       x = "Relative Error", y = "Density") +
  facet_wrap(~species, scales = "free_y")

# Boxplot by year
re_year_bio <- ggplot(summary_df_bio, aes(x = factor(year), y = rel_bias)) +
  geom_boxplot() +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Relative Errors in Biomass by Year",
       x = "Year", y = "Relative Error")

# Save the "Relative Errors in Biomass by Year" plot
# ggsave("plots/re_year_bio.png", plot = re_year_bio, width = 10, height = 6, dpi = 300)

# Q-Q plot
qqplot_bio <- ggplot(summary_df_bio, aes(sample = rel_bias)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(title = "Q-Q Plot of Relative Errors by Species (Biomass)",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  facet_wrap(~species, scales = "free")

#######
#CATCH

# 1. OM CATCH (truth from indexfits)
# ============================================
indexfits_catch <- indexfits[[2]] %>%
  mutate(species = hydraDataList$speciesList[species]) %>%
  select(fishery, year, species, pred_catch) %>%
  mutate(sim = "OM",
         source = "OM")

# ============================================
# 2. Simulation CATCH (from 100 fits)
# ============================================

catch_colnames <- c("fishery", "area", "year", "species",
                    "catch", "cv", "pred_catch", "residual", "nll")

catch_fits_all <- map_dfr(1:30, function(i) {
  df <- as.data.frame(sim_outputs[[i]][["catch"]])
  colnames(df) <- catch_colnames
  
  df %>%
    mutate(sim = paste0("sim", i)) %>%       # tag simulation number
    relocate(sim)                            # put sim column first
})

catch_fits_all <- catch_fits_all %>%
  mutate(
    species = hydraDataList$speciesList[species],  # replace numeric ID with species name
    source  = "Simulation"                         # add new column
  )

# ============================================
# 3. Combine OM + Sims
# ============================================
catch_fits_all_clean <- catch_fits_all %>%
  select(fishery, year, species, pred_catch, sim, source)

catch_all_df <- bind_rows(indexfits_catch, catch_fits_all_clean)

# ==========================
#   PLOT OM vs Simulation
# ==========================
sim_OM_catch <- ggplot(catch_all_df, aes(x = year, y = pred_catch)) +
  geom_line(data = filter(catch_all_df, source == "Simulation"),
            aes(group = interaction(sim, species)),
            color = "red", alpha = 0.2) +
  geom_line(data = filter(catch_all_df, source == "OM"),
            color = "black", linewidth = 1) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(x = "Year",
       y = "Catch (t)",
       title = "Time series of OM catch (black) and simulation fits (red)") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1))

sim_OM_catch
# ggsave("plots/sim_OM_catch.png", plot = sim_OM_catch, width = 10, height = 6, dpi = 300)

# ==========================
#   METRICS
# ==========================

# Separar OM y Simulaciones
om_df_catch <- catch_all_df %>%
  filter(source == "OM") %>%
  select(year, species, true_catch = pred_catch)

sim_df_catch <- catch_all_df %>%
  filter(source == "Simulation") %>%
  select(year, species, sim, est_catch = pred_catch)

# Combinar (cada sim con el valor verdadero correspondiente)
eval_df_catch <- sim_df_catch %>%
  left_join(om_df_catch, by = c("year", "species"))

# Bias
bias_df_catch <- eval_df_catch %>%
  group_by(year, species) %>%
  summarise(bias = mean(est_catch - true_catch, na.rm = TRUE), .groups = "drop")

# Relative Bias
relbias_df_catch <- eval_df_catch %>%
  group_by(year, species) %>%
  summarise(rel_bias = mean((est_catch - true_catch) / true_catch, na.rm = TRUE), .groups = "drop")

# Precision (SD entre sims)
precision_df_catch <- eval_df_catch %>%
  group_by(year, species) %>%
  summarise(precision = sd(est_catch, na.rm = TRUE), .groups = "drop")

# RMSE
rmse_df_catch <- eval_df_catch %>%
  group_by(year, species) %>%
  summarise(rmse = sqrt(mean((est_catch - true_catch)^2, na.rm = TRUE)), .groups = "drop")

# Combinar métricas
summary_df_catch <- bias_df_catch %>%
  left_join(relbias_df_catch, by = c("year", "species")) %>%
  left_join(precision_df_catch, by = c("year", "species")) %>%
  left_join(rmse_df_catch, by = c("year", "species"))

# ==========================
#   PLOTS DIAGNÓSTICOS
# ==========================

# Bias over time
bias_plot_catch <- ggplot(bias_df_catch, aes(x = year, y = bias, color = species)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Bias in Estimated Catch", y = "Bias (est - true)")

# Histogram of relative errors
ggplot(summary_df_catch, aes(x = rel_bias)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Relative Errors by Species (Catch)",
       x = "Relative Error", y = "Count") +
  facet_wrap(~species, scales = "free_y")

# Density plot
dp_catch <- ggplot(summary_df_catch, aes(x = rel_bias)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Density of Relative Errors by Species (Catch)",
       x = "Relative Error", y = "Density") +
  facet_wrap(~species, scales = "free_y")

# Boxplot of relative errors by year
re_year_catch <- ggplot(summary_df_catch, aes(x = factor(year), y = rel_bias)) +
  geom_boxplot() +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Relative Errors in Catch by Year",
       x = "Year", y = "Relative Error")

# Q-Q plot
qqplot_catch <- ggplot(summary_df_catch, aes(sample = rel_bias)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(title = "Q-Q Plot of Relative Errors by Species (Catch)",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  facet_wrap(~species, scales = "free")

#ggsave("plots/re_year_catch.png", plot = re_year_catch, width = 10, height = 6, dpi = 300)
