source("R/read.report.R")
source("R/gettables.R")
library(ggplot2)
library(dplyr)
library(tidyverse)
library(hydradata)
library(scales)

nsims<-100
# READ SARAHS DATA SETS
hydraDataList <- readRDS("sims/hydra_sim_GBself_5bin_combined.rds")

# READ FITS FROM SARAHS DATA SETS
repfile <- "OM_scenarios/OM_base/hydra_sim.rep"
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



library(purrr)
library(dplyr)

# Function to build file paths 
#file_paths <- paste0("sims/low/rep/hydra_sim", 1:nsims, ".rep") 
# Read each repfile with reptoRlist()
#sim_outputs <- lapply(file_paths, reptoRlist)
#save
#saveRDS(sim_outputs, "sims/base/outputs/sim_outputs_base.rds")
#read
sim_outputs_base <- readRDS("sims/base/outputs/rep_simoutputs_base.rds")

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
fits_all <- purrr::map_dfr(1:nsims, function(i) {
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

sim_OM_biomass <- ggplot(bio_all, aes(x = year, y = pred_bio)) +
  # All simulation fits (light blue transparent)
  geom_line(
    data = filter(bio_all, source == "Simulation"),
    aes(group = interaction(sim, species)),
    color = "skyblue3", alpha = 0.2
  ) +
  # OM reference line (black)
  geom_line(
    data = filter(bio_all, source == "OM"),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "Biomass (t)",
    title = "Time series of OM biomass (black) and simulation fits (blue)"
  ) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Display plot
sim_OM_biomass

ggsave("sims/base/plots/sim_OM_biomass.png", plot = sim_OM_biomass, width = 10, height = 6, dpi = 300)


# Separate OM and Simulation
bio_OM <- bio_all %>% filter(source == "OM") %>% select(species, year, OM_bio = pred_bio)
bio_sim <- bio_all %>% filter(source == "Simulation")

# Join OM to Simulation and calculate relative error
bio_relerr <- bio_sim %>%
  left_join(bio_OM, by = c("species", "year")) %>%
  mutate(
    REE = (pred_bio - OM_bio) / OM_bio  # Relative error
  )

bio_relerr_base <- bio_sim %>%
  left_join(bio_OM, by = c("species", "year")) %>%
  mutate(
    REE = (pred_bio - OM_bio) / OM_bio  # Relative error
  )

bio_summary <- bio_relerr %>%
  group_by(species, year) %>%
  summarise(
    median_REE = median(REE, na.rm = TRUE),
    precision  = 1 / var(REE, na.rm = TRUE),
    .groups = "drop"
  )

bio_summary_base <- bio_relerr_base %>%
  group_by(species, year) %>%
  summarise(
    median_REE = median(REE, na.rm = TRUE),
    precision  = 1 / var(REE, na.rm = TRUE),
    .groups = "drop"
  )


rel_err_bio <- ggplot(bio_relerr, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "skyblue", alpha = 0.5, outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~ species, scales = "fixed") +
  coord_cartesian(ylim = c(-0.8, 1.5)) +
  theme_minimal() +
  scale_x_discrete(
    breaks = function(x) x[as.integer(x) %% 2 == 1]   # etiqueta cada 2 años
  ) +
  labs(
    title = "Relative Errors in Biomass Across Simulations",
    x = "Year",
    y = "Relative Error (Simulation vs OM)"
  )

rel_err_bio

ggsave("sims/base/plots/rel_err_bio.png", plot = rel_err_bio, width = 10, height = 6, dpi = 300)

rel_err_bio_low <- ggplot(bio_relerr_low, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "skyblue", alpha = 0.5, outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~ species, scales = "fixed") +
  coord_cartesian(ylim = c(-0.8, 1.5)) +
  theme_minimal() +
  scale_x_discrete(
    breaks = function(x) x[as.integer(x) %% 2 == 1]   # etiqueta cada 2 años
  ) +
  labs(
    title = "Relative Errors in Biomass Across Simulations",
    x = "Year",
    y = "Relative Error (Simulation vs OM)"
  )

rel_err_bio_low

ggsave("sims/base/plots/rel_err_bio.png", plot = rel_err_bio, width = 10, height = 6, dpi = 300)

# ============================================================
# Biomass REE summary over years (using bio_summary_base/low)
# Each year already contains median REE across sims
# We summarize across years to get overall median + IQR
# ============================================================

bio_relerr_both <- bind_rows(
  bio_relerr_base %>% mutate(scenario = "Base"),
  bio_relerr_low  %>% mutate(scenario = "Low interaction")
)

bio_sum_over_years <- bio_relerr_both %>%
  group_by(species, scenario) %>%
  summarise(
    median_REE = median(REE, na.rm = TRUE),
    q25_REE    = quantile(REE, 0.25, na.rm = TRUE),
    q75_REE    = quantile(REE, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p_bio_sum <- plot_relerr_summary(
  bio_sum_over_years,
  title = "Biomass REE summary over years: Base vs Low interaction",
  ylab  = "Median REE(Biomass)"
)

p_bio_sum

ggsave("sims/compare/Biomass_REE_summary_base_vs_low.png",
       plot = p_bio_sum, width = 12, height = 6, dpi = 300)


# Median REE over time
median_err_bio<-ggplot(bio_summary, aes(x = year, y = median_REE)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Median Relative Error in Biomass", y = "Median REE")

ggsave("sims/base/plots/median_err_bio.png", plot = median_err_bio, width = 10, height = 6, dpi = 300)

# Precision (1/Var) over time
ggplot(bio_summary, aes(x = year, y = precision)) +
  geom_line(color = "darkgreen") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Precision of Relative Errors in Biomass", y = "1 / Var(REE)")

ggplot(filter(bio_relerr, species == "Atlantic_cod"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Relative Errors: Atlantic cod", x = "REE", y = "Count")

ggplot(filter(bio_relerr, species == "Atlantic_herring"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Relative Errors: Atlantic herring", x = "REE", y = "Count")

ggplot(filter(bio_relerr, species == "Atlantic_mackerel"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Relative Errors: Atlantic mackerel", x = "REE", y = "Count")

ggplot(filter(bio_relerr, species == "Spiny_dogfish"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Relative Errors: Spiny dogfish", x = "REE", y = "Count")


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

catch_fits_all <- map_dfr(1:nsims, function(i) {
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


# Separate OM and simulation
catch_OM <- catch_all_df %>%
  filter(source == "OM") %>%
  select(fishery, year, species, OM_catch = pred_catch)

catch_sim <- catch_all_df %>%
  filter(source == "Simulation")

# Join OM to Simulation and compute REE
catch_relerr <- catch_sim %>%
  left_join(catch_OM, by = c("fishery", "year", "species")) %>%
  mutate(
    REE = (pred_catch - OM_catch) / OM_catch
  )

catch_relerr_low <- catch_sim %>%
  left_join(catch_OM, by = c("fishery", "year", "species")) %>%
  mutate(
    REE = (pred_catch - OM_catch) / OM_catch
  )

catch_summary <- catch_relerr %>%
  group_by(species, year) %>%
  summarise(
    median_REE = median(REE, na.rm = TRUE),
    precision  = 1 / var(REE, na.rm = TRUE),
    .groups = "drop"
  )

catch_summary_low <- catch_relerr_low %>%
  group_by(species, year) %>%
  summarise(
    median_REE = median(REE, na.rm = TRUE),
    precision  = 1 / var(REE, na.rm = TRUE),
    .groups = "drop"
  )

rel_err_catch<-ggplot(catch_relerr, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "orange", alpha = 0.5, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Relative Errors in Predicted Catch Across Simulations",
    x = "Year",
    y = "Relative Error (Simulation vs OM)"
  )


ggsave("sims/base/plots/rel_err_catch.png", plot = rel_err_catch, width = 10, height = 6, dpi = 300)

rel_err_catch_low<-ggplot(catch_relerr_low, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "orange", alpha = 0.5, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Relative Errors in Predicted Catch Across Simulations",
    x = "Year",
    y = "Relative Error (Simulation vs OM)"
  )


ggsave("sims/base/plots/rel_err_catch.png", plot = rel_err_catch_low, width = 10, height = 6, dpi = 300)

median_err_catch<-ggplot(catch_summary, aes(x = year, y = median_REE)) +
  geom_line(color = "orange", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Median Relative Error in Catch", y = "Median REE")

ggsave("sims/base/plots/median_err_catch.png", plot = median_err_catch, width = 10, height = 6, dpi = 300)

ggplot(catch_summary, aes(x = year, y = precision)) +
  geom_line(color = "darkred", linewidth = 1) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Precision of Relative Errors in Catch", y = "1 / Var(REE)")

ggplot(filter(catch_relerr, species == "Atlantic_cod"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "orange", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Catch – Atlantic cod",
    x = "Relative Error (REE)",
    y = "Count"
  )

ggplot(filter(catch_relerr, species == "Atlantic_herring"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "orange", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Catch – Atlantic herring",
    x = "Relative Error (REE)",
    y = "Count"
  )

ggplot(filter(catch_relerr, species == "Atlantic_mackerel"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "orange", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Catch – Atlantic mackerel",
    x = "Relative Error (REE)",
    y = "Count"
  )

ggplot(filter(catch_relerr, species == "Spiny_dogfish"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "orange", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Catch – Spiny dogfish",
    x = "Relative Error (REE)",
    y = "Count"
  )


sim_OM_catch <- ggplot(catch_all_df, aes(x = year, y = pred_catch)) +
  # All simulated runs (light red transparent)
  geom_line(
    data = filter(catch_all_df, source == "Simulation"),
    aes(group = interaction(sim, species)),
    color = "red", alpha = 0.2
  ) +
  # OM reference (black line)
  geom_line(
    data = filter(catch_all_df, source == "OM"),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "Catch (t)",
    title = "Time series of OM catch (black) and simulation fits (red)"
  ) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Display plot
sim_OM_catch

ggsave("sims/base/plots/sim_OM_catch.png", plot = sim_OM_catch, width = 10, height = 6, dpi = 300)


####################################
######## SURVEY INDEX
###############################

#READ FITS FROM RUNS WITH 100 DATA SETS
rep_simoutputs <- readRDS("sims/base/outputs/rep_simoutputs_base.rds")
#rep_files_sim <- paste0("sims/base/rep/hydra_sim", 1:100, ".rep")
#repfile <- "inputs/initial_run/hydra_sim.rep"

# Read all rep files into a named list
#rep_simoutputs <- lapply(rep_files_sim, function(file) {
#  cat("Reading:", file, "\n")
#  reptoRlist(file)
#})

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


####################################
######## SELECTIVITY
###############################

# ---------------------------------------------
# Setup
# ---------------------------------------------
Nspecies   <- length(hydraDataList$speciesList)  # 4
Nlen_sel   <- 5                                  # length bins in fishsel
fleet_labs <- c("Demersal", "Pelagic")

# helper: turn a fishsel matrix into long format
fishsel_to_long <- function(fishsel_mat, sim_label, source_label) {
  fs_df <- as.data.frame(fishsel_mat)
  
  # assume columns: species_id, fleet_id, then 5 length bins
  names(fs_df)[1:2] <- c("species_id", "fleet_id")
  
  fs_df %>%
    pivot_longer(
      cols = (3:(2 + Nlen_sel)),
      names_to  = "len_col",
      values_to = "selectivity"
    ) %>%
    mutate(
      len_bin = as.integer(gsub("V", "", len_col)) - 2,         # 1..5
      species = hydraDataList$speciesList[species_id],
      fleet   = factor(fleet_id, levels = 1:2, labels = fleet_labs),
      sim     = sim_label,
      source  = source_label
    ) %>%
    select(sim, source, species, fleet, len_bin, selectivity)
}

# ---------------------------------------------
# OM selectivity
# ---------------------------------------------
fishsel_OM <- fishsel_to_long(output$fishsel, sim_label = "OM", source_label = "OM")

# ---------------------------------------------
# Selectivity from 100 simulations
# ---------------------------------------------
fishsel_sims <- purrr::map_dfr(1:nsims, function(i) {
  fishsel_to_long(sim_outputs[[i]]$fishsel,
                  sim_label   = paste0("sim", i),
                  source_label = "Simulation")
})

# ---------------------------------------------
# Combine OM + simulations
# ---------------------------------------------
fishsel_all <- dplyr::bind_rows(fishsel_OM, fishsel_sims)

fishsel_OM_vals <- fishsel_all %>%
  filter(source == "OM") %>%
  select(species, fleet, len_bin, sel_OM = selectivity)

fishsel_relerr <- fishsel_all %>%
  filter(source == "Simulation") %>%
  left_join(fishsel_OM_vals,
            by = c("species", "fleet", "len_bin")) %>%
  mutate(REE = (selectivity - sel_OM) / sel_OM)


# ============================================================
# 4. Plot OM vs Simulation fits (by species & fleet)
# ============================================================
valid_pairs <- tribble(
  ~species,           ~fleet,
  "Atlantic_cod",     "Demersal",
  "Spiny_dogfish",    "Demersal",
  "Atlantic_herring", "Pelagic",
  "Atlantic_mackerel","Pelagic"
)

fishsel_filtered <- fishsel_all %>%
  inner_join(valid_pairs, by = c("species", "fleet"))

sim_OM_sel <- ggplot(fishsel_filtered,
                     aes(x = len_bin, y = selectivity)) +
  geom_line(
    data = filter(fishsel_filtered, source == "Simulation"),
    aes(group = interaction(sim, species, fleet)),
    color = "seagreen3", alpha = 0.25
  ) +
  geom_line(
    data = filter(fishsel_filtered, source == "OM"),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~species + fleet, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Length bin",
    y = "Selectivity",
    title = "OM (black) vs Simulation Selectivity (green) – Valid Species/Fleet Only"
  ) +
  scale_x_continuous(breaks = 1:5)

ggsave("sims/base/plots/sim_OM_selectivity_filtered.png",
     plot = sim_OM_sel, width = 10, height = 6, dpi = 300)


sel_OM <- fishsel_filtered %>%
  filter(source == "OM") %>%
  select(species, fleet, len_bin, OM_sel = selectivity)

sel_sim <- fishsel_filtered %>%
  filter(source == "Simulation")

sel_relerr <- sel_sim %>%
  left_join(sel_OM, by = c("species", "fleet", "len_bin")) %>%
  mutate(REE = (selectivity - OM_sel) / OM_sel)

rel_err_sel <- ggplot(sel_relerr,
                      aes(x = factor(len_bin), y = REE)) +
  geom_boxplot(fill = "seagreen3", alpha = 0.5, outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species + fleet, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(
    title = "Relative Errors in Selectivity – Valid Fleet–Species",
    x = "Length bin",
    y = "Relative Error"
  )

ggsave("sims/base/plots/rel_err_selectivity_filtered.png",
       plot = rel_err_sel, width = 10, height = 6, dpi = 300)


median_err_sel <- ggplot(sel_relerr %>% group_by(species, fleet, len_bin) %>%
                           summarise(median_REE = median(REE, na.rm = TRUE)),
                         aes(x = len_bin, y = median_REE)) +
  geom_line(color = "seagreen4") +
  geom_point(color = "seagreen4") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species + fleet, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(
    title = "Median REE in Selectivity – Valid Fleet–Species",
    x = "Length bin",
    y = "Median REE"
  )

ggsave("sims/base/plots/median_err_selectivity_filtered.png",
       plot = median_err_sel, width = 10, height = 6, dpi = 300)


precision_sel <- ggplot(sel_relerr %>% group_by(species, fleet, len_bin) %>%
                          summarise(precision = 1/var(REE, na.rm = TRUE)),
                        aes(x = len_bin, y = precision)) +
  geom_line(color = "darkgreen") +
  geom_point(color = "darkgreen") +
  facet_wrap(~species + fleet, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(
    title = "Precision of Selectivity – Valid Fleet–Species",
    x = "Length bin",
    y = "1 / Var(REE)"
  )

ggsave("sims/base/plots/precision_selectivity_filtered.png",
       plot = precision_sel, width = 10, height = 6, dpi = 300)

####################################
######## RECRUITMENT
###############################

# ============================================================
# 1. OM recruitment (from the deterministic or base run)
# ============================================================
est_recruits_OM <- output$EstRec %>%
  as.data.frame() %>%
  mutate(species = hydraDataList$speciesList) %>%
  select(species, everything()) %>%
  pivot_longer(cols = -species, names_to = "year") %>%
  mutate(
    year = as.integer(str_remove(year, "V")),
    log_rec = ifelse(value > 0, log(value), NA),
    pred_rec = value,
    sim = "OM",
    source = "OM"
  ) %>%
  select(year, species, pred_rec, sim, source)

# ============================================================
# 2. Simulation recruitments (from 74 replicated fits)
# ============================================================
recruit_fits_all <- map_dfr(1:nsims, function(i) {
  sim_outputs[[i]]$EstRec %>%
    as.data.frame() %>%
    mutate(species = hydraDataList$speciesList) %>%
    select(species, everything()) %>%
    pivot_longer(cols = -species, names_to = "year") %>%
    mutate(
      year = as.integer(str_remove(year, "V")),
      log_rec = ifelse(value > 0, log(value), NA),
      pred_rec = value,
      sim = paste0("sim", i),
      source = "Simulation"
    ) %>%
    select(year, species, pred_rec, sim, source)
})

# ============================================================
# 3. Combine OM + Simulation
# ============================================================
rec_all <- bind_rows(est_recruits_OM, recruit_fits_all)

# ============================================================
# 4. Plot OM vs Simulation fits
# ============================================================
sim_OM_recruitment <- ggplot(rec_all, aes(x = year, y = pred_rec)) +
  geom_line(
    data = filter(rec_all, source == "Simulation"),
    aes(group = interaction(sim, species)),
    color = "darkorange", alpha = 0.25
  ) +
  geom_line(
    data = filter(rec_all, source == "OM"),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "Recruitment (thousands)",
    title = "Time series of OM recruitment (black) and simulation fits (orange)"
  ) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

ggsave("sims/base/plots/sim_OM_recruitment.png",
       plot = sim_OM_recruitment, width = 10, height = 6, dpi = 300)


# ============================================================
# 5. Compute Relative Errors and Summaries
# ============================================================
rec_OM <- rec_all %>%
  filter(source == "OM") %>%
  select(species, year, OM_rec = pred_rec)

rec_sim <- rec_all %>%
  filter(source == "Simulation")

rec_relerr <- rec_sim %>%
  left_join(rec_OM, by = c("species", "year")) %>%
  mutate(REE = (pred_rec - OM_rec) / OM_rec)

rec_summary <- rec_relerr %>%
  group_by(species, year) %>%
  summarise(
    median_REE = median(REE, na.rm = TRUE),
    precision  = 1 / var(REE, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 6. Plot Relative Error Diagnostics
# ============================================================

# (a) Boxplots
rel_err_rec <- ggplot(rec_relerr, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "darkorange", alpha = 0.5, outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Relative Errors in Recruitment Across Simulations",
    x = "Year",
    y = "Relative Error (Simulation vs OM)"
  )

ggsave("sims/base/plots/rel_err_rec.png",
       plot = rel_err_rec, width = 10, height = 6, dpi = 300)

# (b) Median REE
median_err_rec <- ggplot(rec_summary, aes(x = year, y = median_REE)) +
  geom_line(color = "orange") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Median Relative Error in Recruitment", y = "Median REE")

ggsave("sims/base/plots/median_err_rec.png",
       plot = median_err_rec, width = 10, height = 6, dpi = 300)

# (c) Precision (1/Var)
precision_rec <- ggplot(rec_summary, aes(x = year, y = precision)) +
  geom_line(color = "darkred") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Precision of Relative Errors in Recruitment", y = "1 / Var(REE)")

# Atlantic cod
ggplot(filter(rec_relerr, species == "Atlantic_cod"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "darkorange", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Recruitment – Atlantic cod",
    x = "Relative Error (REE)",
    y = "Count"
  )

# Atlantic herring
ggplot(filter(rec_relerr, species == "Atlantic_herring"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "darkorange", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Recruitment – Atlantic herring",
    x = "Relative Error (REE)",
    y = "Count"
  )

# Atlantic mackerel
ggplot(filter(rec_relerr, species == "Atlantic_mackerel"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "darkorange", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Recruitment – Atlantic mackerel",
    x = "Relative Error (REE)",
    y = "Count"
  )

# Spiny dogfish
ggplot(filter(rec_relerr, species == "Spiny_dogfish"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "darkorange", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Recruitment – Spiny dogfish",
    x = "Relative Error (REE)",
    y = "Count"
  )

#######################################
###### F realized
###################################


#######################################
###### Helpers to parse .par files
#######################################

extract_vec_after_header <- function(lines, header_regex, n) {
  i <- which(str_detect(lines, header_regex))[1]
  if (is.na(i)) stop("Header not found: ", header_regex)
  
  vals <- numeric(0)
  j <- i + 1
  while (length(vals) < n && j <= length(lines)) {
    nums <- str_extract_all(
      lines[j],
      "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
    )[[1]]
    if (length(nums) > 0) vals <- c(vals, as.numeric(nums))
    j <- j + 1
  }
  
  if (length(vals) < n) {
    stop("Not enough values after header: ", header_regex,
         " (needed ", n, ", got ", length(vals), ")")
  }
  vals[1:n]
}

extract_mat_after_header <- function(lines, header_regex, nrow, ncol, byrow = TRUE) {
  n <- nrow * ncol
  vec <- extract_vec_after_header(lines, header_regex, n)
  matrix(vec, nrow = nrow, ncol = ncol, byrow = byrow)
}

# Parse one .par file -> returns list(F_base_df, q_df)
parse_par_for_F_and_q <- function(par_path, sim_id, source_id,
                                  nfleets = 2, nyrs = 42, nq = 2,
                                  q_species_order = c("Spiny_dogfish", "Atlantic_mackerel")) {
  
  lines <- readLines(par_path, warn = FALSE)
  
  # avg_F (length = nfleets)
  avg_F <- extract_vec_after_header(lines, "^\\s*#\\s*avg_F\\s*:", n = nfleets)
  avgF_df <- tibble(
    fleet = paste0("fleet", seq_len(nfleets)),
    avg_F = avg_F
  )
  
  # F_devs (nfleets x nyrs)
  Fdev_mat <- extract_mat_after_header(lines, "^\\s*#\\s*F_devs\\s*:", nrow = nfleets, ncol = nyrs, byrow = TRUE)
  
  Fdev_df <- as.data.frame(Fdev_mat) %>%
    mutate(fleet = paste0("fleet", seq_len(nfleets))) %>%
    pivot_longer(cols = -fleet, names_to = "year", values_to = "Fdev") %>%
    mutate(year = as.integer(str_remove(year, "^V"))) %>%
    arrange(fleet, year)
  
  # F_base = exp(avg_F + Fdev)
  F_base_df <- Fdev_df %>%
    left_join(avgF_df, by = "fleet") %>%
    mutate(
      sim = sim_id,
      source = source_id,
      F_base = exp(avg_F + Fdev)
    ) %>%
    select(sim, source, fleet, year, F_base)
  
  # ln_fishery_q (length = nq)
  ln_q <- extract_vec_after_header(lines, "^\\s*#\\s*ln_fishery_q\\s*:", n = nq)
  
  q_df <- tibble(
    sim = sim_id,
    source = source_id,
    species = q_species_order,
    ln_q = ln_q,
    q = exp(ln_q)
  )
  
  list(F_base_df = F_base_df, q_df = q_df)
}

#######################################
###### 1) Build list of .par files (OM + sims)
#######################################
OM_par <- "OM_scenarios/OM_base/hydra_sim.par"
sim_par_files <- list.files("sims/base/par", pattern = "\\.par$", full.names = TRUE)

par_files <- c(OM = OM_par)
par_files <- c(par_files, setNames(sim_par_files, paste0("sim", seq_along(sim_par_files))))

#######################################
###### 2) Parse all .par files -> F_base_all and q_all
#######################################
nfleets <- 2
nyrs    <- 42
nq      <- 2

parsed <- imap(par_files, ~{
  sim_id <- .y
  source_id <- ifelse(sim_id == "OM", "OM", "Simulation")
  parse_par_for_F_and_q(.x, sim_id = sim_id, source_id = source_id,
                        nfleets = nfleets, nyrs = nyrs, nq = nq)
})

F_base_all <- bind_rows(map(parsed, "F_base_df"))
q_all      <- bind_rows(map(parsed, "q_df"))

#######################################
###### 3) Selectivity at size_bin_4 for Fleet 2 (Pelagic)
###### fishsel_all must already exist
#######################################
sel_bin4_fleet2 <- fishsel_all %>%
  filter(fleet == "Pelagic",
         len_bin == 4,
         species %in% c("Atlantic_herring", "Atlantic_mackerel")) %>%
  transmute(sim, source, species, sel_bin4 = selectivity)

# Safety check: should be 1 row per (sim, source, species)
stopifnot(nrow(sel_bin4_fleet2 %>% count(sim, source, species) %>% filter(n > 1)) == 0)

#######################################
###### 4) Compute Realized_F (correct formulas)
#######################################

# Mapping: which species gets which fleet's F_base
# IMPORTANT: fleet names must match F_base_all$fleet (fleet1/fleet2)
species_fleet_map <- tibble(
  fleet   = c("fleet1","fleet1","fleet2","fleet2"),
  species = c("Atlantic_cod","Spiny_dogfish","Atlantic_herring","Atlantic_mackerel")
)

# q lookup: only dogfish & mackerel have q from par; cod/herring default to 1
q_lookup_all <- q_all %>%
  select(sim, source, species, q_mult = q)

RealizedF_all <- F_base_all %>%
  left_join(species_fleet_map, by = "fleet") %>%
  left_join(q_lookup_all, by = c("sim", "source", "species")) %>%
  mutate(q_mult = coalesce(q_mult, 1)) %>%
  left_join(sel_bin4_fleet2, by = c("sim", "source", "species")) %>%
  mutate(sel_bin4 = coalesce(sel_bin4, 1)) %>%
  mutate(
    Realized_F = case_when(
      fleet == "fleet1" & species == "Atlantic_cod"      ~ F_base,
      fleet == "fleet1" & species == "Spiny_dogfish"     ~ q_mult * F_base,
      fleet == "fleet2" & species == "Atlantic_herring"  ~ sel_bin4 * F_base,
      fleet == "fleet2" & species == "Atlantic_mackerel" ~ q_mult * sel_bin4 * F_base,
      TRUE ~ NA_real_
    )
  ) %>%
  select(sim, source, fleet, species, year, Realized_F)

#######################################
###### 5) Quick sanity checks
#######################################

# Should have 4 species for each sim-year (OM and each sim)
RealizedF_all %>%
  group_by(sim, source) %>%
  summarise(
    n_species = n_distinct(species),
    n_years   = n_distinct(year),
    n_rows    = n(),
    .groups = "drop"
  ) %>%
  arrange(source, sim) %>%
  print(n = 50)

# Ratio check: dogfish/cod should be constant in time within each sim (equals q_dogfish)
ratio_check <- RealizedF_all %>%
  filter(species %in% c("Atlantic_cod","Spiny_dogfish")) %>%
  select(sim, source, year, species, Realized_F) %>%
  pivot_wider(names_from = species, values_from = Realized_F) %>%
  mutate(ratio_dogfish_cod = Spiny_dogfish / Atlantic_cod)

ratio_check %>%
  group_by(sim, source) %>%
  summarise(
    ratio_min = min(ratio_dogfish_cod, na.rm = TRUE),
    ratio_max = max(ratio_dogfish_cod, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = 20)

##### save scenario base 
RealizedF_all_base <- RealizedF_all
RealizedF_all_low <- RealizedF_all

#rm(RealizedF_all, Freal_relerr, Freal_summary)

#######################################
###### 6) Plot OM vs simulations (time series)
#######################################
p_Freal <- ggplot(RealizedF_all, aes(x = year, y = Realized_F)) +
  geom_line(
    data = filter(RealizedF_all, source == "Simulation"),
    aes(group = interaction(sim, species)),
    color = "red", alpha = 0.25
  ) +
  geom_line(
    data = filter(RealizedF_all, source == "OM"),
    aes(group = species),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "Realized fishing mortality (F)",
    title = "OM realized F (black) and EM fits (red)"
  )

p_Freal

 ggsave("sims/base/plots/sim_OM_Frealized_timeseries.png",
        plot = p_Freal, width = 10, height = 6, dpi = 300)


# ============================================================
# 5. Compute Relative Errors and Summaries (Realized F)
# ============================================================

Freal_OM <- RealizedF_all %>%
  filter(source == "OM") %>%
  select(species, year, OM_F = Realized_F)

Freal_sim <- RealizedF_all %>%
  filter(source == "Simulation")

Freal_relerr <- Freal_sim %>%
  left_join(Freal_OM, by = c("species", "year")) %>%
  filter(!is.na(OM_F), OM_F > 0) %>%
  mutate(REE = (Realized_F - OM_F) / OM_F)

# Boxplots of REE by year
rel_err_Freal <- ggplot(Freal_relerr, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "red", alpha = 0.5, outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  scale_x_discrete(
    breaks = function(x) x[as.integer(x) %% 2 == 1]  # show every 2 years
  ) +
  labs(
    title = "Relative Errors in Realized F (EM vs OM)",
    x = "Year",
    y = "Relative Error (REE)"
  )

rel_err_Freal


# Optional save
ggsave("sims/base/plots/rel_err_Freal.png",
       plot = rel_err_Freal, width = 10, height = 6, dpi = 300)



# Freal_relerr must already exist from the realized-F script:
# columns: sim, source, fleet, species, year, Realized_F, OM_F, REE

Freal_summary <- Freal_relerr %>%
  group_by(species, year) %>%
  summarise(
    median_REE = median(REE, na.rm = TRUE),
    precision  = {
      v <- var(REE, na.rm = TRUE)
      ifelse(is.na(v) | v == 0, NA_real_, 1 / v)
    },
    .groups = "drop"
  )

# ============================================================
# 2) Plot (b) Median REE time series
# ============================================================

median_err_Freal <- ggplot(Freal_summary, aes(x = year, y = median_REE)) +
  geom_line(color = "red") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Median Relative Error in Realized Fishing Mortality (F)",
    x = "Year",
    y = "Median REE"
  )

median_err_Freal

 ggsave("sims/base/plots/median_err_Freal.png",
        plot = median_err_Freal, width = 10, height = 6, dpi = 300)

# ============================================================
# 3) Plot (c) Precision (1/Var) time series
# ============================================================

precision_Freal <- ggplot(Freal_summary, aes(x = year, y = precision)) +
  geom_line(color = "red") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Precision of Relative Errors in Realized Fishing Mortality (F)",
    x = "Year",
    y = "1 / Var(REE)"
  )

precision_Freal

 ggsave("sims/base/plots/precision_Freal.png",
        plot = precision_Freal, width = 10, height = 6, dpi = 300)

# ============================================================
# 4) Histograms by species (Realized F REE)
# ============================================================

hist_F_cod <- ggplot(filter(Freal_relerr, species == "Atlantic_cod"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "red", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Realized F – Atlantic cod",
    x = "Relative Error (REE)",
    y = "Count"
  )

hist_F_dogfish <- ggplot(filter(Freal_relerr, species == "Spiny_dogfish"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "red", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Realized F – Spiny dogfish",
    x = "Relative Error (REE)",
    y = "Count"
  )

hist_F_herring <- ggplot(filter(Freal_relerr, species == "Atlantic_herring"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "red", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Realized F – Atlantic herring",
    x = "Relative Error (REE)",
    y = "Count"
  )

hist_F_mackerel <- ggplot(filter(Freal_relerr, species == "Atlantic_mackerel"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "red", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Realized F – Atlantic mackerel",
    x = "Relative Error (REE)",
    y = "Count"
  )

hist_F_cod
hist_F_dogfish
hist_F_herring
hist_F_mackerel

# Optional saves:
# ggsave("sims/base/plots/hist_Freal_cod.png",      hist_F_cod,      8, 5, dpi = 300)
# ggsave("sims/base/plots/hist_Freal_dogfish.png",  hist_F_dogfish,  8, 5, dpi = 300)
# ggsave("sims/base/plots/hist_Freal_herring.png",  hist_F_herring,  8, 5, dpi = 300)
# ggsave("sims/compare/hist_Freal_mackerel.png", hist_F_mackerel, 8, 5, dpi = 300)

# ============================================================
# 5) Base vs Low scenario comparison 
# ============================================================

#   RealizedF_all_base
#   RealizedF_all_low
# with columns: sim, source, fleet, species, year, Realized_F

# ---------- BASE ----------
#Freal_base_OM <- RealizedF_all_base %>%
 # filter(source == "OM") %>%
  #select(species, year, OM_F = Realized_F)

# ---------- BASE ----------
Freal_base_OM <- RealizedF_all %>%
  filter(source == "OM") %>%
  select(species, year, OM_F = Realized_F)

# ---------- LOW INTERACTION ----------
Freal_low_OM <- RealizedF_all %>%
  filter(source == "OM") %>%
  select(species, year, OM_F = Realized_F)


# ---------- BASE: ---------
Freal_base_std <- RealizedF_all_base %>%
  filter(source == "Simulation") %>%
  left_join(Freal_base_OM, by = c("species","year")) %>%
  filter(!is.na(OM_F), OM_F > 0) %>%
  mutate(
    REE = (Realized_F - OM_F) / OM_F,
    metric = "Realized F (formula)",
    scenario = "Base"
  ) %>%
  select(year, species, sim, source,
         predicted = Realized_F, OM_value = OM_F, REE, metric, scenario)

# ---------- LOW INTERACTION: ----------
Freal_low_std <- RealizedF_all %>%
  filter(source == "Simulation") %>%
  left_join(Freal_low_OM, by = c("species","year")) %>%
  filter(!is.na(OM_F), OM_F > 0) %>%
  mutate(
    REE = (Realized_F - OM_F) / OM_F,
    metric = "Realized F (formula)",
    scenario = "Low interaction"
  ) %>%
  select(year, species, sim, source,
         predicted = Realized_F, OM_value = OM_F, REE, metric, scenario)

# ---------- Bind both scenarios ----------
Freal_std_both <- bind_rows(Freal_base_std, Freal_low_std) %>%
  filter(source == "Simulation")

summarize_relerr_over_years <- function(df) {
  df %>%
    group_by(species, scenario) %>%
    summarise(
      median_REE = median(REE, na.rm = TRUE),
      q25_REE    = quantile(REE, 0.25, na.rm = TRUE),
      q75_REE    = quantile(REE, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_relerr_summary <- function(df, title, ylab) {
  ggplot(df, aes(x = species, y = median_REE, color = scenario)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(
      aes(ymin = q25_REE, ymax = q75_REE),
      width = 0.2,
      position = position_dodge(width = 0.5)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal(base_size = 12) +
    labs(
      title = title,
      x = "Species",
      y = ylab,
      color = "Scenario"
    )
}

# --------------------
# summarize_relerr_over_years() summarize across years and sims
Freal_sum_over_years <- summarize_relerr_over_years(Freal_std_both)

# plot_relerr_summary() should produce your "median + IQR" style plot
p_F_sum_new <- plot_relerr_summary(
  Freal_sum_over_years,
  title = "Realized F (formula) REE summary over years: Base vs Low interaction",
  ylab  = "Median REE(F_realized)"
)

p_F_sum_new

ggsave("sims/compare/Freal_formula_REE_summary_base_vs_low.png",
       plot = p_F_sum_new, width = 12, height = 6, dpi = 300)


#######################################
###### M2 (Predation Mortality)
###################################

# ============================================================
# 1. OM M2 (from OM)
# ============================================================
nlen <- ncol(output$EstM2size)

est_M2_OM <- output$EstM2size %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.), names_to = "ilen", names_prefix = "V") %>%
  mutate(
    species = rep(hydraDataList$speciesList,
                  each = hydraDataList$Nyrs * nlen * stepperyr),
    year = rep(rep(1:(hydraDataList$Nyrs * stepperyr), each = nlen),
               length(hydraDataList$speciesList)),
    year = (1 - 1 / stepperyr) + year / stepperyr,
    sim = "OM",
    source = "OM",
    pred_M2 = value
  ) %>%
  # Aggregate across length bins to get total M2 per species-year
  group_by(species, year, sim, source) %>%
  summarise(pred_M2 = sum(pred_M2, na.rm = TRUE), .groups = "drop")

# ============================================================
# 2. Simulation M2 (from 74 simulation fits)
# ============================================================
M2_fits_all <- map_dfr(1:nsims, function(i) {
  sim_outputs[[i]]$EstM2size %>%
    as.data.frame() %>%
    pivot_longer(cols = 1:ncol(.), names_to = "ilen", names_prefix = "V") %>%
    mutate(
      species = rep(hydraDataList$speciesList,
                    each = hydraDataList$Nyrs * nlen * stepperyr),
      year = rep(rep(1:(hydraDataList$Nyrs * stepperyr), each = nlen),
                 length(hydraDataList$speciesList)),
      year = (1 - 1 / stepperyr) + year / stepperyr,
      sim = paste0("sim", i),
      source = "Simulation",
      pred_M2 = value
    ) %>%
    # Aggregate across length bins to get total M2 per species-year
    group_by(species, year, sim, source) %>%
    summarise(pred_M2 = sum(pred_M2, na.rm = TRUE), .groups = "drop")
})

# ============================================================
# 3. Combine OM + Simulation
# ============================================================
#M2_all_low <- bind_rows(est_M2_OM, M2_fits_all)
M2_all <- bind_rows(est_M2_OM, M2_fits_all)
# ============================================================
# 4. Plot OM (black) and simulation fits (green)
# ============================================================
sim_OM_M2 <- ggplot(M2_all, aes(x = year, y = pred_M2)) +
  geom_line(
    data = filter(M2_all, source == "Simulation"),
    aes(group = interaction(sim, species)),
    color = "seagreen3", alpha = 0.25
  ) +
  geom_line(
    data = filter(M2_all, source == "OM"),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "Predation mortality (M2)",
    title = "Time series of OM predation mortality (black) and simulation fits (green)"
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

#sim_OM_M2 <- ggplot(M2_all_low, aes(x = year, y = pred_M2)) +
#  geom_line(
#    data = filter(M2_all_low, source == "Simulation"),
#    aes(group = interaction(sim, species)),
#    color = "seagreen3", alpha = 0.25
#  ) +
#  geom_line(
#    data = filter(M2_all_low, source == "OM"),
#    color = "black", linewidth = 1
#  ) +
#  facet_wrap(~species, scales = "free_y") +
#  theme_minimal(base_size = 12) +
#  labs(
#    x = "Year",
#    y = "Predation mortality (M2)",
#    title = "Time series of OM predation mortality (black) and simulation fits (green)"
#  ) +
#  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
#  theme(
#    plot.title = element_text(face = "bold", hjust = 0.5),
#    panel.grid.minor = element_blank()
#  )


#######################
#### PAY ATTENTION TO WHAT YOU ARE DOING!!!#####

ggsave("sims/base/plots/sim_OM_M2.png",
       plot = sim_OM_M2, width = 10, height = 6, dpi = 300)

# ============================================================
# 5. Compute Relative Errors and Summaries
# ============================================================
M2_OM <- M2_all %>%
  filter(source == "OM") %>%
  select(species, year, OM_M2 = pred_M2)

M2_sim <- M2_all %>%
  filter(source == "Simulation")

M2_relerr <- M2_sim %>%
  left_join(M2_OM, by = c("species", "year")) %>%
  mutate(REE = (pred_M2 - OM_M2) / OM_M2)


M2_summary <- M2_relerr %>%
  group_by(species, year) %>%
  summarise(
    median_REE = median(REE, na.rm = TRUE),
    precision  = 1 / var(REE, na.rm = TRUE),
    .groups = "drop"
  )

M2_summary_base <- M2_relerr

# ============================================================
# 6. Plot Relative Error Diagnostics
# ============================================================

# (a) Boxplots
rel_err_M2 <- ggplot(M2_relerr, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "seagreen3", alpha = 0.5, outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Relative Errors in Predation Mortality (M2) Across Simulations",
    x = "Year",
    y = "Relative Error (Simulation vs OM)"
  )

ggsave("sims/base/plots/rel_err_M2.png",
       plot = rel_err_M2, width = 10, height = 6, dpi = 300)

# (b) Median REE
median_err_M2 <- ggplot(M2_summary, aes(x = year, y = median_REE)) +
  geom_line(color = "seagreen4") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Median Relative Error in Predation Mortality (M2)", y = "Median REE")

ggsave("sims/base/plots/median_err_M2.png",
       plot = median_err_M2, width = 10, height = 6, dpi = 300)

# (c) Precision (1/Var)
precision_M2 <- ggplot(M2_summary, aes(x = year, y = precision)) +
  geom_line(color = "darkgreen") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Precision of Relative Errors in Predation Mortality (M2)", y = "1 / Var(REE)")

ggsave("sims/base/plots/precision_M2.png",
        plot = precision_M2, width = 10, height = 6, dpi = 300)


# ============================================================
# 7. Histograms by Species
# ============================================================

# Atlantic herring
ggplot(filter(M2_relerr, species == "Atlantic_herring"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "seagreen3", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Predation Mortality – Atlantic herring",
    x = "Relative Error (REE)",
    y = "Count"
  )

# Atlantic mackerel
ggplot(filter(M2_relerr, species == "Atlantic_mackerel"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "seagreen3", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Predation Mortality – Atlantic mackerel",
    x = "Relative Error (REE)",
    y = "Count"
  )

#######################################
###### Z (Total Mortality)
###################################

# ============================================================
# 1. OM M1 (from OM)
# ============================================================
nlen_M1 <- ncol(output$EstM1size)   

est_M1_OM <- output$EstM1size %>%
  as.data.frame() %>%
  transmute(pred_M1 = V1) %>%     # KEEP ONLY FIRST COLUMN
  mutate(
    species = rep(hydraDataList$speciesList,
                  each = hydraDataList$Nyrs * stepperyr),
    year = rep(1:(hydraDataList$Nyrs * stepperyr),
               times = hydraDataList$Nspecies),
    year = (1 - 1/stepperyr) + year / stepperyr,
    sim = "OM",
    source = "OM"
  )

# ============================================================
# 2. Simulation M1 (from 100 simulation fits)
# ============================================================
M1_fits_all <- map_dfr(1:nsims, function(i) {
  sim_outputs[[i]]$EstM1size[, 1] %>%      # ⬅️ FIRST COLUMN ONLY
    as.data.frame() %>%
    rename(pred_M1 = ".") %>%             # name the column
    mutate(
      species = rep(hydraDataList$speciesList,
                    each  = hydraDataList$Nyrs * stepperyr),
      year    = rep(1:(hydraDataList$Nyrs * stepperyr),
                    times = hydraDataList$Nspecies),
      year    = (1 - 1/stepperyr) + year / stepperyr,
      sim     = paste0("sim", i),
      source  = "Simulation"
    )
})

# ============================================================
# 3. Combine OM + Simulation M1
# ============================================================
M1_all <- bind_rows(est_M1_OM, M1_fits_all)

#######################################
###### Compute Z = F + M1 + M2
###################################

Z_all <- RealizedF_all %>%
  left_join(M1_all, by = c("species", "year", "sim", "source")) %>%
  left_join(M2_all, by = c("species", "year", "sim", "source")) %>%
  mutate(
    pred_F = Realized_F, #pred_F,
    pred_M = pred_M1 + pred_M2,
    pred_Z = pred_F + pred_M
  )

# ============================================================
# 4. Plot OM (black) and simulation fits (purple)
# ============================================================
sim_OM_Z <- ggplot(Z_all, aes(x = year, y = pred_Z)) +
  geom_line(
    data = filter(Z_all, source == "Simulation"),
    aes(group = interaction(sim, species)),
    color = "purple3", alpha = 0.25
  ) +
  geom_line(
    data = filter(Z_all, source == "OM"),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "Total mortality (Z)",
    title = "Time series of OM total mortality Z (black) and simulation fits (purple)"
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

ggsave("sims/base/plots/sim_OM_Z.png",
       plot = sim_OM_Z, width = 10, height = 6, dpi = 300)

###################################
####with realized F

#######################################
###### Compute Z = Realized F + M1 + M2
###################################

Z_all <- RealizedF_all %>%
  left_join(M1_all, by = c("species", "year", "sim", "source")) %>%
  left_join(M2_all, by = c("species", "year", "sim", "source")) %>%
  mutate(
    pred_M = pred_M1 + pred_M2,
    pred_Z = Realized_F + pred_M
  )

# ============================================================
# Plot OM (black) and simulation fits (purple)
# ============================================================
sim_OM_Z <- ggplot(Z_all, aes(x = year, y = pred_Z)) +
  geom_line(
    data = dplyr::filter(Z_all, source == "Simulation"),
    aes(group = interaction(sim, species)),
    color = "purple3", alpha = 0.25
  ) +
  geom_line(
    data = dplyr::filter(Z_all, source == "OM"),
    aes(group = species),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "Total mortality (Z)",
    title = "Time series of OM total mortality Z (black) and simulation fits (purple)"
  ) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

ggsave(
  "sims/base/plots/sim_OM_Z_realizedF.png",
  plot = sim_OM_Z,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# 5. Compute Relative Errors and Summaries for Z
# ============================================================
Z_OM <- Z_all %>%
  filter(source == "OM") %>%
  select(species, year, OM_Z = pred_Z)

Z_sim <- Z_all %>%
  filter(source == "Simulation")

Z_relerr <- Z_sim %>%
  left_join(Z_OM, by = c("species", "year")) %>%
  mutate(REE = (pred_Z - OM_Z) / OM_Z)

Z_summary <- Z_relerr %>%
  group_by(species, year) %>%
  summarise(
    median_REE = median(REE, na.rm = TRUE),
    precision  = 1 / var(REE, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 6. Plot Relative Error Diagnostics for Z
# ============================================================

# (a) Boxplots
rel_err_Z <- ggplot(Z_relerr, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "plum3", alpha = 0.5, outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Relative Errors in Total Mortality (Z) Across Simulations",
    x = "Year",
    y = "Relative Error (Simulation vs OM)"
  )

ggsave("sims/base/plots/rel_err_Z.png",
       plot = rel_err_Z, width = 10, height = 6, dpi = 300)

# (b) Median REE
median_err_Z <- ggplot(Z_summary, aes(x = year, y = median_REE)) +
  geom_line(color = "purple4") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Median Relative Error in Total Mortality (Z)", y = "Median REE")

ggsave("sims/base/plots/median_err_Z.png",
       plot = median_err_Z, width = 10, height = 6, dpi = 300)

# (c) Precision (1/Var)
precision_Z <- ggplot(Z_summary, aes(x = year, y = precision)) +
  geom_line(color = "darkmagenta") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Precision of Relative Errors in Total Mortality (Z)", y = "1 / Var(REE)")

ggsave("sims/base/plots/precision_Z.png",
       plot = precision_Z, width = 10, height = 6, dpi = 300)


# ============================================================
# 7. Histograms by Species – Z
# ============================================================

# Atlantic cod
ggplot(filter(Z_relerr, species == "Atlantic_cod"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "plum3", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Total Mortality (Z) – Atlantic cod",
    x = "Relative Error (REE)",
    y = "Count"
  )

# Atlantic herring
ggplot(filter(Z_relerr, species == "Atlantic_herring"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "plum3", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Total Mortality (Z) – Atlantic herring",
    x = "Relative Error (REE)",
    y = "Count"
  )

# Atlantic mackerel
ggplot(filter(Z_relerr, species == "Atlantic_mackerel"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "plum3", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Total Mortality (Z) – Atlantic mackerel",
    x = "Relative Error (REE)",
    y = "Count"
  )

# Spiny dogfish
ggplot(filter(Z_relerr, species == "Spiny_dogfish"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "plum3", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Total Mortality (Z) – Spiny dogfish",
    x = "Relative Error (REE)",
    y = "Count"
  )

######################################################
######## summary

# Add a column indicating variable type for bookkeeping
bio_summary  <- bio_summary  %>% mutate(metric = "Biomass")
catch_summary <- catch_summary %>% mutate(metric = "Catch")
rec_summary  <- rec_summary  %>% mutate(metric = "Recruitment")
F_summary    <- F_summary    %>% mutate(metric = "Fishing mortality (F)")
M2_summary  <- M2_summary  %>% mutate(metric = "M2")
Z_summary  <- Z_summary  %>% mutate(metric = "Z")

# ============================================================
# 2. Save individual CSVs
# ============================================================
write_csv(bio_summary,   "sims/base/plots/summary_median_precision_biomass.csv")
write_csv(catch_summary, "sims/base/plotssummary_median_precision_catch.csv")
write_csv(rec_summary,   "sims/base/plots/summary_median_precision_recruitment.csv")
write_csv(F_summary,     "sims/base/plots/summary_median_precision_F.csv")
write_csv(M2_summary,     "sims/base/plots/summary_median_precision_M2.csv")
write_csv(Z_summary,     "sims/base/plots/summary_median_precision_Z.csv")
# ============================================================
# 3. (Optional) Combine all into one master summary
# ============================================================
summary_all <- bind_rows(bio_summary, catch_summary, rec_summary, F_summary, M2_summary, Z_summary)
write_csv(summary_all, "sims/base/plots/summary_median_precision_all.csv")


###############################################
############## PARAMETERS

# OM .par file
om_par_path <- "OM_scenarios/OM_base/hydra_sim.par"              # <- edit if different

# Simulation .par files (sim1 ... sim100)
sim_par_paths <- sprintf("sims/base/par/hydra_sim%d.par", 1:nsims)

# Parameters of interest (exact base names *before* any [index])
params_keep <- c(
  "ln_yr1N",
  "ln_avg_recruitment", "recruitment_devs", "F_devs", "avg_F", "fishsel_pars", "ln_fishery_q",
  "ln_survey_q")

# ------------------------------------------------------------------------------
# 1) Reader that is resilient to whitespace and vector-style names like foo[12]
# ------------------------------------------------------------------------------
read_par_block <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  # Keep only lines that aren't completely empty
  lines <- lines[nzchar(trimws(lines))]
  
  current_param <- NULL
  out <- list()
  
  for (ln in lines) {
    # Detect new parameter header
    if (grepl("^# ", ln)) {
      current_param <- gsub("[:#]", "", trimws(ln))
    } else if (!is.null(current_param)) {
      # Read numeric values under the current parameter
      nums <- as.numeric(unlist(strsplit(trimws(ln), "\\s+")))
      if (length(nums) > 0) {
        out[[length(out) + 1]] <- tibble(
          name = current_param,
          value = nums
        )
      }
    }
  }
  
  df <- bind_rows(out) %>%
    mutate(
      param_base = trimws(name),
      idx = ave(value, param_base, FUN = seq_along)
    )
  
  df
}

# ============================================================
# 2. READ OM .par FILE
# ============================================================
par_om <- read_par_block(om_par_path) %>%
  filter(param_base %in% params_keep) %>%
  mutate(sim = "OM", source = "OM")

# ============================================================
# 3. READ SIMULATION .par FILES (74 runs)
# ============================================================
read_sim_par_block <- function(path, sim_id) {
  tryCatch(
    read_par_block(path) %>%
      filter(param_base %in% params_keep) %>%
      mutate(sim = paste0("sim", sim_id), source = "Simulation"),
    error = function(e) {
      warning("Failed to read: ", path)
      NULL
    }
  )
}

par_sims <- map2_dfr(sim_par_paths, 1:nsims, read_sim_par_block)
par_all <- bind_rows(par_om, par_sims)

par_all <- par_all %>%
  mutate(value = exp(value))

# ============================================================
# 5. COMPUTE RELATIVE ERRORS
# ============================================================
par_OM <- par_all %>%
  filter(source == "OM") %>%
  select(param_base, idx, OM_value = value)

par_sim <- par_all %>%
  filter(source == "Simulation")

par_relerr <- par_sim %>%
  left_join(par_OM, by = c("param_base", "idx")) %>%
  mutate(
    REE  = (value - OM_value) / OM_value,
    diff = value - OM_value
  )

# ============================================================
# 6. SUMMARIES
# ============================================================
par_summary_detail <- par_relerr %>%
  group_by(param_base, idx) %>%
  summarise(
    bias_REE      = mean(REE, na.rm = TRUE),                 # mean relative error
    median_REE    = median(REE, na.rm = TRUE),
    median_abs_REE= median(abs(REE), na.rm = TRUE),
    rmse_REE      = sqrt(mean(REE^2, na.rm = TRUE)),         # RMSE in relative error
    emp_sd        = sd(diff, na.rm = TRUE),                  # empirical SD of (hat - true)
    emp_sd_rel    = sd(REE, na.rm = TRUE),                   # SD of relative error
    precision_emp = ifelse(emp_sd > 0, 1 / emp_sd^2, NA_real_),
    .groups = "drop"
  )

# Save summary
write.csv(par_summary_detail, "sims/base/plots/summary_par_REE.csv", row.names = FALSE)

# ============================================================
# 7. PLOTS
# ============================================================

# (a) Boxplots of REE
# Directory where plots will be saved
out_dir <- "sims/base/plots/params_by_idx/"
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

Nyrs      <- 42                      # years
Nyrs_rec  <- Nyrs - 1 
Nspecies  <- length(hydraDataList$speciesList)  # 4
Nlen      <- 5                       # length bins
fleet_labs <- c("Demersal", "Pelagic")

special_params <- c("F_devs", "ln_yr1N", "recruitment_devs")

param_simple <- setdiff(unique(par_relerr$param_base), special_params)

for (p in param_simple) {
  df_p <- par_relerr %>% filter(param_base == p)
  
  plt <- ggplot(df_p, aes(x = factor(idx), y = REE)) +
    geom_boxplot(fill = "gray70", alpha = 0.6, outlier.size = 0.5) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Relative Errors for", p),
      x = "Index within parameter (idx)",
      y = "Relative Error (Simulation vs OM)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = paste0(out_dir, "rel_err_", p, ".png"),
    plot = plt,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# F_devs: 2 fleets x 42 years
par_Fdevs <- par_relerr %>%
  filter(param_base == "F_devs") %>%
  mutate(
    fleet_id = ((idx - 1) %/% Nyrs) + 1,           # 1 or 2
    year     = ((idx - 1) %%  Nyrs) + 1,
    fishery  = factor(fleet_id, levels = 1:2, labels = fleet_labs)
  )

# Demersal-only plot
df_dem <- par_Fdevs %>% filter(fishery == fleet_labs[1])

plt_dem <- ggplot(df_dem, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "gray70", alpha = 0.6, outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(
    title = paste("Relative Errors for F_devs –", fleet_labs[1], "fishery"),
    x = "Year",
    y = "Relative Error (Simulation vs OM)"
  )

ggsave(
  filename = paste0(out_dir, "rel_err_F_devs_", fleet_labs[1], ".png"),
  plot = plt_dem,
  width = 8, height = 5, dpi = 300
)

# Pelagic-only plot
df_pel <- par_Fdevs %>% filter(fishery == fleet_labs[2])

plt_pel <- ggplot(df_pel, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "gray70", alpha = 0.6, outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(
    title = paste("Relative Errors for F_devs –", fleet_labs[2], "fishery"),
    x = "Year",
    y = "Relative Error (Simulation vs OM)"
  )

ggsave(
  filename = paste0(out_dir, "rel_err_F_devs_", fleet_labs[2], ".png"),
  plot = plt_pel,
  width = 8, height = 5, dpi = 300
)


par_lnN <- par_relerr %>%
  filter(param_base == "ln_yr1N") %>%
  mutate(
    species_id = ((idx - 1) %/% Nlen) + 1,           # 1..4
    len_bin    = ((idx - 1) %%  Nlen) + 1,
    species    = hydraDataList$speciesList[species_id]
  )

species_list <- unique(par_lnN$species)

for (sp in species_list) {
  df_sp <- par_lnN %>% filter(species == sp)
  
  plt_sp <- ggplot(df_sp, aes(x = factor(len_bin), y = REE)) +
    geom_boxplot(fill = "gray70", alpha = 0.6, outlier.size = 0.5) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Relative Errors for Initial Abundance –", sp),
      x = "Length bin",
      y = "Relative Error (Simulation vs OM)"
    )
  
  ggsave(
    filename = paste0(out_dir, "rel_err_ln_yr1N_", sp, ".png"),
    plot = plt_sp,
    width = 8, height = 5, dpi = 300
  )
}

par_recdevs <- par_relerr %>%
  filter(param_base == "recruitment_devs") %>%
  mutate(
    species_id = ((idx - 1) %/% Nyrs_rec) + 1,   # 1..4
    dev_idx    = ((idx - 1) %%  Nyrs_rec) + 1,   # 1..41 (index of deviation)
    year       = dev_idx + 1,                    # map devs to years 2..42
    species    = hydraDataList$speciesList[species_id]
  )

species_list_rec <- unique(par_recdevs$species)

for (sp in species_list_rec) {
  df_sp <- par_recdevs %>% filter(species == sp)
  
  plt_sp <- ggplot(df_sp, aes(x = factor(year), y = REE)) +
    geom_boxplot(fill = "gray70", alpha = 0.6, outlier.size = 0.5) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Relative Errors for recruitment_devs –", sp),
      x = "Year",
      y = "Relative Error (Simulation vs OM)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = paste0(out_dir, "rel_err_recruitment_devs_", sp, ".png"),
    plot = plt_sp,
    width = 9, height = 5, dpi = 300
  )
}

# Make a copy to modify
par_relerr_tagged <- par_relerr %>%
  mutate(
    species = NA_character_,
    year     = NA_integer_,
    len_bin  = NA_integer_,
    fishery  = NA_character_
  )

par_relerr_tagged <- par_relerr_tagged %>%
  mutate(
    species = if_else(
      param_base == "recruitment_devs",
      hydraDataList$speciesList[((idx - 1) %/% Nyrs_rec) + 1],
      species
    ),
    year = if_else(
      param_base == "recruitment_devs",
      ((idx - 1) %% Nyrs_rec) + 2,  # years 2..42
      year
    )
  )

par_relerr_tagged <- par_relerr_tagged %>%
  mutate(
    species = if_else(
      param_base == "ln_yr1N",
      hydraDataList$speciesList[((idx - 1) %/% Nlen) + 1],
      species
    ),
    len_bin = if_else(param_base == "ln_yr1N",
                      ((idx - 1) %% Nlen) + 1,
                      len_bin)
  )

par_relerr_tagged <- par_relerr_tagged %>%
  mutate(
    species = if_else(
      param_base == "recruitment_devs",
      hydraDataList$speciesList[((idx - 1) %/% Nyrs) + 1],
      species
    ),
    year = if_else(
      param_base == "recruitment_devs",
      ((idx - 1) %% Nyrs) + 1,
      year
    )
  )

par_relerr_out <- par_relerr_tagged %>%
  select(sim, param_base, idx,
         species, year, len_bin, fishery,
         value, OM_value, REE, source) %>%
  arrange(param_base, idx, sim)


write.csv(par_relerr_out,
          "sims/base/plots/relative_errors_parameters_tagged.csv",
          row.names = FALSE)



###############################################
### OVERALL GOODNESS OF FIT SUMMARY PER SIM ###
###############################################

library(dplyr)

# Combine all relative-error tables (biomass, catch, recruitment, F, M2)
bio_relerr   <- bio_relerr  %>% mutate(metric = "Biomass")
catch_relerr <- catch_relerr %>% mutate(metric = "Catch")
rec_relerr   <- rec_relerr  %>% mutate(metric = "Recruitment")
F_relerr     <- F_relerr    %>% mutate(metric = "Fishing mortality (F)")
M2_relerr    <- M2_relerr   %>% mutate(metric = "M2")
Z_relerr    <- Z_relerr   %>% mutate(metric = "Z")

# 1. Bind them all together
relerr_all <- bind_rows(bio_relerr, catch_relerr, rec_relerr, M2_relerr, Z_relerr) %>%
  filter(!is.na(REE)) # F_relerr,

# 2. Compute GOF metrics by simulation
gof_per_sim <- relerr_all %>%
  group_by(sim, metric) %>%
  summarise(
    mean_REE   = mean(REE, na.rm = TRUE),        # bias
    abs_bias   = mean(abs(REE), na.rm = TRUE),   # absolute bias
    sd_REE     = sd(REE, na.rm = TRUE),          # precision
    rmse       = sqrt(mean((REE)^2, na.rm = TRUE)), # RMSE
    .groups = "drop"
  ) %>%
  group_by(sim) %>%
  summarise(
    mean_bias  = mean(mean_REE, na.rm = TRUE),
    mean_absbias = mean(abs_bias, na.rm = TRUE),
    mean_sd     = mean(sd_REE, na.rm = TRUE),
    mean_rmse   = mean(rmse, na.rm = TRUE),
    n_metrics   = n(),
    .groups = "drop"
  ) %>%
  mutate(
    GOF_score = 100 * (1 - mean_rmse / max(mean_rmse, na.rm = TRUE))  # 0–100 scale
  )

# 3. Rank simulations by overall GOF
gof_ranked <- gof_per_sim %>%
  arrange(desc(GOF_score)) %>%
  mutate(rank = row_number())

# 4. Save output tables
write.csv(gof_ranked, "sims/base/plots/overall_GOF_by_sim.csv", row.names = FALSE)

# 5. Optional: simple visualization
library(ggplot2)

ggplot(gof_ranked, aes(x = reorder(sim, -GOF_score), y = GOF_score)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Overall Goodness of Fit per Simulation",
    x = "Simulation",
    y = "Composite GOF Score (0–100)"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("sims/base/plots/overall_GOF_by_sim.png",
       width = 10, height = 5, dpi = 300)


library(ggplot2)

ggplot(par_relerr, aes(x = OM_value, y = value)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~param_base, scales = "free") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Operating Model (OM) vs Estimated Parameters (EM)",
    x = "OM True Value",
    y = "Estimated Value (EM)"
  )

ggsave("sims/base/plots/par_recovery.png",
       width = 10, height = 5, dpi = 300)





# (b) Median REE per parameter
median_err_par <- ggplot(par_summary, aes(x = param_base, y = median_REE)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Median Relative Error (REE) per Estimated Parameter",
    x = "Parameter",
    y = "Median REE"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("sims/low/plots/median_err_params.png",
       plot = median_err_par, width = 10, height = 6, dpi = 300)

# (c) Precision (1/Var) per parameter
precision_par <- ggplot(par_summary, aes(x = param_base, y = precision)) +
  geom_col(fill = "darkgreen") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Precision (1 / Var(REE)) per Estimated Parameter",
    x = "Parameter",
    y = "Precision (1 / Var)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("sims/low/plots/precision_params.png",
       plot = precision_par, width = 10, height = 6, dpi = 300)


# ============================================================
# 5B. SAVE FULL RELATIVE ERRORS TABLE
# ============================================================
par_relerr_out <- par_relerr %>%
  select(sim, param_base, idx, value, OM_value, REE, source) %>%
  arrange(param_base, idx, sim)

# Save full REE dataset
write.csv(par_relerr_out, "sims/low/plots/relative_errors_parameters.csv", row.names = FALSE)


##############################################################################

# COMPARISION

library(dplyr)
library(ggplot2)
library(readr)

# ============================================================
# START: already have these 4 dataframes:
# bio_base_std, bio_low_std, F_base_std, F_low_std
# with columns: year, species, sim, source, predicted, OM_value, REE, metric
# ============================================================

# ---------- BASE ----------
F_base_OM <- RealizedF_all_base %>%
  filter(source == "OM") %>%
  select(species, year, OM_value = Realized_F)

F_base_std <- RealizedF_all_base %>%
  filter(source == "Simulation") %>%
  left_join(F_base_OM, by = c("species", "year")) %>%
  filter(!is.na(OM_value), OM_value > 0) %>%
  mutate(
    REE = (Realized_F - OM_value) / OM_value,
    metric = "Realized F (formula)",
    predicted = Realized_F
  ) %>%
  select(year, species, sim, source, predicted, OM_value, REE, metric)

# ---------- LOW ----------
F_low_OM <- RealizedF_all %>%
  filter(source == "OM") %>%
  select(species, year, OM_value = Realized_F)

F_low_std <- RealizedF_all %>%
  filter(source == "Simulation") %>%
  left_join(F_low_OM, by = c("species", "year")) %>%
  filter(!is.na(OM_value), OM_value > 0) %>%
  mutate(
    REE = (Realized_F - OM_value) / OM_value,
    metric = "Realized F (formula)",
    predicted = Realized_F
  ) %>%
  select(year, species, sim, source, predicted, OM_value, REE, metric)


# 1) Add scenario labels (do NOT touch 'source')
bio_base_std2 <- bio_base_std %>% mutate(scenario = "Base")
bio_low_std2  <- bio_low_std  %>% mutate(scenario = "Low interaction")

F_base_std2   <- F_base_std   %>% mutate(scenario = "Base")
F_low_std2    <- F_low_std    %>% mutate(scenario = "Low interaction")

# 2) Bind scenarios and keep ONLY simulation estimates
bio_std_both <- bind_rows(bio_base_std2, bio_low_std2) %>%
  filter(source == "Simulation")

F_std_both <- bind_rows(F_base_std2, F_low_std2) %>%
  filter(source == "Simulation")

# 3) Summaries pooled across years (and sims)
summarize_relerr_over_years <- function(df) {
  df %>%
    group_by(metric, scenario, species) %>%
    summarise(
      n = sum(!is.na(REE)),
      bias_mean = mean(REE, na.rm = TRUE),
      median_REE = median(REE, na.rm = TRUE),
      q25 = quantile(REE, 0.25, na.rm = TRUE),
      q75 = quantile(REE, 0.75, na.rm = TRUE),
      IQR = IQR(REE, na.rm = TRUE),
      sd_REE = sd(REE, na.rm = TRUE),
      rmse_REE = sqrt(mean((REE)^2, na.rm = TRUE)),
      .groups = "drop"
    )
}

bio_sum_over_years <- summarize_relerr_over_years(bio_std_both)
F_sum_over_years   <- summarize_relerr_over_years(F_std_both)

# 4) Plots
plot_relerr_dist <- function(df, title, ylab) {
  ggplot(df, aes(x = scenario, y = REE)) +
    geom_violin(trim = TRUE, alpha = 0.25) +
    geom_boxplot(width = 0.18, outlier.size = 0.4, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~species, scales = "free_y") +
    theme_minimal() +
    labs(title = title, x = NULL, y = ylab)
}

plot_relerr_summary <- function(sumdf, title, ylab) {
  ggplot(sumdf, aes(x = scenario, y = median_REE)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.15) +
    facet_wrap(~species, scales = "free_y") +
    theme_minimal() +
    labs(title = title, x = NULL, y = ylab,
         caption = "Point = median; bars = IQR (25–75%)")
}

p_bio_dist <- plot_relerr_dist(
  bio_std_both,
  title = "Biomass REE pooled across years: Base vs Low interaction",
  ylab  = "REE(B)"
)

p_bio_sum <- plot_relerr_summary(
  bio_sum_over_years,
  title = "Biomass REE summary over years: median and IQR",
  ylab  = "Median REE(B)"
)

p_F_dist <- plot_relerr_dist(
  F_std_both,
  title = "Realized F REE pooled across years: Base vs Low interaction",
  ylab  = "REE(F_realized)"
)

p_F_sum <- plot_relerr_summary(
  F_sum_over_years,
  title = "Realized F REE summary over years: median and IQR",
  ylab  = "Median REE(F_realized)"
)

# 5) Optional saves
dir.create("sims/compare", recursive = TRUE, showWarnings = FALSE)

ggsave("sims/compare/biomass_REE_dist_base_vs_low.png", p_bio_dist, width = 12, height = 6, dpi = 300)
ggsave("sims/compare/biomass_REE_summary_base_vs_low.png", p_bio_sum, width = 12, height = 6, dpi = 300)
write_csv(bio_sum_over_years, "sims/compare/biomass_REE_summary_over_years.csv")

ggsave("sims/compare/Frealized_REE_dist_base_vs_low.png", p_F_dist, width = 12, height = 6, dpi = 300)
ggsave("sims/compare/Frealized_REE_summary_base_vs_low.png", p_F_sum, width = 12, height = 6, dpi = 300)
write_csv(F_sum_over_years, "sims/compare/Frealized_REE_summary_over_years.csv")

# Print plots
p_bio_dist; p_bio_sum
p_F_dist; p_F_sum


# ============================================================
# 1) Add scenario label + keep only Simulation rows
# ============================================================
M2_base_std <- M2_relerr_base %>%
  mutate(
    scenario = "Base",
    metric = "M2",
    source = trimws(as.character(source))
  ) %>%
  filter(tolower(source) == "simulation")

M2_low_std <- M2_relerr_low %>%
  mutate(
    scenario = "Low interaction",
    metric = "M2",
    source = trimws(as.character(source))
  ) %>%
  filter(tolower(source) == "simulation")

M2_std_both <- bind_rows(M2_base_std, M2_low_std)

# Optional: set scenario order for nicer plots
M2_std_both$scenario <- factor(M2_std_both$scenario, levels = c("Base", "Low interaction"))

# ============================================================
# 2) Summaries pooled across YEARS (and sims)
# ============================================================
M2_sum_over_years <- M2_std_both %>%
  group_by(metric, scenario, species) %>%
  summarise(
    n = sum(!is.na(REE)),
    bias_mean  = mean(REE, na.rm = TRUE),
    median_REE = median(REE, na.rm = TRUE),
    q25 = quantile(REE, 0.25, na.rm = TRUE),
    q75 = quantile(REE, 0.75, na.rm = TRUE),
    IQR = IQR(REE, na.rm = TRUE),
    sd_REE = sd(REE, na.rm = TRUE),
    rmse_REE = sqrt(mean((REE)^2, na.rm = TRUE)),
    .groups = "drop"
  )

# ============================================================
# 3) Plots: distribution + median/IQR summary (same style as B & F)
# ============================================================

# (a) Distribution plot (violin + box) pooled across years
p_M2_dist <- ggplot(M2_std_both, aes(x = scenario, y = REE)) +
  geom_violin(trim = TRUE, alpha = 0.25) +
  geom_boxplot(width = 0.18, outlier.size = 0.4, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "M2 relative error pooled across years: Base vs Low interaction",
    x = NULL,
    y = "REE(M2)"
  )

# (b) Main slide plot: Median ± IQR summary over years
p_M2_sum <- ggplot(M2_sum_over_years, aes(x = scenario, y = median_REE)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.15) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "M2 REE summary over years: median and IQR",
    subtitle = "Simulation estimates only",
    x = NULL,
    y = "Median REE(M2)"
  )

p_M2_dist
p_M2_sum

# ============================================================
# 4) Save (match your compare folder pattern)
# ============================================================
dir.create("sims/compare", recursive = TRUE, showWarnings = FALSE)

ggsave("sims/compare/M2_REE_dist_base_vs_low.png",
       plot = p_M2_dist, width = 12, height = 6, dpi = 300)

ggsave("sims/compare/M2_REE_summary_base_vs_low.png",
       plot = p_M2_sum, width = 12, height = 6, dpi = 300)

write_csv(M2_sum_over_years, "sims/compare/M2_REE_summary_over_years.csv")







































#######################################
###### F
###################################

# ============================================================
# 1. OM F (from OM)
# ============================================================
est_F_OM <- output$Fyr %>%
  as.data.frame() %>%
  pivot_longer(cols = 3:ncol(.), names_to = "year", names_prefix = "V") %>%
  rename(species = "V1", fleet = "V2") %>%
  mutate(
    year = as.numeric(year) - 2,
    species = hydraDataList$speciesList[species],
    sim = "OM",
    source = "OM",
    pred_F = value
  ) %>%
  # Keep only target species–fleet combinations
  filter(
    (species == "Atlantic_cod"       & fleet == 1) |
      (species == "Atlantic_herring"   & fleet == 2) |
      (species == "Atlantic_mackerel"  & fleet == 2) |
      (species == "Spiny_dogfish"      & fleet == 1)
  ) %>%
  select(fleet, year, species, pred_F, sim, source)

# ============================================================
# 2. Simulation F (from 74 simulation fits)
# ============================================================

F_fits_all <- map_dfr(1:nsims, function(i) {
  sim_outputs[[i]]$Fyr %>%
    as.data.frame() %>%
    pivot_longer(cols = 3:ncol(.), names_to = "year", names_prefix = "V") %>%
    rename(species = "V1", fleet = "V2") %>%
    mutate(
      year = as.numeric(year) - 2,
      species = hydraDataList$speciesList[species],
      sim = paste0("sim", i),
      source = "Simulation",
      pred_F = value
    ) %>%
    # keep only the target (species, fleet) combinations
    filter(
      (species == "Atlantic_cod"       & fleet == 1) |
        (species == "Atlantic_herring"   & fleet == 2) |
        (species == "Atlantic_mackerel"  & fleet == 2) |
        (species == "Spiny_dogfish"      & fleet == 1)
    ) %>%
    select(fleet, year, species, pred_F, sim, source)
})

# ============================================================
# 3. Combine OM + Simulation
# ============================================================
F_all <- bind_rows(est_F_OM, F_fits_all)


sim_OM_F <- ggplot(F_all, aes(x = year, y = pred_F)) +
  geom_line(
    data = filter(F_all, source == "Simulation"),
    aes(group = interaction(sim, species)),
    color = "firebrick2", alpha = 0.25
  ) +
  geom_line(
    data = filter(F_all, source == "OM"),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "Fishing mortality (F)",
    title = "Time series of OM fishing mortality (black) and simulation fits (red)"
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

sim_OM_F

ggsave("sims/low/plots/sim_OM_F.png",
       plot = sim_OM_F, width = 10, height = 6, dpi = 300)


om_F_only <- ggplot(filter(F_all, source == "OM"),
                    aes(x = year, y = pred_F)) +
  geom_line(
    aes(group = interaction(sim, species)),
    color = "firebrick2", alpha = 0.25
  ) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "Fishing mortality (F)",
    title = "Simulation fishing mortality (74 replicates)"
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Save the plot
ggsave("sims/low/plots/om_F_only.png",
       plot = om_F_only, width = 10, height = 6, dpi = 300)



# ============================================================
# 5. Compute Relative Errors and Summaries
# ============================================================

F_OM <- F_all %>%
  filter(source == "OM") %>%
  select(species, year, OM_F = pred_F)

F_sim <- F_all %>%
  filter(source == "Simulation")

F_relerr <- F_sim %>%
  left_join(F_OM, by = c("species", "year")) %>%
  mutate(REE = (pred_F - OM_F) / OM_F)

F_summary <- F_relerr %>%
  group_by(species, year) %>%
  summarise(
    median_REE = median(REE, na.rm = TRUE),
    precision  = 1 / var(REE, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 6. Plot Relative Error Diagnostics
# ============================================================

# (a) Boxplots
rel_err_F <- ggplot(F_relerr, aes(x = factor(year), y = REE)) +
  geom_boxplot(fill = "firebrick2", alpha = 0.5, outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Relative Errors in Fishing Mortality Across Simulations",
    x = "Year",
    y = "Relative Error (Simulation vs OM)"
  )

ggsave("sims/low/plots/rel_err_F.png",
       plot = rel_err_F, width = 10, height = 6, dpi = 300)

# (b) Median REE
median_err_F <- ggplot(F_summary, aes(x = year, y = median_REE)) +
  geom_line(color = "red3") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Median Relative Error in Fishing Mortality", y = "Median REE")

ggsave("sims/low/plots/median_err_F.png",
       plot = median_err_F, width = 10, height = 6, dpi = 300)

# (c) Precision (1/Var)
precision_F <- ggplot(F_summary, aes(x = year, y = precision)) +
  geom_line(color = "darkred") +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  labs(title = "Precision of Relative Errors in Fishing Mortality", y = "1 / Var(REE)")

ggsave("sims/low/plots/precision_F.png",
       plot = precision_F, width = 10, height = 6, dpi = 300)


# ============================================================
# 7. Histograms by Species
# ============================================================

# Atlantic cod
ggplot(filter(F_relerr, species == "Atlantic_cod"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "firebrick2", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Fishing Mortality – Atlantic cod",
    x = "Relative Error (REE)",
    y = "Count"
  )

# Atlantic herring
ggplot(filter(F_relerr, species == "Atlantic_herring"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "firebrick2", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Fishing Mortality – Atlantic herring",
    x = "Relative Error (REE)",
    y = "Count"
  )

# Atlantic mackerel
ggplot(filter(F_relerr, species == "Atlantic_mackerel"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "firebrick2", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Fishing Mortality – Atlantic mackerel",
    x = "Relative Error (REE)",
    y = "Count"
  )

# Spiny dogfish
ggplot(filter(F_relerr, species == "Spiny_dogfish"), aes(x = REE)) +
  geom_histogram(bins = 30, fill = "firebrick2", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Relative Errors in Fishing Mortality – Spiny dogfish",
    x = "Relative Error (REE)",
    y = "Count"
  )


