source("R/read.report.R")
source("R/gettables.R")
library(tidyverse)
#library(magrittr)
#library(FSA)
library(hydradata)


### Read observed and estimated values, Hydra data list from Sarahs 4 species scenario
hydraDataList <- readRDS("inputs/hydra_sim_GBself_5bin.rds")
#rep file model (check latest version) 
repfile <- "inputs/initial_run/hydra_sim.rep"
output<-reptoRlist(repfile)


#### READ CATCH AND SURVEY OBSERVED BIOMASS ####
obs_surveyB <- hydraDataList$observedBiomass
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


##### SIMULATE CATCH DATA ######
# create simulation object
sim_data <- NULL
isim <- 1

for (isim in 1:100) { # 5 just to try 
  
  # replace index with simulated data
  
  
  hydraDataList$observedCatch <- obs_catchB %>%
    mutate(catch = rnorm(nrow(.), indexfits[[2]]$pred_catch,indexfits[[2]]$cv))
  
  sim_data[[isim]] <- hydraDataList
  
  
  ##### SIMULATE SURVEY BIOMASS DATA ######
  hydraDataList$observedBiomass <- obs_surveyB %>%
    mutate(biomass = (rnorm(nrow(.), indexfits[[1]]$pred_bio, indexfits[[1]]$cv)))
  # store simulated object
  sim_data[[isim]] <- hydraDataList
  
  
  #### SIMULATE CATCH SIZE COMPOSITION DATA ####
  
  obscatch_size <- hydraDataList$observedCatchSize %>% tibble()
  obscatch_size<-obscatch_size %>% pivot_longer(cols=7:ncol(.), names_to = "lenbin") %>%
    mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")),
           label = rep("catch",nrow(.)))#,
  #species = hydraDataList$speciesList[species])# %>% filter(value != -999)
  #obscatch_size<- obscatch_size %>% filter(value != -999)
  obscatch_size$value[which(obscatch_size$value == -999)] = 0.000001
  obscatch_size <- select(obscatch_size, -label)
  
  pred_catchsize<-output$pred_catch_size
  nll_catch<-output$nll_catch_size
  
  #hydraDataList$observedCatchSize <- obscatch_size %>% 
  #  mutate(value = (rmultinom(1, obscatch_size$inpN, pred_catchsize)))
  # store simulated object
  #sim_data[[isim]] <- hydraDataList
  
  temporal1 = numeric()
  especie = numeric(); especie = sort(unique(obscatch_size$species)) # especies
  for(e in 1:length(especie)){
    pos1 = numeric(); pos1 = which(obscatch_size$species == especie[e])
    year = numeric(); year = sort(unique(obscatch_size$year[pos1]))
    for(y in 1:length(year)){
      pos2 = numeric(); pos2 = which(obscatch_size$year[pos1]== year[y])
      
      temp = numeric(); temp = rmultinom(1, unique(obscatch_size$inpN[pos1][pos2]), pred_catchsize[pos1][pos2])
      temporal1 = c(temporal1, temp)
    }
  }
  
  obscatch_size$value = temporal1
  hydraDataList$observedCatchSize<-obscatch_size
  
  hydraDataList$observedCatchSize["value"]<-hydraDataList$observedCatchSize["value"]/hydraDataList$observedCatchSize["inpN"]
  
  hydraDataList$observedCatchSize<-hydraDataList$observedCatchSize %>% pivot_wider(names_from = "lenbin") %>%
    rename(sizebin1 = `1`, sizebin2 = `2`, sizebin3 = `3`, sizebin4 = `4`, sizebin5 = `5`)
  
  
  #### SIMULATE SURVEY SIZE COMPOSITION DATA ####
  
  obssurv_size <- hydraDataList$observedSurvSize %>% tibble()
  obssurv_size <- obssurv_size %>% pivot_longer(cols=6:ncol(.), names_to = "lenbin") %>% #filter(value != -999)%>%
    
    mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")),
           label = rep("survey",nrow(.)))#,
  #species = hydraDataList$speciesList[species])
  #surv_size<- surv_size %>% filter(value != -999)
  obssurv_size$value[which(obssurv_size$value == -999)] = 0.000001
  obssurv_size <- select(obssurv_size, -label)
  
  pred_survsize<-output$pred_survey_size
  nll_survey<-output$nll_survey_size
  
  #hydraDataList$observedSurvSize <- obssurv_size %>% 
  #  mutate(value = (rmultinom(1, obssurv_size$inpN, pred_survsize)))
  # store simulated object
  #sim_data[[isim]] <- hydraDataList
  
  temporal1 = numeric()
  number = numeric(); number = sort(unique(obssurv_size$survey))
  for (n in 1:length(number)) {
    pos0 = numeric(); pos0 = which(obssurv_size$survey == number[n])
    
    especie = numeric(); especie = sort(unique(obssurv_size$species[pos0])) # especies
    for(e in 1:length(especie)){
      pos1 = numeric(); pos1 = which(obssurv_size$species[pos0] == especie[e])
      
      year = numeric(); year = sort(unique(obssurv_size$year[pos0][pos1]))
      for(y in 1:length(year)){
        pos2 = numeric(); pos2 = which(obssurv_size$year[pos0][pos1]== year[y])
        
        temp = numeric(); temp = rmultinom(1, unique(obssurv_size$inpN[pos0][pos1][pos2]), pred_survsize[pos0][pos1][pos2])
        temporal1 = c(temporal1, temp)
      }
    }
  }
  
  obssurv_size$value = temporal1
  hydraDataList$observedSurvSize<-obssurv_size
  
  hydraDataList$observedSurvSize["value"]<-hydraDataList$observedSurvSize["value"]/hydraDataList$observedSurvSize["inpN"]
  
  hydraDataList$observedSurvSize<-hydraDataList$observedSurvSize %>% pivot_wider(names_from = "lenbin") %>%
    rename(sizebin1 = `1`, sizebin2 = `2`, sizebin3 = `3`, sizebin4 = `4`, sizebin5 = `5`)
  
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
  number = numeric(); number = sort(unique(obsdiet_comp$survey))
  for (n in 1:length(number)) {
    pos0 = numeric(); pos0 = which(obsdiet_comp$survey == number[n])
    
    species = numeric(); species = sort(unique(obsdiet_comp$species[pos0])) # especies
    for(e in 1:length(species)){
      pos1 = numeric(); pos1 = which(obsdiet_comp$species[pos0] == species[e])
      
      year = numeric(); year = sort(unique(obsdiet_comp$year[pos0][pos1]))
      for(y in 1:length(year)){
        pos2 = numeric(); pos2 = which(obsdiet_comp$year[pos0][pos1]== year[y])
        
        lenbin = numeric(); lenbin = sort(unique(obsdiet_comp$sizebin[pos0][pos1][pos2]))
        for(l in 1:length(lenbin)){
          pos3 = numeric(); pos3 = which(obsdiet_comp$sizebin[pos0][pos1][pos2] == lenbin[l])
          
          temp = numeric(); temp = rmultinom(1, unique(obsdiet_comp$inpN[pos0][pos1][pos2][pos3]), obsdiet_comp$value[pos0][pos1][pos2][pos3])
          temporal1 = c(temporal1, temp)
        }
      }
    }
  }
  
  # replace data with simulated data
  
  obsdiet_comp$value = temporal1
  hydraDataList$observedSurvDiet<-obsdiet_comp
  hydraDataList$observedSurvDiet["value"]<-hydraDataList$observedSurvDiet["value"]/hydraDataList$observedSurvDiet["inpN"]
  
  hydraDataList$observedSurvDiet<-hydraDataList$observedSurvDiet %>% pivot_wider(names_from = "prey")
  
  #to check if I am getting the correct values
  #write.csv(diet_comp, file = "diet_comp.csv", row.names = T)
  
  # use write_tsDatFile function to save the new ts file with simulated data 
  #n sim equal to scenario identifier from 1 to 30 simulations 
  
}

#### DIAGNOSTICS ####

source("R/read.report.R")
source("R/gettables.R")
library(ggforce)
library(tidyverse)

hydraDataList <- readRDS("inputs/hydra_sim_GBself_5bin.rds")
#sim_data

write_rds(sim_data, "sim_data.rds")
hydraDataList2 <- readRDS("sim_data.rds")

## OLD WAY TO READ MY SIM DATA
repfiles <- c("simulated data/sim1/hydra_sim.rep","simulated data/sim2/hydra_sim.rep",
              "simulated data/sim3/hydra_sim.rep","simulated data/sim4/hydra_sim.rep" ,"simulated data/sim5/hydra_sim.rep",
              "simulated data/sim6/hydra_sim.rep", "simulated data/sim7/hydra_sim.rep" ,"simulated data/sim8/hydra_sim.rep",
              "simulated data/sim9/hydra_sim.rep", "simulated data/sim10/hydra_sim.rep" ,"simulated data/sim11/hydra_sim.rep",
              "simulated data/sim12/hydra_sim.rep", "simulated data/sim13/hydra_sim.rep" ,"simulated data/sim14/hydra_sim.rep",
              "simulated data/sim15/hydra_sim.rep", "simulated data/sim16/hydra_sim.rep" ,"simulated data/sim17/hydra_sim.rep",
              "simulated data/sim18/hydra_sim.rep", "simulated data/sim19/hydra_sim.rep" ,"simulated data/sim20/hydra_sim.rep",
              "simulated data/sim21/hydra_sim.rep", "simulated data/sim22/hydra_sim.rep" ,"simulated data/sim23/hydra_sim.rep",
              "simulated data/sim24/hydra_sim.rep", "simulated data/sim25/hydra_sim.rep" ,"simulated data/sim26/hydra_sim.rep",
              "simulated data/sim27/hydra_sim.rep", "simulated data/sim28/hydra_sim.rep" ,"simulated data/sim29/hydra_sim.rep",
              "simulated data/sim30/hydra_sim.rep")

## I need to read the 100 sims in one line and then select one variable 
## for the 100 sims, but I have a data list object with 1oo data sets and xxxx variables on each data set



output_est<-purrr::map(repfiles, reptoRlist)
nmodel <- length(repfiles)

F_true<-output$EstFsize

F_true<- as.data.frame(F_true)  %>%
  rename(sizebin1 = 'V1', sizebin2 = 'V2', sizebin3 = 'V3', sizebin4 = 'V4', sizebin5 = 'V5')

F_true <- F_true %>%  #output$EstBsize %>%
  # rowSums() %>%
  tibble() %>%
  mutate(Ftot = rowSums(across(where(is.numeric)))/5) %>%
  mutate(species = (rep(hydraDataList$speciesList, each = hydraDataList$Nyrs)),
         year  = (rep(1:(hydraDataList$Nyrs),4)),
         #         #year = 0.8 + year / 5,  #5 time steps per year
         log_F = ifelse(Ftot>0,log(Ftot),NA))

FF<-rep(F_true$Ftot,times=30)


est_survey <- map(output_est, "survey") %>%
  #walk(as.data.frame) %>%
  map_dfr(as.data.frame, .id = "model")

write.csv(est_survey, file = "biomass.csv", row.names = T)


est_survey<- est_survey %>%
  rename(n_sim = 'model', survey = 'V1', year = 'V2', species = 'V3', obs_value = 'V4', cv = 'V5', pred_value = 'V6', res = 'V7', nll = 'V8') %>%
  mutate(species = hydraDataList$speciesList[species])

surv1plot<-est_survey %>% filter(survey==1)%>%
  ggplot() +
  aes(x = year, y = (pred_value), col = n_sim) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  #geom_point(aes(x=year, y=log(obs_value)), col = "red")+
  #geom_errorbar(aes(ymin = (log(obs_value)-1.96*cv), ymax = (log(obs_value)+1.96*cv)), col="gray")+
  theme_minimal() +
  labs(x = "Year",
       y = "Biomass (t)",
       title = "Time series of estimated LN(biomass)")

print(surv1plot)


surv2plot<-est_survey %>% filter(survey==2)%>%
  ggplot() +
  aes(x = year, y = log(pred_value), col = n_sim) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  geom_point(aes(x=year, y=log(obs_value)), col = "red")+
  geom_errorbar(aes(ymin = (log(obs_value)-1.96*cv), ymax = (log(obs_value)+1.96*cv)), col="gray")+
  theme_minimal() +
  labs(x = "Year",
       y = "Biomass (t)",
       title = "Time series of estimated LN(biomass)")

print(surv2plot)


est_catch <- map(output_est, "catch") %>%
  #walk(as.data.frame) %>%
  map_dfr(as.data.frame, .id = "model")

est_catch<- est_catch %>%
  rename(n_sim = 'model', fleet = 'V1', area = 'V2', year = 'V3', species = 'V4', obs_value = 'V5', cv = 'V6', pred_value = 'V7', res = 'V8', nll = 'V9') %>%
  mutate(species = hydraDataList$speciesList[species])

fleet1plot<-est_catch %>% filter(fleet==1)%>%
  ggplot() +
  aes(x = year, y = log(pred_value), col = n_sim) +
  geom_line() +
  facet_wrap(~species, scales = "free",dir="v") +
  #geom_line() +
  geom_point(aes(x=year, y=log(obs_value)), col = "red")+
  geom_errorbar(aes(ymin = (log(obs_value)-1.96*cv), ymax = (log(obs_value)+1.96*cv)), col="gray")+
  theme_minimal() +
  labs(x = "Year",
       y = "Catch (t)",
       title = "Time series of estimated LN(catch)")

print(fleet1plot)

fleet2plot<-est_catch %>% filter(fleet==2)%>%
  ggplot() +
  aes(x = year, y = log(pred_value), col = n_sim) +
  geom_line() +
  facet_wrap(~species, scales = "free",dir="v") +
  geom_point(aes(x=year, y=log(obs_value)), col = "red")+
  geom_errorbar(aes(ymin = (log(obs_value)-1.96*cv), ymax = (log(obs_value)+1.96*cv)), col="gray")+
  theme_minimal() +
  labs(x = "Year",
       y = "Catch (t)",
       title = "Time series of estimated LN(catch)")

print(fleet2plot)


est_F <- map(output_est, "EstFsize") %>%
  #walk(as.data.frame) %>%
  map_dfr(as.data.frame, .id = "model")

est_F<- est_F %>%
  rename(model = 'model', sizebin1 = 'V1', sizebin2 = 'V2', sizebin3 = 'V3', sizebin4 = 'V4', sizebin5 = 'V5')

est_F <- est_F %>%  #output$EstBsize %>%
  # rowSums() %>%
  tibble() %>%
  mutate(Ftot = rowSums(across(where(is.numeric)))/5) %>%
  mutate(species = rep(rep(hydraDataList$speciesList, each = hydraDataList$Nyrs),nmodel),
         year  = rep(rep(1:(hydraDataList$Nyrs),4),nmodel),
         #year = 0.8 + year / 5,  #5 time steps per year
         log_F = ifelse(Ftot>0,log(Ftot),NA))
#model = as.factor(model)) %>%


ftot1plot<-est_F %>%
  ggplot() +
  aes(x = year, y = Ftot/mean(Ftot), col = model) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  geom_line() +
  geom_point(aes(x=year, y=FF/mean(FF)), col = "red")+
  #geom_errorbar(aes(ymin = (log(obs_value)-1.96*cv), ymax = (log(obs_value)+1.96*cv)), col="gray")+
  theme_minimal() +
  labs(x = "Year",
       y = "Ftot (year-1)",
       title = "Time series of estimated fishing mortality")

print(ftot1plot)


est_M2 <- map(output_est, "EstM2size") %>%
  #walk(as.data.frame) %>%
  map_dfr(as.data.frame, .id = "model")

est_M2<- est_M2 %>%
  rename(model = 'model', sizebin1 = 'V1', sizebin2 = 'V2', sizebin3 = 'V3', sizebin4 = 'V4', sizebin5 = 'V5')

est_M2 <- est_M2 %>%  #output$EstBsize %>%
  # rowSums() %>%
  tibble() %>%
  mutate(M2 = rowSums(across(where(is.numeric)))) %>%
  mutate(species = rep(rep(hydraDataList$speciesList, each = hydraDataList$Nyrs),nmodel),
         year  = rep(rep(1:(hydraDataList$Nyrs),4),nmodel),
         #year = 0.8 + year / 5,  #5 time steps per year
         log_M2 = ifelse(M2>0,log(M2),NA))
#model = as.factor(model)) %>%

M2_true<-output$EstM2size

M2_true<- as.data.frame(M2_true)  %>%
  rename(sizebin1 = 'V1', sizebin2 = 'V2', sizebin3 = 'V3', sizebin4 = 'V4', sizebin5 = 'V5')

M2_true <- M2_true %>%  #output$EstBsize %>%
  # rowSums() %>%
  tibble() %>%
  mutate(M2 = rowSums(across(where(is.numeric)))) %>%
  mutate(species = (rep(hydraDataList$speciesList, each = hydraDataList$Nyrs)),
         year  = (rep(1:(hydraDataList$Nyrs),4)),
         #         #year = 0.8 + year / 5,  #5 time steps per year
         log_M2 = ifelse(M2>0,log(M2),NA))

pred_mort<-rep(M2_true$M2,times=30)

M2plot<-est_M2 %>% #select(year, species, model, sizebin5, M2) %>%
  ggplot() +
  aes(x = year, y = M2/mean(M2), col = model) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  geom_line() +
  geom_point(aes(x=year, y=pred_mort/mean(pred_mort)), col = "red")+
  theme_minimal() +
  labs(x = "Year",
       y = "M2 (year-1)",
       title = "Time series of estimated natural mortality")

print(M2plot)


survey_obspred<-indexfits[1][[1]]%>%
  mutate(obs = biomass + 1e-07,
         pred = pred_bio + 1e-07,
         log_obs = log(obs),
         log_pred = log(pred),
         log_lo = log_obs - 1.96*cv,
         log_hi = log_obs + 1.96*cv,
         obs_lo = exp(log_lo),
         obs_hi = exp(log_hi),
         species = hydraDataList$speciesList[species])

p1 <- survey_obspred %>%
  filter(survey == 1) %>%
  ggplot() +
  aes(x= year, y = log_obs, group = species, col=factor(survey)) +
  geom_errorbar(aes(ymin = log_lo, ymax = log_hi)) +
  #geom_point() +
  geom_line(aes(x=year, y=log_pred), col = "blue") +
  facet_wrap(~species, scales = "free_y") +
  theme_bw() +
  guides(species = "None")
print(p1)

p2 <- survey_obspred %>%
  filter(survey == 2) %>%
  ggplot() +
  aes(x= year, y = log_obs, group = species, col=factor(survey)) +
  geom_errorbar(aes(ymin = log_lo, ymax = log_hi)) +
  #geom_point() +
  geom_line(aes(x=year, y=log_pred), col = "blue") +
  facet_wrap(~species, scales = "free_y") +
  theme_bw() +
  guides(species = "None")
print(p2)

catch_obspred<-indexfits[2][[1]]%>%
  mutate(obs = catch + 1e-07,
         pred = pred_catch + 1e-07,
         log_obs = log(obs),
         log_pred = log(pred),
         log_lo = log_obs - 1.96*cv,
         log_hi = log_obs + 1.96*cv,
         obs_lo = exp(log_lo),
         obs_hi = exp(log_hi),
         species = hydraDataList$speciesList[species])

p3 <- catch_obspred %>%
  filter(fishery == 1) %>%
  ggplot() +
  aes(x= year, y = log_obs, group = species, col=factor(fishery)) +
  geom_errorbar(aes(ymin = log_lo, ymax = log_hi)) +
  geom_point() +
  geom_line(aes(x=year, y=log_pred), col = "blue") +
  facet_wrap(~species, scales = "free_y") +
  theme_bw() +
  guides(species = "None")
print(p3)

p4 <- catch_obspred %>%
  filter(fishery == 2) %>%
  ggplot() +
  aes(x= year, y = log_obs, group = species, col=factor(fishery)) +
  geom_errorbar(aes(ymin = log_lo, ymax = log_hi)) +
  geom_point() +
  geom_line(aes(x=year, y=log_pred), col = "blue") +
  facet_wrap(~species, scales = "free_y") +
  theme_bw() +
  guides(species = "None")
print(p4)













