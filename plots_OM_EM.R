#### RUN THE MODEL WITH 100 SIMS ####
library(here)
dir<-here()
dir<-paste0(dir,"/","sims","/","initial")
setwd(dir)

#multiple calls to 'system()' given different folders/filenames.

#isim<-1
#system("cp sims/hydra_sim1-ts.dat sims/hydra_sim_GBself_5bin-ts.dat")
#file.copy(from="sims/hydra_sim1-ts", to="sims/hydra_sim_GBself_5bin-ts")

for (nsim in 1:100)
{
  file.copy(from=paste0("hydra_sim",nsim,"-ts.dat"), to= "hydra_GBself_5bin_simdata-ts.dat", overwrite = TRUE)
  system("./hydra_sim -ind hydra_sim_GBself_5bin.dat -ainp hydra_sim_GBself_5bin.pin")
  file.copy(from = "hydra_sim.rep", to = paste0("rep/hydra_sim",nsim,".rep"))
  file.copy(from = "hydra_sim.par", to = paste0("par/hydra_sim",nsim,".par"))
}



#### DIAGNOSTICS ####
#browseVignettes("hydradata")

source("R/read.report.R")
source("R/gettables.R")

library(ggforce)
library(tidyverse)
hydraDataList2 <- readRDS("sim_data.rds")

#### PLOT SIM CATCH ####

sim_obs_catch<-purrr::map_dfr(hydraDataList2,"observedCatch",.id = "isim")
sim_obs_catch<- sim_obs_catch %>%
  mutate(species = hydraDataList$speciesList[species])

fleet1plot<-sim_obs_catch %>% filter(fishery==1)%>%
  ggplot() +
  aes(x = year, y = log(catch), col = isim) +
  geom_line() +
  facet_wrap(~species, scales = "free",dir="v") +
  #theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "Catch (t)",
       title = "Time series of estimated LN(catch)")

print(fleet1plot)

fleet2plot<-sim_obs_catch %>% filter(fishery==2)%>%
  ggplot() +
  aes(x = year, y = (catch), col = isim) +
  geom_line() +
  facet_wrap(~species, scales = "free",dir="v") +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "Catch (t)",
       title = "Time series of estimated LN(catch)")

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
       title = "Time series of estimated LN(biomass)")

print(surv1plot)


surv2plot<-sim_obs_bio %>% filter(survey==2)%>%
  ggplot() +
  aes(x = year, y = (biomass), col = isim) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "Biomass (t)",
       title = "Time series of estimated LN(biomass)")

print(surv2plot)

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

sp<-1
plot_surv <- list()
especies<-unique(sim_surv_lenght$species)
for (sp in especies) {
  
  temp_size<-sim_surv_lenght %>% filter(species == sp & survey==2) %>%
    group_by(year) %>%
    summarize(mu_ss=mean(inpN))
  
  plot_surv[[sp]] <- sim_surv_lenght %>% filter (species==sp & survey==2) %>%
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
plot_catch$Atlantic_herring
plot_catch$Atlantic_mackerel
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

plot_catch$Atlantic_cod
plot_catch$Atlantic_herring
plot_catch$Atlantic_mackerel
plot_catch$Spiny_dogfish


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




