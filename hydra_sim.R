source("R/read.report.R")
source("R/gettables.R")
library(tidyverse)
#library(magrittr)
#library(FSA)
library(hydradata)


### Read observed and estimated values, Hydra data list from Sarahs 4 species scenario
hydraDataList <- readRDS("inputs/Sarah_files/hydra_sim_GBself_5bin.rds")
#rep file model (check latest version) 
repfile <- "inputs/mod01_1_survey/hydra_sim.rep"
output_OM<-reptoRlist(repfile)

est_bio_OM <- output_OM[["EstBsize"]] 
est_bio_OM <- est_bio_OM %>%  
  #rowSums() %>% 
  tibble() %>% 
  mutate(bio = rowSums(across(where(is.numeric)))) %>% 
  mutate(species = rep(rep(hydraDataList$speciesList, each = hydraDataList$Nyrs)),
         year  = rep((1:hydraDataList$Nyrs),length(hydraDataList$speciesList)),
         log_bio = ifelse(bio>0,log(bio),NA)) 
  
plot_bio_OM<-est_bio_OM %>% 
              ggplot() +
              geom_line(aes(x = year, y = log(bio))) +
              facet_wrap(~species, scales = "free") +
              theme_minimal() +
              labs(x = "Year",
                   y = "Biomass (t)",
                   title = "Time series of estimated biomass")

#ggsave("biomassOM.jpeg", plot_bio_OM, width=10, height=5, dpi = 300)#, width=3000, height=2500, res=250)


est_F_OM <- output_OM[["EstFsize"]] 
est_F_OM <- est_F_OM %>%  
  #rowSums() %>% 
  tibble() %>% 
  mutate(fmort = rowSums(across(where(is.numeric)))) %>% 
  mutate(species = rep(rep(hydraDataList$speciesList, each = hydraDataList$Nyrs)),
         year  = rep((1:hydraDataList$Nyrs),length(hydraDataList$speciesList)))

plot_F_OM<-est_F_OM %>% 
  ggplot() +
  geom_line(aes(x = year, y = fmort)) +
  facet_wrap(~species, scales = "free") +
  theme_minimal() +
  labs(x = "Year",
       y = "Fishing mortality",
       title = "Time series of estimated fishing mortality")

est_R_OM <- output_OM[["EstRec"]] 
est_R_OM <- est_R_OM %>%  
  tibble () %>% mutate(species = rep(rep(hydraDataList$speciesList, each = hydraDataList$Nyrs-1)),
       year  = rep((1:hydraDataList$Nyrs),length(hydraDataList$speciesList)))

#ggsave("F_OM.jpeg", plot_F_OM, width=10, height=5, dpi = 300)#, width=3000, height=2500, res=250)
#### READ CATCH AND SURVEY OBSERVED BIOMASS ####
# if you have one survey use this line
obs_surveyB <- hydraDataList$observedBiomass %>% 
  filter(survey == 1)
# if you have 2 surveys use this line
#obs_surveyB <- hydraDataList$observedBiomass
obs_catchB <- hydraDataList$observedCatch

biorows <- dim(obs_surveyB)[1]
catchrows <- dim(obs_catchB)[1]

#create a table with estimated and observed values
#indexfits <- gettables(repfile, biorows, catchrows)
  
# select survey biomass

biomass<-indexfits[[1]] %>%
  mutate(species = hydraDataList$speciesList[species])
# select catch
catch<-indexfits[[2]]%>%
  mutate(species = hydraDataList$speciesList[species])


### READ FILES 100 .REP ###
library(here)
library(tidyverse)
hydraDataList2 <- readRDS("sim_data_1survey.rds")
dir<-here()
dir<-paste0(dir,"/","sims","/","1_survey","/","rep")
#dir<-paste0(dir,"/","sims","/","initial","/","rep")
setwd(dir)


filelist = list.files(pattern = ".*.rep")
output2 = purrr::map(filelist, function(x)reptoRlist(x)) 

# check biorows from one survey 
biorows <- dim(obs_surveyB)[1]
catchrows <- dim(obs_catchB)[1]

fit_data <- NULL
#nsim <- 100
#"sims","/","1_survey","/","rep",/,filelist[nsim]
for (nsim in 1:100)
{
  repp<-paste0(filelist[nsim])
  ind <- gettables(repp, biorows, catchrows)
  
  fit_data[nsim] <- ind
  
}

# Extraer las variables de cada dataframe y añadir una columna de simulación
sim_bio_pred <- do.call(rbind, lapply(seq_along(fit_data), function(i) {
  df <- fit_data[[i]]
  df$simulacion <- i
  return(df)
}))

sim_bio_pred$species[sim_bio_pred$species == "1"] <- "Atlantic cod"
sim_bio_pred$species[sim_bio_pred$species == "2"] <- "Atlantic herring"
sim_bio_pred$species[sim_bio_pred$species == "3"] <- "Atlantic mackerel"
sim_bio_pred$species[sim_bio_pred$species == "4"] <- "Spinny dogfish"

surv1plot<-ggplot(sim_bio_pred, aes(x = year, y = log(pred_bio), color = as.factor(simulacion))) +
  geom_line(alpha = 0.5) + # Transparencia para las líneas de simulación
  geom_line(data = pred_bioOM, aes(x = year, y = log(pred_bio)), color = "black", size = 1) +
  facet_wrap(~species, scales = "free") +
      labs(title = "Time series of estimated LN(biomass)", 
       x = "Years", 
       y = "Predicted Biomass", 
       color = "Simulación") +
  theme_minimal()+
  theme(legend.position = "none")

print(surv1plot)


hydraDataList <- readRDS("inputs/Sarah_files/hydra_sim_GBself_5bin.rds")
repfile <- "inputs/mod01_1_survey/hydra_sim.rep"
output_OM<-reptoRlist(repfile)

obs_surveyB <- hydraDataList$observedBiomass %>% 
  filter(survey==1)
obs_catchB <- hydraDataList$observedCatch

biorows <- dim(obs_surveyB)[1]
catchrows <- dim(obs_catchB)[1]

#create a table with estimated and observed values
indexfits <- gettables(repfile, biorows, catchrows)

pred_bioOM<-data.frame(indexfits[1])
pred_bioOM$species[pred_bioOM$species == "1"] <- "Atlantic cod"
pred_bioOM$species[pred_bioOM$species == "2"] <- "Atlantic herring"
pred_bioOM$species[pred_bioOM$species == "3"] <- "Atlantic mackerel"
pred_bioOM$species[pred_bioOM$species == "4"] <- "Spinny dogfish"




































































































est_survey<-fit_data %>% 
  ggplot() +
  aes(x = year, y = log(biomass), col = isim) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "Biomass (t)",
       title = "Time series of estimated LN(biomass)")

print(surv1plot)


surv1plot<-sim_obs_bio %>% filter(survey==1)%>%
  ggplot() +
  aes(x = year, y = log(biomass), col = isim) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "Biomass (t)",
       title = "Time series of estimated LN(biomass)")

print(surv1plot)


est_bio <- map(output, "EstBsize") %>%  
  #walk(as.data.frame) %>% 
  map_dfr(as.data.frame, .id = "model") 

nmodel <- length(output)
est_bio <- est_bio %>%  #output$EstBsize %>%
  #rowSums() %>% 
  tibble() %>% 
  mutate(bio = rowSums(across(where(is.numeric)))) %>% 
  mutate(species = rep(rep(hydraDataList$speciesList, each = hydraDataList$Nyrs),nmodel),
         year  =   rep(rep(1:hydraDataList$Nyrs),(nmodel*length(hydraDataList$speciesList))),
         #year =  year / 5,  #5 time steps per year
         log_bio = ifelse(bio>0,log(bio),NA)) 

#para chequear  
#write.csv(est_bio, file = "biomass.csv", row.names = T)

survey_obspred <- indexfits[1][[1]] %>%
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
      geom_point() +
      geom_line(data= est_bio, aes(x = year, y = log_bio, group= model, col = model), cex=0.2) +
      facet_wrap(~species, scales = "free_y") +
      theme(legend.position = "none") 
      
p1  

p1 <- survey_obspred %>% 
  filter(survey == 1) %>% 
  ggplot() +
  aes(x= year, y = log_pred, group = species, col=factor(survey)) +
  #geom_errorbar(aes(ymin = log_lo, ymax = log_hi)) +
  geom_point() +
  geom_line(data= est_bio, aes(x = year, y = log_bio, group= model, col = model), cex=0.2) +
  facet_wrap(~species, scales = "free_y") +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "log Biomass",
       title = "Time series of estimated LN(biomass)")

p1  

#ggsave("plots/biomass_pred_OM.jpeg", p1, width=10, height=5, dpi = 300)#, width=3000, height=2500, res=250)

plot_bio<-est_bio %>% 
          ggplot() +
          aes(x = year, y = bio, col = model) +
          #geom_point() +x
          geom_line() +
          facet_wrap(~species, scales = "free") +
          theme(legend.position = "none") +
          labs(x = "Year",
          y = "Biomass (t)",
          title = "Time series of estimated biomass")
plot_bio


plot_logBio<- est_bio %>% 
              ggplot() +
              aes(x = year, y = log_bio, col = model) +
              geom_line() +
              facet_wrap(~species, scales = "free") +
              theme(legend.position = "none") +
              labs(x = "Year",
              y = "log Biomass (t)",
              title = "Time series of estimated LN(biomass)")
plot_logBio


#ggsave("plots/biomass.jpeg", plot_logBio, width=10, height=5, dpi = 300)#, width=3000, height=2500, res=250)


#### PLOT 100 FITS CATCH ####

est_catch <- map(output, "catch") %>%  
  #walk(as.data.frame) %>% 
  map_dfr(as.data.frame, .id = "model") 

nmodel <- length(output)
est_catch <- est_catch %>%  #output$EstBsize %>%
  #rowSums() %>% 
  tibble() %>% 
  	     rename(fleet = "V1",
         area = "V2",
         year = "V3",
         species = "V4",
         observed = "V5",
         cv = "V6",
         pred = "V7",
         residual = "V8",
         nll = "V9") %>% 
  mutate(species = hydraDataList$speciesList[species])

#para chequear  
#write.csv(est_catch, file = "est_catch.csv", row.names = T)

catch_obspred <- indexfits[2][[1]] %>% 
  mutate(obs = catch + 1e-07,
         pred = pred_catch + 1e-07,
         log_obs = log(obs),
         log_pred = log(pred),
         log_lo = log_obs - 1.96*cv,
         log_hi = log_obs + 1.96*cv,
         obs_lo = exp(log_lo),
         obs_hi = exp(log_hi),
         species = hydraDataList$speciesList[species])

plot_pred_catch<- est_catch %>% 
  ggplot() +
  aes(x = year, y = pred, col = model) +
  geom_line() +
  geom_point(aes(x = year, y = observed, col = model), cex=0.5) +
  facet_wrap(~species, scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "Pred catch (t)",
       title = "Time series of estimated catch")

plot_pred_catch

#ggsave("plots/pred_catch.jpeg", plot_pred_catch, width=10, height=5, dpi = 300)#, width=3000, height=2500, res=250)

p100 <- catch_obspred %>% 
          ggplot() +
          aes(x= year, y = log_pred, group = species, col="red") +
          #geom_errorbar(aes(ymin = log_lo, ymax = log_hi)) +
          geom_point() +
          geom_line(data= est_catch, aes(x = year, y = pred, group =model, col = model), cex=0.5) +
          facet_wrap(~species, scales = "free_y") +
          theme(legend.position = "none") +
          labs(x = "Year",
          y = "Pred log catch",
          title = "Time series of log estimated catch")


p100

#ggsave("plots/pred_catch2.jpeg", p100, width=10, height=5, dpi = 300)#, width=3000, height=2500, res=250)

#### PLOT 100 FITS F ####

est_F <- map(output, "Fyr") %>%  
         map_dfr(as.data.frame, .id = "model")  
        
est_F <- est_F %>% 
  as.data.frame() %>% 
  pivot_longer(cols=4:ncol(.), names_to = "year", names_prefix = "V") %>% 
  rename(species = "V1",
         fleet = "V2") %>% 
  mutate(year = as.numeric(year)-2,
         species = hydraDataList$speciesList[species])
  

#para chequear  
#write.csv(est_F, file = "F.csv", row.names = T)

plot_F1<- est_F %>% filter(fleet == 1) %>%
          ggplot() +
          aes(x = year, y = value, col = model) +
          geom_line() +
          facet_wrap(~species, scales = "free") +
          theme(legend.position = "none") +
          labs(x = "Year",
          y = "F",
          col = "model",
          title = "Time series of estimated fishing mortality")

plot_F1


plot_F2<- est_F %>% filter(fleet == 2) %>%
  ggplot() +
  aes(x = year, y = value, col = model) +
  #geom_point() +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = "F",
       col = "model",
       title = "Time series of estimated fishing mortality")

plot_F2

#ggsave("plots/F.jpeg", plot_F1, dpi = 300)#, width=3000, height=2500, res=250)

#### PLOT 100 FITS M2 ####

est_M2 <- map(output, "EstM2size") %>%  
  map_dfr(as.data.frame, .id = "model")  

est_M2 <- est_M2 %>%  
  tibble() %>% 
  mutate(M2 = rowSums(across(where(is.numeric)))) %>% 
  mutate(species = rep(rep(hydraDataList$speciesList, each = hydraDataList$Nyrs),nmodel),
         year  =   rep(rep(1:hydraDataList$Nyrs),(nmodel*length(hydraDataList$speciesList))))

#para chequear  
#write.csv(est_M2, file = "M2.csv", row.names = T)

plot_M2<- est_M2 %>% 
          ggplot() +
          aes(x = year, y = M2, col = model) +
          geom_line() +
          facet_wrap(~species, scales = "free") +
          theme(legend.position = "none") +
          labs(x = "Year",
          y = "M2",
          col = "model",
          title = "Time series of estimated predation mortality") 

plot_M2  

#ggsave("plots/M2.jpeg", plot_M2, dpi = 300)#, width=3000, height=2500, res=250)

#### PLOT 100 FITS RECRUITS ####

est_Rec <- map(output, "EstRec") %>%  
  map_dfr(as.data.frame, .id = "model")  

est_Rec <- est_Rec %>% 
  as.data.frame() %>% 
  pivot_longer(cols=2:ncol(.), names_to = "year", names_prefix = "V") %>% 
  mutate(species = rep(rep(hydraDataList$speciesList, each = hydraDataList$Nyrs-1),nmodel),
         log_rec = ifelse(value > 0,log(value),NA))
         
#para chequear  
#write.csv(est_Rec, file = "est_Rec.csv", row.names = T)
level_order <- seq(1,41,1) 
plot_Rec<-est_Rec %>% 
    ggplot() +
    geom_line(aes(x = year, y = value, col= model, group= model)) +
    scale_x_discrete(limits= level_order)+
    facet_wrap(~species, scales = "free") +
    theme(legend.position = "none") +
    labs(x = "Year",
         y = "Recruitment (thousands)",
         col= "model",
         title = "Time series of estimated recruitment")

plot_Rec


plot_logRec<-est_Rec %>% 
  ggplot() +
  geom_line(aes(x = year, y = log_rec, col= model, group= model), size=0.5) +
  scale_x_discrete(limits= level_order)+
  facet_wrap(~species, scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Year",
       y = " log Recruitment",
       col= "model",
       title = "Time series of estimated log recruitment")

plot_logRec

#ggsave("plots/plot_logRec.jpeg", plot_logRec,  width=10, height=5, dpi = 300)#, width=3000, height=2500, res=250)






stepperyr <- output$Nstepsyr
if (length(stepperyr)==0) stepperyr <- nrow(output$EstBsize)/hydraDataList$Nyrs/length(hydraDataList$speciesList)

est_bio <- output$EstBsize %>%
  rowSums() %>% 
  as.data.frame() %>% 
  rename(bio = ".") %>% 
  mutate(species = rep(hydraDataList$speciesList, each = hydraDataList$Nyrs*stepperyr),
         year  = rep(1:(hydraDataList$Nyrs*stepperyr),hydraDataList$Nspecies),
         year = (1-1/stepperyr) + year / stepperyr,  #5 time steps per year
         log_bio = ifelse(bio>0,log(bio),NA)) %>%
  I()

if(skill){
  truebio <- omlist_ss$truetotbio_ss %>%
    dplyr::filter(time %in% fittimes.days) %>%
    dplyr::mutate(year = (time/365),
                  year = year-fitstartyr,
                  biomass = atoutput,
                  logbio = log(biomass)) %>%
    I()
  
}
#est_bio

plotB <-ggplot() +
  geom_line(data=truebio, aes(x=year,y=biomass, color="True B"), 
            alpha = 10/10) +
  geom_point(data=est_bio, aes(x=year,y=bio, color="Hydra Est B"),
             alpha = 10/10) + 
  theme_minimal() +
  theme(legend.position = "top") +
  labs(x = "Year",
       y = "Biomass (t)",
       title = "Total biomass skill, Hydra estimated vs Atlantis",
       colour="")

plotB +
  facet_wrap(~species, scales="free") 


plotlogB <-ggplot() +
  geom_line(data=truebio, aes(x=year,y=logbio, color="True log B"), 
            alpha = 10/10) +
  geom_point(data=est_bio, aes(x=year,y=log_bio, color="Hydra Est log B"),
             alpha = 10/10) + 
  theme_minimal() +
  theme(legend.position = "top") +
  labs(x = "Year",
       y = "Biomass (t)",
       title = "Total biomass skill (log scale), Hydra estimated vs Atlantis",
       colour="")

plotlogB +
  facet_wrap(~species, scales="free") 















survey_sel <- map(output, "survey_sel") %>%  
  map_dfr(as.data.frame, .id = "model")  

survey_sel <- survey_sel %>% 
  as.data.frame() %>% 
  pivot_longer(cols=4:ncol(.), names_to = "lenbin", names_prefix = "V") %>% 
  rename(species = "V1",
         survey = "V2") %>% 
  mutate(year = as.numeric(year)-2,
         species = hydraDataList$speciesList[species])

#para chequear  
write.csv(survey_sel, file = "survey_sel.csv", row.names = T)


survey_sel <- output$survey_sel %>% 
  as.data.frame() %>% 
  pivot_longer(cols = 3:ncol(.), names_prefix = "V", names_to = "ilen") %>% 
  rename(species = "V1",
         survey = "V2") %>% 
  mutate(ilen = as.numeric(ilen) - 2,
         species = hydraDataList$speciesList[species])  
  

survey_sel %>% filter(model==1)
  ggplot() +
  aes(x=ilen, y = value, col = factor(survey)) +
  geom_line() +
  facet_wrap(~species) +
  labs(x = "length bin",
       y = "selectivity",
       col = "survey") +
 

survey_sel <- output$fishsel %>% 
  as.data.frame() %>% 
  pivot_longer(cols = 3:ncol(.), names_prefix = "V", names_to = "ilen") %>% 
  rename(species = "V1",
         fleet = "V2") %>% 
  mutate(ilen = as.numeric(ilen) - 2,
         species = hydraDataList$speciesList[species]) %>% 
  I()

survey_sel %>% 
  ggplot() +
  aes(x=ilen, y = value, col = factor(fleet)) +
  geom_line() +
  facet_wrap(~species) +
  labs(x = "length bin",
       y = "selectivity",
       col = "fleet") +
  NULL

