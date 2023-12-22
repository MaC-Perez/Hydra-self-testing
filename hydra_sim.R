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
              geom_line(aes(x = year, y = bio)) +
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
  #rowSums() %>% 
  tibble() %>% 
  mutate(species = rep(rep(hydraDataList$speciesList, each = hydraDataList$Nyrs)),
         year  = rep((1:hydraDataList$Nyrs),length(hydraDataList$speciesList)))

#ggsave("F_OM.jpeg", plot_F_OM, width=10, height=5, dpi = 300)#, width=3000, height=2500, res=250)
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


### READ FILES 100 .REP ###
library(here)
library(tidyverse)
dir<-here()
dir<-paste0(dir,"/","sims","/","rep")
setwd(dir)

filelist = list.files(pattern = ".*.rep")
output = purrr::map(filelist, function(x)reptoRlist(x)) 

#### PLOT 100 FITS BIOMASS ####

setwd("C:/Users/macristina.perez/Documents/GitHub/Hydra-self-testing")

fit_data <- NULL
#nsim <- 1

for (nsim in 1:100)
{
  repp<-paste0("sims","/","rep","/", filelist[nsim])
  ind <- gettables(repp, biorows, catchrows)
  
  fit_data[nsim] <- ind
  
}


for (nsim in 1:100) {
  bioma_pred<- fit_data[[nsim]] %>% 
  mutate(species = hydraDataList$speciesList[species])
}

#fit_data<-purrr::map_dfr(fit_data,nsim,.id = "isim")
#fit_data

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



















































































temporal2 = numeric()
especie = numeric(); especie = sort(unique(temp.catch$species)) # especies
for(e in 1:length(especie)){
  pos1 = numeric(); pos1 = which(temp.catch$species == especie[e])
  
  temporal1 = numeric()
  year = numeric(); year = sort(unique(temp.catch$year[pos1])) # ano
  for(y in 1:length(year)){
    pos2 = numeric(); pos2 = which(temp.catch$year[pos1] == year[y])
    
    N = numeric(); N = unique(temp.catch$InpN[pos1][pos2]) # sample size para el estrato especie ano
    n_ejemplares_obs = numeric(); n_ejemplares_obs = round(N * temp.catch$obs_value[pos1][pos2] , 0) # crea vector de proporcion por el sample size
    n_ejemplares_est = numeric(); n_ejemplares_est = round(N * temp.catch$pred_value[pos1][pos2], 0)
    
    tt1 = numeric()
    tt1 = data.frame(YEAR = year[y], BIN = temp.catch$lenbin[pos1][pos2], N = N, N_EJEM_OBS = n_ejemplares_obs, N_EJEM_EST = n_ejemplares_est)
    
    temporal1 = rbind(temporal1, tt1) # guarda la proporcion * sample size para la especie seleccionada y ano
  }
  
  bin1_obs = numeric(); bin1_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 1)])
  bin2_obs = numeric(); bin2_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 2)])
  bin3_obs = numeric(); bin3_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 3)])
  bin4_obs = numeric(); bin4_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 4)])
  bin5_obs = numeric(); bin5_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 5)])
  
  total_obs = sum(bin1_obs, bin2_obs, bin3_obs, bin4_obs, bin5_obs)
  prop_new_obs = numeric(); prop_new_obs = c(bin1_obs, bin2_obs, bin3_obs, bin4_obs, bin5_obs) / total_obs
  
  bin1_pred = numeric(); bin1_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 1)])
  bin2_pred = numeric(); bin2_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 2)])
  bin3_pred = numeric(); bin3_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 3)])
  bin4_pred = numeric(); bin4_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 4)])
  bin5_pred = numeric(); bin5_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 5)])
  
  total_pred = sum(bin1_pred, bin2_pred, bin3_pred, bin4_pred, bin5_pred)
  prop_new_est = numeric(); prop_new_est = c(bin1_pred, bin2_pred, bin3_pred, bin4_pred, bin5_pred) / total_pred
  
  tt2 = numeric()
  tt2 = data.frame(ESPECIE = especie[e], BIN = seq(1, 5, 1), PROP_OBS = prop_new_obs, PROP_EST = prop_new_est)
  
  temporal2 = rbind(temporal2, tt2) # almacena por cada especie
}


long_temp2 <- temporal2 %>%
  pivot_longer(cols = c("PROP_OBS","PROP_EST"),
               names_to = c("kind","junk"),names_sep = " ") %>%
  select(-junk)

complot_year<- long_temp2 %>%
  ggplot() +
  aes(x = BIN, y = value) +
  geom_line(aes(col=kind)) +
  theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  labs(x="length bin", y="proportion value", title="Length composition aggregated by year catch data") +
  facet_wrap(.~ESPECIE) +
  theme(legend.position = "bottom") +
  labs(col="") +
  guides(col = guide_legend(nrow = 1))

#ggsave(paste0(plotdir,"complot_year_catch.jpeg"), complot_year)#, width=3000, height=2500, res=250)
print(complot_year)





temporal2 = numeric()
especie = numeric(); especie = sort(unique(temp.surv$species)) # especies
for(e in 1:length(especie)){
  pos1 = numeric(); pos1 = which(temp.surv$species == especie[e])
  
  temporal1 = numeric()
  year = numeric(); year = sort(unique(temp.surv$year[pos1])) # ano
  for(y in 1:length(year)){
    pos2 = numeric(); pos2 = which(temp.surv$year[pos1] == year[y])
    
    N = numeric(); N = unique(temp.surv$InpN[pos1][pos2]) # sample size para el estrato especie ano
    n_ejemplares_obs = numeric(); n_ejemplares_obs = round(N * temp.surv$obs_value[pos1][pos2] , 0) # crea vector de proporcion por el sample size
    n_ejemplares_est = numeric(); n_ejemplares_est = round(N * temp.surv$pred_value[pos1][pos2], 0)
    
    tt1 = numeric()
    tt1 = data.frame(YEAR = year[y], BIN = temp.surv$lenbin[pos1][pos2], N = N, N_EJEM_OBS = n_ejemplares_obs, N_EJEM_EST = n_ejemplares_est)
    
    temporal1 = rbind(temporal1, tt1) # guarda la proporcion * sample size para la especie seleccionada y ano
  }
  
  bin1_obs = numeric(); bin1_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 1)])
  bin2_obs = numeric(); bin2_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 2)])
  bin3_obs = numeric(); bin3_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 3)])
  bin4_obs = numeric(); bin4_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 4)])
  bin5_obs = numeric(); bin5_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 5)])
  
  total_obs = sum(bin1_obs, bin2_obs, bin3_obs, bin4_obs, bin5_obs)
  prop_new_obs = numeric(); prop_new_obs = c(bin1_obs, bin2_obs, bin3_obs, bin4_obs, bin5_obs) / total_obs
  
  bin1_pred = numeric(); bin1_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 1)])
  bin2_pred = numeric(); bin2_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 2)])
  bin3_pred = numeric(); bin3_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 3)])
  bin4_pred = numeric(); bin4_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 4)])
  bin5_pred = numeric(); bin5_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 5)])
  
  total_pred = sum(bin1_pred, bin2_pred, bin3_pred, bin4_pred, bin5_pred)
  prop_new_est = numeric(); prop_new_est = c(bin1_pred, bin2_pred, bin3_pred, bin4_pred, bin5_pred) / total_pred
  
  tt2 = numeric()
  tt2 = data.frame(ESPECIE = especie[e], BIN = seq(1, 5, 1), PROP_OBS = prop_new_obs, PROP_EST = prop_new_est)
  
  temporal2 = rbind(temporal2, tt2) # almacena por cada especie
}

long_temp2 <- temporal2 %>%
  pivot_longer(cols = c("PROP_OBS","PROP_EST"),
               names_to = c("kind","junk"),names_sep = " ") %>%
  select(-junk)

complot_year<- long_temp2 %>%
  ggplot() +
  aes(x = BIN, y = value) +
  geom_line(aes(col=kind)) +
  theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  labs(x="length bin", y="proportion value", title="Length composition aggregated by year survey data") +
  facet_wrap(.~ESPECIE) +
  theme(legend.position = "bottom") +
  labs(col="") +
  guides(col = guide_legend(nrow = 1))

#ggsave(paste0(plotdir,"complot_year_survey.jpeg"), complot_year)#, width=3000, height=2500, res=250)
print(complot_year)



temp.surv <- temp.surv %>% dplyr::filter(number ==2)

### survey, species = 1 to 11
#tiff(paste0(plotdir,"bubbleplot_survey.jpeg"), width=3000, height=2500, res=250)
temp.surv<-ggplot(temp.surv, aes(x=year, y=lenbin, size = res_abs, color=factor(residual))) +
  geom_point(alpha=0.7) + theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  facet_wrap(.~species) + labs(x="year", y="length bin", title="Pearson residuals")
#dev.off()

ggsave("plots/pearson_survey2.jpeg",temp.surv, width = 10, height = 7, dpi = 300)#, width=3000, height=2500, res=250)


temp.catch <- temp.catch %>% dplyr::filter(number ==2)

temp.catch1<-ggplot(temp.catch, aes(x=year, y=lenbin, size = res_abs, color=factor(residual))) +
  geom_point(alpha=0.7) + theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  facet_wrap(.~species) + labs(x="year", y="length bin", title="Pearson residuals")
#dev.off()

ggsave("plots/pearson_catch2.jpeg",temp.catch1, width = 10, height = 7, dpi = 300)#, width=3000, height=2500, res=250)



plot_osa <- list()
especies<-unique(temp.surv$species)
for (sp in especies) {
  
  survdat <- temp.surv %>%
    filter(number == 2, species == sp)
  obs <- survdat %>% 
    select(year, lenbin, obs_value) %>% 
    mutate(obs_value = obs_value + 0.00001) %>% 
    pivot_wider(names_from = "year",
                values_from = "obs_value") %>% 
    column_to_rownames(var = "lenbin") %>% 
    as.matrix()
  pred <- survdat %>% 
    select(year, lenbin, pred_value) %>% 
    pivot_wider(names_from = "year",
                values_from = "pred_value") %>% 
    column_to_rownames(var = "lenbin") %>% 
    as.matrix()
  res<-resMulti(obs,pred)
  plot_osa[[sp]] <- res
}


plot_osa$Atlantic_mackerel

tmp <- plot_osa$Atlantic_mackerel

plot(tmp)

ggsave("plots/osa_S1_atlantic_cod.jpeg",tmp, width = 10, height = 7, dpi = 300)#, width=3000, height=2500, res=250)



plot_osa <- list()
especies<-unique(temp.catch$species)
for (sp in especies) {
  
  catchdat <- temp.catch %>%
    filter(number == 2, species == sp)
  obs <- catchdat %>% 
    select(year, lenbin, obs_value) %>% 
    mutate(obs_value = obs_value + 0.00001) %>% 
    pivot_wider(names_from = "year",
                values_from = "obs_value") %>% 
    column_to_rownames(var = "lenbin") %>% 
    as.matrix()
  pred <- catchdat %>% 
    select(year, lenbin, pred_value) %>% 
    pivot_wider(names_from = "year",
                values_from = "pred_value") %>% 
    column_to_rownames(var = "lenbin") %>% 
    as.matrix()
  res<-resMulti(obs,pred)
  plot_osa[[sp]] <- res
}


plot_osa$Atlantic_mackerel

tmp <- plot_osa$Atlantic_mackerel

plot(tmp)

ggsave("plots/osa_S1_atlantic_cod.jpeg",tmp, width = 10, height = 7, dpi = 300)#, width=3000, height=2500, res=250)


