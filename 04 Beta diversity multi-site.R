# Calculate multi-site beta diversity metrics for each taxonomic group

##Input community matrices are the same as used for the joint-species modeling

library(tidyverse)
library(betapart)



##birds ####
##using list Y_birds
rl<- as.data.frame(do.call(rbind, lapply(Y_birds, dim)))  #to get dimensions of each zone x decade data subset


for (u in 1:dim(rl)[1]){
  
  aux <- betapart.core(Y_birds[[u]][,-ncol(Y_birds[[u]])]) ##remove last col=model name
  
  assign(paste("betacore", unique(Y_birds[[u]][,ncol(Y_birds[[u]])]), sep="_"), aux)
  
  u<-u+1
}



##calculate beta metrics values using the min number of sites so that the estimates are comparable
multi_samp_bird_1_HB_SB <- beta.sample(betacore_bird_1_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_bird_1_MB <- beta.sample(betacore_bird_1_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_bird_1_NB <- beta.sample(betacore_bird_1_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_bird_2_HB_SB <- beta.sample(betacore_bird_2_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_bird_2_MB <- beta.sample(betacore_bird_2_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_bird_2_NB <- beta.sample(betacore_bird_2_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_bird_3_HB_SB <- beta.sample(betacore_bird_3_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_bird_3_MB <- beta.sample(betacore_bird_3_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_bird_3_NB <- beta.sample(betacore_bird_3_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_bird_4_HB_SB <- beta.sample(betacore_bird_4_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_bird_4_MB <- beta.sample(betacore_bird_4_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_bird_4_NB <- beta.sample(betacore_bird_4_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)


##combine and re-arrange the data to long format
multi_bird <- rbind(
  
  pivot_longer(multi_samp_bird_1_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "HB_SB"),
  
  pivot_longer(multi_samp_bird_1_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "MB"),
  
  pivot_longer(multi_samp_bird_1_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "NB"),
  
  pivot_longer(multi_samp_bird_2_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "HB_SB"),
  
  pivot_longer(multi_samp_bird_2_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "MB"),
  
  pivot_longer(multi_samp_bird_2_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "NB"),
  
  pivot_longer(multi_samp_bird_3_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "HB_SB"),
  
  pivot_longer(multi_samp_bird_3_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "MB"),
  
  pivot_longer(multi_samp_bird_3_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "NB"),
  
  pivot_longer(multi_samp_bird_4_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "HB_SB"),
  
  pivot_longer(multi_samp_bird_4_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "MB"),
  
  pivot_longer(multi_samp_bird_4_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "NB")) %>%
  
  mutate(Metric = ifelse(Beta_metric=="beta.SIM", "Turnover",
                             ifelse(Beta_metric=="beta.SNE", "Nestedness","Total Beta")))  ##better names for the plots


##to order the plots
multi_bird$Metric = factor(multi_bird$Metric, levels=c("Total Beta", "Turnover", "Nestedness"))



###
##moths ####
##using list Y_moth
rl<- as.data.frame(do.call(rbind, lapply(Y_moth, dim)))  #to get dimensions of each zone x decade data subset


for (u in 1:dim(rl)[1]){
  
  aux <- betapart.core(Y_moth[[u]][,-ncol(Y_moth[[u]])]) ##remove last col=model name
  
  assign(paste("betacore", unique(Y_moth[[u]][,ncol(Y_moth[[u]])]), sep="_"), aux)
  
  u<-u+1
}



##calculate beta metrics values using the min number of sites so that the estimates are comparable
multi_samp_moth_2_HB_SB <- beta.sample(betacore_Moth_2_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_moth_2_MB <- beta.sample(betacore_Moth_2_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_moth_2_NB <- beta.sample(betacore_Moth_2_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_moth_3_HB_SB <- beta.sample(betacore_Moth_3_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_moth_3_MB <- beta.sample(betacore_Moth_3_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_moth_3_NB <- beta.sample(betacore_Moth_3_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_moth_4_HB_SB <- beta.sample(betacore_Moth_4_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_moth_4_MB <- beta.sample(betacore_Moth_4_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_moth_4_NB <- beta.sample(betacore_Moth_4_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)


##combine and re-arrange the data to long format
multi_moth <- rbind(
  
  pivot_longer(multi_samp_moth_2_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "HB_SB"),
  
  pivot_longer(multi_samp_moth_2_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "MB"),
  
  pivot_longer(multi_samp_moth_2_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "NB"),
  
  pivot_longer(multi_samp_moth_3_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "HB_SB"),
  
  pivot_longer(multi_samp_moth_3_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "MB"),
  
  pivot_longer(multi_samp_moth_3_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "NB"),
  
  pivot_longer(multi_samp_moth_4_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "HB_SB"),
  
  pivot_longer(multi_samp_moth_4_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "MB"),
  
  pivot_longer(multi_samp_moth_4_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "NB")) %>%
  
  mutate(Metric = ifelse(Beta_metric=="beta.SIM", "Turnover",
                         ifelse(Beta_metric=="beta.SNE", "Nestedness","Total Beta")))  ##better names for the plots


##to order the plots
multi_moth$Metric = factor(multi_moth$Metric, levels=c("Total Beta", "Turnover", "Nestedness"))



###
##butterflies ####
##using list Y_butter
rl<- as.data.frame(do.call(rbind, lapply(Y_butter, dim)))  #to get dimensions of each zone x decade data subset


for (u in 1:dim(rl)[1]){
  
  aux <- betapart.core(Y_butter[[u]][,-ncol(Y_butter[[u]])]) ##remove last col=model name
  
  assign(paste("betacore", unique(Y_butter[[u]][,ncol(Y_butter[[u]])]), sep="_"), aux)
  
  u<-u+1
}



##calculate beta metrics values using the min number of sites so that the estimates are comparable
multi_samp_butter_3_HB_SB <- beta.sample(betacore_butterfly_3_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_butter_3_MB <- beta.sample(betacore_butterfly_3_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_butter_4_HB_SB <- beta.sample(betacore_butterfly_4_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_butter_4_MB <- beta.sample(betacore_butterfly_4_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)



##combine and re-arrange the data to long format
multi_butter <- rbind(
  
  pivot_longer(multi_samp_butter_3_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "HB_SB"),
  
  pivot_longer(multi_samp_butter_3_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "MB"),
  
  pivot_longer(multi_samp_butter_4_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "HB_SB"),
  
  pivot_longer(multi_samp_butter_4_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "MB")) %>%
  
  mutate(Metric = ifelse(Beta_metric=="beta.SIM", "Turnover",
                         ifelse(Beta_metric=="beta.SNE", "Nestedness","Total Beta")))  ##better names for the plots


##to order the plots
multi_butter$Metric = factor(multi_butter$Metric, levels=c("Total Beta", "Turnover", "Nestedness"))



###
##forest understory plants ####
##using list Y_forest
rl<- as.data.frame(do.call(rbind, lapply(Y_forest, dim)))  #to get dimensions of each zone x decade data subset


for (u in 1:dim(rl)[1]){
  
  aux <- betapart.core(Y_forest[[u]][,-ncol(Y_forest[[u]])]) ##remove last col=model name
  
  assign(paste("betacore", unique(Y_forest[[u]][,ncol(Y_forest[[u]])]), sep="_"), aux)
  
  u<-u+1
}



##calculate beta metrics values using the min number of sites so that the estimates are comparable
multi_samp_forest_1_HB_SB <- beta.sample(betacore_forest_1_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_forest_1_MB <- beta.sample(betacore_forest_1_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_forest_1_NB <- beta.sample(betacore_forest_1_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_forest_2_HB_SB <- beta.sample(betacore_forest_2_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_forest_2_MB <- beta.sample(betacore_forest_2_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_forest_2_NB <- beta.sample(betacore_forest_2_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_forest_3_HB_SB <- beta.sample(betacore_forest_3_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_forest_3_MB <- beta.sample(betacore_forest_3_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_forest_3_NB <- beta.sample(betacore_forest_3_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)


##combine and re-arrange the data to long format
multi_forest <- rbind(
  
  pivot_longer(multi_samp_forest_1_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "HB_SB"),
  
  pivot_longer(multi_samp_forest_1_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "MB"),
  
  pivot_longer(multi_samp_forest_1_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "NB"),
  
  pivot_longer(multi_samp_forest_2_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "HB_SB"),
  
  pivot_longer(multi_samp_forest_2_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "MB"),
  
  pivot_longer(multi_samp_forest_2_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "NB"),
  
  pivot_longer(multi_samp_forest_3_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "HB_SB"),
  
  pivot_longer(multi_samp_forest_3_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "MB"),
  
  pivot_longer(multi_samp_forest_3_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "NB")) %>%
  
  mutate(Metric = ifelse(Beta_metric=="beta.SIM", "Turnover",
                         ifelse(Beta_metric=="beta.SNE", "Nestedness","Total Beta")))  ##better names for the plots


##to order the plots
multi_forest$Metric = factor(multi_forest$Metric, levels=c("Total Beta", "Turnover", "Nestedness"))



###
##mammals ####
##using list Y_mamm
rl<- as.data.frame(do.call(rbind, lapply(Y_mamm, dim)))  #to get dimensions of each zone x decade data subset


for (u in 1:dim(rl)[1]){
  
  aux <- betapart.core(Y_mamm[[u]][,-ncol(Y_mamm[[u]])]) ##remove last col=model name
  
  assign(paste("betacore", unique(Y_mamm[[u]][,ncol(Y_mamm[[u]])]), sep="_"), aux)
  
  u<-u+1
}



##calculate beta metrics values using the min number of sites so that the estimates are comparable
multi_samp_mamm_2_HB_SB <- beta.sample(betacore_mammal_2_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_mamm_2_MB <- beta.sample(betacore_mammal_2_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_mamm_2_NB <- beta.sample(betacore_mammal_2_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_mamm_3_HB_SB <- beta.sample(betacore_mammal_3_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_mamm_3_MB <- beta.sample(betacore_mammal_3_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_mamm_3_NB <- beta.sample(betacore_mammal_3_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_mamm_4_HB_SB <- beta.sample(betacore_mammal_4_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_mamm_4_MB <- beta.sample(betacore_mammal_4_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_mamm_4_NB <- beta.sample(betacore_mammal_4_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)


##combine and re-arrange the data to long format
multi_mammal <- rbind(
  
  pivot_longer(multi_samp_mamm_2_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "HB_SB"),
  
  pivot_longer(multi_samp_mamm_2_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "MB"),
  
  pivot_longer(multi_samp_mamm_2_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "NB"),
  
  pivot_longer(multi_samp_mamm_3_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "HB_SB"),
  
  pivot_longer(multi_samp_mamm_3_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "MB"),
  
  pivot_longer(multi_samp_mamm_3_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "NB"),
  
  pivot_longer(multi_samp_mamm_4_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "HB_SB"),
  
  pivot_longer(multi_samp_mamm_4_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "MB"),
  
  pivot_longer(multi_samp_mamm_4_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "NB")) %>%
  
  mutate(Metric = ifelse(Beta_metric=="beta.SIM", "Turnover",
                         ifelse(Beta_metric=="beta.SNE", "Nestedness","Total Beta")))  ##better names for the plots


##to order the plots
multi_mammal$Metric = factor(multi_mammal$Metric, levels=c("Total Beta", "Turnover", "Nestedness"))



###
##phytoplankton ####
##using list Y_phyto
rl<- as.data.frame(do.call(rbind, lapply(Y_phyto, dim)))  #to get dimensions of each zone x decade data subset


for (u in 1:dim(rl)[1]){
  
  aux <- betapart.core(Y_phyto[[u]][,-ncol(Y_phyto[[u]])]) ##remove last col=model name
  
  assign(paste("betacore", unique(Y_phyto[[u]][,ncol(Y_phyto[[u]])]), sep="_"), aux)
  
  u<-u+1
}



##calculate beta metrics values using the min number of sites so that the estimates are comparable
multi_samp_phyto_1_HB_SB <- beta.sample(betacore_phytoplankton_1_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_phyto_1_MB <- beta.sample(betacore_phytoplankton_1_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_phyto_1_NB <- beta.sample(betacore_phytoplankton_1_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_phyto_2_HB_SB <- beta.sample(betacore_phytoplankton_2_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_phyto_2_MB <- beta.sample(betacore_phytoplankton_2_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_phyto_2_NB <- beta.sample(betacore_phytoplankton_2_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_phyto_3_HB_SB <- beta.sample(betacore_phytoplankton_3_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_phyto_3_MB <- beta.sample(betacore_phytoplankton_3_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_phyto_3_NB <- beta.sample(betacore_phytoplankton_3_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_phyto_4_HB_SB <- beta.sample(betacore_phytoplankton_4_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_phyto_4_MB <- beta.sample(betacore_phytoplankton_4_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_phyto_4_NB <- beta.sample(betacore_phytoplankton_4_NB, index.family="sorensen", sites=min(rl$V1), samples=1000)



##combine and re-arrange the data to long format
multi_phyto <- rbind(
  
  pivot_longer(multi_samp_phyto_1_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "HB_SB"),
  
  pivot_longer(multi_samp_phyto_1_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "MB"),
  
  pivot_longer(multi_samp_phyto_1_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "NB"),
  
  pivot_longer(multi_samp_phyto_2_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "HB_SB"),
  
  pivot_longer(multi_samp_phyto_2_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "MB"),
  
  pivot_longer(multi_samp_phyto_2_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "NB"),
  
  pivot_longer(multi_samp_phyto_3_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "HB_SB"),
  
  pivot_longer(multi_samp_phyto_3_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "MB"),
  
  pivot_longer(multi_samp_phyto_3_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "NB"),
  
  pivot_longer(multi_samp_phyto_4_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "HB_SB"),
  
  pivot_longer(multi_samp_phyto_4_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "MB"),
  
  pivot_longer(multi_samp_phyto_4_NB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "NB")) %>%
  
  mutate(Metric = ifelse(Beta_metric=="beta.SIM", "Turnover",
                         ifelse(Beta_metric=="beta.SNE", "Nestedness","Total Beta")))  ##better names for the plots


##to order the plots
multi_phyto$Metric = factor(multi_phyto$Metric, levels=c("Total Beta", "Turnover", "Nestedness"))




###
##small mammals ####
##using list Y_vole
rl<- as.data.frame(do.call(rbind, lapply(Y_vole, dim)))  #to get dimensions of each zone x decade data subset


for (u in 1:dim(rl)[1]){
  
  aux <- betapart.core(Y_vole[[u]][,-ncol(Y_vole[[u]])]) ##remove last col=model name
  
  assign(paste("betacore", unique(Y_vole[[u]][,ncol(Y_vole[[u]])]), sep="_"), aux)
  
  u<-u+1
}



##calculate beta metrics values using the min number of sites so that the estimates are comparable
multi_samp_vole_1_HB_SB <- beta.sample(betacore_Vole_1_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_vole_1_MB <- beta.sample(betacore_Vole_1_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_vole_2_HB_SB <- beta.sample(betacore_Vole_2_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_vole_2_MB <- beta.sample(betacore_Vole_2_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_vole_3_HB_SB <- beta.sample(betacore_Vole_3_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_vole_3_MB <- beta.sample(betacore_Vole_3_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)

multi_samp_vole_4_HB_SB <- beta.sample(betacore_Vole_4_HB_SB, index.family="sorensen", sites=min(rl$V1), samples=1000)
multi_samp_vole_4_MB <- beta.sample(betacore_Vole_4_MB, index.family="sorensen", sites=min(rl$V1), samples=1000)


##combine and re-arrange the data to long format
multi_vole <- rbind(
  
  pivot_longer(multi_samp_vole_1_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "HB_SB"),
  
  pivot_longer(multi_samp_vole_1_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=1, Region= "MB"),
  
  pivot_longer(multi_samp_vole_2_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "HB_SB"),
  
  pivot_longer(multi_samp_vole_2_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=2, Region= "MB"),
  
  pivot_longer(multi_samp_vole_3_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "HB_SB"),
  
  pivot_longer(multi_samp_vole_3_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=3, Region= "MB"),
  
  pivot_longer(multi_samp_vole_4_HB_SB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "HB_SB"),
  
  pivot_longer(multi_samp_vole_4_MB$sampled.values, cols = c(1:3), names_to = "Beta_metric", values_to = "values") %>%
    mutate(Decade=4, Region= "MB")) %>%
  
  mutate(Metric = ifelse(Beta_metric=="beta.SIM", "Turnover",
                         ifelse(Beta_metric=="beta.SNE", "Nestedness","Total Beta")))  ##better names for the plots


##to order the plots
multi_vole$Metric = factor(multi_vole$Metric, levels=c("Total Beta", "Turnover", "Nestedness"))




##end


