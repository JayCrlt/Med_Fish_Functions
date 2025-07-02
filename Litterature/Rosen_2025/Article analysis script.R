# Clear R, if needed
rm(list=ls())

# load required packages
Packages <- c("ggplot2","patchwork","lmerTest","lmtest","dplyr","Rmisc","ggpubr", "performance","MuMIn","multcomp","marginaleffects","tidyverse","jtools","broom.mixed","huxtable","ggeffects","MCMCglmm","rptR", "scales","grid","ggh4x")
lapply(Packages, library, character.only = TRUE)

# load data
#Data_raw <- read.csv("C:/Users/tnor/Dropbox/Work/DTU/Metabolic scaling project/Scaling data combined/Data_all.csv", sep=";", encoding ="latin1")
Data_raw <- read.csv("C:/Users/alero/My Drive/phd/Artikels/Big article/Data_analysis/Data_all.csv", sep=";", encoding ="latin1")

names(Data_raw) <- c("Species","Date","FishID","Measurement","Treatment","Mass","SMR","MMR","AS","Comment", "Chambersize","Temp", "Sex")
Data_raw$Measurement <- as.numeric(Data_raw$Measurement) #Makiong sure messurement nr is numeric
Data_raw$Date <- dmy(Data_raw$Date) # Making sure date is in a date formate

# grouping Species and FishID to create unique identifier, as some individuals have same ID across species
Data_raw$Group <- paste(Data_raw$Species, Data_raw$FishID, sep="_")

#Ordering
Data_raw <- Data_raw[order(Data_raw$Species,Data_raw$FishID,Data_raw$Measurement),] #making the dataframe more ordered to look at

#determine if RR is always used in analysis, RI is alwasy also compared but not used
use_RR=TRUE
use_RRM=TRUE
use_treatments=TRUE
# remove fish with less than four repeated Trials
Data_raw <- Data_raw %>% group_by(FishID) %>% filter(n() > 3)




# Adding new collums to calculate scaling of growth rate
Data_raw$Day_since_start <- 0
Data_raw$SGR <- NA
Data_raw$AGR <- NA
Data_raw$MassGrowth <- NA
Data_raw$MassMean <- NA
Data_raw$MassCent <- NA
Data_raw$MassMeanGrowth <- NA
Data_raw$MassCentGrowth <- NA
Data_raw$logMassCent <- NA
Data_raw$logMassCentGrowth <- NA
Data_raw$SMRGrowth <- NA
Data_raw$GE <- NA


# Finding days since start of first messuremnt for all subsequnt messurements for each individual fish, SGR and AGR is the standart growthrate and abselut growth rate bewteen the last messurement and the given messurement,
# MassGrowth is the awege mass between the last and the current messurement, MassMean is the avgee Mass for a given individual, MassCent is the centralized mass for each individual
for (x in 1:length(Data_raw$Date)) {
  Data_ID <- Data_raw[Data_raw$Species==Data_raw$Species[x] & Data_raw$FishID==Data_raw$FishID[x],]
  Data_raw$Day_since_start[x] <- difftime(Data_raw$Date[x] ,min(Data_ID$Date) , units = c("days"))
  last_mes <- Data_ID[Data_ID$Measurement==(Data_raw$Measurement[x]-1),]
  if(length(last_mes$Mass)>0){
    Data_raw$SGR[x] <- (log(Data_raw$Mass[x])-log(last_mes$Mass))/(Data_raw$Day_since_start[x]-last_mes$Day_since_start)*100
    Data_raw$AGR[x] <- (Data_raw$Mass[x]-last_mes$Mass)/(Data_raw$Day_since_start[x]-last_mes$Day_since_start)
    Data_raw$MassGrowth[x] <- (Data_raw$Mass[x]+last_mes$Mass)/2
    }
  Data_raw$MassMean[x] <- mean(Data_ID$Mass)
  Data_raw$MassCent[x] <- Data_raw$Mass[x] - Data_raw$MassMean[x]
  Data_raw$logMassCent[x] <- log10(Data_raw$Mass[x]) - log10(Data_raw$MassMean[x])
}


# Using MassGrowth to calculate centrylized values since that is used to calculate scaling of growth, while using SMR is used so SMR and growth values can be direktly compared
for (x in 1:length(Data_raw$Date)) {
  Data_ID <- Data_raw[Data_raw$Species==Data_raw$Species[x] & Data_raw$FishID==Data_raw$FishID[x],]
  last_mes <- Data_ID[Data_ID$Measurement==(Data_raw$Measurement[x]-1),]
  if(length(last_mes$Mass)>0){
    Data_raw$MassMeanGrowth[x] <- mean(Data_ID[!is.na(Data_ID$MassGrowth),]$MassGrowth)
    Data_raw$MassCentGrowth[x] <- Data_raw$MassGrowth[x] - Data_raw$MassMeanGrowth[x]
    Data_raw$logMassCentGrowth[x] <- log10(Data_raw$MassGrowth[x]) - log10(Data_raw$MassMeanGrowth[x])
    Data_raw$SMRGrowth[x] <- (Data_raw$SMR[x]+last_mes$SMR)/2
    Data_raw$GE[x] <- Data_raw$AGR[x]/Data_raw$SMRGrowth[x]*24
    }
  }


# replace negative absolute growth rates for (zebra)fish with NAs
Data_raw <- Data_raw %>% mutate(AGR = ifelse(AGR > 0, AGR, NA))
Data_raw <- Data_raw[(!Data_raw$Sex=="n"&!Data_raw$Sex=="d")|is.na(Data_raw$Sex),]
 
# log transformation
Data_raw$logMass <- log10(Data_raw$Mass)
Data_raw$logSMR <- log10(Data_raw$SMR)
Data_raw$logMMR <- log10(Data_raw$MMR)
Data_raw$logAS <- log10(Data_raw$AS)
Data_raw$logMassGrowth <- log10(Data_raw$MassGrowth)
Data_raw$logSGR <- log10(Data_raw$SGR)
Data_raw$logAGR <- log10(Data_raw$AGR)
Data_raw$logMassMean <- log10(Data_raw$MassMean)
Data_raw$logMassMeanGrowth <- log10(Data_raw$MassMeanGrowth)
Data_raw$logSMRGrowth <- log10(Data_raw$SMRGrowth)
Data_raw$logGE <- log10(Data_raw$GE)





##### MCMCglmm models #####
prior_1t_G1 <- list(R = list(V = 1, nu = 0.002),
                    G = list(G1 = list(V = 1,
                                       nu = 1,
                                       alpha.mu = 0,
                                       alpha.V = 25^2)))

prior_1t_RR_G1 <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1 = list(V = diag(2),
                                          nu = 2,
                                          alpha.mu = c(0,0),
                                          alpha.V = diag(2)*25^2)))








#Handy table to keep track of which model is best:
Modeluse <- data.frame(matrix(ncol = 13,nrow = 0))
names(Modeluse) <- c("Species","SMR","SMR_dif","MMR","MMR_dif","AS","AS_dif","AGR","AGR_dif","SGR","SGR_dif","GE","GE_dif")
ModeluseM <- data.frame(matrix(ncol = 13,nrow = 0))
names(ModeluseM) <- c("Species","SMR","SMR_dif","MMR","MMR_dif","AS","AS_dif","AGR","AGR_dif","SGR","SGR_dif","GE","GE_dif")




### SMR ###
#SMR dataframe
Data_raw_SMR <- Data_raw
Data_raw_SMR$SMR_staticM <- NA
Data_raw_SMR$SMR_staticT <- NA

#Subsetting data
Data_raw_SMR_Zebrafish <- subset(Data_raw_SMR, Species=="Zebrafish")
Data_raw_SMR_Rainbow_trout <- subset(Data_raw_SMR, Species=="Rainbow_trout")
Data_raw_SMR_Brown_trout <- subset(Data_raw_SMR, Species=="Brown_trout")
Data_raw_SMR_Cunner <- subset(Data_raw_SMR, Species=="Cunner")
Data_raw_SMR_Guppy <- subset(Data_raw_SMR, Species=="Guppy", Sex != "n" & Sex != "d")
Data_raw_SMR_Chromis <- subset(Data_raw_SMR, Species=="Chromis")
Data_raw_SMR_Clownfish <- subset(Data_raw_SMR, Species=="Clownfish")


## Zebrafish ##
#With treatments
mcmc_SMR_RR_Zebrafish  <- MCMCglmm(logSMR ~ logMassCent * Treatment + logMassMean * Treatment + Temp,
                                                  random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                                  prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Zebrafish))

# mcmc_SMR_RI_Zebrafish  <- MCMCglmm(logSMR ~ logMassCent * Treatment + logMassMean * Treatment + Temp,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Zebrafish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SMR_Zebrafish_use <- mcmc_SMR_RR_Zebrafish,ifelse(mcmc_SMR_RR_Zebrafish$DIC<mcmc_SMR_RI_Zebrafish$DIC,mcmc_SMR_Zebrafish_use <- mcmc_SMR_RR_Zebrafish,mcmc_SMR_Zebrafish_use <- mcmc_SMR_RI_Zebrafish))
# Modeluse[1,1] <- "Zebrafish"
# Modeluse[1,2] <-ifelse(mcmc_SMR_RR_Zebrafish$DIC<mcmc_SMR_RI_Zebrafish$DIC, "RR","RI")
# Modeluse[1,3] <- abs(mcmc_SMR_RR_Zebrafish$DIC-mcmc_SMR_RI_Zebrafish$DIC)

#Without treatments
mcmc_SMR_RR_ZebrafishM  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean  + Temp,
                                   random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Zebrafish))

# mcmc_SMR_RI_ZebrafishM  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean  + Temp,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Zebrafish))
# 
ifelse(use_RRM==TRUE,mcmc_SMR_Zebrafish_useM <- mcmc_SMR_RR_ZebrafishM,ifelse(mcmc_SMR_RR_ZebrafishM$DIC<mcmc_SMR_RI_ZebrafishM$DIC,mcmc_SMR_Zebrafish_useM <- mcmc_SMR_RR_ZebrafishM,mcmc_SMR_Zebrafish_useM <- mcmc_SMR_RI_ZebrafishM))
# ModeluseM[1,1] <- "Zebrafish"
# ModeluseM[1,2] <-ifelse(mcmc_SMR_RR_ZebrafishM$DIC<mcmc_SMR_RI_ZebrafishM$DIC, "RR","RI")
# ModeluseM[1,3] <- abs(mcmc_SMR_RR_ZebrafishM$DIC-mcmc_SMR_RI_ZebrafishM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_SMR_pred_full_Zebrafish <- cbind(Data_raw_SMR_Zebrafish, predict(mcmc_SMR_Zebrafish_use, marginal=NULL, interval = "prediction")),
                                                                df_SMR_pred_full_Zebrafish <- cbind(Data_raw_SMR_Zebrafish, predict(mcmc_SMR_Zebrafish_useM, marginal=NULL, interval = "prediction"))) 


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SMR_pred_full_Zebrafish$SMR_staticM <- df_SMR_pred_full_Zebrafish$logMass*summary(mcmc_SMR_Zebrafish_useM)$solutions[3,1]+summary(mcmc_SMR_Zebrafish_useM)$solutions[1,1]+summary(mcmc_SMR_Zebrafish_useM)$solutions[4,1]*df_SMR_pred_full_Zebrafish$Temp

# finding the statick scaling value for a given individual at a given point for the given treatment
df_SMR_pred_full_Zebrafish[df_SMR_pred_full_Zebrafish$Treatment=="H",]$SMR_staticT <- df_SMR_pred_full_Zebrafish[df_SMR_pred_full_Zebrafish$Treatment=="H",]$logMass*summary(mcmc_SMR_Zebrafish_use)$solutions[5,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[6,1]*df_SMR_pred_full_Zebrafish[df_SMR_pred_full_Zebrafish$Treatment=="H",]$Temp+summary(mcmc_SMR_Zebrafish_use)$solutions[1,1]
df_SMR_pred_full_Zebrafish[df_SMR_pred_full_Zebrafish$Treatment=="L",]$SMR_staticT <- df_SMR_pred_full_Zebrafish[df_SMR_pred_full_Zebrafish$Treatment=="L",]$logMass*(summary(mcmc_SMR_Zebrafish_use)$solutions[5,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[9,1])+summary(mcmc_SMR_Zebrafish_use)$solutions[6,1]*df_SMR_pred_full_Zebrafish[df_SMR_pred_full_Zebrafish$Treatment=="L",]$Temp+summary(mcmc_SMR_Zebrafish_use)$solutions[1,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[3,1]
df_SMR_pred_full_Zebrafish[df_SMR_pred_full_Zebrafish$Treatment=="M",]$SMR_staticT <- df_SMR_pred_full_Zebrafish[df_SMR_pred_full_Zebrafish$Treatment=="M",]$logMass*(summary(mcmc_SMR_Zebrafish_use)$solutions[5,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[10,1])+summary(mcmc_SMR_Zebrafish_use)$solutions[6,1]*df_SMR_pred_full_Zebrafish[df_SMR_pred_full_Zebrafish$Treatment=="M",]$Temp+summary(mcmc_SMR_Zebrafish_use)$solutions[1,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[4,1]



## Rainbow trout ##
#With treatments
mcmc_SMR_RR_Rainbow_trout  <- MCMCglmm(logSMR ~ logMassCent * Treatment + logMassMean * Treatment + Temp,
                                   random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Rainbow_trout))

# mcmc_SMR_RI_Rainbow_trout  <- MCMCglmm(logSMR ~ logMassCent * Treatment + logMassMean * Treatment + Temp,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Rainbow_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SMR_Rainbow_trout_use <- mcmc_SMR_RR_Rainbow_trout,ifelse(mcmc_SMR_RR_Rainbow_trout$DIC<mcmc_SMR_RI_Rainbow_trout$DIC,mcmc_SMR_Rainbow_trout_use <- mcmc_SMR_RR_Rainbow_trout,mcmc_SMR_Rainbow_trout_use <- mcmc_SMR_RI_Rainbow_trout))
# Modeluse[2,1] <- "Rainbow_trout"
# Modeluse[2,2] <-ifelse(mcmc_SMR_RR_Rainbow_trout$DIC<mcmc_SMR_RI_Rainbow_trout$DIC, "RR","RI")
# Modeluse[2,3] <- abs(mcmc_SMR_RR_Rainbow_trout$DIC-mcmc_SMR_RI_Rainbow_trout$DIC)

#Without treatments
mcmc_SMR_RR_Rainbow_troutM  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean  + Temp,
                                    random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                    prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Rainbow_trout))

# mcmc_SMR_RI_Rainbow_troutM  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean  + Temp,
#                                     random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                     prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Rainbow_trout))

ifelse(use_RRM==TRUE,mcmc_SMR_Rainbow_trout_useM <- mcmc_SMR_RR_Rainbow_troutM,ifelse(mcmc_SMR_RR_Rainbow_troutM$DIC<mcmc_SMR_RI_Rainbow_troutM$DIC,mcmc_SMR_Rainbow_trout_useM <- mcmc_SMR_RR_Rainbow_troutM,mcmc_SMR_Rainbow_trout_useM <- mcmc_SMR_RI_Rainbow_troutM))
# ModeluseM[2,1] <- "Rainbow_trout"
# ModeluseM[2,2] <-ifelse(mcmc_SMR_RR_Rainbow_troutM$DIC<mcmc_SMR_RI_Rainbow_troutM$DIC, "RR","RI")
# ModeluseM[2,3] <- abs(mcmc_SMR_RR_Rainbow_troutM$DIC-mcmc_SMR_RI_Rainbow_troutM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_SMR_pred_full_Rainbow_trout <- cbind(Data_raw_SMR_Rainbow_trout, predict(mcmc_SMR_Rainbow_trout_use, marginal=NULL, interval = "prediction")),
       df_SMR_pred_full_Rainbow_trout <- cbind(Data_raw_SMR_Rainbow_trout, predict(mcmc_SMR_Rainbow_trout_useM, marginal=NULL, interval = "prediction"))) 


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SMR_pred_full_Rainbow_trout$SMR_staticM <- df_SMR_pred_full_Rainbow_trout$logMass*summary(mcmc_SMR_Rainbow_trout_useM)$solutions[3,1]+summary(mcmc_SMR_Rainbow_trout_useM)$solutions[1,1]+summary(mcmc_SMR_Rainbow_trout_useM)$solutions[4,1]*df_SMR_pred_full_Rainbow_trout$Temp

# finding the statick scaling value for a given individual at a given point for the given treatment
df_SMR_pred_full_Rainbow_trout[df_SMR_pred_full_Rainbow_trout$Treatment=="10Cegg",]$SMR_staticT <- df_SMR_pred_full_Rainbow_trout[df_SMR_pred_full_Rainbow_trout$Treatment=="10Cegg",]$logMass*summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[6,1]*df_SMR_pred_full_Rainbow_trout[df_SMR_pred_full_Rainbow_trout$Treatment=="10Cegg",]$Temp+summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,1]
df_SMR_pred_full_Rainbow_trout[df_SMR_pred_full_Rainbow_trout$Treatment=="14Cegg",]$SMR_staticT <- df_SMR_pred_full_Rainbow_trout[df_SMR_pred_full_Rainbow_trout$Treatment=="14Cegg",]$logMass*(summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[9,1])+summary(mcmc_SMR_Rainbow_trout_use)$solutions[6,1]*df_SMR_pred_full_Rainbow_trout[df_SMR_pred_full_Rainbow_trout$Treatment=="14Cegg",]$Temp+summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[3,1]
df_SMR_pred_full_Rainbow_trout[df_SMR_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$SMR_staticT <- df_SMR_pred_full_Rainbow_trout[df_SMR_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$logMass*(summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[10,1])+summary(mcmc_SMR_Rainbow_trout_use)$solutions[6,1]*df_SMR_pred_full_Rainbow_trout[df_SMR_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$Temp+summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[4,1]


## Brown trout ##
mcmc_SMR_RR_Brown_trout  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean ,
                                       random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                       prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                       verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Brown_trout))

mcmc_SMR_RI_Brown_trout  <- MCMCglmm(logSMR ~ logMassCent + logMassMean ,
                                       random = ~ FishID, rcov = ~ units, family = "gaussian",
                                       prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
                                       verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Brown_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SMR_Brown_trout_use <- mcmc_SMR_RR_Brown_trout,ifelse(mcmc_SMR_RR_Brown_trout$DIC<mcmc_SMR_RI_Brown_trout$DIC,mcmc_SMR_Brown_trout_use <- mcmc_SMR_RR_Brown_trout,mcmc_SMR_Brown_trout_use <- mcmc_SMR_RI_Brown_trout))
# Modeluse[3,1] <- "Brown_trout"
# Modeluse[3,2] <-ifelse(mcmc_SMR_RR_Brown_trout$DIC<mcmc_SMR_RI_Brown_trout$DIC, "RR","RI")
# Modeluse[3,3] <- abs(mcmc_SMR_RR_Brown_trout$DIC-mcmc_SMR_RI_Brown_trout$DIC)

df_SMR_pred_full_Brown_trout <- cbind(Data_raw_SMR_Brown_trout, predict(mcmc_SMR_Brown_trout_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_SMR_pred_full_Brown_trout$SMR_staticT <- summary(mcmc_SMR_Brown_trout_use)$solutions[1,1] + summary(mcmc_SMR_Brown_trout_use)$solutions[3,1]*df_SMR_pred_full_Brown_trout$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SMR_pred_full_Brown_trout$SMR_staticM <- summary(mcmc_SMR_Brown_trout_use)$solutions[1,1] + summary(mcmc_SMR_Brown_trout_use)$solutions[3,1]*df_SMR_pred_full_Brown_trout$logMass


## Cunner ##
mcmc_SMR_RR_Cunner  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean ,
                                     random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                     prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Cunner))

mcmc_SMR_RI_Cunner  <- MCMCglmm(logSMR ~ logMassCent + logMassMean ,
                                     random = ~ FishID, rcov = ~ units, family = "gaussian",
                                     prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Cunner))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SMR_Cunner_use <- mcmc_SMR_RR_Cunner,ifelse(mcmc_SMR_RR_Cunner$DIC<mcmc_SMR_RI_Cunner$DIC,mcmc_SMR_Cunner_use <- mcmc_SMR_RR_Cunner,mcmc_SMR_Cunner_use <- mcmc_SMR_RI_Cunner))
# Modeluse[4,1] <- "Cunner"
# Modeluse[4,2] <-ifelse(mcmc_SMR_RR_Cunner$DIC<mcmc_SMR_RI_Cunner$DIC, "RR","RI")
# Modeluse[4,3] <- abs(mcmc_SMR_RR_Cunner$DIC-mcmc_SMR_RI_Cunner$DIC)

df_SMR_pred_full_Cunner <- cbind(Data_raw_SMR_Cunner, predict(mcmc_SMR_Cunner_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_SMR_pred_full_Cunner$SMR_staticT <- summary(mcmc_SMR_Cunner_use)$solutions[1,1] + summary(mcmc_SMR_Cunner_use)$solutions[3,1]*df_SMR_pred_full_Cunner$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SMR_pred_full_Cunner$SMR_staticM <- summary(mcmc_SMR_Cunner_use)$solutions[1,1] + summary(mcmc_SMR_Cunner_use)$solutions[3,1]*df_SMR_pred_full_Cunner$logMass


## Guppy ##
#with sex
mcmc_SMR_RR_Guppy  <- MCMCglmm(logSMR ~ logMassCent * Sex + logMassMean * Sex,
                                     random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                     prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Guppy))

# mcmc_SMR_RI_Guppy  <- MCMCglmm(logSMR ~ logMassCent * Sex + logMassMean * Sex,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Guppy))

#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SMR_Guppy_use <- mcmc_SMR_RR_Guppy,ifelse(mcmc_SMR_RR_Guppy$DIC<mcmc_SMR_RI_Guppy$DIC,mcmc_SMR_Guppy_use <- mcmc_SMR_RR_Guppy,mcmc_SMR_Guppy_use <- mcmc_SMR_RI_Guppy))
# Modeluse[5,1] <- "Guppy"
# Modeluse[5,2] <-ifelse(mcmc_SMR_RR_Guppy$DIC<mcmc_SMR_RI_Guppy$DIC, "RR","RI")
# Modeluse[5,3] <- abs(mcmc_SMR_RR_Guppy$DIC-mcmc_SMR_RI_Guppy$DIC)


#without sex
mcmc_SMR_RR_GuppyM  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean ,
                               random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                               prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Guppy))

# mcmc_SMR_RI_GuppyM  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean ,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Guppy))

ifelse(use_RRM==TRUE,mcmc_SMR_Guppy_useM <- mcmc_SMR_RR_GuppyM,ifelse(mcmc_SMR_RR_GuppyM$DIC<mcmc_SMR_RI_GuppyM$DIC,mcmc_SMR_Guppy_useM <- mcmc_SMR_RR_GuppyM,mcmc_SMR_Guppy_useM <- mcmc_SMR_RI_GuppyM))
# ModeluseM[5,1] <- "Guppy"
# ModeluseM[5,2] <-ifelse(mcmc_SMR_RR_GuppyM$DIC<mcmc_SMR_RI_GuppyM$DIC, "RR","RI")
# ModeluseM[5,3] <- abs(mcmc_SMR_RR_GuppyM$DIC-mcmc_SMR_RI_GuppyM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_SMR_pred_full_Guppy <- cbind(Data_raw_SMR_Guppy, predict(mcmc_SMR_Guppy_use, marginal=NULL, interval = "prediction")),
       df_SMR_pred_full_Guppy <- cbind(Data_raw_SMR_Guppy, predict(mcmc_SMR_Guppy_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_SMR_pred_full_Guppy[df_SMR_pred_full_Guppy$Sex=="f",]$SMR_staticT <- summary(mcmc_SMR_Guppy_use)$solutions[1,1] + summary(mcmc_SMR_Guppy_use)$solutions[4,1]*df_SMR_pred_full_Guppy[df_SMR_pred_full_Guppy$Sex=="f",]$logMass
df_SMR_pred_full_Guppy[df_SMR_pred_full_Guppy$Sex=="m",]$SMR_staticT <- summary(mcmc_SMR_Guppy_use)$solutions[1,1]+summary(mcmc_SMR_Guppy_use)$solutions[3,1] + (summary(mcmc_SMR_Guppy_use)$solutions[4,1]+summary(mcmc_SMR_Guppy_use)$solutions[6,1])*df_SMR_pred_full_Guppy[df_SMR_pred_full_Guppy$Sex=="m",]$logMass

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SMR_pred_full_Guppy$SMR_staticM <- summary(mcmc_SMR_Guppy_useM)$solutions[1,1] + summary(mcmc_SMR_Guppy_useM)$solutions[3,1]*df_SMR_pred_full_Guppy$logMass


## Chromis ##
#With treatments
mcmc_SMR_RR_Chromis  <- MCMCglmm(logSMR ~ logMassCent * Treatment + logMassMean * Treatment ,
                                   random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Chromis))

# mcmc_SMR_RI_Chromis  <- MCMCglmm(logSMR ~ logMassCent * Treatment + logMassMean * Treatment,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Chromis))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SMR_Chromis_use <- mcmc_SMR_RR_Chromis,ifelse(mcmc_SMR_RR_Chromis$DIC<mcmc_SMR_RI_Chromis$DIC,mcmc_SMR_Chromis_use <- mcmc_SMR_RR_Chromis,mcmc_SMR_Chromis_use <- mcmc_SMR_RI_Chromis))
# Modeluse[6,1] <- "Chromis"
# Modeluse[6,2] <-ifelse(mcmc_SMR_RR_Chromis$DIC<mcmc_SMR_RI_Chromis$DIC, "RR","RI")
# Modeluse[6,3] <- abs(mcmc_SMR_RR_Chromis$DIC-mcmc_SMR_RI_Chromis$DIC)

#Without treatments
mcmc_SMR_RR_ChromisM  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean ,
                                    random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                    prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Chromis))

# mcmc_SMR_RI_ChromisM  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean,
#                                     random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                     prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Chromis))

ifelse(use_RRM==TRUE,mcmc_SMR_Chromis_useM <- mcmc_SMR_RR_ChromisM,ifelse(mcmc_SMR_RR_ChromisM$DIC<mcmc_SMR_RI_ChromisM$DIC,mcmc_SMR_Chromis_useM <- mcmc_SMR_RR_ChromisM,mcmc_SMR_Chromis_useM <- mcmc_SMR_RI_ChromisM))
# ModeluseM[6,1] <- "Chromis"
# ModeluseM[6,2] <-ifelse(mcmc_SMR_RR_ChromisM$DIC<mcmc_SMR_RI_ChromisM$DIC, "RR","RI")
# ModeluseM[6,3] <- abs(mcmc_SMR_RR_ChromisM$DIC-mcmc_SMR_RI_ChromisM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_SMR_pred_full_Chromis <- cbind(Data_raw_SMR_Chromis, predict(mcmc_SMR_Chromis_use, marginal=NULL, interval = "prediction")),
       df_SMR_pred_full_Chromis <- cbind(Data_raw_SMR_Chromis, predict(mcmc_SMR_Chromis_useM, marginal=NULL, interval = "prediction"))) 



# finding the statick scaling value for a given individual at a given point for the given treatment
df_SMR_pred_full_Chromis[df_SMR_pred_full_Chromis$Treatment=="high",]$SMR_staticT <- summary(mcmc_SMR_Chromis_use)$solutions[1,1] + summary(mcmc_SMR_Chromis_use)$solutions[4,1]*df_SMR_pred_full_Chromis[df_SMR_pred_full_Chromis$Treatment=="high",]$logMass
df_SMR_pred_full_Chromis[df_SMR_pred_full_Chromis$Treatment=="low",]$SMR_staticT <- summary(mcmc_SMR_Chromis_use)$solutions[1,1]+summary(mcmc_SMR_Chromis_use)$solutions[3,1] + (summary(mcmc_SMR_Chromis_use)$solutions[4,1]+summary(mcmc_SMR_Chromis_use)$solutions[6,1])*df_SMR_pred_full_Chromis[df_SMR_pred_full_Chromis$Treatment=="low",]$logMass
#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SMR_pred_full_Chromis$SMR_staticM <- summary(mcmc_SMR_Chromis_useM)$solutions[1,1] + summary(mcmc_SMR_Chromis_useM)$solutions[3,1]*df_SMR_pred_full_Chromis$logMass



## Clownfish ##
mcmc_SMR_RR_Clownfish  <- MCMCglmm(logSMR ~ logMassCent  + logMassMean ,
                                random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Clownfish))

# mcmc_SMR_RI_Clownfish  <- MCMCglmm(logSMR ~ logMassCent + logMassMean ,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SMR_Clownfish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SMR_Clownfish_use <- mcmc_SMR_RR_Clownfish,ifelse(mcmc_SMR_RR_Clownfish$DIC<mcmc_SMR_RI_Clownfish$DIC,mcmc_SMR_Clownfish_use <- mcmc_SMR_RR_Clownfish,mcmc_SMR_Clownfish_use <- mcmc_SMR_RI_Clownfish))
# Modeluse[7,1] <- "Clownfish"
# Modeluse[7,2] <-ifelse(mcmc_SMR_RR_Clownfish$DIC<mcmc_SMR_RI_Clownfish$DIC, "RR","RI")
# Modeluse[7,3] <- abs(mcmc_SMR_RR_Clownfish$DIC-mcmc_SMR_RI_Clownfish$DIC)

df_SMR_pred_full_Clownfish <- cbind(Data_raw_SMR_Clownfish, predict(mcmc_SMR_Clownfish_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_SMR_pred_full_Clownfish$SMR_staticT <- summary(mcmc_SMR_Clownfish_use)$solutions[1,1] + summary(mcmc_SMR_Clownfish_use)$solutions[3,1]*df_SMR_pred_full_Clownfish$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SMR_pred_full_Clownfish$SMR_staticM <- summary(mcmc_SMR_Clownfish_use)$solutions[1,1] + summary(mcmc_SMR_Clownfish_use)$solutions[3,1]*df_SMR_pred_full_Clownfish$logMass


# Recombining #
Data_SMR <- bind_rows(df_SMR_pred_full_Zebrafish, df_SMR_pred_full_Rainbow_trout, df_SMR_pred_full_Brown_trout, df_SMR_pred_full_Cunner, df_SMR_pred_full_Guppy, df_SMR_pred_full_Chromis, df_SMR_pred_full_Clownfish)

names(Data_SMR)[names(Data_SMR) == "fit" ] <- "SMR_fit"
names(Data_SMR)[names(Data_SMR) == "lwr" ] <- "SMR_fit_low"
names(Data_SMR)[names(Data_SMR) == "upr" ] <- "SMR_fit_high"


### MMR ###
#MMR dataframe
Data_raw_MMR <- Data_raw
Data_raw_MMR$MMR_staticM <- NA
Data_raw_MMR$MMR_staticT <- NA

#Subsetting data
Data_raw_MMR_Rainbow_trout <- subset(Data_raw_MMR, Species=="Rainbow_trout")
Data_raw_MMR_Brown_trout <- subset(Data_raw_MMR, Species=="Brown_trout")
Data_raw_MMR_Cunner <- subset(Data_raw_MMR, Species=="Cunner")
Data_raw_MMR_Guppy <- subset(Data_raw_MMR, Species=="Guppy", Sex != "n" & Sex != "d")
Data_raw_MMR_Chromis <- subset(Data_raw_MMR, Species=="Chromis")
Data_raw_MMR_Clownfish <- subset(Data_raw_MMR, Species=="Clownfish")


## Rainbow trout ##
#With treatments
mcmc_MMR_RR_Rainbow_trout  <- MCMCglmm(logMMR ~ logMassCent * Treatment + logMassMean * Treatment + Temp,
                                       random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                       prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                       verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Rainbow_trout))

# mcmc_MMR_RI_Rainbow_trout  <- MCMCglmm(logMMR ~ logMassCent * Treatment + logMassMean * Treatment + Temp,
#                                        random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                        prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                        verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Rainbow_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_MMR_Rainbow_trout_use <- mcmc_MMR_RR_Rainbow_trout,ifelse(mcmc_MMR_RR_Rainbow_trout$DIC<mcmc_MMR_RI_Rainbow_trout$DIC,mcmc_MMR_Rainbow_trout_use <- mcmc_MMR_RR_Rainbow_trout,mcmc_MMR_Rainbow_trout_use <- mcmc_MMR_RI_Rainbow_trout))
# Modeluse[2,4] <-ifelse(mcmc_MMR_RR_Rainbow_trout$DIC<mcmc_MMR_RI_Rainbow_trout$DIC, "RR","RI")
# Modeluse[2,5] <- abs(mcmc_MMR_RR_Rainbow_trout$DIC-mcmc_MMR_RI_Rainbow_trout$DIC)

#Without treatments
mcmc_MMR_RR_Rainbow_troutM  <- MCMCglmm(logMMR ~ logMassCent  + logMassMean  + Temp,
                                        random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                        prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                        verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Rainbow_trout))

# mcmc_MMR_RI_Rainbow_troutM  <- MCMCglmm(logMMR ~ logMassCent  + logMassMean  + Temp,
#                                         random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                         prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                         verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Rainbow_trout))

ifelse(use_RRM==TRUE,mcmc_MMR_Rainbow_trout_useM <- mcmc_MMR_RR_Rainbow_troutM,ifelse(mcmc_MMR_RR_Rainbow_troutM$DIC<mcmc_MMR_RI_Rainbow_troutM$DIC,mcmc_MMR_Rainbow_trout_useM <- mcmc_MMR_RR_Rainbow_troutM,mcmc_MMR_Rainbow_trout_useM <- mcmc_MMR_RI_Rainbow_troutM))
# ModeluseM[2,4] <-ifelse(mcmc_MMR_RR_Rainbow_troutM$DIC<mcmc_MMR_RI_Rainbow_troutM$DIC, "RR","RI")
# ModeluseM[2,5] <- abs(mcmc_MMR_RR_Rainbow_troutM$DIC-mcmc_MMR_RI_Rainbow_troutM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_MMR_pred_full_Rainbow_trout <- cbind(Data_raw_MMR_Rainbow_trout, predict(mcmc_MMR_Rainbow_trout_use, marginal=NULL, interval = "prediction")),
       df_MMR_pred_full_Rainbow_trout <- cbind(Data_raw_MMR_Rainbow_trout, predict(mcmc_MMR_Rainbow_trout_useM, marginal=NULL, interval = "prediction"))) 


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_MMR_pred_full_Rainbow_trout$MMR_staticM <- df_MMR_pred_full_Rainbow_trout$logMass*summary(mcmc_MMR_Rainbow_trout_useM)$solutions[3,1]+summary(mcmc_MMR_Rainbow_trout_useM)$solutions[1,1]+summary(mcmc_MMR_Rainbow_trout_useM)$solutions[4,1]*df_MMR_pred_full_Rainbow_trout$Temp

# finding the statick scaling value for a given individual at a given point for the given treatment
df_MMR_pred_full_Rainbow_trout[df_MMR_pred_full_Rainbow_trout$Treatment=="10Cegg",]$MMR_staticT <- df_MMR_pred_full_Rainbow_trout[df_MMR_pred_full_Rainbow_trout$Treatment=="10Cegg",]$logMass*summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[6,1]*df_MMR_pred_full_Rainbow_trout[df_MMR_pred_full_Rainbow_trout$Treatment=="10Cegg",]$Temp+summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,1]
df_MMR_pred_full_Rainbow_trout[df_MMR_pred_full_Rainbow_trout$Treatment=="14Cegg",]$MMR_staticT <- df_MMR_pred_full_Rainbow_trout[df_MMR_pred_full_Rainbow_trout$Treatment=="14Cegg",]$logMass*(summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[9,1])+summary(mcmc_MMR_Rainbow_trout_use)$solutions[6,1]*df_MMR_pred_full_Rainbow_trout[df_MMR_pred_full_Rainbow_trout$Treatment=="14Cegg",]$Temp+summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[3,1]
df_MMR_pred_full_Rainbow_trout[df_MMR_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$MMR_staticT <- df_MMR_pred_full_Rainbow_trout[df_MMR_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$logMass*(summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[10,1])+summary(mcmc_MMR_Rainbow_trout_use)$solutions[6,1]*df_MMR_pred_full_Rainbow_trout[df_MMR_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$Temp+summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[4,1]


## Brown trout ##
mcmc_MMR_RR_Brown_trout  <- MCMCglmm(logMMR ~ logMassCent  + logMassMean ,
                                     random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                     prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Brown_trout))

# mcmc_MMR_RI_Brown_trout  <- MCMCglmm(logMMR ~ logMassCent + logMassMean ,
#                                      random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                      prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                      verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Brown_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_MMR_Brown_trout_use <- mcmc_MMR_RR_Brown_trout,ifelse(mcmc_MMR_RR_Brown_trout$DIC<mcmc_MMR_RI_Brown_trout$DIC,mcmc_MMR_Brown_trout_use <- mcmc_MMR_RR_Brown_trout,mcmc_MMR_Brown_trout_use <- mcmc_MMR_RI_Brown_trout))

# Modeluse[3,4] <-ifelse(mcmc_MMR_RR_Brown_trout$DIC<mcmc_MMR_RI_Brown_trout$DIC, "RR","RI")
# Modeluse[3,5] <- abs(mcmc_MMR_RR_Brown_trout$DIC-mcmc_MMR_RI_Brown_trout$DIC)

df_MMR_pred_full_Brown_trout <- cbind(Data_raw_MMR_Brown_trout, predict(mcmc_MMR_Brown_trout_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_MMR_pred_full_Brown_trout$MMR_staticT <- summary(mcmc_MMR_Brown_trout_use)$solutions[1,1] + summary(mcmc_MMR_Brown_trout_use)$solutions[3,1]*df_MMR_pred_full_Brown_trout$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_MMR_pred_full_Brown_trout$MMR_staticM <- summary(mcmc_MMR_Brown_trout_use)$solutions[1,1] + summary(mcmc_MMR_Brown_trout_use)$solutions[3,1]*df_MMR_pred_full_Brown_trout$logMass


## Cunner ##
mcmc_MMR_RR_Cunner  <- MCMCglmm(logMMR ~ logMassCent  + logMassMean ,
                                random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Cunner))

# mcmc_MMR_RI_Cunner  <- MCMCglmm(logMMR ~ logMassCent + logMassMean ,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Cunner))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_MMR_Cunner_use <- mcmc_MMR_RR_Cunner,ifelse(mcmc_MMR_RR_Cunner$DIC<mcmc_MMR_RI_Cunner$DIC,mcmc_MMR_Cunner_use <- mcmc_MMR_RR_Cunner,mcmc_MMR_Cunner_use <- mcmc_MMR_RI_Cunner))
# Modeluse[4,4] <-ifelse(mcmc_MMR_RR_Cunner$DIC<mcmc_MMR_RI_Cunner$DIC, "RR","RI")
# Modeluse[4,5] <- abs(mcmc_MMR_RR_Cunner$DIC-mcmc_MMR_RI_Cunner$DIC)

df_MMR_pred_full_Cunner <- cbind(Data_raw_MMR_Cunner, predict(mcmc_MMR_Cunner_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_MMR_pred_full_Cunner$MMR_staticT <- summary(mcmc_MMR_Cunner_use)$solutions[1,1] + summary(mcmc_MMR_Cunner_use)$solutions[3,1]*df_MMR_pred_full_Cunner$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_MMR_pred_full_Cunner$MMR_staticM <- summary(mcmc_MMR_Cunner_use)$solutions[1,1] + summary(mcmc_MMR_Cunner_use)$solutions[3,1]*df_MMR_pred_full_Cunner$logMass


## Guppy ##
#with sex
mcmc_MMR_RR_Guppy  <- MCMCglmm(logMMR ~ logMassCent * Sex + logMassMean * Sex,
                               random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                               prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Guppy))

# mcmc_MMR_RI_Guppy  <- MCMCglmm(logMMR ~ logMassCent * Sex + logMassMean * Sex,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Guppy))

#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_MMR_Guppy_use <- mcmc_MMR_RR_Guppy,ifelse(mcmc_MMR_RR_Guppy$DIC<mcmc_MMR_RI_Guppy$DIC,mcmc_MMR_Guppy_use <- mcmc_MMR_RR_Guppy,mcmc_MMR_Guppy_use <- mcmc_MMR_RI_Guppy))
# Modeluse[5,4] <-ifelse(mcmc_MMR_RR_Guppy$DIC<mcmc_MMR_RI_Guppy$DIC, "RR","RI")
# Modeluse[5,5] <- abs(mcmc_MMR_RR_Guppy$DIC-mcmc_MMR_RI_Guppy$DIC)


#without sex
mcmc_MMR_RR_GuppyM  <- MCMCglmm(logMMR ~ logMassCent  + logMassMean ,
                                random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Guppy))

# mcmc_MMR_RI_GuppyM  <- MCMCglmm(logMMR ~ logMassCent  + logMassMean ,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Guppy))

ifelse(use_RRM==TRUE,mcmc_MMR_Guppy_useM <- mcmc_MMR_RR_GuppyM,ifelse(mcmc_MMR_RR_GuppyM$DIC<mcmc_MMR_RI_GuppyM$DIC,mcmc_MMR_Guppy_useM <- mcmc_MMR_RR_GuppyM,mcmc_MMR_Guppy_useM <- mcmc_MMR_RI_GuppyM))
# ModeluseM[5,4] <-ifelse(mcmc_MMR_RR_GuppyM$DIC<mcmc_MMR_RI_GuppyM$DIC, "RR","RI")
# ModeluseM[5,5] <- abs(mcmc_MMR_RR_GuppyM$DIC-mcmc_MMR_RI_GuppyM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_MMR_pred_full_Guppy <- cbind(Data_raw_MMR_Guppy, predict(mcmc_MMR_Guppy_use, marginal=NULL, interval = "prediction")),
       df_MMR_pred_full_Guppy <- cbind(Data_raw_MMR_Guppy, predict(mcmc_MMR_Guppy_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_MMR_pred_full_Guppy[df_MMR_pred_full_Guppy$Sex=="f",]$MMR_staticT <- summary(mcmc_MMR_Guppy_use)$solutions[1,1] + summary(mcmc_MMR_Guppy_use)$solutions[4,1]*df_MMR_pred_full_Guppy[df_MMR_pred_full_Guppy$Sex=="f",]$logMass
df_MMR_pred_full_Guppy[df_MMR_pred_full_Guppy$Sex=="m",]$MMR_staticT <- summary(mcmc_MMR_Guppy_use)$solutions[1,1]+summary(mcmc_MMR_Guppy_use)$solutions[3,1] + (summary(mcmc_MMR_Guppy_use)$solutions[4,1]+summary(mcmc_MMR_Guppy_use)$solutions[6,1])*df_MMR_pred_full_Guppy[df_MMR_pred_full_Guppy$Sex=="m",]$logMass

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_MMR_pred_full_Guppy$MMR_staticM <- summary(mcmc_MMR_Guppy_useM)$solutions[1,1] + summary(mcmc_MMR_Guppy_useM)$solutions[3,1]*df_MMR_pred_full_Guppy$logMass


## Chromis ##
#With treatments
mcmc_MMR_RR_Chromis  <- MCMCglmm(logMMR ~ logMassCent * Treatment + logMassMean * Treatment ,
                                 random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                 prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Chromis))

# mcmc_MMR_RI_Chromis  <- MCMCglmm(logMMR ~ logMassCent * Treatment + logMassMean * Treatment,
#                                  random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                  prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Chromis))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_MMR_Chromis_use <- mcmc_MMR_RR_Chromis,ifelse(mcmc_MMR_RR_Chromis$DIC<mcmc_MMR_RI_Chromis$DIC,mcmc_MMR_Chromis_use <- mcmc_MMR_RR_Chromis,mcmc_MMR_Chromis_use <- mcmc_MMR_RI_Chromis))
# Modeluse[6,4] <-ifelse(mcmc_MMR_RR_Chromis$DIC<mcmc_MMR_RI_Chromis$DIC, "RR","RI")
# Modeluse[6,5] <- abs(mcmc_MMR_RR_Chromis$DIC-mcmc_MMR_RI_Chromis$DIC)

#Without treatments
mcmc_MMR_RR_ChromisM  <- MCMCglmm(logMMR ~ logMassCent  + logMassMean ,
                                  random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                  prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Chromis))

# mcmc_MMR_RI_ChromisM  <- MCMCglmm(logMMR ~ logMassCent  + logMassMean,
#                                   random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                   prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Chromis))

ifelse(use_RRM==TRUE,mcmc_MMR_Chromis_useM <- mcmc_MMR_RR_ChromisM,ifelse(mcmc_MMR_RR_ChromisM$DIC<mcmc_MMR_RI_ChromisM$DIC,mcmc_MMR_Chromis_useM <- mcmc_MMR_RR_ChromisM,mcmc_MMR_Chromis_useM <- mcmc_MMR_RI_ChromisM))
# ModeluseM[6,4] <-ifelse(mcmc_MMR_RR_ChromisM$DIC<mcmc_MMR_RI_ChromisM$DIC, "RR","RI")
# ModeluseM[6,5] <- abs(mcmc_MMR_RR_ChromisM$DIC-mcmc_MMR_RI_ChromisM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_MMR_pred_full_Chromis <- cbind(Data_raw_MMR_Chromis, predict(mcmc_MMR_Chromis_use, marginal=NULL, interval = "prediction")),
       df_MMR_pred_full_Chromis <- cbind(Data_raw_MMR_Chromis, predict(mcmc_MMR_Chromis_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_MMR_pred_full_Chromis[df_MMR_pred_full_Chromis$Treatment=="high",]$MMR_staticT <- summary(mcmc_MMR_Chromis_use)$solutions[1,1] + summary(mcmc_MMR_Chromis_use)$solutions[4,1]*df_MMR_pred_full_Chromis[df_MMR_pred_full_Chromis$Treatment=="high",]$logMass
df_MMR_pred_full_Chromis[df_MMR_pred_full_Chromis$Treatment=="low",]$MMR_staticT <- summary(mcmc_MMR_Chromis_use)$solutions[1,1]+summary(mcmc_MMR_Chromis_use)$solutions[3,1] + (summary(mcmc_MMR_Chromis_use)$solutions[4,1]+summary(mcmc_MMR_Chromis_use)$solutions[6,1])*df_MMR_pred_full_Chromis[df_MMR_pred_full_Chromis$Treatment=="low",]$logMass
#Finding the Avg statick scaling for a given species, disregaringing treatment
df_MMR_pred_full_Chromis$MMR_staticM <- summary(mcmc_MMR_Chromis_useM)$solutions[1,1] + summary(mcmc_MMR_Chromis_useM)$solutions[3,1]*df_MMR_pred_full_Chromis$logMass



## Clownfish ##
mcmc_MMR_RR_Clownfish  <- MCMCglmm(logMMR ~ logMassCent  + logMassMean ,
                                   random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Clownfish))

# mcmc_MMR_RI_Clownfish  <- MCMCglmm(logMMR ~ logMassCent + logMassMean ,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_MMR_Clownfish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_MMR_Clownfish_use <- mcmc_MMR_RR_Clownfish,ifelse(mcmc_MMR_RR_Clownfish$DIC<mcmc_MMR_RI_Clownfish$DIC,mcmc_MMR_Clownfish_use <- mcmc_MMR_RR_Clownfish,mcmc_MMR_Clownfish_use <- mcmc_MMR_RI_Clownfish))
# Modeluse[7,4] <-ifelse(mcmc_MMR_RR_Clownfish$DIC<mcmc_MMR_RI_Clownfish$DIC, "RR","RI")
# Modeluse[7,5] <- abs(mcmc_MMR_RR_Clownfish$DIC-mcmc_MMR_RI_Clownfish$DIC)

df_MMR_pred_full_Clownfish <- cbind(Data_raw_MMR_Clownfish, predict(mcmc_MMR_Clownfish_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_MMR_pred_full_Clownfish$MMR_staticT <- summary(mcmc_MMR_Clownfish_use)$solutions[1,1] + summary(mcmc_MMR_Clownfish_use)$solutions[3,1]*df_MMR_pred_full_Clownfish$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_MMR_pred_full_Clownfish$MMR_staticM <- summary(mcmc_MMR_Clownfish_use)$solutions[1,1] + summary(mcmc_MMR_Clownfish_use)$solutions[3,1]*df_MMR_pred_full_Clownfish$logMass


# Recombining #
Data_MMR <- bind_rows( df_MMR_pred_full_Rainbow_trout, df_MMR_pred_full_Brown_trout, df_MMR_pred_full_Cunner, df_MMR_pred_full_Guppy, df_MMR_pred_full_Chromis, df_MMR_pred_full_Clownfish)

names(Data_MMR)[names(Data_MMR) == "fit" ] <- "MMR_fit"
names(Data_MMR)[names(Data_MMR) == "lwr" ] <- "MMR_fit_low"
names(Data_MMR)[names(Data_MMR) == "upr" ] <- "MMR_fit_high"

### FAS ###
#FAS dataframe
Data_raw_FAS <- Data_raw
Data_raw_FAS$FAS <- Data_raw_FAS$MMR/Data_raw_FAS$SMR
Data_raw_FAS$logFAS <- log10(Data_raw_FAS$FAS)
Data_raw_FAS$FAS_staticM <- NA
Data_raw_FAS$FAS_staticT <- NA

#Subsetting data
Data_raw_FAS_Rainbow_trout <- subset(Data_raw_FAS, Species=="Rainbow_trout")
Data_raw_FAS_Brown_trout <- subset(Data_raw_FAS, Species=="Brown_trout")
Data_raw_FAS_Cunner <- subset(Data_raw_FAS, Species=="Cunner")
Data_raw_FAS_Guppy <- subset(Data_raw_FAS, Species=="Guppy", Sex != "n" & Sex != "d")
Data_raw_FAS_Chromis <- subset(Data_raw_FAS, Species=="Chromis")
Data_raw_FAS_Clownfish <- subset(Data_raw_FAS, Species=="Clownfish")

## Rainbow trout ##
#With treatments
mcmc_FAS_RR_Rainbow_trout  <- MCMCglmm(logFAS ~ logMassCent * Treatment + logMassMean * Treatment + Temp,
                                      random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                      prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                      verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Rainbow_trout))

# mcmc_FAS_RI_Rainbow_trout  <- MCMCglmm(logFAS ~ logMassCent * Treatment + logMassMean * Treatment + Temp,
#                                       random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                       prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                       verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Rainbow_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_FAS_Rainbow_trout_use <- mcmc_FAS_RR_Rainbow_trout,ifelse(mcmc_FAS_RR_Rainbow_trout$DIC<mcmc_FAS_RI_Rainbow_trout$DIC,mcmc_FAS_Rainbow_trout_use <- mcmc_FAS_RR_Rainbow_trout,mcmc_FAS_Rainbow_trout_use <- mcmc_FAS_RI_Rainbow_trout))

# Modeluse[2,6] <-ifelse(mcmc_FAS_RR_Rainbow_trout$DIC<mcmc_FAS_RI_Rainbow_trout$DIC, "RR","RI")
# Modeluse[2,7] <- abs(mcmc_FAS_RR_Rainbow_trout$DIC-mcmc_FAS_RI_Rainbow_trout$DIC)

#Without treatments
mcmc_FAS_RR_Rainbow_troutM  <- MCMCglmm(logFAS ~ logMassCent  + logMassMean  + Temp,
                                       random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                       prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                       verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Rainbow_trout))

# mcmc_FAS_RI_Rainbow_troutM  <- MCMCglmm(logFAS ~ logMassCent  + logMassMean  + Temp,
#                                        random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                        prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                        verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Rainbow_trout))

ifelse(use_RRM==TRUE,mcmc_FAS_Rainbow_trout_useM <- mcmc_FAS_RR_Rainbow_troutM,ifelse(mcmc_FAS_RR_Rainbow_troutM$DIC<mcmc_FAS_RI_Rainbow_troutM$DIC,mcmc_FAS_Rainbow_trout_useM <- mcmc_FAS_RR_Rainbow_troutM,mcmc_FAS_Rainbow_trout_useM <- mcmc_FAS_RI_Rainbow_troutM))
# ModeluseM[2,6] <-ifelse(mcmc_FAS_RR_Rainbow_troutM$DIC<mcmc_FAS_RI_Rainbow_troutM$DIC, "RR","RI")
# ModeluseM[2,7] <- abs(mcmc_FAS_RR_Rainbow_troutM$DIC-mcmc_FAS_RI_Rainbow_troutM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_FAS_pred_full_Rainbow_trout <- cbind(Data_raw_FAS_Rainbow_trout, predict(mcmc_FAS_Rainbow_trout_use, marginal=NULL, interval = "prediction")),
       df_FAS_pred_full_Rainbow_trout <- cbind(Data_raw_FAS_Rainbow_trout, predict(mcmc_FAS_Rainbow_trout_useM, marginal=NULL, interval = "prediction"))) 


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_FAS_pred_full_Rainbow_trout$FAS_staticM <- df_FAS_pred_full_Rainbow_trout$logMass*summary(mcmc_FAS_Rainbow_trout_useM)$solutions[3,1]+summary(mcmc_FAS_Rainbow_trout_useM)$solutions[1,1]+summary(mcmc_FAS_Rainbow_trout_useM)$solutions[4,1]*df_FAS_pred_full_Rainbow_trout$Temp

# finding the statick scaling value for a given individual at a given point for the given treatment
df_FAS_pred_full_Rainbow_trout[df_FAS_pred_full_Rainbow_trout$Treatment=="10Cegg",]$FAS_staticT <- df_FAS_pred_full_Rainbow_trout[df_FAS_pred_full_Rainbow_trout$Treatment=="10Cegg",]$logMass*summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[6,1]*df_FAS_pred_full_Rainbow_trout[df_FAS_pred_full_Rainbow_trout$Treatment=="10Cegg",]$Temp+summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,1]
df_FAS_pred_full_Rainbow_trout[df_FAS_pred_full_Rainbow_trout$Treatment=="14Cegg",]$FAS_staticT <- df_FAS_pred_full_Rainbow_trout[df_FAS_pred_full_Rainbow_trout$Treatment=="14Cegg",]$logMass*(summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[9,1])+summary(mcmc_FAS_Rainbow_trout_use)$solutions[6,1]*df_FAS_pred_full_Rainbow_trout[df_FAS_pred_full_Rainbow_trout$Treatment=="14Cegg",]$Temp+summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[3,1]
df_FAS_pred_full_Rainbow_trout[df_FAS_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$FAS_staticT <- df_FAS_pred_full_Rainbow_trout[df_FAS_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$logMass*(summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[10,1])+summary(mcmc_FAS_Rainbow_trout_use)$solutions[6,1]*df_FAS_pred_full_Rainbow_trout[df_FAS_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$Temp+summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[4,1]


## Brown trout ##
mcmc_FAS_RR_Brown_trout  <- MCMCglmm(logFAS ~ logMassCent  + logMassMean ,
                                    random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                    prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Brown_trout))

# mcmc_FAS_RI_Brown_trout  <- MCMCglmm(logFAS ~ logMassCent + logMassMean ,
#                                     random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                     prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Brown_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_FAS_Brown_trout_use <- mcmc_FAS_RR_Brown_trout,ifelse(mcmc_FAS_RR_Brown_trout$DIC<mcmc_FAS_RI_Brown_trout$DIC,mcmc_FAS_Brown_trout_use <- mcmc_FAS_RR_Brown_trout,mcmc_FAS_Brown_trout_use <- mcmc_FAS_RI_Brown_trout))
# Modeluse[3,6] <-ifelse(mcmc_FAS_RR_Brown_trout$DIC<mcmc_FAS_RI_Brown_trout$DIC, "RR","RI")
# Modeluse[3,7] <- abs(mcmc_FAS_RR_Brown_trout$DIC-mcmc_FAS_RI_Brown_trout$DIC)

df_FAS_pred_full_Brown_trout <- cbind(Data_raw_FAS_Brown_trout, predict(mcmc_FAS_Brown_trout_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_FAS_pred_full_Brown_trout$FAS_staticT <- summary(mcmc_FAS_Brown_trout_use)$solutions[1,1] + summary(mcmc_FAS_Brown_trout_use)$solutions[3,1]*df_FAS_pred_full_Brown_trout$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_FAS_pred_full_Brown_trout$FAS_staticM <- summary(mcmc_FAS_Brown_trout_use)$solutions[1,1] + summary(mcmc_FAS_Brown_trout_use)$solutions[3,1]*df_FAS_pred_full_Brown_trout$logMass


# Cunner #
mcmc_FAS_RR_Cunner  <- MCMCglmm(logFAS ~ logMassCent  + logMassMean ,
                               random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                               prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Cunner))

# mcmc_FAS_RI_Cunner  <- MCMCglmm(logFAS ~ logMassCent + logMassMean ,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Cunner))
# #Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_FAS_Cunner_use <- mcmc_FAS_RR_Cunner,ifelse(mcmc_FAS_RR_Cunner$DIC<mcmc_FAS_RI_Cunner$DIC,mcmc_FAS_Cunner_use <- mcmc_FAS_RR_Cunner,mcmc_FAS_Cunner_use <- mcmc_FAS_RI_Cunner))
# Modeluse[4,6] <-ifelse(mcmc_FAS_RR_Cunner$DIC<mcmc_FAS_RI_Cunner$DIC, "RR","RI")
# Modeluse[4,7] <- abs(mcmc_FAS_RR_Cunner$DIC-mcmc_FAS_RI_Cunner$DIC)

df_FAS_pred_full_Cunner <- cbind(Data_raw_FAS_Cunner, predict(mcmc_FAS_Cunner_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_FAS_pred_full_Cunner$FAS_staticT <- summary(mcmc_FAS_Cunner_use)$solutions[1,1] + summary(mcmc_FAS_Cunner_use)$solutions[3,1]*df_FAS_pred_full_Cunner$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_FAS_pred_full_Cunner$FAS_staticM <- summary(mcmc_FAS_Cunner_use)$solutions[1,1] + summary(mcmc_FAS_Cunner_use)$solutions[3,1]*df_FAS_pred_full_Cunner$logMass


## Guppy ##
#with sex
mcmc_FAS_RR_Guppy  <- MCMCglmm(logFAS ~ logMassCent * Sex + logMassMean * Sex,
                              random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                              prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                              verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Guppy))

# mcmc_FAS_RI_Guppy  <- MCMCglmm(logFAS ~ logMassCent * Sex + logMassMean * Sex,
#                               random = ~ FishID, rcov = ~ units, family = "gaussian",
#                               prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Guppy))

#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_FAS_Guppy_use <- mcmc_FAS_RR_Guppy,ifelse(mcmc_FAS_RR_Guppy$DIC<mcmc_FAS_RI_Guppy$DIC,mcmc_FAS_Guppy_use <- mcmc_FAS_RR_Guppy,mcmc_FAS_Guppy_use <- mcmc_FAS_RI_Guppy))
# Modeluse[5,6] <-ifelse(mcmc_FAS_RR_Guppy$DIC<mcmc_FAS_RI_Guppy$DIC, "RR","RI")
# Modeluse[5,7] <- abs(mcmc_FAS_RR_Guppy$DIC-mcmc_FAS_RI_Guppy$DIC)


#without sex
mcmc_FAS_RR_GuppyM  <- MCMCglmm(logFAS ~ logMassCent  + logMassMean ,
                               random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                               prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Guppy))

# mcmc_FAS_RI_GuppyM  <- MCMCglmm(logFAS ~ logMassCent  + logMassMean ,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Guppy))

ifelse(use_RRM==TRUE,mcmc_FAS_Guppy_useM <- mcmc_FAS_RR_GuppyM,ifelse(mcmc_FAS_RR_GuppyM$DIC<mcmc_FAS_RI_GuppyM$DIC,mcmc_FAS_Guppy_useM <- mcmc_FAS_RR_GuppyM,mcmc_FAS_Guppy_useM <- mcmc_FAS_RI_GuppyM))
# ModeluseM[5,6] <-ifelse(mcmc_FAS_RR_GuppyM$DIC<mcmc_FAS_RI_GuppyM$DIC, "RR","RI")
# ModeluseM[5,7] <- abs(mcmc_FAS_RR_GuppyM$DIC-mcmc_FAS_RI_GuppyM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_FAS_pred_full_Guppy <- cbind(Data_raw_FAS_Guppy, predict(mcmc_FAS_Guppy_use, marginal=NULL, interval = "prediction")),
       df_FAS_pred_full_Guppy <- cbind(Data_raw_FAS_Guppy, predict(mcmc_FAS_Guppy_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_FAS_pred_full_Guppy[df_FAS_pred_full_Guppy$Sex=="f",]$FAS_staticT <- summary(mcmc_FAS_Guppy_use)$solutions[1,1] + summary(mcmc_FAS_Guppy_use)$solutions[4,1]*df_FAS_pred_full_Guppy[df_FAS_pred_full_Guppy$Sex=="f",]$logMass
df_FAS_pred_full_Guppy[df_FAS_pred_full_Guppy$Sex=="m",]$FAS_staticT <- summary(mcmc_FAS_Guppy_use)$solutions[1,1]+summary(mcmc_FAS_Guppy_use)$solutions[3,1] + (summary(mcmc_FAS_Guppy_use)$solutions[4,1]+summary(mcmc_FAS_Guppy_use)$solutions[6,1])*df_FAS_pred_full_Guppy[df_FAS_pred_full_Guppy$Sex=="m",]$logMass

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_FAS_pred_full_Guppy$FAS_staticM <- summary(mcmc_FAS_Guppy_useM)$solutions[1,1] + summary(mcmc_FAS_Guppy_useM)$solutions[3,1]*df_FAS_pred_full_Guppy$logMass


## Chromis ##
#With treatments
mcmc_FAS_RR_Chromis  <- MCMCglmm(logFAS ~ logMassCent * Treatment + logMassMean * Treatment ,
                                random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Chromis))

# mcmc_FAS_RI_Chromis  <- MCMCglmm(logFAS ~ logMassCent * Treatment + logMassMean * Treatment,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Chromis))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_FAS_Chromis_use <- mcmc_FAS_RR_Chromis,ifelse(mcmc_FAS_RR_Chromis$DIC<mcmc_FAS_RI_Chromis$DIC,mcmc_FAS_Chromis_use <- mcmc_FAS_RR_Chromis,mcmc_FAS_Chromis_use <- mcmc_FAS_RI_Chromis))
# Modeluse[6,6] <-ifelse(mcmc_FAS_RR_Chromis$DIC<mcmc_FAS_RI_Chromis$DIC, "RR","RI")
# Modeluse[6,7] <- abs(mcmc_FAS_RR_Chromis$DIC-mcmc_FAS_RI_Chromis$DIC)

#Without treatments
mcmc_FAS_RR_ChromisM  <- MCMCglmm(logFAS ~ logMassCent  + logMassMean ,
                                 random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                 prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Chromis))

# mcmc_FAS_RI_ChromisM  <- MCMCglmm(logFAS ~ logMassCent  + logMassMean,
#                                  random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                  prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Chromis))

ifelse(use_RRM==TRUE,mcmc_FAS_Chromis_useM <- mcmc_FAS_RR_ChromisM,ifelse(mcmc_FAS_RR_ChromisM$DIC<mcmc_FAS_RI_ChromisM$DIC,mcmc_FAS_Chromis_useM <- mcmc_FAS_RR_ChromisM,mcmc_FAS_Chromis_useM <- mcmc_FAS_RI_ChromisM))
# ModeluseM[6,6] <-ifelse(mcmc_FAS_RR_ChromisM$DIC<mcmc_FAS_RI_ChromisM$DIC, "RR","RI")
# ModeluseM[6,7] <- abs(mcmc_FAS_RR_ChromisM$DIC-mcmc_FAS_RI_ChromisM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_FAS_pred_full_Chromis <- cbind(Data_raw_FAS_Chromis, predict(mcmc_FAS_Chromis_use, marginal=NULL, interval = "prediction")),
       df_FAS_pred_full_Chromis <- cbind(Data_raw_FAS_Chromis, predict(mcmc_FAS_Chromis_useM, marginal=NULL, interval = "prediction"))) 



# finding the statick scaling value for a given individual at a given point for the given treatment
df_FAS_pred_full_Chromis[df_FAS_pred_full_Chromis$Treatment=="high",]$FAS_staticT <- summary(mcmc_FAS_Chromis_use)$solutions[1,1] + summary(mcmc_FAS_Chromis_use)$solutions[4,1]*df_FAS_pred_full_Chromis[df_FAS_pred_full_Chromis$Treatment=="high",]$logMass
df_FAS_pred_full_Chromis[df_FAS_pred_full_Chromis$Treatment=="low",]$FAS_staticT <- summary(mcmc_FAS_Chromis_use)$solutions[1,1]+summary(mcmc_FAS_Chromis_use)$solutions[3,1] + (summary(mcmc_FAS_Chromis_use)$solutions[4,1]+summary(mcmc_FAS_Chromis_use)$solutions[6,1])*df_FAS_pred_full_Chromis[df_FAS_pred_full_Chromis$Treatment=="low",]$logMass
#Finding the Avg statick scaling for a given species, disregaringing treatment
df_FAS_pred_full_Chromis$FAS_staticM <- summary(mcmc_FAS_Chromis_useM)$solutions[1,1] + summary(mcmc_FAS_Chromis_useM)$solutions[3,1]*df_FAS_pred_full_Chromis$logMass



## Clownfish ##
mcmc_FAS_RR_Clownfish  <- MCMCglmm(logFAS ~ logMassCent  + logMassMean ,
                                  random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                  prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Clownfish))

# mcmc_FAS_RI_Clownfish  <- MCMCglmm(logFAS ~ logMassCent + logMassMean ,
#                                   random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                   prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_FAS_Clownfish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_FAS_Clownfish_use <- mcmc_FAS_RR_Clownfish,ifelse(mcmc_FAS_RR_Clownfish$DIC<mcmc_FAS_RI_Clownfish$DIC,mcmc_FAS_Clownfish_use <- mcmc_FAS_RR_Clownfish,mcmc_FAS_Clownfish_use <- mcmc_FAS_RI_Clownfish))
# Modeluse[7,6] <-ifelse(mcmc_FAS_RR_Clownfish$DIC<mcmc_FAS_RI_Clownfish$DIC, "RR","RI")
# Modeluse[7,7] <- abs(mcmc_FAS_RR_Clownfish$DIC-mcmc_FAS_RI_Clownfish$DIC)

df_FAS_pred_full_Clownfish <- cbind(Data_raw_FAS_Clownfish, predict(mcmc_FAS_Clownfish_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_FAS_pred_full_Clownfish$FAS_staticT <- summary(mcmc_FAS_Clownfish_use)$solutions[1,1] + summary(mcmc_FAS_Clownfish_use)$solutions[3,1]*df_FAS_pred_full_Clownfish$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_FAS_pred_full_Clownfish$FAS_staticM <- summary(mcmc_FAS_Clownfish_use)$solutions[1,1] + summary(mcmc_FAS_Clownfish_use)$solutions[3,1]*df_FAS_pred_full_Clownfish$logMass


# Recombining #
Data_FAS <- bind_rows(df_FAS_pred_full_Rainbow_trout, df_FAS_pred_full_Brown_trout, df_FAS_pred_full_Cunner, df_FAS_pred_full_Guppy, df_FAS_pred_full_Chromis, df_FAS_pred_full_Clownfish)

names(Data_FAS)[names(Data_FAS) == "fit" ] <- "FAS_fit"
names(Data_FAS)[names(Data_FAS) == "lwr" ] <- "FAS_fit_low"
names(Data_FAS)[names(Data_FAS) == "upr" ] <- "FAS_fit_high"



### AS ###
#AS dataframe
Data_raw_AS <- Data_raw
Data_raw_AS$AS_staticM <- NA
Data_raw_AS$AS_staticT <- NA

#Subsetting data
Data_raw_AS_Rainbow_trout <- subset(Data_raw_AS, Species=="Rainbow_trout")
Data_raw_AS_Brown_trout <- subset(Data_raw_AS, Species=="Brown_trout")
Data_raw_AS_Cunner <- subset(Data_raw_AS, Species=="Cunner")
Data_raw_AS_Guppy <- subset(Data_raw_AS, Species=="Guppy", Sex != "n" & Sex != "d")
Data_raw_AS_Chromis <- subset(Data_raw_AS, Species=="Chromis")
Data_raw_AS_Clownfish <- subset(Data_raw_AS, Species=="Clownfish")

## Rainbow trout ##
#With treatments
mcmc_AS_RR_Rainbow_trout  <- MCMCglmm(logAS ~ logMassCent * Treatment + logMassMean * Treatment + Temp,
                                       random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                       prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                       verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Rainbow_trout))

# mcmc_AS_RI_Rainbow_trout  <- MCMCglmm(logAS ~ logMassCent * Treatment + logMassMean * Treatment + Temp,
#                                        random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                        prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                        verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Rainbow_trout))
# #Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AS_Rainbow_trout_use <- mcmc_AS_RR_Rainbow_trout,ifelse(mcmc_AS_RR_Rainbow_trout$DIC<mcmc_AS_RI_Rainbow_trout$DIC,mcmc_AS_Rainbow_trout_use <- mcmc_AS_RR_Rainbow_trout,mcmc_AS_Rainbow_trout_use <- mcmc_AS_RI_Rainbow_trout))

# Modeluse[2,6] <-ifelse(mcmc_AS_RR_Rainbow_trout$DIC<mcmc_AS_RI_Rainbow_trout$DIC, "RR","RI")
# Modeluse[2,7] <- abs(mcmc_AS_RR_Rainbow_trout$DIC-mcmc_AS_RI_Rainbow_trout$DIC)

#Without treatments
mcmc_AS_RR_Rainbow_troutM  <- MCMCglmm(logAS ~ logMassCent  + logMassMean  + Temp,
                                        random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                        prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                        verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Rainbow_trout))

# mcmc_AS_RI_Rainbow_troutM  <- MCMCglmm(logAS ~ logMassCent  + logMassMean  + Temp,
#                                         random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                         prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                         verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Rainbow_trout))

ifelse(use_RRM==TRUE,mcmc_AS_Rainbow_trout_useM <- mcmc_AS_RR_Rainbow_troutM,ifelse(mcmc_AS_RR_Rainbow_troutM$DIC<mcmc_AS_RI_Rainbow_troutM$DIC,mcmc_AS_Rainbow_trout_useM <- mcmc_AS_RR_Rainbow_troutM,mcmc_AS_Rainbow_trout_useM <- mcmc_AS_RI_Rainbow_troutM))
# ModeluseM[2,6] <-ifelse(mcmc_AS_RR_Rainbow_troutM$DIC<mcmc_AS_RI_Rainbow_troutM$DIC, "RR","RI")
# ModeluseM[2,7] <- abs(mcmc_AS_RR_Rainbow_troutM$DIC-mcmc_AS_RI_Rainbow_troutM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_AS_pred_full_Rainbow_trout <- cbind(Data_raw_AS_Rainbow_trout, predict(mcmc_AS_Rainbow_trout_use, marginal=NULL, interval = "prediction")),
       df_AS_pred_full_Rainbow_trout <- cbind(Data_raw_AS_Rainbow_trout, predict(mcmc_AS_Rainbow_trout_useM, marginal=NULL, interval = "prediction"))) 


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AS_pred_full_Rainbow_trout$AS_staticM <- df_AS_pred_full_Rainbow_trout$logMass*summary(mcmc_AS_Rainbow_trout_useM)$solutions[3,1]+summary(mcmc_AS_Rainbow_trout_useM)$solutions[1,1]+summary(mcmc_AS_Rainbow_trout_useM)$solutions[4,1]*df_AS_pred_full_Rainbow_trout$Temp

# finding the statick scaling value for a given individual at a given point for the given treatment
df_AS_pred_full_Rainbow_trout[df_AS_pred_full_Rainbow_trout$Treatment=="10Cegg",]$AS_staticT <- df_AS_pred_full_Rainbow_trout[df_AS_pred_full_Rainbow_trout$Treatment=="10Cegg",]$logMass*summary(mcmc_AS_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[6,1]*df_AS_pred_full_Rainbow_trout[df_AS_pred_full_Rainbow_trout$Treatment=="10Cegg",]$Temp+summary(mcmc_AS_Rainbow_trout_use)$solutions[1,1]
df_AS_pred_full_Rainbow_trout[df_AS_pred_full_Rainbow_trout$Treatment=="14Cegg",]$AS_staticT <- df_AS_pred_full_Rainbow_trout[df_AS_pred_full_Rainbow_trout$Treatment=="14Cegg",]$logMass*(summary(mcmc_AS_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[9,1])+summary(mcmc_AS_Rainbow_trout_use)$solutions[6,1]*df_AS_pred_full_Rainbow_trout[df_AS_pred_full_Rainbow_trout$Treatment=="14Cegg",]$Temp+summary(mcmc_AS_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[3,1]
df_AS_pred_full_Rainbow_trout[df_AS_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$AS_staticT <- df_AS_pred_full_Rainbow_trout[df_AS_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$logMass*(summary(mcmc_AS_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[10,1])+summary(mcmc_AS_Rainbow_trout_use)$solutions[6,1]*df_AS_pred_full_Rainbow_trout[df_AS_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$Temp+summary(mcmc_AS_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[4,1]


## Brown trout ##
mcmc_AS_RR_Brown_trout  <- MCMCglmm(logAS ~ logMassCent  + logMassMean ,
                                     random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                     prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Brown_trout))

# mcmc_AS_RI_Brown_trout  <- MCMCglmm(logAS ~ logMassCent + logMassMean ,
#                                      random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                      prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                      verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Brown_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AS_Brown_trout_use <- mcmc_AS_RR_Brown_trout,ifelse(mcmc_AS_RR_Brown_trout$DIC<mcmc_AS_RI_Brown_trout$DIC,mcmc_AS_Brown_trout_use <- mcmc_AS_RR_Brown_trout,mcmc_AS_Brown_trout_use <- mcmc_AS_RI_Brown_trout))
# Modeluse[3,6] <-ifelse(mcmc_AS_RR_Brown_trout$DIC<mcmc_AS_RI_Brown_trout$DIC, "RR","RI")
# Modeluse[3,7] <- abs(mcmc_AS_RR_Brown_trout$DIC-mcmc_AS_RI_Brown_trout$DIC)

df_AS_pred_full_Brown_trout <- cbind(Data_raw_AS_Brown_trout, predict(mcmc_AS_Brown_trout_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_AS_pred_full_Brown_trout$AS_staticT <- summary(mcmc_AS_Brown_trout_use)$solutions[1,1] + summary(mcmc_AS_Brown_trout_use)$solutions[3,1]*df_AS_pred_full_Brown_trout$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AS_pred_full_Brown_trout$AS_staticM <- summary(mcmc_AS_Brown_trout_use)$solutions[1,1] + summary(mcmc_AS_Brown_trout_use)$solutions[3,1]*df_AS_pred_full_Brown_trout$logMass


# Cunner #
mcmc_AS_RR_Cunner  <- MCMCglmm(logAS ~ logMassCent  + logMassMean ,
                                random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Cunner))

# mcmc_AS_RI_Cunner  <- MCMCglmm(logAS ~ logMassCent + logMassMean ,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Cunner))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AS_Cunner_use <- mcmc_AS_RR_Cunner,ifelse(mcmc_AS_RR_Cunner$DIC<mcmc_AS_RI_Cunner$DIC,mcmc_AS_Cunner_use <- mcmc_AS_RR_Cunner,mcmc_AS_Cunner_use <- mcmc_AS_RI_Cunner))
# Modeluse[4,6] <-ifelse(mcmc_AS_RR_Cunner$DIC<mcmc_AS_RI_Cunner$DIC, "RR","RI")
# Modeluse[4,7] <- abs(mcmc_AS_RR_Cunner$DIC-mcmc_AS_RI_Cunner$DIC)

df_AS_pred_full_Cunner <- cbind(Data_raw_AS_Cunner, predict(mcmc_AS_Cunner_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_AS_pred_full_Cunner$AS_staticT <- summary(mcmc_AS_Cunner_use)$solutions[1,1] + summary(mcmc_AS_Cunner_use)$solutions[3,1]*df_AS_pred_full_Cunner$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AS_pred_full_Cunner$AS_staticM <- summary(mcmc_AS_Cunner_use)$solutions[1,1] + summary(mcmc_AS_Cunner_use)$solutions[3,1]*df_AS_pred_full_Cunner$logMass


## Guppy ##
#with sex
mcmc_AS_RR_Guppy  <- MCMCglmm(logAS ~ logMassCent * Sex + logMassMean * Sex,
                               random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                               prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Guppy))

# mcmc_AS_RI_Guppy  <- MCMCglmm(logAS ~ logMassCent * Sex + logMassMean * Sex,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Guppy))

#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AS_Guppy_use <- mcmc_AS_RR_Guppy,ifelse(mcmc_AS_RR_Guppy$DIC<mcmc_AS_RI_Guppy$DIC,mcmc_AS_Guppy_use <- mcmc_AS_RR_Guppy,mcmc_AS_Guppy_use <- mcmc_AS_RI_Guppy))
# Modeluse[5,6] <-ifelse(mcmc_AS_RR_Guppy$DIC<mcmc_AS_RI_Guppy$DIC, "RR","RI")
# Modeluse[5,7] <- abs(mcmc_AS_RR_Guppy$DIC-mcmc_AS_RI_Guppy$DIC)


#without sex
mcmc_AS_RR_GuppyM  <- MCMCglmm(logAS ~ logMassCent  + logMassMean ,
                                random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Guppy))

# mcmc_AS_RI_GuppyM  <- MCMCglmm(logAS ~ logMassCent  + logMassMean ,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Guppy))

ifelse(use_RRM==TRUE,mcmc_AS_Guppy_useM <- mcmc_AS_RR_GuppyM,ifelse(mcmc_AS_RR_GuppyM$DIC<mcmc_AS_RI_GuppyM$DIC,mcmc_AS_Guppy_useM <- mcmc_AS_RR_GuppyM,mcmc_AS_Guppy_useM <- mcmc_AS_RI_GuppyM))
# ModeluseM[5,6] <-ifelse(mcmc_AS_RR_GuppyM$DIC<mcmc_AS_RI_GuppyM$DIC, "RR","RI")
# ModeluseM[5,7] <- abs(mcmc_AS_RR_GuppyM$DIC-mcmc_AS_RI_GuppyM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_AS_pred_full_Guppy <- cbind(Data_raw_AS_Guppy, predict(mcmc_AS_Guppy_use, marginal=NULL, interval = "prediction")),
       df_AS_pred_full_Guppy <- cbind(Data_raw_AS_Guppy, predict(mcmc_AS_Guppy_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_AS_pred_full_Guppy[df_AS_pred_full_Guppy$Sex=="f",]$AS_staticT <- summary(mcmc_AS_Guppy_use)$solutions[1,1] + summary(mcmc_AS_Guppy_use)$solutions[4,1]*df_AS_pred_full_Guppy[df_AS_pred_full_Guppy$Sex=="f",]$logMass
df_AS_pred_full_Guppy[df_AS_pred_full_Guppy$Sex=="m",]$AS_staticT <- summary(mcmc_AS_Guppy_use)$solutions[1,1]+summary(mcmc_AS_Guppy_use)$solutions[3,1] + (summary(mcmc_AS_Guppy_use)$solutions[4,1]+summary(mcmc_AS_Guppy_use)$solutions[6,1])*df_AS_pred_full_Guppy[df_AS_pred_full_Guppy$Sex=="m",]$logMass

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AS_pred_full_Guppy$AS_staticM <- summary(mcmc_AS_Guppy_useM)$solutions[1,1] + summary(mcmc_AS_Guppy_useM)$solutions[3,1]*df_AS_pred_full_Guppy$logMass


## Chromis ##
#With treatments
mcmc_AS_RR_Chromis  <- MCMCglmm(logAS ~ logMassCent * Treatment + logMassMean * Treatment ,
                                 random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                 prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Chromis))

# mcmc_AS_RI_Chromis  <- MCMCglmm(logAS ~ logMassCent * Treatment + logMassMean * Treatment,
#                                  random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                  prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Chromis))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AS_Chromis_use <- mcmc_AS_RR_Chromis,ifelse(mcmc_AS_RR_Chromis$DIC<mcmc_AS_RI_Chromis$DIC,mcmc_AS_Chromis_use <- mcmc_AS_RR_Chromis,mcmc_AS_Chromis_use <- mcmc_AS_RI_Chromis))
# Modeluse[6,6] <-ifelse(mcmc_AS_RR_Chromis$DIC<mcmc_AS_RI_Chromis$DIC, "RR","RI")
# Modeluse[6,7] <- abs(mcmc_AS_RR_Chromis$DIC-mcmc_AS_RI_Chromis$DIC)

#Without treatments
mcmc_AS_RR_ChromisM  <- MCMCglmm(logAS ~ logMassCent  + logMassMean ,
                                  random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                  prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Chromis))

# mcmc_AS_RI_ChromisM  <- MCMCglmm(logAS ~ logMassCent  + logMassMean,
#                                   random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                   prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Chromis))

ifelse(use_RRM==TRUE,mcmc_AS_Chromis_useM <- mcmc_AS_RR_ChromisM,ifelse(mcmc_AS_RR_ChromisM$DIC<mcmc_AS_RI_ChromisM$DIC,mcmc_AS_Chromis_useM <- mcmc_AS_RR_ChromisM,mcmc_AS_Chromis_useM <- mcmc_AS_RI_ChromisM))
# ModeluseM[6,6] <-ifelse(mcmc_AS_RR_ChromisM$DIC<mcmc_AS_RI_ChromisM$DIC, "RR","RI")
# ModeluseM[6,7] <- abs(mcmc_AS_RR_ChromisM$DIC-mcmc_AS_RI_ChromisM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_AS_pred_full_Chromis <- cbind(Data_raw_AS_Chromis, predict(mcmc_AS_Chromis_use, marginal=NULL, interval = "prediction")),
       df_AS_pred_full_Chromis <- cbind(Data_raw_AS_Chromis, predict(mcmc_AS_Chromis_useM, marginal=NULL, interval = "prediction"))) 



# finding the statick scaling value for a given individual at a given point for the given treatment
df_AS_pred_full_Chromis[df_AS_pred_full_Chromis$Treatment=="high",]$AS_staticT <- summary(mcmc_AS_Chromis_use)$solutions[1,1] + summary(mcmc_AS_Chromis_use)$solutions[4,1]*df_AS_pred_full_Chromis[df_AS_pred_full_Chromis$Treatment=="high",]$logMass
df_AS_pred_full_Chromis[df_AS_pred_full_Chromis$Treatment=="low",]$AS_staticT <- summary(mcmc_AS_Chromis_use)$solutions[1,1]+summary(mcmc_AS_Chromis_use)$solutions[3,1] + (summary(mcmc_AS_Chromis_use)$solutions[4,1]+summary(mcmc_AS_Chromis_use)$solutions[6,1])*df_AS_pred_full_Chromis[df_AS_pred_full_Chromis$Treatment=="low",]$logMass
#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AS_pred_full_Chromis$AS_staticM <- summary(mcmc_AS_Chromis_useM)$solutions[1,1] + summary(mcmc_AS_Chromis_useM)$solutions[3,1]*df_AS_pred_full_Chromis$logMass



## Clownfish ##
mcmc_AS_RR_Clownfish  <- MCMCglmm(logAS ~ logMassCent  + logMassMean ,
                                   random = ~ us(1 + logMassCent):FishID, rcov = ~ units, family = "gaussian",
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Clownfish))

# mcmc_AS_RI_Clownfish  <- MCMCglmm(logAS ~ logMassCent + logMassMean ,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AS_Clownfish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AS_Clownfish_use <- mcmc_AS_RR_Clownfish,ifelse(mcmc_AS_RR_Clownfish$DIC<mcmc_AS_RI_Clownfish$DIC,mcmc_AS_Clownfish_use <- mcmc_AS_RR_Clownfish,mcmc_AS_Clownfish_use <- mcmc_AS_RI_Clownfish))
# Modeluse[7,6] <-ifelse(mcmc_AS_RR_Clownfish$DIC<mcmc_AS_RI_Clownfish$DIC, "RR","RI")
# Modeluse[7,7] <- abs(mcmc_AS_RR_Clownfish$DIC-mcmc_AS_RI_Clownfish$DIC)

df_AS_pred_full_Clownfish <- cbind(Data_raw_AS_Clownfish, predict(mcmc_AS_Clownfish_use, marginal=NULL, interval = "prediction")) #finding predicted varible, using the best model
# finding the statick scaling value for a given individual at a given point for the given treatment
df_AS_pred_full_Clownfish$AS_staticT <- summary(mcmc_AS_Clownfish_use)$solutions[1,1] + summary(mcmc_AS_Clownfish_use)$solutions[3,1]*df_AS_pred_full_Clownfish$logMass


#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AS_pred_full_Clownfish$AS_staticM <- summary(mcmc_AS_Clownfish_use)$solutions[1,1] + summary(mcmc_AS_Clownfish_use)$solutions[3,1]*df_AS_pred_full_Clownfish$logMass


# Recombining #
Data_AS <- bind_rows(df_AS_pred_full_Rainbow_trout, df_AS_pred_full_Brown_trout, df_AS_pred_full_Cunner, df_AS_pred_full_Guppy, df_AS_pred_full_Chromis, df_AS_pred_full_Clownfish)

names(Data_AS)[names(Data_AS) == "fit" ] <- "AS_fit"
names(Data_AS)[names(Data_AS) == "lwr" ] <- "AS_fit_low"
names(Data_AS)[names(Data_AS) == "upr" ] <- "AS_fit_high"


### AGR ###
#AGR dataframe
Data_raw_AGR <- Data_raw[Data_raw$AGR>0,]
Data_raw_AGR$AGR_staticM <- NA
Data_raw_AGR$AGR_staticT <- NA

#Subsetting data
Data_raw_AGR_Zebrafish <- subset(Data_raw_AGR, Species=="Zebrafish")
Data_raw_AGR_Rainbow_trout <- subset(Data_raw_AGR, Species=="Rainbow_trout")
Data_raw_AGR_Brown_trout <- subset(Data_raw_AGR, Species=="Brown_trout")
Data_raw_AGR_Cunner <- subset(Data_raw_AGR, Species=="Cunner")
Data_raw_AGR_Guppy <- subset(Data_raw_AGR, Species=="Guppy")
Data_raw_AGR_Chromis <- subset(Data_raw_AGR, Species=="Chromis")
Data_raw_AGR_Clownfish <- subset(Data_raw_AGR, Species=="Clownfish")

## Zebrafish ##
#With treatments
mcmc_AGR_RR_Zebrafish  <- MCMCglmm(logAGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
                                   random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Zebrafish))

# mcmc_AGR_RI_Zebrafish  <- MCMCglmm(logAGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Zebrafish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AGR_Zebrafish_use <- mcmc_AGR_RR_Zebrafish,ifelse(mcmc_AGR_RR_Zebrafish$DIC<mcmc_AGR_RI_Zebrafish$DIC,mcmc_AGR_Zebrafish_use <- mcmc_AGR_RR_Zebrafish,mcmc_AGR_Zebrafish_use <- mcmc_AGR_RI_Zebrafish))

# Modeluse[1,8] <-ifelse(mcmc_AGR_RR_Zebrafish$DIC<mcmc_AGR_RI_Zebrafish$DIC, "RR","RI")
# Modeluse[1,9] <- abs(mcmc_AGR_RR_Zebrafish$DIC-mcmc_AGR_RI_Zebrafish$DIC)
#Without treatments
mcmc_AGR_RR_ZebrafishM  <- MCMCglmm(logAGR ~ logMassCentGrowth +logMassMeanGrowth,
                                    random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                    prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Zebrafish))

# mcmc_AGR_RI_ZebrafishM  <- MCMCglmm(logAGR ~ logMassCentGrowth +logMassMeanGrowth,
#                                     random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                     prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Zebrafish))
#Finding best model and noting it
ifelse(use_RRM==TRUE,mcmc_AGR_Zebrafish_useM <- mcmc_AGR_RR_ZebrafishM,ifelse(mcmc_AGR_RR_ZebrafishM$DIC<mcmc_AGR_RI_ZebrafishM$DIC,mcmc_AGR_Zebrafish_useM <- mcmc_AGR_RR_ZebrafishM,mcmc_AGR_Zebrafish_useM <- mcmc_AGR_RI_ZebrafishM))

# ModeluseM[1,8] <-ifelse(mcmc_AGR_RR_ZebrafishM$DIC<mcmc_AGR_RI_ZebrafishM$DIC, "RR","RI")
# ModeluseM[1,9] <- abs(mcmc_AGR_RR_ZebrafishM$DIC-mcmc_AGR_RI_ZebrafishM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_AGR_pred_full_Zebrafish <- cbind(Data_raw_AGR_Zebrafish, predict(mcmc_AGR_Zebrafish_use, marginal=NULL, interval = "prediction")),
       df_AGR_pred_full_Zebrafish <- cbind(Data_raw_AGR_Zebrafish, predict(mcmc_AGR_Zebrafish_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_AGR_pred_full_Zebrafish[df_AGR_pred_full_Zebrafish$Treatment=="H",]$AGR_staticT <- df_AGR_pred_full_Zebrafish[df_AGR_pred_full_Zebrafish$Treatment=="H",]$logMassGrowth*summary(mcmc_AGR_Zebrafish_use)$solutions[5,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[1,1]
df_AGR_pred_full_Zebrafish[df_AGR_pred_full_Zebrafish$Treatment=="L",]$AGR_staticT <- df_AGR_pred_full_Zebrafish[df_AGR_pred_full_Zebrafish$Treatment=="L",]$logMassGrowth*(summary(mcmc_AGR_Zebrafish_use)$solutions[5,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[8,1])+summary(mcmc_AGR_Zebrafish_use)$solutions[1,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[3,1]
df_AGR_pred_full_Zebrafish[df_AGR_pred_full_Zebrafish$Treatment=="M",]$AGR_staticT <- df_AGR_pred_full_Zebrafish[df_AGR_pred_full_Zebrafish$Treatment=="M",]$logMassGrowth*(summary(mcmc_AGR_Zebrafish_use)$solutions[5,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[9,1])+summary(mcmc_AGR_Zebrafish_use)$solutions[1,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[4,1]

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AGR_pred_full_Zebrafish$AGR_staticM <- summary(mcmc_AGR_Zebrafish_useM)$solutions[1,1] + summary(mcmc_AGR_Zebrafish_useM)$solutions[3,1]*df_AGR_pred_full_Zebrafish$logMassGrowth


## Rainbow trout ##
#with treatments
mcmc_AGR_RR_Rainbow_trout  <- MCMCglmm(logAGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
                                   random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Rainbow_trout))

# mcmc_AGR_RI_Rainbow_trout  <- MCMCglmm(logAGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Rainbow_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AGR_Rainbow_trout_use <- mcmc_AGR_RR_Rainbow_trout,ifelse(mcmc_AGR_RR_Rainbow_trout$DIC<mcmc_AGR_RI_Rainbow_trout$DIC,mcmc_AGR_Rainbow_trout_use <- mcmc_AGR_RR_Rainbow_trout,mcmc_AGR_Rainbow_trout_use <- mcmc_AGR_RI_Rainbow_trout))

# Modeluse[2,8] <-ifelse(mcmc_AGR_RR_Rainbow_trout$DIC<mcmc_AGR_RI_Rainbow_trout$DIC, "RR","RI")
# Modeluse[2,9] <- abs(mcmc_AGR_RR_Rainbow_trout$DIC-mcmc_AGR_RI_Rainbow_trout$DIC)
#Without treatments
mcmc_AGR_RR_Rainbow_troutM  <- MCMCglmm(logAGR ~ logMassCentGrowth  + logMassMeanGrowth ,
                                    random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                    prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Rainbow_trout))

# mcmc_AGR_RI_Rainbow_troutM  <- MCMCglmm(logAGR ~ logMassCentGrowth  + logMassMeanGrowth,
#                                     random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                     prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Rainbow_trout))

ifelse(use_RRM==TRUE,mcmc_AGR_Rainbow_trout_useM <- mcmc_AGR_RR_Rainbow_troutM,ifelse(mcmc_AGR_RR_Rainbow_troutM$DIC<mcmc_AGR_RI_Rainbow_troutM$DIC,mcmc_AGR_Rainbow_trout_useM <- mcmc_AGR_RR_Rainbow_troutM,mcmc_AGR_Rainbow_trout_useM <- mcmc_AGR_RI_Rainbow_troutM))

# ModeluseM[2,8] <-ifelse(mcmc_AGR_RR_Rainbow_troutM$DIC<mcmc_AGR_RI_Rainbow_troutM$DIC, "RR","RI")
# ModeluseM[2,9] <- abs(mcmc_AGR_RR_Rainbow_troutM$DIC-mcmc_AGR_RI_Rainbow_troutM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_AGR_pred_full_Rainbow_trout <- cbind(Data_raw_AGR_Rainbow_trout, predict(mcmc_AGR_Rainbow_trout_use, marginal=NULL, interval = "prediction")),
       df_AGR_pred_full_Rainbow_trout <- cbind(Data_raw_AGR_Rainbow_trout, predict(mcmc_AGR_Rainbow_trout_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_AGR_pred_full_Rainbow_trout[df_AGR_pred_full_Rainbow_trout$Treatment=="10Cegg",]$AGR_staticT <- df_AGR_pred_full_Rainbow_trout[df_AGR_pred_full_Rainbow_trout$Treatment=="10Cegg",]$logMassGrowth*summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,1]
df_AGR_pred_full_Rainbow_trout[df_AGR_pred_full_Rainbow_trout$Treatment=="14Cegg",]$AGR_staticT <- df_AGR_pred_full_Rainbow_trout[df_AGR_pred_full_Rainbow_trout$Treatment=="14Cegg",]$logMassGrowth*(summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[8,1])+summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[3,1]
df_AGR_pred_full_Rainbow_trout[df_AGR_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$AGR_staticT <- df_AGR_pred_full_Rainbow_trout[df_AGR_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$logMassGrowth*(summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[9,1])+summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[4,1]

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AGR_pred_full_Rainbow_trout$AGR_staticM <- summary(mcmc_AGR_Rainbow_trout_useM)$solutions[1,1] + summary(mcmc_AGR_Rainbow_trout_useM)$solutions[3,1]*df_AGR_pred_full_Rainbow_trout$logMassGrowth


## Brown trout ##
mcmc_AGR_RR_Brown_trout  <- MCMCglmm(logAGR ~ logMassCentGrowth + logMassMeanGrowth,
                                     random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian", 
                                     prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Brown_trout))
# mcmc_AGR_RI_Brown_trout  <- MCMCglmm(logAGR ~ logMassCentGrowth  + logMassMeanGrowth ,
#                                        random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                        prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                        verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Brown_trout))
# #Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AGR_Brown_trout_use <- mcmc_AGR_RR_Brown_trout,ifelse(mcmc_AGR_RR_Brown_trout$DIC<mcmc_AGR_RI_Brown_trout$DIC,mcmc_AGR_Brown_trout_use <- mcmc_AGR_RR_Brown_trout,mcmc_AGR_Brown_trout_use <- mcmc_AGR_RI_Brown_trout))

# Modeluse[3,8] <-ifelse(mcmc_AGR_RR_Brown_trout$DIC<mcmc_AGR_RI_Brown_trout$DIC, "RR","RI")
# Modeluse[3,9] <- abs(mcmc_AGR_RR_Brown_trout$DIC-mcmc_AGR_RI_Brown_trout$DIC)

df_AGR_pred_full_Brown_trout <- cbind(Data_raw_AGR_Brown_trout, predict(mcmc_AGR_Brown_trout_use, marginal=NULL, interval = "prediction"))#finding predicted varible, using the best model

# finding the statick scaling value for a given individual at a given point for the given treatment
df_AGR_pred_full_Brown_trout$AGR_staticT <- summary(mcmc_AGR_Brown_trout_use)$solutions[1,1] + summary(mcmc_AGR_Brown_trout_use)$solutions[3,1]*df_AGR_pred_full_Brown_trout$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AGR_pred_full_Brown_trout$AGR_staticM <- summary(mcmc_AGR_Brown_trout_use)$solutions[1,1] + summary(mcmc_AGR_Brown_trout_use)$solutions[3,1]*df_AGR_pred_full_Brown_trout$logMassGrowth


## Cunner ##
mcmc_AGR_RR_Cunner  <- MCMCglmm(logAGR ~ logMassCentGrowth + logMassMeanGrowth,
                                     random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian", 
                                     prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Cunner))
# mcmc_AGR_RI_Cunner  <- MCMCglmm(logAGR ~ logMassCentGrowth  + logMassMeanGrowth ,
#                                      random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                      prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                      verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Cunner))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AGR_Cunner_use <- mcmc_AGR_RR_Cunner,ifelse(mcmc_AGR_RR_Cunner$DIC<mcmc_AGR_RI_Cunner$DIC,mcmc_AGR_Cunner_use <- mcmc_AGR_RR_Cunner,mcmc_AGR_Cunner_use <- mcmc_AGR_RI_Cunner))

# Modeluse[4,8] <-ifelse(mcmc_AGR_RR_Cunner$DIC<mcmc_AGR_RI_Cunner$DIC, "RR","RI")
# Modeluse[4,9] <- abs(mcmc_AGR_RR_Cunner$DIC-mcmc_AGR_RI_Cunner$DIC)

df_AGR_pred_full_Cunner <- cbind(Data_raw_AGR_Cunner, predict(mcmc_AGR_Cunner_use, marginal=NULL, interval = "prediction"))#finding predicted varible, using the best model

# finding the statick scaling value for a given individual at a given point for the given treatment
df_AGR_pred_full_Cunner$AGR_staticT <- summary(mcmc_AGR_Cunner_use)$solutions[1,1] + summary(mcmc_AGR_Cunner_use)$solutions[3,1]*df_AGR_pred_full_Cunner$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AGR_pred_full_Cunner$AGR_staticM <- summary(mcmc_AGR_Cunner_use)$solutions[1,1] + summary(mcmc_AGR_Cunner_use)$solutions[3,1]*df_AGR_pred_full_Cunner$logMassGrowth

## Guppy ##

#with sex
mcmc_AGR_RR_Guppy  <- MCMCglmm(logAGR ~ logMassCentGrowth * Sex + logMassMeanGrowth * Sex,
                               random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                               prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Guppy))
# mcmc_AGR_RI_Guppy  <- MCMCglmm(logAGR ~ logMassCentGrowth * Sex + logMassMeanGrowth * Sex,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Guppy))
# #Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AGR_Guppy_use <- mcmc_AGR_RR_Guppy,ifelse(mcmc_AGR_RR_Guppy$DIC<mcmc_AGR_RI_Guppy$DIC,mcmc_AGR_Guppy_use <- mcmc_AGR_RR_Guppy,mcmc_AGR_Guppy_use <- mcmc_AGR_RI_Guppy))

# Modeluse[5,8] <-ifelse(mcmc_AGR_RR_Guppy$DIC<mcmc_AGR_RI_Guppy$DIC, "RR","RI")
# Modeluse[5,9] <- abs(mcmc_AGR_RR_Guppy$DIC-mcmc_AGR_RI_Guppy$DIC)
# 
#without sex
mcmc_AGR_RR_GuppyM  <- MCMCglmm(logAGR ~ logMassCentGrowth + logMassMeanGrowth,
                               random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                               prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Guppy))

# mcmc_AGR_RI_GuppyM  <- MCMCglmm(logAGR ~ logMassCentGrowth + logMassMeanGrowth,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Guppy))
# #Finding best model and noting it
ifelse(use_RRM==TRUE,mcmc_AGR_Guppy_useM <- mcmc_AGR_RR_GuppyM,ifelse(mcmc_AGR_RR_GuppyM$DIC<mcmc_AGR_RI_GuppyM$DIC,mcmc_AGR_Guppy_useM <- mcmc_AGR_RR_GuppyM,mcmc_AGR_Guppy_useM <- mcmc_AGR_RI_GuppyM))
# 
# ModeluseM[5,8] <-ifelse(mcmc_AGR_RR_GuppyM$DIC<mcmc_AGR_RI_GuppyM$DIC, "RR","RI")
# ModeluseM[5,9] <- abs(mcmc_AGR_RR_GuppyM$DIC-mcmc_AGR_RI_GuppyM$DIC)



ifelse(use_treatments==TRUE,df_AGR_pred_full_Guppy <- cbind(Data_raw_AGR_Guppy, predict(mcmc_AGR_Guppy_use, marginal=NULL, interval = "prediction")),
       df_AGR_pred_full_Guppy <- cbind(Data_raw_AGR_Guppy, predict(mcmc_AGR_Guppy_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_AGR_pred_full_Guppy[df_AGR_pred_full_Guppy$Sex=="f",]$AGR_staticT <- summary(mcmc_AGR_Guppy_use)$solutions[1,1] + summary(mcmc_AGR_Guppy_use)$solutions[4,1]*df_AGR_pred_full_Guppy[df_AGR_pred_full_Guppy$Sex=="f",]$logMassGrowth
df_AGR_pred_full_Guppy[df_AGR_pred_full_Guppy$Sex=="m",]$AGR_staticT <- summary(mcmc_AGR_Guppy_use)$solutions[1,1]+summary(mcmc_AGR_Guppy_use)$solutions[3,1] + (summary(mcmc_AGR_Guppy_use)$solutions[4,1]+summary(mcmc_AGR_Guppy_use)$solutions[6,1])*df_AGR_pred_full_Guppy[df_AGR_pred_full_Guppy$Sex=="m",]$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AGR_pred_full_Guppy$AGR_staticM <- summary(mcmc_AGR_Guppy_useM)$solutions[1,1] + summary(mcmc_AGR_Guppy_useM)$solutions[3,1]*df_AGR_pred_full_Guppy$logMassGrowth


## Chromis ##
#with treatment
mcmc_AGR_RR_Chromis  <- MCMCglmm(logAGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
                               random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                               prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Chromis))

# mcmc_AGR_RI_Chromis  <- MCMCglmm(logAGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Chromis))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AGR_Chromis_use <- mcmc_AGR_RR_Chromis,ifelse(mcmc_AGR_RR_Chromis$DIC<mcmc_AGR_RI_Chromis$DIC,mcmc_AGR_Chromis_use <- mcmc_AGR_RR_Chromis,mcmc_AGR_Chromis_use <- mcmc_AGR_RI_Chromis))
# 
# Modeluse[6,8] <-ifelse(mcmc_AGR_RR_Chromis$DIC<mcmc_AGR_RI_Chromis$DIC, "RR","RI")
# Modeluse[6,9] <- abs(mcmc_AGR_RR_Chromis$DIC-mcmc_AGR_RI_Chromis$DIC)

#without treatment
mcmc_AGR_RR_ChromisM  <- MCMCglmm(logAGR ~ logMassCentGrowth + logMassMeanGrowth,
                                random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Chromis))

# mcmc_AGR_RI_ChromisM  <- MCMCglmm(logAGR ~ logMassCentGrowth + logMassMeanGrowth,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Chromis))
# #Finding best model and noting it
ifelse(use_RRM==TRUE,mcmc_AGR_Chromis_useM <- mcmc_AGR_RR_ChromisM,ifelse(mcmc_AGR_RR_ChromisM$DIC<mcmc_AGR_RI_ChromisM$DIC,mcmc_AGR_Chromis_useM <- mcmc_AGR_RR_ChromisM,mcmc_AGR_Chromis_useM <- mcmc_AGR_RI_ChromisM))

# ModeluseM[6,8] <-ifelse(mcmc_AGR_RR_ChromisM$DIC<mcmc_AGR_RI_ChromisM$DIC, "RR","RI")
# ModeluseM[6,9] <- abs(mcmc_AGR_RR_ChromisM$DIC-mcmc_AGR_RI_ChromisM$DIC)



ifelse(use_treatments==TRUE,df_AGR_pred_full_Chromis <- cbind(Data_raw_AGR_Chromis, predict(mcmc_AGR_Chromis_use, marginal=NULL, interval = "prediction")),
       df_AGR_pred_full_Chromis <- cbind(Data_raw_AGR_Chromis, predict(mcmc_AGR_Chromis_useM, marginal=NULL, interval = "prediction"))) 


# finding the statick scaling value for a given individual at a given point for the given treatment
df_AGR_pred_full_Chromis[df_AGR_pred_full_Chromis$Treatment=="high",]$AGR_staticT <- summary(mcmc_AGR_Chromis_use)$solutions[1,1] + summary(mcmc_AGR_Chromis_use)$solutions[4,1]*df_AGR_pred_full_Chromis[df_AGR_pred_full_Chromis$Treatment=="high",]$logMassGrowth
df_AGR_pred_full_Chromis[df_AGR_pred_full_Chromis$Treatment=="low",]$AGR_staticT <- summary(mcmc_AGR_Chromis_use)$solutions[1,1]+summary(mcmc_AGR_Chromis_use)$solutions[3,1] + (summary(mcmc_AGR_Chromis_use)$solutions[4,1]+summary(mcmc_AGR_Chromis_use)$solutions[6,1])*df_AGR_pred_full_Chromis[df_AGR_pred_full_Chromis$Treatment=="low",]$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AGR_pred_full_Chromis$AGR_staticM <- summary(mcmc_AGR_Chromis_useM)$solutions[1,1] + summary(mcmc_AGR_Chromis_useM)$solutions[3,1]*df_AGR_pred_full_Chromis$logMassGrowth


## Clownfish ##
mcmc_AGR_RR_Clownfish  <- MCMCglmm(logAGR ~ logMassCentGrowth + logMassMeanGrowth,
                                     random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian", 
                                     prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Clownfish))
# mcmc_AGR_RI_Clownfish  <- MCMCglmm(logAGR ~ logMassCentGrowth  + logMassMeanGrowth ,
#                                      random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                      prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                      verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_AGR_Clownfish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_AGR_Clownfish_use <- mcmc_AGR_RR_Clownfish,ifelse(mcmc_AGR_RR_Clownfish$DIC<mcmc_AGR_RI_Clownfish$DIC,mcmc_AGR_Clownfish_use <- mcmc_AGR_RR_Clownfish,mcmc_AGR_Clownfish_use <- mcmc_AGR_RI_Clownfish))

# Modeluse[7,8] <-ifelse(mcmc_AGR_Clownfish_use$DIC<mcmc_AGR_RI_Clownfish$DIC, "RR","RI")
# Modeluse[7,9] <- abs(mcmc_AGR_Clownfish_use$DIC-mcmc_AGR_RI_Clownfish$DIC)

df_AGR_pred_full_Clownfish <- cbind(Data_raw_AGR_Clownfish, predict(mcmc_AGR_Clownfish_use, marginal=NULL, interval = "prediction"))#finding predicted varible, using the best model

# finding the statick scaling value for a given individual at a given point for the given treatment
df_AGR_pred_full_Clownfish$AGR_staticT <- summary(mcmc_AGR_Clownfish_use)$solutions[1,1] + summary(mcmc_AGR_Clownfish_use)$solutions[3,1]*df_AGR_pred_full_Clownfish$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_AGR_pred_full_Clownfish$AGR_staticM <- summary(mcmc_AGR_Clownfish_use)$solutions[1,1] + summary(mcmc_AGR_Clownfish_use)$solutions[3,1]*df_AGR_pred_full_Clownfish$logMassGrowth


# Recombining #
Data_AGR <- bind_rows(df_AGR_pred_full_Zebrafish, df_AGR_pred_full_Rainbow_trout, df_AGR_pred_full_Brown_trout, df_AGR_pred_full_Cunner, df_AGR_pred_full_Guppy, df_AGR_pred_full_Chromis, df_AGR_pred_full_Clownfish)

names(Data_AGR)[names(Data_AGR) == "fit" ] <- "AGR_fit"
names(Data_AGR)[names(Data_AGR) == "lwr" ] <- "AGR_fit_low"
names(Data_AGR)[names(Data_AGR) == "upr" ] <- "AGR_fit_high"


### SGR ###
#SGR dataframe
Data_raw_SGR <- Data_raw[Data_raw$SGR>0,]
Data_raw_SGR$SGR_staticM <- NA
Data_raw_SGR$SGR_staticT <- NA

#Subsetting data
Data_raw_SGR_Zebrafish <- subset(Data_raw_SGR, Species=="Zebrafish")
Data_raw_SGR_Rainbow_trout <- subset(Data_raw_SGR, Species=="Rainbow_trout")
Data_raw_SGR_Brown_trout <- subset(Data_raw_SGR, Species=="Brown_trout")
Data_raw_SGR_Cunner <- subset(Data_raw_SGR, Species=="Cunner")
Data_raw_SGR_Guppy <- subset(Data_raw_SGR, Species=="Guppy",Sex != "n" & Sex != "d")
Data_raw_SGR_Chromis <- subset(Data_raw_SGR, Species=="Chromis")
Data_raw_SGR_Clownfish <- subset(Data_raw_SGR, Species=="Clownfish")

## Zebrafish ##
#With treatments
mcmc_SGR_RR_Zebrafish  <- MCMCglmm(logSGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
                                   random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Zebrafish))

# mcmc_SGR_RI_Zebrafish  <- MCMCglmm(logSGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Zebrafish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SGR_Zebrafish_use <- mcmc_SGR_RR_Zebrafish,ifelse(mcmc_SGR_RR_Zebrafish$DIC<mcmc_SGR_RI_Zebrafish$DIC,mcmc_SGR_Zebrafish_use <- mcmc_SGR_RR_Zebrafish,mcmc_SGR_Zebrafish_use <- mcmc_SGR_RI_Zebrafish))

# Modeluse[1,10] <-ifelse(mcmc_SGR_RR_Zebrafish$DIC<mcmc_SGR_RI_Zebrafish$DIC, "RR","RI")
# Modeluse[1,11] <- abs(mcmc_SGR_RR_Zebrafish$DIC-mcmc_SGR_RI_Zebrafish$DIC)
#Without treatments
mcmc_SGR_RR_ZebrafishM  <- MCMCglmm(logSGR ~ logMassCentGrowth +logMassMeanGrowth,
                                    random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                    prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Zebrafish))

# mcmc_SGR_RI_ZebrafishM  <- MCMCglmm(logSGR ~ logMassCentGrowth +logMassMeanGrowth,
#                                     random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                     prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Zebrafish))
#Finding best model and noting it
ifelse(use_RRM==TRUE,mcmc_SGR_Zebrafish_useM <- mcmc_SGR_RR_ZebrafishM,ifelse(mcmc_SGR_RR_ZebrafishM$DIC<mcmc_SGR_RI_ZebrafishM$DIC,mcmc_SGR_Zebrafish_useM <- mcmc_SGR_RR_ZebrafishM,mcmc_SGR_Zebrafish_useM <- mcmc_SGR_RI_ZebrafishM))

# ModeluseM[1,10] <-ifelse(mcmc_SGR_RR_ZebrafishM$DIC<mcmc_SGR_RI_ZebrafishM$DIC, "RR","RI")
# ModeluseM[1,11] <- abs(mcmc_SGR_RR_ZebrafishM$DIC-mcmc_SGR_RI_ZebrafishM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_SGR_pred_full_Zebrafish <- cbind(Data_raw_SGR_Zebrafish, predict(mcmc_SGR_Zebrafish_use, marginal=NULL, interval = "prediction")),
       df_SGR_pred_full_Zebrafish <- cbind(Data_raw_SGR_Zebrafish, predict(mcmc_SGR_Zebrafish_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_SGR_pred_full_Zebrafish[df_SGR_pred_full_Zebrafish$Treatment=="H",]$SGR_staticT <- df_SGR_pred_full_Zebrafish[df_SGR_pred_full_Zebrafish$Treatment=="H",]$logMassGrowth*summary(mcmc_SGR_Zebrafish_use)$solutions[5,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[1,1]
df_SGR_pred_full_Zebrafish[df_SGR_pred_full_Zebrafish$Treatment=="L",]$SGR_staticT <- df_SGR_pred_full_Zebrafish[df_SGR_pred_full_Zebrafish$Treatment=="L",]$logMassGrowth*(summary(mcmc_SGR_Zebrafish_use)$solutions[5,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[8,1])+summary(mcmc_SGR_Zebrafish_use)$solutions[1,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[3,1]
df_SGR_pred_full_Zebrafish[df_SGR_pred_full_Zebrafish$Treatment=="M",]$SGR_staticT <- df_SGR_pred_full_Zebrafish[df_SGR_pred_full_Zebrafish$Treatment=="M",]$logMassGrowth*(summary(mcmc_SGR_Zebrafish_use)$solutions[5,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[9,1])+summary(mcmc_SGR_Zebrafish_use)$solutions[1,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[4,1]

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SGR_pred_full_Zebrafish$SGR_staticM <- summary(mcmc_SGR_Zebrafish_useM)$solutions[1,1] + summary(mcmc_SGR_Zebrafish_useM)$solutions[3,1]*df_SGR_pred_full_Zebrafish$logMassGrowth


## Rainbow trout ##
#with treatments
mcmc_SGR_RR_Rainbow_trout  <- MCMCglmm(logSGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
                                       random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                       prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                       verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Rainbow_trout))

# mcmc_SGR_RI_Rainbow_trout  <- MCMCglmm(logSGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
#                                        random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                        prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                        verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Rainbow_trout))
# #Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SGR_Rainbow_trout_use <- mcmc_SGR_RR_Rainbow_trout,ifelse(mcmc_SGR_RR_Rainbow_trout$DIC<mcmc_SGR_RI_Rainbow_trout$DIC,mcmc_SGR_Rainbow_trout_use <- mcmc_SGR_RR_Rainbow_trout,mcmc_SGR_Rainbow_trout_use <- mcmc_SGR_RI_Rainbow_trout))

# Modeluse[2,10] <-ifelse(mcmc_SGR_RR_Rainbow_trout$DIC<mcmc_SGR_RI_Rainbow_trout$DIC, "RR","RI")
# Modeluse[2,11] <- abs(mcmc_SGR_RR_Rainbow_trout$DIC-mcmc_SGR_RI_Rainbow_trout$DIC)
#Without treatments
mcmc_SGR_RR_Rainbow_troutM  <- MCMCglmm(logSGR ~ logMassCentGrowth  + logMassMeanGrowth  + Temp,
                                        random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                        prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                        verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Rainbow_trout))

# mcmc_SGR_RI_Rainbow_troutM  <- MCMCglmm(logSGR ~ logMassCentGrowth  + logMassMeanGrowth  + Temp,
#                                         random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                         prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                         verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Rainbow_trout))

ifelse(use_RRM==TRUE,mcmc_SGR_Rainbow_trout_useM <- mcmc_SGR_RR_Rainbow_troutM,ifelse(mcmc_SGR_RR_Rainbow_troutM$DIC<mcmc_SGR_RI_Rainbow_troutM$DIC,mcmc_SGR_Rainbow_trout_useM <- mcmc_SGR_RR_Rainbow_troutM,mcmc_SGR_Rainbow_trout_useM <- mcmc_SGR_RI_Rainbow_troutM))

# ModeluseM[2,10] <-ifelse(mcmc_SGR_RR_Rainbow_troutM$DIC<mcmc_SGR_RI_Rainbow_troutM$DIC, "RR","RI")
# ModeluseM[2,11] <- abs(mcmc_SGR_RR_Rainbow_troutM$DIC-mcmc_SGR_RI_Rainbow_troutM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_SGR_pred_full_Rainbow_trout <- cbind(Data_raw_SGR_Rainbow_trout, predict(mcmc_SGR_Rainbow_trout_use, marginal=NULL, interval = "prediction")),
       df_SGR_pred_full_Rainbow_trout <- cbind(Data_raw_SGR_Rainbow_trout, predict(mcmc_SGR_Rainbow_trout_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_SGR_pred_full_Rainbow_trout[df_SGR_pred_full_Rainbow_trout$Treatment=="10Cegg",]$SGR_staticT <- df_SGR_pred_full_Rainbow_trout[df_SGR_pred_full_Rainbow_trout$Treatment=="10Cegg",]$logMassGrowth*summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,1]
df_SGR_pred_full_Rainbow_trout[df_SGR_pred_full_Rainbow_trout$Treatment=="14Cegg",]$SGR_staticT <- df_SGR_pred_full_Rainbow_trout[df_SGR_pred_full_Rainbow_trout$Treatment=="14Cegg",]$logMassGrowth*(summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[8,1])+summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[3,1]
df_SGR_pred_full_Rainbow_trout[df_SGR_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$SGR_staticT <- df_SGR_pred_full_Rainbow_trout[df_SGR_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$logMassGrowth*(summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[9,1])+summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[4,1]

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SGR_pred_full_Rainbow_trout$SGR_staticM <- summary(mcmc_SGR_Rainbow_trout_useM)$solutions[1,1] + summary(mcmc_SGR_Rainbow_trout_useM)$solutions[3,1]*df_SGR_pred_full_Rainbow_trout$logMassGrowth


## Brown trout ##
mcmc_SGR_RR_Brown_trout  <- MCMCglmm(logSGR ~ logMassCentGrowth + logMassMeanGrowth,
                                     random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian", 
                                     prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Brown_trout))
# mcmc_SGR_RI_Brown_trout  <- MCMCglmm(logSGR ~ logMassCentGrowth  + logMassMeanGrowth ,
#                                      random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                      prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                      verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Brown_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SGR_Brown_trout_use <- mcmc_SGR_RR_Brown_trout,ifelse(mcmc_SGR_RR_Brown_trout$DIC<mcmc_SGR_RI_Brown_trout$DIC,mcmc_SGR_Brown_trout_use <- mcmc_SGR_RR_Brown_trout,mcmc_SGR_Brown_trout_use <- mcmc_SGR_RI_Brown_trout))

# Modeluse[3,10] <-ifelse(mcmc_SGR_RR_Brown_trout$DIC<mcmc_SGR_RI_Brown_trout$DIC, "RR","RI")
# Modeluse[3,11] <- abs(mcmc_SGR_RR_Brown_trout$DIC-mcmc_SGR_RI_Brown_trout$DIC)

df_SGR_pred_full_Brown_trout <- cbind(Data_raw_SGR_Brown_trout, predict(mcmc_SGR_Brown_trout_use, marginal=NULL, interval = "prediction"))#finding predicted varible, using the best model

# finding the statick scaling value for a given individual at a given point for the given treatment
df_SGR_pred_full_Brown_trout$SGR_staticT <- summary(mcmc_SGR_Brown_trout_use)$solutions[1,1] + summary(mcmc_SGR_Brown_trout_use)$solutions[3,1]*df_SGR_pred_full_Brown_trout$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SGR_pred_full_Brown_trout$SGR_staticM <- summary(mcmc_SGR_Brown_trout_use)$solutions[1,1] + summary(mcmc_SGR_Brown_trout_use)$solutions[3,1]*df_SGR_pred_full_Brown_trout$logMassGrowth


## Cunner ##
mcmc_SGR_RR_Cunner  <- MCMCglmm(logSGR ~ logMassCentGrowth + logMassMeanGrowth,
                                random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian", 
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Cunner))
# mcmc_SGR_RI_Cunner  <- MCMCglmm(logSGR ~ logMassCentGrowth  + logMassMeanGrowth ,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Cunner))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SGR_Cunner_use <- mcmc_SGR_RR_Cunner,ifelse(mcmc_SGR_RR_Cunner$DIC<mcmc_SGR_RI_Cunner$DIC,mcmc_SGR_Cunner_use <- mcmc_SGR_RR_Cunner,mcmc_SGR_Cunner_use <- mcmc_SGR_RI_Cunner))

# Modeluse[4,10] <-ifelse(mcmc_SGR_RR_Cunner$DIC<mcmc_SGR_RI_Cunner$DIC, "RR","RI")
# Modeluse[4,11] <- abs(mcmc_SGR_RR_Cunner$DIC-mcmc_SGR_RI_Cunner$DIC)

df_SGR_pred_full_Cunner <- cbind(Data_raw_SGR_Cunner, predict(mcmc_SGR_Cunner_use, marginal=NULL, interval = "prediction"))#finding predicted varible, using the best model

# finding the statick scaling value for a given individual at a given point for the given treatment
df_SGR_pred_full_Cunner$SGR_staticT <- summary(mcmc_SGR_Cunner_use)$solutions[1,1] + summary(mcmc_SGR_Cunner_use)$solutions[3,1]*df_SGR_pred_full_Cunner$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SGR_pred_full_Cunner$SGR_staticM <- summary(mcmc_SGR_Cunner_use)$solutions[1,1] + summary(mcmc_SGR_Cunner_use)$solutions[3,1]*df_SGR_pred_full_Cunner$logMassGrowth

## Guppy ##
#with sex
mcmc_SGR_RR_Guppy  <- MCMCglmm(logSGR ~ logMassCentGrowth * Sex + logMassMeanGrowth * Sex,
                               random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                               prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Guppy))

# mcmc_SGR_RI_Guppy  <- MCMCglmm(logSGR ~ logMassCentGrowth * Sex + logMassMeanGrowth * Sex,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Guppy))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SGR_Guppy_use <- mcmc_SGR_RR_Guppy,ifelse(mcmc_SGR_RR_Guppy$DIC<mcmc_SGR_RI_Guppy$DIC,mcmc_SGR_Guppy_use <- mcmc_SGR_RR_Guppy,mcmc_SGR_Guppy_use <- mcmc_SGR_RI_Guppy))

# Modeluse[5,10] <-ifelse(mcmc_SGR_RR_Guppy$DIC<mcmc_SGR_RI_Guppy$DIC, "RR","RI")
# Modeluse[5,11] <- abs(mcmc_SGR_RR_Guppy$DIC-mcmc_SGR_RI_Guppy$DIC)


#without sex
mcmc_SGR_RR_GuppyM  <- MCMCglmm(logSGR ~ logMassCentGrowth + logMassMeanGrowth,
                                random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Guppy))

# mcmc_SGR_RI_GuppyM  <- MCMCglmm(logSGR ~ logMassCentGrowth + logMassMeanGrowth,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Guppy))
#Finding best model and noting it
ifelse(use_RRM==TRUE,mcmc_SGR_Guppy_useM <- mcmc_SGR_RR_GuppyM,ifelse(mcmc_SGR_RR_GuppyM$DIC<mcmc_SGR_RI_GuppyM$DIC,mcmc_SGR_Guppy_useM <- mcmc_SGR_RR_GuppyM,mcmc_SGR_Guppy_useM <- mcmc_SGR_RI_GuppyM))

# ModeluseM[5,10] <-ifelse(mcmc_SGR_RR_GuppyM$DIC<mcmc_SGR_RI_GuppyM$DIC, "RR","RI")
# ModeluseM[5,11] <- abs(mcmc_SGR_RR_GuppyM$DIC-mcmc_SGR_RI_GuppyM$DIC)



ifelse(use_treatments==TRUE,df_SGR_pred_full_Guppy <- cbind(Data_raw_SGR_Guppy, predict(mcmc_SGR_Guppy_use, marginal=NULL, interval = "prediction")),
       df_SGR_pred_full_Guppy <- cbind(Data_raw_SGR_Guppy, predict(mcmc_SGR_Guppy_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_SGR_pred_full_Guppy[df_SGR_pred_full_Guppy$Sex=="f",]$SGR_staticT <- summary(mcmc_SGR_Guppy_use)$solutions[1,1] + summary(mcmc_SGR_Guppy_use)$solutions[4,1]*df_SGR_pred_full_Guppy[df_SGR_pred_full_Guppy$Sex=="f",]$logMassGrowth
df_SGR_pred_full_Guppy[df_SGR_pred_full_Guppy$Sex=="m",]$SGR_staticT <- summary(mcmc_SGR_Guppy_use)$solutions[1,1]+summary(mcmc_SGR_Guppy_use)$solutions[3,1] + (summary(mcmc_SGR_Guppy_use)$solutions[4,1]+summary(mcmc_SGR_Guppy_use)$solutions[6,1])*df_SGR_pred_full_Guppy[df_SGR_pred_full_Guppy$Sex=="m",]$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SGR_pred_full_Guppy$SGR_staticM <- summary(mcmc_SGR_Guppy_useM)$solutions[1,1] + summary(mcmc_SGR_Guppy_useM)$solutions[3,1]*df_SGR_pred_full_Guppy$logMassGrowth


## Chromis ##
#with treatment
mcmc_SGR_RR_Chromis  <- MCMCglmm(logSGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
                                 random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                 prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Chromis))

# mcmc_SGR_RI_Chromis  <- MCMCglmm(logSGR ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
#                                  random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                  prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Chromis))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SGR_Chromis_use <- mcmc_SGR_RR_Chromis,ifelse(mcmc_SGR_RR_Chromis$DIC<mcmc_SGR_RI_Chromis$DIC,mcmc_SGR_Chromis_use <- mcmc_SGR_RR_Chromis,mcmc_SGR_Chromis_use <- mcmc_SGR_RI_Chromis))

# Modeluse[6,10] <-ifelse(mcmc_SGR_RR_Chromis$DIC<mcmc_SGR_RI_Chromis$DIC, "RR","RI")
# Modeluse[6,11] <- abs(mcmc_SGR_RR_Chromis$DIC-mcmc_SGR_RI_Chromis$DIC)

#without treatment
mcmc_SGR_RR_ChromisM  <- MCMCglmm(logSGR ~ logMassCentGrowth + logMassMeanGrowth,
                                  random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                  prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Chromis))

# mcmc_SGR_RI_ChromisM  <- MCMCglmm(logSGR ~ logMassCentGrowth + logMassMeanGrowth,
#                                   random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                   prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Chromis))
#Finding best model and noting it
ifelse(use_RRM==TRUE,mcmc_SGR_Chromis_useM <- mcmc_SGR_RR_ChromisM,ifelse(mcmc_SGR_RR_ChromisM$DIC<mcmc_SGR_RI_ChromisM$DIC,mcmc_SGR_Chromis_useM <- mcmc_SGR_RR_ChromisM,mcmc_SGR_Chromis_useM <- mcmc_SGR_RI_ChromisM))

# ModeluseM[6,10] <-ifelse(mcmc_SGR_RR_ChromisM$DIC<mcmc_SGR_RI_ChromisM$DIC, "RR","RI")
# ModeluseM[6,11] <- abs(mcmc_SGR_RR_ChromisM$DIC-mcmc_SGR_RI_ChromisM$DIC)



ifelse(use_treatments==TRUE,df_SGR_pred_full_Chromis <- cbind(Data_raw_SGR_Chromis, predict(mcmc_SGR_Chromis_use, marginal=NULL, interval = "prediction")),
       df_SGR_pred_full_Chromis <- cbind(Data_raw_SGR_Chromis, predict(mcmc_SGR_Chromis_useM, marginal=NULL, interval = "prediction"))) 


# finding the statick scaling value for a given individual at a given point for the given treatment
df_SGR_pred_full_Chromis[df_SGR_pred_full_Chromis$Treatment=="high",]$SGR_staticT <- summary(mcmc_SGR_Chromis_use)$solutions[1,1] + summary(mcmc_SGR_Chromis_use)$solutions[4,1]*df_SGR_pred_full_Chromis[df_SGR_pred_full_Chromis$Treatment=="high",]$logMassGrowth
df_SGR_pred_full_Chromis[df_SGR_pred_full_Chromis$Treatment=="low",]$SGR_staticT <- summary(mcmc_SGR_Chromis_use)$solutions[1,1]+summary(mcmc_SGR_Chromis_use)$solutions[3,1] + (summary(mcmc_SGR_Chromis_use)$solutions[4,1]+summary(mcmc_SGR_Chromis_use)$solutions[6,1])*df_SGR_pred_full_Chromis[df_SGR_pred_full_Chromis$Treatment=="low",]$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SGR_pred_full_Chromis$SGR_staticM <- summary(mcmc_SGR_Chromis_useM)$solutions[1,1] + summary(mcmc_SGR_Chromis_useM)$solutions[3,1]*df_SGR_pred_full_Chromis$logMassGrowth


## Clownfish ##
mcmc_SGR_RR_Clownfish  <- MCMCglmm(logSGR ~ logMassCentGrowth + logMassMeanGrowth,
                                   random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian", 
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Clownfish))
# mcmc_SGR_RI_Clownfish  <- MCMCglmm(logSGR ~ logMassCentGrowth  + logMassMeanGrowth ,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_SGR_Clownfish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_SGR_Clownfish_use <- mcmc_SGR_RR_Clownfish,ifelse(mcmc_SGR_RR_Clownfish$DIC<mcmc_SGR_RI_Clownfish$DIC,mcmc_SGR_Clownfish_use <- mcmc_SGR_RR_Clownfish,mcmc_SGR_Clownfish_use <- mcmc_SGR_RI_Clownfish))

# Modeluse[7,10] <-ifelse(mcmc_SGR_Clownfish_use$DIC<mcmc_SGR_RI_Clownfish$DIC, "RR","RI")
# Modeluse[7,11] <- abs(mcmc_SGR_Clownfish_use$DIC-mcmc_SGR_RI_Clownfish$DIC)

df_SGR_pred_full_Clownfish <- cbind(Data_raw_SGR_Clownfish, predict(mcmc_SGR_Clownfish_use, marginal=NULL, interval = "prediction"))#finding predicted varible, using the best model

# finding the statick scaling value for a given individual at a given point for the given treatment
df_SGR_pred_full_Clownfish$SGR_staticT <- summary(mcmc_SGR_Clownfish_use)$solutions[1,1] + summary(mcmc_SGR_Clownfish_use)$solutions[3,1]*df_SGR_pred_full_Clownfish$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_SGR_pred_full_Clownfish$SGR_staticM <- summary(mcmc_SGR_Clownfish_use)$solutions[1,1] + summary(mcmc_SGR_Clownfish_use)$solutions[3,1]*df_SGR_pred_full_Clownfish$logMassGrowth


# Recombining #
Data_SGR <- bind_rows(df_SGR_pred_full_Zebrafish, df_SGR_pred_full_Rainbow_trout, df_SGR_pred_full_Brown_trout, df_SGR_pred_full_Cunner, df_SGR_pred_full_Guppy, df_SGR_pred_full_Chromis, df_SGR_pred_full_Clownfish)

names(Data_SGR)[names(Data_SGR) == "fit" ] <- "SGR_fit"
names(Data_SGR)[names(Data_SGR) == "lwr" ] <- "SGR_fit_low"
names(Data_SGR)[names(Data_SGR) == "upr" ] <- "SGR_fit_high"


### GE ###
#GE dataframe
Data_raw_GE <- Data_raw[Data_raw$GE>0,]
Data_raw_GE$GE_staticM <- NA
Data_raw_GE$GE_staticT <- NA

#Subsetting data
Data_raw_GE_Zebrafish <- subset(Data_raw_GE, Species=="Zebrafish")
Data_raw_GE_Rainbow_trout <- subset(Data_raw_GE, Species=="Rainbow_trout")
Data_raw_GE_Brown_trout <- subset(Data_raw_GE, Species=="Brown_trout")
Data_raw_GE_Cunner <- subset(Data_raw_GE, Species=="Cunner")
Data_raw_GE_Guppy <- subset(Data_raw_GE, Species=="Guppy",Sex != "n" & Sex != "d")
Data_raw_GE_Chromis <- subset(Data_raw_GE, Species=="Chromis")
Data_raw_GE_Clownfish <- subset(Data_raw_GE, Species=="Clownfish")

## Zebrafish ##
#With treatments
mcmc_GE_RR_Zebrafish  <- MCMCglmm(logGE ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
                                   random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Zebrafish))

# mcmc_GE_RI_Zebrafish  <- MCMCglmm(logGE ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Zebrafish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_GE_Zebrafish_use <- mcmc_GE_RR_Zebrafish,ifelse(mcmc_GE_RR_Zebrafish$DIC<mcmc_GE_RI_Zebrafish$DIC,mcmc_GE_Zebrafish_use <- mcmc_GE_RR_Zebrafish,mcmc_GE_Zebrafish_use <- mcmc_GE_RI_Zebrafish))

# Modeluse[1,12] <-ifelse(mcmc_GE_RR_Zebrafish$DIC<mcmc_GE_RI_Zebrafish$DIC, "RR","RI")
# Modeluse[1,13] <- abs(mcmc_GE_RR_Zebrafish$DIC-mcmc_GE_RI_Zebrafish$DIC)
#Without treatments
mcmc_GE_RR_ZebrafishM  <- MCMCglmm(logGE ~ logMassCentGrowth +logMassMeanGrowth,
                                    random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                    prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Zebrafish))

# mcmc_GE_RI_ZebrafishM  <- MCMCglmm(logGE ~ logMassCentGrowth +logMassMeanGrowth,
#                                     random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                     prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Zebrafish))
#Finding best model and noting it
ifelse(use_RRM==TRUE,mcmc_GE_Zebrafish_useM <- mcmc_GE_RR_ZebrafishM,ifelse(mcmc_GE_RR_ZebrafishM$DIC<mcmc_GE_RI_ZebrafishM$DIC,mcmc_GE_Zebrafish_useM <- mcmc_GE_RR_ZebrafishM,mcmc_GE_Zebrafish_useM <- mcmc_GE_RI_ZebrafishM))

# ModeluseM[1,12] <-ifelse(mcmc_GE_RR_ZebrafishM$DIC<mcmc_GE_RI_ZebrafishM$DIC, "RR","RI")
# ModeluseM[1,13] <- abs(mcmc_GE_RR_ZebrafishM$DIC-mcmc_GE_RI_ZebrafishM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_GE_pred_full_Zebrafish <- cbind(Data_raw_GE_Zebrafish, predict(mcmc_GE_Zebrafish_use, marginal=NULL, interval = "prediction")),
       df_GE_pred_full_Zebrafish <- cbind(Data_raw_GE_Zebrafish, predict(mcmc_GE_Zebrafish_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_GE_pred_full_Zebrafish[df_GE_pred_full_Zebrafish$Treatment=="H",]$GE_staticT <- df_GE_pred_full_Zebrafish[df_GE_pred_full_Zebrafish$Treatment=="H",]$logMassGrowth*summary(mcmc_GE_Zebrafish_use)$solutions[5,1]+summary(mcmc_GE_Zebrafish_use)$solutions[1,1]
df_GE_pred_full_Zebrafish[df_GE_pred_full_Zebrafish$Treatment=="L",]$GE_staticT <- df_GE_pred_full_Zebrafish[df_GE_pred_full_Zebrafish$Treatment=="L",]$logMassGrowth*(summary(mcmc_GE_Zebrafish_use)$solutions[5,1]+summary(mcmc_GE_Zebrafish_use)$solutions[8,1])+summary(mcmc_GE_Zebrafish_use)$solutions[1,1]+summary(mcmc_GE_Zebrafish_use)$solutions[3,1]
df_GE_pred_full_Zebrafish[df_GE_pred_full_Zebrafish$Treatment=="M",]$GE_staticT <- df_GE_pred_full_Zebrafish[df_GE_pred_full_Zebrafish$Treatment=="M",]$logMassGrowth*(summary(mcmc_GE_Zebrafish_use)$solutions[5,1]+summary(mcmc_GE_Zebrafish_use)$solutions[9,1])+summary(mcmc_GE_Zebrafish_use)$solutions[1,1]+summary(mcmc_GE_Zebrafish_use)$solutions[4,1]

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_GE_pred_full_Zebrafish$GE_staticM <- summary(mcmc_GE_Zebrafish_useM)$solutions[1,1] + summary(mcmc_GE_Zebrafish_useM)$solutions[3,1]*df_GE_pred_full_Zebrafish$logMassGrowth


## Rainbow trout ##
#with treatments
mcmc_GE_RR_Rainbow_trout  <- MCMCglmm(logGE ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
                                       random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                       prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                       verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Rainbow_trout))

# mcmc_GE_RI_Rainbow_trout  <- MCMCglmm(logGE ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
#                                        random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                        prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                        verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Rainbow_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_GE_Rainbow_trout_use <- mcmc_GE_RR_Rainbow_trout,ifelse(mcmc_GE_RR_Rainbow_trout$DIC<mcmc_GE_RI_Rainbow_trout$DIC,mcmc_GE_Rainbow_trout_use <- mcmc_GE_RR_Rainbow_trout,mcmc_GE_Rainbow_trout_use <- mcmc_GE_RI_Rainbow_trout))

# Modeluse[2,12] <-ifelse(mcmc_GE_RR_Rainbow_trout$DIC<mcmc_GE_RI_Rainbow_trout$DIC, "RR","RI")
# Modeluse[2,13] <- abs(mcmc_GE_RR_Rainbow_trout$DIC-mcmc_GE_RI_Rainbow_trout$DIC)
#Without treatments
mcmc_GE_RR_Rainbow_troutM  <- MCMCglmm(logGE ~ logMassCentGrowth  + logMassMeanGrowth  + Temp,
                                        random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                        prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                        verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Rainbow_trout))

# mcmc_GE_RI_Rainbow_troutM  <- MCMCglmm(logGE ~ logMassCentGrowth  + logMassMeanGrowth  + Temp,
#                                         random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                         prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                         verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Rainbow_trout))

ifelse(use_RRM==TRUE,mcmc_GE_Rainbow_trout_useM <- mcmc_GE_RR_Rainbow_troutM,ifelse(mcmc_GE_RR_Rainbow_troutM$DIC<mcmc_GE_RI_Rainbow_troutM$DIC,mcmc_GE_Rainbow_trout_useM <- mcmc_GE_RR_Rainbow_troutM,mcmc_GE_Rainbow_trout_useM <- mcmc_GE_RI_Rainbow_troutM))

# ModeluseM[2,12] <-ifelse(mcmc_GE_RR_Rainbow_troutM$DIC<mcmc_GE_RI_Rainbow_troutM$DIC, "RR","RI")
# ModeluseM[2,13] <- abs(mcmc_GE_RR_Rainbow_troutM$DIC-mcmc_GE_RI_Rainbow_troutM$DIC)

#finding predicted varible, using the best model with or without treatmetns, depending on use_treatment
ifelse(use_treatments==TRUE,df_GE_pred_full_Rainbow_trout <- cbind(Data_raw_GE_Rainbow_trout, predict(mcmc_GE_Rainbow_trout_use, marginal=NULL, interval = "prediction")),
       df_GE_pred_full_Rainbow_trout <- cbind(Data_raw_GE_Rainbow_trout, predict(mcmc_GE_Rainbow_trout_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_GE_pred_full_Rainbow_trout[df_GE_pred_full_Rainbow_trout$Treatment=="10Cegg",]$GE_staticT <- df_GE_pred_full_Rainbow_trout[df_GE_pred_full_Rainbow_trout$Treatment=="10Cegg",]$logMassGrowth*summary(mcmc_GE_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[1,1]
df_GE_pred_full_Rainbow_trout[df_GE_pred_full_Rainbow_trout$Treatment=="14Cegg",]$GE_staticT <- df_GE_pred_full_Rainbow_trout[df_GE_pred_full_Rainbow_trout$Treatment=="14Cegg",]$logMassGrowth*(summary(mcmc_GE_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[8,1])+summary(mcmc_GE_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[3,1]
df_GE_pred_full_Rainbow_trout[df_GE_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$GE_staticT <- df_GE_pred_full_Rainbow_trout[df_GE_pred_full_Rainbow_trout$Treatment=="14Cyolk",]$logMassGrowth*(summary(mcmc_GE_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[9,1])+summary(mcmc_GE_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[4,1]

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_GE_pred_full_Rainbow_trout$GE_staticM <- summary(mcmc_GE_Rainbow_trout_useM)$solutions[1,1] + summary(mcmc_GE_Rainbow_trout_useM)$solutions[3,1]*df_GE_pred_full_Rainbow_trout$logMassGrowth


## Brown trout ##
mcmc_GE_RR_Brown_trout  <- MCMCglmm(logGE ~ logMassCentGrowth + logMassMeanGrowth,
                                     random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian", 
                                     prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                     verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Brown_trout))
# mcmc_GE_RI_Brown_trout  <- MCMCglmm(logGE ~ logMassCentGrowth  + logMassMeanGrowth ,
#                                      random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                      prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                      verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Brown_trout))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_GE_Brown_trout_use <- mcmc_GE_RR_Brown_trout,ifelse(mcmc_GE_RR_Brown_trout$DIC<mcmc_GE_RI_Brown_trout$DIC,mcmc_GE_Brown_trout_use <- mcmc_GE_RR_Brown_trout,mcmc_GE_Brown_trout_use <- mcmc_GE_RI_Brown_trout))

# Modeluse[3,12] <-ifelse(mcmc_GE_RR_Brown_trout$DIC<mcmc_GE_RI_Brown_trout$DIC, "RR","RI")
# Modeluse[3,13] <- abs(mcmc_GE_RR_Brown_trout$DIC-mcmc_GE_RI_Brown_trout$DIC)

df_GE_pred_full_Brown_trout <- cbind(Data_raw_GE_Brown_trout, predict(mcmc_GE_Brown_trout_use, marginal=NULL, interval = "prediction"))#finding predicted varible, using the best model

# finding the statick scaling value for a given individual at a given point for the given treatment
df_GE_pred_full_Brown_trout$GE_staticT <- summary(mcmc_GE_Brown_trout_use)$solutions[1,1] + summary(mcmc_GE_Brown_trout_use)$solutions[3,1]*df_GE_pred_full_Brown_trout$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_GE_pred_full_Brown_trout$GE_staticM <- summary(mcmc_GE_Brown_trout_use)$solutions[1,1] + summary(mcmc_GE_Brown_trout_use)$solutions[3,1]*df_GE_pred_full_Brown_trout$logMassGrowth


## Cunner ##
mcmc_GE_RR_Cunner  <- MCMCglmm(logGE ~ logMassCentGrowth + logMassMeanGrowth,
                                random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian", 
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Cunner))
# mcmc_GE_RI_Cunner  <- MCMCglmm(logGE ~ logMassCentGrowth  + logMassMeanGrowth ,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Cunner))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_GE_Cunner_use <- mcmc_GE_RR_Cunner,ifelse(mcmc_GE_RR_Cunner$DIC<mcmc_GE_RI_Cunner$DIC,mcmc_GE_Cunner_use <- mcmc_GE_RR_Cunner,mcmc_GE_Cunner_use <- mcmc_GE_RI_Cunner))

# Modeluse[4,12] <-ifelse(mcmc_GE_RR_Cunner$DIC<mcmc_GE_RI_Cunner$DIC, "RR","RI")
# Modeluse[4,13] <- abs(mcmc_GE_RR_Cunner$DIC-mcmc_GE_RI_Cunner$DIC)

df_GE_pred_full_Cunner <- cbind(Data_raw_GE_Cunner, predict(mcmc_GE_Cunner_use, marginal=NULL, interval = "prediction"))#finding predicted varible, using the best model

# finding the statick scaling value for a given individual at a given point for the given treatment
df_GE_pred_full_Cunner$GE_staticT <- summary(mcmc_GE_Cunner_use)$solutions[1,1] + summary(mcmc_GE_Cunner_use)$solutions[3,1]*df_GE_pred_full_Cunner$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_GE_pred_full_Cunner$GE_staticM <- summary(mcmc_GE_Cunner_use)$solutions[1,1] + summary(mcmc_GE_Cunner_use)$solutions[3,1]*df_GE_pred_full_Cunner$logMassGrowth

## Guppy ##
#with sex
mcmc_GE_RR_Guppy  <- MCMCglmm(logGE ~ logMassCentGrowth * Sex + logMassMeanGrowth * Sex,
                               random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                               prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                               verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Guppy))

# mcmc_GE_RI_Guppy  <- MCMCglmm(logGE ~ logMassCentGrowth * Sex + logMassMeanGrowth * Sex,
#                                random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Guppy))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_GE_Guppy_use <- mcmc_GE_RR_Guppy,ifelse(mcmc_GE_RR_Guppy$DIC<mcmc_GE_RI_Guppy$DIC,mcmc_GE_Guppy_use <- mcmc_GE_RR_Guppy,mcmc_GE_Guppy_use <- mcmc_GE_RI_Guppy))

# Modeluse[5,12] <-ifelse(mcmc_GE_RR_Guppy$DIC<mcmc_GE_RI_Guppy$DIC, "RR","RI")
# Modeluse[5,13] <- abs(mcmc_GE_RR_Guppy$DIC-mcmc_GE_RI_Guppy$DIC)


#without sex
mcmc_GE_RR_GuppyM  <- MCMCglmm(logGE ~ logMassCentGrowth + logMassMeanGrowth,
                                random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Guppy))

# mcmc_GE_RI_GuppyM  <- MCMCglmm(logGE ~ logMassCentGrowth + logMassMeanGrowth,
#                                 random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                 prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Guppy))
#Finding best model and noting it
ifelse(use_RRM==TRUE,mcmc_GE_Guppy_useM <- mcmc_GE_RR_GuppyM,ifelse(mcmc_GE_RR_GuppyM$DIC<mcmc_GE_RI_GuppyM$DIC,mcmc_GE_Guppy_useM <- mcmc_GE_RR_GuppyM,mcmc_GE_Guppy_useM <- mcmc_GE_RI_GuppyM))

# ModeluseM[5,12] <-ifelse(mcmc_GE_RR_GuppyM$DIC<mcmc_GE_RI_GuppyM$DIC, "RR","RI")
# ModeluseM[5,13] <- abs(mcmc_GE_RR_GuppyM$DIC-mcmc_GE_RI_GuppyM$DIC)



ifelse(use_treatments==TRUE,df_GE_pred_full_Guppy <- cbind(Data_raw_GE_Guppy, predict(mcmc_GE_Guppy_use, marginal=NULL, interval = "prediction")),
       df_GE_pred_full_Guppy <- cbind(Data_raw_GE_Guppy, predict(mcmc_GE_Guppy_useM, marginal=NULL, interval = "prediction"))) 

# finding the statick scaling value for a given individual at a given point for the given treatment
df_GE_pred_full_Guppy[df_GE_pred_full_Guppy$Sex=="f",]$GE_staticT <- summary(mcmc_GE_Guppy_use)$solutions[1,1] + summary(mcmc_GE_Guppy_use)$solutions[4,1]*df_GE_pred_full_Guppy[df_GE_pred_full_Guppy$Sex=="f",]$logMassGrowth
df_GE_pred_full_Guppy[df_GE_pred_full_Guppy$Sex=="m",]$GE_staticT <- summary(mcmc_GE_Guppy_use)$solutions[1,1]+summary(mcmc_GE_Guppy_use)$solutions[3,1] + (summary(mcmc_GE_Guppy_use)$solutions[4,1]+summary(mcmc_GE_Guppy_use)$solutions[6,1])*df_GE_pred_full_Guppy[df_GE_pred_full_Guppy$Sex=="m",]$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_GE_pred_full_Guppy$GE_staticM <- summary(mcmc_GE_Guppy_useM)$solutions[1,1] + summary(mcmc_GE_Guppy_useM)$solutions[3,1]*df_GE_pred_full_Guppy$logMassGrowth


## Chromis ##
#with treatment
mcmc_GE_RR_Chromis  <- MCMCglmm(logGE ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
                                 random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                 prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                 verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Chromis))

# mcmc_GE_RI_Chromis  <- MCMCglmm(logGE ~ logMassCentGrowth * Treatment + logMassMeanGrowth * Treatment,
#                                  random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                  prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Chromis))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_GE_Chromis_use <- mcmc_GE_RR_Chromis,ifelse(mcmc_GE_RR_Chromis$DIC<mcmc_GE_RI_Chromis$DIC,mcmc_GE_Chromis_use <- mcmc_GE_RR_Chromis,mcmc_GE_Chromis_use <- mcmc_GE_RI_Chromis))

# Modeluse[6,12] <-ifelse(mcmc_GE_RR_Chromis$DIC<mcmc_GE_RI_Chromis$DIC, "RR","RI")
# Modeluse[6,13] <- abs(mcmc_GE_RR_Chromis$DIC-mcmc_GE_RI_Chromis$DIC)

#without treatment
mcmc_GE_RR_ChromisM  <- MCMCglmm(logGE ~ logMassCentGrowth + logMassMeanGrowth,
                                  random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian",
                                  prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                  verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Chromis))

# mcmc_GE_RI_ChromisM  <- MCMCglmm(logGE ~ logMassCentGrowth + logMassMeanGrowth,
#                                   random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                   prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Chromis))
#Finding best model and noting it
ifelse(use_RRM==TRUE,mcmc_GE_Chromis_useM <- mcmc_GE_RR_ChromisM,ifelse(mcmc_GE_RR_ChromisM$DIC<mcmc_GE_RI_ChromisM$DIC,mcmc_GE_Chromis_useM <- mcmc_GE_RR_ChromisM,mcmc_GE_Chromis_useM <- mcmc_GE_RI_ChromisM))

# ModeluseM[6,12] <-ifelse(mcmc_GE_RR_ChromisM$DIC<mcmc_GE_RI_ChromisM$DIC, "RR","RI")
# ModeluseM[6,13] <- abs(mcmc_GE_RR_ChromisM$DIC-mcmc_GE_RI_ChromisM$DIC)



ifelse(use_treatments==TRUE,df_GE_pred_full_Chromis <- cbind(Data_raw_GE_Chromis, predict(mcmc_GE_Chromis_use, marginal=NULL, interval = "prediction")),
       df_GE_pred_full_Chromis <- cbind(Data_raw_GE_Chromis, predict(mcmc_GE_Chromis_useM, marginal=NULL, interval = "prediction"))) 


# finding the statick scaling value for a given individual at a given point for the given treatment
df_GE_pred_full_Chromis[df_GE_pred_full_Chromis$Treatment=="high",]$GE_staticT <- summary(mcmc_GE_Chromis_use)$solutions[1,1] + summary(mcmc_GE_Chromis_use)$solutions[4,1]*df_GE_pred_full_Chromis[df_GE_pred_full_Chromis$Treatment=="high",]$logMassGrowth
df_GE_pred_full_Chromis[df_GE_pred_full_Chromis$Treatment=="low",]$GE_staticT <- summary(mcmc_GE_Chromis_use)$solutions[1,1]+summary(mcmc_GE_Chromis_use)$solutions[3,1] + (summary(mcmc_GE_Chromis_use)$solutions[4,1]+summary(mcmc_GE_Chromis_use)$solutions[6,1])*df_GE_pred_full_Chromis[df_GE_pred_full_Chromis$Treatment=="low",]$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_GE_pred_full_Chromis$GE_staticM <- summary(mcmc_GE_Chromis_useM)$solutions[1,1] + summary(mcmc_GE_Chromis_useM)$solutions[3,1]*df_GE_pred_full_Chromis$logMassGrowth


## Clownfish ##
mcmc_GE_RR_Clownfish  <- MCMCglmm(logGE ~ logMassCentGrowth + logMassMeanGrowth,
                                   random = ~ us(1 + logMassCentGrowth):FishID, rcov = ~ units, family = "gaussian", 
                                   prior = prior_1t_RR_G1, nitt = 350000, burnin = 50000, thin = 40,
                                   verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Clownfish))
# mcmc_GE_RI_Clownfish  <- MCMCglmm(logGE ~ logMassCentGrowth  + logMassMeanGrowth ,
#                                    random = ~ FishID, rcov = ~ units, family = "gaussian",
#                                    prior = prior_1t_G1, nitt = 350000, burnin = 50000, thin = 40,
#                                    verbose = TRUE, pr = TRUE, data = as.data.frame(Data_raw_GE_Clownfish))
#Finding best model and noting it
ifelse(use_RR==TRUE,mcmc_GE_Clownfish_use <- mcmc_GE_RR_Clownfish,ifelse(mcmc_GE_RR_Clownfish$DIC<mcmc_GE_RI_Clownfish$DIC,mcmc_GE_Clownfish_use <- mcmc_GE_RR_Clownfish,mcmc_GE_Clownfish_use <- mcmc_GE_RI_Clownfish))

# Modeluse[7,12] <-ifelse(mcmc_GE_Clownfish_use$DIC<mcmc_GE_RI_Clownfish$DIC, "RR","RI")
# Modeluse[7,13] <- abs(mcmc_GE_Clownfish_use$DIC-mcmc_GE_RI_Clownfish$DIC)

df_GE_pred_full_Clownfish <- cbind(Data_raw_GE_Clownfish, predict(mcmc_GE_Clownfish_use, marginal=NULL, interval = "prediction"))#finding predicted varible, using the best model

# finding the statick scaling value for a given individual at a given point for the given treatment
df_GE_pred_full_Clownfish$GE_staticT <- summary(mcmc_GE_Clownfish_use)$solutions[1,1] + summary(mcmc_GE_Clownfish_use)$solutions[3,1]*df_GE_pred_full_Clownfish$logMassGrowth

#Finding the Avg statick scaling for a given species, disregaringing treatment
df_GE_pred_full_Clownfish$GE_staticM <- summary(mcmc_GE_Clownfish_use)$solutions[1,1] + summary(mcmc_GE_Clownfish_use)$solutions[3,1]*df_GE_pred_full_Clownfish$logMassGrowth


# Recombining #
Data_GE <- bind_rows(df_GE_pred_full_Zebrafish, df_GE_pred_full_Rainbow_trout, df_GE_pred_full_Brown_trout, df_GE_pred_full_Cunner, df_GE_pred_full_Guppy, df_GE_pred_full_Chromis, df_GE_pred_full_Clownfish)

names(Data_GE)[names(Data_GE) == "fit" ] <- "GE_fit"
names(Data_GE)[names(Data_GE) == "lwr" ] <- "GE_fit_low"
names(Data_GE)[names(Data_GE) == "upr" ] <- "GE_fit_high"


#Merging all parameter dataframes
Data_all <- merge(Data_raw,Data_SMR, all=TRUE)
Data_all <- merge(Data_all,Data_MMR, all=TRUE)
Data_all <- merge(Data_all,Data_AS, all=TRUE)
Data_all <- merge(Data_all,Data_AGR, all=TRUE)
Data_all <- merge(Data_all,Data_SGR, all=TRUE)
Data_all <- merge(Data_all,Data_GE, all=TRUE)
Data_all <- merge(Data_all, Data_FAS, all=TRUE)

Data_fit <- Data_all




## finding temperatur adjustet values and evolutionary scaling of parameters
# AGR model with all data to analyse across-species (evolutionary) scaling
mcmc_AGR_RR_Evol  <- MCMCglmm(logAGR ~ logMassGrowth + Temp,
                              random = ~ us(1 + logMassGrowth):Group,
                              rcov = ~ units,
                              family = "gaussian",
                              prior = prior_1t_RR_G1,
                              nitt = 350000,
                              burnin = 50000,
                              thin = 40,
                              verbose = TRUE,
                              pr = TRUE,
                              data = as.data.frame(Data_fit[!is.na(Data_fit$logMassGrowth) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMassGrowth) & !is.infinite(Data_fit$logMassGrowth) & !is.na(Data_fit$logAGR)&!is.nan(Data_fit$logAGR) & !is.infinite(Data_fit$logAGR),]))

#summary(mcmc_AGR_RR_Evol)

Data_fit <- Data_fit[!is.na(Data_fit$Species),]


# Q10 adjust raw data using Q10 from model
Data_fit$AGR_TempAdj <- Data_fit$AGR * ((10^(summary(mcmc_AGR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust model fits using Q10 from model
Data_fit$AGR_fit_TempAdj <- 10^(Data_fit$AGR_fit) * ((10^(summary(mcmc_AGR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust static scaling fits using Q10 from model
Data_fit$AGR_staticM_TempAdj <- 10^(Data_fit$AGR_staticM) * ((10^(summary(mcmc_AGR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))



# model with temperature-adjusted data to extract temperature-independent intercept for plotting
mcmc_AGR_RR_Evol_TempAdj  <- MCMCglmm(log10(AGR_TempAdj) ~ logMassGrowth,
                                      random = ~ us(1 + logMassGrowth):Group,
                                      rcov = ~ units,
                                      family = "gaussian",
                                      prior = prior_1t_RR_G1,
                                      nitt = 350000,
                                      burnin = 50000,
                                      thin = 40,
                                      verbose = TRUE,
                                      pr = TRUE,
                                      data = as.data.frame(Data_fit[!is.na(Data_fit$logMassGrowth) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMassGrowth) & !is.infinite(Data_fit$logMassGrowth) & !is.na(Data_fit$logAGR)&!is.nan(Data_fit$logAGR) & !is.infinite(Data_fit$logAGR),]))

# SGR model with all data to analyse across-species (evolutionary) scaling
mcmc_SGR_RR_Evol  <- MCMCglmm(logSGR ~ logMassGrowth + Temp,
                              random = ~ us(1 + logMassGrowth):Group,
                              rcov = ~ units,
                              family = "gaussian",
                              prior = prior_1t_RR_G1,
                              nitt = 350000,
                              burnin = 50000,
                              thin = 40,
                              verbose = TRUE,
                              pr = TRUE,
                              data = as.data.frame(Data_fit[!is.na(Data_fit$logMassGrowth) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMassGrowth) & !is.infinite(Data_fit$logMassGrowth) & !is.na(Data_fit$logSGR)&!is.nan(Data_fit$logSGR) & !is.infinite(Data_fit$logSGR),]))



# Q10 adjust raw data using Q10 from model
Data_fit$SGR_TempAdj <- Data_fit$SGR * ((10^(summary(mcmc_SGR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust model fits using Q10 from model
Data_fit$SGR_fit_TempAdj <- 10^(Data_fit$SGR_fit) * ((10^(summary(mcmc_SGR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust static scaling fits using Q10 from model
Data_fit$SGR_staticM_TempAdj <- 10^(Data_fit$SGR_staticM) * ((10^(summary(mcmc_SGR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))



# model with temperature-adjusted data to extract temperature-independent intercept for plotting
mcmc_SGR_RR_Evol_TempAdj  <- MCMCglmm(log10(SGR_TempAdj) ~ logMassGrowth,
                                      random = ~ us(1 + logMassGrowth):Group,
                                      rcov = ~ units,
                                      family = "gaussian",
                                      prior = prior_1t_RR_G1,
                                      nitt = 350000,
                                      burnin = 50000,
                                      thin = 40,
                                      verbose = TRUE,
                                      pr = TRUE,
                                      data = as.data.frame(Data_fit[!is.na(Data_fit$logMassGrowth) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMassGrowth) & !is.infinite(Data_fit$logMassGrowth) & !is.na(Data_fit$logSGR)&!is.nan(Data_fit$logSGR) & !is.infinite(Data_fit$logSGR),]))

# GE model with all data to analyse across-species (evolutionary) scaling
mcmc_GE_RR_Evol  <- MCMCglmm(logGE ~ logMassGrowth + Temp,
                              random = ~ us(1 + logMassGrowth):Group,
                              rcov = ~ units,
                              family = "gaussian",
                              prior = prior_1t_RR_G1,
                              nitt = 350000,
                              burnin = 50000,
                              thin = 40,
                              verbose = TRUE,
                              pr = TRUE,
                              data = as.data.frame(Data_fit[!is.na(Data_fit$logMassGrowth) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMassGrowth) & !is.infinite(Data_fit$logMassGrowth) & !is.na(Data_fit$logGE)&!is.nan(Data_fit$logGE) & !is.infinite(Data_fit$logGE),]))

#summary(mcmc_GE_RR_Evol)

Data_fit <- Data_fit[!is.na(Data_fit$Species),]


# Q10 adjust raw data using Q10 from model
Data_fit$GE_TempAdj <- Data_fit$GE * ((10^(summary(mcmc_GE_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust model fits using Q10 from model
Data_fit$GE_fit_TempAdj <- 10^(Data_fit$GE_fit) * ((10^(summary(mcmc_GE_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust static scaling fits using Q10 from model
Data_fit$GE_staticM_TempAdj <- 10^(Data_fit$GE_staticM) * ((10^(summary(mcmc_GE_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))



# model with temperature-adjusted data to extract temperature-independent intercept for plotting
mcmc_GE_RR_Evol_TempAdj  <- MCMCglmm(log10(GE_TempAdj) ~ logMassGrowth,
                                      random = ~ us(1 + logMassGrowth):Group,
                                      rcov = ~ units,
                                      family = "gaussian",
                                      prior = prior_1t_RR_G1,
                                      nitt = 350000,
                                      burnin = 50000,
                                      thin = 40,
                                      verbose = TRUE,
                                      pr = TRUE,
                                      data = as.data.frame(Data_fit[!is.na(Data_fit$logMassGrowth) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMassGrowth) & !is.infinite(Data_fit$logMassGrowth) & !is.na(Data_fit$logGE)&!is.nan(Data_fit$logGE) & !is.infinite(Data_fit$logGE),]))

# SMR model with all data to analyse across-species (evolutionary) scaling
mcmc_SMR_RR_Evol  <- MCMCglmm(logSMR ~ logMass + Temp,
                              random = ~ us(1 + logMass):Group,
                              rcov = ~ units,
                              family = "gaussian",
                              prior = prior_1t_RR_G1,
                              nitt = 350000,
                              burnin = 50000,
                              thin = 40,
                              verbose = TRUE,
                              pr = TRUE,
                              data = as.data.frame(Data_fit[!is.na(Data_fit$logMass) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMass) & !is.infinite(Data_fit$logMass) & !is.na(Data_fit$logSMR)&!is.nan(Data_fit$logSMR) & !is.infinite(Data_fit$logSMR),]))

#summary(mcmc_SMR_RR_Evol_TempAdj)

Data_fit <- Data_fit[!is.na(Data_fit$Species),]


# Q10 adjust raw data using Q10 from model
Data_fit$SMR_TempAdj <- Data_fit$SMR * ((10^(summary(mcmc_SMR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust model fits using Q10 from model
Data_fit$SMR_fit_TempAdj <- 10^(Data_fit$SMR_fit) * ((10^(summary(mcmc_SMR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust static scaling fits using Q10 from model
Data_fit$SMR_staticM_TempAdj <- 10^(Data_fit$SMR_staticM) * ((10^(summary(mcmc_SMR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))



# model with temperature-adjusted data to extract temperature-independent intercept for plotting
mcmc_SMR_RR_Evol_TempAdj  <- MCMCglmm(log10(SMR_TempAdj) ~ logMass,
                                      random = ~ us(1 + logMass):Group,
                                      rcov = ~ units,
                                      family = "gaussian",
                                      prior = prior_1t_RR_G1,
                                      nitt = 350000,
                                      burnin = 50000,
                                      thin = 40,
                                      verbose = TRUE,
                                      pr = TRUE,
                                      data = as.data.frame(Data_fit[!is.na(Data_fit$logMass) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMass) & !is.infinite(Data_fit$logMass) & !is.na(Data_fit$logSMR)&!is.nan(Data_fit$logSMR) & !is.infinite(Data_fit$logSMR),]))

# MMR model with all data to analyse across-species (evolutionary) scaling
mcmc_MMR_RR_Evol  <- MCMCglmm(logMMR ~ logMass + Temp,
                              random = ~ us(1 + logMass):Group,
                              rcov = ~ units,
                              family = "gaussian",
                              prior = prior_1t_RR_G1,
                              nitt = 350000,
                              burnin = 50000,
                              thin = 40,
                              verbose = TRUE,
                              pr = TRUE,
                              data = as.data.frame(Data_fit[!is.na(Data_fit$logMass) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMass) & !is.infinite(Data_fit$logMass) & !is.na(Data_fit$logMMR)&!is.nan(Data_fit$logMMR) & !is.infinite(Data_fit$logMMR),]))




# Q10 adjust raw data using Q10 from model
Data_fit$MMR_TempAdj <- Data_fit$MMR * ((10^(summary(mcmc_MMR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust model fits using Q10 from model
Data_fit$MMR_fit_TempAdj <- 10^(Data_fit$MMR_fit) * ((10^(summary(mcmc_MMR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust static scaling fits using Q10 from model
Data_fit$MMR_staticM_TempAdj <- 10^(Data_fit$MMR_staticM) * ((10^(summary(mcmc_MMR_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))



# model with temperature-adjusted data to extract temperature-independent intercept for plotting
mcmc_MMR_RR_Evol_TempAdj  <- MCMCglmm(log10(MMR_TempAdj) ~ logMass,
                                      random = ~ us(1 + logMass):Group,
                                      rcov = ~ units,
                                      family = "gaussian",
                                      prior = prior_1t_RR_G1,
                                      nitt = 350000,
                                      burnin = 50000,
                                      thin = 40,
                                      verbose = TRUE,
                                      pr = TRUE,
                                      data = as.data.frame(Data_fit[!is.na(Data_fit$logMass) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMass) & !is.infinite(Data_fit$logMass) & !is.na(Data_fit$logMMR)&!is.nan(Data_fit$logMMR) & !is.infinite(Data_fit$logMMR),]))

# FAS model with all data to analyse across-species (evolutionary) scaling
mcmc_FAS_RR_Evol  <- MCMCglmm(logFAS ~ logMass + Temp,
                              random = ~ us(1 + logMass):Group,
                              rcov = ~ units,
                              family = "gaussian",
                              prior = prior_1t_RR_G1,
                              nitt = 350000,
                              burnin = 50000,
                              thin = 40,
                              verbose = TRUE,
                              pr = TRUE,
                              data = as.data.frame(Data_fit[!is.na(Data_fit$logMass) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMass) & !is.infinite(Data_fit$logMass) & !is.na(Data_fit$logFAS)&!is.nan(Data_fit$logFAS) & !is.infinite(Data_fit$logFAS),]))




# Q10 adjust raw data using Q10 from model
Data_fit$FAS_TempAdj <- Data_fit$FAS * ((10^(summary(mcmc_FAS_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust model fits using Q10 from model
Data_fit$FAS_fit_TempAdj <- 10^(Data_fit$FAS_fit) * ((10^(summary(mcmc_FAS_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust static scaling fits using Q10 from model
Data_fit$FAS_staticM_TempAdj <- 10^(Data_fit$FAS_staticM) * ((10^(summary(mcmc_FAS_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))



# model with temperature-adjusted data to extract temperature-independent intercept for plotting
mcmc_FAS_RR_Evol_TempAdj  <- MCMCglmm(log10(FAS_TempAdj) ~ logMass,
                                      random = ~ us(1 + logMass):Group,
                                      rcov = ~ units,
                                      family = "gaussian",
                                      prior = prior_1t_RR_G1,
                                      nitt = 350000,
                                      burnin = 50000,
                                      thin = 40,
                                      verbose = TRUE,
                                      pr = TRUE,
                                      data = as.data.frame(Data_fit[!is.na(Data_fit$logMass) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMass) & !is.infinite(Data_fit$logMass) & !is.na(Data_fit$logFAS)&!is.nan(Data_fit$logFAS) & !is.infinite(Data_fit$logFAS),]))

# AS model with all data to analyse across-species (evolutionary) scaling
mcmc_AS_RR_Evol  <- MCMCglmm(logAS ~ logMass + Temp,
                             random = ~ us(1 + logMass):Group,
                             rcov = ~ units,
                             family = "gaussian",
                             prior = prior_1t_RR_G1,
                             nitt = 350000,
                             burnin = 50000,
                             thin = 40,
                             verbose = TRUE,
                             pr = TRUE,
                             data = as.data.frame(Data_fit[!is.na(Data_fit$logMass) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMass) & !is.infinite(Data_fit$logMass) & !is.na(Data_fit$logAS)&!is.nan(Data_fit$logAS) & !is.infinite(Data_fit$logAS),]))




# Q10 adjust raw data using Q10 from model
Data_fit$AS_TempAdj <- Data_fit$AS * ((10^(summary(mcmc_AS_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust model fits using Q10 from model
Data_fit$AS_fit_TempAdj <- 10^(Data_fit$AS_fit) * ((10^(summary(mcmc_AS_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))

# Q10 adjust static scaling fits using Q10 from model
Data_fit$AS_staticM_TempAdj <- 10^(Data_fit$AS_staticM) * ((10^(summary(mcmc_AS_RR_Evol)$solutions[3,1] * 10))^((20-Data_fit$Temp)/10))



# model with temperature-adjusted data to extract temperature-independent intercept for plotting
mcmc_AS_RR_Evol_TempAdj  <- MCMCglmm(log10(AS_TempAdj) ~ logMass,
                                     random = ~ us(1 + logMass):Group,
                                     rcov = ~ units,
                                     family = "gaussian",
                                     prior = prior_1t_RR_G1,
                                     nitt = 350000,
                                     burnin = 50000,
                                     thin = 40,
                                     verbose = TRUE,
                                     pr = TRUE,
                                     data = as.data.frame(Data_fit[!is.na(Data_fit$logMass) & !is.na(Data_fit$Temp) & !is.nan(Data_fit$logMass) & !is.infinite(Data_fit$logMass) & !is.na(Data_fit$logAS)&!is.nan(Data_fit$logAS) & !is.infinite(Data_fit$logAS),]))





####### finding individual scaling for SMR, MMR, AS, AGR, SGR and GE and FAS #########
Results_fit <- data.frame(matrix(nrow=0,ncol=36))
names(Results_fit) <- c("Species","FishID","Treatment", "Start_weight","End_weight", "Fold_change","Number_of_messurements","SGR", "AGR","Slope_SMR_fit","p_vaule_SMR_fit","Intersept_SMR_fit",
                        "Slope_MMR_fit","p_vaule_MMR_fit","Intersept_MMR_fit","Slope_AS_fit","p_vaule_AS_fit","Intersept_AS_fit", "Slope_SGR_fit","p_vaule_SGR_fit",
                        "Intersept_SGR_fit","Slope_AGR_fit","p_vaule_AGR_fit","Intersept_AGR_fit","Slope_GE_fit","p_vaule_GE_fit","Intersept_GE_fit","Slope_FAS_fit","p_value_FAS_fit","Intersept_FAS_fit","Plot_treatment","Plot_group", "SMR_end","Mean_w",
                        "SMR_mean_w","Temp")
Data_fit <-Data_fit[!duplicated(Data_fit[,c(1:4)]),]
#ID vector
for (x in 1:length(unique(Data_fit$Group))) {
  Data_ID <- Data_fit[Data_fit$Group==unique(Data_fit$Group)[x],] #a dataframe that contains all messurements of a given individual
  if( max(Data_ID$Mass)/min(Data_ID$Mass)>1 ){ #only uses individuals that has grown in size
    Mass_start <- Data_ID[which.min(Data_ID$Measurement),6]
    Mass_end <- Data_ID[which.max(Data_ID$Measurement),6]
    SGR <- (log(Mass_end)-log(Mass_start))/(max(Data_ID$Day)-min(Data_ID$Day))*100 #overall SGR
    AGR <- (Mass_end-Mass_start)/(max(Data_ID$Day)-min(Data_ID$Day)) #overall AGR
    if(length(Data_ID[!is.na(Data_ID$SMR_fit),1])>2){slope_SMR <-lm(log10(Data_ID$SMR_fit_TempAdj) ~Data_ID$logMassCent) }
    if(length(Data_ID[!is.na(Data_ID$MMR_fit),1])>2){slope_MMR <-lm(log10(Data_ID$MMR_fit_TempAdj) ~Data_ID$logMassCent) }
    if(length(Data_ID[!is.na(Data_ID$AS_fit),1])>2){slope_AS <-lm(log10(Data_ID$AS_fit_TempAdj) ~Data_ID$logMassCent) }
    Data_ID_1 <- Data_ID[!is.na(Data_ID$logMassGrowth),] #removingt datapoints with NA growth, ia fisrt messurement as growth is always messured from the previus messurement to the current
    #only finds the slope of growth paremeters for individuals that has more then 2 datapoints
    if(length(Data_ID_1[!is.na(Data_ID_1$AGR_fit),1])>2){slope_AGR <-lm(data=Data_ID_1[!is.na(Data_ID_1$AGR_fit),], log10(AGR_fit_TempAdj) ~logMassCentGrowth) }
    if(length(Data_ID_1[!is.na(Data_ID_1$SGR_fit),1])>2){slope_SGR <-lm(data=Data_ID_1[!is.na(Data_ID_1$SGR_fit),], log10(SGR_fit_TempAdj) ~logMassCentGrowth) }
    if(length(Data_ID_1[!is.na(Data_ID_1$GE_fit),1])>2){slope_GE <-lm(data=Data_ID_1[!is.na(Data_ID_1$GE_fit),], log10(GE_fit_TempAdj) ~logMassCentGrowth) }
    if(length(Data_ID_1[!is.na(Data_ID_1$FAS_fit),1])>2){slope_FAS <-lm(data=Data_ID_1[!is.na(Data_ID_1$FAS_fit),], log10(FAS_fit_TempAdj) ~logMassCent) }
    Results_fit[x,1] <- Data_ID$Species[1]
    Results_fit[x,2] <- Data_ID$FishID[1]
    Results_fit[x,3] <- ifelse(!Data_ID$Species[1]=="Guppy",Results_fit[x,3] <- Data_ID$Treatment[1],Results_fit[x,3] <- Data_ID$Sex[1])
    Results_fit[x,4] <- Mass_start
    Results_fit[x,5] <- Mass_end
    Results_fit[x,6] <- Mass_end/Mass_start
    Results_fit[x,7] <- length(Data_ID$Measurement)
    Results_fit[x,8] <- SGR
    Results_fit[x,9] <- AGR
    Results_fit[x,10] <- slope_SMR$coefficients[2]
    Results_fit[x,11] <- summary(slope_SMR)$coefficients[2,4]
    Results_fit[x,12] <- slope_SMR$coefficients[1]
    ifelse(length(Data_ID[!is.na(Data_ID$MMR_fit),1])>2, Results_fit[x,13] <- slope_MMR$coefficients[2], Results_fit[x,13] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$MMR_fit),1])>2, Results_fit[x,14] <- summary(slope_MMR)$coefficients[2,4], Results_fit[x,14] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$MMR_fit),1])>2, Results_fit[x,15] <- slope_MMR$coefficients[1], Results_fit[x,15] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$AS_fit),1])>2, Results_fit[x,16] <- slope_AS$coefficients[2], Results_fit[x,16] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$AS_fit),1])>2, Results_fit[x,17] <- summary(slope_AS)$coefficients[2,4], Results_fit[x,17] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$AS_fit),1])>2, Results_fit[x,18] <- slope_AS$coefficients[1], Results_fit[x,18] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$SGR_fit),1])>2, Results_fit[x,19] <- slope_SGR$coefficients[2], Results_fit[x,19] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$SGR_fit),1])>2, Results_fit[x,20] <- summary(slope_SGR)$coefficients[2,4], Results_fit[x,20] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$SGR_fit),1])>2, Results_fit[x,21] <- slope_SGR$coefficients[1], Results_fit[x,21] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$AGR_fit),1])>2, Results_fit[x,22] <- slope_AGR$coefficients[2], Results_fit[x,22] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$AGR_fit),1])>2, Results_fit[x,23] <- summary(slope_AGR)$coefficients[2,4], Results_fit[x,23] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$AGR_fit),1])>2, Results_fit[x,24] <- slope_AGR$coefficients[1], Results_fit[x,24] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$GE_fit),1])>2, Results_fit[x,25] <- slope_GE$coefficients[2], Results_fit[x,25] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$GE_fit),1])>2, Results_fit[x,26] <- summary(slope_GE)$coefficients[2,4], Results_fit[x,26] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$GE_fit),1])>2, Results_fit[x,27] <- slope_GE$coefficients[1], Results_fit[x,27] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$FAS_fit),1])>2, Results_fit[x,28] <- slope_FAS$coefficients[2], Results_fit[x,28] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$FAS_fit),1])>2, Results_fit[x,29] <- summary(slope_FAS)$coefficients[2,4], Results_fit[x,29] <- NA )
    ifelse(length(Data_ID[!is.na(Data_ID$FAS_fit),1])>2, Results_fit[x,30] <- slope_FAS$coefficients[1], Results_fit[x,30] <- NA )

    # Adding simular treatment ID's, so each species only has the treatments 1, 2 and 3. It is used later for plots
    if(Results_fit[x,1]=="Cunner"|Results_fit[x,1]=="Brown_trout"|Results_fit[x,1]=="Clownfish"){Results_fit[x,31] <- 1}
    if(Results_fit[x,1]=="Rainbow_trout"){ifelse(Results_fit[x,3]=="10Cegg",Results_fit[x,31] <- 1,ifelse(Results_fit[x,3]=="14Cegg",Results_fit[x,31] <- 2,Results_fit[x,31] <- 3))}
    if(Results_fit[x,1]=="Zebrafish"){ifelse(Results_fit[x,3]=="H",Results_fit[x,31] <- 1,ifelse(Results_fit[x,3]=="L",Results_fit[x,31] <- 2,Results_fit[x,31] <- 3))}
    if(Results_fit[x,1]=="Chromis"){ifelse(Results_fit[x,3]=="high",Results_fit[x,31] <- 1,Results_fit[x,31] <- 2)}
    if(Results_fit[x,1]=="Guppy"){ifelse(Results_fit[x,3]=="f",Results_fit[x,31] <- 1,Results_fit[x,31] <- 2)}
    Results_fit[x,32] <- paste(Results_fit[x,1],Results_fit[x,3],sep="")
    Results_fit[x,33] <- Data_ID[which.max(Data_ID$Measurement),7]
    Results_fit[x,34] <- mean(na.omit(Data_ID$Mass))
    Results_fit[x,35] <- log10(Results_fit[x,34])*slope_SMR$coefficients[2]+slope_SMR$coefficients[1]
    Results_fit[x,36] <- mean(na.omit(Data_ID$Temp))
      }
}
# small adjustments to DF
Results_fit <- Results_fit[!is.na(Results_fit$Start_weight),]
Results_fit$Plot_treatment <- as.factor(Results_fit$Plot_treatment)



## restructuring Results, varibles in one culum
Results_fit_p <- data.frame(matrix(nrow=0,ncol=20))
names(Results_fit_p) <- c("Species","FishID","Treatment", "Start_weight","End_weight", "Fold_change","Number_of_messurements","SGR", "AGR","Parameter",
                        "Plot_treatment","Plot_group", "SMR_end","Mean_w","SMR_mean_w","Temp","Slope_fit","p_vaule_fit","Intersept_fit","Int_at_mean")
Para <- c("SMR","MMR","AS","SGR","AGR","GE","FAS")
for (x in 1:length(Results_fit$Species)) {
  for (y in 1:7) {
    ifelse(x==1&y==1,nr <-1, nr<-nr+1)
    Results_fit_p[nr,1] <- Results_fit[x,1]
    Results_fit_p[nr,2] <- Results_fit[x,2]
    Results_fit_p[nr,3] <- Results_fit[x,3]
    Results_fit_p[nr,4] <- Results_fit[x,4]
    Results_fit_p[nr,5] <- Results_fit[x,5]
    Results_fit_p[nr,6] <- Results_fit[x,6]
    Results_fit_p[nr,7] <- Results_fit[x,7]
    Results_fit_p[nr,8] <- Results_fit[x,8]
    Results_fit_p[nr,9] <- Results_fit[x,9]
    Results_fit_p[nr,10] <- Para[y]
    Results_fit_p[nr,11] <- Results_fit[x,31]
    Results_fit_p[nr,12] <- Results_fit[x,32]
    Results_fit_p[nr,13] <- Results_fit[x,33]
    Results_fit_p[nr,14] <- Results_fit[x,34]
    Results_fit_p[nr,15] <- Results_fit[x,35]
    Results_fit_p[nr,16] <- Results_fit[x,36]
    Results_fit_p[nr,17] <- Results_fit[x,10+y*3-3]
    Results_fit_p[nr,18] <- Results_fit[x,11+y*3-3]
    Results_fit_p[nr,19] <- Results_fit[x,12+y*3-3]
    Results_fit_p[nr,20] <- Results_fit_p[nr,17]*LogOverAllMeanMass+Results_fit_p[nr,19]
  }
}


### Finding the scaling exponents for each species for each paremeter
##Zebrafish##
Zebrafish_scal <- data.frame("Species"=rep("Zebrafish",4),"Treatment"=c("Mean","High_fed","Low_fed","Moderate_fed") ,
                             "bSMR_onto"=c(summary(mcmc_SMR_Zebrafish_useM)$solutions[2,1],summary(mcmc_SMR_Zebrafish_use)$solutions[2,1],(summary(mcmc_SMR_Zebrafish_use)$solutions[2,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[7,1]),(summary(mcmc_SMR_Zebrafish_use)$solutions[2,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[8,1])),
                             "bSMR_onto_low"=c(summary(mcmc_SMR_Zebrafish_useM)$solutions[2,2],summary(mcmc_SMR_Zebrafish_use)$solutions[2,2],summary(mcmc_SMR_Zebrafish_use)$solutions[2,2]+summary(mcmc_SMR_Zebrafish_use)$solutions[7,2],summary(mcmc_SMR_Zebrafish_use)$solutions[2,2]+summary(mcmc_SMR_Zebrafish_use)$solutions[8,2]),
                             "bSMR_onto_high"=c(summary(mcmc_SMR_Zebrafish_useM)$solutions[2,3],summary(mcmc_SMR_Zebrafish_use)$solutions[2,3],summary(mcmc_SMR_Zebrafish_use)$solutions[2,3]+summary(mcmc_SMR_Zebrafish_use)$solutions[7,3],summary(mcmc_SMR_Zebrafish_use)$solutions[2,3]+summary(mcmc_SMR_Zebrafish_use)$solutions[8,3]),
                             "bMMR_onto"=c(NA,NA,NA,NA),"bMMR_onto_low"=c(NA,NA,NA,NA),"bMMR_onto_high"=c(NA,NA,NA,NA) ,
                             "bAS_onto"=c(NA,NA,NA,NA),"bAS_onto_low"=c(NA,NA,NA,NA),"bAS_onto_high"=c(NA,NA,NA,NA),
                             "bFAS_onto"=c(rep(NA,4)),
                             "bFAS_onto_low"=rep(NA,4),
                             "bFAS_onto_high"=rep(NA,4),
                             "bAGR_onto"=c(summary(mcmc_AGR_Zebrafish_useM)$solutions[2,1],summary(mcmc_AGR_Zebrafish_use)$solutions[2,1],(summary(mcmc_AGR_Zebrafish_use)$solutions[2,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[6,1]),(summary(mcmc_AGR_Zebrafish_use)$solutions[2,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[7,1])),
                             "bAGR_onto_low"=c(summary(mcmc_AGR_Zebrafish_useM)$solutions[2,2],summary(mcmc_AGR_Zebrafish_use)$solutions[2,2],summary(mcmc_AGR_Zebrafish_use)$solutions[2,2]+summary(mcmc_AGR_Zebrafish_use)$solutions[6,2],summary(mcmc_AGR_Zebrafish_use)$solutions[2,2]+summary(mcmc_AGR_Zebrafish_use)$solutions[7,2]),
                             "bAGR_onto_high"=c(summary(mcmc_AGR_Zebrafish_useM)$solutions[2,3],summary(mcmc_AGR_Zebrafish_use)$solutions[2,3],summary(mcmc_AGR_Zebrafish_use)$solutions[2,3]+summary(mcmc_AGR_Zebrafish_use)$solutions[6,3],summary(mcmc_AGR_Zebrafish_use)$solutions[2,3]+summary(mcmc_AGR_Zebrafish_use)$solutions[7,3]),
                             "bSGR_onto"=c(summary(mcmc_SGR_Zebrafish_useM)$solutions[2,1],summary(mcmc_SGR_Zebrafish_use)$solutions[2,1],(summary(mcmc_SGR_Zebrafish_use)$solutions[2,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[6,1]),(summary(mcmc_SGR_Zebrafish_use)$solutions[2,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[7,1])),
                             "bSGR_onto_low"=c(summary(mcmc_SGR_Zebrafish_useM)$solutions[2,2],summary(mcmc_SGR_Zebrafish_use)$solutions[2,2],summary(mcmc_SGR_Zebrafish_use)$solutions[2,2]+summary(mcmc_SGR_Zebrafish_use)$solutions[6,2],summary(mcmc_SGR_Zebrafish_use)$solutions[2,2]+summary(mcmc_SGR_Zebrafish_use)$solutions[7,2]),
                             "bSGR_onto_high"=c(summary(mcmc_SGR_Zebrafish_useM)$solutions[2,3],summary(mcmc_SGR_Zebrafish_use)$solutions[2,3],summary(mcmc_SGR_Zebrafish_use)$solutions[2,3]+summary(mcmc_SGR_Zebrafish_use)$solutions[6,3],summary(mcmc_SGR_Zebrafish_use)$solutions[2,3]+summary(mcmc_SGR_Zebrafish_use)$solutions[7,3]) ,
                             "bGE_onto"=c(summary(mcmc_GE_Zebrafish_useM)$solutions[2,1],summary(mcmc_GE_Zebrafish_use)$solutions[2,1],(summary(mcmc_GE_Zebrafish_use)$solutions[2,1]+summary(mcmc_GE_Zebrafish_use)$solutions[6,1]),(summary(mcmc_GE_Zebrafish_use)$solutions[2,1]+summary(mcmc_GE_Zebrafish_use)$solutions[7,1])),
                             "bGE_onto_low"=c(summary(mcmc_GE_Zebrafish_useM)$solutions[2,2],summary(mcmc_GE_Zebrafish_use)$solutions[2,2],summary(mcmc_GE_Zebrafish_use)$solutions[2,2]+summary(mcmc_GE_Zebrafish_use)$solutions[6,2],summary(mcmc_GE_Zebrafish_use)$solutions[2,2]+summary(mcmc_GE_Zebrafish_use)$solutions[7,2]),
                             "bGE_onto_high"=c(summary(mcmc_GE_Zebrafish_useM)$solutions[2,3],summary(mcmc_GE_Zebrafish_use)$solutions[2,3],summary(mcmc_GE_Zebrafish_use)$solutions[2,3]+summary(mcmc_GE_Zebrafish_use)$solutions[6,3],summary(mcmc_GE_Zebrafish_use)$solutions[2,3]+summary(mcmc_GE_Zebrafish_use)$solutions[7,3]),
                            
                             
                                               "bSMR_static"=c(summary(mcmc_SMR_Zebrafish_useM)$solutions[3,1],summary(mcmc_SMR_Zebrafish_use)$solutions[5,1],(summary(mcmc_SMR_Zebrafish_use)$solutions[5,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[9,1]),(summary(mcmc_SMR_Zebrafish_use)$solutions[5,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[10,1])),
                                               "bSMR_static_low"=c(summary(mcmc_SMR_Zebrafish_useM)$solutions[3,2],summary(mcmc_SMR_Zebrafish_use)$solutions[5,2],summary(mcmc_SMR_Zebrafish_use)$solutions[5,2]+summary(mcmc_SMR_Zebrafish_use)$solutions[9,2],summary(mcmc_SMR_Zebrafish_use)$solutions[5,2]+summary(mcmc_SMR_Zebrafish_use)$solutions[10,2]),
                                               "bSMR_static_high"=c(summary(mcmc_SMR_Zebrafish_useM)$solutions[3,3],summary(mcmc_SMR_Zebrafish_use)$solutions[5,3],summary(mcmc_SMR_Zebrafish_use)$solutions[5,3]+summary(mcmc_SMR_Zebrafish_use)$solutions[9,3],summary(mcmc_SMR_Zebrafish_use)$solutions[5,3]+summary(mcmc_SMR_Zebrafish_use)$solutions[10,3]),
                                               "bMMR_static"=c(NA,NA,NA,NA),"bMMR_static_low"=c(NA,NA,NA,NA),"bMMR_static_high"=c(NA,NA,NA,NA) ,
                                               "bAS_static"=c(NA,NA,NA,NA),"bAS_static_low"=c(NA,NA,NA,NA),"bAS_static_high"=c(NA,NA,NA,NA),
                             "bFAS_static"=c(NA,NA,NA,NA),
                             "bFAS_static_low"=c(NA,NA,NA,NA),
                             "bFAS_static_high"=c(NA,NA,NA,NA),
                                               "bAGR_static"=c(summary(mcmc_AGR_Zebrafish_useM)$solutions[3,1],summary(mcmc_AGR_Zebrafish_use)$solutions[5,1],(summary(mcmc_AGR_Zebrafish_use)$solutions[5,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[8,1]),(summary(mcmc_AGR_Zebrafish_use)$solutions[5,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[9,1])),
                                               "bAGR_static_low"=c(summary(mcmc_AGR_Zebrafish_useM)$solutions[3,2],summary(mcmc_AGR_Zebrafish_use)$solutions[5,2],summary(mcmc_AGR_Zebrafish_use)$solutions[5,2]+summary(mcmc_AGR_Zebrafish_use)$solutions[8,2],summary(mcmc_AGR_Zebrafish_use)$solutions[5,2]+summary(mcmc_AGR_Zebrafish_use)$solutions[9,2]),
                                               "bAGR_static_high"=c(summary(mcmc_AGR_Zebrafish_useM)$solutions[3,3],summary(mcmc_AGR_Zebrafish_use)$solutions[5,3],summary(mcmc_AGR_Zebrafish_use)$solutions[5,3]+summary(mcmc_AGR_Zebrafish_use)$solutions[8,3],summary(mcmc_AGR_Zebrafish_use)$solutions[5,3]+summary(mcmc_AGR_Zebrafish_use)$solutions[9,3]),
                                               "bSGR_static"=c(summary(mcmc_SGR_Zebrafish_useM)$solutions[3,1],summary(mcmc_SGR_Zebrafish_use)$solutions[5,1],(summary(mcmc_SGR_Zebrafish_use)$solutions[5,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[8,1]),(summary(mcmc_SGR_Zebrafish_use)$solutions[5,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[9,1])),
                                               "bSGR_static_low"=c(summary(mcmc_SGR_Zebrafish_useM)$solutions[3,2],summary(mcmc_SGR_Zebrafish_use)$solutions[5,2],summary(mcmc_SGR_Zebrafish_use)$solutions[5,2]+summary(mcmc_SGR_Zebrafish_use)$solutions[8,2],summary(mcmc_SGR_Zebrafish_use)$solutions[5,2]+summary(mcmc_SGR_Zebrafish_use)$solutions[9,2]),
                                               "bSGR_static_high"=c(summary(mcmc_SGR_Zebrafish_useM)$solutions[3,3],summary(mcmc_SGR_Zebrafish_use)$solutions[5,3],summary(mcmc_SGR_Zebrafish_use)$solutions[5,3]+summary(mcmc_SGR_Zebrafish_use)$solutions[8,3],summary(mcmc_SGR_Zebrafish_use)$solutions[5,3]+summary(mcmc_SGR_Zebrafish_use)$solutions[9,3]) ,
                                               "bGE_static"=c(summary(mcmc_GE_Zebrafish_useM)$solutions[3,1],summary(mcmc_GE_Zebrafish_use)$solutions[5,1],(summary(mcmc_GE_Zebrafish_use)$solutions[5,1]+summary(mcmc_GE_Zebrafish_use)$solutions[8,1]),(summary(mcmc_GE_Zebrafish_use)$solutions[5,1]+summary(mcmc_GE_Zebrafish_use)$solutions[9,1])),
                                               "bGE_static_low"=c(summary(mcmc_GE_Zebrafish_useM)$solutions[3,2],summary(mcmc_GE_Zebrafish_use)$solutions[5,2],summary(mcmc_GE_Zebrafish_use)$solutions[5,2]+summary(mcmc_GE_Zebrafish_use)$solutions[8,2],summary(mcmc_GE_Zebrafish_use)$solutions[5,2]+summary(mcmc_GE_Zebrafish_use)$solutions[9,2]),
                                               "bGE_static_high"=c(summary(mcmc_GE_Zebrafish_useM)$solutions[3,3],summary(mcmc_GE_Zebrafish_use)$solutions[5,3],summary(mcmc_GE_Zebrafish_use)$solutions[5,3]+summary(mcmc_GE_Zebrafish_use)$solutions[8,3],summary(mcmc_GE_Zebrafish_use)$solutions[5,3]+summary(mcmc_GE_Zebrafish_use)$solutions[9,3])
                              
                             
                             )
                             

##Rainbow trout##
Rainbow_trout_scal <- data.frame("Species"=rep("Rainbow_trout",4),"Treatment"=c("Mean","10Cegg","14Cegg","14Cyolk") ,
                             "bSMR_onto"=c(summary(mcmc_SMR_Rainbow_trout_useM)$solutions[2,1],summary(mcmc_SMR_Rainbow_trout_use)$solutions[2,1],(summary(mcmc_SMR_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[7,1]),(summary(mcmc_SMR_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[8,1])),
                             "bSMR_onto_low"=c(summary(mcmc_SMR_Rainbow_trout_useM)$solutions[2,2],summary(mcmc_SMR_Rainbow_trout_use)$solutions[2,2],summary(mcmc_SMR_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[7,2],summary(mcmc_SMR_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[8,2]),
                             "bSMR_onto_high"=c(summary(mcmc_SMR_Rainbow_trout_useM)$solutions[2,3],summary(mcmc_SMR_Rainbow_trout_use)$solutions[2,3],summary(mcmc_SMR_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[7,3],summary(mcmc_SMR_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[8,3]),
                             "bMMR_onto"=c(summary(mcmc_MMR_Rainbow_trout_useM)$solutions[2,1],summary(mcmc_MMR_Rainbow_trout_use)$solutions[2,1],(summary(mcmc_MMR_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[7,1]),(summary(mcmc_MMR_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[8,1])),
                             "bMMR_onto_low"=c(summary(mcmc_MMR_Rainbow_trout_useM)$solutions[2,2],summary(mcmc_MMR_Rainbow_trout_use)$solutions[2,2],summary(mcmc_MMR_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[7,2],summary(mcmc_MMR_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[8,2]),
                             "bMMR_onto_high"=c(summary(mcmc_MMR_Rainbow_trout_useM)$solutions[2,3],summary(mcmc_MMR_Rainbow_trout_use)$solutions[2,3],summary(mcmc_MMR_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[7,3],summary(mcmc_MMR_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[8,3]),
                             "bAS_onto"=c(summary(mcmc_AS_Rainbow_trout_useM)$solutions[2,1],summary(mcmc_AS_Rainbow_trout_use)$solutions[2,1],(summary(mcmc_AS_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[7,1]),(summary(mcmc_AS_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[8,1])),
                             "bAS_onto_low"=c(summary(mcmc_AS_Rainbow_trout_useM)$solutions[2,2],summary(mcmc_AS_Rainbow_trout_use)$solutions[2,2],summary(mcmc_AS_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_AS_Rainbow_trout_use)$solutions[7,2],summary(mcmc_AS_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_AS_Rainbow_trout_use)$solutions[8,2]),
                             "bAS_onto_high"=c(summary(mcmc_AS_Rainbow_trout_useM)$solutions[2,3],summary(mcmc_AS_Rainbow_trout_use)$solutions[2,3],summary(mcmc_AS_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_AS_Rainbow_trout_use)$solutions[7,3],summary(mcmc_AS_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_AS_Rainbow_trout_use)$solutions[8,3]),
                             "bFAS_onto"=c(summary(mcmc_FAS_Rainbow_trout_useM)$solutions[2,1],summary(mcmc_FAS_Rainbow_trout_use)$solutions[2,1],(summary(mcmc_FAS_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[7,1]),(summary(mcmc_FAS_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[8,1])),
                             "bFAS_onto_low"=c(summary(mcmc_FAS_Rainbow_trout_useM)$solutions[2,2],summary(mcmc_FAS_Rainbow_trout_use)$solutions[2,2],summary(mcmc_FAS_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[7,2],summary(mcmc_FAS_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[8,2]),
                             "bFAS_onto_high"=c(summary(mcmc_FAS_Rainbow_trout_useM)$solutions[2,3],summary(mcmc_FAS_Rainbow_trout_use)$solutions[2,3],summary(mcmc_FAS_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[7,3],summary(mcmc_FAS_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[8,3]),
                             
                             "bAGR_onto"=c(summary(mcmc_AGR_Rainbow_trout_useM)$solutions[2,1],summary(mcmc_AGR_Rainbow_trout_use)$solutions[2,1],(summary(mcmc_AGR_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[6,1]),(summary(mcmc_AGR_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[7,1])),
                             "bAGR_onto_low"=c(summary(mcmc_AGR_Rainbow_trout_useM)$solutions[2,2],summary(mcmc_AGR_Rainbow_trout_use)$solutions[2,2],summary(mcmc_AGR_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[6,2],summary(mcmc_AGR_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[7,2]),
                             "bAGR_onto_high"=c(summary(mcmc_AGR_Rainbow_trout_useM)$solutions[2,3],summary(mcmc_AGR_Rainbow_trout_use)$solutions[2,3],summary(mcmc_AGR_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[6,3],summary(mcmc_AGR_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[7,3]),
                             "bSGR_onto"=c(summary(mcmc_SGR_Rainbow_trout_useM)$solutions[2,1],summary(mcmc_SGR_Rainbow_trout_use)$solutions[2,1],(summary(mcmc_SGR_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[6,1]),(summary(mcmc_SGR_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[7,1])),
                             "bSGR_onto_low"=c(summary(mcmc_SGR_Rainbow_trout_useM)$solutions[2,2],summary(mcmc_SGR_Rainbow_trout_use)$solutions[2,2],summary(mcmc_SGR_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[6,2],summary(mcmc_SGR_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[7,2]),
                             "bSGR_onto_high"=c(summary(mcmc_SGR_Rainbow_trout_useM)$solutions[2,3],summary(mcmc_SGR_Rainbow_trout_use)$solutions[2,3],summary(mcmc_SGR_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[6,3],summary(mcmc_SGR_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[7,3]) ,
                             "bGE_onto"=c(summary(mcmc_GE_Rainbow_trout_useM)$solutions[2,1],summary(mcmc_GE_Rainbow_trout_use)$solutions[2,1],(summary(mcmc_GE_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[6,1]),(summary(mcmc_GE_Rainbow_trout_use)$solutions[2,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[7,1])),
                             "bGE_onto_low"=c(summary(mcmc_GE_Rainbow_trout_useM)$solutions[2,2],summary(mcmc_GE_Rainbow_trout_use)$solutions[2,2],summary(mcmc_GE_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_GE_Rainbow_trout_use)$solutions[6,2],summary(mcmc_GE_Rainbow_trout_use)$solutions[2,2]+summary(mcmc_GE_Rainbow_trout_use)$solutions[7,2]),
                             "bGE_onto_high"=c(summary(mcmc_GE_Rainbow_trout_useM)$solutions[2,3],summary(mcmc_GE_Rainbow_trout_use)$solutions[2,3],summary(mcmc_GE_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_GE_Rainbow_trout_use)$solutions[6,3],summary(mcmc_GE_Rainbow_trout_use)$solutions[2,3]+summary(mcmc_GE_Rainbow_trout_use)$solutions[7,3]),
                                               "bSMR_static"=c(summary(mcmc_SMR_Rainbow_trout_useM)$solutions[3,1],summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,1],(summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[9,1]),(summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[10,1])),
                                               "bSMR_static_low"=c(summary(mcmc_SMR_Rainbow_trout_useM)$solutions[3,2],summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,2],summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[9,2],summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[10,2]),
                                               "bSMR_static_high"=c(summary(mcmc_SMR_Rainbow_trout_useM)$solutions[3,3],summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,3],summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[9,3],summary(mcmc_SMR_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[10,3]),
                                               "bMMR_static"=c(summary(mcmc_MMR_Rainbow_trout_useM)$solutions[3,1],summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,1],(summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[9,1]),(summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[10,1])),
                                               "bMMR_static_low"=c(summary(mcmc_MMR_Rainbow_trout_useM)$solutions[3,2],summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,2],summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[9,2],summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[10,2]),
                                               "bMMR_static_high"=c(summary(mcmc_MMR_Rainbow_trout_useM)$solutions[3,3],summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,3],summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[9,3],summary(mcmc_MMR_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[10,3]),
                                               "bAS_static"=c(summary(mcmc_AS_Rainbow_trout_useM)$solutions[3,1],summary(mcmc_AS_Rainbow_trout_use)$solutions[5,1],(summary(mcmc_AS_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[9,1]),(summary(mcmc_AS_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[10,1])),
                                               "bAS_static_low"=c(summary(mcmc_AS_Rainbow_trout_useM)$solutions[3,2],summary(mcmc_AS_Rainbow_trout_use)$solutions[5,2],summary(mcmc_AS_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_AS_Rainbow_trout_use)$solutions[9,2],summary(mcmc_AS_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_AS_Rainbow_trout_use)$solutions[10,2]),
                                               "bAS_static_high"=c(summary(mcmc_AS_Rainbow_trout_useM)$solutions[3,3],summary(mcmc_AS_Rainbow_trout_use)$solutions[5,3],summary(mcmc_AS_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_AS_Rainbow_trout_use)$solutions[9,3],summary(mcmc_AS_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_AS_Rainbow_trout_use)$solutions[10,3]),
                             "bFAS_static"=c(summary(mcmc_FAS_Rainbow_trout_useM)$solutions[3,1],summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,1],(summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[9,1]),(summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[10,1])),
                             "bFAS_static_low"=c(summary(mcmc_FAS_Rainbow_trout_useM)$solutions[3,2],summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,2],summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[9,2],summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[10,2]),
                             "bFAS_static_high"=c(summary(mcmc_FAS_Rainbow_trout_useM)$solutions[3,3],summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,3],summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[9,3],summary(mcmc_FAS_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[10,3]),
                             
                                               "bAGR_static"=c(summary(mcmc_AGR_Rainbow_trout_useM)$solutions[3,1],summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,1],(summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[8,1]),(summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[9,1])),
                                               "bAGR_static_low"=c(summary(mcmc_AGR_Rainbow_trout_useM)$solutions[3,2],summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,2],summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[8,2],summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[9,2]),
                                               "bAGR_static_high"=c(summary(mcmc_AGR_Rainbow_trout_useM)$solutions[3,3],summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,3],summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[8,3],summary(mcmc_AGR_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[9,3]),
                                               "bSGR_static"=c(summary(mcmc_SGR_Rainbow_trout_useM)$solutions[3,1],summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,1],(summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[8,1]),(summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[9,1])),
                                               "bSGR_static_low"=c(summary(mcmc_SGR_Rainbow_trout_useM)$solutions[3,2],summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,2],summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[8,2],summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[9,2]),
                                               "bSGR_static_high"=c(summary(mcmc_SGR_Rainbow_trout_useM)$solutions[3,3],summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,3],summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[8,3],summary(mcmc_SGR_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[9,3]) ,
                                               "bGE_static"=c(summary(mcmc_GE_Rainbow_trout_useM)$solutions[3,1],summary(mcmc_GE_Rainbow_trout_use)$solutions[5,1],(summary(mcmc_GE_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[8,1]),(summary(mcmc_GE_Rainbow_trout_use)$solutions[5,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[9,1])),
                                               "bGE_static_low"=c(summary(mcmc_GE_Rainbow_trout_useM)$solutions[3,2],summary(mcmc_GE_Rainbow_trout_use)$solutions[5,2],summary(mcmc_GE_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_GE_Rainbow_trout_use)$solutions[8,2],summary(mcmc_GE_Rainbow_trout_use)$solutions[5,2]+summary(mcmc_GE_Rainbow_trout_use)$solutions[9,2]),
                                               "bGE_static_high"=c(summary(mcmc_GE_Rainbow_trout_useM)$solutions[3,3],summary(mcmc_GE_Rainbow_trout_use)$solutions[5,3],summary(mcmc_GE_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_GE_Rainbow_trout_use)$solutions[8,3],summary(mcmc_GE_Rainbow_trout_use)$solutions[5,3]+summary(mcmc_GE_Rainbow_trout_use)$solutions[9,3])
                             )
Chromis_scal <- data.frame("Species"=rep("Chromis",3),"Treatment"=c("Mean","High_fed","Low_fed") ,
                                 "bSMR_onto"=c(summary(mcmc_SMR_Chromis_useM)$solutions[2,1],summary(mcmc_SMR_RR_Chromis)$solutions[2,1],(summary(mcmc_SMR_RR_Chromis)$solutions[2,1]+summary(mcmc_SMR_RR_Chromis)$solutions[5,1])),
                                 "bSMR_onto_low"=c(summary(mcmc_SMR_Chromis_useM)$solutions[2,2],summary(mcmc_SMR_RR_Chromis)$solutions[2,2],summary(mcmc_SMR_RR_Chromis)$solutions[2,2]+summary(mcmc_SMR_RR_Chromis)$solutions[5,2]),
                                 "bSMR_onto_high"=c(summary(mcmc_SMR_Chromis_useM)$solutions[2,3],summary(mcmc_SMR_RR_Chromis)$solutions[2,3],summary(mcmc_SMR_RR_Chromis)$solutions[2,3]+summary(mcmc_SMR_RR_Chromis)$solutions[5,3]),
                                 "bMMR_onto"=c(summary(mcmc_MMR_Chromis_useM)$solutions[2,1],summary(mcmc_MMR_RR_Chromis)$solutions[2,1],(summary(mcmc_MMR_RR_Chromis)$solutions[2,1]+summary(mcmc_MMR_RR_Chromis)$solutions[5,1])),
                                 "bMMR_onto_low"=c(summary(mcmc_MMR_Chromis_useM)$solutions[2,2],summary(mcmc_MMR_RR_Chromis)$solutions[2,2],summary(mcmc_MMR_RR_Chromis)$solutions[2,2]+summary(mcmc_MMR_RR_Chromis)$solutions[5,2]),
                                 "bMMR_onto_high"=c(summary(mcmc_MMR_Chromis_useM)$solutions[2,3],summary(mcmc_MMR_RR_Chromis)$solutions[2,3],summary(mcmc_MMR_RR_Chromis)$solutions[2,3]+summary(mcmc_MMR_RR_Chromis)$solutions[5,3]),
                                 "bAS_onto"=c(summary(mcmc_AS_Chromis_useM)$solutions[2,1],summary(mcmc_AS_RR_Chromis)$solutions[2,1],(summary(mcmc_AS_RR_Chromis)$solutions[2,1]+summary(mcmc_AS_RR_Chromis)$solutions[5,1])),
                                 "bAS_onto_low"=c(summary(mcmc_AS_Chromis_useM)$solutions[2,2],summary(mcmc_AS_RR_Chromis)$solutions[2,2],summary(mcmc_AS_RR_Chromis)$solutions[2,2]+summary(mcmc_AS_RR_Chromis)$solutions[5,2]),
                                 "bAS_onto_high"=c(summary(mcmc_AS_Chromis_useM)$solutions[2,3],summary(mcmc_AS_RR_Chromis)$solutions[2,3],summary(mcmc_AS_RR_Chromis)$solutions[2,3]+summary(mcmc_AS_RR_Chromis)$solutions[5,3]),
                           "bFAS_onto"=c(summary(mcmc_FAS_Chromis_useM)$solutions[2,1],summary(mcmc_FAS_RR_Chromis)$solutions[2,1],(summary(mcmc_FAS_RR_Chromis)$solutions[2,1]+summary(mcmc_FAS_RR_Chromis)$solutions[5,1])),
                           "bFAS_onto_low"=c(summary(mcmc_FAS_Chromis_useM)$solutions[2,2],summary(mcmc_FAS_RR_Chromis)$solutions[2,2],summary(mcmc_FAS_RR_Chromis)$solutions[2,2]+summary(mcmc_FAS_RR_Chromis)$solutions[5,2]),
                           "bFAS_onto_high"=c(summary(mcmc_FAS_Chromis_useM)$solutions[2,3],summary(mcmc_FAS_RR_Chromis)$solutions[2,3],summary(mcmc_FAS_RR_Chromis)$solutions[2,3]+summary(mcmc_FAS_RR_Chromis)$solutions[5,3]),
                           
                                 "bAGR_onto"=c(summary(mcmc_AGR_Chromis_useM)$solutions[2,1],summary(mcmc_AGR_RR_Chromis)$solutions[2,1],(summary(mcmc_AGR_RR_Chromis)$solutions[2,1]+summary(mcmc_AGR_RR_Chromis)$solutions[5,1])),
                                 "bAGR_onto_low"=c(summary(mcmc_AGR_Chromis_useM)$solutions[2,2],summary(mcmc_AGR_RR_Chromis)$solutions[2,2],summary(mcmc_AGR_RR_Chromis)$solutions[2,2]+summary(mcmc_AGR_RR_Chromis)$solutions[5,2]),
                                 "bAGR_onto_high"=c(summary(mcmc_AGR_Chromis_useM)$solutions[2,3],summary(mcmc_AGR_RR_Chromis)$solutions[2,3],summary(mcmc_AGR_RR_Chromis)$solutions[2,3]+summary(mcmc_AGR_RR_Chromis)$solutions[5,3]),
                                 "bSGR_onto"=c(summary(mcmc_SGR_Chromis_useM)$solutions[2,1],summary(mcmc_SGR_RR_Chromis)$solutions[2,1],(summary(mcmc_SGR_RR_Chromis)$solutions[2,1]+summary(mcmc_SGR_RR_Chromis)$solutions[5,1])),
                                 "bSGR_onto_low"=c(summary(mcmc_SGR_Chromis_useM)$solutions[2,2],summary(mcmc_SGR_RR_Chromis)$solutions[2,2],summary(mcmc_SGR_RR_Chromis)$solutions[2,2]+summary(mcmc_SGR_RR_Chromis)$solutions[5,2]),
                                 "bSGR_onto_high"=c(summary(mcmc_SGR_Chromis_useM)$solutions[2,3],summary(mcmc_SGR_RR_Chromis)$solutions[2,3],summary(mcmc_SGR_RR_Chromis)$solutions[2,3]+summary(mcmc_SGR_RR_Chromis)$solutions[5,3]) ,
                                 "bGE_onto"=c(summary(mcmc_GE_Chromis_useM)$solutions[2,1],summary(mcmc_GE_RR_Chromis)$solutions[2,1],(summary(mcmc_GE_RR_Chromis)$solutions[2,1]+summary(mcmc_GE_RR_Chromis)$solutions[5,1])),
                                 "bGE_onto_low"=c(summary(mcmc_GE_Chromis_useM)$solutions[2,2],summary(mcmc_GE_RR_Chromis)$solutions[2,2],summary(mcmc_GE_RR_Chromis)$solutions[2,2]+summary(mcmc_GE_RR_Chromis)$solutions[5,2]),
                                 "bGE_onto_high"=c(summary(mcmc_GE_Chromis_useM)$solutions[2,3],summary(mcmc_GE_RR_Chromis)$solutions[2,3],summary(mcmc_GE_RR_Chromis)$solutions[2,3]+summary(mcmc_GE_RR_Chromis)$solutions[5,3]),
                                                   "bSMR_static"=c(summary(mcmc_SMR_Chromis_useM)$solutions[3,1],summary(mcmc_SMR_RR_Chromis)$solutions[4,1],(summary(mcmc_SMR_RR_Chromis)$solutions[4,1]+summary(mcmc_SMR_RR_Chromis)$solutions[6,1])),
                                                   "bSMR_static_low"=c(summary(mcmc_SMR_Chromis_useM)$solutions[3,2],summary(mcmc_SMR_RR_Chromis)$solutions[4,2],summary(mcmc_SMR_RR_Chromis)$solutions[4,2]+summary(mcmc_SMR_RR_Chromis)$solutions[6,2]),
                                                   "bSMR_static_high"=c(summary(mcmc_SMR_Chromis_useM)$solutions[3,3],summary(mcmc_SMR_RR_Chromis)$solutions[4,3],summary(mcmc_SMR_RR_Chromis)$solutions[4,3]+summary(mcmc_SMR_RR_Chromis)$solutions[6,3]),
                                                   "bMMR_static"=c(summary(mcmc_MMR_Chromis_useM)$solutions[3,1],summary(mcmc_MMR_RR_Chromis)$solutions[4,1],(summary(mcmc_MMR_RR_Chromis)$solutions[4,1]+summary(mcmc_MMR_RR_Chromis)$solutions[6,1])),
                                                   "bMMR_static_low"=c(summary(mcmc_MMR_Chromis_useM)$solutions[3,2],summary(mcmc_MMR_RR_Chromis)$solutions[4,2],summary(mcmc_MMR_RR_Chromis)$solutions[4,2]+summary(mcmc_MMR_RR_Chromis)$solutions[6,2]),
                                                   "bMMR_static_high"=c(summary(mcmc_MMR_Chromis_useM)$solutions[3,3],summary(mcmc_MMR_RR_Chromis)$solutions[4,3],summary(mcmc_MMR_RR_Chromis)$solutions[4,3]+summary(mcmc_MMR_RR_Chromis)$solutions[6,3]),
                                                   "bAS_static"=c(summary(mcmc_AS_Chromis_useM)$solutions[3,1],summary(mcmc_AS_RR_Chromis)$solutions[4,1],(summary(mcmc_AS_RR_Chromis)$solutions[4,1]+summary(mcmc_AS_RR_Chromis)$solutions[6,1])),
                                                   "bAS_static_low"=c(summary(mcmc_AS_Chromis_useM)$solutions[3,2],summary(mcmc_AS_RR_Chromis)$solutions[4,2],summary(mcmc_AS_RR_Chromis)$solutions[4,2]+summary(mcmc_AS_RR_Chromis)$solutions[6,2]),
                                                   "bAS_static_high"=c(summary(mcmc_AS_Chromis_useM)$solutions[3,3],summary(mcmc_AS_RR_Chromis)$solutions[4,3],summary(mcmc_AS_RR_Chromis)$solutions[4,3]+summary(mcmc_AS_RR_Chromis)$solutions[6,3]),
                           "bFAS_static"=c(summary(mcmc_FAS_Chromis_useM)$solutions[3,1],summary(mcmc_FAS_RR_Chromis)$solutions[4,1],(summary(mcmc_FAS_RR_Chromis)$solutions[4,1]+summary(mcmc_FAS_RR_Chromis)$solutions[6,1])),
                           "bFAS_static_low"=c(summary(mcmc_FAS_Chromis_useM)$solutions[3,2],summary(mcmc_FAS_RR_Chromis)$solutions[4,2],summary(mcmc_FAS_RR_Chromis)$solutions[4,2]+summary(mcmc_FAS_RR_Chromis)$solutions[6,2]),
                           "bFAS_static_high"=c(summary(mcmc_FAS_Chromis_useM)$solutions[3,3],summary(mcmc_FAS_RR_Chromis)$solutions[4,3],summary(mcmc_FAS_RR_Chromis)$solutions[4,3]+summary(mcmc_FAS_RR_Chromis)$solutions[6,3]),
                           
                                                   "bAGR_static"=c(summary(mcmc_AGR_Chromis_useM)$solutions[3,1],summary(mcmc_AGR_RR_Chromis)$solutions[4,1],(summary(mcmc_AGR_RR_Chromis)$solutions[4,1]+summary(mcmc_AGR_RR_Chromis)$solutions[6,1])),
                                                   "bAGR_static_low"=c(summary(mcmc_AGR_Chromis_useM)$solutions[3,2],summary(mcmc_AGR_RR_Chromis)$solutions[4,2],summary(mcmc_AGR_RR_Chromis)$solutions[4,2]+summary(mcmc_AGR_RR_Chromis)$solutions[6,2]),
                                                   "bAGR_static_high"=c(summary(mcmc_AGR_Chromis_useM)$solutions[3,3],summary(mcmc_AGR_RR_Chromis)$solutions[4,3],summary(mcmc_AGR_RR_Chromis)$solutions[4,3]+summary(mcmc_AGR_RR_Chromis)$solutions[6,3]),
                                                   "bSGR_static"=c(summary(mcmc_SGR_Chromis_useM)$solutions[3,1],summary(mcmc_SGR_RR_Chromis)$solutions[4,1],(summary(mcmc_SGR_RR_Chromis)$solutions[4,1]+summary(mcmc_SGR_RR_Chromis)$solutions[6,1])),
                                                   "bSGR_static_low"=c(summary(mcmc_SGR_Chromis_useM)$solutions[3,2],summary(mcmc_SGR_RR_Chromis)$solutions[4,2],summary(mcmc_SGR_RR_Chromis)$solutions[4,2]+summary(mcmc_SGR_RR_Chromis)$solutions[6,2]),
                                                   "bSGR_static_high"=c(summary(mcmc_SGR_Chromis_useM)$solutions[3,3],summary(mcmc_SGR_RR_Chromis)$solutions[4,3],summary(mcmc_SGR_RR_Chromis)$solutions[4,3]+summary(mcmc_SGR_RR_Chromis)$solutions[6,3]) ,
                                                   "bGE_static"=c(summary(mcmc_GE_Chromis_useM)$solutions[3,1],summary(mcmc_GE_RR_Chromis)$solutions[4,1],(summary(mcmc_GE_RR_Chromis)$solutions[4,1]+summary(mcmc_GE_RR_Chromis)$solutions[6,1])),
                                                   "bGE_static_low"=c(summary(mcmc_GE_Chromis_useM)$solutions[3,2],summary(mcmc_GE_RR_Chromis)$solutions[4,2],summary(mcmc_GE_RR_Chromis)$solutions[4,2]+summary(mcmc_GE_RR_Chromis)$solutions[6,2]),
                                                   "bGE_static_high"=c(summary(mcmc_GE_Chromis_useM)$solutions[3,3],summary(mcmc_GE_RR_Chromis)$solutions[4,3],summary(mcmc_GE_RR_Chromis)$solutions[4,3]+summary(mcmc_GE_RR_Chromis)$solutions[6,3])
                                 )

Guppy_scal <- data.frame("Species"=rep("Guppy",3),"Treatment"=c("Mean","Female","Male") ,
                           "bSMR_onto"=c(summary(mcmc_SMR_Guppy_useM)$solutions[2,1],summary(mcmc_SMR_RR_Guppy)$solutions[2,1],(summary(mcmc_SMR_RR_Guppy)$solutions[2,1]+summary(mcmc_SMR_RR_Guppy)$solutions[5,1])),
                           "bSMR_onto_low"=c(summary(mcmc_SMR_Guppy_useM)$solutions[2,2],summary(mcmc_SMR_RR_Guppy)$solutions[2,2],summary(mcmc_SMR_RR_Guppy)$solutions[2,2]+summary(mcmc_SMR_RR_Guppy)$solutions[5,2]),
                           "bSMR_onto_high"=c(summary(mcmc_SMR_Guppy_useM)$solutions[2,3],summary(mcmc_SMR_RR_Guppy)$solutions[2,3],summary(mcmc_SMR_RR_Guppy)$solutions[2,3]+summary(mcmc_SMR_RR_Guppy)$solutions[5,3]),
                           "bMMR_onto"=c(summary(mcmc_MMR_Guppy_useM)$solutions[2,1],summary(mcmc_MMR_RR_Guppy)$solutions[2,1],(summary(mcmc_MMR_RR_Guppy)$solutions[2,1]+summary(mcmc_MMR_RR_Guppy)$solutions[5,1])),
                           "bMMR_onto_low"=c(summary(mcmc_MMR_Guppy_useM)$solutions[2,2],summary(mcmc_MMR_RR_Guppy)$solutions[2,2],summary(mcmc_MMR_RR_Guppy)$solutions[2,2]+summary(mcmc_MMR_RR_Guppy)$solutions[5,2]),
                           "bMMR_onto_high"=c(summary(mcmc_MMR_Guppy_useM)$solutions[2,3],summary(mcmc_MMR_RR_Guppy)$solutions[2,3],summary(mcmc_MMR_RR_Guppy)$solutions[2,3]+summary(mcmc_MMR_RR_Guppy)$solutions[5,3]),
                           "bAS_onto"=c(summary(mcmc_AS_Guppy_useM)$solutions[2,1],summary(mcmc_AS_RR_Guppy)$solutions[2,1],(summary(mcmc_AS_RR_Guppy)$solutions[2,1]+summary(mcmc_AS_RR_Guppy)$solutions[5,1])),
                           "bAS_onto_low"=c(summary(mcmc_AS_Guppy_useM)$solutions[2,2],summary(mcmc_AS_RR_Guppy)$solutions[2,2],summary(mcmc_AS_RR_Guppy)$solutions[2,2]+summary(mcmc_AS_RR_Guppy)$solutions[5,2]),
                           "bAS_onto_high"=c(summary(mcmc_AS_Guppy_useM)$solutions[2,3],summary(mcmc_AS_RR_Guppy)$solutions[2,3],summary(mcmc_AS_RR_Guppy)$solutions[2,3]+summary(mcmc_AS_RR_Guppy)$solutions[5,3]),
                         "bFAS_onto"=c(summary(mcmc_FAS_Guppy_useM)$solutions[2,1],summary(mcmc_FAS_RR_Guppy)$solutions[2,1],(summary(mcmc_FAS_RR_Guppy)$solutions[2,1]+summary(mcmc_FAS_RR_Guppy)$solutions[5,1])),
                         "bFAS_onto_low"=c(summary(mcmc_FAS_Guppy_useM)$solutions[2,2],summary(mcmc_FAS_RR_Guppy)$solutions[2,2],summary(mcmc_FAS_RR_Guppy)$solutions[2,2]+summary(mcmc_FAS_RR_Guppy)$solutions[5,2]),
                         "bFAS_onto_high"=c(summary(mcmc_FAS_Guppy_useM)$solutions[2,3],summary(mcmc_FAS_RR_Guppy)$solutions[2,3],summary(mcmc_FAS_RR_Guppy)$solutions[2,3]+summary(mcmc_FAS_RR_Guppy)$solutions[5,3]),
                         
                           "bAGR_onto"=c(summary(mcmc_AGR_Guppy_useM)$solutions[2,1],summary(mcmc_AGR_RR_Guppy)$solutions[2,1],(summary(mcmc_AGR_RR_Guppy)$solutions[2,1]+summary(mcmc_AGR_RR_Guppy)$solutions[5,1])),
                           "bAGR_onto_low"=c(summary(mcmc_AGR_Guppy_useM)$solutions[2,2],summary(mcmc_AGR_RR_Guppy)$solutions[2,2],summary(mcmc_AGR_RR_Guppy)$solutions[2,2]+summary(mcmc_AGR_RR_Guppy)$solutions[5,2]),
                           "bAGR_onto_high"=c(summary(mcmc_AGR_Guppy_useM)$solutions[2,3],summary(mcmc_AGR_RR_Guppy)$solutions[2,3],summary(mcmc_AGR_RR_Guppy)$solutions[2,3]+summary(mcmc_AGR_RR_Guppy)$solutions[5,3]),
                           "bSGR_onto"=c(summary(mcmc_SGR_Guppy_useM)$solutions[2,1],summary(mcmc_SGR_RR_Guppy)$solutions[2,1],(summary(mcmc_SGR_RR_Guppy)$solutions[2,1]+summary(mcmc_SGR_RR_Guppy)$solutions[5,1])),
                           "bSGR_onto_low"=c(summary(mcmc_SGR_Guppy_useM)$solutions[2,2],summary(mcmc_SGR_RR_Guppy)$solutions[2,2],summary(mcmc_SGR_RR_Guppy)$solutions[2,2]+summary(mcmc_SGR_RR_Guppy)$solutions[5,2]),
                           "bSGR_onto_high"=c(summary(mcmc_SGR_Guppy_useM)$solutions[2,3],summary(mcmc_SGR_RR_Guppy)$solutions[2,3],summary(mcmc_SGR_RR_Guppy)$solutions[2,3]+summary(mcmc_SGR_RR_Guppy)$solutions[5,3]) ,
                           "bGE_onto"=c(summary(mcmc_GE_Guppy_useM)$solutions[2,1],summary(mcmc_GE_RR_Guppy)$solutions[2,1],(summary(mcmc_GE_RR_Guppy)$solutions[2,1]+summary(mcmc_GE_RR_Guppy)$solutions[5,1])),
                           "bGE_onto_low"=c(summary(mcmc_GE_Guppy_useM)$solutions[2,2],summary(mcmc_GE_RR_Guppy)$solutions[2,2],summary(mcmc_GE_RR_Guppy)$solutions[2,2]+summary(mcmc_GE_RR_Guppy)$solutions[5,2]),
                           "bGE_onto_high"=c(summary(mcmc_GE_Guppy_useM)$solutions[2,3],summary(mcmc_GE_RR_Guppy)$solutions[2,3],summary(mcmc_GE_RR_Guppy)$solutions[2,3]+summary(mcmc_GE_RR_Guppy)$solutions[5,3]),
                           
                         "bSMR_static"=c(summary(mcmc_SMR_Guppy_useM)$solutions[3,1],summary(mcmc_SMR_RR_Guppy)$solutions[4,1],(summary(mcmc_SMR_RR_Guppy)$solutions[4,1]+summary(mcmc_SMR_RR_Guppy)$solutions[6,1])),
                           "bSMR_static_low"=c(summary(mcmc_SMR_Guppy_useM)$solutions[3,2],summary(mcmc_SMR_RR_Guppy)$solutions[4,2],summary(mcmc_SMR_RR_Guppy)$solutions[4,2]+summary(mcmc_SMR_RR_Guppy)$solutions[6,2]),
                           "bSMR_static_high"=c(summary(mcmc_SMR_Guppy_useM)$solutions[3,3],summary(mcmc_SMR_RR_Guppy)$solutions[4,3],summary(mcmc_SMR_RR_Guppy)$solutions[4,3]+summary(mcmc_SMR_RR_Guppy)$solutions[6,3]),
                           "bMMR_static"=c(summary(mcmc_MMR_Guppy_useM)$solutions[3,1],summary(mcmc_MMR_RR_Guppy)$solutions[4,1],(summary(mcmc_MMR_RR_Guppy)$solutions[4,1]+summary(mcmc_MMR_RR_Guppy)$solutions[6,1])),
                           "bMMR_static_low"=c(summary(mcmc_MMR_Guppy_useM)$solutions[3,2],summary(mcmc_MMR_RR_Guppy)$solutions[4,2],summary(mcmc_MMR_RR_Guppy)$solutions[4,2]+summary(mcmc_MMR_RR_Guppy)$solutions[6,2]),
                           "bMMR_static_high"=c(summary(mcmc_MMR_Guppy_useM)$solutions[3,3],summary(mcmc_MMR_RR_Guppy)$solutions[4,3],summary(mcmc_MMR_RR_Guppy)$solutions[4,3]+summary(mcmc_MMR_RR_Guppy)$solutions[6,3]),
                           "bAS_static"=c(summary(mcmc_AS_Guppy_useM)$solutions[3,1],summary(mcmc_AS_RR_Guppy)$solutions[4,1],(summary(mcmc_AS_RR_Guppy)$solutions[4,1]+summary(mcmc_AS_RR_Guppy)$solutions[6,1])),
                           "bAS_static_low"=c(summary(mcmc_AS_Guppy_useM)$solutions[3,2],summary(mcmc_AS_RR_Guppy)$solutions[4,2],summary(mcmc_AS_RR_Guppy)$solutions[4,2]+summary(mcmc_AS_RR_Guppy)$solutions[6,2]),
                           "bAS_static_high"=c(summary(mcmc_AS_Guppy_useM)$solutions[3,3],summary(mcmc_AS_RR_Guppy)$solutions[4,3],summary(mcmc_AS_RR_Guppy)$solutions[4,3]+summary(mcmc_AS_RR_Guppy)$solutions[6,3]),
                         "bFAS_static"=c(summary(mcmc_FAS_Guppy_useM)$solutions[3,1],summary(mcmc_FAS_RR_Guppy)$solutions[4,1],(summary(mcmc_FAS_RR_Guppy)$solutions[4,1]+summary(mcmc_FAS_RR_Guppy)$solutions[6,1])),
                         "bFAS_static_low"=c(summary(mcmc_FAS_Guppy_useM)$solutions[3,2],summary(mcmc_FAS_RR_Guppy)$solutions[4,2],summary(mcmc_FAS_RR_Guppy)$solutions[4,2]+summary(mcmc_FAS_RR_Guppy)$solutions[6,2]),
                         "bFAS_static_high"=c(summary(mcmc_FAS_Guppy_useM)$solutions[3,3],summary(mcmc_FAS_RR_Guppy)$solutions[4,3],summary(mcmc_FAS_RR_Guppy)$solutions[4,3]+summary(mcmc_FAS_RR_Guppy)$solutions[6,3]),
                         
                           "bAGR_static"=c(summary(mcmc_AGR_Guppy_useM)$solutions[3,1],summary(mcmc_AGR_RR_Guppy)$solutions[4,1],(summary(mcmc_AGR_RR_Guppy)$solutions[4,1]+summary(mcmc_AGR_RR_Guppy)$solutions[6,1])),
                           "bAGR_static_low"=c(summary(mcmc_AGR_Guppy_useM)$solutions[3,2],summary(mcmc_AGR_RR_Guppy)$solutions[4,2],summary(mcmc_AGR_RR_Guppy)$solutions[4,2]+summary(mcmc_AGR_RR_Guppy)$solutions[6,2]),
                           "bAGR_static_high"=c(summary(mcmc_AGR_Guppy_useM)$solutions[3,3],summary(mcmc_AGR_RR_Guppy)$solutions[4,3],summary(mcmc_AGR_RR_Guppy)$solutions[4,3]+summary(mcmc_AGR_RR_Guppy)$solutions[6,3]),
                           "bSGR_static"=c(summary(mcmc_SGR_Guppy_useM)$solutions[3,1],summary(mcmc_SGR_RR_Guppy)$solutions[4,1],(summary(mcmc_SGR_RR_Guppy)$solutions[4,1]+summary(mcmc_SGR_RR_Guppy)$solutions[6,1])),
                           "bSGR_static_low"=c(summary(mcmc_SGR_Guppy_useM)$solutions[3,2],summary(mcmc_SGR_RR_Guppy)$solutions[4,2],summary(mcmc_SGR_RR_Guppy)$solutions[4,2]+summary(mcmc_SGR_RR_Guppy)$solutions[6,2]),
                           "bSGR_static_high"=c(summary(mcmc_SGR_Guppy_useM)$solutions[3,3],summary(mcmc_SGR_RR_Guppy)$solutions[4,3],summary(mcmc_SGR_RR_Guppy)$solutions[4,3]+summary(mcmc_SGR_RR_Guppy)$solutions[6,3]) ,
                           "bGE_static"=c(summary(mcmc_GE_Guppy_useM)$solutions[3,1],summary(mcmc_GE_RR_Guppy)$solutions[4,1],(summary(mcmc_GE_RR_Guppy)$solutions[4,1]+summary(mcmc_GE_RR_Guppy)$solutions[6,1])),
                           "bGE_static_low"=c(summary(mcmc_GE_Guppy_useM)$solutions[3,2],summary(mcmc_GE_RR_Guppy)$solutions[4,2],summary(mcmc_GE_RR_Guppy)$solutions[4,2]+summary(mcmc_GE_RR_Guppy)$solutions[6,2]),
                           "bGE_static_high"=c(summary(mcmc_GE_Guppy_useM)$solutions[3,3],summary(mcmc_GE_RR_Guppy)$solutions[4,3],summary(mcmc_GE_RR_Guppy)$solutions[4,3]+summary(mcmc_GE_RR_Guppy)$solutions[6,3])
)

Brown_trout_scal <- data.frame("Species"=rep("Brown_trout",1),"Treatment"=c("Mean") ,
                         "bSMR_onto"=c(summary(mcmc_SMR_RR_Brown_trout)$solutions[2,1]),
                         "bSMR_onto_low"=c(summary(mcmc_SMR_RR_Brown_trout)$solutions[2,2]),
                         "bSMR_onto_high"=c(summary(mcmc_SMR_RR_Brown_trout)$solutions[2,3]),
                         "bMMR_onto"=c(summary(mcmc_MMR_RR_Brown_trout)$solutions[2,1]),
                         "bMMR_onto_low"=c(summary(mcmc_MMR_RR_Brown_trout)$solutions[2,2]),
                         "bMMR_onto_high"=c(summary(mcmc_MMR_RR_Brown_trout)$solutions[2,3]),
                         "bAS_onto"=c(summary(mcmc_AS_RR_Brown_trout)$solutions[2,1]),
                         "bAS_onto_low"=c(summary(mcmc_AS_RR_Brown_trout)$solutions[2,2]),
                         "bAS_onto_high"=c(summary(mcmc_AS_RR_Brown_trout)$solutions[2,3]),
                         "bFAS_onto"=c(summary(mcmc_FAS_RR_Brown_trout)$solutions[2,1]),
                         "bFAS_onto_low"=c(summary(mcmc_FAS_RR_Brown_trout)$solutions[2,2]),
                         "bFAS_onto_high"=c(summary(mcmc_FAS_RR_Brown_trout)$solutions[2,3]),
                         
                         "bAGR_onto"=c(summary(mcmc_AGR_RR_Brown_trout)$solutions[2,1]),
                         "bAGR_onto_low"=c(summary(mcmc_AGR_RR_Brown_trout)$solutions[2,2]),
                         "bAGR_onto_high"=c(summary(mcmc_AGR_RR_Brown_trout)$solutions[2,3]),
                         "bSGR_onto"=c(summary(mcmc_SGR_RR_Brown_trout)$solutions[2,1]),
                         "bSGR_onto_low"=c(summary(mcmc_SGR_RR_Brown_trout)$solutions[2,2]),
                         "bSGR_onto_high"=c(summary(mcmc_SGR_RR_Brown_trout)$solutions[2,3]) ,
                         "bGE_onto"=c(summary(mcmc_GE_RR_Brown_trout)$solutions[2,1]),
                         "bGE_onto_low"=c(summary(mcmc_GE_RR_Brown_trout)$solutions[2,2]),
                         "bGE_onto_high"=c(summary(mcmc_GE_RR_Brown_trout)$solutions[2,3]),
                         
                         "bSMR_static"=c(summary(mcmc_SMR_RR_Brown_trout)$solutions[3,1]),
                         "bSMR_static_low"=c(summary(mcmc_SMR_RR_Brown_trout)$solutions[3,2]),
                         "bSMR_static_high"=c(summary(mcmc_SMR_RR_Brown_trout)$solutions[3,3]),
                         "bMMR_static"=c(summary(mcmc_MMR_RR_Brown_trout)$solutions[3,1]),
                         "bMMR_static_low"=c(summary(mcmc_MMR_RR_Brown_trout)$solutions[3,2]),
                         "bMMR_static_high"=c(summary(mcmc_MMR_RR_Brown_trout)$solutions[3,3]),
                         "bAS_static"=c(summary(mcmc_AS_RR_Brown_trout)$solutions[3,1]),
                         "bAS_static_low"=c(summary(mcmc_AS_RR_Brown_trout)$solutions[3,2]),
                         "bAS_static_high"=c(summary(mcmc_AS_RR_Brown_trout)$solutions[3,3]),
                         "bFAS_static"=c(summary(mcmc_FAS_RR_Brown_trout)$solutions[3,1]),
                         "bFAS_static_low"=c(summary(mcmc_FAS_RR_Brown_trout)$solutions[3,2]),
                         "bFAS_static_high"=c(summary(mcmc_FAS_RR_Brown_trout)$solutions[3,3]),
                         "bAGR_static"=c(summary(mcmc_AGR_RR_Brown_trout)$solutions[3,1]),
                         "bAGR_static_low"=c(summary(mcmc_AGR_RR_Brown_trout)$solutions[3,2]),
                         "bAGR_static_high"=c(summary(mcmc_AGR_RR_Brown_trout)$solutions[3,3]),
                         "bSGR_static"=c(summary(mcmc_SGR_RR_Brown_trout)$solutions[3,1]),
                         "bSGR_static_low"=c(summary(mcmc_SGR_RR_Brown_trout)$solutions[3,2]),
                         "bSGR_static_high"=c(summary(mcmc_SGR_RR_Brown_trout)$solutions[3,3]) ,
                         "bGE_static"=c(summary(mcmc_GE_RR_Brown_trout)$solutions[3,1]),
                         "bGE_static_low"=c(summary(mcmc_GE_RR_Brown_trout)$solutions[3,2]),
                         "bGE_static_high"=c(summary(mcmc_GE_RR_Brown_trout)$solutions[3,3])
)
Cunner_scal <- data.frame("Species"=rep("Cunner",1),"Treatment"=c("Mean") ,
                               "bSMR_onto"=c(summary(mcmc_SMR_RR_Cunner)$solutions[2,1]),
                               "bSMR_onto_low"=c(summary(mcmc_SMR_RR_Cunner)$solutions[2,2]),
                               "bSMR_onto_high"=c(summary(mcmc_SMR_RR_Cunner)$solutions[2,3]),
                               "bMMR_onto"=c(summary(mcmc_MMR_RR_Cunner)$solutions[2,1]),
                               "bMMR_onto_low"=c(summary(mcmc_MMR_RR_Cunner)$solutions[2,2]),
                               "bMMR_onto_high"=c(summary(mcmc_MMR_RR_Cunner)$solutions[2,3]),
                               "bAS_onto"=c(summary(mcmc_AS_RR_Cunner)$solutions[2,1]),
                               "bAS_onto_low"=c(summary(mcmc_AS_RR_Cunner)$solutions[2,2]),
                               "bAS_onto_high"=c(summary(mcmc_AS_RR_Cunner)$solutions[2,3]),
                          "bFAS_onto"=c(summary(mcmc_FAS_RR_Cunner)$solutions[2,1]),
                          "bFAS_onto_low"=c(summary(mcmc_FAS_RR_Cunner)$solutions[2,2]),
                          "bFAS_onto_high"=c(summary(mcmc_FAS_RR_Cunner)$solutions[2,3]),
                          
                               "bAGR_onto"=c(summary(mcmc_AGR_RR_Cunner)$solutions[2,1]),
                               "bAGR_onto_low"=c(summary(mcmc_AGR_RR_Cunner)$solutions[2,2]),
                               "bAGR_onto_high"=c(summary(mcmc_AGR_RR_Cunner)$solutions[2,3]),
                               "bSGR_onto"=c(summary(mcmc_SGR_RR_Cunner)$solutions[2,1]),
                               "bSGR_onto_low"=c(summary(mcmc_SGR_RR_Cunner)$solutions[2,2]),
                               "bSGR_onto_high"=c(summary(mcmc_SGR_RR_Cunner)$solutions[2,3]) ,
                               "bGE_onto"=c(summary(mcmc_GE_RR_Cunner)$solutions[2,1]),
                               "bGE_onto_low"=c(summary(mcmc_GE_RR_Cunner)$solutions[2,2]),
                               "bGE_onto_high"=c(summary(mcmc_GE_RR_Cunner)$solutions[2,3]),
                               
                               "bSMR_static"=c(summary(mcmc_SMR_RR_Cunner)$solutions[3,1]),
                               "bSMR_static_low"=c(summary(mcmc_SMR_RR_Cunner)$solutions[3,2]),
                               "bSMR_static_high"=c(summary(mcmc_SMR_RR_Cunner)$solutions[3,3]),
                               "bMMR_static"=c(summary(mcmc_MMR_RR_Cunner)$solutions[3,1]),
                               "bMMR_static_low"=c(summary(mcmc_MMR_RR_Cunner)$solutions[3,2]),
                               "bMMR_static_high"=c(summary(mcmc_MMR_RR_Cunner)$solutions[3,3]),
                               "bAS_static"=c(summary(mcmc_AS_RR_Cunner)$solutions[3,1]),
                               "bAS_static_low"=c(summary(mcmc_AS_RR_Cunner)$solutions[3,2]),
                               "bAS_static_high"=c(summary(mcmc_AS_RR_Cunner)$solutions[3,3]),
                          "bFAS_static"=c(summary(mcmc_FAS_RR_Cunner)$solutions[3,1]),
                          "bFAS_static_low"=c(summary(mcmc_FAS_RR_Cunner)$solutions[3,2]),
                          "bFAS_static_high"=c(summary(mcmc_FAS_RR_Cunner)$solutions[3,3]),
                          
                               "bAGR_static"=c(summary(mcmc_AGR_RR_Cunner)$solutions[3,1]),
                               "bAGR_static_low"=c(summary(mcmc_AGR_RR_Cunner)$solutions[3,2]),
                               "bAGR_static_high"=c(summary(mcmc_AGR_RR_Cunner)$solutions[3,3]),
                               "bSGR_static"=c(summary(mcmc_SGR_RR_Cunner)$solutions[3,1]),
                               "bSGR_static_low"=c(summary(mcmc_SGR_RR_Cunner)$solutions[3,2]),
                               "bSGR_static_high"=c(summary(mcmc_SGR_RR_Cunner)$solutions[3,3]) ,
                               "bGE_static"=c(summary(mcmc_GE_RR_Cunner)$solutions[3,1]),
                               "bGE_static_low"=c(summary(mcmc_GE_RR_Cunner)$solutions[3,2]),
                               "bGE_static_high"=c(summary(mcmc_GE_RR_Cunner)$solutions[3,3])
)
Clownfish_scal <- data.frame("Species"=rep("Clownfish",1),"Treatment"=c("Mean") ,
                          "bSMR_onto"=c(summary(mcmc_SMR_Clownfish_use)$solutions[2,1]),
                          "bSMR_onto_low"=c(summary(mcmc_SMR_Clownfish_use)$solutions[2,2]),
                          "bSMR_onto_high"=c(summary(mcmc_SMR_Clownfish_use)$solutions[2,3]),
                          "bMMR_onto"=c(summary(mcmc_MMR_Clownfish_use)$solutions[2,1]),
                          "bMMR_onto_low"=c(summary(mcmc_MMR_Clownfish_use)$solutions[2,2]),
                          "bMMR_onto_high"=c(summary(mcmc_MMR_Clownfish_use)$solutions[2,3]),
                          "bAS_onto"=c(summary(mcmc_AS_Clownfish_use)$solutions[2,1]),
                          "bAS_onto_low"=c(summary(mcmc_AS_Clownfish_use)$solutions[2,2]),
                          "bAS_onto_high"=c(summary(mcmc_AS_Clownfish_use)$solutions[2,3]),
                          "bFAS_onto"=c(summary(mcmc_FAS_Clownfish_use)$solutions[2,1]),
                          "bFAS_onto_low"=c(summary(mcmc_FAS_Clownfish_use)$solutions[2,2]),
                          "bFAS_onto_high"=c(summary(mcmc_FAS_Clownfish_use)$solutions[2,3]),
                          
                          "bAGR_onto"=c(summary(mcmc_AGR_Clownfish_use)$solutions[2,1]),
                          "bAGR_onto_low"=c(summary(mcmc_AGR_Clownfish_use)$solutions[2,2]),
                          "bAGR_onto_high"=c(summary(mcmc_AGR_Clownfish_use)$solutions[2,3]),
                          "bSGR_onto"=c(summary(mcmc_SGR_Clownfish_use)$solutions[2,1]),
                          "bSGR_onto_low"=c(summary(mcmc_SGR_Clownfish_use)$solutions[2,2]),
                          "bSGR_onto_high"=c(summary(mcmc_SGR_Clownfish_use)$solutions[2,3]) ,
                          "bGE_onto"=c(summary(mcmc_GE_Clownfish_use)$solutions[2,1]),
                          "bGE_onto_low"=c(summary(mcmc_GE_Clownfish_use)$solutions[2,2]),
                          "bGE_onto_high"=c(summary(mcmc_GE_Clownfish_use)$solutions[2,3]),
                          
                          "bSMR_static"=c(summary(mcmc_SMR_Clownfish_use)$solutions[3,1]),
                          "bSMR_static_low"=c(summary(mcmc_SMR_Clownfish_use)$solutions[3,2]),
                          "bSMR_static_high"=c(summary(mcmc_SMR_Clownfish_use)$solutions[3,3]),
                          "bMMR_static"=c(summary(mcmc_MMR_Clownfish_use)$solutions[3,1]),
                          "bMMR_static_low"=c(summary(mcmc_MMR_Clownfish_use)$solutions[3,2]),
                          "bMMR_static_high"=c(summary(mcmc_MMR_Clownfish_use)$solutions[3,3]),
                          "bAS_static"=c(summary(mcmc_AS_Clownfish_use)$solutions[3,1]),
                          "bAS_static_low"=c(summary(mcmc_AS_Clownfish_use)$solutions[3,2]),
                          "bAS_static_high"=c(summary(mcmc_AS_Clownfish_use)$solutions[3,3]),
                          "bFAS_static"=c(summary(mcmc_FAS_Clownfish_use)$solutions[3,1]),
                          "bFAS_static_low"=c(summary(mcmc_FAS_Clownfish_use)$solutions[3,2]),
                          "bFAS_static_high"=c(summary(mcmc_FAS_Clownfish_use)$solutions[3,3]),
                          
                          "bAGR_static"=c(summary(mcmc_AGR_Clownfish_use)$solutions[3,1]),
                          "bAGR_static_low"=c(summary(mcmc_AGR_Clownfish_use)$solutions[3,2]),
                          "bAGR_static_high"=c(summary(mcmc_AGR_Clownfish_use)$solutions[3,3]),
                          "bSGR_static"=c(summary(mcmc_SGR_Clownfish_use)$solutions[3,1]),
                          "bSGR_static_low"=c(summary(mcmc_SGR_Clownfish_use)$solutions[3,2]),
                          "bSGR_static_high"=c(summary(mcmc_SGR_Clownfish_use)$solutions[3,3]) ,
                          "bGE_static"=c(summary(mcmc_GE_Clownfish_use)$solutions[3,1]),
                          "bGE_static_low"=c(summary(mcmc_GE_Clownfish_use)$solutions[3,2]),
                          "bGE_static_high"=c(summary(mcmc_GE_Clownfish_use)$solutions[3,3])
)

Fish_scal <- rbind(Zebrafish_scal,rbind(Rainbow_trout_scal,rbind(Brown_trout_scal,rbind(Cunner_scal,rbind(Guppy_scal,rbind(Clownfish_scal, Chromis_scal))))))

### Finding the inting exponents for each species for each paremeter
##Zebrafish##
Zebrafish_int <- data.frame("Species"=rep("Zebrafish",4),"Treatment"=c("Mean","High_fed","Low_fed","Moderate_fed") ,
                             "aSMR"=c(summary(mcmc_SMR_Zebrafish_useM)$solutions[1,1],summary(mcmc_SMR_Zebrafish_use)$solutions[1,1],(summary(mcmc_SMR_Zebrafish_use)$solutions[1,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[3,1]),(summary(mcmc_SMR_Zebrafish_use)$solutions[1,1]+summary(mcmc_SMR_Zebrafish_use)$solutions[4,1])),
                             "aSMR_low"=c(summary(mcmc_SMR_Zebrafish_useM)$solutions[1,2],summary(mcmc_SMR_Zebrafish_use)$solutions[1,2],summary(mcmc_SMR_Zebrafish_use)$solutions[1,2]+summary(mcmc_SMR_Zebrafish_use)$solutions[3,2],summary(mcmc_SMR_Zebrafish_use)$solutions[1,2]+summary(mcmc_SMR_Zebrafish_use)$solutions[4,2]),
                             "aSMR_high"=c(summary(mcmc_SMR_Zebrafish_useM)$solutions[1,3],summary(mcmc_SMR_Zebrafish_use)$solutions[1,3],summary(mcmc_SMR_Zebrafish_use)$solutions[1,3]+summary(mcmc_SMR_Zebrafish_use)$solutions[3,3],summary(mcmc_SMR_Zebrafish_use)$solutions[1,3]+summary(mcmc_SMR_Zebrafish_use)$solutions[4,3]),
                             "aMMR"=c(NA,NA,NA,NA),"aMMR_low"=c(NA,NA,NA,NA),"aMMR_high"=c(NA,NA,NA,NA) ,
                             "aAS"=c(NA,NA,NA,NA),"aAS_low"=c(NA,NA,NA,NA),"aAS_high"=c(NA,NA,NA,NA),
                            "aFAS"=c(NA,NA,NA,NA),"aFAS_low"=c(NA,NA,NA,NA),"aFAS_high"=c(NA,NA,NA,NA),
                             "aAGR"=c(summary(mcmc_AGR_Zebrafish_useM)$solutions[1,1],summary(mcmc_AGR_Zebrafish_use)$solutions[1,1],(summary(mcmc_AGR_Zebrafish_use)$solutions[1,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[3,1]),(summary(mcmc_AGR_Zebrafish_use)$solutions[1,1]+summary(mcmc_AGR_Zebrafish_use)$solutions[4,1])),
                             "aAGR_low"=c(summary(mcmc_AGR_Zebrafish_useM)$solutions[1,2],summary(mcmc_AGR_Zebrafish_use)$solutions[1,2],summary(mcmc_AGR_Zebrafish_use)$solutions[1,2]+summary(mcmc_AGR_Zebrafish_use)$solutions[3,2],summary(mcmc_AGR_Zebrafish_use)$solutions[1,2]+summary(mcmc_AGR_Zebrafish_use)$solutions[4,2]),
                             "aAGR_high"=c(summary(mcmc_AGR_Zebrafish_useM)$solutions[1,3],summary(mcmc_AGR_Zebrafish_use)$solutions[1,3],summary(mcmc_AGR_Zebrafish_use)$solutions[1,3]+summary(mcmc_AGR_Zebrafish_use)$solutions[3,3],summary(mcmc_AGR_Zebrafish_use)$solutions[1,3]+summary(mcmc_AGR_Zebrafish_use)$solutions[4,3]),
                             "aSGR"=c(summary(mcmc_SGR_Zebrafish_useM)$solutions[1,1],summary(mcmc_SGR_Zebrafish_use)$solutions[1,1],(summary(mcmc_SGR_Zebrafish_use)$solutions[1,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[3,1]),(summary(mcmc_SGR_Zebrafish_use)$solutions[1,1]+summary(mcmc_SGR_Zebrafish_use)$solutions[4,1])),
                             "aSGR_low"=c(summary(mcmc_SGR_Zebrafish_useM)$solutions[1,2],summary(mcmc_SGR_Zebrafish_use)$solutions[1,2],summary(mcmc_SGR_Zebrafish_use)$solutions[1,2]+summary(mcmc_SGR_Zebrafish_use)$solutions[3,2],summary(mcmc_SGR_Zebrafish_use)$solutions[1,2]+summary(mcmc_SGR_Zebrafish_use)$solutions[4,2]),
                             "aSGR_high"=c(summary(mcmc_SGR_Zebrafish_useM)$solutions[1,3],summary(mcmc_SGR_Zebrafish_use)$solutions[1,3],summary(mcmc_SGR_Zebrafish_use)$solutions[1,3]+summary(mcmc_SGR_Zebrafish_use)$solutions[3,3],summary(mcmc_SGR_Zebrafish_use)$solutions[1,3]+summary(mcmc_SGR_Zebrafish_use)$solutions[4,3]) ,
                             "aGE"=c(summary(mcmc_GE_Zebrafish_useM)$solutions[1,1],summary(mcmc_GE_Zebrafish_use)$solutions[1,1],(summary(mcmc_GE_Zebrafish_use)$solutions[1,1]+summary(mcmc_GE_Zebrafish_use)$solutions[3,1]),(summary(mcmc_GE_Zebrafish_use)$solutions[1,1]+summary(mcmc_GE_Zebrafish_use)$solutions[4,1])),
                             "aGE_low"=c(summary(mcmc_GE_Zebrafish_useM)$solutions[1,2],summary(mcmc_GE_Zebrafish_use)$solutions[1,2],summary(mcmc_GE_Zebrafish_use)$solutions[1,2]+summary(mcmc_GE_Zebrafish_use)$solutions[3,2],summary(mcmc_GE_Zebrafish_use)$solutions[1,2]+summary(mcmc_GE_Zebrafish_use)$solutions[4,2]),
                             "aGE_high"=c(summary(mcmc_GE_Zebrafish_useM)$solutions[1,3],summary(mcmc_GE_Zebrafish_use)$solutions[1,3],summary(mcmc_GE_Zebrafish_use)$solutions[1,3]+summary(mcmc_GE_Zebrafish_use)$solutions[3,3],summary(mcmc_GE_Zebrafish_use)$solutions[1,3]+summary(mcmc_GE_Zebrafish_use)$solutions[4,3]))


##Rainbow trout##
Rainbow_trout_int <- data.frame("Species"=rep("Rainbow_trout",4),"Treatment"=c("Mean","10Cegg","14Cegg","14Cyolk") ,
                                 "aSMR"=c(summary(mcmc_SMR_Rainbow_trout_useM)$solutions[1,1],summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,1],(summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[3,1]),(summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[4,1])),
                                 "aSMR_low"=c(summary(mcmc_SMR_Rainbow_trout_useM)$solutions[1,2],summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,2],summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[3,2],summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[4,2]),
                                 "aSMR_high"=c(summary(mcmc_SMR_Rainbow_trout_useM)$solutions[1,3],summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,3],summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[3,3],summary(mcmc_SMR_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_SMR_Rainbow_trout_use)$solutions[4,3]),
                                 "aMMR"=c(summary(mcmc_MMR_Rainbow_trout_useM)$solutions[1,1],summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,1],(summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[3,1]),(summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[4,1])),
                                 "aMMR_low"=c(summary(mcmc_MMR_Rainbow_trout_useM)$solutions[1,2],summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,2],summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[3,2],summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[4,2]),
                                 "aMMR_high"=c(summary(mcmc_MMR_Rainbow_trout_useM)$solutions[1,3],summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,3],summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[3,3],summary(mcmc_MMR_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_MMR_Rainbow_trout_use)$solutions[4,3]),
                                 "aAS"=c(summary(mcmc_AS_Rainbow_trout_useM)$solutions[1,1],summary(mcmc_AS_Rainbow_trout_use)$solutions[1,1],(summary(mcmc_AS_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[3,1]),(summary(mcmc_AS_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_AS_Rainbow_trout_use)$solutions[4,1])),
                                 "aAS_low"=c(summary(mcmc_AS_Rainbow_trout_useM)$solutions[1,2],summary(mcmc_AS_Rainbow_trout_use)$solutions[1,2],summary(mcmc_AS_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_AS_Rainbow_trout_use)$solutions[3,2],summary(mcmc_AS_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_AS_Rainbow_trout_use)$solutions[4,2]),
                                 "aAS_high"=c(summary(mcmc_AS_Rainbow_trout_useM)$solutions[1,3],summary(mcmc_AS_Rainbow_trout_use)$solutions[1,3],summary(mcmc_AS_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_AS_Rainbow_trout_use)$solutions[3,3],summary(mcmc_AS_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_AS_Rainbow_trout_use)$solutions[4,3]),
                                "aFAS"=c(summary(mcmc_FAS_Rainbow_trout_useM)$solutions[1,1],summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,1],(summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[3,1]),(summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[4,1])),
                                "aFAS_low"=c(summary(mcmc_FAS_Rainbow_trout_useM)$solutions[1,2],summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,2],summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[3,2],summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[4,2]),
                                "aFAS_high"=c(summary(mcmc_FAS_Rainbow_trout_useM)$solutions[1,3],summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,3],summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[3,3],summary(mcmc_FAS_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_FAS_Rainbow_trout_use)$solutions[4,3]),
                                
                                 "aAGR"=c(summary(mcmc_AGR_Rainbow_trout_useM)$solutions[1,1],summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,1],(summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[3,1]),(summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[4,1])),
                                 "aAGR_low"=c(summary(mcmc_AGR_Rainbow_trout_useM)$solutions[1,2],summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,2],summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[3,2],summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[4,2]),
                                 "aAGR_high"=c(summary(mcmc_AGR_Rainbow_trout_useM)$solutions[1,3],summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,3],summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[3,3],summary(mcmc_AGR_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_AGR_Rainbow_trout_use)$solutions[4,3]),
                                 "aSGR"=c(summary(mcmc_SGR_Rainbow_trout_useM)$solutions[1,1],summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,1],(summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[3,1]),(summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[4,1])),
                                 "aSGR_low"=c(summary(mcmc_SGR_Rainbow_trout_useM)$solutions[1,2],summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,2],summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[3,2],summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[4,2]),
                                 "aSGR_high"=c(summary(mcmc_SGR_Rainbow_trout_useM)$solutions[1,3],summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,3],summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[3,3],summary(mcmc_SGR_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_SGR_Rainbow_trout_use)$solutions[4,3]) ,
                                 "aGE"=c(summary(mcmc_GE_Rainbow_trout_useM)$solutions[1,1],summary(mcmc_GE_Rainbow_trout_use)$solutions[1,1],(summary(mcmc_GE_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[3,1]),(summary(mcmc_GE_Rainbow_trout_use)$solutions[1,1]+summary(mcmc_GE_Rainbow_trout_use)$solutions[4,1])),
                                 "aGE_low"=c(summary(mcmc_GE_Rainbow_trout_useM)$solutions[1,2],summary(mcmc_GE_Rainbow_trout_use)$solutions[1,2],summary(mcmc_GE_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_GE_Rainbow_trout_use)$solutions[3,2],summary(mcmc_GE_Rainbow_trout_use)$solutions[1,2]+summary(mcmc_GE_Rainbow_trout_use)$solutions[4,2]),
                                 "aGE_high"=c(summary(mcmc_GE_Rainbow_trout_useM)$solutions[1,3],summary(mcmc_GE_Rainbow_trout_use)$solutions[1,3],summary(mcmc_GE_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_GE_Rainbow_trout_use)$solutions[3,3],summary(mcmc_GE_Rainbow_trout_use)$solutions[1,3]+summary(mcmc_GE_Rainbow_trout_use)$solutions[4,3]))
Chromis_int <- data.frame("Species"=rep("Chromis",3),"Treatment"=c("Mean","High_fed","Low_fed") ,
                           "aSMR"=c(summary(mcmc_SMR_Chromis_useM)$solutions[1,1],summary(mcmc_SMR_RR_Chromis)$solutions[1,1],(summary(mcmc_SMR_RR_Chromis)$solutions[1,1]+summary(mcmc_SMR_RR_Chromis)$solutions[3,1])),
                           "aSMR_low"=c(summary(mcmc_SMR_Chromis_useM)$solutions[1,2],summary(mcmc_SMR_RR_Chromis)$solutions[1,2],summary(mcmc_SMR_RR_Chromis)$solutions[1,2]+summary(mcmc_SMR_RR_Chromis)$solutions[3,2]),
                           "aSMR_high"=c(summary(mcmc_SMR_Chromis_useM)$solutions[1,3],summary(mcmc_SMR_RR_Chromis)$solutions[1,3],summary(mcmc_SMR_RR_Chromis)$solutions[1,3]+summary(mcmc_SMR_RR_Chromis)$solutions[3,3]),
                           "aMMR"=c(summary(mcmc_MMR_Chromis_useM)$solutions[1,1],summary(mcmc_MMR_RR_Chromis)$solutions[1,1],(summary(mcmc_MMR_RR_Chromis)$solutions[1,1]+summary(mcmc_MMR_RR_Chromis)$solutions[3,1])),
                           "aMMR_low"=c(summary(mcmc_MMR_Chromis_useM)$solutions[1,2],summary(mcmc_MMR_RR_Chromis)$solutions[1,2],summary(mcmc_MMR_RR_Chromis)$solutions[1,2]+summary(mcmc_MMR_RR_Chromis)$solutions[3,2]),
                           "aMMR_high"=c(summary(mcmc_MMR_Chromis_useM)$solutions[1,3],summary(mcmc_MMR_RR_Chromis)$solutions[1,3],summary(mcmc_MMR_RR_Chromis)$solutions[1,3]+summary(mcmc_MMR_RR_Chromis)$solutions[3,3]),
                           "aAS"=c(summary(mcmc_AS_Chromis_useM)$solutions[1,1],summary(mcmc_AS_RR_Chromis)$solutions[1,1],(summary(mcmc_AS_RR_Chromis)$solutions[1,1]+summary(mcmc_AS_RR_Chromis)$solutions[3,1])),
                           "aAS_low"=c(summary(mcmc_AS_Chromis_useM)$solutions[1,2],summary(mcmc_AS_RR_Chromis)$solutions[1,2],summary(mcmc_AS_RR_Chromis)$solutions[1,2]+summary(mcmc_AS_RR_Chromis)$solutions[3,2]),
                           "aAS_high"=c(summary(mcmc_AS_Chromis_useM)$solutions[1,3],summary(mcmc_AS_RR_Chromis)$solutions[1,3],summary(mcmc_AS_RR_Chromis)$solutions[1,3]+summary(mcmc_AS_RR_Chromis)$solutions[3,3]),
                          "aFAS"=c(summary(mcmc_FAS_Chromis_useM)$solutions[1,1],summary(mcmc_FAS_RR_Chromis)$solutions[1,1],(summary(mcmc_FAS_RR_Chromis)$solutions[1,1]+summary(mcmc_FAS_RR_Chromis)$solutions[3,1])),
                          "aFAS_low"=c(summary(mcmc_FAS_Chromis_useM)$solutions[1,2],summary(mcmc_FAS_RR_Chromis)$solutions[1,2],summary(mcmc_FAS_RR_Chromis)$solutions[1,2]+summary(mcmc_FAS_RR_Chromis)$solutions[3,2]),
                          "aFAS_high"=c(summary(mcmc_FAS_Chromis_useM)$solutions[1,3],summary(mcmc_FAS_RR_Chromis)$solutions[1,3],summary(mcmc_FAS_RR_Chromis)$solutions[1,3]+summary(mcmc_FAS_RR_Chromis)$solutions[3,3]),
                          
                           "aAGR"=c(summary(mcmc_AGR_Chromis_useM)$solutions[1,1],summary(mcmc_AGR_RR_Chromis)$solutions[1,1],(summary(mcmc_AGR_RR_Chromis)$solutions[1,1]+summary(mcmc_AGR_RR_Chromis)$solutions[3,1])),
                           "aAGR_low"=c(summary(mcmc_AGR_Chromis_useM)$solutions[1,2],summary(mcmc_AGR_RR_Chromis)$solutions[1,2],summary(mcmc_AGR_RR_Chromis)$solutions[1,2]+summary(mcmc_AGR_RR_Chromis)$solutions[3,2]),
                           "aAGR_high"=c(summary(mcmc_AGR_Chromis_useM)$solutions[1,3],summary(mcmc_AGR_RR_Chromis)$solutions[1,3],summary(mcmc_AGR_RR_Chromis)$solutions[1,3]+summary(mcmc_AGR_RR_Chromis)$solutions[3,3]),
                           "aSGR"=c(summary(mcmc_SGR_Chromis_useM)$solutions[1,1],summary(mcmc_SGR_RR_Chromis)$solutions[1,1],(summary(mcmc_SGR_RR_Chromis)$solutions[1,1]+summary(mcmc_SGR_RR_Chromis)$solutions[3,1])),
                           "aSGR_low"=c(summary(mcmc_SGR_Chromis_useM)$solutions[1,2],summary(mcmc_SGR_RR_Chromis)$solutions[1,2],summary(mcmc_SGR_RR_Chromis)$solutions[1,2]+summary(mcmc_SGR_RR_Chromis)$solutions[3,2]),
                           "aSGR_high"=c(summary(mcmc_SGR_Chromis_useM)$solutions[1,3],summary(mcmc_SGR_RR_Chromis)$solutions[1,3],summary(mcmc_SGR_RR_Chromis)$solutions[2,3]+summary(mcmc_SGR_RR_Chromis)$solutions[3,3]) ,
                           "aGE"=c(summary(mcmc_GE_Chromis_useM)$solutions[1,1],summary(mcmc_GE_RR_Chromis)$solutions[1,1],(summary(mcmc_GE_RR_Chromis)$solutions[1,1]+summary(mcmc_GE_RR_Chromis)$solutions[3,1])),
                           "aGE_low"=c(summary(mcmc_GE_Chromis_useM)$solutions[1,2],summary(mcmc_GE_RR_Chromis)$solutions[1,2],summary(mcmc_GE_RR_Chromis)$solutions[1,2]+summary(mcmc_GE_RR_Chromis)$solutions[3,2]),
                           "aGE_high"=c(summary(mcmc_GE_Chromis_useM)$solutions[1,3],summary(mcmc_GE_RR_Chromis)$solutions[1,3],summary(mcmc_GE_RR_Chromis)$solutions[1,3]+summary(mcmc_GE_RR_Chromis)$solutions[3,3]))

Guppy_int <- data.frame("Species"=rep("Guppy",3),"Treatment"=c("Mean","Female","Male") ,
                         "aSMR"=c(summary(mcmc_SMR_Guppy_useM)$solutions[1,1],summary(mcmc_SMR_RR_Guppy)$solutions[1,1],(summary(mcmc_SMR_RR_Guppy)$solutions[1,1]+summary(mcmc_SMR_RR_Guppy)$solutions[3,1])),
                         "aSMR_low"=c(summary(mcmc_SMR_Guppy_useM)$solutions[1,2],summary(mcmc_SMR_RR_Guppy)$solutions[1,2],summary(mcmc_SMR_RR_Guppy)$solutions[1,2]+summary(mcmc_SMR_RR_Guppy)$solutions[3,2]),
                         "aSMR_high"=c(summary(mcmc_SMR_Guppy_useM)$solutions[1,3],summary(mcmc_SMR_RR_Guppy)$solutions[1,3],summary(mcmc_SMR_RR_Guppy)$solutions[1,3]+summary(mcmc_SMR_RR_Guppy)$solutions[3,3]),
                         "aMMR"=c(summary(mcmc_MMR_Guppy_useM)$solutions[1,1],summary(mcmc_MMR_RR_Guppy)$solutions[1,1],(summary(mcmc_MMR_RR_Guppy)$solutions[1,1]+summary(mcmc_MMR_RR_Guppy)$solutions[3,1])),
                         "aMMR_low"=c(summary(mcmc_MMR_Guppy_useM)$solutions[1,2],summary(mcmc_MMR_RR_Guppy)$solutions[1,2],summary(mcmc_MMR_RR_Guppy)$solutions[1,2]+summary(mcmc_MMR_RR_Guppy)$solutions[3,2]),
                         "aMMR_high"=c(summary(mcmc_MMR_Guppy_useM)$solutions[1,3],summary(mcmc_MMR_RR_Guppy)$solutions[1,3],summary(mcmc_MMR_RR_Guppy)$solutions[1,3]+summary(mcmc_MMR_RR_Guppy)$solutions[3,3]),
                         "aAS"=c(summary(mcmc_AS_Guppy_useM)$solutions[1,1],summary(mcmc_AS_RR_Guppy)$solutions[1,1],(summary(mcmc_AS_RR_Guppy)$solutions[1,1]+summary(mcmc_AS_RR_Guppy)$solutions[3,1])),
                         "aAS_low"=c(summary(mcmc_AS_Guppy_useM)$solutions[1,2],summary(mcmc_AS_RR_Guppy)$solutions[1,2],summary(mcmc_AS_RR_Guppy)$solutions[1,2]+summary(mcmc_AS_RR_Guppy)$solutions[3,2]),
                         "aAS_high"=c(summary(mcmc_AS_Guppy_useM)$solutions[1,3],summary(mcmc_AS_RR_Guppy)$solutions[1,3],summary(mcmc_AS_RR_Guppy)$solutions[1,3]+summary(mcmc_AS_RR_Guppy)$solutions[3,3]),
                        "aFAS"=c(summary(mcmc_FAS_Guppy_useM)$solutions[1,1],summary(mcmc_FAS_RR_Guppy)$solutions[1,1],(summary(mcmc_FAS_RR_Guppy)$solutions[1,1]+summary(mcmc_FAS_RR_Guppy)$solutions[3,1])),
                        "aFAS_low"=c(summary(mcmc_FAS_Guppy_useM)$solutions[1,2],summary(mcmc_FAS_RR_Guppy)$solutions[1,2],summary(mcmc_FAS_RR_Guppy)$solutions[1,2]+summary(mcmc_FAS_RR_Guppy)$solutions[3,2]),
                        "aFAS_high"=c(summary(mcmc_FAS_Guppy_useM)$solutions[1,3],summary(mcmc_FAS_RR_Guppy)$solutions[1,3],summary(mcmc_FAS_RR_Guppy)$solutions[1,3]+summary(mcmc_FAS_RR_Guppy)$solutions[3,3]),
                        
                         "aAGR"=c(summary(mcmc_AGR_Guppy_useM)$solutions[1,1],summary(mcmc_AGR_RR_Guppy)$solutions[1,1],(summary(mcmc_AGR_RR_Guppy)$solutions[1,1]+summary(mcmc_AGR_RR_Guppy)$solutions[3,1])),
                         "aAGR_low"=c(summary(mcmc_AGR_Guppy_useM)$solutions[1,2],summary(mcmc_AGR_RR_Guppy)$solutions[1,2],summary(mcmc_AGR_RR_Guppy)$solutions[1,2]+summary(mcmc_AGR_RR_Guppy)$solutions[3,2]),
                         "aAGR_high"=c(summary(mcmc_AGR_Guppy_useM)$solutions[1,3],summary(mcmc_AGR_RR_Guppy)$solutions[1,3],summary(mcmc_AGR_RR_Guppy)$solutions[1,3]+summary(mcmc_AGR_RR_Guppy)$solutions[3,3]),
                         "aSGR"=c(summary(mcmc_SGR_Guppy_useM)$solutions[1,1],summary(mcmc_SGR_RR_Guppy)$solutions[1,1],(summary(mcmc_SGR_RR_Guppy)$solutions[1,1]+summary(mcmc_SGR_RR_Guppy)$solutions[3,1])),
                         "aSGR_low"=c(summary(mcmc_SGR_Guppy_useM)$solutions[1,2],summary(mcmc_SGR_RR_Guppy)$solutions[1,2],summary(mcmc_SGR_RR_Guppy)$solutions[1,2]+summary(mcmc_SGR_RR_Guppy)$solutions[3,2]),
                         "aSGR_high"=c(summary(mcmc_SGR_Guppy_useM)$solutions[1,3],summary(mcmc_SGR_RR_Guppy)$solutions[1,3],summary(mcmc_SGR_RR_Guppy)$solutions[1,3]+summary(mcmc_SGR_RR_Guppy)$solutions[3,3]) ,
                         "aGE"=c(summary(mcmc_GE_Guppy_useM)$solutions[1,1],summary(mcmc_GE_RR_Guppy)$solutions[1,1],(summary(mcmc_GE_RR_Guppy)$solutions[1,1]+summary(mcmc_GE_RR_Guppy)$solutions[3,1])),
                         "aGE_low"=c(summary(mcmc_GE_Guppy_useM)$solutions[1,2],summary(mcmc_GE_RR_Guppy)$solutions[1,2],summary(mcmc_GE_RR_Guppy)$solutions[1,2]+summary(mcmc_GE_RR_Guppy)$solutions[3,1]),
                         "aGE_high"=c(summary(mcmc_GE_Guppy_useM)$solutions[1,3],summary(mcmc_GE_RR_Guppy)$solutions[1,3],summary(mcmc_GE_RR_Guppy)$solutions[1,3]+summary(mcmc_GE_RR_Guppy)$solutions[3,3]))

Brown_trout_int <- data.frame("Species"=rep("Brown_trout",1),"Treatment"=c("Mean") ,
                               "aSMR"=c(summary(mcmc_SMR_RR_Brown_trout)$solutions[1,1]),
                               "aSMR_low"=c(summary(mcmc_SMR_RR_Brown_trout)$solutions[1,2]),
                               "aSMR_high"=c(summary(mcmc_SMR_RR_Brown_trout)$solutions[1,3]),
                               "aMMR"=c(summary(mcmc_MMR_RR_Brown_trout)$solutions[1,1]),
                               "aMMR_low"=c(summary(mcmc_MMR_RR_Brown_trout)$solutions[1,2]),
                               "aMMR_high"=c(summary(mcmc_MMR_RR_Brown_trout)$solutions[1,3]),
                               "aAS"=c(summary(mcmc_AS_RR_Brown_trout)$solutions[1,1]),
                               "aAS_low"=c(summary(mcmc_AS_RR_Brown_trout)$solutions[1,2]),
                               "aAS_high"=c(summary(mcmc_AS_RR_Brown_trout)$solutions[1,3]),
                              "aFAS"=c(summary(mcmc_FAS_RR_Brown_trout)$solutions[1,1]),
                              "aFAS_low"=c(summary(mcmc_FAS_RR_Brown_trout)$solutions[1,2]),
                              "aFAS_high"=c(summary(mcmc_FAS_RR_Brown_trout)$solutions[1,3]),
                              
                               "aAGR"=c(summary(mcmc_AGR_RR_Brown_trout)$solutions[1,1]),
                               "aAGR_low"=c(summary(mcmc_AGR_RR_Brown_trout)$solutions[1,2]),
                               "aAGR_high"=c(summary(mcmc_AGR_RR_Brown_trout)$solutions[1,3]),
                               "aSGR"=c(summary(mcmc_SGR_RR_Brown_trout)$solutions[1,1]),
                               "aSGR_low"=c(summary(mcmc_SGR_RR_Brown_trout)$solutions[1,2]),
                               "aSGR_high"=c(summary(mcmc_SGR_RR_Brown_trout)$solutions[1,3]) ,
                               "aGE"=c(summary(mcmc_GE_RR_Brown_trout)$solutions[1,1]),
                               "aGE_low"=c(summary(mcmc_GE_RR_Brown_trout)$solutions[1,2]),
                               "aGE_high"=c(summary(mcmc_GE_RR_Brown_trout)$solutions[1,3]))

Cunner_int <- data.frame("Species"=rep("Cunner",1),"Treatment"=c("Mean") ,
                          "aSMR"=c(summary(mcmc_SMR_RR_Cunner)$solutions[1,1]),
                          "aSMR_low"=c(summary(mcmc_SMR_RR_Cunner)$solutions[1,2]),
                          "aSMR_high"=c(summary(mcmc_SMR_RR_Cunner)$solutions[1,3]),
                          "aMMR"=c(summary(mcmc_MMR_RR_Cunner)$solutions[1,1]),
                          "aMMR_low"=c(summary(mcmc_MMR_RR_Cunner)$solutions[1,2]),
                          "aMMR_high"=c(summary(mcmc_MMR_RR_Cunner)$solutions[1,3]),
                          "aAS"=c(summary(mcmc_AS_RR_Cunner)$solutions[1,1]),
                          "aAS_low"=c(summary(mcmc_AS_RR_Cunner)$solutions[1,2]),
                          "aAS_high"=c(summary(mcmc_AS_RR_Cunner)$solutions[1,3]),
                         "aFAS"=c(summary(mcmc_FAS_RR_Cunner)$solutions[1,1]),
                         "aFAS_low"=c(summary(mcmc_FAS_RR_Cunner)$solutions[1,2]),
                         "aFAS_high"=c(summary(mcmc_FAS_RR_Cunner)$solutions[1,3]),
                         
                          "aAGR"=c(summary(mcmc_AGR_RR_Cunner)$solutions[1,1]),
                          "aAGR_low"=c(summary(mcmc_AGR_RR_Cunner)$solutions[1,2]),
                          "aAGR_high"=c(summary(mcmc_AGR_RR_Cunner)$solutions[1,3]),
                          "aSGR"=c(summary(mcmc_SGR_RR_Cunner)$solutions[1,1]),
                          "aSGR_low"=c(summary(mcmc_SGR_RR_Cunner)$solutions[1,2]),
                          "aSGR_high"=c(summary(mcmc_SGR_RR_Cunner)$solutions[1,3]) ,
                          "aGE"=c(summary(mcmc_GE_RR_Cunner)$solutions[1,1]),
                          "aGE_low"=c(summary(mcmc_GE_RR_Cunner)$solutions[1,2]),
                          "aGE_high"=c(summary(mcmc_GE_RR_Cunner)$solutions[1,3]))

Clownfish_int <- data.frame("Species"=rep("Clownfish",1),"Treatment"=c("Mean") ,
                             "aSMR"=c(summary(mcmc_SMR_Clownfish_use)$solutions[1,1]),
                             "aSMR_low"=c(summary(mcmc_SMR_Clownfish_use)$solutions[1,2]),
                             "aSMR_high"=c(summary(mcmc_SMR_Clownfish_use)$solutions[1,3]),
                             "aMMR"=c(summary(mcmc_MMR_Clownfish_use)$solutions[1,1]),
                             "aMMR_low"=c(summary(mcmc_MMR_Clownfish_use)$solutions[1,2]),
                             "aMMR_high"=c(summary(mcmc_MMR_Clownfish_use)$solutions[1,3]),
                             "aAS"=c(summary(mcmc_AS_Clownfish_use)$solutions[1,1]),
                             "aAS_low"=c(summary(mcmc_AS_Clownfish_use)$solutions[1,2]),
                             "aAS_high"=c(summary(mcmc_AS_Clownfish_use)$solutions[1,3]),
                            "aFAS"=c(summary(mcmc_FAS_Clownfish_use)$solutions[1,1]),
                            "aFAS_low"=c(summary(mcmc_FAS_Clownfish_use)$solutions[1,2]),
                            "aFAS_high"=c(summary(mcmc_FAS_Clownfish_use)$solutions[1,3]),
                            
                             "aAGR"=c(summary(mcmc_AGR_Clownfish_use)$solutions[1,1]),
                             "aAGR_low"=c(summary(mcmc_AGR_Clownfish_use)$solutions[1,2]),
                             "aAGR_high"=c(summary(mcmc_AGR_Clownfish_use)$solutions[1,3]),
                             "aSGR"=c(summary(mcmc_SGR_Clownfish_use)$solutions[1,1]),
                             "aSGR_low"=c(summary(mcmc_SGR_Clownfish_use)$solutions[1,2]),
                             "aSGR_high"=c(summary(mcmc_SGR_Clownfish_use)$solutions[1,3]) ,
                             "aGE"=c(summary(mcmc_GE_Clownfish_use)$solutions[1,1]),
                             "aGE_low"=c(summary(mcmc_GE_Clownfish_use)$solutions[1,2]),
                             "aGE_high"=c(summary(mcmc_GE_Clownfish_use)$solutions[1,3]))

Fish_int <- rbind(Zebrafish_int,rbind(Rainbow_trout_int,rbind(Brown_trout_int,rbind(Cunner_int,rbind(Guppy_int,rbind(Clownfish_int, Chromis_int))))))


### Finding between treatment adjusted means###
## Scaling exponent
# rescucturing dataframe 
Results_species2 <- data.frame(matrix(nrow=0,ncol=11))
names(Results_species2) <- c("Species","Treatment","Parameter","Level","bValue","bLow","bHigh","aValue","aLow","aHigh","CV")
para <- c("SMR","MMR","AS","FAS","AGR","SGR","GE")
Level <- c("Onto","Static")
for (x in 1:length(unique(Fish_scal$Species))) {
  for (b in 1:length(unique(Fish_scal[Fish_scal$Species==unique(Fish_scal$Species)[x],]$Treatment))) {
    for(q in 1:2){
      for (y in 1:7) {
        ifelse(x==1&q==1&y==1&b==1,nr<-1,nr<-nr+1)
        Results_species2[nr,1]<- unique(Fish_scal$Species)[x]
        Results_species2[nr,2]<- unique(Fish_scal[Fish_scal$Species==unique(Fish_scal$Species)[x],]$Treatment)[b]
        Results_species2[nr,3]<- para[y]
        Results_species2[nr,4]<- Level[q]
        row <- which(Fish_scal$Species==Results_species2[nr,1]&Fish_scal$Treatment==Results_species2[nr,2])
        Results_species2[nr,5]<- Fish_scal[row,ifelse(q==1,2+y*3-2,23+y*3-2)]
        Results_species2[nr,6]<- Fish_scal[row,ifelse(q==1,3+y*3-2,24+y*3-2)]
        Results_species2[nr,7]<- Fish_scal[row,ifelse(q==1,4+y*3-2,25+y*3-2)]
        Results_species2[nr,8]<- Fish_int[row,2+y*3-2]
        Results_species2[nr,9]<- Fish_int[row,3+y*3-2]
        Results_species2[nr,10]<- Fish_int[row,4+y*3-2]
        Results_species2[nr,11] <- (((Results_species2[nr,7]-Results_species2[nr,6])/4)^2)/Results_species2[nr,5]*100
        
        }  
    }
    }
}

#weigthing messurement
##Weigthing between treatments 
Results_species3 <- Results_species2[!Results_species2$Treatment=="Mean",]
Results_species4 <- data.frame(matrix(nrow=0,ncol=10))
Fish_with_treatments <- c("Zebrafish","Rainbow_trout","Guppy","Chromis")
names(Results_species4) <- c("Species","Parameter","Level","bValue","bLow","bHigh","aValue","aLow","aHigh","CV")
for (x in 1:length(Fish_with_treatments)){
  for (q in 1:2) {
    for (y in 1:7) {
      ifelse(x==1&q==1&y==1,nr<-1,nr<-nr+1)
      Results_species4[nr,1]<- Fish_with_treatments[x]
      Results_species4[nr,2]<- para[y]
      Results_species4[nr,3]<- Level[q]
      Data_ID <- Results_species3[Results_species3$Species==Fish_with_treatments[x]&Results_species3$Parameter==para[y]& Results_species3$Level==Level[q],]
      Data_ID$bSD <- ((Data_ID$bHigh-Data_ID$bLow)/4)^2
      Data_ID$bSDI <-1/Data_ID$bSD
      Data_ID$aSD <- ((Data_ID$aHigh-Data_ID$aLow)/4)^2
      Data_ID$aSDI <-1/Data_ID$aSD
      total_bSDI <- sum(na.omit(Data_ID$bSDI))
      total_aSDI <- sum(na.omit(Data_ID$aSDI))
      Data_ID$bAD_adjust <- Data_ID$bValue*Data_ID$bSDI
      total_bvalue <- sum(na.omit(Data_ID$bAD_adjust))
      total_bSD <- sum((Data_ID$bSDI/total_bSDI)^2*Data_ID$bSD)
      Data_ID$aAD_adjust <- Data_ID$aValue*Data_ID$aSDI
      total_avalue <- sum(na.omit(Data_ID$aAD_adjust))
      total_aSD <- sum((Data_ID$bSDI/total_bSDI)^2*Data_ID$bSD)
      Results_species4[nr,4]<- total_bvalue/total_bSDI
      Results_species4[nr,5] <- Results_species4[nr,4]-2*sqrt(total_bSD)
      Results_species4[nr,6] <- Results_species4[nr,4]+2*sqrt(total_bSD)
      Results_species4[nr,7]<- total_avalue/total_aSDI
      Results_species4[nr,8] <- Results_species4[nr,7]-2*sqrt(total_aSD)
      Results_species4[nr,9] <- Results_species4[nr,7]+2*sqrt(total_aSD)
      Results_species4[nr,10] <- sqrt(total_bSD)/Results_species4[nr,4]*100
    }
  }
}




Scaling_exponents_adj <- rbind(Results_species4,Results_species2[Results_species2$Species %in% c("Clownfish","Cunner","Brown_trout"),c(1,3,4,5,6,7,8,9,10,11)])
Evo_scal <- data.frame("Species"=rep("Mean",7),"Parameter"=c("SMR","MMR","AS","FAS","AGR","SGR","GE"),"Level"=rep("Evo",7),
                       "bValue"=c(summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[2,1],summary(mcmc_MMR_RR_Evol_TempAdj)$solutions[2,1],summary(mcmc_AS_RR_Evol_TempAdj)$solutions[2,1],summary(mcmc_FAS_RR_Evol_TempAdj)$solutions[2,1],summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[2,1],summary(mcmc_SGR_RR_Evol_TempAdj)$solutions[2,1],summary(mcmc_GE_RR_Evol_TempAdj)$solutions[2,1]),
                       "bLow"=c(summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[2,2],summary(mcmc_MMR_RR_Evol_TempAdj)$solutions[2,2],summary(mcmc_AS_RR_Evol_TempAdj)$solutions[2,2],summary(mcmc_FAS_RR_Evol_TempAdj)$solutions[2,2],summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[2,2],summary(mcmc_SGR_RR_Evol_TempAdj)$solutions[2,2],summary(mcmc_GE_RR_Evol_TempAdj)$solutions[2,2]),
                       "bHigh"=c(summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[2,3],summary(mcmc_MMR_RR_Evol_TempAdj)$solutions[2,3],summary(mcmc_AS_RR_Evol_TempAdj)$solutions[2,3],summary(mcmc_FAS_RR_Evol_TempAdj)$solutions[2,3],summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[2,3],summary(mcmc_SGR_RR_Evol_TempAdj)$solutions[2,3],summary(mcmc_GE_RR_Evol_TempAdj)$solutions[2,3]),
                       "aValue"=c(summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[1,1],summary(mcmc_MMR_RR_Evol_TempAdj)$solutions[1,1],summary(mcmc_AS_RR_Evol_TempAdj)$solutions[1,1],summary(mcmc_FAS_RR_Evol_TempAdj)$solutions[1,1],summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[1,1],summary(mcmc_SGR_RR_Evol_TempAdj)$solutions[1,1],summary(mcmc_GE_RR_Evol_TempAdj)$solutions[1,1]),
                       "aLow"=c(summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[1,2],summary(mcmc_MMR_RR_Evol_TempAdj)$solutions[1,2],summary(mcmc_AS_RR_Evol_TempAdj)$solutions[1,2],summary(mcmc_FAS_RR_Evol_TempAdj)$solutions[1,2],summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[1,2],summary(mcmc_SGR_RR_Evol_TempAdj)$solutions[1,2],summary(mcmc_GE_RR_Evol_TempAdj)$solutions[1,2]),
                       "aHigh"=c(summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[1,3],summary(mcmc_MMR_RR_Evol_TempAdj)$solutions[1,3],summary(mcmc_AS_RR_Evol_TempAdj)$solutions[1,3],summary(mcmc_FAS_RR_Evol_TempAdj)$solutions[1,3],summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[1,3],summary(mcmc_SGR_RR_Evol_TempAdj)$solutions[1,3],summary(mcmc_GE_RR_Evol_TempAdj)$solutions[1,3]))
Evo_scal$CV <- (((Evo_scal$bHigh-Evo_scal$bLow)/4)^2)/Evo_scal$bValue*100
Mean_onto_static <- data.frame(matrix(nrow=0,ncol=10))
names(Mean_onto_static) <- c("Species","Parameter","Level","bValue","bLow","bHigh","aValue","aLow","aHigh","CV")

for (q in 1:2) {
    for (y in 1:7) {
      ifelse(q==1&y==1,nr<-1,nr<-nr+1)
      Mean_onto_static[nr,1]<- "Mean"
      Mean_onto_static[nr,2]<- para[y]
      Mean_onto_static[nr,3]<- Level[q]
      Data_ID <- Scaling_exponents_adj[Scaling_exponents_adj$Parameter==para[y]& Scaling_exponents_adj$Level==Level[q],]
      Data_ID$bSD <- ((Data_ID$bHigh-Data_ID$bLow)/4)^2
      Data_ID$bSDI <-1/Data_ID$bSD
      Data_ID$aSD <- ((Data_ID$aHigh-Data_ID$aLow)/4)^2
      Data_ID$aSDI <-1/Data_ID$aSD
      total_bSDI <- sum(na.omit(Data_ID$bSDI))
      total_aSDI <- sum(na.omit(Data_ID$aSDI))
      Data_ID$bAD_adjust <- Data_ID$bValue*Data_ID$bSDI
      total_bvalue <- sum(na.omit(Data_ID$bAD_adjust))
      total_bSD <- sum(na.omit((Data_ID$bSDI/total_bSDI)^2*Data_ID$bSD))
      Data_ID$aAD_adjust <- Data_ID$aValue*Data_ID$aSDI
      total_avalue <- sum(na.omit(Data_ID$aAD_adjust))
      total_aSD <- sum(na.omit((Data_ID$bSDI/total_bSDI)^2*Data_ID$bSD))
      Mean_onto_static[nr,4]<- total_bvalue/total_bSDI
      Mean_onto_static[nr,5] <- Mean_onto_static[nr,4]-2*sqrt(total_bSD)
      Mean_onto_static[nr,6] <- Mean_onto_static[nr,4]+2*sqrt(total_bSD)
      Mean_onto_static[nr,7]<- total_avalue/total_aSDI
      Mean_onto_static[nr,8] <- Mean_onto_static[nr,7]-2*sqrt(total_aSD)
      Mean_onto_static[nr,9] <- Mean_onto_static[nr,7]+2*sqrt(total_aSD)
      Mean_onto_static[nr,10] <- sqrt(total_bSD)/Mean_onto_static[nr,4]*100
    }
}
Scaling_exponents_all <- rbind(Scaling_exponents_adj,rbind(Mean_onto_static,Evo_scal))
Scaling_exponents_all$CV <- format(Scaling_exponents_all$CV,scientific=F)
Scaling_exponents_all$Parameter <- factor(Scaling_exponents_all$Parameter,levels=c("SMR","MMR","AS","FAS","AGR","SGR","GE"))
Results_fit[Results_fit$Species=="Brown_trout",]$Slope_SMR_fit


Dif_onto_statik <- data.frame("Species"=c(1),"Parameter"=c(1),"Comparason"=c(1),"Value"=c(1),"Significans"=c(1))
para <- c("SMR","MMR","AS","FAS","AGR","SGR","GE")
Levels <- c("Onto","Static","Evo")
for (x in 1:8) {
  for (y in 1:6) {
    Data_ID <- Scaling_exponents_all[(Scaling_exponents_all$Parameter==para[y]&Scaling_exponents_all$Species==unique(Scaling_exponents_all$Species)[x])|(Scaling_exponents_all$Level=="Evo"&Scaling_exponents_all$Parameter==para[y]),]
    for (q in 1:3) {
      ifelse(q==1&y==1&x==1,nr<-1,nr<-nr+1)
      Dif_onto_statik[nr,1]<- unique(Scaling_exponents_all$Species)[x]
      Dif_onto_statik[nr,2]<- para[y]
      ifelse(q==1,Dif_onto_statik[nr,3]<-"Onto vs Static",ifelse(q==2,Dif_onto_statik[nr,3]<-"Onto vs Evo",Dif_onto_statik[nr,3]<-"Static vs Evo"))
      ifelse(q==1,Dif_onto_statik[nr,4]<-Data_ID$bValue[1]-Data_ID$bValue[2],ifelse(q==2,Dif_onto_statik[nr,4]<-Data_ID$bValue[1]-Data_ID$bValue[3],Dif_onto_statik[nr,4]<-Data_ID$bValue[2]-Data_ID$bValue[3]))
      ifelse(q==1,ifelse((Data_ID$bLow[1]>Data_ID$bHigh[2]|Data_ID$bLow[2]>Data_ID$bHigh[1]),Dif_onto_statik[nr,5]<-"yes",Dif_onto_statik[nr,5]<-"No"),
             ifelse(q==2,ifelse((Data_ID$bLow[1]>Data_ID$bHigh[3]|Data_ID$bLow[3]>Data_ID$bHigh[1]),Dif_onto_statik[nr,5]<-"yes",Dif_onto_statik[nr,5]<-"No"),
                    ifelse((Data_ID$bLow[2]>Data_ID$bHigh[3]|Data_ID$bLow[3]>Data_ID$bHigh[2]),Dif_onto_statik[nr,5]<-"yes",Dif_onto_statik[nr,5]<-"No")))
      
    }
  }
}

# Static scaling 
Data_fit$SMR_staticMT <- NA
Data_fit$MMR_staticMT <- NA
Data_fit$AS_staticMT <- NA
Data_fit$FAS_staticMT <- NA
Data_fit$AGR_staticMT <- NA
Data_fit$SGR_staticMT <- NA
Data_fit$GE_staticMT <- NA

Z_SMR_f <- summary(mcmc_SMR_Zebrafish_use)$solutions[6,1]
R_SMR_f <- summary(mcmc_SMR_Rainbow_trout_use)$solutions[6,1]
R_MMR_f <- summary(mcmc_MMR_Rainbow_trout_use)$solutions[6,1]
R_AS_f <- summary(mcmc_AS_Rainbow_trout_use)$solutions[6,1]
R_FAS_f <- summary(mcmc_FAS_Rainbow_trout_use)$solutions[6,1]
T_SMR <- summary(mcmc_SMR_RR_Evol)$solutions[3,1]
T_MMR <- summary(mcmc_MMR_RR_Evol)$solutions[3,1]
T_AS <- summary(mcmc_AS_RR_Evol)$solutions[3,1]
T_FAS <- summary(mcmc_FAS_RR_Evol)$solutions[3,1]
T_AGR <- summary(mcmc_AGR_RR_Evol)$solutions[3,1]
T_SGR <- summary(mcmc_SGR_RR_Evol)$solutions[3,1]
T_GE <- summary(mcmc_GE_RR_Evol)$solutions[3,1]



for (x in 1:length(Data_fit$Species)) {

  Species <- Data_fit$Species[x]
  Data_ID <- Scaling_exponents_all[Scaling_exponents_all$Species==Species & Scaling_exponents_all$Level=="Static",]
  if(Species=="Zebrafish"){
    Data_fit$SMR_staticMT[x] <- 10^(Data_fit$logMass[x]*Data_ID[Data_ID$Parameter=="SMR",4]+Data_ID[Data_ID$Parameter=="SMR",7]+Z_SMR_f*Data_fit$Temp[x])* ((10^(T_SMR * 10))^((20-Data_fit$Temp[x])/10))
    Data_fit$MMR_staticMT[x] <- NA
    Data_fit$AS_staticMT[x] <- NA
    Data_fit$FAS_staticMT[x] <- NA
    Data_fit$AGR_staticMT[x] <- 10^(Data_fit$logMassGrowth[x]*Data_ID[Data_ID$Parameter=="AGR",4]+Data_ID[Data_ID$Parameter=="AGR",7])* (10^(T_AGR * 10))^((20-Data_fit$Temp[x])/10)
    Data_fit$SGR_staticMT[x] <- 10^(Data_fit$logMassGrowth[x]*Data_ID[Data_ID$Parameter=="SGR",4]+Data_ID[Data_ID$Parameter=="SGR",7])* (10^(T_SGR * 10))^((20-Data_fit$Temp[x])/10)
    Data_fit$GE_staticMT[x] <- 10^(Data_fit$logMassGrowth[x]*Data_ID[Data_ID$Parameter=="GE",4]+Data_ID[Data_ID$Parameter=="GE",7])* (10^(T_GE * 10))^((20-Data_fit$Temp[x])/10)}
  if(Species=="Rainbow_trout"){
    Data_fit$SMR_staticMT[x] <- 10^(Data_fit$logMass[x]*Data_ID[Data_ID$Parameter=="SMR",4]+Data_ID[Data_ID$Parameter=="SMR",7]+R_SMR_f*Data_fit$Temp[x])* (10^(T_SMR * 10))^((20-Data_fit$Temp[x])/10)
    Data_fit$MMR_staticMT[x] <- 10^(Data_fit$logMass[x]*Data_ID[Data_ID$Parameter=="MMR",4]+Data_ID[Data_ID$Parameter=="MMR",7]+R_MMR_f*Data_fit$Temp[x])* (10^(T_MMR * 10))^((20-Data_fit$Temp[x])/10)
    Data_fit$AS_staticMT[x] <- 10^(Data_fit$logMass[x]*Data_ID[Data_ID$Parameter=="AS",4]+Data_ID[Data_ID$Parameter=="AS",7]+R_AS_f*Data_fit$Temp[x])* (10^(T_AS * 10))^((20-Data_fit$Temp[x])/10)
    Data_fit$FAS_staticMT[x] <- 10^(Data_fit$logMass[x]*Data_ID[Data_ID$Parameter=="FAS",4]+Data_ID[Data_ID$Parameter=="FAS",7]+R_FAS_f*Data_fit$Temp[x])* (10^(T_FAS * 10))^((20-Data_fit$Temp[x])/10)
    Data_fit$AGR_staticMT[x] <- 10^(Data_fit$logMassGrowth[x]*Data_ID[Data_ID$Parameter=="AGR",4]+Data_ID[Data_ID$Parameter=="AGR",7])* (10^(T_AGR * 10))^((20-Data_fit$Temp[x])/10)
    Data_fit$SGR_staticMT[x] <- 10^(Data_fit$logMassGrowth[x]*Data_ID[Data_ID$Parameter=="SGR",4]+Data_ID[Data_ID$Parameter=="SGR",7])* (10^(T_SGR * 10))^((20-Data_fit$Temp[x])/10)
    Data_fit$GE_staticMT[x] <- 10^(Data_fit$logMassGrowth[x]*Data_ID[Data_ID$Parameter=="GE",4]+Data_ID[Data_ID$Parameter=="GE",7])* (10^(T_GE * 10))^((20-Data_fit$Temp[x])/10)}
  if(!Species=="Rainbow_trout"&!Species=="Zebrafish"){
    Data_fit$SMR_staticMT[x] <- 10^(Data_fit$logMass[x]*Data_ID[Data_ID$Parameter=="SMR",4]+Data_ID[Data_ID$Parameter=="SMR",7])* ((10^(T_SMR * 10))^((20-Data_fit$Temp[x])/10))
    Data_fit$MMR_staticMT[x] <- 10^(Data_fit$logMass[x]*Data_ID[Data_ID$Parameter=="MMR",4]+Data_ID[Data_ID$Parameter=="MMR",7])* ((10^(T_MMR * 10))^((20-Data_fit$Temp[x])/10))
    Data_fit$AS_staticMT[x] <- 10^(Data_fit$logMass[x]*Data_ID[Data_ID$Parameter=="AS",4]+Data_ID[Data_ID$Parameter=="AS",7])* ((10^(T_AS * 10))^((20-Data_fit$Temp[x])/10))
    Data_fit$FAS_staticMT[x] <- 10^(Data_fit$logMass[x]*Data_ID[Data_ID$Parameter=="FAS",4]+Data_ID[Data_ID$Parameter=="FAS",7])* ((10^(T_FAS * 10))^((20-Data_fit$Temp[x])/10))
    Data_fit$AGR_staticMT[x] <- 10^(Data_fit$logMassGrowth[x]*Data_ID[Data_ID$Parameter=="AGR",4]+Data_ID[Data_ID$Parameter=="AGR",7])* ((10^(T_AGR * 10))^((20-Data_fit$Temp[x])/10))
    Data_fit$SGR_staticMT[x] <- 10^(Data_fit$logMassGrowth[x]*Data_ID[Data_ID$Parameter=="SGR",4]+Data_ID[Data_ID$Parameter=="SGR",7])* ((10^(T_SGR * 10))^((20-Data_fit$Temp[x])/10))
    Data_fit$GE_staticMT[x] <- 10^(Data_fit$logMassGrowth[x]*Data_ID[Data_ID$Parameter=="GE",4]+Data_ID[Data_ID$Parameter=="GE",7])* ((10^(T_GE * 10))^((20-Data_fit$Temp[x])/10))}
}


#### Quick corelations ####

#Statestik
bcor_table <- data.frame(matrix(nrow=0,ncol=9))
names(bcor_table)<-c("Par_1","Par_2","Cor","p_val","level", "Compereson","Col","top","but")
for (q in 1:2) {
  for (x in 1:(length(unique(Scaling_exponents_all$Parameter) )-1)) {
    for (y in (x+1):length(unique(Scaling_exponents_all$Parameter))) {
      ifelse(x==1&y==2&q==1,nr <-1,nr<-nr+1)
      Para_1 <- unique(Scaling_exponents_all$Parameter)[x]
      Para_2 <- unique(Scaling_exponents_all$Parameter)[y]
      ifelse(q==2,Data_1_v <- Scaling_exponents_all[Scaling_exponents_all$Level=="Static"& !Scaling_exponents_all$Species=="Mean" &Scaling_exponents_all$Parameter==Para_1,]$bValue,Data_1_v <-Results_fit_p[Results_fit_p$Parameter==Para_1,]$Slope_fit)
      ifelse(q==2,Data_2_v <- Scaling_exponents_all[Scaling_exponents_all$Level=="Static"& !Scaling_exponents_all$Species=="Mean" &Scaling_exponents_all$Parameter==Para_2,]$bValue,Data_2_v <-Results_fit_p[Results_fit_p$Parameter==Para_2,]$Slope_fit)
      Cor_ID <- cor.test(Data_1_v,Data_2_v)
      bcor_table[nr,1] <- Para_1
      bcor_table[nr,2] <- Para_2
      bcor_table[nr,3] <- Cor_ID$estimate
      bcor_table[nr,4] <- round(Cor_ID$p.value,3)
      bcor_table[nr,5] <- unique(Scaling_exponents_all$Level)[q]
      bcor_table[nr,6] <- paste(unique(Scaling_exponents_all$Level)[q]," ",Para_1," vs ",Para_2," p = ",round(Cor_ID$p.value,3),sep = "")
      ifelse(bcor_table[nr,4]>0.05,bcor_table[nr,7] <-"red",bcor_table[nr,7] <-"green")
      bcor_table[nr,8]<- max(Data_2_v)
      bcor_table[nr,9]<- max(Data_2_v)*1.1
      
    }
  }
}

bcor_table1 <- data.frame(matrix(nrow=0,ncol=5))
names(bcor_table1)<-c("Par_1","Par_2","Cor","p_val", "Compereson")

  for (x in 1:(length(unique(Scaling_exponents_all[!Scaling_exponents_all$Parameter%in%c("SGR","AS"),]$Parameter) )-1)) {
    for (y in (x+1):length(unique(Scaling_exponents_all[!Scaling_exponents_all$Parameter%in%c("SGR","AS"),]$Parameter))) {
      ifelse(x==1&y==2,nr <-1,nr<-nr+1)
      Para_1 <- unique(Scaling_exponents_all[!Scaling_exponents_all$Parameter%in%c("SGR","AS"),]$Parameter)[x]
      Para_2 <- unique(Scaling_exponents_all[!Scaling_exponents_all$Parameter%in%c("SGR","AS"),]$Parameter)[y]
      Data_ID <-data.frame("Species"=Results_fit_p[Results_fit_p$Parameter==Para_1&Results_fit$Species=="Brown_trout",]$Species,"Para_1"=Results_fit_p[Results_fit_p$Parameter==Para_1&Results_fit$Species=="Brown_trout",]$Slope_fit,"Para_2"=Results_fit_p[Results_fit_p$Parameter==Para_2&Results_fit$Species=="Brown_trout",]$Slope_fit) 
      
      reg <- summary(lmer(data = na.omit(Data_ID),Para_1~Para_2+(1|Species)))
      
      bcor_table1[nr,1] <- Para_1
      bcor_table1[nr,2] <- Para_2
      bcor_table1[nr,3] <- round(reg$coefficients[2,1],2)
      bcor_table1[nr,4] <- round(reg$coefficients[2,5],3)
      
      bcor_table1[nr,5] <- paste(Para_1," vs ",Para_2," p = ",round(reg$coefficients[2,5],3),sep = "")

      
    }
  }


## tranfering data to multible plotting graf
Plot_scaling_variblels <- Scaling_exponents_all[!Scaling_exponents_all$Parameter%in%c("SGR","AS")&!Scaling_exponents_all$Species%in%c("Mean","Brown_trout")& Scaling_exponents_all$Level=="Onto",]
for (q in 1:2) {
  for (x in 1:(length(unique(Plot_scaling_variblels$Parameter) )-1)) {
    for (y in (x+1):length(unique(Plot_scaling_variblels$Parameter))) {
      Para_1 <- unique(Plot_scaling_variblels$Parameter)[x]
      Para_2 <- unique(Plot_scaling_variblels$Parameter)[y]
      ifelse(q==2,Data_1 <- Plot_scaling_variblels[Plot_scaling_variblels$Level=="Static"& Plot_scaling_variblels$Parameter==Para_1,],Data_1 <-Results_fit_p[Results_fit_p$Parameter==Para_1,])
      ifelse(q==2,Data_2 <- Plot_scaling_variblels[Plot_scaling_variblels$Level=="Static"& Plot_scaling_variblels$Parameter==Para_2,],Data_2 <-Results_fit_p[Results_fit_p$Parameter==Para_2,])
      ifelse(q==1,Cor_ID <- cor.test(Data_1$Slope_fit,Data_2$Slope_fit),Cor_ID <- cor.test(Data_1$bValue,Data_2$bValue))
      
      bCor_plot_table_ID <- data.frame(matrix(nrow=ifelse(q==1,length(Results_fit_p[Results_fit_p$Parameter==Para_1,]$Species),7),ncol=9))
      names(bCor_plot_table_ID)<-c("Par_1","Par_2","Level","Compereson", "Par_1_low","Par_1_high","Par_2_low","Par_2_high","Species")
      ifelse(q==1,bCor_plot_table_ID[,1] <- Data_1$Slope_fit,bCor_plot_table_ID[,1] <-Data_1$bValue)
      ifelse(q==1,bCor_plot_table_ID[,2] <- Data_2$Slope_fit,bCor_plot_table_ID[,2] <-Data_2$bValue)
      bCor_plot_table_ID[,3] <- unique(Plot_scaling_variblels$Level)[q]
      bCor_plot_table_ID[,4] <- paste(unique(Plot_scaling_variblels$Level)[q]," ","b",Para_1," vs ","b",Para_2," p ", ifelse(round(Cor_ID$p.value,3)>0.000,paste("= ",round(Cor_ID$p.value,3),sep=""),"< 0.001")," cor = ",round(Cor_ID$estimate,3),sep = "")
      ifelse(q==1,bCor_plot_table_ID[,5] <- NA,bCor_plot_table_ID[,5] <- Data_1$bLow)
      ifelse(q==1,bCor_plot_table_ID[,6] <- NA,bCor_plot_table_ID[,6] <- Data_1$bHigh)
      ifelse(q==1,bCor_plot_table_ID[,7] <- NA,bCor_plot_table_ID[,7] <- Data_2$bLow)
      ifelse(q==1,bCor_plot_table_ID[,8] <- NA,bCor_plot_table_ID[,8] <- Data_2$bHigh)
      bCor_plot_table_ID[,9]<- Data_1$Species
      ifelse(x==1&y==2&q==1,bCor_plot_table <-bCor_plot_table_ID,bCor_plot_table<-rbind(bCor_plot_table,bCor_plot_table_ID))
      }
  }
}

All_cor_plot <- ggplot(data = bCor_plot_table)+
  geom_point(aes(y=Par_1,x=Par_2,col=Species))+
  geom_line(aes(y = Par_1, x = Par_2),col="black", stat = "smooth", alpha = 1, size = 0.8, method = "lm", se=F, linetype="dotted",show.legend = FALSE)+
  scale_colour_manual(name = NULL, values = c( "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF"))+
  facet_wrap(~Compereson,scales = "free")






#SMR~AGR
cor.test(Results_fit$Slope_SMR_fit,Results_fit$Slope_MMR_fit)

lm(data=Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_SMR_fit~Slope_AGR_fit )
summary(lm(data=Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_SMR_fit~Slope_SGR_fit ))
summary(lm(data=Results_fit[!is.na(Results_fit$Slope_AGR_fit)& Results_fit$Species=="Zebrafish",],Slope_SMR_fit~Slope_SGR_fit ))



cor.test(Results_fit[Results_fit$Species=="Brown_trout",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Brown_trout",]$Slope_AGR_fit)
cor.test(Results_fit[Results_fit$Species=="Chromis",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Chromis",]$Slope_AGR_fit)
cor.test(Results_fit[Results_fit$Species=="Clownfish",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Clownfish",]$Slope_AGR_fit)
cor.test(Results_fit[Results_fit$Species=="Cunner",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Cunner",]$Slope_AGR_fit)
cor.test(Results_fit[Results_fit$Species=="Guppy",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Guppy",]$Slope_AGR_fit)
cor.test(Results_fit[Results_fit$Species=="Rainbow_trout",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Rainbow_trout",]$Slope_AGR_fit)
cor.test(Results_fit[Results_fit$Species=="Zebrafish",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Zebrafish",]$Slope_AGR_fit)

#SMR~MMR
cor.test(Results_fit$Slope_SMR_fit,Results_fit$Slope_MMR_fit)

cor.test(Results_fit[Results_fit$Species=="Brown_trout",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Brown_trout",]$Slope_MMR_fit)
cor.test(Results_fit[Results_fit$Species=="Chromis",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Chromis",]$Slope_MMR_fit)
cor.test(Results_fit[Results_fit$Species=="Clownfish",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Clownfish",]$Slope_MMR_fit)
cor.test(Results_fit[Results_fit$Species=="Cunner",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Cunner",]$Slope_MMR_fit)
cor.test(Results_fit[Results_fit$Species=="Guppy",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Guppy",]$Slope_MMR_fit)
cor.test(Results_fit[Results_fit$Species=="Rainbow_trout",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Rainbow_trout",]$Slope_MMR_fit)
cor.test(Results_fit[Results_fit$Species=="Zebrafish",]$Slope_SMR_fit,Results_fit[Results_fit$Species=="Zebrafish",]$Slope_MMR_fit)


bSMR_bAGR_data <- Results_scaling[Results_scaling$Treatment=="Mean"&(Results_scaling$Parameter=="SMR"|Results_scaling$Parameter=="AGR"),]



bSMR_bAGR_data_evo <- data.frame("Parameter"=c("SMR","AGR"),"Level"=c("Evolutionary","Evolutionary"),"Value"=c(summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[2,1],summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[2,1]),
                                 "Low"=c(summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[2,2],summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[2,2]),"High"=c(summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[2,3],summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[2,3]))
bSMR_bAGR_data_mean <- data.frame("Value"=c(mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="SMR"& bSMR_bAGR_data$Level=="Onto",]$Value),mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="AGR"& bSMR_bAGR_data$Level=="Onto",]$Value),mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="SMR"& bSMR_bAGR_data$Level=="Static",]$Value),mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="AGR"& bSMR_bAGR_data$Level=="Static",]$Value)),
                                  "Level"=c("Onto","Onto","Static","Static"),"Parameter"=c("SMR","AGR","SMR","AGR"),
                                  "Low"=c(mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="SMR"& bSMR_bAGR_data$Level=="Onto",]$Low),mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="AGR"& bSMR_bAGR_data$Level=="Onto",]$Low),mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="SMR"& bSMR_bAGR_data$Level=="Static",]$Low),mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="AGR"& bSMR_bAGR_data$Level=="Static",]$Low)),
                                  "High"=c(mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="SMR"& bSMR_bAGR_data$Level=="Onto",]$High),mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="AGR"& bSMR_bAGR_data$Level=="Onto",]$High),mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="SMR"& bSMR_bAGR_data$Level=="Static",]$High),mean(bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="AGR"& bSMR_bAGR_data$Level=="Static",]$High)))
bSMR_bAGR_data_evo$Species <- "Mean"
bSMR_bAGR_data[bSMR_bAGR_data$Parameter=="SMR",]

bSMR_bAGR_data_mean_adj <- data.frame(matrix(nrow = 0,ncol = 5))
names(bSMR_bAGR_data_mean_adj) <- c("Parameter","Level","Value","Low","High")
for (x in 1:length(unique(bSMR_bAGR_data$Parameter))) {
  Data_ID <- bSMR_bAGR_data[bSMR_bAGR_data$Parameter==unique(bSMR_bAGR_data$Parameter)[x],]
  for (q in 1:length(unique(Data_ID$Level))) {
    Data_ID_q <- Data_ID[Data_ID$Level==unique(Data_ID$Level)[q],]
    ifelse(x==1&q==1, nr <-1,nr<-1+nr)
    bSMR_bAGR_data_mean_adj[nr,1]<- unique(bSMR_bAGR_data$Parameter)[x]
    bSMR_bAGR_data_mean_adj[nr,2]<- unique(Data_ID$Level)[q]
    total_SDI <- sum(na.omit(Data_ID_q$SDI))
    Data_ID_q$AD_adjust <- Data_ID_q$SDI*Data_ID_q$Value
    total_value <- sum(na.omit(Data_ID_q$AD_adjust))
    bSMR_bAGR_data_mean_adj[nr,3]<- total_value/total_SDI
    total_SD <- sum((Data_ID_q$W)^2*Data_ID_q$SD)
    bSMR_bAGR_data_mean_adj[nr,5] <- bSMR_bAGR_data_mean_adj[nr,3]+2*sqrt(total_SD)
    bSMR_bAGR_data_mean_adj[nr,4] <- bSMR_bAGR_data_mean_adj[nr,3]-2*sqrt(total_SD)
  }
}

  

#####Figure 1

nA <- grobTree(textGrob("A", x=0.03,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nB <- grobTree(textGrob("B", x=0.03,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nC <- grobTree(textGrob("C", x=0.03,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nD <- grobTree(textGrob("D", x=0.03,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nE <- grobTree(textGrob("E", x=0.03,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nF <- grobTree(textGrob("F", x=0.03,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))

nSMR <- grobTree(textGrob("Standard metabolic rate", x=0.1,  y=0.94, hjust=0,gp=gpar(col="black", fontsize=8)))
nMMR <- grobTree(textGrob("Maximum metabolic rate", x=0.1,  y=0.94, hjust=0,gp=gpar(col="black", fontsize=8)))
nAS <- grobTree(textGrob("Aerobic scope", x=0.1,  y=0.94, hjust=0,gp=gpar(col="black", fontsize=8)))
nFAS <- grobTree(textGrob("Factorial aerobic scope", x=0.13,  y=0.94, hjust=0,gp=gpar(col="black", fontsize=8)))
nGE <- grobTree(textGrob("Growth efficiency", x=0.1,  y=0.94, hjust=0,gp=gpar(col="black", fontsize=8)))
nAGR <- grobTree(textGrob("Abselute growth rate", x=0.1,  y=0.94, hjust=0,gp=gpar(col="black", fontsize=8)))

nSMRu <- grobTree(textGrob(expression((mg~O[2]~h^-1)), x=0.1,  y=0.87, hjust=0,gp=gpar(col="black", fontsize=8)))
nMMRu <- grobTree(textGrob(expression((mg~O[2]~h^-1)), x=0.1,  y=0.87, hjust=0,gp=gpar(col="black", fontsize=8)))
nASu <- grobTree(textGrob(expression((mg~O[2]~h^-1)), x=0.1,  y=0.87, hjust=0,gp=gpar(col="black", fontsize=8)))
nGEu <- grobTree(textGrob(expression((mg~mg~ O[2]^-1)), x=0.1,  y=0.87, hjust=0,gp=gpar(col="black", fontsize=8)))
nAGRu <- grobTree(textGrob(expression((mg ~ day^{-1})), x=0.1,  y=0.87, hjust=0,gp=gpar(col="black", fontsize=8)))


Data_fit$Species1 <- paste(Data_fit$Species,"1",sep = "")


Legend <- ggplot(data = Results_fit, aes(x=Slope_AGR_fit,y=Slope_SMR_fit))+
  geom_point(aes(col=Species))+
  geom_line(aes(col="Mean"),stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F )+
  geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  #geom_line(aes(col=Species), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  #geom_line(aes(col=Species, group=Plot_group), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  #theme(legend.position="none")+
  theme(legend.position="bottom",legend.key.size = unit(7,"mm"),legend.text = element_text(size=10,margin = margin(r = 3,l=0)),
        legend.box.background = element_rect(colour = "black"),legend.spacing.x = unit(0.0, 'mm'))+
  force_panelsizes(rows = unit(3, "in"), cols = unit(3, "in"))+
  theme(axis.ticks.length=unit(2,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=14, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=14,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(size=16, margin = margin(t = 8, r = 0, b = 0, l = 0))) +
  #scale_colour_manual(values = c("red3", "dodgerblue2")) +
  
  labs(y=expression(Scaling ~ exponent ~ (italic(b)) ~ "for" ~ SMR), x=expression(Scaling ~ exponent ~ (italic(b)) ~ "for" ~ AGR))+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,1,1,1, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "#00B6EB", "Brown_trout" = "black",  "Chromis" = "#DA1F26", "Clownfish" = "#F26739" , "Cunner" = "#00BE00","Guppy"="#FFBE00","Rainbow_trout"="#0000FF","Zebrafish"="black","ref"="#B400FF"),
                      guide = guide_legend(title = "", nrow = 1,override.aes = list(linetype = c("solid",rep("blank", 7),"dotted"),size=2,shape = c(NA,rep(16, 7),NA))), 
                      labels =c("Overall", "Brown trout","Damselfish","Clownfish","Cunner","Guppy","Rainbow trout","Zebrafish", "1:1 Reference line")) 




SMR_plot <- ggplot(Data_fit, aes(x = Mass, y = 10^SMR_fit_TempAdj, group = Group)) +
  geom_abline(aes(slope=0.75, intercept=-0.75, colour="ref"), linetype="dotted") +
  geom_abline(aes(slope=summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[2,1], intercept=summary(mcmc_SMR_RR_Evol_TempAdj)$solutions[1,1],col="Mean"), size=0.8) +
      #geom_point(aes(x = Mass, y = SMR_TempAdj), shape = 16, colour = "grey90", alpha = 1, size = 0.6) +
  #geom_line(aes(x = Mass, y = SMR_TempAdj), colour = "grey90", size = 0.2, alpha = 1) +
  geom_line(data=Data_fit, aes(x = Mass, y = SMR_fit_TempAdj, group = Group, colour=Species), stat = "smooth", alpha = 0.3, size = 0.3, method = "lm", se=F) +
  geom_line(data=Data_fit, aes(x = Mass, y = SMR_staticMT, group = Species, colour=Species1), stat = "smooth", alpha = 1, size = 0.8, method = "lm", se=F, linetype="dotted", show.legend = FALSE) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  annotation_custom(nA)+ annotation_custom(nSMR)+annotation_custom(nSMRu)+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1.0,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_blank())+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  #scale_colour_manual(values = c("red3", "dodgerblue2")) +
  scale_x_log10(labels = label_number(c(0.01,0.1,1,10,100)),breaks=c(0.01,0.1,1,10,100), limits=c(0.0009,68.2))+
  scale_y_log10(labels = label_number(c(0.001,0.01,0.1,1,10)),breaks=c(0.001,0.01,0.1,1,10))+
  labs( x=expression(Body ~ mass ~ (g)))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF",
                                              "ref"="grey24", "Brown_trout1" = "#b91a20",  "Chromis1" = "#ee4710", "Clownfish1" = "#00a100" , "Cunner1" = "#d9a100","Guppy1"="#00a9d9","Rainbow_trout1"="#0000d9","Zebrafish1"="#9900d9")) +
  theme(plot.title = element_text(hjust = 0,vjust = -1.7,size = 12))

MMR_plot <- ggplot(Data_fit, aes(x = Mass, y = 10^MMR_fit_TempAdj, group = Group)) +
  geom_abline(aes(slope=0.75, intercept=-0.1, colour="ref"), linetype="dotted") +
  geom_abline(aes(slope=summary(mcmc_MMR_RR_Evol_TempAdj)$solutions[2,1], intercept=summary(mcmc_MMR_RR_Evol_TempAdj)$solutions[1,1],col="Mean"), size=0.8) +
  #geom_point(aes(x = Mass, y = SMR_TempAdj), shape = 16, colour = "grey90", alpha = 1, size = 0.6) +
  #geom_line(aes(x = Mass, y = SMR_TempAdj), colour = "grey90", size = 0.2, alpha = 1) +
  geom_line(data=Data_fit, aes(x = Mass, y = MMR_fit_TempAdj, group = Group, colour=Species), stat = "smooth", alpha = 0.3, size = 0.3, method = "lm", se=F) +
  geom_line(data=Data_fit, aes(x = Mass, y = MMR_staticMT, group = Species, colour=Species1), stat = "smooth", alpha = 1, size = 0.8, method = "lm", se=F, linetype="dotted", show.legend = FALSE) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  annotation_custom(nB)+annotation_custom(nMMR)+annotation_custom(nMMRu)+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1.0,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_blank())+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  #scale_colour_manual(values = c("red3", "dodgerblue2")) +
  scale_x_log10(labels = label_number(c(0.01,0.1,1,10,100)),breaks=c(0.01,0.1,1,10,100), limits=c(0.0009,68.2))+
  scale_y_log10(labels = label_number(c(0.001,0.01,0.1,1,10)),breaks=c(0.001,0.01,0.1,1,10))+
  labs( x=expression(Body ~ mass ~ (g)))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF",
                                              "ref"="grey24", "Brown_trout1" = "#b91a20",  "Chromis1" = "#ee4710", "Clownfish1" = "#00a100" , "Cunner1" = "#d9a100","Guppy1"="#00a9d9","Rainbow_trout1"="#0000d9","Zebrafish1"="#9900d9")) +
  theme(plot.title = element_text(hjust = 0,vjust = -1.7,size = 12))

AGR_plot <- ggplot(Data_fit[Data_fit$AGR>0,], aes(x = MassGrowth, y = 10^AGR_fit_TempAdj, group = Group)) +
  geom_abline(aes(slope=0.75, intercept=1, colour="ref"),linetype="dotted",show.legend = FALSE) +
  geom_abline(aes(slope=summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[2,1], intercept=summary(mcmc_AGR_RR_Evol_TempAdj)$solutions[1,1]*(-0.45),col="Mean"),show.legend = FALSE) +
  #geom_point(aes(x = MassGrowth, y = AGR_TempAdj*1000), shape = 16, colour = "grey90", alpha = 1, size = 0.6) +
  #geom_line(aes(x = MassGrowth, y = AGR_TempAdj*1000), colour = "grey90", size = 0.2, alpha = 1) +
  geom_line(data=Data_fit, aes(x = MassGrowth, y = AGR_fit_TempAdj*1000, group = Group, colour=Species), stat = "smooth", alpha = 0.3, size = 0.3, method = "lm", se=F) +
  geom_line(data=Data_fit, aes(x = MassGrowth, y = AGR_staticMT*1000, group = Species, colour=Species1), stat = "smooth", alpha = 1, size = 0.8, method = "lm", se=F, linetype="dotted",show.legend = FALSE) +
   theme_bw()+
  annotation_custom(nE)+annotation_custom(nAGR)+annotation_custom(nAGRu)+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  #annotate("text", x = -Inf, y = Inf, label = paste("a)"), vjust = 6, hjust = 6)+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1.0,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=14, margin = margin(t = 8, r = 0, b = 0, l = 0))) +
  #scale_colour_manual(values = c("red3", "dodgerblue2")) +
  scale_x_log10(labels = label_number(c(0.01,0.1,1,10,100)),breaks=c(0.01,0.1,1,10,100), limits=c(0.0009,68.2))+
  scale_y_log10(labels = label_number(c(0.01,0.1,1,10,100)),breaks=c(0.01,0.1,1,10,100))+
  labs(x=expression(Body ~ mass ~ (g)))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF",
                                              "ref"="grey24", "Brown_trout1" = "#b91a20",  "Chromis1" = "#ee4710", "Clownfish1" = "#00a100" , "Cunner1" = "#d9a100","Guppy1"="#00a9d9","Rainbow_trout1"="#0000d9","Zebrafish1"="#9900d9")) +
  theme(plot.title = element_text(hjust = 0,vjust = -1.7,size = 12))

GE_plot <- ggplot(Data_fit[Data_fit$GE>0,], aes(x = MassGrowth, y = 10^GE_fit_TempAdj, group = Group)) +
  geom_abline(aes(slope=0.75, intercept=0.6, colour="ref"),linetype="dotted",show.legend = FALSE) +
  geom_abline(aes(slope=summary(mcmc_GE_RR_Evol_TempAdj)$solutions[2,1], intercept=summary(mcmc_GE_RR_Evol_TempAdj)$solutions[1,1],col="Mean"),show.legend = FALSE) +
    #geom_point(aes(x = MassGrowth, y = GE_TempAdj), shape = 16, colour = "grey90", alpha = 1, size = 0.6) +
  #geom_line(aes(x = MassGrowth, y = GE_TempAdj), colour = "grey90", size = 0.2, alpha = 1) +
  geom_line(data=Data_fit, aes(x = MassGrowth, y = GE_fit_TempAdj, group = Group, colour=Species), stat = "smooth", alpha = 0.3, size = 0.3, method = "lm", se=F) +
  geom_line(data=Data_fit, aes(x = MassGrowth, y = GE_staticMT, group = Species, colour=Species1), stat = "smooth", alpha = 1, size = 0.8, method = "lm", se=F, linetype="dotted",show.legend = FALSE) +
  geom_abline(aes(slope=0.75, intercept=0.6, colour="ref"),linetype="dotted",show.legend = FALSE) +
  geom_abline(aes(slope=summary(mcmc_GE_RR_Evol_TempAdj)$solutions[2,1], intercept=summary(mcmc_GE_RR_Evol_TempAdj)$solutions[1,1],col="Mean"),show.legend = FALSE) +
  theme_bw()+
  annotation_custom(nF)+annotation_custom(nGE)+annotation_custom(nGEu)+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  #annotate("text", x = -Inf, y = Inf, label = paste("a)"), vjust = 6, hjust = 6)+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  #scale_colour_manual(values = c("red3", "dodgerblue2")) +
  scale_x_log10(labels = label_number(c(0.01,0.1,1,10,100)),breaks=c(0.01,0.1,1,10,100))+
  scale_y_log10(labels = label_number(c(0.2,0.5,2,10)),breaks=c(0.2,0.5,2,10),limits=c(0.2,20))+
  labs( x=expression(Body ~ mass ~ (g)))+
    scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF",
                                              "ref"="grey24", "Brown_trout1" = "#b91a20",  "Chromis1" = "#ee4710", "Clownfish1" = "#00a100" , "Cunner1" = "#d9a100","Guppy1"="#00a9d9","Rainbow_trout1"="#0000d9","Zebrafish1"="#9900d9")) +
  theme(plot.title = element_text(hjust = 0,vjust = -1.7,size = 12))


FAS_plot <- ggplot(Data_fit[Data_fit$FAS>0,], aes(x = Mass, y = 10^FAS_fit_TempAdj, group = Group)) +
  geom_abline(aes(slope=0.75, intercept=0.8, colour="ref"),linetype="dotted",show.legend = FALSE) +
  geom_abline(aes(slope=summary(mcmc_FAS_RR_Evol_TempAdj)$solutions[2,1], intercept=summary(mcmc_FAS_RR_Evol_TempAdj)$solutions[1,1],col="Mean"),show.legend = FALSE) +
  #geom_point(aes(x = Mass, y = FAS_TempAdj), shape = 16, colour = "grey90", alpha = 1, size = 0.6) +
  #geom_line(aes(x = Mass, y = FAS_TempAdj), colour = "grey90", size = 0.2, alpha = 1) +
  geom_line(data=Data_fit, aes(x = Mass, y = FAS_fit_TempAdj, group = Group, colour=Species), stat = "smooth", alpha = 0.3, size = 0.3, method = "lm", se=F) +
  geom_line(data=Data_fit, aes(x = Mass, y = FAS_staticM_TempAdj, group = Species, colour=Species1), stat = "smooth", alpha = 1, size = 0.8, method = "lm", se=F, linetype="dotted",show.legend = FALSE) +
  geom_abline(aes(slope=0.75, intercept=0.8, colour="ref"),linetype="dotted",show.legend = FALSE) +
  geom_abline(aes(slope=summary(mcmc_FAS_RR_Evol_TempAdj)$solutions[2,1], intercept=summary(mcmc_FAS_RR_Evol_TempAdj)$solutions[1,1],col="Mean"),show.legend = FALSE) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  #annotate("text", x = -Inf, y = Inf, label = paste("a)"), vjust = 6, hjust = 6)+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  #scale_colour_manual(values = c("red3", "dodgerblue2")) +
  scale_x_log10(labels = label_number(c(0.01,0.1,1,10,100)),breaks=c(0.01,0.1,1,10,100), limits=c(0.0009,68.2))+
  scale_y_log10(labels = label_number(c(2,4,7,10)),breaks=c(2,4,7,10),limits=c(2,13.5))+
  labs(x=expression(Body ~ mass ~ (g)))+
  annotate("rect",xmin =1 ,xmax=10,ymin=12,ymax = 13.5,fill="white")+
  annotation_custom(nD)+annotation_custom(nFAS)+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF",
                                              "ref"="grey24", "Brown_trout1" = "#b91a20",  "Chromis1" = "#ee4710", "Clownfish1" = "#00a100" , "Cunner1" = "#d9a100","Guppy1"="#00a9d9","Rainbow_trout1"="#0000d9","Zebrafish1"="#9900d9")) +
  theme(plot.title = element_text(hjust = 0,vjust = -1.2,size = 12))

AS_plot <- ggplot(Data_fit[Data_fit$AS>0,], aes(x = Mass, y = 10^AS_fit_TempAdj, group = Group)) +
  geom_abline(aes(slope=0.75, intercept=-0.2, colour="ref"),linetype="dotted",show.legend = FALSE) +
  geom_abline(aes(slope=summary(mcmc_AS_RR_Evol_TempAdj)$solutions[2,1], intercept=summary(mcmc_AS_RR_Evol_TempAdj)$solutions[1,1],col="Mean"),show.legend = FALSE) +
    #geom_point(aes(x = Mass, y = FAS_TempAdj), shape = 16, colour = "grey90", alpha = 1, size = 0.6) +
  #geom_line(aes(x = Mass, y = FAS_TempAdj), colour = "grey90", size = 0.2, alpha = 1) +
  geom_line(data=Data_fit, aes(x = Mass, y = AS_fit_TempAdj, group = Group, colour=Species), stat = "smooth", alpha = 0.3, size = 0.3, method = "lm", se=F) +
  geom_line(data=Data_fit, aes(x = Mass, y = AS_staticM_TempAdj, group = Species, colour=Species1), stat = "smooth", alpha = 1, size = 0.8, method = "lm", se=F, linetype="dotted",show.legend = FALSE) +
  theme_bw()+
  annotation_custom(nC)+annotation_custom(nAS)+annotation_custom(nASu)+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  #annotate("text", x = -Inf, y = Inf, label = paste("a)"), vjust = 6, hjust = 6)+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  #scale_colour_manual(values = c("red3", "dodgerblue2")) +
  scale_x_log10(labels = label_number(c(0.01,0.1,1,10,100)),breaks=c(0.01,0.1,1,10,100), limits=c(0.0009,68.2))+
  scale_y_log10(labels = label_number(c(0.01,0.1,1,10)),breaks=c(0.01,0.1,1,10),limits=c(0.01,50))+
  labs( x=expression(Body ~ mass ~ (g)))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF",
                                              "ref"="grey24", "Brown_trout1" = "#b91a20",  "Chromis1" = "#ee4710", "Clownfish1" = "#00a100" , "Cunner1" = "#d9a100","Guppy1"="#00a9d9","Rainbow_trout1"="#0000d9","Zebrafish1"="#9900d9")) +
  theme(plot.title = element_text(hjust = 0,vjust = -1.2,size = 12))
fig1 <- ggarrange(SMR_plot,MMR_plot,AS_plot,FAS_plot,AGR_plot,ncol = 3,nrow = 2, hjust = -2)

#Figure 2
#Data_raw <- read.csv("C:/Users/alero/My Drive/phd/Artikels/Big article/Data_analysis/Results_fit.csv", sep=";", encoding ="latin1")

nA <- grobTree(textGrob("A", x=0.03,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nB <- grobTree(textGrob("B", x=0.03,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nC <- grobTree(textGrob("C", x=0.03,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nD <- grobTree(textGrob("D", x=0.03,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nE <- grobTree(textGrob("E", x=0.03,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nF <- grobTree(textGrob("F", x=0.03,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))

nSMR <- grobTree(textGrob("Standart metabolic rate", x=0.13,  y=0.89, hjust=0,gp=gpar(col="black", fontsize=8)))
nMMR <- grobTree(textGrob("Maximun metabolic rate", x=0.13,  y=0.89, hjust=0,gp=gpar(col="black", fontsize=8)))
nAS <- grobTree(textGrob("Aerobic scope", x=0.13,  y=0.89, hjust=0,gp=gpar(col="black", fontsize=8)))
nFAS <- grobTree(textGrob("Factorial aerobic scope", x=0.13,  y=0.89, hjust=0,gp=gpar(col="black", fontsize=8)))
nGE <- grobTree(textGrob("Growth efficiency", x=0.13,  y=0.89, hjust=0,gp=gpar(col="black", fontsize=8)))
nAGR <- grobTree(textGrob("Abselute growth rate", x=0.13,  y=0.89, hjust=0,gp=gpar(col="black", fontsize=8)))

Scaling_exponents_all$Species <- factor(Scaling_exponents_all$Species, levels=c( "Brown_trout",  "Chromis","Clownfish","Mean","Cunner","Guppy","Rainbow_trout","Zebrafish"))
bSMR_plot <- ggplot(data = Scaling_exponents_all[Scaling_exponents_all$Parameter=="SMR",],aes(x=Level,y=bValue, col=Species))+ 
  geom_point(position = position_dodge(width = 0.7),size=2,aes(shape=Species))+
  geom_errorbar(aes(ymin=bLow, ymax=bHigh,),position = position_dodge(width = 0.7),width=0, size=0.7)+
  theme_bw()+theme(panel.spacing = unit(0.2, "lines"), strip.background = element_blank(),strip.placement = "outside", strip.text.x = element_text(size=16), strip.background.x = element_blank())+
  labs(title = expression(italic(b)[SMR]))+
  geom_hline(aes(yintercept=0.75,col="ref"), linetype="dotted")+
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(), axis.line.y = element_blank())+
  theme(legend.position="none")+
  annotation_custom(nA)+annotation_custom(nSMR)+
  force_panelsizes(rows = unit(1.0, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  #theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  scale_x_discrete(labels=c("Onto","Static"),limits = c("Onto","Static"))+
  scale_y_continuous(labels = label_number(c(0.2,0.6,1,1.4)),breaks =c(0.2,0.6,1,1.4), limits = c(0.2,1.62))+
  scale_colour_manual(name = NULL, values = c("Mean" = "grey", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="grey24"))+
  theme(plot.title = element_blank())+
  scale_shape_manual(name = NULL, values = c(16,  16,16, 15 , 16,16,16,16))

bMMR_plot <- ggplot(data = Scaling_exponents_all[Scaling_exponents_all$Parameter=="MMR",],aes(x=Level,y=bValue, col=Species))+ 
  geom_point(position = position_dodge(width = 0.7),size=2,aes(shape=Species))+
  geom_errorbar(aes(ymin=bLow, ymax=bHigh,),position = position_dodge(width = 0.7),width=0, size=0.7)+
  theme_bw()+theme(panel.spacing = unit(0.2, "lines"), strip.background = element_blank(),strip.placement = "outside", strip.text.x = element_text(size=16), strip.background.x = element_blank())+
  labs(title = expression(italic(b)[MMR]))+
  geom_hline(aes(yintercept=0.75,col="ref"), linetype="dotted")+
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(), axis.line.y = element_blank())+
  theme(legend.position="none")+
  annotation_custom(nB)+annotation_custom(nMMR)+
  force_panelsizes(rows = unit(1.0, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  #theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  scale_x_discrete(labels=c("Onto","Static"),limits = c("Onto","Static"))+
  scale_y_continuous(labels = label_number(c(0.6,0.8,1,1.2)),breaks =c(0.6,0.8,1,1.2), limits = c(0.45,1.4))+
  scale_colour_manual(name = NULL, values = c("Mean" = "grey", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="grey24"))+
  theme(plot.title = element_blank())+
  scale_shape_manual(name = NULL, values = c(16,  16,16, 15 , 16,16,16,16))

bAGR_plot <- ggplot(data = Scaling_exponents_all[Scaling_exponents_all$Parameter=="AGR",],aes(x=Level,y=bValue, col=Species))+ 
  geom_point(position = position_dodge(width = 0.7),size=2,aes(shape=Species))+
  geom_errorbar(aes(ymin=bLow, ymax=bHigh,),position = position_dodge(width = 0.7),width=0, size=0.7)+
  theme_bw()+theme(panel.spacing = unit(0.2, "lines"), strip.background = element_blank(),strip.placement = "outside", strip.text.x = element_text(size=16), strip.background.x = element_blank())+
  labs(title = expression(italic(b)[AGR]))+
  geom_hline(aes(yintercept=0.75,col="ref"), linetype="dotted")+
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(), axis.line.y = element_blank())+
  theme(legend.position="none")+
  annotation_custom(nE)+annotation_custom(nAGR)+
  force_panelsizes(rows = unit(1.0, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  #theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=10,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  scale_x_discrete(labels=c("Ontogenetic","Static"),limits = c("Onto","Static"))+
  scale_y_continuous(labels =(c(-3,-1.5,0.5,2)),breaks=c(-3,-1.5,0.5,2),limits = c(-3,2.9))+
  scale_colour_manual(name = NULL, values = c("Mean" = "grey", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="grey24"))+
  theme(plot.title = element_blank())+
  scale_shape_manual(name = NULL, values = c(16,  16,16, 15 , 16,16,16,16))

bFAS_plot <- ggplot(data = Scaling_exponents_all[Scaling_exponents_all$Parameter=="FAS",],aes(x=Level,y=bValue, col=Species))+ 
  geom_point(position = position_dodge(width = 0.7),size=2,aes(shape=Species))+
  geom_errorbar(aes(ymin=bLow, ymax=bHigh,),position = position_dodge(width = 0.7),width=0, size=0.7)+
  theme_bw()+theme(panel.spacing = unit(0.2, "lines"), strip.background = element_blank(),strip.placement = "outside", strip.text.x = element_text(size=16), strip.background.x = element_blank())+
  labs(title = expression(italic(b)[FAS]))+
  #geom_hline(aes(yintercept=0.75,col="ref"), linetype="dotted")+
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(), axis.line.y = element_blank())+
  theme(legend.position="none")+
  annotation_custom(nD)+annotation_custom(nFAS)+
  force_panelsizes(rows = unit(1.0, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  #theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  scale_x_discrete(labels=c("Ontogenetic","Static"),limits = c("Onto","Static"))+
  scale_y_continuous(labels =(c(-0.6,-0.3,0,0.4)),breaks=c(-0.6,-0.3,0,0.4),limits = c(-0.67,0.9))+
  scale_colour_manual(name = NULL, values = c("Mean" = "grey", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="grey24"))+
  theme(plot.title = element_blank())+
  scale_shape_manual(name = NULL, values = c(16,  16,16, 15 , 16,16,16,16))

bGE_plot <- ggplot(data = Scaling_exponents_all[Scaling_exponents_all$Parameter=="GE",],aes(x=Level,y=bValue, col=Species))+ 
  geom_point(position = position_dodge(width = 0.7),size=2,aes(shape=Species))+
  geom_errorbar(aes(ymin=bLow, ymax=bHigh,),position = position_dodge(width = 0.7),width=0, size=0.7)+
  theme_bw()+theme(panel.spacing = unit(0.2, "lines"), strip.background = element_blank(),strip.placement = "outside", strip.text.x = element_text(size=16), strip.background.x = element_blank())+
  labs(title = expression(italic(b)[GE]))+
  #geom_hline(aes(yintercept=0.75,col="ref"), linetype="dotted")+
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(), axis.line.y = element_blank())+
  theme(legend.position="none")+
  annotation_custom(nF)+annotation_custom(nGE)+
  force_panelsizes(rows = unit(1.0, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  #theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=10,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  scale_x_discrete(labels=c("Ontogenetic","Static"),limits = c("Onto","Static"))+
  scale_y_continuous(labels = (c(-3.0,-1.7,-0.5,0.7)),breaks=c(-3.0,-1.7,-0.5,0.7), limits = c(-3.28,1.6))+
  scale_colour_manual(name = NULL, values = c("Mean" = "grey", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="grey24"))+
  theme(plot.title = element_blank())+
  scale_shape_manual(name = NULL, values = c(16,  16,16, 15 , 16,16,16,16))

bAS_plot <- ggplot(data = Scaling_exponents_all[Scaling_exponents_all$Parameter=="AS",],aes(x=Level,y=bValue, col=Species))+ 
  geom_point(position = position_dodge(width = 0.7),size=2,aes(shape=Species))+
  geom_errorbar(aes(ymin=bLow, ymax=bHigh,),position = position_dodge(width = 0.7),width=0, size=0.7)+
  theme_bw()+theme(panel.spacing = unit(0.2, "lines"), strip.background = element_blank(),strip.placement = "outside", strip.text.x = element_text(size=16), strip.background.x = element_blank())+
  labs(title = expression(italic(b)[AS]))+
  geom_hline(aes(yintercept=0.75,col="ref"), linetype="dotted")+
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(), axis.line.y = element_blank())+
  theme(legend.position="none")+
  annotation_custom(nC)+annotation_custom(nAS)+
  force_panelsizes(rows = unit(1.0, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  #theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  scale_x_discrete(labels=c("Onto","Static"),limits = c("Onto","Static"))+
  scale_y_continuous(labels = (c(0.4,1.0,1.3,0.7)),breaks=c(0.4,1.0,1.3,0.7),limits = c(0.4,1.6))+
  scale_colour_manual(name = NULL, values = c("Mean" = "grey", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="grey24"))+
  theme(plot.title = element_blank())+
  scale_shape_manual(name = NULL, values = c(16,  16,16, 15 , 16,16,16,16))
fig2 <- ggarrange(bSMR_plot,bMMR_plot,bAS_plot,bFAS_plot,bAGR_plot,bGE_plot,ncol = 2,nrow = 3, hjust = -2)


#Figure 3
nA <- grobTree(textGrob("A", x=0.02,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nB <- grobTree(textGrob("B", x=0.02,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nC <- grobTree(textGrob("C", x=0.02,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nD <- grobTree(textGrob("D", x=0.02,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))

nSMR <- grobTree(textGrob("Standart metabolic rate", x=0.09,  y=0.94, hjust=0,gp=gpar(col="black", fontsize=8)))
nMMR <- grobTree(textGrob("Maximun meatabolic rate", x=0.09,  y=0.94, hjust=0,gp=gpar(col="black", fontsize=8)))
nAS <- grobTree(textGrob("Aerobic scope", x=0.09,  y=0.94, hjust=0,gp=gpar(col="black", fontsize=8)))
nFAS <- grobTree(textGrob("Factorial aerobic scope", x=0.09,  y=0.94, hjust=0,gp=gpar(col="black", fontsize=8)))


cv_bSMR <- sd(Results_fit$Slope_SMR_fit)/mean(Results_fit$Slope_SMR_fit)*100
cv_bMMR <- sd(na.omit(Results_fit$Slope_MMR_fit))/mean(na.omit(Results_fit$Slope_MMR_fit))*100
cv_bAS <- sd(na.omit(Results_fit$Slope_AS_fit))/mean(na.omit(Results_fit$Slope_AS_fit))*100
cv_bFAS <- sd(na.omit(Results_fit$Slope_FAS_fit))/mean(na.omit(Results_fit$Slope_FAS_fit))*100
cv_bAGR <- sd(na.omit(Results_fit$Slope_AGR_fit))/mean(na.omit(Results_fit$Slope_AGR_fit))*100

Results_fit$Plot_treatment2 <- "1"
Results_fit[Results_fit$Treatment=="f",]$Plot_treatment2 <- "2"

bSMRvbAGR_plot <- ggplot(data = Results_fit, aes(x=Slope_AGR_fit,y=Slope_SMR_fit))+
  geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  geom_point(aes(col=Species,shape=Plot_treatment2),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.8, method = "lm", se=F )+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  annotation_custom(nA)+  annotation_custom(nSMR)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  labs(y=expression(SMR~scaling ~ exponent ~ (italic(b)) ), x=expression(AGR~scaling ~ exponent ~ (italic(b)) ))+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,0,1,0, "mm")))+ 
  scale_shape_manual(name=NULL, values=c("1"=16,"2"=17))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="grey24")) +
  annotate("rect",xmin = c(-1),xmax=c(1),ymin=c(-0.3),ymax = c(-0.22),fill=c("white"))+
  scale_y_continuous(limits = c(-0.3,1.4),breaks = c(0,0.4,0.8,1.2))+
  annotate(geom = "text",y=-0.29,x=0.2, label=paste("r","==","\"0.82, p < 0.001\"", sep=""), parse=T)+theme(plot.margin = unit(c(0,0,0,0), 'lines'))


bAGR_dens <- ggplot(Results_fit, aes(x=Slope_AGR_fit,fill=Species))+
  geom_density(alpha = 0.8,linewidth = 0.2)+
  theme_void() + 
  #xlim(-2.3,1.4) +
  theme(legend.position = "none") +theme(plot.margin = unit(c(0,-1.8,-20.4,0), 'lines'))+
  force_panelsizes(rows = unit(0.5, "in"), cols = unit(2, "in"))+
  scale_fill_manual(name = NULL, values = c( "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 

bSMR_dens <- ggplot(Results_fit, aes(x=Slope_SMR_fit,fill=Species))+
  geom_density(alpha = 0.8,linewidth = 0.2) + 
  theme_void() + 
  theme(legend.position = "none") +
  coord_flip()+theme(plot.margin = unit(c(0,0,0,-37), 'lines'))+
  force_panelsizes(rows = unit(2, "in"), cols = unit(0.5, "in"))+
  scale_fill_manual(name = NULL, values = c( "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 

ggarrange(bAGR_dens,NULL,bSMRvbAGR_plot,bSMR_dens,ncol = 2,nrow = 2)







summary(lmer(data = Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_SMR_fit~Slope_AGR_fit+(1|Species)))
r.squaredGLMM(lmer(data = Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_SMR_fit~Slope_AGR_fit+(1|Species)))
cor.test(Results_fit$Slope_SMR_fit,Results_fit$Slope_AGR_fit)
summary(lm(data = Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_SMR_fit~Slope_AGR_fit))
confint(lm(data = Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_SMR_fit~Slope_AGR_fit))
bMMRvbAGR_plot <- ggplot(data = Results_fit, aes(x=Slope_AGR_fit,y=Slope_MMR_fit))+
  geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
    geom_point(aes(col=Species, shape=Plot_treatment2),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.8, method = "lm", se=F )+
   theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  annotation_custom(nB)+annotation_custom(nMMR)+
  scale_shape_manual(name=NULL, values=c("1"=16,"2"=17))+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  #scale_colour_manual(values = c("red3", "dodgerblue2")) +
  labs(y=expression(MMR~scaling ~ exponent ~ (italic(b)) ), x=expression(AGR~scaling ~ exponent ~ (italic(b)) ))+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,0,1,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="grey24")) +
  #annotate("rect",xmin = c(-1),xmax=c(1),ymin=c(-0.3),ymax = c(-0.22),fill=c("white"))+
  scale_y_continuous(limits = c(0.25,1.2),breaks = c(0.4,0.6,0.8,1))+
  annotate(geom = "text",y=0.26,x=-1.2, label=paste("r","==","\"0.21, p < 0.001\"", sep=""), parse=T)

bMMR_dens <- ggplot(Results_fit, aes(x=Slope_MMR_fit,fill=Species))+
  geom_density(alpha = 0.8,linewidth = 0.2) + 
  theme_void() + 
  theme(legend.position = "none") +
  coord_flip()+
  force_panelsizes(rows = unit(2, "in"), cols = unit(0.5, "in"))+
  scale_fill_manual(name = NULL, values = c( "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 



summary(lmer(data = Results_fit[!is.na(Results_fit$Slope_MMR_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_MMR_fit~Slope_AGR_fit+(1|Species)))
r.squaredGLMM(lmer(data = Results_fit[!is.na(Results_fit$Slope_MMR_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_MMR_fit~Slope_AGR_fit+(1|Species)))
summary(lm(data = Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_MMR_fit~Slope_AGR_fit))
cor.test(Results_fit$Slope_MMR_fit,Results_fit$Slope_AGR_fit)


bFASvbAGR_plot <- ggplot(data = Results_fit, aes(x=Slope_AGR_fit,y=Slope_FAS_fit))+
  geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  geom_point(aes(col=Species, shape=Plot_treatment2),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.8, method = "lm", se=F )+
   theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  annotation_custom(nD)+annotation_custom(nFAS)+
  scale_shape_manual(name=NULL, values=c("1"=16,"2"=17))+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_text(size=12, margin = margin(t = 8, r = 0, b = 0, l = 0)), axis.title.x = element_text(size=12, margin = margin(t = 8, r = 0, b = 0, l = 0))) +
  xlim(-2.3,1.3)+
  labs(y=expression(Scaling ~ exponent ~ (italic(b))~"for"~metabolic~trait ), x=expression(Scaling ~ exponent ~ (italic(b))~"for"~growth~rate) )+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,0,1,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="grey24"))+
  annotate("rect",xmin = c(-0.9),xmax=c(-0),ymin=c(-0.53),ymax = c(-0.46),fill=c("white"))+
  scale_y_continuous(breaks = c(-0.5,-0.2,0.2,0.5),limits = c(-0.53,0.55))+
  annotate(geom = "text",y=-0.52,x=-1.2, label=paste("r","==","\"-0.72, p < 0.001\"", sep=""), parse=T)

bFAS_dens <- ggplot(Results_fit, aes(x=Slope_FAS_fit,fill=Species))+
  geom_density(alpha = 0.8,linewidth = 0.2) + 
  theme_void() + 
  theme(legend.position = "none") +
  xlim(-0.73,0.55) +
  coord_flip()+
  force_panelsizes(rows = unit(2, "in"), cols = unit(0.5, "in"))+
  scale_fill_manual(name = NULL, values = c( "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 

bAGR_dens <- ggplot(Results_fit, aes(x=Slope_AGR_fit,fill=Species))+
  geom_density(alpha = 0.8,linewidth = 0.2)+
  theme_void() + 
  xlim(-2.3,1.4) +
  theme(legend.position = "none") +
  force_panelsizes(rows = unit(0.5, "in"), cols = unit(2, "in"))+
  scale_fill_manual(name = NULL, values = c( "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 




summary(lmer(data = Results_fit[!is.na(Results_fit$Slope_FAS_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_FAS_fit~Slope_AGR_fit+(1|Species)))
r.squaredGLMM(lmer(data = Results_fit[!is.na(Results_fit$Slope_FAS_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_FAS_fit~Slope_AGR_fit+(1|Species)))
summary(lm(data = Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_FAS_fit~Slope_AGR_fit))
cor.test(Results_fit$Slope_AS_fit,Results_fit$Slope_AGR_fit)



bASvbAGR_plot <- ggplot(data = Results_fit, aes(x=Slope_AGR_fit,y=Slope_AS_fit))+
  geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  geom_point(aes(col=Species, shape=Plot_treatment2),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.8, method = "lm", se=F )+
   theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  annotation_custom(nC)+annotation_custom(nAS)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  #ylim(-0.73,0.55) +
  labs(y=expression(Metabolic ~scaling ~ exponent ~ (italic(b)) ), x=expression(AGR~scaling ~ exponent ~ (italic(b)) ))+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,0,1,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="grey24"))+
  scale_shape_manual(name=NULL, values=c("1"=16,"2"=17))+
  #annotate("rect",xmin = c(-0.9),xmax=c(0.4),ymin=c(0.05),ymax = c(0.15),fill=c("white"))+
  scale_y_continuous(breaks = c(0,0.4,0.8,1.2),limits = c(0,1.4))+
  annotate(geom = "text",y=0.01,x=-1.2, label=paste("r","==","\"0.06, p = 0.289\"", sep=""), parse=T)

bAS_dens <- ggplot(Results_fit, aes(x=Slope_AS_fit,fill=Species))+
  geom_density(alpha = 0.8,linewidth = 0.2) + 
  theme_void() + 
  theme(legend.position = "none") +
  coord_flip()+
  force_panelsizes(rows = unit(2, "in"), cols = unit(0.5, "in"))+
  scale_fill_manual(name = NULL, values = c( "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 
ggarrange(bASvbAGR_plot,bAS_dens, ncol=1,nrow=1)

summary(lmer(data = Results_fit[!is.na(Results_fit$Slope_FAS_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_FAS_fit~Slope_AGR_fit+(1|Species)))
r.squaredGLMM(lmer(data = Results_fit[!is.na(Results_fit$Slope_FAS_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_FAS_fit~Slope_AGR_fit+(1|Species)))
summary(lm(data = Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_AS_fit~Slope_AGR_fit))
cor.test(Results_fit$Slope_FAS_fit,Results_fit$Slope_AGR_fit)

check_model(lmer(data = Results_fit[!is.na(Results_fit$Slope_SMR_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_SMR_fit~Slope_AGR_fit+(1|Species)))
Trout_bSMR_bAGR <- lm(data = Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_AS_fit~Slope_AGR_fit)
check_model(Trout_bSMR_bAGR)
check_outliers(Trout_bSMR_bAGR)
cor.test(Results_fit[!Results_fit$Species=="Brown_trout",]$Slope_SMR_fit,Results_fit[!Results_fit$Species=="Brown_trout",]$Slope_AGR_fit)
cor.test(Results_fit[!Results_fit$Species=="Brown_trout",]$Slope_MMR_fit,Results_fit[!Results_fit$Species=="Brown_trout",]$Slope_AGR_fit)
cor.test(Results_fit[!Results_fit$Species=="Brown_trout",]$Slope_AS_fit,Results_fit[!Results_fit$Species=="Brown_trout",]$Slope_AGR_fit)
cor.test(Results_fit[!Results_fit$Species=="Brown_trout",]$Slope_FAS_fit,Results_fit[!Results_fit$Species=="Brown_trout",]$Slope_AGR_fit)

Results_fit[!Results_fit$Species=="Brown_trout",]


ggarrange(bSMRvbAGR_plot,bMMRvbAGR_plot,bASvbAGR_plot,bFASvbAGR_plot,ncol = 2,nrow = 2)


summary(lmer(data = Results_fit[!is.na(Results_fit$Slope_AS_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_AS_fit~Slope_AGR_fit+(1|Species)))
r.squaredGLMM(lmer(data = Results_fit[!is.na(Results_fit$Slope_AS_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_AS_fit~Slope_AGR_fit+(1|Species)))
summary(lm(data = Results_fit[!is.na(Results_fit$Slope_AGR_fit),],Slope_AS_fit~Slope_AGR_fit))



Legend <- ggplot(data = Results_fit, aes(x=Slope_AGR_fit,y=Slope_SMR_fit))+
  geom_point(aes(col=Species))+
  geom_line(aes(col="Mean"),stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F )+
  geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  theme_bw()+
   theme(legend.position="bottom",legend.key.size = unit(7,"mm"),legend.text = element_text(size=10,margin = margin(r = 3,l=0)),
        legend.box.background = element_rect(colour = "black"),legend.spacing.x = unit(0.0, 'mm'))+
   labs(y=expression(Scaling ~ exponent ~ (italic(b)) ~ "for" ~ SMR), x=expression(Scaling ~ exponent ~ (italic(b)) ~ "for" ~ AGR))+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,1,1,1, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "#00B6EB", "Brown_trout" = "black",  "Chromis" = "#DA1F26", "Clownfish" = "#F26739" , "Cunner" = "#00BE00","Guppy"="#FFBE00","Rainbow_trout"="#0000FF","Zebrafish"="black","ref"="#B400FF"),
                      guide = guide_legend(title = "", nrow = 1,override.aes = list(linetype = c("solid",rep("blank", 7),"dotted"),size=2,shape = c(NA,rep(16, 7),NA))), 
                      labels =c("Regression", "Brown trout","Damselfish","Clownfish","Cunner","Guppy","Rainbow trout","Zebrafish", "1:1 Reference line")) 

#Figure 4
nA <- grobTree(textGrob("A", x=0.05,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nB <- grobTree(textGrob("B", x=0.05,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nC <- grobTree(textGrob("C", x=0.05,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nD <- grobTree(textGrob("D", x=0.05,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))
nE <- grobTree(textGrob("E", x=0.05,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12,fontface="bold")))

SMR_cor_plot <- ggplot(data = Results_fit, aes(x=Intersept_SMR_fit,y=Slope_SMR_fit))+
  geom_point(aes(col=Species),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.8, method = "lm", se=F )+
  #geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  #geom_line(aes(col=Species), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  #geom_line(aes(col=Species, group=Plot_group), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(1, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  annotation_custom(nA)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=12, margin = margin(t = 8, r = 0, b = 0, l = 0))) +
  ylim(-0.3,1.5) +
  labs(title =expression(italic(b)[SMR]) , x=expression(italic(a)[SMR]) )+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,0,1,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="black"))+
  theme(plot.title = element_text(hjust = 0,vjust = -1.7,size = 12))+
  #annotate(geom = "text",y=0,x=-1.6, label=paste("italic(b)[SMR] == 0.75-0.04 *italic(a)[SMR]"), parse=T,size=3.5)+
  annotate(geom = "text",y=-0.2,x=-1.6, label=paste("Cor","==","\"-0.4, p < 0.001\"", sep=""), parse=T,size=3.5)

cor.test(Results_fit$Slope_SMR_fit,Results_fit$Intersept_SMR_fit)
summary(lmer(data = Results_fit[!is.na(Results_fit$Intersept_SMR_fit)&!is.na(Results_fit$Slope_SMR_fit),],Slope_SMR_fit~Intersept_SMR_fit+(1|Species)))
r.squaredGLMM(lmer(data = Results_fit[!is.na(Results_fit$Intersept_SMR_fit)&!is.na(Results_fit$Slope_SMR_fit),],Slope_SMR_fit~Intersept_SMR_fit+(1|Species)))

MMR_cor_plot <- ggplot(data = Results_fit, aes(x=Intersept_MMR_fit,y=Slope_MMR_fit))+
  geom_point(aes(col=Species),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.8, method = "lm", se=F )+
  #geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  #geom_line(aes(col=Species), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  #geom_line(aes(col=Species, group=Plot_group), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(1, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  annotation_custom(nB)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=12, margin = margin(t = 8, r = 0, b = 0, l = 0))) +
  #ylim(-0.7,0.55) +
  labs(title=expression(italic(b)[MMR] ), x=expression(italic(a)[MMR] ))+
  theme(plot.title = element_text(hjust = 0,vjust = -1.7,size = 12))+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,0,1,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="black"))+
  #annotate(geom = "text",y=0.4,x=0.4, label=paste("italic(b)[MMR] == 0.8+0.02 *italic(a)[MMR]"), parse=T)+
  annotate(geom = "text",y=0.3,x=0.4, label=paste("Cor","==","\"0.25, p < 0.001\"", sep=""), parse=T,size=3.5)

cor.test(Results_fit$Slope_MMR_fit,Results_fit$Intersept_MMR_fit)
summary(lmer(data = Results_fit[!is.na(Results_fit$Intersept_MMR_fit)&!is.na(Results_fit$Slope_MMR_fit),],Slope_MMR_fit~Intersept_MMR_fit+(1|Species)))
r.squaredGLMM(lmer(data = Results_fit[!is.na(Results_fit$Intersept_MMR_fit)&!is.na(Results_fit$Slope_MMR_fit),],Slope_MMR_fit~Intersept_MMR_fit+(1|Species)))

AGR_cor_plot <- ggplot(data = Results_fit, aes(x=Intersept_AGR_fit,y=Slope_AGR_fit))+
  geom_point(aes(col=Species),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.8, method = "lm", se=F )+
  #geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  #geom_line(aes(col=Species), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  #geom_line(aes(col=Species, group=Plot_group), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(1, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  annotation_custom(nC)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=12, margin = margin(t = 8, r = 0, b = 0, l = 0))) +
  ylim(-2.2,1.6) +
  labs(title=expression(italic(b)[AGR] ), x=expression(italic(a)[AGR] ))+
  theme(plot.title = element_text(hjust = 0,vjust = -1.7,size = 12))+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,0,1,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="black"))+
  #annotate(geom = "text",y=-1.5,x=-2.6, label=paste("italic(b)[AGR] == 0.4+0.2 *italic(a)[AGR]"), parse=T)+
  annotate(geom = "text",y=-1.5,x=-2.6, label=paste("Cor","==","\"-0.46, p < 0.001\"", sep=""), parse=T, size=3.5)

cor.test(Results_fit$Slope_AGR_fit,Results_fit$Intersept_AGR_fit)
summary(lmer(data = Results_fit[!is.na(Results_fit$Intersept_AGR_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_AGR_fit~Intersept_AGR_fit+(1|Species)))
r.squaredGLMM(lmer(data = Results_fit[!is.na(Results_fit$Intersept_AGR_fit)&!is.na(Results_fit$Slope_AGR_fit),],Slope_AGR_fit~Intersept_AGR_fit+(1|Species)))

FAS_cor_plot <- ggplot(data = Results_fit, aes(x=Intersept_FAS_fit,y=Slope_FAS_fit))+
  geom_point(aes(col=Species),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.8, method = "lm", se=F )+
  #geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  #geom_line(aes(col=Species), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  #geom_line(aes(col=Species, group=Plot_group), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(1.5, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  annotation_custom(nE)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=12, margin = margin(t = 8, r = 0, b = 0, l = 0))) +
  ylim(-0.7,0.4) +
  #xlim(0,0.9) +
  labs(title=expression(italic(b)[FAS]) , x=expression(italic(a)[FAS] ))+
  theme(plot.title = element_text(hjust = 0,vjust = -1.7,size = 12))+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,0,1,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="black"))+
  #annotate(geom = "text",y=-0.62,x=0.7, label=paste("italic(b)[FAS] == 0.09-0.04 *italic(a)[FAS]"), parse=T)+
  annotate(geom = "text",y=-0.65,x=0.7, label=paste("Cor","==","\"0.15, p = 0.011\"", sep=""), parse=T,size=3.5)

cor.test(Results_fit$Slope_FAS_fit,Results_fit$Intersept_FAS_fit)
summary(lmer(data = Results_fit[!is.na(Results_fit$Intersept_FAS_fit)&!is.na(Results_fit$Slope_FAS_fit),],Slope_FAS_fit~Intersept_FAS_fit+(1|Species)))
r.squaredGLMM(lmer(data = Results_fit[!is.na(Results_fit$Intersept_FAS_fit)&!is.na(Results_fit$Slope_FAS_fit),],Slope_FAS_fit~Intersept_FAS_fit+(1|Species)))

GE_cor_plot <- ggplot(data = Results_fit, aes(x=Intersept_GE_fit,y=Slope_GE_fit))+
  geom_point(aes(col=Species),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.8, method = "lm", se=F )+
  #geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  #geom_line(aes(col=Species), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  #geom_line(aes(col=Species, group=Plot_group), stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(1.5, "in"), cols = unit(1.5, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  annotation_custom(nD)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=12, margin = margin(t = 8, r = 0, b = 0, l = 0))) +
  ylim(-2.9,0.9) +
  #xlim(0,0.9) +
  labs(title=expression(italic(b)[GE] ), x=expression(italic(a)[GE] ))+
  theme(plot.title = element_text(hjust = 0,vjust = -1.7,size = 12))+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,0,1,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="black"))+
  #annotate(geom = "text",y=-2.3,x=0.21, label=paste("italic(b)[GE] == -0.9+1.3 *italic(a)[GE]"), parse=T)+
  annotate(geom = "text",y=-2.8,x=0.14, label=paste("Cor","==","\"0.51, p < 0.001\"", sep=""), parse=T,size=3.5)

cor.test(Results_fit$Slope_GE_fit,Results_fit$Intersept_GE_fit)
summary(lmer(data = Results_fit[!is.na(Results_fit$Intersept_GE_fit)&!is.na(Results_fit$Slope_GE_fit),],Slope_GE_fit~Intersept_GE_fit+(1|Species)))
r.squaredGLMM(lmer(data = Results_fit[!is.na(Results_fit$Intersept_GE_fit)&!is.na(Results_fit$Slope_GE_fit),],Slope_GE_fit~Intersept_GE_fit+(1|Species)))

Legend <- ggplot(data = Results_fit, aes(x=Slope_AGR_fit,y=Slope_SMR_fit))+
  geom_point(aes(col=Species))+
  geom_line(aes(col="Mean"),stat = "smooth", alpha = 0.5, linewidth = 0.8, method = "lm", se=F )+
  #geom_abline(aes(slope=1, intercept=0, colour="ref"),linetype="dotted",show.legend = FALSE) +
  theme_bw()+
  theme(legend.position="bottom",legend.key.size = unit(7,"mm"),legend.text = element_text(size=10,margin = margin(r = 3,l=0)),
        legend.box.background = element_rect(colour = "black"),legend.spacing.x = unit(0.0, 'mm'))+
  labs(y=expression(Scaling ~ exponent ~ (italic(b)) ~ "for" ~ SMR), x=expression(Scaling ~ exponent ~ (italic(b)) ~ "for" ~ AGR))+
  theme(strip.background = element_rect(linewidth=0.3), strip.text.x = element_text(size=10, colour="black", margin=margin(1,1,1,1, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "#00B6EB", "Brown_trout" = "black",  "Chromis" = "#DA1F26", "Clownfish" = "#F26739" , "Cunner" = "#00BE00","Guppy"="#FFBE00","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF","ref"="#B400FF"),
                      guide = guide_legend(title = "", nrow = 1,override.aes = list(linetype = c("solid",rep("blank", 7)),size=2,shape = c(NA,rep(16, 7)))), 
                      labels =c("Regression", "Brown trout","Damselfish","Clownfish","Cunner","Guppy","Rainbow trout","Zebrafish")) 

#Table 1

Scaling_exponents_m <-Scaling_exponents_all[Scaling_exponents_all$Species=="Mean",]
Scaling_exponents_adj_m
Table_Scaling <- data.frame("Parameter"=unique(Scaling_exponents_m$Parameter),
                            "Ontogenetic"=c(paste(round(Scaling_exponents_m[Scaling_exponents_m$Level=="Onto",]$bValue,2)," [",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Onto",]$bHigh,2),";",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Onto",]$bLow,2),"]",sep="")),
                            "Static"=c(paste(round(Scaling_exponents_m[Scaling_exponents_m$Level=="Static",]$bValue,2)," [",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Static",]$bHigh,2),";",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Static",]$bLow,2),"]",sep="")),
                            "Evolutionary"=c(paste(round(Scaling_exponents_m[Scaling_exponents_m$Level=="Evo",]$bValue,2)," [",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Evo",]$bHigh,2),";",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Evo",]$bLow,2),"]",sep="")),
                            "Ontogenetic"=c(paste(round(Scaling_exponents_m[Scaling_exponents_m$Level=="Onto",]$aValue,2)," [",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Onto",]$aHigh,2),";",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Onto",]$aLow,2),"]",sep="")),
                            "Static"=c(paste(round(Scaling_exponents_m[Scaling_exponents_m$Level=="Static",]$aValue,2)," [",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Static",]$aHigh,2),";",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Static",]$aLow,2),"]",sep="")),
                            "Evolutionary"=c(paste(round(Scaling_exponents_m[Scaling_exponents_m$Level=="Evo",]$aValue,2)," [",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Evo",]$aHigh,2),";",round(Scaling_exponents_m[Scaling_exponents_m$Level=="Evo",]$aLow,2),"]",sep="")))
Scaling_exponents_all[(Scaling_exponents_all$Parameter=="SMR"|Scaling_exponents_all$Level=="Evo"),]
# Sub table 2
Sub_table <- data.frame(matrix(nrow=0,ncol = 10))
names(Sub_table)<- c("Trait","Level", "Zebrafish", "Rainbow trout", "Guppy","Damselfish","Brown trout", "Cunner", "Clownfish",   "Overall")
Trait <- c("SMR","MMR","AS","FAS","AGR","GE")
Level <- c("Onto","Static")
for (x in 1:6) {
  for (y in 1:2) {
    ifelse(x==1&y==1,nr<-1,nr<-nr+1)
    Sub_table[nr,1] <- Trait[x]
    Sub_table[nr,2] <- Level[y]
    Data_ID <- Scaling_exponents_all[Scaling_exponents_all$Parameter==Trait[x]&Scaling_exponents_all$Level==Level[y],]
    for (q in 1:8) {
      Sub_table[nr,2+q] <- paste(round(Data_ID$bValue[q],2)," [",round(Data_ID$bLow[q],2),"; ",round(Data_ID$bHigh[q],2),"]",sep="")
    }
  }
}
write.csv(Sub_table,"C:/Users/alero/My Drive/phd/Artikels/Big article/Data_analysis/Sub_table.csv", sep="T")

### SUB figures

#Figure S3


cv_bSMR <- sd(na.omit(Results_fit$Slope_SMR_fit))/mean(na.omit(Results_fit$Slope_SMR_fit))*100
cv_bMMR <- sd(na.omit(Results_fit$Slope_MMR_fit))/mean(na.omit(Results_fit$Slope_MMR_fit))*100
cv_bAS <- sd(na.omit(Results_fit$Slope_AS_fit))/mean(na.omit(Results_fit$Slope_AS_fit))*100
cv_bFAS <- sd(na.omit(Results_fit$Slope_FAS_fit))/mean(na.omit(Results_fit$Slope_FAS_fit))*100
cv_bAGR <- sd(na.omit(Results_fit$Slope_AGR_fit))/mean(na.omit(Results_fit$Slope_AGR_fit))*100
cv_bGE <- sd(na.omit(Results_fit$Slope_GE_fit))/mean(na.omit(Results_fit$Slope_GE_fit))*100


Results_fit$Species2 <- Results_fit$Species
Results_fit[Results_fit$Species=="Brown_trout",]$Species2 <- "Brown trout"
Results_fit[Results_fit$Species=="Rainbow_trout",]$Species2 <- "Rainbow trout"
Results_fit[Results_fit$Species=="Chromis",]$Species2 <- "Damselfish"
Results_fit[Results_fit$Species=="Clownfish",]$Species2 <- "Anemonefish"
Results_fit2 <- Results_fit
Results_fit2$Species2 <- "All species"
Results_fit3 <- rbind(Results_fit,Results_fit2)

#Interspecies  corelations
bcor_table_species <- data.frame(matrix(nrow=0,ncol=6))
names(bcor_table_species)<-c("Species","Par_1","Par_2","Cor","p_val", "Compereson")
Parameters_1 <- c("SMR","MMR","AS","FAS")
Parameters_2 <- c("AGR")
Species_1 <- c("Mean","Brown_trout","Chromis","Clownfish","Cunner","Guppy","Rainbow_trout","Zebrafish")
  for (x in 1:length(Parameters_1)) {
    for (y in 1:length(Parameters_2)) {
      for (q in 1:length(Species_1)) {
        Para_1 <- Parameters_1[x]
        Para_2 <- Parameters_2[y]
        Species_q <- Species_1[q]
        if(!(Species_q=="Zebrafish"&(Para_1%in%c("MMR","AS","FAS")|Para_2%in%c("MMR","AS","FAS")))){
          ifelse(x==1&y==1&q==1,nr <-1,nr<-nr+1)
          ifelse(Species_q=="Mean", Data_1_v <- Results_fit_p[Results_fit_p$Parameter==Para_1,]$Slope_fit,Data_1_v <-Results_fit_p[Results_fit_p$Parameter==Para_1&Results_fit_p$Species==Species_q,]$Slope_fit) 
          ifelse(Species_q=="Mean",Data_2_v <- Results_fit_p[Results_fit_p$Parameter==Para_2,]$Slope_fit,Data_2_v <-Results_fit_p[Results_fit_p$Parameter==Para_2&Results_fit_p$Species==Species_q,]$Slope_fit)
          Cor_ID <- cor.test(Data_1_v,Data_2_v)
          bcor_table_species[nr,1] <- Species_q
          bcor_table_species[nr,2] <- Para_1
          bcor_table_species[nr,3] <- Para_2
          bcor_table_species[nr,4] <- round(Cor_ID$estimate,3)
          bcor_table_species[nr,5] <- round(Cor_ID$p.value,3)
          bcor_table_species[nr,6] <- paste(Species_q," ",Para_1," vs ",Para_2," p = ",round(Cor_ID$p.value,3),sep = "")
        }
      }
    }
  }

Results_fit3$Species2 <- factor(Results_fit3$Species2,levels=c("Brown trout","Anemonefish","Cunner","Damselfish","Guppy","Rainbow trout","Zebrafish", "All species"))
Results_fit3$Plot_treatment <- as.character(Results_fit3$Plot_treatment)


bSMRvbAGR_plot <- ggplot(data = Results_fit3, aes(x=Slope_AGR_fit,y=Slope_SMR_fit))+
  geom_point(aes(col=Species, shape=Treatment ),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.5, method = "lm", se=F )+
   theme_bw()+
  theme(plot.background = element_blank(),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(size=14, margin = margin(t = 8, r = 0, b = -9, l = 0))) +
  labs(y=expression(SMR~scaling ~ exponent ~ (italic(b)) ), x=expression(GR~scaling ~ exponent ~ (italic(b)) ))+
  facet_wrap(~Species2,scales = "free")+
  scale_shape_manual(name=NULL, values = c("Control"=16, "m"=16,"f"=17,"low"=16,"L"=16,"10Cegg"=16,"M"=15,"14Cegg"=15,"H"=17,"high"=17,"14Cyolk"=17))+
  theme(strip.background = element_blank(), strip.text.x = element_text(size=12, colour="black", margin=margin(1,0,0.5,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 

bcor_table_species[bcor_table_species$Par_1=="SMR",]

unique(Results_fit3$Treatment)

bMMRvbAGR_plot <- ggplot(data = Results_fit3[!Results_fit$Species=="Zebrafish",], aes(x=Slope_AGR_fit,y=Slope_MMR_fit))+
  geom_point(aes(col=Species, shape=Treatment),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.5, method = "lm", se=F )+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(size=14, margin = margin(t = 8, r = 0, b = -9, l = 0))) +
  labs(y=expression(MMR~scaling ~ exponent ~ (italic(b)) ), x=expression(AGR~scaling ~ exponent ~ (italic(b)) ))+
  facet_wrap(~Species2,scales = "free")+
  scale_shape_manual(name=NULL, values = c("Control"=16, "m"=16,"f"=17,"low"=16,"L"=16,"10Cegg"=16,"M"=15,"14Cegg"=15,"H"=17,"high"=17,"14Cyolk"=17))+
  theme(strip.background = element_blank(), strip.text.x = element_text(size=12, colour="black", margin=margin(1,0,0.5,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 
bcor_table_species[bcor_table_species$Par_1=="MMR",]


bFASvbAGR_plot <- ggplot(data = Results_fit3[!Results_fit$Species=="Zebrafish",], aes(x=Slope_AGR_fit,y=Slope_FAS_fit))+
  geom_point(aes(col=Species, shape=Treatment),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.5, method = "lm", se=F )+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(size=14, margin = margin(t = 8, r = 0, b = -9, l = 0))) +
  labs(y=expression(FAS~scaling ~ exponent ~ (italic(b)) ), x=expression(GR~scaling ~ exponent ~ (italic(b)) ))+
  facet_wrap(~Species2,scales = "free")+
  scale_shape_manual(name=NULL, values = c("Control"=16, "m"=16,"f"=17,"low"=16,"L"=16,"10Cegg"=16,"M"=15,"14Cegg"=15,"H"=17,"high"=17,"14Cyolk"=17))+
    theme(strip.background = element_blank(), strip.text.x = element_text(size=12, colour="black", margin=margin(1,0,0.5,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 
bcor_table_species[bcor_table_species$Par_1=="FAS",]




bASvbAGR_plot <- ggplot(data = Results_fit3[!Results_fit$Species=="Zebrafish",], aes(x=Slope_AGR_fit,y=Slope_AS_fit))+
  geom_point(aes(col=Species, shape=Treatment),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.5, method = "lm", se=F )+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(size=14, margin = margin(t = 8, r = 0, b = -9, l = 0))) +
  labs(y=expression(AS~scaling ~ exponent ~ (italic(b)) ), x=expression(GR~scaling ~ exponent ~ (italic(b)) ))+
  facet_wrap(~Species2,scales = "free")+
  scale_shape_manual(name=NULL, values = c("Control"=16, "m"=16,"f"=17,"low"=16,"L"=16,"10Cegg"=16,"M"=15,"14Cegg"=15,"H"=17,"high"=17,"14Cyolk"=17))+
  theme(strip.background = element_blank(), strip.text.x = element_text(size=12, colour="black", margin=margin(1,0,0.5,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 
bcor_table_species[bcor_table_species$Par_1=="AS",]
C_AS_AGR <- lm(data=Results_fit[Results_fit$Species=="Clownfish",][-15,],Slope_AS_fit~Slope_AGR_fit)
check_outliers(C_AS_AGR)
check_model(C_AS_AGR)
summary(C_AS_AGR)
cor.test(Results_fit[Results_fit$Species=="Clownfish",][-15,]$Slope_AS_fit,Results_fit[Results_fit$Species=="Clownfish",][-15,]$Slope_AGR_fit)
Results_fit[Results_fit$Species=="Clownfish",][15,]$Slope_AS_fit
Scaling_exponents_all[Scaling_exponents_all$Parameter=="AGR",]
lm(data = Data_fit, log10(GE_fit_TempAdj)~logMass)
ggplot()+geom_point(data = Data_fit, aes(x=logMass, y=log10(GE_fit_TempAdj),col=Species))+geom_abline(intercept=-0.0478,slope=-0.238)
min(na.omit(Results_fit[Results_fit$Species=="Brown_trout",]$Slope_SMR_fit))


bASvbAGR_plot <- ggplot(data = Results_fit3, aes(x=Slope_AGR_fit,y=Slope_GE_fit))+
  geom_point(aes(col=Species, shape=Plot_treatment),alpha=0.8)+
  geom_line(aes(col="Mean"), stat = "smooth", method = "lm", alpha = 1, linewidth = 0.5, method = "lm", se=F )+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(color="black", linewidth = 0.3))+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position="none")+
  force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))+
  #theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.background=element_blank())+
  theme(axis.ticks.length=unit(1,"mm"), axis.ticks.y = element_line(colour = "black", size = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3))+
  theme(aspect.ratio=1/1)+
  theme(axis.text.y = element_text(size=12, margin = margin(r = 2), colour="black"), axis.text.x = element_text(size=12,margin = margin(t = 3), colour="black"))+
  theme(axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(size=14, margin = margin(t = 8, r = 0, b = -9, l = 0))) +
  labs(y=expression(GE~scaling ~ exponent ~ (italic(b)) ), x=expression(AGR~scaling ~ exponent ~ (italic(b)) ))+
  facet_wrap(~Species2,scales = "free")+
  theme(strip.background = element_blank(), strip.text.x = element_text(size=12, colour="black", margin=margin(1,0,0.5,0, "mm")))+
  scale_colour_manual(name = NULL, values = c("Mean" = "black", "Brown_trout" = "#DA1F26",  "Chromis" = "#F26739", "Clownfish" = "#00BE00" , "Cunner" = "#FFBE00","Guppy"="#00C7FF","Rainbow_trout"="#0000FF","Zebrafish"="#B400FF")) 
