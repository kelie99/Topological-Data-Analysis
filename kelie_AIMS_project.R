###############################################################################
##############################################################################
########### PERSISTENT HOMOLOGY FOR BREAST CANCER(WISCONSIN) DATASET ########
###############################################################################
##############################################################################

#Import necessary modules.
library(remotes)
library(ripserr)
library(TDA)
library(ggplot2)
library(gplots)
library(tidyverse)
library(rdist)
library(TDApplied)

#Load the breast cancer data.
data <- read.csv("data.csv")

#Split into Benign and Malignant dataset
data_mali <- subset(data, diagnosis =='M')
data_beng <- subset(data, diagnosis == 'B')

#Extract relevant columns in both dataset. 
data_Mali<- data_mali[,c(3:32)]
data_Beng <- data_beng[,c(3:32)]

#Implement the persistence algorithm on the datasets
#For the malignant dataset
RipsDiag_m<- ripserr::vietoris_rips(dataset = data_Mali,
                                    threshold = max(dist(data_Mali)),
                                    dim=2, standardize = TRUE)
Diag_m <- as.matrix(RipsDiag_m)

#For the benign dataset
RipsDiag_b<- ripserr::vietoris_rips(dataset = data_Beng,
                                    threshold = max(dist(data_Beng)), 
                                    dim=2, standardize = TRUE)
Diag_b <- as.matrix(RipsDiag_b)

#Plots the persistent barcode of both groups.
par(mfrow=c(1,2))
plot.diagram(Diag_m, barcode = TRUE, main='Malignant persisence barcode')
plot.diagram(Diag_b, barcode = TRUE, main = 'Benign Persistence barcode')


#Plots the 1-dimensional persistence barcode of malignant versus benign
par(mfrow=c(1,2))
plot.diagram(Diag_b[c(357:620),], barcode = TRUE, main="", col='red')
plot.diagram(Diag_m[c(212:340),], barcode = TRUE, main="", col='red')

par(mfrow=c(1,2))
plot.diagram(Diag_b[c(357:620),], barcode = TRUE, main="", col='blue')
plot.diagram(Diag_m[c(341:381),], barcode = TRUE, main="", col='blue')

#Boostrapping on both groups of data to generate 3 subgroups
set.seed(15)
subgroup_m <- list()
for (i in 1:3){
  subgroup_m[[i]] <- data_Mali[sample(nrow(data_Mali), size = 210),]
}

#Compute the persistence on subgroup list of malignant
RipsDiag_sbgm <- list()
for (i in 1:3){
  RipsDiag_sbgm[[i]] <- ripserr::vietoris_rips(subgroup_m[[i]], dim=2,standardize = TRUE ,threshold = max(dist(subgroup_m[[i]])))
}

#plot persistence barcode of all subgroups
par(mfrow = c(1,3))
for (i in 1:3){
  plot.diagram(as.matrix(RipsDiag_sbgm[[i]]), barcode = TRUE)
}


set.seed(15)
#similar code for benign
subgroup_b <- list()
for (i in 1:3){
  subgroup_b[[i]] <- data_Beng[sample(nrow(data_Beng), size = 210),]
  }
RipsDiag_sbgb <- list()
for (i in 1:3){ 
  RipsDiag_sbgb[[i]] <- ripserr::vietoris_rips(subgroup_b[[i]], dim=2,standardize = TRUE ,threshold = max(dist(subgroup_b[[i]])))
}

#plot persistence barcode of all subgroups
par(mfrow = c(1,3))
for (i in 1:3){
  plot.diagram(as.matrix(RipsDiag_sbgb[[i]]), barcode = TRUE)
}



#############################################################
######################### Some statistics ####################
##############################################################

#Boxplot of 2 dimensional holes of benign and malignant
pers_begnin_b2 <- subset(RipsDiag_sbgb[[1]], dimension == 2)
pers_begnin_b2$diagnosis <- 'B'
pers_begnin_b2$length <- (pers_begnin_b2$death - pers_begnin_b2$birth)


pers_mali_b2 <- subset(RipsDiag_sbgm[[1]], dimension ==2)
pers_mali_b2$diagnosis <- 'M'
pers_mali_b2$length <- (pers_mali_b2$death - pers_mali_b2$birth)

pers_data_b2 <- rbind(pers_begnin_b2, pers_mali_b2)

#Draw a boxplot of group/start times of barcode 
ggplot(data= pers_data_b2, aes(x=diagnosis , y= birth, fill =diagnosis))+stat_boxplot(geom="errorbar")+ geom_boxplot()
ggplot(data= pers_data_b2, aes(x=diagnosis , y= length, fill =diagnosis))+stat_boxplot(geom="errorbar")+ geom_boxplot()


#some statistics on 1 homology group of malignant and Begnin
pers_begnin_b1 <- subset(RipsDiag_sbgb[[1]], dimension == 1)
pers_begnin_b1$diagnosis <- 'B'
pers_begnin_b1$length <- (pers_begnin_b1$death - pers_begnin_b1$birth)

pers_mali_b1 <- subset(RipsDiag_sbgm[[1]], dimension ==1)
pers_mali_b1$diagnosis <- 'M'
pers_mali_b1$length <- (pers_mali_b1$death - pers_mali_b1$birth)

pers_data_b1 <- rbind(pers_begnin_b1, pers_mali_b1)

#Draw a boxplot of group/start times of barcode 
ggplot(data= pers_data_b1, aes(x=diagnosis , y= birth, fill =diagnosis))+stat_boxplot(geom="errorbar")+ geom_boxplot()
ggplot(data= pers_data_b1, aes(x=diagnosis , y= length, fill =diagnosis ))+stat_boxplot(geom="errorbar")+ geom_boxplot()



#compute stats each group on birth(x) and death(y) for B
names <- c("mean", "median", "sd")

stat_birth_b<- c(mean(pers_mali_b2$birth),median(pers_mali_b2$birth),sd(pers_mali_b2$birth))
stat_death_b<- c(mean(pers_mali_b2$death),median(pers_mali_b2$death),sd(pers_mali_b2$death))
stat_length_b<- c(mean(pers_mali_b2$length), median(pers_mali_b2$length), sd(pers_mali_b2$length))
stat_diffdeath_b <- c(mean(max(pers_mali_b2$death)-(pers_mali_b2$death)), median(max(pers_mali_b2$death)-(pers_mali_b2$death)), sd(max(pers_mali_b2$death)-(pers_mali_b2$death))  )

stat_mali_B2 <- data.frame(names,stat_birth_b, stat_death_b, stat_length_b, stat_diffdeath_b)


stat_birth_m<- c(mean(pers_begnin_b2$birth),median(pers_begnin_b2$birth),sd(pers_begnin_b2$birth))
stat_death_m<- c(mean(pers_begnin_b2$death),median(pers_begnin_b2$death),sd(pers_begnin_b2$death))
stat_length_m<- c(mean(pers_begnin_b2$length), median(pers_begnin_b2$length), sd(pers_begnin_b2$length))
stat_diffdeath_m<- c(mean(max(pers_begnin_b2$death)-(pers_begnin_b2$death)), median(max(pers_begnin_b2$death)-(pers_begnin_b2$death)), sd(max(pers_begnin_b2$death)-(pers_begnin_b2$death))  )

stat_begnin_b2 <- data.frame(names,stat_birth_m, stat_death_m, stat_length_m, stat_diffdeath_m)

##############################
# END ANALYSIS
