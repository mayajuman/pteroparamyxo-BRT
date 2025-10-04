## pteropodidae-paramyxovirus associations: virus data collection and assembly
## maya juman + dan becker
## updated 01/23/24

## clean up space
rm(list=ls()) 
graphics.off()

## set up

## libraries
library(vroom)
library(tidyverse)
library(stringr)
library(plyr)
library(R.utils)
library(ape)
library(caper)
library(dplyr)

## lit review raw virus data

setwd("~/Documents/PhD/1. pteropodidae")
lit <- read.csv("pteroparamyxo_simple.csv") #simplified file from Juman et al. 2025 pteroparamyxo database
lit$host.species[which(lit$host.species == "labiatus minor")] <- "labiatus" #ignore subspecies
lit <- lit %>% drop_na(host.genus, host.species) %>% filter(!host.species == "sp.") #remove NAs and "sp."s
lit$host <- with(lit, paste0(host.genus," ",host.species)) #create host name variable

length(unique(lit$host)) #60 species

## virion data (accessed 3 august 2023)

#virion <- vroom("Virion.csv.gz")
#pteroparamyxos <- virion %>% filter(VirusFamily == "paramyxoviridae", HostFamily == "pteropodidae")

#write.csv(pteroparamyxos, file="pteroparamyxos.csv", row.names = FALSE)
virion <- read.csv("virion.csv")

## dump GLOBI b/c no prep information or citation info for cross-ref, also redundant w/ VIRION
## only pteropus rayneri is in GLOBI but none of the other databases... cannot find citation for this host-virus association
## remove questionable observations in VIRION: 
## eonycteris spelaea NiV isolation that was mis-cited from Calisher et al. 2006
## cynopterus brachyotis NiV isolation that was attributed to WHO but not findable anywhere online

virion <- virion %>% drop_na(Host) %>% 
  filter(!Database == "GLOBI", !DetectionMethod == "Not specified") %>% 
  filter((AssocID != 174639) %>% replace_na(TRUE)) %>% 
  filter((AssocID != 174030) %>% replace_na(TRUE))

## constructing host-virus association file

####creating host-virus association file from literature search + VIRION

data <- data.frame(matrix(ncol = 0, nrow = 60))
data$host <- unique(lit$host)
data$paramyxovirus.PCR.positive <- NA
data$paramyxovirus.seropositive <- NA
data$paramyxovirus.isolated <- NA
data$orthoparamyxovirus.PCR.positive <- NA
data$orthoparamyxovirus.seropositive <- NA
data$orthoparamyxovirus.isolated <- NA
data$rubulavirus.PCR.positive <- NA
data$rubulavirus.seropositive <- NA
data$rubulavirus.isolated <- NA
data$notes <- NA

for (i in 1:nrow(data)) {
  data$paramyxovirus.PCR.positive[i] <- 
    sum(lit[which(lit$host == data$host[i]),]$PCR.positive,na.rm = TRUE)
  data$paramyxovirus.seropositive[i] <- 
    sum(lit[which(lit$host == data$host[i]),]$sero.positive,na.rm = TRUE)
  data$paramyxovirus.isolated[i] <- 
    sum(lit[which(lit$host == data$host[i]),]$isolation.positive,na.rm = TRUE)
  data$orthoparamyxovirus.PCR.positive[i] <-
    sum(lit[which(lit$host == data$host[i] & 
                    lit$virus.subfamily == "orthoparamyxovirinae"),]$PCR.positive, na.rm = TRUE)
  data$orthoparamyxovirus.seropositive[i] <- 
    sum(lit[which(lit$host == data$host[i] & 
                    lit$virus.subfamily == "orthoparamyxovirinae"),]$sero.positive, na.rm = TRUE)
  data$orthoparamyxovirus.isolated[i] <- 
    sum(lit[which(lit$host == data$host[i] & 
                    lit$virus.subfamily == "orthoparamyxovirinae"),]$isolation.positive, na.rm = TRUE)
  data$rubulavirus.PCR.positive[i] <- 
    sum(lit[which(lit$host == data$host[i] & 
                    lit$virus.subfamily == "rubulavirinae"),]$PCR.positive, na.rm = TRUE)
  data$rubulavirus.seropositive[i] <- 
    sum(lit[which(lit$host == data$host[i] & 
                    lit$virus.subfamily == "rubulavirinae"),]$sero.positive, na.rm = TRUE)
  data$rubulavirus.isolated[i] <- 
    sum(lit[which(lit$host == data$host[i] & 
                    lit$virus.subfamily == "rubulavirinae"),]$isolation.positive, na.rm = TRUE)
}

data$paramyxovirus.PCR.positive[which(data$paramyxovirus.PCR.positive>1)] <- 1
data$paramyxovirus.seropositive[which(data$paramyxovirus.seropositive>1)] <- 1
data$paramyxovirus.isolated[which(data$paramyxovirus.isolated>1)] <- 1
data$orthoparamyxovirus.PCR.positive[which(data$orthoparamyxovirus.PCR.positive>1)] <- 1
data$orthoparamyxovirus.seropositive[which(data$orthoparamyxovirus.seropositive>1)] <- 1
data$orthoparamyxovirus.isolated[which(data$orthoparamyxovirus.isolated>1)] <- 1
data$rubulavirus.PCR.positive[which(data$rubulavirus.PCR.positive>1)] <- 1
data$rubulavirus.seropositive[which(data$rubulavirus.seropositive>1)] <- 1
data$rubulavirus.isolated[which(data$rubulavirus.isolated>1)] <- 1

extras <- virion %>% filter(!Host %in% data$host) %>% 
  filter(!Host == "pteropus melanotus") %>% distinct(Host)
extras <- list(extras)
data[61:66, 1] <- extras[[1]]
data[61:66,c(2:10)] <- 0

virionchecks <- data.frame(matrix(ncol = ncol(virion), nrow = 0))

for (i in 1:nrow(data)) {
  if (data[i,2] == 0 & data$host[i] %in% virion$Host & 
      any(virion$DetectionMethod[which(virion$Host == data$host[i])] == "PCR/Sequencing")) {
    data[i,2] <- 1
    data$notes[i] <- "supplemented from VIRION"
    virionchecks <- rbind(virionchecks, 
                          virion[which(virion$Host == data$host[i] & 
                                         virion$DetectionMethod == "PCR/Sequencing"),])
  } 
  if (data[i,3] == 0 & data$host[i] %in% virion$Host & 
      any(virion$DetectionMethod[which(virion$Host == data$host[i])] == "Antibodies")) {
    data[i,3] <- 1
    data$notes[i] <- "supplemented from VIRION"
    virionchecks <- rbind(virionchecks, 
                          virion[which(virion$Host == data$host[i] & 
                                         virion$DetectionMethod == "Antibodies"),])
  } 
  if (data[i,4] == 0 & data$host[i] %in% virion$Host & 
      any(virion$DetectionMethod[which(virion$Host == data$host[i])] == "Isolation/Observation")) {
    data[i,4] <- 1
    data$notes[i] <- "supplemented from VIRION"
    virionchecks <- rbind(virionchecks, 
                          virion[which(virion$Host == data$host[i] & 
                                         virion$DetectionMethod == "Isolation/Observation"),])
  }
  if (data[i,5] == 0 & data$host[i] %in% virion$Host & 
      any(virion$DetectionMethod[which(virion$Host == data$host[i] & 
                                       virion$VirusGenus == "henipavirus")] == "PCR/Sequencing")) {
    data[i,5] <- 1
    data$notes[i] <- "supplemented from VIRION"
    virionchecks <- rbind(virionchecks, 
                          virion[which(virion$Host == data$host[i] & 
                                         virion$DetectionMethod == "PCR/Sequencing" & 
                                         virion$VirusGenus == "henipavirus"),])
  } 
  if (data[i,6] == 0 & data$host[i] %in% virion$Host & 
      any(virion$DetectionMethod[which(virion$Host == data$host[i] & 
                                       virion$VirusGenus == "henipavirus")] == "Antibodies")) {
    data[i,6] <- 1
    data$notes[i] <- "supplemented from VIRION"
    virionchecks <- rbind(virionchecks, 
                          virion[which(virion$Host == data$host[i] & 
                                         virion$DetectionMethod == "Antibodies" & 
                                         virion$VirusGenus == "henipavirus"),])
  } 
  if (data[i,7] == 0 & data$host[i] %in% virion$Host & 
      any(virion$DetectionMethod[which(virion$Host == data$host[i] & 
                                       virion$VirusGenus == "henipavirus")] == "Isolation/Observation")) {
    data[i,7] <- 1
    data$notes[i] <- "supplemented from VIRION"
    virionchecks <- rbind(virionchecks, 
                          virion[which(virion$Host == data$host[i] & 
                                         virion$DetectionMethod == "Isolation/Observation" & 
                                         virion$VirusGenus == "henipavirus"),])
  }
  if (data[i,8] == 0 & data$host[i] %in% virion$Host & 
      any(virion$DetectionMethod[which(virion$Host == data$host[i] & 
                                       virion$VirusGenus == "pararubulavirus")] == "PCR/Sequencing")) {
    data[i,8] <- 1
    data$notes[i] <- "supplemented from VIRION"
    virionchecks <- rbind(virionchecks, 
                          virion[which(virion$Host == data$host[i] & 
                                         virion$DetectionMethod == "PCR/Sequencing" & 
                                         virion$VirusGenus == "pararubulavirus"),])
  } 
  if (data[i,9] == 0 & data$host[i] %in% virion$Host & 
      any(virion$DetectionMethod[which(virion$Host == data$host[i] & 
                                       virion$VirusGenus == "pararubulavirus")] == "Antibodies")) {
    data[i,9] <- 1
    data$notes[i] <- "supplemented from VIRION"
    virionchecks <- rbind(virionchecks, 
                          virion[which(virion$Host == data$host[i] & 
                                         virion$DetectionMethod == "Antibodies" & 
                                         virion$VirusGenus == "pararubulavirus"),])
  } 
  if (data[i,10] == 0 & data$host[i] %in% virion$Host & 
      any(virion$DetectionMethod[which(virion$Host == data$host[i] & 
                                       virion$VirusGenus == "pararubulavirus")] == "Isolation/Observation")) {
    data[i,10] <- 1
    data$notes[i] <- "supplemented from VIRION"
    virionchecks <- rbind(virionchecks, 
                          virion[which(virion$Host == data$host[i] & 
                                         virion$DetectionMethod == "Isolation/Observation" & 
                                         virion$VirusGenus == "pararubulavirus"),])
  }
}

## 17 checks to make
virionchecks <- virionchecks %>% group_by_all %>% distinct

## 11 changes by VIRION
virionchecks %>% group_by(Host) %>% distinct(Host)

#write.csv(set, file="pteropodidae paramyxo associations.csv", row.names = FALSE)

## matching to taxonomy

## load in bat phylogeny and taxonomy
setwd("~/Documents/PhD/1. pteropodidae/phylo")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)

## trim to family
taxa=taxa[taxa$fam=="PTEROPODIDAE",]

## trim tree
tree=keep.tip(tree,taxa$tiplabel)

## tip
taxa$tip=taxa$Species_Name
taxa=taxa[c("tip","clade","fam","gen","tiplabel")]

## make tip for data
data$tip=gsub(" ","_",capitalize(data$host))

## are all data in taxa?
setdiff(data$tip,taxa$tip)

## fix data
data$tip=revalue(data$tip,
                 c("Dobsonia_magna"="Dobsonia_moluccensis",
                   "Pteropus_medius"="Pteropus_giganteus",
                   "Lissonycteris_angolensis"="Myonycteris_angolensis",
                   "Pteropus_melanotus_natalis"="Pteropus_melanotus",
                   "Pteropus_seychellensis_comorensis"="Pteropus_seychellensis",
                   "Dobsonia_andersoni"="Dobsonia_anderseni"))

## are all data in taxa?
setdiff(data$tip,taxa$tip)

## re-aggregate data
set=data
set$notes=NULL
set$host=NULL

## aggregate all
set=set %>%
  group_by(tip) %>% 
  summarise_each(funs(sum))

## replace 2s with 1s for doubled species
set[which(set$tip == "Pteropus_giganteus"),c(2:7)] <- 1
set$orthoparamyxovirus.seropositive[which(set$tip == "Dobsonia_moluccensis")] <- 1
set$paramyxovirus.seropositive[which(set$tip == "Dobsonia_moluccensis")] <- 1

## drop missing species
drop=setdiff(data$tip,taxa$tip)
set=set[!set$tip%in%drop,]

## clean tree tips
x=strsplit(tree$tip.label,"_")
tree$tip.label=sapply(x,function(x) paste(x[1],"_",x[2],sep=""))

## merge set with all pteropids
set=merge(set,taxa,by="tip",all.y=T)

## save
set$label=set$tip

## replace NAs with 0 (pseudoabsences)
set[is.na(set)]=0

## merge into single comparative data frame
cdata=comparative.data(phy=tree,data=set,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

setwd("~/Documents/PhD/pteropodidae")
#write.csv(set, file="virus data.csv", row.names = FALSE)
