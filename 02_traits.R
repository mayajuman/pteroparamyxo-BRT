## pteropodidae-paramyxovirus: trait merging and cleaning
## maya juman + dan becker
## updated 01/23/24

## clean up space
rm(list=ls()) 
graphics.off()

## set up

library(tidyverse)
library(dplyr)
library(readxl)
library(rredlist)
library(R.utils)
library(plyr)
library(terra)
library(reshape2)
library(easyPubMed)
library(phyloregion)
library(fastDummies)
library(ape)

## load host-virus associations

setwd("~/Documents/PhD/1. pteropodidae")
set <- read.csv("virus data.csv")

## load databases

setwd("~/Documents/PhD/1. pteropodidae/databases")

## guy et al. 2020
guy <- read.csv("Guy et al. 2021.csv")
guy <- guy[which(guy$Family == "Pteropodidae"),]

## crane et al. 2022
wings <- read_excel("mam12270-sup-0004-appendixs7-s9.xlsx", sheet = "Appendix S8_species level data")
wings <- wings %>% filter(familyName == "PTEROPODIDAE") %>% 
  dplyr::select(scientificName.iucn.sub, wingspan_cm,
                liftingSurfaceArea_cm2, totalAspectRatio, totalRelativeLoading)
wings[wings == "NA"] <- NA
wings$tip=gsub(" ","_",capitalize(wings$scientificName.iucn.sub))
wings$scientificName.iucn.sub <- NULL

## Cosentino et al. 2023
setwd("~/Documents/PhD/1. pteropodidae/databases/AfroBaT pre-imputation dataset")
ABTmorph <- read.csv("AfroBaT_morphology.csv")
ABTrepro <- read.csv("AfroBaT_reproduction.csv")
ABTlh <- read.csv("AfroBaT_life_history.csv")
ABT <- merge(ABTmorph,ABTrepro, all.x = TRUE, by="SpeciesEpithet")
ABT <- merge(ABT, ABTlh, by="SpeciesEpithet")
rm(ABTmorph,ABTrepro, ABTlh)

## COMBINE (soria et al. 2021)
setwd("~/Documents/PhD/1. pteropodidae/databases/COMBINE_archives")
combine <- read.csv("trait_data_reported.csv")
combine <- combine[which(combine$family == "Pteropodidae"),]

guy$Binomial <- str_to_title(gsub(" ", "_", guy$Binomial))
combine$Binomial <- with(combine, paste0(genus,"_",species))
combine <- combine %>% dplyr::select(Binomial, everything())

guy$Binomial=revalue(guy$Binomial,
                     c("Boneia_bidens"="Rousettus_bidens",
                       "Stenonycteris_lanosus"="Rousettus_lanosus"))

guy %>% filter(!Binomial %in% set$tip)

combine$Binomial=revalue(combine$Binomial,
                         c("Boneia_bidens"="Rousettus_bidens",
                           "Lissonycteris_angolensis"="Myonycteris_angolensis",
                           "Pteropus_leucopterus"="Desmalopex_leucopterus",
                           "Scotonycteris_ophiodon"="Casinycteris_ophiodon"))

combine %>% filter(!Binomial %in% set$tip)
set %>% filter(!tip %in% combine$Binomial)
setdiff(combine$Binomial, set$tip)
setdiff(set$tip, combine$Binomial)

guy <- guy[,c(4,12:38)]
combine <- combine[,-c(2:7)]

traits <- merge(set, guy, all.x = TRUE, by.x="tip", by.y = "Binomial")
traits <- merge(traits, combine, all.x = TRUE, by.x="tip", by.y = "Binomial")

setwd("~/Documents/PhD/1. pteropodidae")
IUCN <- read.csv("IUCNchecks.csv")
IUCN <- IUCN %>% dplyr::select(tip, StatusM, TrendM)

traits <- merge(IUCN, traits)
traits <- traits %>% dplyr::select(-Status, -Trend)
traits <- traits %>% relocate(StatusM, TrendM, .after = last_col())

## remove redundant or logged traits
traits <- traits %>% 
  dplyr::select(-log_AR, -log_Forearm, -log_Mass, -sqrt_Rarea, 
                -log_Median, -Mass, -Forearm, -Pulse, -Lifespan,
                -sin_Citations, -Torpor)

## supplementing from AfroBaT

ABT <- ABT %>% filter(Family.x == "Pteropodidae") %>% dplyr::select(SpeciesEpithet,
                                                                    W_loadingU_mean,
                                                                    W_aspect_ratioU_mean,
                                                                    Colony_SizeN_mean,
                                                                    Migratory,
                                                                    BodyMass_Neonatal_mean,
                                                                    BodyMass_Weaning_mean,
                                                                    Lifespan_Captivity_max,
                                                                    Lifespan_Wild_max,
                                                                    SexMaturity_AgeM_mean,
                                                                    SexMaturityF_Age_mean,
                                                                    SexMaturity_AgeU_mean,
                                                                    Active_Gestation_Length_mean,
                                                                    Litter_Size_mean,
                                                                    Num_Litter_year_mean,
                                                                    Births_Interval_mean,
                                                                    Weaning_Length_mean,
                                                                    Dispersal_DistU_mean,
                                                                    Hibernation) %>%
  dplyr::rename(tip = SpeciesEpithet)
ABT <- ABT[colSums(!is.na(ABT)) > 0]

## converting units
traits$Migration[which(traits$Migration == "N")] <- 0
traits$Migration[which(traits$Migration == "Y")] <- 1
traits$Migration <- as.integer(traits$Migration)
ABT$Lifespan_Captivity_max <- ABT$Lifespan_Captivity_max*365

## coalesce

traits <- left_join(traits, ABT, by = "tip")

traits <- traits %>%
  mutate(Group.Size = coalesce(Group.Size, Colony_SizeN_mean),
         Migration = coalesce(Migration, Migratory),
         neonate_mass_g = coalesce(neonate_mass_g, BodyMass_Neonatal_mean),
         max_longevity_d = coalesce(max_longevity_d, Lifespan_Captivity_max),
         maturity_d = coalesce(maturity_d, SexMaturity_AgeU_mean),
         male_maturity_d = coalesce(male_maturity_d, SexMaturity_AgeM_mean),
         female_maturity_d = coalesce(female_maturity_d, SexMaturityF_Age_mean),
         gestation_length_d = coalesce(gestation_length_d, Active_Gestation_Length_mean),
         litter_size_n = coalesce(litter_size_n, Litter_Size_mean),
         litters_per_year_n = coalesce(litters_per_year_n, Num_Litter_year_mean),
         interbirth_interval_d = coalesce(interbirth_interval_d, Births_Interval_mean),
         weaning_age_d = coalesce(weaning_age_d, Weaning_Length_mean),
         hibernation_torpor = coalesce(hibernation_torpor, Hibernation))

traits$Colony_SizeN_mean=NULL
traits$Migratory=NULL
traits$BodyMass_Neonatal_mean=NULL
traits$Lifespan_Captivity_max=NULL
traits$SexMaturity_AgeM_mean=NULL
traits$SexMaturity_AgeU_mean=NULL
traits$SexMaturityF_Age_mean=NULL
traits$Active_Gestation_Length_mean=NULL
traits$Litter_Size_mean=NULL
traits$Num_Litter_year_mean=NULL
traits$Births_Interval_mean=NULL
traits$Weaning_Length_mean=NULL
traits$Hibernation=NULL

#minitraits <- traits %>% select(tip, adult_mass_g, adult_forearm_length_mm, adult_body_length_mm)
#write.csv(minitraits, file="HMWgapsblank.csv", row.names = FALSE)

setwd("~/Documents/PhD/1. pteropodidae")
HMW <- read.csv("HMWgaps.csv")
HMW <- HMW %>% mutate(HMWlength = rowMeans(cbind(HMWlengthmin, HMWlengthmax), na.rm=T),
                      HMWforearm = rowMeans(cbind(HMWforearmmin, HMWforearmmax), na.rm=T),
                      HMWmass = rowMeans(cbind(HMWmassmin, HMWmassmax), na.rm=T),
                      ear_mm = rowMeans(cbind(ear_min, ear_max), na.rm=T),
                      hindfoot_mm = rowMeans(cbind(hindfoot_min, hindfoot_max), na.rm=T),
                      tail_mm = rowMeans(cbind(tail_min, tail_max), na.rm=T)) %>%
  dplyr::select(tip, HMWlength, adult_body_length_mm, HMWforearm, adult_forearm_length_mm, 
         HMWmass, adult_mass_g, ear_mm, hindfoot_mm, tail_mm)

HMW[HMW == "NaN"] <- NA

HMW <- HMW %>%
  mutate(Adult_length_mm = coalesce(HMWlength, adult_body_length_mm),
         Adult_forearm_mm = coalesce(HMWforearm, adult_forearm_length_mm),
         Adult_mass_g = coalesce(HMWmass, adult_mass_g)) %>%
  dplyr::select(tip, ear_mm, hindfoot_mm, tail_mm, Adult_length_mm,
         Adult_forearm_mm, Adult_mass_g)

traits <- left_join(traits, HMW, by = "tip")
traits$adult_mass_g=NULL
traits$adult_forearm_length_mm=NULL
traits$adult_body_length_mm=NULL

## add wing traits

traits <- left_join(traits, wings, by = "tip")

## remove less complete wing and IUCN traits
traits <- traits %>% dplyr::select(-RWL, -AR, -Range_Area, -Median_Y_Lat, -Median_X_Long)

## add IUCN range data

setwd("~/Documents/PhD/1. pteropodidae/databases/IUCN")
IUCN <- vect('data_0.shp')

traits$MeanLat <- NA
traits$MeanLong <- NA
traits$RangeArea <- NA

IUCNsp <- unique(IUCN$SCI_NAME)

for (i in 1:nrow(traits)) {
  name <- gsub("_"," ",traits$tip[i])
  if (name%in%IUCNsp) {
  traits$MeanLat[i] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == name]))[c(3,4)])
  traits$MeanLong[i] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == name]))[c(1,2)])
  traits$RangeArea[i] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == name]), unit="km")
  print(paste(i,"/",nrow(traits)))
  }
}

## specific cases of IUCN taxonomic mismatches

traits$MeanLat[which(traits$tip == "Pteropus_giganteus")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Pteropus medius"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Pteropus_giganteus")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Pteropus medius"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Pteropus_giganteus")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Pteropus medius"]), unit="km")

traits$MeanLat[which(traits$tip == "Rousettus_bidens")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Boneia bidens"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Rousettus_bidens")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Boneia bidens"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Rousettus_bidens")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Boneia bidens"]), unit="km")

traits$MeanLat[which(traits$tip == "Rousettus_celebensis")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Pilonycteris celebensis"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Rousettus_celebensis")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Pilonycteris celebensis"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Rousettus_celebensis")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Pilonycteris celebensis"]), unit="km")

traits$MeanLat[which(traits$tip == "Rousettus_lanosus")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Stenonycteris lanosus"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Rousettus_lanosus")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Stenonycteris lanosus"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Rousettus_lanosus")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Stenonycteris lanosus"]), unit="km")

traits$MeanLat[which(traits$tip == "Epomops_dobsonii")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Epomophorus dobsonii"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Epomops_dobsonii")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Epomophorus dobsonii"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Epomops_dobsonii")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Epomophorus dobsonii"]), unit="km")

traits$MeanLat[which(traits$tip == "Melonycteris_fardoulisi")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Nesonycteris fardoulisi"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Melonycteris_fardoulisi")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Nesonycteris fardoulisi"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Melonycteris_fardoulisi")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Nesonycteris fardoulisi"]), unit="km")

traits$MeanLat[which(traits$tip == "Melonycteris_woodfordi")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Nesonycteris woodfordi"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Melonycteris_woodfordi")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Nesonycteris woodfordi"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Melonycteris_woodfordi")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Nesonycteris woodfordi"]), unit="km")

traits$MeanLat[which(traits$tip == "Micropteropus_intermedius")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Epomophorus intermedius"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Micropteropus_intermedius")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Epomophorus intermedius"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Micropteropus_intermedius")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Epomophorus intermedius"]), unit="km")

traits$MeanLat[which(traits$tip == "Micropteropus_pusillus")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Epomophorus pusillus"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Micropteropus_pusillus")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Epomophorus pusillus"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Micropteropus_pusillus")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Epomophorus pusillus"]), unit="km")

traits$MeanLat[which(traits$tip == "Notopteris_macdonaldi")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Notopteris macdonaldii"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Notopteris_macdonaldi")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Notopteris macdonaldii"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Notopteris_macdonaldi")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Notopteris macdonaldii"]), unit="km")

traits$MeanLat[which(traits$tip == "Notopteris_neocaledonica")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Notopteris neocaledonicus"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Notopteris_neocaledonica")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Notopteris neocaledonicus"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Notopteris_neocaledonica")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Notopteris neocaledonicus"]), unit="km")

traits$MeanLat[which(traits$tip == "Pteropus_vetulus")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Pteropus vetula"]))[c(3,4)])
traits$MeanLong[which(traits$tip == "Pteropus_vetulus")] <- mean(ext(centroids(IUCN[IUCN$SCI_NAME == "Pteropus vetula"]))[c(1,2)])
traits$RangeArea[which(traits$tip == "Pteropus_vetulus")] <- expanse(aggregate(IUCN[IUCN$SCI_NAME == "Pteropus vetula"]), unit="km")

## add citations

## citation compilation, extracted february 22, 2024
## collect any citations per bat species
cites=c()
for(i in 1:length(traits$tip)) {
  
  counts=as.numeric(as.character(get_pubmed_ids(gsub('_','-',traits$tip[i]))$Count))
  cites[i]=counts
  print(paste(i,"/",nrow(traits)))
}

## virus related citation counts, could also make paramyxo specific
vcites=c()
for(i in 1:length(traits$tip)) {
  
  x=gsub('_','-',traits$tip[i])
  x=paste("(",x,")",sep="")
  x=paste(x,"AND (virus OR viral)")
  counts=as.numeric(as.character(get_pubmed_ids(x)$Count))
  vcites[i]=counts
  print(paste(i,"/",nrow(traits)))
}

## compile all citations
cdata=data.frame(tip=traits$tip,
                 cites=cites,
                 vcites=vcites) ## add vcites if desired

## merge
traits=merge(traits,cdata,by="tip")

## clean
rm(cites,i,counts,vcites,IUCN,combine,guy,set,wings,ABT,cdata,HMW,IUCNsp,name,x)

## add dummy variables for genus
dums=dummy_cols(traits["gen"])

## unique
dums=dums[!duplicated(dums$gen),]

## ensure all factor
for(i in 1:ncol(dums)){
  
  ## column as factor
  dums[,i]=factor(dums[,i])
  
}

traits=merge(traits,dums,by="gen",all.x=T)
rm(dums)

## simplify
temp=traits[c("tip","biogeographical_realm")]

## convert bioregion
x=temp$biogeographical_realm
x=strsplit(x,", ")

## get all unique regions
ureg=unique((melt(x)$value))
ureg=ureg[!is.na(ureg)]

## loop through and save
for(i in 1:length(ureg)){
  
  ## ifelse
  col=ifelse(x==ureg[i],1,0)
  temp[,2+i]=factor(col)
  
}

## assign names
names(temp)[3:(2+length(ureg))]=ureg
temp$biogeographical_realm=NULL

## merge
traits$biogeographical_realm=NULL
traits=merge(traits,temp,by="tip",all.x=T)
rm(temp,ureg,col,i,x)

## add evolutionary distance (upham et al. 2019)

## load Upham phylogeny
#setwd("/Users/danielbecker/Library/CloudStorage/OneDrive-UniversityofOklahoma/Becker Lab/Datasets/Upham phylo")
setwd("~/Documents/PhD/1. pteropodidae/phylo")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## fix tree tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep='_'))

## drop tips
tree=keep.tip(tree,traits$tip)
length(tree$tip.label)==nrow(traits)

## evolutionary distinctiveness
a=phyloregion::evol_distinct(tree,type="fair.proportion")
b=phyloregion::evol_distinct(tree,type="equal.splits")

## data
pdata=data.frame(ed_fp=a,
                 ed_es=b)
pdata$tip=rownames(pdata)
rm(a,b)

## merge with data
traits=merge(traits,pdata,by="tip")
rm(pdata, tree)

## save

setwd("~/Documents/PhD/1. pteropodidae")
#write.csv(traits, file="traits and response.csv", row.names = FALSE)
traits <- read.csv("traits and response.csv")

## mode function
mode.prop <- function(x) {
  ux <- unique(x[is.na(x)==FALSE])
  tab <- tabulate(match(na.omit(x), ux))
  max(tab)/length(x[is.na(x)==FALSE])
}

## assess variation across columns
vars=data.frame(apply(traits,2,function(x) mode.prop(x)),
                apply(traits,2,function(x) length(unique(x))))

vars$variables=rownames(vars)
names(vars)=c("var","uniq","column")
vars$var=round(vars$var,2)

## if homogenous (100%)
vars$keep=ifelse(vars$var<1,"keep","cut")
#vars$keep=ifelse(vars$column%in%c('hPCR','competence',"fam"),'keep',vars$keep)
vars=vars[order(vars$keep),]

## trim
keeps=vars[-which(vars$keep=="cut"),]$column

## drop if no variation
traits=traits[keeps]
rm(keeps,vars)

## assess missing values
mval=data.frame(apply(traits,2,function(x) length(x[!is.na(x)])/nrow(traits)))

## get names
mval$variables=rownames(mval)
names(mval)=c("comp","column")

#visualize NA distribution
ggplot(mval[!mval$column%in%c("gen","paramyxovirus.PCR.positive",
                              "paramyxovirus.seropositive",
                              "paramyxovirus.isolated",
                              "orthoparamyxovirus.PCR.positive",
                              "orthoparamyxovirus.seropositive",
                              "orthoparamyxovirus.isolated",
                              "rubulavirus.PCR.positive",
                              "rubulavirus.seropositive",
                              "rubulavirus.isolated","tip"),],
       aes(comp))+
  geom_histogram(bins=50)+
  geom_vline(xintercept=0.40,linetype=2,linewidth=0.5)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  labs(y="frequency",
       x="trait coverage across pteropodids")+
  scale_x_continuous(labels = scales::percent)

ggsave("Fig. S1.png", width=4, height=4)

## 25% cutoff for now
## try for a 40% here
mval$comp=round(mval$comp,2)
#mval$keep=ifelse(mval$comp>=0.25,"keep","cut")
mval$keep=ifelse(mval$comp>=0.40,"keep","cut")
table(mval$keep)
mval=mval[order(mval$keep),]
keeps=mval[-which(mval$keep=="cut"),]$column

## order
mval=mval[order(mval$comp),]

## drop if not well represented
traits=traits[keeps]
rm(mval,keeps)

## drop unnecessary columns
traits$clade=NULL
traits$fam=NULL
traits$tiplabel=NULL
traits$label=NULL

## transforming skewed variables

traits$upper_elevation_m <- log(traits$upper_elevation_m)
traits$lower_elevation_m <- log(traits$lower_elevation_m+1)
traits$altitude_breadth_m <- log(traits$altitude_breadth_m+1)
traits$Adult_mass_g <- log(traits$Adult_mass_g)
traits$Adult_forearm_mm <- log(traits$Adult_forearm_mm)
traits$Adult_length_mm <- log(traits$Adult_length_mm)
traits$RangeArea <- log(traits$RangeArea)
traits$cites <- log(traits$cites+1)
traits$vcites <- log(traits$vcites+1)

####residuals on body mass

##resid(lm(PMV_sero_2$Adult_mass_g ~ PMV_sero_2$Adult_length_mm, na.action = na.exclude))

setwd("~/Documents/PhD/1. pteropodidae")
write.csv(traits, file="cleaned traits and response.csv", row.names = FALSE)



#####################MAKING PHYLOGENY FIGURE

library(diversitree)

setwd("~/Documents/PhD/1. pteropodidae/phylo")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## fix tree tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep='_'))

## drop tips
tree=keep.tip(tree,data$tip)
length(tree$tip.label)==nrow(data)

treeorder <- data.frame(matrix(nrow=194))
treeorder$tip <- tree$tip.label
treeorder$id  <- 1:nrow(treeorder)
treeorder  <- merge(treeorder,data[,c(1,3,4,5)], by = "tip")
colnames(treeorder)[4:6] <- c("PCR", "serology", "isolation")
treeorder <- treeorder[order(treeorder$id), ]
tree$state <- treeorder[,c(4:6)]
rownames(tree$state) <- tree$tip.label

palette <- c('black','red')[tree$PCRpos]
plot(tree, type = 'fan', tip.color = palette, cex = 0.2)

png(file="~/Documents/PhD/pteropodidae/phylo/phylo_plot.png", 
    width=500, height=500, res=100, pointsize = 8)

jpeg(file="~/Documents/PhD/pteropodidae/phylo/phylo_plot.jpg", 
     width=500, height=500, res=100, pointsize = 8, quality=500)

pdf(file="~/Documents/PhD/pteropodidae/phylo/phylo_plot.pdf")
trait.plot(tree, tree$state, cols = list(serology = c("grey","blue"),
                                         PCR = c("grey","red"), 
                                         isolation = c("grey","green")),
           cex.lab = 0.5, margin=0.4, w=0.08, cex.legend = 1,
           lab = c("Serology","PCR","Isolation"))
dev.off()

tree$tip.label <- c(1:194)
rownames(tree$state) <- tree$tip.label
pdf(file="~/Documents/PhD/pteropodidae/phylo/phylo_plot_mini.pdf")
trait.plot(tree, tree$state, cols = list(serology = c("grey","blue"),
                                         PCR = c("grey","red"), 
                                         isolation = c("grey","green")),
           cex.lab = 0.01, margin=0, w=0.1, legend = FALSE, cex.legend = 1,
           lab = c("Serology","PCR","Isolation"))
dev.off()

######making IUCN ranges figure

library(terra)
library(letsR)

setwd("~/Documents/PhD/1. pteropodidae/databases/IUCN")
IUCN <- vect('data_0.shp')

setwd("~/Documents/PhD/1. pteropodidae")
predictedsp <- read.csv("all response predictions.csv")

maps <- letsR::lets.presab(IUCN)

# Define axis for colour breaks - this should be symmetric around zero and encompass the magnitude of species richness changes
brks.diff <- seq(-0,20,by=1)

# Define red-blue diverging colour palette
cols_diff_palette <- colorRampPalette(
  rev(c('#650A13','#b2182b','#d6604d','#f4a582','grey90','#92c5de','#4393c3','#2166ac','#0B2F52')))

# Now we define a vector of colours from our new palette that holds as many colours as we have break points above
cols.diff = cols_diff_palette(length(brks.diff)) 

plot(maps, xlim = c(-20,175))

# Assign range map of the Eurasian lynx to separate object
#top10 <- gsub("_"," ",predictedsp[order(-predictedsp$pred.x),]$tip[1:10])
known <- predictedsp[which(predictedsp$paramyxovirus.PCR.positive == 1),]
known <- gsub("_"," ",known$tip)

newsp <- predictedsp[which(predictedsp$paramyxovirus.PCR.positive == 0),]
top10novel <- gsub("_"," ",newsp[order(-newsp$pred.x),]$tip[1:10])
knownplus10 <- c(known,top10novel)

IUCNknown <- IUCN[IUCN$SCI_NAME %in% known,]
knownmap <- letsR::lets.presab(IUCNknown)

IUCNknownplus10 <- IUCN[IUCN$SCI_NAME %in% knownplus10,]
knownplus10map <- letsR::lets.presab(IUCNknownplus10)

#brks <- c(0,1,2,3,4,5,6,7,8,9,10)
cols1 <- colorRampPalette((c('#FBE9E7', '#FFCCBC', '#FFAB91', '#FF8A65', '#FF7043', '#FF5722')))
cols2 <- colorRampPalette((c('#FBE9E7', '#FFCCBC', '#FFAB91', '#FF8A65', '#FF7043', '#FF5722',
                            '#D84315', '#BF360C', '#B71c1c', '#880E4F')))

knownmap[["name"]] <- knownmap[["Species_name"]]

#par(mfrow=c(2,1))
plot(knownmap, xlim = c(-20,175), col_rich=cols1, breaks=brks, legend=NULL)
plot(knownplus10map, xlim = c(-20,175), col_rich=cols2, breaks=brks, legend=NULL)

plot(knownmap, xlim = c(110,160), ylim = c(-45,0), col_rich=cols2, breaks=brks, legend=NULL)


##########################


predictedsp$mean <- rowMeans(predictedsp[,c("predictedsp$pred.x", 
                                            "predictedsp$pred.y",
                                            "predictedsp$pred.x.x",
                                            "predictedsp$pred.y.y",
                                            "predictedsp$pred.x.x.x",
                                            "predictedsp$pred.y.y.y",
                                            "predictedsp$pred.x.x.x.x",
                                            "predictedsp$pred.y.y.y.y",
                                            "predictedsp$pred")], na.rm=TRUE)

predictedsp$mean <- rowMeans(predictedsp[,c(2,6,10,14,18,22,26,30,34)], na.rm=TRUE)

