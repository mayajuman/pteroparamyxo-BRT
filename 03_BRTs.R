## pteropodidae-paramyxovirus associations
## maya juman + dan becker
## updated 01/12/24
## updated/edited by DB 1/15/24

## clean up space
rm(list=ls()) 
graphics.off()

## set up

library(dplyr)
library(gbm)
library(fastDummies)
library(rsample)
library(ROCR)
library(sciplot)
library(ggplot2)
library(pdp)
library(PresenceAbsence)
library(tidyr)
library(viridis)
library(caper)
library(caret) 
library(InformationValue)
library(mgcv)
library(ggpubr)
library(plotrix)

setwd("~/Documents/PhD/1. pteropodidae")
#setwd("/Users/danielbecker/Library/CloudStorage/OneDrive-UniversityofOklahoma/Becker only/Students/Maya Juman")
data <- read.csv("cleaned traits and response.csv")

#tabulate PCR, sero, isolation for all paramyxos
table(data$paramyxovirus.PCR.positive)
table(data$paramyxovirus.seropositive)
table(data$paramyxovirus.isolated)

#tabulate PCR, sero, isolation for orthoparamyxovirinae
table(data$orthoparamyxovirus.PCR.positive)
table(data$orthoparamyxovirus.seropositive)
table(data$orthoparamyxovirus.isolated)

#tabulate PCR, sero, isolation for rubulavirinae
table(data$rubulavirus.PCR.positive)
table(data$rubulavirus.seropositive)
table(data$rubulavirus.isolated)

#column for ANY positives (n = 35)
data <- data %>% mutate(allpos = paramyxovirus.PCR.positive 
                        + paramyxovirus.seropositive
                        + paramyxovirus.isolated)
data$allpos[which(data$allpos == 2 | data$allpos == 3)] <- 1
table(data$allpos)

## simplify

set=data
set$gen=NULL
set$tip=NULL
set$StatusM=factor(set$StatusM)
set$TrendM=factor(set$TrendM)
set$activity_cycle=factor(set$activity_cycle)
set$island_dwelling=factor(set$island_dwelling)
set$dissected_by_mountains=factor(set$dissected_by_mountains)
set$hibernation_torpor=factor(set$hibernation_torpor)
set$PhyloClust41=factor(set$PhyloClust41)
set$sin_Citations=NULL

## covariates
ncol(set)-10

## coverage table s1
ts1=data.frame(apply(set,2,function(x) length(x[!is.na(x)])/nrow(set)))
set2 <- set[which(set$paramyxovirus.PCR.positive == 1),]
ts2=data.frame(apply(set2,2,function(x) length(x[!is.na(x)])/nrow(set2)))
ts1 <- cbind(ts1,ts2)

## get names
ts1$variables=rownames(ts1)
names(ts1)=c("comptot","comppos","column")
rownames(ts1)=NULL

## sort
ts1=ts1[order(ts1$column),]
ts1$comptot=round(ts1$comptot,2)

## trim disease data
ts1=ts1[!ts1$column%in%c("paramyxovirus.PCR.positive",
                         "paramyxovirus.seropositive",
                         "paramyxovirus.isolated",
                         "orthoparamyxovirus.PCR.positive",
                         "orthoparamyxovirus.seropositive",
                         "orthoparamyxovirus.isolated",
                         "rubulavirus.PCR.positive",
                         "rubulavirus.seropositive",
                         "rubulavirus.isolated",
                         "allpos"),]
ts1$feature=ts1$column
ts1$column=NULL
ts1$coverage_tot=ts1$comptot
ts1$comptot=NULL
ts1$coverage_pos=ts1$comppos
ts1$comppos=NULL

## create specific set for preliminary BRT: all-paramyxovirus PCR as response
PMV_PCR <- set %>% dplyr::select(paramyxovirus.PCR.positive, 
                              everything(), 
                              -paramyxovirus.seropositive, 
                              -paramyxovirus.isolated, 
                              -orthoparamyxovirus.PCR.positive, 
                              -orthoparamyxovirus.seropositive, 
                              -orthoparamyxovirus.isolated, 
                              -rubulavirus.PCR.positive, 
                              -rubulavirus.seropositive, 
                              -rubulavirus.isolated,
                              -allpos)

PMV_sero <- set %>% dplyr::select(paramyxovirus.seropositive, 
                              everything(), 
                              -paramyxovirus.PCR.positive, 
                              -paramyxovirus.isolated, 
                              -orthoparamyxovirus.PCR.positive, 
                              -orthoparamyxovirus.seropositive, 
                              -orthoparamyxovirus.isolated, 
                              -rubulavirus.PCR.positive, 
                              -rubulavirus.seropositive, 
                              -rubulavirus.isolated,
                              -allpos)

PMV_iso <- set %>% dplyr::select(paramyxovirus.isolated,
                              everything(), 
                              -paramyxovirus.PCR.positive, 
                              -paramyxovirus.seropositive, 
                              -orthoparamyxovirus.PCR.positive, 
                              -orthoparamyxovirus.seropositive, 
                              -orthoparamyxovirus.isolated, 
                              -rubulavirus.PCR.positive, 
                              -rubulavirus.seropositive, 
                              -rubulavirus.isolated,
                              -allpos)

OPMV_PCR <- set %>% dplyr::select(orthoparamyxovirus.PCR.positive,
                              everything(), 
                              -paramyxovirus.PCR.positive, 
                              -paramyxovirus.seropositive, 
                              -paramyxovirus.isolated, 
                              -orthoparamyxovirus.seropositive, 
                              -orthoparamyxovirus.isolated, 
                              -rubulavirus.PCR.positive, 
                              -rubulavirus.seropositive, 
                              -rubulavirus.isolated,
                              -allpos)

OPMV_sero <- set %>% dplyr::select(orthoparamyxovirus.seropositive,
                              everything(), 
                              -paramyxovirus.PCR.positive, 
                              -paramyxovirus.seropositive, 
                              -paramyxovirus.isolated, 
                              -orthoparamyxovirus.PCR.positive, 
                              -orthoparamyxovirus.isolated, 
                              -rubulavirus.PCR.positive, 
                              -rubulavirus.seropositive, 
                              -rubulavirus.isolated,
                              -allpos)

OPMV_iso <- set %>% dplyr::select(orthoparamyxovirus.isolated,
                              everything(), 
                              -paramyxovirus.PCR.positive, 
                              -paramyxovirus.seropositive, 
                              -paramyxovirus.isolated, 
                              -orthoparamyxovirus.PCR.positive, 
                              -orthoparamyxovirus.seropositive, 
                              -rubulavirus.PCR.positive, 
                              -rubulavirus.seropositive, 
                              -rubulavirus.isolated,
                              -allpos)

RPMV_PCR <- set %>% dplyr::select(rubulavirus.PCR.positive,
                              everything(), 
                              -paramyxovirus.PCR.positive, 
                              -paramyxovirus.seropositive, 
                              -paramyxovirus.isolated, 
                              -orthoparamyxovirus.PCR.positive, 
                              -orthoparamyxovirus.seropositive, 
                              -orthoparamyxovirus.isolated, 
                              -rubulavirus.seropositive, 
                              -rubulavirus.isolated,
                              -allpos)

RPMV_sero <- set %>% dplyr::select(rubulavirus.seropositive,
                              everything(), 
                              -paramyxovirus.PCR.positive, 
                              -paramyxovirus.seropositive, 
                              -paramyxovirus.isolated, 
                              -orthoparamyxovirus.PCR.positive, 
                              -orthoparamyxovirus.seropositive, 
                              -orthoparamyxovirus.isolated, 
                              -rubulavirus.PCR.positive, 
                              -rubulavirus.isolated,
                              -allpos)

RPMV_iso <- set %>% dplyr::select(rubulavirus.isolated,
                               everything(), 
                               -paramyxovirus.PCR.positive, 
                               -paramyxovirus.seropositive, 
                               -paramyxovirus.isolated, 
                               -orthoparamyxovirus.PCR.positive, 
                               -orthoparamyxovirus.seropositive, 
                               -orthoparamyxovirus.isolated, 
                               -rubulavirus.PCR.positive, 
                               -rubulavirus.seropositive,
                               -allpos)

allpos <- set %>%
  dplyr::select(allpos, everything(),
                -paramyxovirus.PCR.positive,
                -paramyxovirus.seropositive, 
                -paramyxovirus.isolated, 
                -orthoparamyxovirus.PCR.positive, 
                -orthoparamyxovirus.seropositive, 
                -orthoparamyxovirus.isolated, 
                -rubulavirus.PCR.positive, 
                -rubulavirus.seropositive, 
                -rubulavirus.isolated)

## hyperparameter tuning ifelse

## hyperparameter grid
hgrid=expand.grid(n.trees=3000,
                  interaction.depth=c(2,3,4),
                  shrinkage=c(0.01,0.001,0.0005),
                  prop=c(0.7,0.8,0.9), ## added
                  n.minobsinnode=4,
                  seed=seq(1,5,by=1))

## fix trees
hgrid$n.trees=ifelse(hgrid$shrinkage<0.001,hgrid$n.trees*3,hgrid$n.trees)

## trees, depth, shrink, min, prop
hgrid$id=with(hgrid,paste(n.trees,interaction.depth,shrinkage,n.minobsinnode))

## sort by id then seed
hgrid=hgrid[order(hgrid$id,hgrid$seed),]

## now add rows
hgrid$row=1:nrow(hgrid)

## factor id
hgrid$id2=factor(as.numeric(factor(hgrid$id)))

## function to assess each hyperpar combination
hfit=function(row,resp,data){
  
  ## make new data
  ndata=data
  
  ## correct response
  ndata$response=ndata[resp][,1]
  
  ## remove raw
  ndata[resp][,1]=NULL
  
  ## use rsample to split
  set.seed(hgrid$seed[row])
  split=initial_split(ndata,prop=hgrid$prop[row],strata="response")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response
  
  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=hgrid$n.trees[row],
             distribution="bernoulli",
             shrinkage=hgrid$shrinkage[row],
             interaction.depth=hgrid$interaction.depth[row],
             n.minobsinnode=hgrid$n.minobsinnode[row],
             cv.folds=5,class.stratify.cv=TRUE,
             bag.fraction=0.5,train.fraction=1,
             n.cores=1,
             verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter=gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## known
  result=dataTest$response
  
  ## sensitiviy and specificity
  sen=InformationValue::sensitivity(result,preds)
  spec=InformationValue::specificity(result,preds)
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## print
  print(paste("hpar row ",row," done; test AUC is ",auc_test,sep=""))
  
  ## save outputs
  return(list(best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              wrow=row))
}

## run the function
hpars=lapply(1:nrow(hgrid),function(x) hfit(x,resp="paramyxovirus.PCR.positive",data=PMV_PCR))

## get results
hresults=data.frame(sapply(hpars,function(x) x$trainAUC),
                    sapply(hpars,function(x) x$testAUC),
                    sapply(hpars,function(x) x$spec),
                    sapply(hpars,function(x) x$sen),
                    sapply(hpars,function(x) x$wrow),
                    sapply(hpars,function(x) x$best))
names(hresults)=c("trainAUC","testAUC",
                  "spec","sen","row","best")

## combine and save
search=merge(hresults,hgrid,by="row")

## export
setwd("~/Documents/PhD/1. pteropodidae")
write.csv(search,"par tuning data summary 2.csv")

search <- read.csv("par tuning data summary PMV PCR.csv")

## factor parameters
search$shrinkage=factor(search$shrinkage)
lvl=rev(sort(unique(search$shrinkage)))
search$shrinkage=factor(search$shrinkage,levels=lvl); rm(lvl)

## factor other
search$interaction.depth=factor(search$interaction.depth)
search$prop=factor(search$prop)

## PCR beta regression for AUC
mod=gam(testAUC~interaction.depth*shrinkage*prop,
        data=search,method="REML",family=betar)
anova(mod)

## PCR beta regression for sensitivity
mod=gam(sen~interaction.depth*shrinkage*prop,
        data=search,method="REML",family=betar)
anova(mod)

## PCR beta regression for specificity
mod=gam(spec~interaction.depth*shrinkage*prop,
        data=search,method="REML",family=betar)
anova(mod)

## recast from wide to long
search2=gather(search,measure,value,testAUC:sen)

## revalue and factor
search2$measure=plyr::revalue(search2$measure,
                              c("sen"="sensitivity",
                                "spec"="specificity",
                                "testAUC"="test AUC"))
search2$measure=factor(search2$measure,
                       levels=c("test AUC","sensitivity","specificity"))

## visualize
png("Figure S2.png",width=5,height=8,units="in",res=600)
set.seed(1)
ggplot(search2,aes(shrinkage,value,
                   colour=interaction.depth,fill=interaction.depth))+
  geom_boxplot(alpha=0.25)+
  geom_point(alpha=0.75,
             position = position_jitterdodge(dodge.width=0.75))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  #facet_wrap(~measure)+
  facet_grid(measure~prop,scales="free_y",switch="y")+
  theme(strip.placement="outside",
        strip.background=element_blank())+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  scale_color_brewer(palette="Pastel2")+
  scale_fill_brewer(palette="Pastel2")+
  guides(colour=guide_legend(title="interaction depth"),
         fill=guide_legend(title="interaction depth"))+
  labs(y=NULL,
       x="learning rate")+
  scale_y_continuous(n.breaks=4)
dev.off()

## clean
rm(search,search2,mod)

## brt function to use different data partitions
brt_part=function(seed,response){
  
  ## make new data
  ndata=PMV_PCR
  
  ## correct response
  ndata$response=ndata[response][,1]
  
  ## remove raw
  ndata$paramyxovirus.PCR.positive=NULL
  
  ## fix cites if response
  if(response=="cites"){
    
    ## plus 1 for 0
    ndata$cites=ifelse(ndata$cites==0,1,ndata$cites)
    
  }else{
    
    ndata=ndata
    
  }
  
  ## use rsample to split, 80/20
  set.seed(seed)
  split=initial_split(ndata,prop=0.8,strata="response")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response
  
  ## dist
  dist=ifelse(response=="cites","poisson","bernoulli")
  
  ## n.trees
  nt=ifelse(response=="cites",5000,3000)
  
  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=nt,
             distribution=dist,
             shrinkage=0.001,
             interaction.depth=3,
             n.minobsinnode=4,
             cv.folds=5,class.stratify.cv=TRUE,
             bag.fraction=0.5,train.fraction=1,
             n.cores=1,
             verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter=gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## known
  result=dataTest$response
  
  ## sensitiviy and specificity
  sen=InformationValue::sensitivity(result,preds)
  spec=InformationValue::specificity(result,preds)
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## skip if poisson
  if(response=="cites"){
    
    perf=NA
    
  }else{
    
    ## inner loop if yTest is all 0
    if(var(yTest)==0){
      
      perf=NA
    }else{
      
      ## ROC
      pr=prediction(preds,dataTest$response)
      perf=performance(pr,measure="tpr",x.measure="fpr")
      perf=data.frame(perf@x.values,perf@y.values)
      names(perf)=c("fpr","tpr")
      
      ## add seed
      perf$seed=seed
      
    }
  }
  
  ## relative importance
  bars=summary(gbmOut,n.trees=best.iter,plotit=F)
  bars$rel.inf=round(bars$rel.inf,2)
  
  ## predict with cites as is
  preds=predict(gbmOut,data,n.trees=best.iter,type="response")
  pred_data=data[c("tip","gen","paramyxovirus.PCR.positive")]
  pred_data$pred=preds
  pred_data$type=response
  
  ## predict with mean cites
  pdata=data
  pdata$cites=mean(pdata$cites)
  pred_data$cpred=predict(gbmOut,pdata,n.trees=best.iter,type="response")
  
  ## predict with mean viral cites
  vdata=data
  vdata$vcites=mean(vdata$vcites)
  pred_data$vcpred=predict(gbmOut,vdata,n.trees=best.iter,type="response")
  
  ## predict with mean cites and viral cites
  vcdata=data
  vcdata$cites=mean(vcdata$cites)
  vcdata$vcites=mean(vcdata$vcites)
  pred_data$vcpred2=predict(gbmOut,vcdata,n.trees=best.iter,type="response")
  
  ## sort
  pred_data=pred_data[order(pred_data$pred,decreasing=T),]
  
  ## print
  print(paste("BRT ",seed," done; test AUC = ",auc_test,sep=""))
  
  ## save outputs
  return(list(mod=gbmOut,
              best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              roc=perf,
              rinf=bars,
              predict=pred_data,
              traindata=dataTrain,
              testdata=dataTest,
              seed=seed))
}

## apply across X splits each (100 better)
smax=25
total_brts=lapply(1:smax,function(x) brt_part(seed=x,response="paramyxovirus.PCR.positive"))

## write to files
#setwd("~/Documents/PhD/pteropodidae")
#saveRDS(pcr_brts,"pcr brts.rds")

## run wos brts
#pm_brts=lapply(1:(smax-1),function(x) brt_part(seed=x,response="cites"))

## write
#saveRDS(pm_brts,"pm brts.rds")

## mean AUC
round(mean(sapply(total_brts,function(x) x$testAUC),na.rm=TRUE),4) #test
round(se(sapply(total_brts,function(x) x$testAUC))*100,2) ## 92% accuracy
round(mean(sapply(total_brts,function(x) x$trainAUC),na.rm=TRUE),3) #train


round(mean(sapply(total_brts,function(x) x$best),na.rm=TRUE),2) #best.iter

## relative importance
vinf=lapply(total_brts,function(x) x$rinf)
vinf=do.call(rbind,vinf)

## aggregate mean
vdata=data.frame(aggregate(rel.inf~var,data=vinf,mean),
                 aggregate(rel.inf~var,data=vinf,se)["rel.inf"])
names(vdata)=c("var","rel.inf","rse")
vdata=vdata[order(vdata$rel.inf,decreasing=T),]
rm(vinf)

## make rmin and rmax
vdata$rmin=vdata$rel.inf-vdata$rse
vdata$rmax=vdata$rel.inf+vdata$rse

## visualize
vset=vdata
vset=vset[vset$rel.inf>1,]
vset=vset[order(vset$rel.inf),]
vset$var=factor(vset$var,levels=unique(vset$var))
ggplot(vset,aes(reorder(var,rel.inf,max),rel.inf))+
  geom_segment(aes(y=0,yend=rel.inf,
                   x=var,xend=var))+
  geom_point(size=1.5)+
  #scale_y_sqrt()+
  coord_flip()+
  theme_bw()+
  labs(x=NULL,
       y="mean feature importance")+
  theme(axis.title=element_text(size=10),
        axis.text=element_text(size=8))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

## average predictions: PCR
total_apreds=lapply(total_brts,function(x) x$predict)
total_apreds=do.call(rbind,total_apreds)

## aggregate
total_apreds=data.frame(aggregate(pred~tip,data=total_apreds,mean),
                      aggregate(cpred~tip,data=total_apreds,mean)['cpred'], ## holding wos constant
                      aggregate(vcpred~tip,data=total_apreds,mean)['vcpred'],
                      aggregate(vcpred2~tip,data=total_apreds,mean)['vcpred2'],
                      aggregate(paramyxovirus.PCR.positive~tip,data=total_apreds,prod)
                      ["paramyxovirus.PCR.positive"])

## type
PMV_pcr_apreds$type='PCR'

#write.csv(vars, file="vars.csv", row.names = FALSE)
#write.csv(pcr_apreds, file="predictions.csv", row.names = FALSE)

df_list <- list(PMV_pcr_apreds, PMV_sero_apreds, PMV_iso_apreds,
                OPMV_pcr_apreds, OPMV_sero_apreds, OPMV_iso_apreds,
                RPMV_pcr_apreds, RPMV_sero_apreds, RPMV_iso_apreds)      

#merge all data frames together
all<- df_list %>% reduce(full_join, by='tip')
write.csv(all, file="all response predictions.csv", row.names = FALSE)

vdata1 <- vdata_PMVPCR %>% dplyr::select(var,rel.inf,rse) %>%
  dplyr::rename(PMVpcr_rel.inf = rel.inf, OMVpcr_rse = rse)

vdata2 <- vdata_PMVsero %>% dplyr::select(var,rel.inf,rse) %>%
  dplyr::rename(PMVsero_rel.inf = rel.inf, PMVsero_rse = rse)

vdata3 <- vdata_PMViso %>% dplyr::select(var,rel.inf,rse) %>%
  dplyr::rename(PMViso_rel.inf = rel.inf, PMViso_rse = rse)

vdata4 <- vdata_OPMVPCR %>% dplyr::select(var,rel.inf,rse) %>%
  dplyr::rename(OPMVpcr_rel.inf = rel.inf, OPMVpcr_rse = rse)

vdata5 <- vdata_OPMVsero %>% dplyr::select(var,rel.inf,rse) %>%
  dplyr::rename(OPMVsero_rel.inf = rel.inf, OPMVsero_rse = rse)

vdata6 <- vdata_OPMViso %>% dplyr::select(var,rel.inf,rse) %>%
  dplyr::rename(OPMViso_rel.inf = rel.inf, OPMViso_rse = rse)

vdata7 <- vdata_RPMVPCR %>% dplyr::select(var,rel.inf,rse) %>%
  dplyr::rename(RPMVpcr_rel.inf = rel.inf, RPMVpcr_rse = rse)

vdata8 <- vdata_RPMVsero %>% dplyr::select(var,rel.inf,rse) %>%
  dplyr::rename(RPMVsero_rel.inf = rel.inf, RPMVsero_rse = rse)

vdata9 <- vdata_RPMViso %>% dplyr::select(var,rel.inf,rse) %>%
  dplyr::rename(RPMViso_rel.inf = rel.inf, RPMViso_rse = rse)

vlist <- list(vdata1, vdata2, vdata3, vdata4, vdata5, vdata6, vdata7,
              vdata8, vdata9)      

#merge all data frames together
allv <- vlist %>% reduce(full_join, by='var')
write.csv(allv, file="all response variables.csv", row.names = FALSE)


#####target shuffling: modified from barbara's code
library(caTools)

permutedAUC<-c()

i=1
while (i <= 10) {
  # for permutation loop
  
  ## random permutation of Label
  randomLabel<-sample(data$paramyxovirus.PCR.positive)
  data2<-cbind(randomLabel,data)
  data2 <- data2[,-c(2:12,85)]
  #data2 <- data2 %>% dplyr::rename(response = paramyxovirus.PCR.positive)
  #data2[,1]<-sapply(data2[,1],as.character)
  data2$PhyloClust41 <- as.factor(data2$PhyloClust41)
  data2$StatusM <- as.factor(data2$StatusM)
  data2$TrendM <- as.factor(data2$TrendM)

  ## create training and test sets
  intrain2<-createDataPartition(y=data2$randomLabel,
                                times=1,
                                p=0.8,
                                list=FALSE)
  
  test2<-data2[-intrain2,]
  training2<-data2[intrain2,]
  
  check<-1-is.na(training2)*1
  checksum<-apply(check,2,sum)
  if(length(which(checksum>=2))==84){
    
    
    ## random permutation of Labels ~ families + traits
    model2<-as.formula(paste(colnames(data2)[1], "~",
                             paste(colnames(data2)[c(2:84)], collapse="+"), #traits
                             sep = ""))
    
    batgbm2<- gbm(model2,
                  data=training2, 
                  distribution="bernoulli",
                  n.trees=3000,
                  shrinkage=0.001,
                  interaction.depth=3,
                  bag.fraction=0.50,
                  train.fraction=1,
                  n.minobsinnode=4,
                  cv.folds=5,
                  keep.data=TRUE,
                  class.stratify.cv=TRUE)
    
    #check performance using 5-fold cross-validation
    best.iter2 <- gbm.perf(batgbm2,method="cv",plot.it=FALSE) #OOB method under predicts
    #   batsum2<-summary.gbm(batgbm2,n.trees=best.iter,method=relative.influence,plotit=FALSE)
    
    ## LABEL
    ## predictions on the TRAINING SET
    output2<-predict(batgbm2, newdata=training2, n.trees=best.iter2, type="response") 
    output2<-cbind(output2,as.numeric(as.factor(training2$randomLabel)))
    #   colnames(output2)<-c("output","label")
    #   output2<-output2[order(-as.numeric(output2[,1])),]
    
    # # training AUC for Bernoulli distributed responses
    auc2=colAUC(output2[,1],output2[,2])
    
    # Predictions on the TEST set
    output.test2<-predict(batgbm2, newdata=test2, n.trees=best.iter2, type="response") 
    output.test2<-cbind(output.test2,as.numeric(as.factor(test2$randomLabel)))
    # colnames(output.test2)<-c("output","label")
    # output.test2<-output.test2[order(-output.test2[,1]),]
    # plot(output.test)
    
    ## test AUC for Bernoulli distributed responses
    auctest2=colAUC(output.test2[,1],output.test2[,2])
    
    permutedAUC[i]<-auctest2
    print(auctest2)
    i=i+1
  } else i=i
}

sum(is.na(permutedAUC)*1) #how many NAs
permutedAUC2<-na.omit(permutedAUC)
mean(permutedAUC2)
se(permutedAUC2)



#######partial dependence plots

pdp_agg=function(mod,feature){
  
  ## just the plot function
  pdep=plot(mod$mod,feature,
            return.grid=T,
            n.trees=mod$best,
            plot=F,
            continuous.resolution=200,
            type="response")
  
  ## add seed
  pdep$seed=unique(mod$roc$seed)
  
  ## save predictor
  pdep$predictor=pdep[feature][,1]
  
  ## order
  pdep=pdep[order(pdep$predictor),]
  
  ## get rank
  pdep$rank=1:nrow(pdep)
  
  ## save yhat
  pdep$yhat=pdep$y
  
  ## return
  return(pdep)
  
}


pdp_plot=function(bmods,feature,xlab){
  
  ## pdp_agg
  agg=do.call(rbind,lapply(bmods,function(x) pdp_agg(x,feature)))
  
  ## get class of the feature
  cl=class(data[feature][,1])
  
  ## if else based on type
  if(cl%in%c("numeric","integer")){
    
    ## get element-wise means
    x=with(agg,tapply(predictor,rank,mean))
    y=with(agg,tapply(yhat,rank,mean))
    
    ## save as mean
    pmean=data.frame(predictor=x,yhat=y)
    
    ## get yrange
    yrange=range(agg$yhat,pmean$yhat,na.rm=T)
    
    ## get histogram
    hi=hist(data[feature][,1],breaks=30,plot=F)
    hi=with(hi,data.frame(breaks[1:(length(breaks)-1)],counts))
    names(hi)=c("mids","counts")
    
    ## ggplot it
    ggplot(agg,aes(predictor,yhat,group=seed))+
      
      ## add histogram
      geom_segment(data=hi,inherit.aes=F,
                   aes(x=mids,xend=mids,
                       y=yrange[1],yend=plotrix::rescale(counts,yrange)),
                   size=1,colour="grey",alpha=0.25)+
      
      ## add lines
      geom_line(linewidth=1,alpha=0.25,colour="grey")+
      
      ## add mean
      geom_line(data=pmean,linewidth=1,inherit.aes=F,
                aes(predictor,yhat))+
      
      ## theme
      theme_bw()+
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=8))+
      theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0)))+
      theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0)))+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      labs(x=xlab,y="marginal effect")+
      scale_y_continuous(labels=scales::number_format(accuracy=0.01))
    
    ## end numeric
  }else{ ## factor-based plot
    
    ## get element-wise means
    y=with(agg,tapply(yhat,predictor,mean))
    
    ## save as mean
    #pmean=data.frame(predictor=x,yhat=y)
    pmean=data.frame(y)
    names(pmean)="yhat"
    pmean$predictor=rownames(pmean)
    rownames(pmean)=NULL
    
    ## make temp data
    temp=data
    temp$predictor=temp[feature][,1]
    
    ## do nothing
    agg=agg
    pmean=pmean
    temp=temp
    
    ## get yrange
    yrange=range(agg$yhat,pmean$yhat,na.rm=T)
    
    ## fix temp to yrange
    temp$yhat=ifelse(temp$comp==1,max(yrange),min(yrange))
    
    ## ggplot with rug
    set.seed(1)
    ggplot(agg,aes(predictor,yhat,group=seed))+
      
      ## add individual BRTs
      geom_jitter(size=1,alpha=0.25,colour="grey",width=0.1)+
      
      ## add mean
      geom_point(data=pmean,size=2,inherit.aes=F,shape=15,
                 aes(predictor,yhat))+
      
      ## add rug
      geom_rug(data=temp,inherit.aes=F,
               aes(predictor,yhat),
               sides="b",position="jitter",
               colour="grey",alpha=0.25,
               na.rm=T)+
      
      ## theme
      theme_bw()+
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7))+
      theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0)))+
      theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0)))+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      labs(x=xlab,y="marginal effect")+
      scale_y_continuous(limits=c(yrange[1]-0.01,yrange[2]+0.01),
                         labels=scales::number_format(accuracy=0.01))
    
  }
  
}

#partial dependence plots
vset=vset[rev(order(vset$rel.inf)),]

#vcites

pdp1 <- pdp_plot(total_brts,as.character(vset$var[1]),"Virus-related citations")
pdp2 <- pdp_plot(total_brts,as.character(vset$var[2]),"Range Area")
pdp3 <- pdp_plot(total_brts,as.character(vset$var[3]),"Citations")
pdp4 <- pdp_plot(total_brts,as.character(vset$var[4]),"Adult length (mm)")
pdp5 <- pdp_plot(total_brts,as.character(vset$var[5]),"Ear length (mm)")
pdp6 <- pdp_plot(total_brts,as.character(vset$var[6]), "Evolutionary distance (es)")
#pdp7 <- pdp_plot(total_brts,as.character(vset$var[7]))
#pdp8 <- pdp_plot(total_brts,as.character(vset$var[8]))

pdps <- ggarrange(pdp1,pdp2,pdp3,pdp4,pdp5,pdp6,
          ncol=3,nrow=2,#widths=c(4,4),heights=c(22,1),
          labels=c("A","B","C","D","E","F"),
          #label.x=c(0,-0.1), label.y=0.001,
          font.label=list(face="plain",size=12))

ggsave("partial dependence plots.png", width = 6, height = 4)


####threholding predictions

setwd("~/Documents/PhD/1. pteropodidae")
data <- read.csv("Table S3.csv")
data <- data[,-1]

library(PresenceAbsence)
optimal.thresholds(DATA = data, which.model = 1, 
                   model.names = c("pred","vcpred2"),
                   opt.methods = 10,
                   req.sens = 0.90)
optimal.thresholds(DATA = data, which.model = 2, 
                   model.names = c("pred","vcpred2"),
                   opt.methods = 10,
                   req.sens = 0.90)
