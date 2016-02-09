###################################
## Replication code for models
## reported in "Legislative Effects 
## of Electoral Mandates" (with M.
## Tavits; forthcoming in BJPS).
## 
## Author: Santiago Olivella
##        (olivella@wustl.edu)
##
## Last update: SO 11/18/12
####################################


library(foreign)
library(car)
library(rjags)
library(MCMCpack)

rm(list=ls())
set.seed(831213)
setwd("~/Dropbox/ElectionsAndDiscipline/Data and Analysis/")



tiersdb  <- read.dta("../BJPSFinal/ReplicationData.dta")
tiersdb$nomchange <- ifelse(tiersdb$nom_lumped==tiersdb$lag_nom_lumped,0,1)
tiersdb$mandatechangeLumped <- ifelse(tiersdb$smd_pr==1|tiersdb$pr_smd==1,1
                                      ,0)
tiersdb$mandatechange1 <- ifelse(tiersdb$smd_rpr==1
                                 |tiersdb$rpr_smd==1
                                 |tiersdb$npr_smd==1,1
                                 ,0)
tiersdb$mandatechange2 <- ifelse(tiersdb$smd_npr==1
                                 |tiersdb$rpr_npr==1
                                 |tiersdb$npr_rpr==1,1
                                 ,0)
tiersdb$fromSMD <- ifelse(tiersdb$smd_smd==1|tiersdb$smd_pr==1,1,0)
tiersdb$fromPR <- ifelse(tiersdb$npr_npr==1
                         |tiersdb$npr_smd==1
                         |tiersdb$npr_rpr==1
                         |tiersdb$rpr_rpr==1
                         |tiersdb$rpr_smd==1  
                         |tiersdb$rpr_npr==1   
                         |tiersdb$pr_smd,1,0)
tiersdb$fromNPR <- ifelse(tiersdb$npr_npr==1
                          |tiersdb$npr_smd==1
                          |tiersdb$npr_rpr==1
                          ,1,0)
tiersdb$fromRPR <- ifelse(tiersdb$rpr_rpr==1
                          |tiersdb$rpr_smd==1  
                          |tiersdb$rpr_npr==1   
                          ,1,0)


tiersdb$originAgg <- ifelse(tiersdb$fromSMD==1,1
                            ,ifelse(tiersdb$fromPR==1,2
                                    ,NA)) 
tiersdb$originDisagg <- ifelse(tiersdb$fromSMD==1,1
                               ,ifelse(tiersdb$fromRPR==1,2
                                       ,ifelse(tiersdb$fromNPR==1,3
                                               ,NA))) 
DicSafeSMD <- ifelse(tiersdb$mar_50>0,1,ifelse(tiersdb$mar_50<0,0,NA))
DicSafePR <- ifelse(tiersdb$intseatpr>1|tiersdb$intseatnat>1,1,ifelse(is.na(tiersdb$intseatnat)|is.na(tiersdb$intseatpr),NA,0))
tiersdb$SafeIndicator <- ifelse(DicSafeSMD==1&DicSafePR==1,"sSMDsPR",
                                ifelse(DicSafeSMD==1&DicSafePR==0,"sSMDnPR",
                                       ifelse(DicSafeSMD==0&DicSafePR==1,"nSMDsPR",
                                              ifelse(DicSafeSMD==0&DicSafePR==0,"nSMDnPR","OneNom"))))
tiersdb$SafeIndicator <- recode(tiersdb$SafeIndicator, "NA='OneNom'")
tiersdb$SafeIndicator <- as.factor(tiersdb$SafeIndicator)

temp <- tiersdb[!duplicated(subset(tiersdb,select=c(year,size))),]
sizeGov <- c(by(temp$size,list(temp$year,temp$government),sum))
sizeGovDiff <- sizeGov[6:10]-sizeGov[1:5]
tiersdb$sizeGovDiff <- ifelse(tiersdb$year==1994,sizeGovDiff[1],
                              ifelse(tiersdb$year==1998,sizeGovDiff[2],
                                     ifelse(tiersdb$year==2002,sizeGovDiff[3],
                                            ifelse(tiersdb$year==2006,sizeGovDiff[4],sizeGovDiff[5]))))


model.data <- subset(tiersdb
                     ,excl_switch==0 & not_consecutive==0
                     ,select=c(percchangew
                               ,wscorediff
                               ,mandatechangeLumped
                               ,mandatechange1
                               ,mandatechange2
                               ,nomchange
                               ,originAgg
                               ,originDisagg
                               ,year
                               ,finalparty
                               ,leadership
                               ,distance
                               ,security
                               ,doublenom
                               ,SafeIndicator
                               ,sizeGovDiff
                               ))

full.model.data <- model.data
write.csv(full.model.data,file="ModelData.csv",row.names=FALSE)

niter <- 10e4

## Aggregated
model.data.agg <- subset(full.model.data
                         ,select=c(wscorediff
                                   ,percchangew
                                   ,mandatechangeLumped
                                   ,originAgg
                                   ,year
                                   ,finalparty
                                   #,leadership
                                   ,distance
                                   ,security
                                   ,doublenom
                                   #,SafeIndicator
                                   #,sizeGovDiff
                                   ))
model.data.agg <- model.data[complete.cases(model.data.agg),]
model.data.agg$year <- as.factor(model.data.agg[,"year"])
model.data.agg$finalparty <- as.factor(model.data.agg[,"finalparty"])
y <- model.data.agg$percchangew
origin <- model.data.agg$originAgg
man.change <- model.data.agg$mandatechangeLumped
X <- model.matrix(y~.-1,subset(model.data.agg, select=c(
                                                        finalparty  
                                                        ,year
                                                        ,doublenom
                                                        #,sizeGovDiff
                                                        #,SafeIndicator
                                                        ,security
                                                        #,leadership
                                                        ,distance
                                                        )))[,-1]#[,-c(2,3,5,8,11)]
J=length(unique(origin))
N=length(y)
K=dim(X)[2]
W=diag(2)
model.data.list.agg <- list(y=y
                            ,man.change=man.change
                            ,origin=origin
                            ,X=X
                            ,N=N
                            ,J=J
                            ,K=K
                            ,W=W
                            )

model.inits.agg <- function(){
                              list(
                                B.raw=array(rnorm(2*J),c(J,2))
                                ,mu.a.raw=rnorm(1)
                                ,mu.t.raw=rnorm(1)
                                ,sigma.y=runif(1)
                                ,Tau.B.raw=rwish(3,diag(2))
                                ,xi.a=runif(1)
                                ,xi.t=runif(1)
                                ,betavec=rnorm(K)
                                )
}

el.model.agg <- jags.model(file="hierModel.jag"
                           ,data=model.data.list.agg
                           ,inits=model.inits.agg
                           ,n.chains=2
                           ,n.adapt=niter)

post.samples.agg <- coda.samples(
                                el.model.agg
                                ,n.iter=niter
                                ,thin=15
                                ,variable.names=c("betavec"
                                                  ,"e.y"
                                                  ,"alpha"
                                                  ,"theta"
                                                  ,"sigma.a"
                                                  ,"sigma.t"
                                                  ,"sigma.y")
                                )
summary.agg <- summary(post.samples.agg,quantiles=c(0.05,0.5,0.95))
residuals.agg <- post.samples.agg[[1]][,grep("e.y",colnames(post.samples.agg[[1]]))]
R2.agg <- 1-((mean(apply(residuals.agg,1,var)))/(var(y)))
summary.agg$statistics[-grep("e.y",rownames(summary.agg$statistics)),1:2]

##Vanilla aggregated model
vanilla.agg <- lmer(y ~ X + man.change + (1+man.change|origin)
                    ,control=list(maxIter=1000)
                    #,REML=FALSE
                    )
summary(vanilla.agg)

## Disaggregated

model.data.dis <- subset(full.model.data
                         ,select=c(percchangew
                                   ,mandatechange1
                                   ,mandatechange2
                                   ,originDisagg
                                   ,year
                                   ,finalparty
                                   #,leadership
                                   ,distance
                                   ,security
                                   ,doublenom))
model.data.dis <- model.data[complete.cases(model.data.dis),]
y <- model.data.dis$percchangew
model.data.dis$year <- as.factor(model.data.dis[,"year"])
model.data.dis$finalparty <- as.factor(model.data.dis[,"finalparty"])
origin <- model.data.dis$originDisagg
man.change1 <- model.data.dis$mandatechange1
man.change2 <- model.data.dis$mandatechange2
X <- model.matrix(y~.-1,subset(model.data.dis, select=c(
                                                        finalparty  
                                                        ,year
                                                        ,doublenom
                                                        ,security
                                                        #,leadership
                                                        ,distance
                                                        )))[,-1]
J=length(unique(origin))
N=length(y)
K=dim(X)[2]
W=diag(3)
model.data.list.dis <- list(y=y
                            ,man.change1=man.change1
                            ,man.change2=man.change2
                            ,origin=origin
                            ,X=X
                            ,N=N
                            ,J=J
                            ,K=K
                            ,W=W
                            )
model.inits.dis <- function(){
                              list(
                                B.raw=array(rnorm(J*3),c(J,3))
                                ,mu.a.raw=rnorm(1)
                                ,mu.t1.raw=rnorm(1)
                                ,mu.t2.raw=rnorm(1)
                                ,sigma.y=runif(1)
                                ,Tau.B.raw=rwish(4,diag(3))
                                ,xi.a=runif(1)
                                ,xi.t1=runif(1)
                                ,xi.t2=runif(1)
                                ,betavec=rnorm(K)
                                )
}

el.model.dis <- jags.model(file="hierModelDisagg.jag"
                           ,data=model.data.list.dis
                           ,inits=model.inits.dis
                           ,n.chains=2
                           ,n.adapt=8e4)

post.samples.dis <- coda.samples(
                                  el.model.dis
                                  ,n.iter=niter
                                  ,thin=15
                                  ,variable.names=c("betavec"
                                                    ,"e.y"
                                                    ,"alpha"
                                                    ,"theta1"
                                                    , "theta2"
                                                    ,"sigma.a"
                                                    ,"sigma.t1"
                                                    ,"sigma.t2"
                                                    ,"sigma.y")
                                  )
summary.dis <- summary(post.samples.dis,quantiles=c(0.05,0.5,0.95))
residuals.dis <- post.samples.dis[[1]][,grep("e.y",colnames(post.samples.dis[[1]]))]
R2.dis <- 1-((mean(apply(residuals.dis,1,var)))/(var(y)))
summary.dis$statistics[-grep("e.y",rownames(summary.dis$statistics)),1:2]


## Interaction
model.data.inter <- subset(full.model.data
                           ,select=c(percchangew
                                     ,mandatechangeLumped
                                     ,originAgg
                                     ,nomchange
                                     ,year
                                     ,finalparty
                                     ,leadership
                                     ,distance
                                     ,security
                                     ,doublenom))
model.data.inter <- model.data.inter[complete.cases(model.data.inter),]
y <- model.data.inter$percchangew
model.data.inter$year <- as.factor(model.data.inter[,"year"])
model.data.inter$finalparty <- as.factor(model.data.inter[,"finalparty"])
origin <- ifelse(model.data.inter$originAgg==1 & model.data.inter$nomchange==0,1
                 ,ifelse(model.data.inter$originAgg==1 & model.data.inter$nomchange==1,2
                         ,ifelse(model.data.inter$originAgg==2 & model.data.inter$nomchange==0,3
                                 ,ifelse(model.data.inter$originAgg==2 & model.data.inter$nomchange==1,4
                                         ,NA))))
man.change <- model.data.inter$mandatechangeLumped
X <- model.matrix(y~.-1
                  ,subset(model.data.inter
                          , select=c(finalparty  
                            ,year
                            ,doublenom
                            ,security
                            #,leadership
                            ,distance
                            )))[,-1]
J=length(unique(origin))
N=length(y)
K=dim(X)[2]
W=diag(2)
model.data.list.inter <- list(y=y
                              ,man.change=man.change
                              ,origin=origin
                              ,X=X
                              ,N=N
                              ,J=J
                              ,K=K
                              ,W=W
                              )

model.inits.inter <- function(){
  list(
    B.raw=array(rnorm(2*J),c(J,2))
    ,mu.a.raw=rnorm(1)
    ,mu.t.raw=rnorm(1)
    ,sigma.y=runif(1)
    ,Tau.B.raw=rwish(3,diag(2))
    ,xi.a=runif(1)
    ,xi.t=runif(1)
    ,betavec=rnorm(K)
    )
}

el.model.inter <- jags.model(file="hierModelInteract.jag"
                             ,data=model.data.list.inter
                             ,inits=model.inits.inter
                             ,n.chains=2
                             ,n.adapt=niter)

post.samples.inter <- coda.samples(
  el.model.inter
  ,n.iter=niter
  ,thin=20
  ,variable.names=c("betavec"
                    ,"e.y"
                    ,"alpha"
                    ,"theta"
                    ,"sigma.a"
                    ,"sigma.t"
                    ,"sigma.y")
  )
summary.inter <- summary(post.samples.inter,quantiles=c(0.05,0.5,0.95))
residuals.inter <- post.samples.inter[[1]][,grep("e.y",colnames(post.samples.inter[[1]]))]
R2.inter <- 1-((mean(apply(residuals.inter,1,var)))/(var(y)))
summary.inter$statistics[-grep("e.y",rownames(summary.inter$statistics)),1:2]
