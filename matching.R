#SEVENTH SET OF MODELS: MATCHING USING GENMATCH

##LUMPED MATCHING
library(Matching)
set.seed(831213)


### SMD
glmSMD <- glm(I(mandatechangeLumped==1)~doublenom
                                                  +finalparty
                                                  +year
                                                  +distance
                                                  +security
                                                  #+leadership
                             ,family=binomial()
                             ,data=subset(full.model.data,originAgg==1&is.na(percchangew)==FALSE))
p.score.smdpr <- predict(glmSMD,type="response")

XSMD <- cbind(subset(glmSMD$model,select=c(doublenom
                                           ,finalparty
                                           ,year
                                           ,distance
                                           ,security
                                           #,leadership
                                           )),p.score.smdpr)
XSMD[,"finalparty"] <- as.numeric(as.factor(XSMD[,"finalparty"]))
XSMD <- data.matrix(XSMD)
                                  
TrSMD <- full.model.data[rownames(full.model.data) %in% rownames(glmSMD$model),"mandatechangeLumped"]
YSMD <- full.model.data[rownames(full.model.data) %in% rownames(glmSMD$model),"percchangew"]
dataSMDpre <- cbind(XSMD,TrSMD)
genMSMD <- GenMatch(Tr=TrSMD,X=XSMD,pop.size=2000,wait.generations=10,max.generations=150,estimand="ATT",M=1,paired=TRUE)
matchSMD <- Match(Y=YSMD, Tr=TrSMD,X=XSMD,Weight.matrix=genMSMD, estimand="ATT")
dataSMDpost <- as.data.frame(cbind(matchSMD$mdata$Y,matchSMD$mdata$X,matchSMD$mdata$Tr))
colnames(dataSMDpost) <- c("percchangew",names(XSMD),"TrSMD")

cat("(",matchSMD$est-matchSMD$se*qnorm(0.95),",",matchSMD$est+matchSMD$se*qnorm(0.95),")",sep="")

####PR
glmPR <- glm(I(mandatechangeLumped==1)~doublenom
                                                  +finalparty
                                                  +year
                                                  +distance
                                                  +security
                                                  #+leadership
                             ,family=binomial()
                             ,data=subset(full.model.data,originAgg==2&is.na(percchangew)==FALSE))
p.score.prsmd <- predict(glmPR,type="response")

XPR <- cbind(subset(glmPR$model,select=c(doublenom
                                         ,finalparty
                                         ,year
                                         ,distance
                                         ,security
                                         #,leadership
                                         )),p.score.prsmd)
XPR$finalparty <- as.numeric(as.factor(XPR$finalparty))
XPR <- data.matrix(XPR)
TrPR <- full.model.data[rownames(full.model.data) %in% rownames(glmPR$model),"mandatechangeLumped"]
YPR <- full.model.data[rownames(full.model.data) %in% rownames(glmPR$model),"percchangew"]
dataPRpre <- cbind(XPR,TrPR)
genMPR <- GenMatch(Tr=TrPR,X=XPR,pop.size=2000,wait.generations=10,max.generations=150,estimand="ATT",M=1,paired=TRUE)
matchPR <- Match(Y=YPR, Tr=TrPR,X=XPR,Weight.matrix=genMPR,estimand="ATT")
dataPRpost <- as.data.frame(cbind(matchPR$mdata$Y,matchPR$mdata$X,matchPR$mdata$Tr))
colnames(dataPRpost) <- c("percchangew",names(XPR),"TrPR")

cat("(",matchPR$est-matchPR$se*qnorm(0.95),",",matchPR$est+matchPR$se*qnorm(0.95),")",sep="")

# ## DISAGGREGATED MATCHING
# 
# # Propensity scores
glmSmdRpr <- glm(I(mandatechange1==1)~
                        doublenom
                        +finalparty
                        +year
                        +distance
                        +security
                        #+leadership
                        ,family=binomial()
                        ,data=subset(full.model.data,originDisagg==1&is.na(percchangew)==FALSE))
p.score.SmdRpr <- predict(glmSmdRpr,type="response")

glmSmdNpr <- glm(I(mandatechange2==1)~
                      doublenom
                        +finalparty
                        +year
                        +distance
                        +security
                        #+leadership
                        ,family=binomial()
                      ,data=subset(full.model.data,originDisagg==1&is.na(percchangew)==FALSE))
p.score.SmdNpr <- predict(glmSmdNpr,type="response")

glmRprSmd <- glm(I(mandatechange1==1)~
                        doublenom
                        +finalparty
                        +year
                        +distance
                        +security
                        #+leadership
                        ,family=binomial()
                        ,data=subset(full.model.data,originDisagg==2&is.na(percchangew)==FALSE))
p.score.RprSmd <- predict(glmRprSmd,type="response")

glmNprSmd <- glm(I(mandatechange1==1)~
                    doublenom
                        +finalparty
                        +year
                        +distance
                        +security
                        #+leadership
                 ,family=binomial()
                    ,data=subset(full.model.data,originDisagg==3&is.na(percchangew)==FALSE))
p.score.NprSmd <- predict(glmNprSmd,type="response")


#SMD to RPR
XSMDRPR <- cbind(subset(glmSmdRpr$model,select=c(doublenom
                                         ,finalparty
                                         ,year
                                         ,distance
                                         ,security
                                         #,leadership
                                         )),p.score.SmdRpr)
XSMDRPR$finalparty <- as.numeric(as.factor(XSMDRPR$finalparty))
XSMDRPR <- data.matrix(XSMDRPR)
TrSMDRPR <- full.model.data[rownames(full.model.data) %in% rownames(glmSmdRpr$model),"mandatechange1"]
YSMDRPR <- full.model.data[rownames(full.model.data) %in% rownames(glmSmdRpr$model),"percchangew"]
dataSMDRPRpre <- cbind(XSMDRPR,TrSMDRPR)
genMSMDRPR <- GenMatch(Tr=TrSMDRPR,X=XSMDRPR,pop.size=2000,wait.generations=10,max.generations=150,estimand="ATT",M=1,paired=TRUE)
matchSMDRPR <- Match(Y=YSMDRPR, Tr=TrSMDRPR,X=XSMDRPR,Weight.matrix=genMSMDRPR,estimand="ATT")
dataSMDRPRpost <- as.data.frame(cbind(matchSMDRPR$mdata$Y,matchSMDRPR$mdata$X,matchSMDRPR$mdata$Tr))
colnames(dataSMDRPRpost) <- c("percchangew",names(XSMDRPR),"TrSMDRPR")
cat("(",matchSMDRPR$est-matchSMDRPR$se*qnorm(0.95),",",matchSMDRPR$est+matchSMDRPR$se*qnorm(0.95),")",sep="")


#SMD to NPR
XSMDNPR <- cbind(subset(glmSmdNpr$model,select=c(doublenom
                                         ,finalparty
                                         ,year
                                         ,distance
                                         ,security
                                         #,leadership
                                         )),p.score.SmdNpr)
XSMDNPR$finalparty <- as.numeric(as.factor(XSMDNPR$finalparty))
XSMDNPR <- data.matrix(XSMDNPR)
TrSMDNPR <- full.model.data[rownames(full.model.data) %in% rownames(glmSmdNpr$model),"mandatechange2"]
YSMDNPR <- full.model.data[rownames(full.model.data) %in% rownames(glmSmdNpr$model),"percchangew"]
dataSMDNPRpre <- cbind(XSMDNPR,TrSMDNPR)
genMSMDNPR <- GenMatch(Tr=TrSMDNPR,X=XSMDNPR,pop.size=2000,wait.generations=10,max.generations=150,estimand="ATT",M=1,paired=TRUE)
matchSMDNPR <- Match(Y=YSMDNPR, Tr=TrSMDNPR,X=XSMDNPR,Weight.matrix=genMSMDNPR,estimand="ATT")
dataSMDNPRpost <- as.data.frame(cbind(matchSMDNPR$mdata$Y,matchSMDNPR$mdata$X,matchSMDNPR$mdata$Tr))
colnames(dataSMDNPRpost) <- c("percchangew",names(XSMDNPR),"TrSMDNPR")
cat("(",matchSMDNPR$est-matchSMDNPR$se*qnorm(0.95),",",matchSMDNPR$est+matchSMDNPR$se*qnorm(0.95),")",sep="")

#RPR to SMD
XRPRSMD <- cbind(subset(glmRprSmd$model,select=c(doublenom
                                         ,finalparty
                                         ,year
                                         ,distance
                                         ,security
                                         #,leadership
                                         )),p.score.RprSmd)
XRPRSMD$finalparty <- as.numeric(as.factor(XRPRSMD$finalparty))
XRPRSMD <- data.matrix(XRPRSMD)
TrRPRSMD <- full.model.data[rownames(full.model.data) %in% rownames(glmRprSmd$model),"mandatechange1"]
YRPRSMD <- full.model.data[rownames(full.model.data) %in% rownames(glmRprSmd$model),"percchangew"]
dataRPRSMDpre <- cbind(XRPRSMD,TrRPRSMD)
genMRPRSMD <- GenMatch(Tr=TrRPRSMD,X=XRPRSMD,pop.size=2000,wait.generations=10,max.generations=150,estimand="ATT",M=1,paired=TRUE)
matchRPRSMD <- Match(Y=YRPRSMD, Tr=TrRPRSMD,X=XRPRSMD,Weight.matrix=genMRPRSMD,estimand="ATT")
dataRPRSMDpost <- as.data.frame(cbind(matchRPRSMD$mdata$Y,matchRPRSMD$mdata$X,matchRPRSMD$mdata$Tr))
colnames(dataRPRSMDpost) <- c("percchangew",names(XRPRSMD),"TrRPRSMD")
cat("(",matchRPRSMD$est-matchRPRSMD$se*qnorm(0.95),",",matchRPRSMD$est+matchRPRSMD$se*qnorm(0.95),")",sep="")



#NPR to SMD
XNPRSMD <- cbind(subset(glmNprSmd$model,select=c(doublenom
                                         ,finalparty
                                         ,year
                                         ,distance
                                         ,security
                                         #,leadership
                                         )),p.score.NprSmd)
XNPRSMD$finalparty <- as.numeric(as.factor(XNPRSMD$finalparty))
XNPRSMD <- data.matrix(XNPRSMD)
TrNPRSMD <- full.model.data[rownames(full.model.data) %in% rownames(glmNprSmd$model),"mandatechange1"]
YNPRSMD <- full.model.data[rownames(full.model.data) %in% rownames(glmNprSmd$model),"percchangew"]
dataNPRSMDpre <- cbind(XNPRSMD,TrNPRSMD)
genMNPRSMD <- GenMatch(Tr=TrNPRSMD,X=XNPRSMD,pop.size=2000,wait.generations=10,max.generations=150,estimand="ATT",M=1,paired=TRUE)
matchNPRSMD <- Match(Y=YNPRSMD, Tr=TrNPRSMD,X=XNPRSMD,Weight.matrix=genMNPRSMD,estimand="ATT")
dataNPRSMDpost <- as.data.frame(cbind(matchNPRSMD$mdata$Y,matchNPRSMD$mdata$X,matchNPRSMD$mdata$Tr))
colnames(dataNPRSMDpost) <- c("percchangew",names(XNPRSMD),"TrNPRSMD")
cat("(",matchNPRSMD$est-matchNPRSMD$se*qnorm(0.95),",",matchNPRSMD$est+matchNPRSMD$se*qnorm(0.95),")",sep="")


### Everything
#Aggregated
cat("(",round((matchSMD$est-matchSMD$se*qnorm(0.95)),3),",",round((matchSMD$est+matchSMD$se*qnorm(0.95)),3),")",sep="")
cat("(",round((matchPR$est-matchPR$se*qnorm(0.9)),3),",",round((matchPR$est+matchPR$se*qnorm(0.9)),3),")",sep="")

#Disaggregated
cat("(",round((matchSMDRPR$est-matchSMDRPR$se*qnorm(0.9)),3),",",round((matchSMDRPR$est+matchSMDRPR$se*qnorm(0.9)),3),")",sep="")
cat("(",round((matchSMDNPR$est-matchSMDNPR$se*qnorm(0.9)),3),",",round((matchSMDNPR$est+matchSMDNPR$se*qnorm(0.9)),3),")",sep="")
cat("(",round((matchRPRSMD$est-matchRPRSMD$se*qnorm(0.9)),3),",",round((matchRPRSMD$est+matchRPRSMD$se*qnorm(0.9)),3),")",sep="")
cat("(",round((matchNPRSMD$est-matchNPRSMD$se*qnorm(0.9)),3),",",round((matchNPRSMD$est+matchNPRSMD$se*qnorm(0.9)),3),")",sep="")

 



