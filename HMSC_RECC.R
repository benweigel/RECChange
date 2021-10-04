library(BayesLogit)
library(Hmsc)
library(tidyr)
library(reshape2)
library(plyr)
library(abind)
library(viridis)
library(ggplot2)

load('input.Rdata')
XData<-input_data[1][[1]]
Y<-input_data[2][[1]]
xycoords<-input_data[3][[1]]


if (dim(XData)[1]>1000){
  sub_s<-order(sample(1:dim(XData)[1],1000))
  XData<-XData[sub_s,]
  Y<-Y[sub_s,]
  xycoords<-xycoords[sub_s,]
}


row.names(XData)<-1:dim(XData)[1]
row.names(Y)<-1:dim(XData)[1]
row.names(xycoords)<-1:dim(XData)[1]


n<-dim(XData)[1]	#number of sites
studyDesign <- data.frame(Site = as.factor(row.names(XData)),Year = as.factor(XData$Year))

nChains= 4
samples = 250
thin = 100
transient = as.integer(.5*(samples*thin))

#Define random levels
if (grepl('forest',getwd())){
    rL1=HmscRandomLevel(sData=xycoords)
    XFormula= ~ poly(annuprec, degree = 2, raw=TRUE) + poly(annutemp, degree = 2, raw=TRUE) + poly(annusnowcount, degree =2, raw=TRUE)
    m <- Hmsc(Y=Y, XData=XData, XFormula=XFormula,studyDesign=studyDesign,ranLevels=list("Site"=rL1),distr="probit")
  } else {
    rL1=HmscRandomLevel(sData=xycoords)
    rL2=HmscRandomLevel(units=unique(studyDesign$Year))
    XFormula= ~ poly(annuprec, degree = 2, raw=TRUE) + poly(annutemp, degree = 2, raw=TRUE) + poly(annusnowcount, degree =2, raw=TRUE) + annuNAO
    m <- Hmsc(Y=Y, XData=XData, XFormula=XFormula,studyDesign=studyDesign,ranLevels=list("Site"=rL1,"Year"=rL2),distr="probit")
  }


a=Sys.time()
m <- sampleMcmc(m, samples = samples, thin = thin, transient = transient,
nChains = nChains, nParallel = nChains, updater=list(GammaEta=FALSE))
print (paste('analysis_time',Sys.time()-a))
post=convertToCodaObject(m)
predY=computePredictedValues(m)
MF = evaluateModelFit(hM = m, predY = predY)

res<-list(XData,Y,xycoords,m,post,predY,MF)
save(res,file='results.Rdata')
