#Example script to fit the Hmsc model for birds, south boreal, decade 1

library(Hmsc)

#Load data for birds decade 1 South Boreal
data<-read.csv("birds_decade_1_SB.csv")

#Site locations
xycoords<-data[,2:3]

#Enviornmental data 
XData<-data[,4:8]

#Community data
Y<-data[,9:ncol(data)]


# For large data subsets, failing to complete with available resources; we randomly subsetted 1000 records
if (dim(XData)[1]>1000){
  sub_s<-order(sample(1:dim(XData)[1],1000))
  XData<-XData[sub_s,]
  Y<-Y[sub_s,]
  xycoords<-xycoords[sub_s,]
}


#Define model study design 
studyDesign <- data.frame(Site = as.factor(row.names(XData)),Year = as.factor(XData$Year))

#Define random levels
  rL1=HmscRandomLevel(sData=xycoords)
  rL2=HmscRandomLevel(units=unique(studyDesign$Year))
  
#Define regression function 
  XFormula= ~ poly(annuprec, degree = 2, raw=TRUE) + poly(annutemp, degree = 2, raw=TRUE) + poly(annusnowcount, degree =2, raw=TRUE) + annuNAO
  
#Define model
  m <- Hmsc(Y=Y, XData=XData, XFormula=XFormula,studyDesign=studyDesign,ranLevels=list("Site"=rL1,"Year"=rL2),distr="probit")

#Define MCMC sampling parameters
  nChains= 4
  samples = 250
  thin = 100
  transient = as.integer(.5*(samples*thin))

#Run model 
a=Sys.time()
m <- sampleMcmc(m, samples = samples, thin = thin, transient = transient,
                nChains = nChains, nParallel = nChains, updater=list(GammaEta=FALSE))
print (paste('analysis_time',Sys.time()-a))

#Convert model to coda object
post=convertToCodaObject(m)

#Compute predicted occurrence probabilities 
predY=computePredictedValues(m)

#Compute model fit
MF = evaluateModelFit(hM = m, predY = predY)

#Safe results
res<-list(XData,Y,xycoords,m,post,predY,MF)
save(res,file='results.Rdata')
