# Extract model parameters for post-processing for each taxonomic group, here exemplified for birds ---------------------------------

##creating objects to save results of all model runs per taxonomic group

# to save variance partitioning (VP) values
Bird_VP <- matrix(NA, 12, 7) 

# to save species level variance partitioning and model fit including AUC and Tjur values
spsVP_birds<- list() 
auxlist<- list()

# to save include species level bete paramters from HMSC object
spsBeta_birds<- list()

# to save model convergence diagnostics
diag_birds<- list()

# to save predicted communities including coordinates of site
predY_birds<-list()

# to save original raw data based communtiy 
Y_birds<-list()

# to save polynomial responses within species niches 
poly_birds<-list()

# to save polynomial responses within species niches including different thresholds (70 and 90 % of gradient) of derivatives (for sensitivity to threshhold)
poly_birds_70<-list()
poly_birds_90<-list()


## run loop over all models for birds extracting the different values
u=1

for (dec_n in c("1","2","3","4")){
  for (reg in c("SB","MB","NB")){
    for (tname in c("bird")){
      
      
      
      #modif foldnames to specifics
      fold_name<-paste(tname,dec_n,reg,sep="_")
      
      {  
        
        tnamefull<-paste0(tname,"_",dec_n,"_",reg)
        
        load(paste0(tnamefull,"/results.Rdata"))
        XData<-res[[1]]
        Y<-res[[2]]
        spn<-dim(Y)[2]
        xy<-res[[3]]
        m<-res[[4]]
        MF<- res[[7]]
        post<- res[[5]]
        
        
        ##for model convergence checks
        diag_birds[[u]]<- cbind(es.beta= effectiveSize(post$Beta),
                                gelman.diag(post$Beta,multivariate=FALSE)$psrf)
        
        print(colMeans(diag_birds[[u]]))
        
        ##1
        VP = computeVariancePartitioning(hM = m, group = group, groupnames = groupnames)
        
        Bird_VP[u,]<- c(round(100*rowMeans(VP$vals), 3) * mean(MF$TjurR2, na.rm=T), tnamefull)
        
        print(Bird_VP[u,])
        
        ##2
        for(i in 1:length(variablegroup)){
          
          ##loop to save for each variable within each model
          auxlist[[i]]<- data.frame(cbind(variable= variablegroup[i],
                                          rawvalues= round(as.numeric(VP$vals[i,]), 4),
                                          weighvalues = round(as.numeric(VP$vals[i,]) * MF$TjurR2, 4),
                                          weighvaluesAUC = round(as.numeric(VP$vals[i,]) * MF$AUC, 4),
                                          Tjur = MF$TjurR2,
                                          AUC = MF$AUC,
                                          model=tnamefull))
          
          ##then combine all the elements and save to list
          spsVP_birds[[u]]<- do.call(rbind,auxlist)
        }
        
        
        ##3
        beta=getPostEstimate(m, "Beta")
        
        ##save mean Beta values and add model name
        spsBeta_birds[[u]]<- data.frame(t(beta$mean),     ##get beta mean values per species
                                        t(beta$support),  ##get beta support values per species
                                        model=tnamefull)  ##add model name
        
        print(tnamefull)
        
        
        #4 
        pred= computePredictedValues(m,expected = FALSE)
        predY = apply(abind(pred,along=3),c(1,2),mean)
        predY_birds[[u]]<-data.frame((predY),
                                     xy[,1],
                                     xy[,2],
                                     model=tnamefull)
        
        #5
        Y_birds[[u]]<-data.frame(Y, 
                                 model=tnamefull)
        
        #6
        mpost = convertToCodaObject(m)
        Beta = as.matrix(mpost$Beta)
        
        prec_pred = as.vector(seq(min(XData$annuprec),max(XData$annuprec),length=100))
        temp_pred = as.vector(seq(min(XData$annutemp),max(XData$annutemp),length=100))
        snow_pred = as.vector(seq(min(XData$annusnowcount),max(XData$annusnowcount),length=100))
        res_prec<-c()
        res_temp<-c()
        res_snow<-c()
        for (sp in 0:(spn-1)){
          sp1<-Beta[,(sp*8+1):(sp*8+8)]
          LpredDeriv = sp1[,2]+ 2*sp1[,3]%*%t(prec_pred)
          res_prec<-rbind(res_prec,check_beta_trend(LpredDeriv))
          LpredDeriv = sp1[,4]+ 2*sp1[,5]%*%t(temp_pred)
          res_temp<-rbind(res_temp,check_beta_trend(LpredDeriv))
          LpredDeriv = sp1[,6]+ 2*sp1[,7]%*%t(snow_pred)
          res_snow<-rbind(res_snow,check_beta_trend(LpredDeriv))}
        
        
        colnames(res_prec)<-c('pos_res_prec','neg_res_prec','change_res_prec')
        colnames(res_temp)<-c('pos_res_temp','neg_res_temp','change_res_temp')
        colnames(res_snow)<-c('pos_res_snow','neg_res_snow','change_res_snow')
        
        poly_birds[[u]]<-data.frame(res_prec, 
                                    res_temp,
                                    res_snow,
                                    
                                    model=tnamefull)
        
        
        ##7 70% threshold
        mpost = convertToCodaObject(m)
        Beta = as.matrix(mpost$Beta)
        
        prec_pred_70 = as.vector(seq(min(XData$annuprec),max(XData$annuprec),length=100))
        temp_pred_70 = as.vector(seq(min(XData$annutemp),max(XData$annutemp),length=100))
        snow_pred_70 = as.vector(seq(min(XData$annusnowcount),max(XData$annusnowcount),length=100))
        res_prec_70<-c()
        res_temp_70<-c()
        res_snow_70<-c()
        for (sp in 0:(spn-1)){
          sp1<-Beta[,(sp*8+1):(sp*8+8)]
          LpredDeriv = sp1[,2]+ 2*sp1[,3]%*%t(prec_pred_70)
          res_prec_70<-rbind(res_prec_70,check_beta_trend_70(LpredDeriv))
          LpredDeriv = sp1[,4]+ 2*sp1[,5]%*%t(temp_pred_70)
          res_temp_70<-rbind(res_temp_70,check_beta_trend_70(LpredDeriv))
          LpredDeriv = sp1[,6]+ 2*sp1[,7]%*%t(snow_pred_70)
          res_snow_70<-rbind(res_snow_70,check_beta_trend_70(LpredDeriv))}
        
        
        colnames(res_prec_70)<-c('pos_res_prec','neg_res_prec','change_res_prec')
        colnames(res_temp_70)<-c('pos_res_temp','neg_res_temp','change_res_temp')
        colnames(res_snow_70)<-c('pos_res_snow','neg_res_snow','change_res_snow')
        
        poly_birds_70[[u]]<-data.frame(res_prec_70,
                                       res_temp_70,
                                       res_snow_70,
                                       
                                       model=tnamefull)
        
        
        ##8 90% threshold
        mpost = convertToCodaObject(m)
        Beta = as.matrix(mpost$Beta)
        
        prec_pred_90 = as.vector(seq(min(XData$annuprec),max(XData$annuprec),length=100))
        temp_pred_90 = as.vector(seq(min(XData$annutemp),max(XData$annutemp),length=100))
        snow_pred_90 = as.vector(seq(min(XData$annusnowcount),max(XData$annusnowcount),length=100))
        res_prec_90<-c()
        res_temp_90<-c()
        res_snow_90<-c()
        for (sp in 0:(spn-1)){
          sp1<-Beta[,(sp*8+1):(sp*8+8)]
          LpredDeriv = sp1[,2]+ 2*sp1[,3]%*%t(prec_pred_90)
          res_prec_90<-rbind(res_prec_90,check_beta_trend_90(LpredDeriv))
          LpredDeriv = sp1[,4]+ 2*sp1[,5]%*%t(temp_pred_90)
          res_temp_90<-rbind(res_temp_90,check_beta_trend_90(LpredDeriv))
          LpredDeriv = sp1[,6]+ 2*sp1[,7]%*%t(snow_pred_90)
          res_snow_90<-rbind(res_snow_90,check_beta_trend_90(LpredDeriv))}
        
        
        colnames(res_prec_90)<-c('pos_res_prec','neg_res_prec','change_res_prec')
        colnames(res_temp_90)<-c('pos_res_temp','neg_res_temp','change_res_temp')
        colnames(res_snow_90)<-c('pos_res_snow','neg_res_snow','change_res_snow')
        
        poly_birds_90[[u]]<-data.frame(res_prec_90,
                                       res_temp_90,
                                       res_snow_90,
                                       
                                       model=tnamefull)
        
        u<-u+1
      }
    }
  }
}

save(list = c("Bird_VP","spsVP_birds", "spsBeta_birds", "diag_birds", "predY_birds", "Y_birds","poly_birds","poly_birds_70", "poly_birds_90"  ), file="Birds.RData")
