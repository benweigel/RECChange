# Niche position functions
# polynomial repsones of HMSC beta parameters

check_beta_trend<-function(LpredDeriv){
  pn<-dim(LpredDeriv)[2]
  response = matrix(NA,nrow = nrow(LpredDeriv),ncol = 1)
  for (i in 1:nrow(LpredDeriv)){
    if (LpredDeriv[i,1]>0 & LpredDeriv[i,pn]>0){   
      # derivative is positive at both ends -> clear case 1
      response[i] = 1 
    } else if (LpredDeriv[i,1]<0 & LpredDeriv[i,pn]<0) {
      # derivative is negative at both ends -> clear case 2
      response[i] = 2
    } else if (LpredDeriv[i,1]>0 & LpredDeriv[i,pn]<0) {
      # derivative changes from postitive to negative -> have to check in more detail
      
      if (sum(LpredDeriv[i,]>0)/pn<0.2){
        # derivative is positive over less than 20% of the x-gradient extent -> classify as case 2
        response[i] = 2
        
      } else if (sum(LpredDeriv[i,]>0)/pn>0.8){
        # derivative is positive over more than 80% of the x-gradient extent -> classify as case 1
        response[i] = 1
        
      } else {
        # derivative is postive or negative only over at most 60% of the x-gradient extent -> classify as case 3
        response[i] = 3
      }
      
    } else if (LpredDeriv[i,1]<0 & LpredDeriv[i,pn]>0) {
      # derivative changes from negative to positive -> have to check in more detail
      
      if (sum(LpredDeriv[i,]<0)/pn<0.2){
        # derivative is negative over less than 20% of the x-gradient extent -> classify as case 1
        response[i] = 1
        
      } else if (sum(LpredDeriv[i,]<0)/pn>0.8){
        # derivative is negative over more than 80% of the x-gradient extent -> classify as case 2
        response[i] = 2
        
      } else {response[i] = 3 
      # derivative is positive or negative only over at most 60% of the x-gradient extent 
      
      }
    }
    
  }
  # calculate the probabilities for the cases
  Pr1 = sum(response==1,na.rm=T)/length(response) #positive; derivative is positive over more than 90%
  Pr2 = sum(response==2,na.rm=T)/length(response) #negative, derivative is negative over more than 90%
  Pr3 = sum(response==3,na.rm=T)/length(response) #change; derivative is positive or negative only over at most 60% 
  
  return(c(Pr1,Pr2,Pr3))}