

L2ErrorSegNorm<-function(EstimateD, OriginalD){
  
  
  l2errorpart<-sqrt(rowSums((EstimateD-OriginalD)^2))/sqrt(rowSums((OriginalD)^2))
  l2error<-mean(l2errorpart)
  
  l2errorInTime<-colMeans(sqrt(((EstimateD-OriginalD)^2))/sqrt(rowSums((OriginalD)^2)))
  #l2errorInTimepart<-sqrt(colSums(rowSums(EstimateD-OriginalD)^2))/sqrt(colSums((OriginalD)^2))
  return(list(c(l2error=l2error),l2errorInTime=l2errorInTime))
}

L2Error<-function(EstimateD, OriginalD){
  
  l2error<-sqrt(sum((EstimateD-OriginalD)^2))/sqrt(sum((OriginalD)^2))
  l2errorInTime<-sqrt(colSums((EstimateD-OriginalD)^2))/sqrt(colSums((OriginalD)^2))
  return(list(l2error=l2error,l2errorInTime=l2errorInTime))
}

L2ErrorAbs<-function(EstimateD, OriginalD){
  
  l2errorAbs<-1/length(EstimateD)*sqrt(sum((EstimateD-OriginalD)^2))
  l2errorInTimeAbs<-1/dim(EstimateD)[1]*sqrt(colSums((EstimateD-OriginalD)^2))
  return(list(l2errorAbs=l2errorAbs,l2errorInTimeAbs=l2errorInTimeAbs))
}

ErrorNormalized<-function(EstimateD, OriginalD){
  errornormalized<-sqrt(sum((EstimateD/sum(abs(EstimateD)) -OriginalD/sum(abs(OriginalD)) )^2))/sqrt(sum((OriginalD/sum(abs(OriginalD)))^2))
  errornormalizedtime<-sqrt(colSums((EstimateD/sum(abs(EstimateD))-OriginalD/sum(abs(OriginalD)))^2))/sqrt(colSums(abs(OriginalD/sum(abs(OriginalD)))^2))
  return(list( errornormalized= errornormalized,errornormalizedtime=errornormalizedtime))  
}

ErrorNormalizedAbs<-function(EstimateD, OriginalD){
  errornormalizedabs<-1/length(EstimateD)*sqrt(sum((EstimateD/sum(abs(EstimateD)) -OriginalD/sum(abs(OriginalD)) )^2))
  errornormalizedtimeabs<-1/dim(EstimateD)[1]*sqrt(colSums((EstimateD/sum(abs(EstimateD))-OriginalD/sum(abs(OriginalD)))^2))
  return(list( errornormalizedabs= errornormalizedabs,errornormalizedtimeabs=errornormalizedtimeabs))  
}

SmoothTheCurrents<-function(Rsmooth, morpho, MembCurrentsLine){
  #Rsmooth=1/(sqrt(2)*R
  #Smoothed membrane currents
  #how far are from eachother the neighbouring segments... smoothing based on this
  length<-numeric()
  comp.place.curr<-0
  for(i in 2:dim(morpho)[1]){
    length<-sqrt(sum((morpho[i,]-morpho[i-1,])^2))
    comp.place.curr<-c(comp.place.curr,comp.place.curr[i-1]+length)
    
  }
  memb.currents.smoothed<-array(0,dim(MembCurrentsLine))
  for (t in 1: dim(MembCurrentsLine)[2]){
    memb.currents.smoothed[,t]<-ksmooth(comp.place.curr,MembCurrentsLine[,t],"normal", bandwidth= Rsmooth,x.points=comp.place.curr)$y
  }  
  return(list(Rsmooth=Rsmooth, comp.place.curr=comp.place.curr, memb.currents.smoothed= memb.currents.smoothed ))
}


###################################################
####################################################

L1Error<-function(EstimateD, OriginalD){
  
  lerror<-sum(abs(EstimateD-OriginalD))/sum(abs(OriginalD))
  lerrorInTime<-(colSums(abs(EstimateD-OriginalD)))/(colSums(abs(OriginalD)))
  return(list(lerror=lerror,lerrorInTime=lerrorInTime))
}

L1ErrorAbs<-function(EstimateD, OriginalD){
  
  lerrorAbs<-1/length(EstimateD)*(sum(abs(EstimateD-OriginalD)))
  lerrorInTimeAbs<-1/dim(EstimateD)[1]*(colSums(abs(EstimateD-OriginalD)))
  return(list(lerrorAbs=lerrorAbs,lerrorInTimeAbs=lerrorInTimeAbs))
}

ErrorNormalizedL1<-function(EstimateD, OriginalD){
  errornormalized<-(sum(abs(EstimateD/sum(abs(EstimateD)) -OriginalD/sum(abs(OriginalD)) )))/(sum(abs(OriginalD/sum(abs(OriginalD)))))
  errornormalizedtime<-(colSums(abs(EstimateD/sum(abs(EstimateD))-OriginalD/sum(abs(OriginalD)))))/(colSums(abs(OriginalD/sum(abs(OriginalD)))))
  return(list( errornormalized= errornormalized,errornormalizedtime=errornormalizedtime))  
}

ErrorNormalizedAbsL1<-function(EstimateD, OriginalD){
  errornormalizedabs<-1/length(EstimateD)*(sum(abs(EstimateD/sum(abs(EstimateD)) -OriginalD/sum(abs(OriginalD)) )))
  errornormalizedtimeabs<-1/dim(EstimateD)[1]*(colSums(abs(EstimateD/sum(abs(EstimateD))-OriginalD/sum(abs(OriginalD)))))
  return(list( errornormalizedabs= errornormalizedabs,errornormalizedtimeabs=errornormalizedtimeabs))  
}

# ErrorinPoints<-function(EstimateD, OriginalD){
#   WhereZeros<- unique(c(which(OriginalD==0),which(EstimateD==0)))
#   error1<-sum(abs(EstimateD[-WhereZeros]- OriginalD[-WhereZeros])/abs(OriginalD[-WhereZeros]))/(dim(EstimateD)[1]*dim(EstimateD)[2]-length(WhereZeros))
#   error1intime<-colSums((EstimateD- OriginalD)^2/OriginalD^2)/(dim(EstimateD)[1])
#   
#   
# }
# 
# 
cross.valid<-function(fi, transfermatrix, compare2, ElCoords){
  N<-dim(fi)[1]
  t<-dim(fi)[2]
  M<-dim(transfermatrix)[1]
  cross.valid.error<-numeric()
  cross.valid.error.time<-numeric()
  cross.WhereisBigDiff<-array(0,c(N,t))
  fi_star<-array(0,c(N,t))
  cat(range(fi))
  cat("\n")
  
  
  
  for(electrode in 1:N){
   # if(ElCoords[electrode,1]< -450 | ElCoords[electrode,1]> 450 | 
    #   ElCoords[electrode,2]< -450 | ElCoords[electrode,2]> 450) next
    
    Epsilon = 1e-03 * max(abs(fi[electrode,])) #1e-08 
    cat(electrode)
    cat("\n")
    #skCSD_star<-array(0,c(M,t))
    skCSD_star<-transfermatrix[,-electrode]%*%fi[-electrode,]
    fi_star<-ginv(transfermatrix[compare2,])%*%skCSD_star[compare2,]
    #cross.valid.error<-c(cross.valid.error,sqrt(sum((fi-fi_star)^2))/sqrt(sum((fi)^2)))
   
    #sum(
     # abs(fi[electrode,]-fi_star[electrode,]))/(abs(fi[electrode])+epsilon)
    #)


#Epsilon = 10^{-8} * max(abs(fi[electrode]))
    
    cross.valid.error<-c(cross.valid.error,sum(abs(fi[electrode,]-fi_star[electrode,]))/(sum(abs(fi[electrode,]))))
                         #/dim(fi)[1])    
    
     #cross.valid.error<-c(cross.valid.error,sum(abs(fi[electrode,]-fi_star[electrode,])/(abs(fi[electrode,]) 
      #                                          + Epsilon ) ))
    #cross.valid.error<-c(cross.valid.error,sum(abs(fi-fi_star))/(rowSums(abs(fi)))/dim(fi)[1])    
    
    cross.valid.error.time<-c(cross.valid.error.time,(colMeans(abs(fi-fi_star)))/(colSums(abs(fi))))
    
    
    #cross.valid.error.time<-c(cross.valid.error.time,sqrt(colMeans((fi-fi_star)^2))/sqrt(colSums((fi)^2)))
    cross.WhereisBigDiff<-cross.WhereisBigDiff + sqrt(((fi-fi_star)^2))/sqrt(sum((fi)^2))
    cat( cross.valid.error[electrode])
    
    cat("\n")
  }#N
  cross.valid.error.all<-mean(cross.valid.error)
  cross.valid.error.time.all<-mean(cross.valid.error.time)
  cross.WhereisBigDiff<-cross.WhereisBigDiff/N
  return(list(cross.valid.error.all,fi_star,cross.valid.error.time.all,cross.WhereisBigDiff))
}
