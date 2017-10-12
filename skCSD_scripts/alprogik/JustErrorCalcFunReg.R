#source(paste0(ScriptLocation,"alprogik/JustErrorCalcFun.R"))
#JustErrorCalcFun(SimPath,locationsKernel,locationsData,SNR, lammbda)
M<-512
lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" ) #c(0.1,0.01,0.001)#
R.all<-2^(3:7) #c(8,16,128)#


JustErrorCalcFunReg<-function(SimPath,locationsKernel,locationsData,SNR){
Reg<-TRUE

ErrorResults <- vector("list", length(locationsKernel))
ErrorResultsL1 <- vector("list", length(locationsKernel))
SmoothR<- vector("list", length(locationsKernel))
for(Dlength in 1:length(locationsKernel)){

inname<-paste0(SimPath, locationsKernel[Dlength],"/")
outname<-paste0(SimPath, locationsData[Dlength],"/")

if(file.exists(paste0(outname,outname1))==FALSE) dir.create(paste0(outname,outname1))


lambda<-0
#SNR<-0#10 #0 means no noise



#readin in data independent from R and M
sigma<-as.matrix(read.table(paste(inname,"elprop",sep="")))[2]
ElCoords<-as.matrix(read.table(paste(inname,"elcoord_x_y_z",sep="")))
ElCoords<-matrix(ElCoords,ncol=3)
DistEl<-as.matrix(dist(ElCoords))
diag(DistEl)<-Inf
SmoothElWidth<-mean(apply(DistEl,1,min))

membcurr <-read.table(paste(outname,"/membcurr",sep=""))
seg.length<-as.matrix(read.table(paste(inname,"seglength",sep="")))
#cat(paste(outname,"/membcurr",sep=""))
seg.nb<-length(seg.length)
funaramvonal<-function(x) x/seg.length
memb.currents.line<-array(0, c(dim(membcurr)))
memb.currents.line<-apply(membcurr,2,funaramvonal) 
morpho<-as.matrix(read.table(paste0(inname,'segcoordinates.txt')))

############# !!!!!!!!!!!!!!!
SmoothPar<- 30 #as.matrix(read.table(paste0(SimPath,"/SmoothinKernelWidth")))
############### !!!!!!!!!!!!!

#BandwidtAll<-c(sqrt(SmoothPar[2]/SmoothPar[3]),2*sqrt(SmoothPar[2]/SmoothPar[3]), SmoothElWidth/2, SmoothElWidth/4, sqrt(SmoothElWidth)/2 , 10,30,50,70,90)
BandwidtAll<-c(R.all)# 10,,20,30,,50,70,90)

BandwidtAllName<-c(paste("Smoothed",c(1:length(BandwidtAll))))




EL2ErrorVarR<-array(0,c((length(BandwidtAll)+1)*length(R.all),length(lambda.all)))
EL1ErrorVarR<-array(0,c((length(BandwidtAll)+1)*length(R.all),length(lambda.all)))

El2ErrorRaw<-array(0,c(length(R.all),length(lambda.all)))
El1ErrorRaw<-array(0,c(length(R.all),length(lambda.all)))
#Boxerror<-array(0,c(length(R.all),length(lambda.all))) #lets try 3 at first
#Reading in the LFP and adding noise
    LFP<-as.matrix(read.table(paste(outname,"/myLFP",sep="")))
    
if (SNR!=0){
      LFPvariance<-var(c(LFP))
      Noisevariance<-LFPvariance/SNR
      NoiseLFP<-array(rnorm(length(c(LFP)), mean = 0, sd = sqrt(Noisevariance)),dim(LFP))
      LFP<-LFP+NoiseLFP
    }
LFPOriginal<-LFP

for(lamb in 1:length(lambda.all)){
  for(r in 1:length(R.all)){
    LFP<-LFPOriginal
    lambda<-lambda.all[lamb]
    R<-R.all[r]
 el2ignorename<-paste(inname,outname1,"/El2Ignore_M",M,"_R",R,sep="")
    if(file.exists(el2ignorename)){
      El2Ignore<-c(as.matrix(read.table(el2ignorename)))
      LFP<-LFP[-El2Ignore,] }
      
       if(file.exists(paste0(outname,outname1,"/skCSDall_M",M,"_R",R, "lambda",lambda))) {
      skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
 skCSD.all.part<-as.matrix(read.table(paste0(outname,outname1,"/skCSDall_M",M,"_R",R, "lambda",lambda)))
    #Reading in tge currents
    }
    if( file.exists(paste(outname,outname1,"/ksCSD_Matr",M,"_R",R,"lambda",
                                                   lambda.all[lamb], "SNR",SNR,sep=""))) {
      skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
 skCSD.all.part<-(4*pi*0.5)^2*as.matrix(read.table( paste(outname,outname1,"/ksCSD_Matr",M,"_R",R,"lambda",
                                                   lambda.all[lamb], "SNR",SNR,sep="")))
    #Reading in tge currents
    }
    
 
 
 SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

 for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
    }
 
 
  

# image.plot(array(rep(seg.length,561),dim(skCSD.all))*skCSD.all)

 # image.plot(as.matrix(membcurr))

 #boxerror
 #CurrNew<-KernSmoothDistance(30,inname, outname, outname1,"10")
  # image.plot(as.matrix(CurrNew),main="Smoothed Ground truth")
  # image.plot(skCSD.all,main="Smoothed Ground truth")
 #  L1Error(skCSD.all,membcurr)[[1]]
   
# Boxerror[r,lamb]<-L1Error(skCSD.all%*%seg.length,membcurr)
 
 
    #skCSD.Smoothed<-KernSmoothDistanceOther(SmoothElWidth/4,inname, outname, outname1,"skCSD",skCSD.all)
    #######################Calculate Currents
    ############################################
    #Calculate Errors
 
 
 
 
for(bw in 1:length(BandwidtAll)){

memb.currents.smoothed<-KernSmoothDistance(BandwidtAll[bw],inname, outname, outname1,BandwidtAllName[bw])
EL2ErrorVarR[r+(bw-1)*length(R.all),lamb]<-L2Error(skCSD.all,memb.currents.smoothed)[[1]]
EL1ErrorVarR[r+(bw-1)*length(R.all),lamb]<-L1Error(skCSD.all,memb.currents.smoothed)[[1]]
#EL2ErrorVarR[r+(bw-1)*length(R),lamb]<-L2Error(skCSD.Smoothed,memb.currents.smoothed)[[1]]
#image(t(memb.currents.smoothed),col=rainbow(40))
}

EL2ErrorVarR[r+length(BandwidtAll)*length(R.all),lamb]<-L2Error(skCSD.all,memb.currents.line)[[1]]
EL1ErrorVarR[r+length(BandwidtAll)*length(R.all),lamb]<-L1Error(skCSD.all,memb.currents.line)[[1]]
El2ErrorRaw[r,lamb]<-L2Error(skCSD.all,memb.currents.line)[[1]]

#EL2ErrorVarR[r,lamb]<-L2Error(skCSD.all,memb.currents.line)[[1]]
El1ErrorRaw[r,lamb]<-L1Error(skCSD.all,memb.currents.line)[[1]]
######################CV error 


#Kmatr<-as.matrix(read.table(paste(inname,outname1,"/K_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
 #   Ktildematr<-as.matrix(read.table(paste(inname,outname1,"/Ktilda_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
      
  #    Tmatr<-try((4*pi*0.5)*t(Ktildematr)%*%ginv(Kmatr)) #!!!!!
      
#EL2ErrorVarR[r+length(BandwidtAll)*length(R.all),lamb]<-try(cross.valid(LFP,Tmatr,1:dim(skCSD.all)[1])[[1]],silent=TRUE)


   }
   }
write.table(El2ErrorRaw, paste(outname,outname1,"/ErrorL2Raw_Noise","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
write.table(El1ErrorRaw, paste(outname,outname1,"/ErrorL1Raw_Noise","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)

write.table(EL2ErrorVarR, paste(outname,outname1,"/EL2ErrorVarR","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
write.table(EL1ErrorVarR, paste(outname,outname1,"/EL1ErrorVarR","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)



    ErrorResults[[Dlength]] <-EL2ErrorVarR
    ErrorResultsL1[[Dlength]] <-EL1ErrorVarR
SmoothR[[Dlength]]<-BandwidtAll
    }
return(list(ErrorResults=ErrorResults,SmoothR=SmoothR, ErrorResultsL1=ErrorResultsL1 ))

}


#####################JUST CV ERROR

sigma<-0.3
JustCVReg<-function(SimPath,locationsKernel,locationsData,SNR){
    Reg<-TRUE
    ErrorResults <- vector("list", length(locationsKernel))
CVConfig<-vector("list", length(R.all)*length(lambda.all))

    for(Dlength in 1:length(locationsKernel)){
        inname<-paste0(SimPath, locationsKernel[Dlength],"/")
        outname<-paste0(SimPath, locationsData[Dlength],"/")
        
        if(file.exists(paste0(outname,outname1))==FALSE) dir.create(paste0(outname,outname1))

        SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
        seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
        LFP<-as.matrix(read.table(paste(outname,"/myLFP",sep="")))
        
if (SNR!=0){
      LFPvariance<-var(c(LFP))
      Noisevariance<-LFPvariance/SNR
      NoiseLFP<-rnorm(length(c(LFP)), mean = 0, sd = sqrt(Noisevariance))
      LFP<-LFP+NoiseLFP
    }
LFPOriginal<-LFP
        
        ElCoords<-as.matrix(read.table(paste(inname,"elcoord_x_y_z",sep="")))
        ElCoords<-matrix(ElCoords,ncol=3)
        Xel<-unique(ElCoords[,2])
        XelNb<-length(Xel)
        Yel<-unique(ElCoords[,3])
        
        ErrorCV<-array(0, c(length(R.all), length(lambda.all)))
        par(mfrow=c(length(lambda.all),length(R.all)))
for(lamb in 1:length(lambda.all)){
  for(r in 1:length(R.all)){
  cvErrorOut<-numeric()
 
  LFP<-LFPOriginal
    lambda<-lambda.all[lamb]
    R<-R.all[r]
#     if(file.exists(paste0(outname,outname1,"/skCSDall_M",M,"_R",R, "lambda",lambda))==FALSE) next
#       skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
#  skCSD.all.part<-as.matrix(read.table(paste0(outname,outname1,"/skCSDall_M",M,"_R",R, "lambda",lambda)))
#     #Reading in tge currents
#  
#  for(i in 1: seg.nb){
#       sameplace<-numeric() 
#       sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
#       #Plot the current related ti the different segments
#     #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
#       if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
#       if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
#     }
#  
  
  
  
  
 el2ignorename<-paste(inname,outname1,"/El2Ignore_M",M,"_R",R,sep="")
    if(file.exists(el2ignorename)){
      El2Ignore<-c(as.matrix(read.table(el2ignorename)))
      LFP<-LFP[-El2Ignore,] }
      


Kmatr<-as.matrix(read.table(paste(inname,outname1,"/K_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
    Ktildematr<-as.matrix(read.table(paste(inname,outname1,"/Ktilda_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
      #try
      Tmatr<-(4*pi*sigma)*t(Ktildematr)%*%ginv(Kmatr) 
      


#skCSD.all.part<-Tmatr%*%LFP
    


# if (SNR!=0) {    write.table(Tmatr%*%LFP, paste(outname,outname1,"/ksCSD_Matr",M,"_R",R,"lambda",
#                                                   lambda.all[lamb], "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)}
      
      
      cvErrorOut<-cross.valid(LFP,Tmatr,1:dim(Tmatr)[1],ElCoords)
ErrorCV[r,lamb]<-cvErrorOut[[1]]
cat("Numb to rteplace:", length(R.all)*(lamb-1)+ r)
#CVConfig[[length(R.all)*(lamb-1)+ r]]<-as.matrix(cvErrorOut[[4]])
#write.table(as.matrix(cvErrorOut[[4]]),paste(outname,outname1,"/CVContribution_M",M,"_R",R,"lambda",lambda.all[lamb], 
#                                             "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
#ErrorContribution
}} #lambda and R
#        write.table(LFP, paste(outname,outname1,"/LFP_Noise","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
        
#ErrorResults[[Dlength]]<-ErrorCV
write.table(ErrorCV, paste(outname,outname1,"/ErrorCV_NoiseElec","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)

#MinIndex<-which(ErrorCV==min(ErrorCV),arr.ind=TRUE)
#write.table(MinIndex, paste(outname,outname1,"/MinIndex_Noise","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)

#CVBest<-CVConfig[[length(R.all)*(MinIndex[2]-1)+ MinIndex[1]]]
# 
# png(paste(inname,outname1,"/CVPlot_SNR",SNR,".png",sep=""), width=600, height=800)
# par(mfrow=c(3,2))
# image(Xel,Yel, matrix(LFP[,60], nrow=XelNb), main="LFP at 60")
# image(Xel,Yel, matrix(CVBest[,60], nrow=XelNb),main="CV at 60")
# image(Xel,Yel, matrix(LFP[,150], nrow=XelNb), main="LFP at 150")
# image(Xel,Yel, matrix(CVBest[,150], nrow=XelNb),main="CV at 150")
# image(Xel,Yel, matrix(LFP[,300], nrow=XelNb), main="LFP at 300")
# image(Xel,Yel, matrix(CVBest[,300], nrow=XelNb),main="CV at 300")
# dev.off()
#write.table(as.matrix(cvErrorOut[[4]]),paste(outname,"CVContribution_M",M,"_R",R,"lambda",lambda.all[lamb],"SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
}#Dlength
#return(list(ErrorResults, CVConfig))
}

break
CVError2<-JustCVReg(SimPath,locationsKernel,locationsData,SNR)
#Cv error for ball & stick + ganglion

Calchere<-0
if(Calchere==1){
SNR<-0
CVError<-JustCVReg(SimPath,locationsKernel,locationsData,SNR)
SNRss<-c(64,16,4,1)
for (i in 1:5){
cat(i)
JustCVReg(SimPath,locationsKernel,locationsData,SNRss[i])
#CVError<-JustCVReg(SimPath,locationsKernel,locationsData,SNR[i])

}
image(1:5, 1:5,CVError[[1]][[24]], col=rainbow(200),xaxt="n",yaxt="n")

matplot(R.all,t(CVError[[1]][[24]]), t='l')
axis(1, at=c(1:5),labels=R.all)
axis(2, at=c(1:5),labels=lambda.all)


matplot(R.all, t(CVError[[1]][[1]]),t='l')
#################JUST CV ERROR END









alma1<-JustErrorCalcFunReg(SimPath,locationsKernel,locationsData,0)



#break

for(iter1 in 1: length(locationsKernel)){
filename1<-paste0(SimPath,locationsKernel[iter1],"/ErrorTable")
filename2<-paste0(SimPath,locationsKernel[iter1],"/SmoothingWidth")
write.table(alma1[[3]][[iter1]], filename1,row.names=FALSE,col.names=FALSE)
write.table(alma1[[2]][[iter1]], filename2,row.names=FALSE,col.names=FALSE)
}


##############
#Plot l2
MinError<-numeric()

for(iter1 in 1: length(locationsData)){
Error2Plot<-as.matrix(read.table(paste0(SimPath,locationsData[iter1],outname1,"/ErrorL1Raw_Noise_SNR0")))
png(paste0(SumPlotPlace,"/",locationsKernel[iter1],"_L1RawRvsLambda.png"),width=500, height=500)
#par(xpd=FALSE)
par(cex=1.5)
matplot(R.all,Error2Plot,t="l",lwd=2, xlab="Basis width", ylab="L2 Error")
legend("bottomright",c(paste(lambda.all)),pch=20,col=1:5,title="lambda")
dev.off()

MinError<-c(MinError, min(Error2Plot))

}

MinErrorArray<-matrix(MinError,nrow=6)
png(paste0(SumPlotPlace, "Y_L1vsElnumb.png"),width=500, height=500)
#par(xpd=FALSE)
par(cex=1.5)
matplot(unique(ElNumb),t(MinErrorArray),pch=20,lwd=2, xlab="Number of Electrodes", ylab="L1 Error")
points(unique(ElNumb),colMeans(MinErrorArray),col="GRAY",pch=20)
legend("topright",c("Grid", "Rand 1","Rand 2","Rand 3","Rand 4","Rand 5","Mean"),pch=20,col=c(1:6,"GRAY"),title="Setup")
dev.off()


#Yrot:
png(paste0(SumPlotPlace, "BS_L1vsElnumb.png"),width=500, height=500)
#par(xpd=FALSE)
par(cex=1.5)
plot(ElNumb,MinError,t="b",lwd=2, xlab="Number of Electrodes", ylab="L1 Error")
dev.off()
#################### Reading in
ErrorResults <- numeric()


for(iter1 in 1: length(locationsKernel)){
png(paste0(SumPlotPlace,"/CampareErrorL1Smoothed_",locationsKernel[iter1],".png"),width=700, height=600)
par(xpd=FALSE, mar=c(4,4,4,7))
Which2Plot<-iter1
SmoothWidth<-as.matrix(read.table(paste0(SimPath,locationsKernel[Which2Plot],"/SmoothingWidth")))

Error2Plot<-as.matrix(read.table(paste0(SimPath,locationsKernel[Which2Plot],"/ErrorTable")))
matplot(Error2Plot[1:length(R.all),],t="l",ylim=range(Error2Plot),xlim=c(1,length(R.all)*(length(SmoothWidth)+1)),xaxt="n",xlab="Basis Width", ylab="Relative Error",lwd=1.5)
for(sw in 2:(length(SmoothWidth)+1)){
matlines(length(R.all)*(sw-1)+ 1:length(R.all) ,Error2Plot[((sw-1)*length(R.all) +1):(sw*length(R.all)),],lwd=1.5)
}
axis(1, 1:dim(Error2Plot)[1], rep(R.all,length(SmoothWidth)+1))

#ErrorResults<-rbind(ErrorResults, Error2Plot[61:66,])
mtext("Smoothing Width", side=3, line=0)
mtext(c(round(SmoothWidth,3), " - "), side=3, at=c(1:(length(SmoothWidth)+1)*length(R.all)-length(R.all)/2),line=-1)
  abline(v=c((1:length(SmoothWidth))*length(R.all)+0.5),col="ORANGE")
  par(xpd=TRUE)
  legend("topright", inset=c(-0.15,0),c(paste(lambda.all)),pch=20,col=1:5,title="lambda",bg="WHITE")
dev.off()
  }


#######################################
#Plot comparable plots with the one of CV Errorss




#Number of electrodes vs minimal error
MinError <- numeric()
MinErrorR<-numeric()
MinErrorL<-numeric()
for(iter1 in 1: length(locationsKernel)){
Error2Plot<-as.matrix(read.table(paste0(SimPath,locationsKernel[iter1],"/ErrorTable")))
#32-es símításos
CurrMin<- min(Error2Plot[(2*length(R.all)+1):(3*length(R.all)),])
IndecesMin<-which(Error2Plot[(2*length(R.all)+1):(3*length(R.all)),]==CurrMin,arr.ind=TRUE)
MinErrorR<-c(MinErrorR,R.all[IndecesMin[1]])
MinErrorL<-c(MinErrorL,lambda.all[IndecesMin[2]])
MinError<-c(MinError,CurrMin)
  } 
png(paste0(SumPlotPlace,"/ElnumbvsMinSmoothed_",CellType,".png"))
par(mfrow=c(1,1), cex=1.2)
#Which2Plot<-iter1
#SmoothWidth<-as.matrix(read.table(paste0(SimPath,locationsKernel[iter1],"/SmoothingWidth")))


if(CellType!="gangNew") {
 plot(1:length(MinError),MinError,pch=4, 
xlab="Number of Electrodes", ylab="L1 Error", cex=2, lwd=2,t="b",xaxt="n")
axis(1, at=c(1:4),labels=c(8,16,32,64))
}
if(CellType=="gangNew"){ 
 MinError1<-matrix(MinError,nrow=2)
 matplot(c(128,256),MinError1,pch=4, 
xlab="Number of Electrodes", ylab="L1 Error", cex=2, lwd=2,t="b",xaxt="n")
axis(1, at=c(128,256),labels=c(128,256))
legend("topright", c("Close", "Sequential", "Random"), lty=1:3,lwd=2,col=1:3)
}
#legend("topright",c("Grid ",paste("Rand",1:5)),pch=20,col=1:6,title="El Geom")
dev.off()

#if celltype Y


MinError <- numeric()
MinErrorR<-numeric()
MinErrorL<-numeric()
for(iter1 in 1: length(locationsKernel)){
Error2Plot<-as.matrix(read.table(paste0(SimPath,locationsKernel[iter1],"/ErrorTable")))
#32-es símításos
CurrMin<- min(Error2Plot[(2*length(R.all)+1):(3*length(R.all)),])
IndecesMin<-which(Error2Plot[(2*length(R.all)+1):(3*length(R.all)),]==CurrMin,arr.ind=TRUE)
MinErrorR<-c(MinErrorR,R.all[IndecesMin[1]])
MinErrorL<-c(MinErrorL,lambda.all[IndecesMin[2]])
MinError<-c(MinError,CurrMin)
}
MinError<-matrix(MinError,nrow=6)
MinErrorR<-matrix(MinErrorR,nrow=6)  
MinErrorL<-matrix(MinErrorL,nrow=6)  

png(paste0(SumPlotPlace,"/ElnumbvsMinSmoothed_",CellType,".png"))
par(mfrow=c(1,1), cex=1.2)

SmoothWidth<-as.matrix(read.table(paste0(SimPath,locationsKernel[iter1],"/SmoothingWidth")))
matplot(unique(ElNumb),t(MinError),pch=20,lwd=2, xlab="Number of Electrodes", ylab="L1 Error")
points(unique(ElNumb),colMeans(MinError),col="GRAY",pch=20)
legend("topright",c("Grid", "Rand 1","Rand 2","Rand 3","Rand 4","Rand 5","Mean"),pch=20,col=c(1:6,"GRAY"),title="Setup")

 
#axis(1, at=c(1:4),labels=c(8,16,32,64))
#legend("topright",c("Grid ",paste("Rand",1:5)),pch=20,col=1:6,title="El Geom")
dev.off()

#plot error and optimal values together


png(paste0(SumPlotPlace,"/ElnumbvsParams_",CellType,".png"),width=800)
par(mfrow=c(1,3), cex=1.2)
SmoothWidth<-as.matrix(read.table(paste0(SimPath,locationsKernel[iter1],"/SmoothingWidth")))
matplot(unique(ElNumb),t(MinError),pch=20,lwd=2, xlab="Number of Electrodes", ylab="L1 Error")
points(unique(ElNumb),colMeans(MinError),col="GRAY",pch=20)
legend("topright",c("Grid", "Rand 1","Rand 2","Rand 3","Rand 4","Rand 5","Mean"),pch=20,col=c(1:6,"GRAY"),title="Setup")

 matplot(unique(ElNumb),t(MinErrorR),pch=20,lwd=2, xlab="Number of Electrodes", ylab="Optimal R")
  matplot(unique(ElNumb),t(MinErrorL),pch=20,lwd=2, xlab="Number of Electrodes", ylab="Optimal lambda")
#axis(1, at=c(1:4),labels=c(8,16,32,64))
#legend("topright",c("Grid ",paste("Rand",1:5)),pch=20,col=1:6,title="El Geom")
dev.off()















#matplot(ErrorResults, pch=10, cex=1.5, ylab="CV Error", xlab="R", xaxt="n", main="Estimation Error of CSD on Y-Rotated")
matplot(ErrorResults, pch=10, cex=1.5, ylab="CV Error", xlab="R", xaxt="n", main="Estimation Error of CSD on Gang")
axis(1, at=1:(4*length(R.all)), labels=rep(R.all,4))
#axis(3, at=(1:4)*length(R.all)-0.5*length(R.all), labels=c("4x2","4x4", "4x8", "4x16"),tick=FALSE,line=-1)
#axis(3, at=(1:6)*length(R.all)-0.5*length(R.all), labels=c("4","8", "16", "32", "64", "128"),tick=FALSE,line=-1)
axis(3, at=(1:4)*length(R.all)-0.5*length(R.all), labels=c("8x8","6x5", "12x10a", "12x10b"),tick=FALSE,line=-1)

abline(v=c((1:(length(locationsKernel)-1))*length(R.all)+0.5),col=1)
legend("topleft",c(paste(lambda.all)),pch=20,col=1:5,title="lambda")


Which2Plot<-5
matplot(alma[[1]][[Which2Plot]],t="l",ylim=c(0,2))
 #abline(v=c((1:3)*length(lambda.all)+0.5))
#axis(side=1, labels=rep(lambda.all,4),at=c(1:(length(lambda.all)*4)))
mtext(round(alma[[2]][[Which2Plot]],3), side=3, at=c(1:length(alma[[2]][[Which2Plot]])*length(lambda.all)-length(lambda.all)/2),line=-1)
  abline(v=c((1:9)*length(lambda.all)+0.5),col="ORANGE")
  #par(xpd=TRUE)
  #legend(4*length(R.all)+2,median(whichError2Plot),c(paste(lambda.all)),pch=20,col=1:6,title="M")

#BandwidtAll

}

#Dicrete inverse problems





#####################JUST L1 ERROR

sigma<-0.3
JustL1Reg<-function(SimPath,locationsKernel,locationsData,SNR){
    Reg<-TRUE
    ErrorResults <- vector("list", length(locationsKernel))

    for(Dlength in 1:length(locationsKernel)){

El1ErrorRaw<-array(0,c(length(R.all),length(lambda.all)))

        inname<-paste0(SimPath, locationsKernel[Dlength],"/")
        outname<-paste0(SimPath, locationsData[Dlength],"/")
        
        if(file.exists(paste0(outname,outname1))==FALSE) dir.create(paste0(outname,outname1))

        SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
        seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
        LFP<-as.matrix(fread(paste(outname,"/myLFP",sep="")))
        
if (SNR!=0){
      LFPvariance<-var(c(LFP))
      Noisevariance<-LFPvariance/SNR
      NoiseLFP<-rnorm(length(c(LFP)), mean = 0, sd = sqrt(Noisevariance))
      LFP<-LFP+NoiseLFP
    }
LFPOriginal<-LFP
        
        ElCoords<-as.matrix(read.table(paste(inname,"elcoord_x_y_z",sep="")))
        ElCoords<-matrix(ElCoords,ncol=3)
        Xel<-unique(ElCoords[,2])
        XelNb<-length(Xel)
        Yel<-unique(ElCoords[,3])
        
       
for(lamb in 1:length(lambda.all)){
  for(r in 1:length(R.all)){
  
  LFP<-LFPOriginal
    lambda<-lambda.all[lamb]
    R<-R.all[r]

  source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))

 CurrNew<-KernSmoothDistance(30,inname, outname, outname1,"30")

 el2ignorename<-paste(inname,outname1,"/El2Ignore_M",M,"_R",R,sep="")
    if(file.exists(el2ignorename)){
      El2Ignore<-c(as.matrix(read.table(el2ignorename)))
      LFP<-LFP[-El2Ignore,] }
      


Kmatr<-as.matrix(read.table(paste(inname,outname1,"/K_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
    Ktildematr<-as.matrix(read.table(paste(inname,outname1,"/Ktilda_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
      
      Tmatr<-(4*pi*sigma)*t(Ktildematr)%*%ginv(Kmatr)
      


skCSD.all.part<-Tmatr%*%LFP
    write.table(skCSD.all.part, paste(outname,outname1,"/ksCSD_Matr",M,"_R",R,"lambda",
                                                   lambda.all[lamb], "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
     


 skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))

 SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

 for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
    }
 

if(CellType=="BSOLD") {skCSD.all<-skCSD.all[c(2,1,3:52),]}
El1ErrorRaw[r,lamb]<-L1Error(skCSD.all,CurrNew)[[1]]
 
      
#      cvErrorOut<-cross.valid(LFP,Tmatr,1:dim(Tmatr)[1])
#ErrorCV[r,lamb]<-cvErrorOut[[1]]
cat("Numb to rteplace:", length(R.all)*(lamb-1)+ r)
#CVConfig[[length(R.all)*(lamb-1)+ r]]<-as.matrix(cvErrorOut[[4]])
#write.table(as.matrix(cvErrorOut[[4]]),paste(outname,outname1,"/CVContribution_M",M,"_R",R,"lambda",lambda.all[lamb], 
#                                             "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
#ErrorContribution
}} #lambda and R

write.table(El1ErrorRaw, paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
        
}#Dlength
#return(list(ErrorResults, CVConfig))
}



korte<-JustL1Reg(SimPath,locationsKernel,locationsData,SNR)


read.table(paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,sep=""))

#####################JUST L1 ERROR within TImeframe


JustL1RegTime<-function(SimPath,locationsKernel,locationsData,SNR, Timeframe){
  Reg<-TRUE
  ErrorResults <- vector("list", length(locationsKernel))
  
  for(Dlength in 1:length(locationsKernel)){
    
    El1ErrorRaw<-array(0,c(length(R.all),length(lambda.all)))
    
    inname<-paste0(SimPath, locationsKernel[Dlength],"/")
    outname<-paste0(SimPath, locationsData[Dlength],"/")
    
    if(file.exists(paste0(outname,outname1))==FALSE) dir.create(paste0(outname,outname1))
    
    SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
    LFP<-as.matrix(fread(paste(outname,"/myLFP",sep="")))
    
    if (SNR!=0){
      LFPvariance<-var(c(LFP))
      Noisevariance<-LFPvariance/SNR
      NoiseLFP<-rnorm(length(c(LFP)), mean = 0, sd = sqrt(Noisevariance))
      LFP<-LFP+NoiseLFP
    }
    LFPOriginal<-LFP
    
    ElCoords<-as.matrix(read.table(paste(inname,"elcoord_x_y_z",sep="")))
    ElCoords<-matrix(ElCoords,ncol=3)
    Xel<-unique(ElCoords[,2])
    XelNb<-length(Xel)
    Yel<-unique(ElCoords[,3])
    
    
    for(lamb in 1:length(lambda.all)){
      for(r in 1:length(R.all)){
        
        LFP<-LFPOriginal
        lambda<-lambda.all[lamb]
        R<-R.all[r]
        
        source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))
        
        CurrNew<-KernSmoothDistance(30,inname, outname, outname1,"30")
        
        el2ignorename<-paste(inname,outname1,"/El2Ignore_M",M,"_R",R,sep="")
        if(file.exists(el2ignorename)){
          El2Ignore<-c(as.matrix(read.table(el2ignorename)))
          LFP<-LFP[-El2Ignore,] }
        
        
        
        Kmatr<-as.matrix(read.table(paste(inname,outname1,"/K_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
        Ktildematr<-as.matrix(read.table(paste(inname,outname1,"/Ktilda_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
        
        Tmatr<-(4*pi*0.5)*t(Ktildematr)%*%ginv(Kmatr)
        
        
        
        skCSD.all.part<-Tmatr%*%LFP
     #   write.table(skCSD.all.part, paste(outname,outname1,"/ksCSD_Matr",M,"_R",R,"lambda",
      #                                    lambda.all[lamb], "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
        
        
        
        skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
        
        SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
        seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
        
        for(i in 1: seg.nb){
          sameplace<-numeric() 
          sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
          #Plot the current related ti the different segments
          #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
          if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
          if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
        }
        
        
        if(CellType=="BSOLD") {skCSD.all<-skCSD.all[c(2,1,3:52),]}
        El1ErrorRaw[r,lamb]<-L1Error(skCSD.all[,Timeframe],CurrNew[,Timeframe])[[1]]
        
        
        #      cvErrorOut<-cross.valid(LFP,Tmatr,1:dim(Tmatr)[1])
        #ErrorCV[r,lamb]<-cvErrorOut[[1]]
        cat("Numb to rteplace:", length(R.all)*(lamb-1)+ r)
        #CVConfig[[length(R.all)*(lamb-1)+ r]]<-as.matrix(cvErrorOut[[4]])
        #write.table(as.matrix(cvErrorOut[[4]]),paste(outname,outname1,"/CVContribution_M",M,"_R",R,"lambda",lambda.all[lamb], 
        #                                             "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
        #ErrorContribution
      }} #lambda and R
    
    write.table(El1ErrorRaw, paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,"_",Timeframe[1],"_",Timeframe[2],sep=""),row.names = FALSE,col.names = FALSE)
    
  }#Dlength
  #return(list(ErrorResults, CVConfig))
}



korte<-JustL1RegTime(SimPath,locationsKernel,locationsData,SNR, c(1,1000))



########################## Calculate L1 just for the best parameters from previous reconstruction for other times




JustL1RegTimeWBestParams<-function(SimPath,locationsKernel,locationsData,SNR, Timeframe,TimeframeEval){
  Reg<-TRUE
  ErrorResults <- vector("list", length(locationsKernel))
  El1ErrorRawAllGangAll<-numeric()
  for(Dlength in 1:length(locationsKernel)){
    
    inname<-paste0(SimPath, locationsKernel[Dlength],"/")
    outname<-paste0(SimPath, locationsData[Dlength],"/")
    
    if(file.exists(paste0(outname,outname1))==FALSE) dir.create(paste0(outname,outname1))
    
    SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
    LFP<-as.matrix(fread(paste(outname,"/myLFP",sep="")))
    
    if (SNR!=0){
      LFPvariance<-var(c(LFP))
      Noisevariance<-LFPvariance/SNR
      NoiseLFP<-rnorm(length(c(LFP)), mean = 0, sd = sqrt(Noisevariance))
      LFP<-LFP+NoiseLFP
    }
    LFPOriginal<-LFP
    
    ElCoords<-as.matrix(read.table(paste(inname,"elcoord_x_y_z",sep="")))
    ElCoords<-matrix(ElCoords,ncol=3)
    Xel<-unique(ElCoords[,2])
    XelNb<-length(Xel)
    Yel<-unique(ElCoords[,3])
    
    El1ErrorRaw<-read.table(paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,"_",Timeframe[1],"_",Timeframe[2],sep=""))
    
    MinIndex<-which(El1ErrorRaw==min(El1ErrorRaw),arr.ind=TRUE)
    
    
      r<-MinIndex[1]
      lamb<-MinIndex[2]
        LFP<-LFPOriginal
        lambda<-lambda.all[lamb]
        R<-R.all[r]
        
        source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))
        
        CurrNew<-KernSmoothDistance(30,inname, outname, outname1,"30")
        
        el2ignorename<-paste(inname,outname1,"/El2Ignore_M",M,"_R",R,sep="")
        if(file.exists(el2ignorename)){
          El2Ignore<-c(as.matrix(read.table(el2ignorename)))
          LFP<-LFP[-El2Ignore,] }
        
        
        
        Kmatr<-as.matrix(read.table(paste(inname,outname1,"/K_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
        Ktildematr<-as.matrix(read.table(paste(inname,outname1,"/Ktilda_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
        
        Tmatr<-(4*pi*0.5)*t(Ktildematr)%*%ginv(Kmatr)
        
        
        
        skCSD.all.part<-Tmatr%*%LFP
       # write.table(skCSD.all.part, paste(outname,outname1,"/ksCSD_Matr",M,"_R",R,"lambda",
         #                                 lambda.all[lamb], "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
        
        
        
        skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
        
        SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
        seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
        
        for(i in 1: seg.nb){
          sameplace<-numeric() 
          sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
          #Plot the current related ti the different segments
          #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
          if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
          if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
        }
        
        
        if(CellType=="BSOLD") {skCSD.all<-skCSD.all[c(2,1,3:52),]}
        El1ErrorRawAllGang<-numeric()
        El1ErrorRawAllGang<-L1Error(skCSD.all[,TimeframeEval],CurrNew[,TimeframeEval])[[1]]
        
        El1ErrorRawAllGangAll<-c(El1ErrorRawAllGangAll,El1ErrorRawAllGang)
        #      cvErrorOut<-cross.valid(LFP,Tmatr,1:dim(Tmatr)[1])
        #ErrorCV[r,lamb]<-cvErrorOut[[1]]
        cat("Numb to rteplace:", length(R.all)*(lamb-1)+ r)
        #CVConfig[[length(R.all)*(lamb-1)+ r]]<-as.matrix(cvErrorOut[[4]])
        #write.table(as.matrix(cvErrorOut[[4]]),paste(outname,outname1,"/CVContribution_M",M,"_R",R,"lambda",lambda.all[lamb], 
        #                                             "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
        #ErrorContribution
      #lambda and R
    
    write.table(El1ErrorRawAllGang, paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,"_",TimeframeEval[1],"_",TimeframeEval[2],sep=""),row.names = FALSE,col.names = FALSE)
    write.table(c(r,lamb), paste(outname,outname1,"/L1Compare_BestParams",sep=""),row.names = FALSE,col.names = FALSE)
    
  }#Dlength
  #return(list(ErrorResults, CVConfig))
}



JustL1RegTimeWBestParams(SimPath,locationsKernel,locationsData,SNR, c(1,1000),c(1001,6800))





###############################################################################


######################for Soike triggered data



JustL1RegSPTrig<-function(SimPath,locationsKernel,locationsData,SNR){
    Reg<-TRUE
    ErrorResults <- vector("list", length(locationsKernel))

    for(Dlength in 1:length(locationsKernel)){

El1ErrorRaw<-array(0,c(length(R.all),length(lambda.all)))

        inname<-paste0(SimPath, locationsKernel[Dlength],"/")
        outname<-paste0(SimPath, locationsData[Dlength],"/")
        
        if(file.exists(paste0(outname,outname1))==FALSE) dir.create(paste0(outname,outname1))

        SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
        seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
        LFP<-as.matrix(fread(paste(outname,"/myLFP",sep="")))
        
if (SNR!=0){
      LFPvariance<-var(c(LFP))
      Noisevariance<-LFPvariance/SNR
      NoiseLFP<-rnorm(length(c(LFP)), mean = 0, sd = sqrt(Noisevariance))
      LFP<-LFP+NoiseLFP
    }
LFPOriginal<-LFP
        
        ElCoords<-as.matrix(read.table(paste(outname,"elcoord_x_y_z",sep="")))
        ElCoords<-matrix(ElCoords,ncol=3)
        Xel<-unique(ElCoords[,2])
        XelNb<-length(Xel)
        Yel<-unique(ElCoords[,3])
        
       
for(lamb in 1:length(lambda.all)){
  for(r in 1:length(R.all)){
  
  LFP<-LFPOriginal
    lambda<-lambda.all[lamb]
    R<-R.all[r]

  source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))

 CurrNew<-KernSmoothDistanceSpTrig(30,inname, outname, outname1,"30")

 el2ignorename<-paste(inname,outname1,"/El2Ignore_M",M,"_R",R,sep="")
    if(file.exists(el2ignorename)){
      El2Ignore<-c(as.matrix(read.table(el2ignorename)))
      LFP<-LFP[-El2Ignore,] }
      


Kmatr<-as.matrix(read.table(paste(inname,outname1,"/K_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
    Ktildematr<-as.matrix(read.table(paste(inname,outname1,"/Ktilda_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
      
      Tmatr<-(4*pi*0.5)*t(Ktildematr)%*%ginv(Kmatr)
      


skCSD.all.part<-Tmatr%*%LFP
    write.table(skCSD.all.part, paste(outname,outname1,"/ksCSD_Matr",M,"_R",R,"lambda",
                                                   lambda.all[lamb], "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
     


 skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))

 SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

 for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
    }
 

if(CellType=="BSOLD") {skCSD.all<-skCSD.all[c(2,1,3:52),]}
El1ErrorRaw[r,lamb]<-L1Error(skCSD.all,CurrNew)[[1]]
 
      
#      cvErrorOut<-cross.valid(LFP,Tmatr,1:dim(Tmatr)[1])
#ErrorCV[r,lamb]<-cvErrorOut[[1]]
cat("Numb to rteplace:", length(R.all)*(lamb-1)+ r)
#CVConfig[[length(R.all)*(lamb-1)+ r]]<-as.matrix(cvErrorOut[[4]])
#write.table(as.matrix(cvErrorOut[[4]]),paste(outname,outname1,"/CVContribution_M",M,"_R",R,"lambda",lambda.all[lamb], 
#                                             "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
#ErrorContribution
}} #lambda and R

write.table(El1ErrorRaw, paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
        
}#Dlength
#return(list(ErrorResults, CVConfig))
}


korte<-JustL1RegSPTrig(SimPath,locationsKernel,locationsData,SNR)







