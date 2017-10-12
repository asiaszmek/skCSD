
#THIS FUNCTION PLOTS THE ERRORS IN CASE OF DIFFERENT VVERSIONS OF THE RECONSTRUCTED CURRENTS
MatPlotError<-function(whichError2Plot, titleofGraph,R.all,M.all){
  par(mar=c(4,3,4,6))
  matplot(whichError2Plot,t="l", main=titleofGraph,xaxt="n",xlab="R", log="y",ylim=c(min(whichError2Plot),median(whichError2Plot)*2))
  axis(side=1, labels=rep(R.all,4),at=c(1:(length(R.all)*4)))
  mtext("Original",side=3, at=length(R.all)/2,line=-1)
  mtext("fix Smooth",side=3, at=3/2*length(R.all),line=-1)
  mtext("R Smooth",side=3, at=5/2*length(R.all),line=-1)
  mtext("R2 Smooth",side=3, at=7/2*length(R.all),line=-1)
  abline(v=c((1:3)*length(R.all)+0.5))
  par(xpd=TRUE)
  legend(4*length(R.all)+2,median(whichError2Plot),c(paste(M.all)),pch=20,col=1:6,title="M")
  par(xpd=FALSE)
}


#################################################
# FUNCTION TO CALCULATE THE CURRENTS, eRRORS


CalcErrorAndPlot<-function(SimPath,locationsKernel,locationsData,SNR, lammbda){
#SimPath<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/simulation/"

for(Dlength in 1:length(locationsKernel)){

inname<-paste0(SimPath, locationsKernel[Dlength],"/")
outname<-paste0(SimPath, locationsData[Dlength])
if(file.exists(paste(inname,outname1,"/params",sep=""))==FALSE) next
if(file.exists(paste0(outname,outname1))==FALSE) dir.create(paste0(outname,outname1))

params<-as.matrix(read.table(paste(inname,outname1,"/params",sep="")))
R.all<-unique(params[1,])#[3]
M.all<-unique(params[2,])#[5]
#lambda<-0#10
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

funaramvonal<-function(x) x/seg.length
memb.currents.line<-array(0, c(dim(membcurr)))
memb.currents.line<-apply(membcurr,2,funaramvonal) 
morpho<-as.matrix(read.table(paste0(inname,'segcoordinates.txt')))
SmoothPar<-as.matrix(read.table(paste0(inname,outname1,"/SmoothinKernelWidth")))

EL1Error<-array(0,c(length(R.all)*4,length(M.all)))
EL1ErrorAbs<-array(0,c(length(R.all)*4,length(M.all)))
E1Normalized<-array(0,c(length(R.all)*4,length(M.all)))
E1NormalizedAbs<-array(0,c(length(R.all)*4,length(M.all)))

EL2Error<-array(0,c(length(R.all)*4,length(M.all)))
EL2ErrorAbs<-array(0,c(length(R.all)*4,length(M.all)))
ENormalized<-array(0,c(length(R.all)*4,length(M.all)))
ENormalizedAbs<-array(0,c(length(R.all)*4,length(M.all)))
E2SegNorm<-array(0,c(length(R.all)*4,length(M.all)))
CrossVal<-array(0,c(length(R.all)*4,length(M.all)))


EL1ErrorSoma<-array(0,c(length(R.all)*4,length(M.all)))
EL1ErrorAbsSoma<-array(0,c(length(R.all)*4,length(M.all)))
E1NormalizedSoma<-array(0,c(length(R.all)*4,length(M.all)))
E1NormalizedAbsSoma<-array(0,c(length(R.all)*4,length(M.all)))

EL2ErrorSoma<-array(0,c(length(R.all)*4,length(M.all)))
EL2ErrorAbsSoma<-array(0,c(length(R.all)*4,length(M.all)))
ENormalizedSoma<-array(0,c(length(R.all)*4,length(M.all)))
ENormalizedAbsSoma<-array(0,c(length(R.all)*4,length(M.all)))
E2SegNormSoma<-array(0,c(length(R.all)*4,length(M.all)))
CrossValSoma<-array(0,c(length(R.all)*4,length(M.all)))

#Reading in the LFP and adding noise
    LFP<-as.matrix(read.table(paste(outname,"/myLFP",sep="")))
if (SNR!=0){
      LFPvariance<-var(c(LFP))
      Noisevariance<-LFPvariance/SNR
      NoiseLFP<-rnorm(length(c(LFP)), mean = 0, sd = sqrt(Noisevariance))
      LFP<-LFP+NoiseLFP
    }
    write.table(LFP,paste0(outname,outname1,"/LFP_SNR",SNR,"_lambda",lambda))

LFPOriginal<-LFP




CVError<-array(0,c(length(R.all),length(M.all)))
Kappaval<-numeric()
for(m in 1:length(M.all)){
  for(r in 1:length(R.all)){
    LFP<-LFPOriginal
    M<-M.all[m]
    R<-R.all[r]
    
    

    if(file.exists(paste(inname,outname1,"/K_M",M,"_R",R,sep=""))==FALSE) {
	Kappaval<-c(Kappaval, NA)
	
} else {
    #Tmatr<-as.matrix(read.table(paste(outname,outname1,"/TransferMatr_M",M,"_R",R,sep="")))
    Kmatr<-as.matrix(read.table(paste(inname,outname1,"/K_M",M,"_R",R,sep="")))
    Ktildematr<-as.matrix(read.table(paste(inname,outname1,"/Ktilda_M",M,"_R",R,sep="")))
    Kappaval<-c(Kappaval, kappa(Kmatr))
	
    
	
    skCSD.all<-array(0,c(dim(membcurr)))
    
    el2ignorename<-paste(inname,outname1,"/El2Ignore_M",M,"_R",R,sep="")
    if(file.exists(el2ignorename)){
      El2Ignore<-c(as.matrix(read.table(el2ignorename)))
      LFP<-LFP[-El2Ignore,]
    }
    
    SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

    RegMatr<-array(0,c(dim(Kmatr)))
    diag(RegMatr)<-lambda
    #C.calc Tikhonov regularization with lambda 
    C.calc<-try((4*pi*sigma)*t(Ktildematr)%*%ginv(Kmatr+RegMatr)%*%LFP ,silent=TRUE)
    Tmatr<-try(1/(4*pi*sigma)*t(Ktildematr)%*%ginv(Kmatr+RegMatr),silent=TRUE)
    skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
    for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<-C.calc[sameplace,]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums(C.calc[c(sameplace),]),nrow=1)
    }
    #Plot C.calc
    #image(t(C.calc),xlab="Time",ylab="seg.nb")
    #Plot the current at the same phisical location -- sameplace
    whIndex<-60
    WhichToplot<-as.integer(unlist(strsplit(SameplaceAll," ")[whIndex]))
    #matplot(cbind(t(C.calc[c(WhichToplot),])),t='l', xlab="Time steps", ylab="ksCSD", main=paste("Values of currents at segment",whIndex))
    #lines(colSums(C.calc[c(WhichToplot),]),col="orange")
    #lines(memb.currents.smoothed[(whIndex),],col="purple")
    #legend("topright",c(as.character(WhichToplot)),col=c(1:length(WhichToplot)), pch=20)
    
    
    #######################################################
    #Contribution of basis function at a point
    Beta<-try((ginv(Kmatr+RegMatr)%*%LFP),silent=TRUE)#[-El2Ignore,]
    if(file.exists(el2ignorename)){Bmatr<-as.matrix(read.table(paste(inname,outname1,"/B_M",M,"_R",R,sep="")))[,-El2Ignore]} else Bmatr<-as.matrix(read.table(paste(inname,outname1,"/B_M",M,"_R",R,sep="")))
    Btildematr<-as.matrix(read.table(paste(inname,outname1,"/Btilda_M",M,"_R",R,sep="")))
    
    j.basis<-1   #1:M
    k.location<-3 # 1:dim(C.calc)[1] 
    timeInst<-which(LFP==min(LFP),arr.ind=TRUE)[2]
    #timeInst<-40
    CBasisCont<-array(0,c(M,dim(C.calc)[1] ))
    for(j.basis in 1:M){
      CBasisCont[j.basis,]<-Beta[,timeInst]%*%Bmatr[j.basis,]%*%Btildematr[j.basis,] #jth source contribution
    }

	write.table(CBasisCont,paste0(outname,outname1,"/CbasisContribute_M",M,"_R",R))
    whichCont<-44
    #plot(CBasisCont[,whichCont], main=paste("contribution of basis function at point",whichCont),xlab="M")
    whichBasis<-1
    #plot(CBasisCont[whichBasis,], main=paste("contribution of one basis",whichBasis,"at diff locations"),xlab="M")
    #plot(C.calc[,timeInst])
    #plot(skCSD.all[,timeInst])


    #dir.create(paste0(outname,outname1))
    write.table(skCSD.all,paste0(outname,outname1,"/skCSDall_M",M,"_R",R,"_SNR", SNR, "_lambda",lambda))
    #Reading in tge currents
    
    ########################Calculate Currents
    ############################################
    #Calculate Errors
    
SmoothPar<-as.matrix(read.table(paste0(inname,outname1,"/SmoothinKernelWidth")))
BandWidth1<-2*sqrt(SmoothPar[2]/SmoothPar[3])

    memb.currents.smoothedR<-KernSmoothDistance(sqrt(SmoothElWidth),inname, outname, outname1,"Smoothed_R1")#(2*SmoothElWidth, morpho, memb.currents.line)
    memb.currents.smoothedsqrR<-KernSmoothDistance(SmoothElWidth,inname, outname, outname1,"Smoothed_R2") #(2*SmoothElWidth, morpho, memb.currents.line,,"Smoothed_R2")
    memb.currents.smoothedfixR<-KernSmoothDistance(BandWidth1,inname, outname, outname1,"Smoothed_fixR")
    #memb.currents.smoothedfixR<-SmoothTheCurrents(SmoothPar[2]/SmoothPar[3], morpho, memb.currents.line)$memb.currents.smoothed
    write.table(memb.currents.smoothedfixR,paste0(outname,outname1,"/FixSmoothed_M",M,"_R",R,"_SNR", SNR, "_lambda",lambda))
    
   
    #Running crossvalidation for M,R and lambda??? ehhh....
    #cv.error<-cross.valid(LFP,Tmatr)
    #which segments should taken into account at compariso?
	 ToCompareSoma<-1:dim(skCSD.all)[1]
    	 ToCompare<-5:dim(skCSD.all)[1]
    Which2CompareVar<-list( memb.currents.line, memb.currents.smoothedfixR, memb.currents.smoothedR,memb.currents.smoothedsqrR)

    for(errvariants in 1: length(Which2CompareVar)){
      Which2Compare<-Which2CompareVar[[errvariants]]
      EL2Error[(errvariants-1)*length(R.all)+r,m]<-L2Error(skCSD.all[ ToCompare,],Which2Compare[ ToCompare,]) [[1]]
      EL2ErrorAbs[(errvariants-1)*length(R.all)+r,m]<-L2ErrorAbs(skCSD.all[ToCompare,],Which2Compare[ToCompare,]) [[1]]
      ENormalized[(errvariants-1)*length(R.all)+r,m]<-ErrorNormalized(skCSD.all[ToCompare,], Which2Compare[ToCompare,]) [[1]]
      ENormalizedAbs[(errvariants-1)*length(R.all)+r,m]<-ErrorNormalizedAbs(skCSD.all[ToCompare,], Which2Compare[ToCompare,]) [[1]]
       E2SegNorm[(errvariants-1)*length(R.all)+r,m]<-L2ErrorSegNorm(skCSD.all[ToCompare,], Which2Compare[ToCompare,]) [[1]]
      CrossVal[(errvariants-1)*length(R.all)+r,m]<-try(cross.valid(LFP,Tmatr,ToCompare)[[1]],silent=TRUE)

    }

for(errvariants in 1: length(Which2CompareVar)){
      Which2Compare<-Which2CompareVar[[errvariants]]
      EL2ErrorSoma[(errvariants-1)*length(R.all)+r,m]<-L2Error(skCSD.all[ ToCompareSoma,],Which2Compare[ ToCompareSoma,]) [[1]]
      EL2ErrorAbsSoma[(errvariants-1)*length(R.all)+r,m]<-L2ErrorAbs(skCSD.all[ToCompareSoma,],Which2Compare[ToCompareSoma,]) [[1]]
      ENormalizedSoma[(errvariants-1)*length(R.all)+r,m]<-ErrorNormalized(skCSD.all[ToCompareSoma,], Which2Compare[ToCompareSoma,]) [[1]]
      ENormalizedAbsSoma[(errvariants-1)*length(R.all)+r,m]<-ErrorNormalizedAbs(skCSD.all[ToCompareSoma,], Which2Compare[ToCompareSoma,]) [[1]]
      CrossValSoma[(errvariants-1)*length(R.all)+r,m]<-try(cross.valid(LFP,Tmatr,ToCompareSoma)[[1]],silent=TRUE)
      
    }
    
    #Plotting the kernel functions
    #memb.currents.line
   
   
    for(errvariants in 1: length(Which2CompareVar)){
      Which2Compare<-Which2CompareVar[[errvariants]]
      EL1Error[(errvariants-1)*length(R.all)+r,m]<-L1Error(skCSD.all[ ToCompare,],Which2Compare[ ToCompare,]) [[1]]
      EL1ErrorAbs[(errvariants-1)*length(R.all)+r,m]<-L1ErrorAbs(skCSD.all[ToCompare,],Which2Compare[ToCompare,]) [[1]]
      E1Normalized[(errvariants-1)*length(R.all)+r,m]<-ErrorNormalizedL1(skCSD.all[ToCompare,], Which2Compare[ToCompare,]) [[1]]
      E1NormalizedAbs[(errvariants-1)*length(R.all)+r,m]< ErrorNormalizedAbsL1(skCSD.all[ToCompare,], Which2Compare[ToCompare,]) [[1]]
      
    }

for(errvariants in 1: length(Which2CompareVar)){
      Which2Compare<-Which2CompareVar[[errvariants]]
      EL1ErrorSoma[(errvariants-1)*length(R.all)+r,m]<-L1Error(skCSD.all[ ToCompareSoma,],Which2Compare[ ToCompareSoma,]) [[1]]
      EL1ErrorAbsSoma[(errvariants-1)*length(R.all)+r,m]<-L1ErrorAbs(skCSD.all[ToCompareSoma,],Which2Compare[ToCompareSoma,]) [[1]]
      E1NormalizedSoma[(errvariants-1)*length(R.all)+r,m]<-ErrorNormalizedL1(skCSD.all[ToCompareSoma,], Which2Compare[ToCompareSoma,]) [[1]]
      E1NormalizedAbsSoma[(errvariants-1)*length(R.all)+r,m]< ErrorNormalizedAbsL1(skCSD.all[ToCompareSoma,], Which2Compare[ToCompareSoma,]) [[1]]
      E2SegNormSoma[(errvariants-1)*length(R.all)+r,m]<-L2ErrorSegNorm(skCSD.all[ToCompareSoma,], Which2Compare[ToCompareSoma,]) [[1]]
      
    }
    CVError[r,m]<-cross.valid(LFP, Tmatr,ToCompareSoma)[[1]]
    
    
  }
}
}
#Writing to file some error values
write.table(EL2Error, paste(outname,outname1,"/EL2ErrorWOsoma","_SNR", SNR, "_lambda",lambda,sep=""))
write.table(EL2ErrorSoma, paste(outname,outname1,"/EL2ErrorSoma","_SNR", SNR, "_lambda",lambda,sep=""))

write.table(E2SegNorm, paste(outname,outname1,"/E2SegNormWOsoma","_SNR", SNR, "_lambda",lambda,sep=""))
write.table(E2SegNormSoma, paste(outname,outname1,"/E2SegNormrSoma","_SNR", SNR, "_lambda",lambda,sep=""))
write.table(CrossVal, paste(outname,outname1,"/ECrossVal","_SNR", SNR, "_lambda",lambda,sep=""))

write.table(CVError, paste(outname,outname1,"/CVError","_SNR", SNR, "_lambda",lambda,sep=""))

CVErrorBest<-c( R.all[which(CVError==min(CVError),arr.ind=TRUE)[1]], 
                        M.all[which(CVError==min(CVError),arr.ind=TRUE)[2]],min(CVError))

write.table(CVErrorBest, paste(outname,outname1,"/CVErrorBest","_SNR", SNR, "_lambda",lambda,sep=""))


#Writing to file the kappa values


KappaMatr<-matrix(Kappaval, nrow=length(R.all))
write.table(KappaMatr, paste(outname,outname1,"/kappaK","_SNR", SNR, "_lambda",lambda,sep=""))
png( paste0(inname,outname1,"/KappaPlot_SNR", SNR, "_sigma",sigma,".png"))
matplot(t(KappaMatr),log="y",t="l", xlab="R", ylab="Kappa", 
main="Kappa value of K",sub=paste("Number of electrodes:", ElNumb[Dlength-(Dlength%/%(length(ElNumb)+1))*length(ElNumb)] ))
legend("bottomright", c(paste(M.all)),col=1:6, pch=20, title="M" )
dev.off()

#E1=memb.currents.line,E2= memb.currents.smoothedfixR, E3=memb.currents.smoothedR,E4=memb.currents.smoothedsqrR



 png(paste0(SumPlotPlace,locationsData[Dlength],"_SNR", SNR, "_lambda",lambda,"withoutSoma.png"),width=800, height=1000)

par(mfrow=c(5,2))
MatPlotError(EL2Error,"L2 Error",R.all,M.all)
MatPlotError(EL2ErrorAbs,"EL2ErrorAbs",R.all,M.all)
MatPlotError(E2SegNorm,"Normalized On Segments L2",R.all,M.all)
MatPlotError(ENormalized,"Normalized L2 Error",R.all,M.all)
MatPlotError(ENormalizedAbs,"Normalized Absolute Error",R.all,M.all)


MatPlotError(EL1Error,"L1 Error",R.all,M.all)
MatPlotError(EL1ErrorAbs,"EL1ErrorAbs",R.all,M.all)
MatPlotError(E1Normalized,"Normalized L1 Error",R.all,M.all)
MatPlotError(E1NormalizedAbs,"Normalized Absolute L1 Error",R.all,M.all)
MatPlotError(CrossVal,"CrossValidated",R.all,M.all)

dev.off()


png(paste0(SumPlotPlace,locationsData[Dlength],"_SNR", SNR, "_lambda",lambda,"wSoma.png"),width=800, height=1000)

par(mfrow=c(5,2))
MatPlotError(EL2ErrorSoma,"L2 Error",R.all,M.all)
MatPlotError(EL2ErrorAbsSoma,"EL2ErrorAbs",R.all,M.all)
MatPlotError(E2SegNormSoma,"Normalized On Segments L2",R.all,M.all)
MatPlotError(ENormalizedSoma,"Normalized L2 Error",R.all,M.all)
MatPlotError(ENormalizedAbsSoma,"Normalized Absolute Error",R.all,M.all)


MatPlotError(EL1ErrorSoma,"L1 Error",R.all,M.all)
MatPlotError(EL1ErrorAbsSoma,"EL1ErrorAbs",R.all,M.all)
MatPlotError(E1NormalizedSoma,"Normalized L1 Error",R.all,M.all)
MatPlotError(E1NormalizedAbsSoma,"Normalized Absolute L1 Error",R.all,M.all)
MatPlotError(CrossValSoma,"CrossValidated",R.all,M.all)
dev.off()

}

} # end funtion



########################################################3
######################################################
#FOR PLOTTING THE RESULTS OF THE PREVIOUS FUNVTION... fOR R AN M THIS PLOTS THE RECONSTRUCTED skCSD

ImageAllParCurrReg<-function(SimPath,locationsKernel,locationsData,SNR,  R2Plot,lambda2Plot){

for(Dlength in 1:length(locationsKernel)){

inname<-paste0(SimPath, locationsKernel[Dlength],"/")
outname<-paste0(SimPath, locationsData[Dlength])
#if(file.exists(paste(inname,outname1,"/params",sep=""))==FALSE) next
#params<-as.matrix(read.table(paste(inname,outname1,"/params",sep="")))
OrigCurr<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/membcurr",sep="")))
 SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

seg.length<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/","seglength",sep="")))
funaramvonal<-function(x) x/seg.length
memb.currents.line1<-array(0, c(dim(OrigCurr)))
memb.currents.line1<-apply(OrigCurr,2,funaramvonal) 
memb.currents.line1[which(memb.currents.line1<(0.5*min(memb.currents.line1)))]<-0.5*min(memb.currents.line1)
#memb.currents.line1[which(memb.currents.line1<(0.7*max(memb.currents.line1)))]<-0.5*max(memb.currents.line1)
#OrigCurrS<-as.matrix(read.table(paste(outname,outname1,"/membcurr_smoothed",sep="")))
sCSD<-as.matrix(read.table(paste0(outname,"/sCSDginv")))
times2Plot<-as.matrix(read.table(paste0(inname,"time")))

#Select just few Rs and Ms for plotting
R.all<-R2Plot#seq(5,130,25)
M<-512



sCSDCoord<-as.matrix(read.table(paste0(outname,"/sCSDCoord")))
if( outname=="/kernelOut_poli") {png(paste0(SumPlotPlace,locationsData[Dlength],"skCSD_all_SNR", SNR, "_lambda",lambda,".png"),width=1000, height=900, pointsize=22)
} else { png(paste0(SumPlotPlace,locationsData[Dlength],"skCSDall_M",M,"SNR", SNR,".png"),width=1000, height=900, pointsize=22)}
plot.new()
par(mfrow=c(length(lambda.all)+2,length(R.all)),oma=c(4,2,4,4))#,mar=c(1,1,1,1))
#par(mfrow=c(length(lambda.all)+2,length(R.all)),oma=c(2,3,4,2))#,mar=c(1,1,1,1))

#layout(rbind(c(1,2,3,4,5)),widths=c(0.5,0.2,1,3,0.2),respect=FALSE) 
#mtext(c("oma[SOUTH<-1]", "oma[WEST<-2]", "oma[NORTH<-3]", "oma[EAST<-4]"),     c(SOUTH<-1, WEST<-2, NORTH<-3, EAST<-4),      line=0.4, cex=1.2, col="green", outer=TRUE)
title(main="skCSD Reconstruction", outer=TRUE, cex=3)
mtext(c("Time (ms)", "Segment #"), c(1,4), outer=TRUE,line=2.5)


 par(mfg=c(1,1))
col2<-ColoursDori(memb.currents.line1)[[1]]
ExtVal<-ColoursDori(memb.currents.line1)[[2]]
par(mar=c(0.1,0.1,0.1,0.1))
SegNb<-dim(memb.currents.line1)[1]
image(times2Plot,1:SegNb,t(memb.currents.line1),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", yaxt="n")
mtext("Ground Truth",3)
axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
axis(4,round(seq(1,SegNb-10,,4)),las=2)
mtext("Segment #",3,line=5)
col2<-ColoursDori(memb.currents.line1)[[1]]
ExtVal<-ColoursDori(memb.currents.line1)[[2]]

par(mfg=c(1,2),mar=c(0.1,0.1,0.1,10)) 
image.plot(times2Plot,1:SegNb,t(memb.currents.line1),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", yaxt="n",axis.args=list( at=c(-ExtVal, ExtVal), labels=c("Sink", "Source")),
 legend.only=TRUE)

if(R2Plot<4) par(mar=c(0,8,0,8)) 
if(R2Plot>3) par(mar=c(0,4,0,4)) 


if(CellType=="BS"){
par(mfg=c(1,3))
par(mar=c(0.1,0.1,0.1,0.1))
image(times2Plot,sCSDCoord[,3],t(sCSD),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n",ylab="z (um)")
mtext("sCSD",3)
axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
#par(mar=c(0,3,0,5)) 
}


#image(c(times2Plot),1:dim(memb.currents.line1)[1],array(0,c(length(times2Plot),dim(memb.currents.line1)[1])),xlab="Time",ylab="Segment #",col="WHITE")
# image(t(OrigCurrS),col=rainbow(500),main="Smoothed Ground Truth")
  par(mfg=c(3,1))

    skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
for(m in 1:length(lambda.all)){
  for(r in 1:length(R.all)){
    
    Lambda<-lambda.all[m]
    R<-R.all[r]
    #Reading in the LFP and adding noise
    #Ktilda_M512_R55lambda1e-04
currName<-paste0(outname,outname1,"/skCSDall_M",M,"_R",R,"lambda",Lambda)
    if(file.exists(currName)) {
      skCSD.all.part<-as.matrix(read.table(currName)) 
      

 for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
    }      
      
      
      
      
                               }else next
col2<-ColoursDori(skCSD.all)[[1]]
ExtVal<-ColoursDori(skCSD.all)[[2]]
par(mar=c(0.1,0.1,0.1,0.1))
image(times2Plot,1:SegNb,t(skCSD.all),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n",yaxt="n")
if(m==1) mtext(paste0("R",R),side=3)
if(r==1) mtext(paste0("Lambda", Lambda),side=2)
if(r==length(R.all)) {axis(4,at=round(seq(1,SegNb-10,,4)),las=2)}
if(m==length(lambda.all)) axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)


#image(t(skCSD.all),main=paste("skCSD Lambda:",Lambda, "R",R),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n",yaxt="n")


}
}
dev.off()
}


}
#ImageAllParCurrReg(SimPath,locationsKernel,locationsData, 0,3,3)
#ImageAllParCurrReg(SimPath,locationsKernel,locationsData, 0,3,3)
######################################








ImageAllParCurrRegBS<-function(SimPath,locationsKernel,locationsData,SNR,  ElPlot){
png(paste0(SumPlotPlace,"/",CellType,"skCSDall.png"),width=1000, height=900, pointsize=22)

plot.new()
par(mfrow=c(length(1)+2,3),oma=c(4,2,4,4))
#Elplot<-c(1,3,5)
for(valami in 1:length(ElPlot)){
Dlength<-valami #ElPlot[valami]
inname<-paste0(SimPath, locationsKernel[Dlength],"/")
outname<-paste0(SimPath, locationsData[Dlength])
#if(file.exists(paste(inname,outname1,"/params",sep=""))==FALSE) next
#params<-as.matrix(read.table(paste(inname,outname1,"/params",sep="")))
OrigCurr<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/membcurr",sep="")))
 SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

seg.length<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/","seglength",sep="")))
funaramvonal<-function(x) x/seg.length
memb.currents.line1<-array(0, c(dim(OrigCurr)))
memb.currents.line1<-apply(OrigCurr,2,funaramvonal) 
memb.currents.line1[which(memb.currents.line1<(0.5*min(memb.currents.line1)))]<-0.5*min(memb.currents.line1)
#memb.currents.line1[which(memb.currents.line1<(0.7*max(memb.currents.line1)))]<-0.5*max(memb.currents.line1)
#OrigCurrS<-as.matrix(read.table(paste(outname,outname1,"/membcurr_smoothed",sep="")))

times2Plot<-as.matrix(read.table(paste0(inname,"time")))

#Select just few Rs and Ms for plotting
R.all<-2^(3:7)
M<-512
lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )



#,mar=c(1,1,1,1))
#par(mfrow=c(length(lambda.all)+2,length(R.all)),oma=c(2,3,4,2))#,mar=c(1,1,1,1))

#layout(rbind(c(1,2,3,4,5)),widths=c(0.5,0.2,1,3,0.2),respect=FALSE) 
#mtext(c("oma[SOUTH<-1]", "oma[WEST<-2]", "oma[NORTH<-3]", "oma[EAST<-4]"),     c(SOUTH<-1, WEST<-2, NORTH<-3, EAST<-4),      line=0.4, cex=1.2, col="green", outer=TRUE)

if(valami==2){
title(main="skCSD Reconstruction", outer=TRUE, cex=3)
#mtext(c("Time (ms)", "Segment #"), c(1,4), outer=TRUE,line=2.5)


 par(mfg=c(1,1))
col2<-ColoursDori(memb.currents.line1)[[1]]
ExtVal<-ColoursDori(memb.currents.line1)[[2]]
par(mar=c(0.1,0.1,0.1,0.1))
SegNb<-dim(memb.currents.line1)[1]
image(times2Plot,1:SegNb,t(memb.currents.line1),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", yaxt="n")#,axis.args=list( at=c(-ExtVal, ExtVal), labels=c("Sink", "Source")) )
mtext("Ground Truth",3)
axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
axis(4,round(seq(1,SegNb-10,,4)),las=2)
mtext("Segment ID",4,cex=0.8,line=2.5)
mtext("Time (ms)",1,cex=0.8,line=2.5)
col2<-ColoursDori(memb.currents.line1)[[1]]
ExtVal<-ColoursDori(memb.currents.line1)[[2]]

par(mfg=c(1,2),mar=c(0.1,0.1,0.1,15)) 
image.plot(times2Plot,1:SegNb,t(memb.currents.line1),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", yaxt="n",axis.args=list( at=c(-ExtVal, ExtVal), labels=c("Sink", "Source")),
 legend.only=TRUE)


#image.scale(t(sCSD),col=col2,axis.pos=4,zlim=c(-ExtVal,ExtVal),add.axis=FALSE)
#axis(4,at=c(-ExtVal,ExtVal),c("Sink", "Source"),las=2)



par(mfg=c(1,3))
par(mar=c(0.1,1,0.1,0.1))
#> BestEMatr[3,]
#[1] "0.125000000000001"  "0.0625000000000028" "0.0312500000000173"
#[4] "0.0156250000432575" "0.0156988326204275"

plot(1:length(data.matrix(BestEMatr[3,])),data.matrix(BestEMatr[3,]),pch=4, 
xlab="Number of Electrodes", ylab="CV Error", cex=2, lwd=2,,t="b",xaxt="n")
mtext("CV Error",3)
if (CellType=="BS") axis(1, at=c(1:5),labels=c(8,16,32,64, 128))
if (CellType=="Yreg") axis(1, at=c(1:4),labels=c(8,16,32,64))
mtext("Number of electrodes",1,cex=0.8,line=2.5)
mtext("CV Error",2,cex=0.8,,line=2.5)
#image(times2Plot,sCSDCoord[,3],t(sCSD),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n",ylab="z (um)")
#mtext("sCSD",3)
#axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
#par(mar=c(0,3,0,5)) 

}

#image(c(times2Plot),1:dim(memb.currents.line1)[1],array(0,c(length(times2Plot),dim(memb.currents.line1)[1])),xlab="Time",ylab="Segment #",col="WHITE")
# image(t(OrigCurrS),col=rainbow(500),main="Smoothed Ground Truth")
  par(mfg=c(3,1))

    skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))


   ETable<-read.table(paste0(outname,outname1,"/ErrorCV_Noise_SNR0"))
    IndexMin<-which(ETable==min(ETable), arr.ind=TRUE)
    Lambda<-lambda.all[IndexMin[2]]
    R<-R.all[IndexMin[1]]
    #Reading in the LFP and adding noise
    #Ktilda_M512_R55lambda1e-04
currName<-paste0(outname,outname1,"/skCSDall_M",M,"_R",R,"lambda",Lambda)
    if(file.exists(currName)) {
      skCSD.all.part<-as.matrix(read.table(currName)) 
      

 for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
    }      
      
      
      
      
                               }else next
col2<-ColoursDori(skCSD.all)[[1]]
ExtVal<-ColoursDori(skCSD.all)[[2]]
par(mfg=c(3,valami))
par(mar=c(0.1,0.1,0.1,0.1))
Els<-c(8, 32, 128)
image(times2Plot,1:seg.nb,t(skCSD.all),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n",yaxt="n")
mtext(paste(Els[valami], "Electrodes"),3)
mtext("Time (ms)",1,cex=0.8,line=2.5)
#if(m==1) mtext(paste0("R",R),side=3)
#if(r==1) mtext(paste0("Lambda", Lambda),side=2)
if(valami==3) {axis(4,at=round(seq(1,seg.nb-10,,4)),las=2)}
axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
mtext("Segment ID",4,cex=0.8,line=2.5)

#image(t(skCSD.all),main=paste("skCSD Lambda:",Lambda, "R",R),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n",yaxt="n")


#}
}

dev.off()



}
#ImageAllParCurrReg(SimPath,locationsKernel,locationsData, 0,3,3)
#ImageAllParCurrReg(SimPath,locationsKernel,locationsData, 0,3,3)
#ImageAllParCurrRegBS(SimPath,locationsKernel,locationsData,SNR,  c(2^(3:6)))
























#####################################

#Plotting the Errors vs the number if electrodes
#ElNum_ErrorPlot(SimPath,locationsKernel,locationsData,0)
ElNum_ErrorPlot<-function(SimPath,locationsKernel,locationsData,SNR){

BesterrorCV<-array(NA,c(3,length(locationsData)))

for(Dlength in 1:length(locationsKernel)){

inname<-paste0(SimPath, locationsKernel[Dlength],"/")
outname<-paste0(SimPath, locationsData[Dlength])
lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
R.all<-2^(3:7)
M.all<-lambda.all#[5]
#if((file.exists(paste(outname,outname1,"/EL2ErrorSoma","_SNR", SNR, "_lambda",lambda,sep="")) | file.exists(paste(outname,outname1,"/EL2ErrorSoma","_SNR",
# SNR, "_lambda",lambda,sep="") ))==FALSE) next
EL2Error<-read.table(paste(outname,outname1,"/ErrorCV_Noise","_SNR",SNR,sep=""))


BesterrorCV[,Dlength]<-c( R.all[which(EL2Error==min(EL2Error),arr.ind=TRUE)[1]], 
                        M.all[which(EL2Error==min(EL2Error),arr.ind=TRUE)[2]],min(EL2Error))

write.table(BesterrorCV, paste(outname,outname1,"/BesterrorCV","_SNR",
                             SNR,sep=""))
png(paste0(outname,"/",locationsData[Dlength],"BesterrorCV_ErrorvsR_lambda",ElNumb,"_SNR", SNR, "_lambda",lambda),width=600, height=600,pointsize=20)
par(mar=c(4,4,4,6))
yLim<-range(BesterrorCV)
matplot(R.all, EL2Error,main="Error dependence on basis width",t="l",
        xlab="Basis Width",ylab="L2 error",lwd=2)#,ylim=yLim)
par(xpd=TRUE)
legend("topleft",c(paste(lambda.all)),pch=20,col=1:6,title="lambda")
par(xpd=FALSE)
dev.off()



ChosenError<-BesterrorCV #E2SegNormSoma

}




BestEMatr<-matrix(c(BesterrorCV),nrow=3)

return(list(BestEMatr=BestEMatr))
}


###############################
ElNum_ErrorPlotL1<-function(SimPath,locationsKernel,locationsData,SNR){

BesterrorCV<-array(NA,c(3,length(locationsData)))

for(Dlength in 1:length(locationsKernel)){

inname<-paste0(SimPath, locationsKernel[Dlength],"/")
outname<-paste0(SimPath, locationsData[Dlength])
lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
R.all<-2^(3:7)
M.all<-lambda.all#[5]
#if((file.exists(paste(outname,outname1,"/EL2ErrorSoma","_SNR", SNR, "_lambda",lambda,sep="")) | file.exists(paste(outname,outname1,"/EL2ErrorSoma","_SNR",
# SNR, "_lambda",lambda,sep="") ))==FALSE) next
EL2Error<-read.table(paste(outname,outname1,"/ErrorL1Smoothed_Noise_SNR0",sep=""))


BesterrorCV[,Dlength]<-c( R.all[which(EL2Error==min(EL2Error),arr.ind=TRUE)[1]], 
                        M.all[which(EL2Error==min(EL2Error),arr.ind=TRUE)[2]],min(EL2Error))

write.table(BesterrorCV, paste(outname,outname1,"/BesterrorL1Smoothed","_SNR",
                             SNR,sep=""))
png(paste0(outname,"/",locationsData[Dlength],"BesterrorL1Smoothed",ElNumb,"_SNR", SNR, "_lambda",lambda, ".png"),width=600, height=600,pointsize=20)
par(mar=c(4,4,4,6))
yLim<-range(BesterrorCV)
matplot(R.all, EL2Error,main="Error dependence on basis width",t="l",
        xlab="Basis Width",ylab="L1 error",lwd=2)#,ylim=yLim)
#par(xpd=TRUE)
legend("bottomright",c(paste(lambda.all)),pch=20,col=1:6,title="lambda")
#par(xpd=FALSE)
dev.off()



ChosenError<-BesterrorCV #E2SegNormSoma

}




BestEMatr<-matrix(c(BesterrorCV),nrow=3)

return(list(BestEMatr=BestEMatr))
}

#ElNum_ErrorPlotL1(SimPath,locationsKernel,locationsData,SNR)



###############################
####



