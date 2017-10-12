

ScriptLocation<-"/home/csdori/ksCSD_2014/trunk/"
  SimPath<-paste0(ScriptLocation,"simulation/")
  ElNumb<-c(441,25,25,25,49,81,81,25)
  #locationsKernel<-c("Ganglion_128Regular","Ganglion_256Regular") 
  locationsKernel<-c("Berd40_21x21_40","gang_5x5_50","gang_5x5_100","gang_9x9_100","gang_5x5_200", "gang_7x7_200","gang_9x9_200", "gang_5x5_400")#"Berd40_21x21"#c("gang_9x9_200","gang_9x9_50","gang_9x9_100","gang_5x5_25","gang_5x5_50","gang_5x5_100","gang_5x5_200","gang_5x5_400")#
  #locationsKernel<-c("gang_10x10_50","gang_10x10_100","gang_20x20_25","gang_20x20_50")#
  locationsData<-locationsKernel

IED<-c(40,50,100,100,200,200,200,400)


if(CellType=="Berdondini"){
  SimPath<-paste0(ScriptLocation,"simulation/")
  ElNumb<-441
  #locationsKernel<-c("Ganglion_128Regular","Ganglion_256Regular") 
  locationsKernel<-"Berd40_21x21"#c("gang_9x9_200","gang_9x9_50","gang_9x9_100","gang_5x5_25","gang_5x5_50","gang_5x5_100","gang_5x5_200","gang_5x5_400")#

  locationsData<-locationsKernel
  outname1<-"/kernelOut_4" ##_kCSD"
}




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
EL2Error<-read.table(paste0(outname,outname1,"/ErrorL1Smoothed_Noise_SNR0_1_1000"))#"/ErrorL1Smoothed_Noise_SNR0",sep=""))


BesterrorCV[,Dlength]<-c( R.all[which(EL2Error==min(EL2Error),arr.ind=TRUE)[1]], 
                        M.all[which(EL2Error==min(EL2Error),arr.ind=TRUE)[2]],min(EL2Error))

write.table(BesterrorCV, paste(outname,outname1,"/BesterrorL1Smoothed","_SNR",
                             SNR,sep=""))
png(paste0(outname,"/",locationsData[Dlength],"BesterrorL1Smoothed",ElNumb,"_SNR", SNR, ".png"),width=600, height=600,pointsize=20)
par(mar=c(4,4,4,6))
yLim<-range(BesterrorCV)
matplot(R.all, EL2Error,main="Error dependence on basis width",t="l",
        xlab="Basis Width",ylab="L1 error",lwd=2)#,ylim=yLim)
#par(xpd=TRUE)
legend("topright",c(paste(lambda.all)),pch=20,col=1:6,title="lambda")
#par(xpd=FALSE)
dev.off()



ChosenError<-BesterrorCV #E2SegNormSoma
}
BestEMatr<-matrix(c(BesterrorCV),nrow=3)

return(list(BestEMatr=BestEMatr))
}
alma<-ElNum_ErrorPlotL1(SimPath,locationsKernel,locationsData,SNR)




par(mar=c(8,3,3,3))
plot(1:length(IED),as.numeric(alma$BestEMatr[1,]),pch=20, ylab="R",xaxt="n", yaxt="n",xlab="",main="Basis width in various setups")#,ylim=c(0.2,0.82))
axis(1,at=c(1:length(IED)),locationsKernel,las=2)
axis(2,at=c(unique(as.numeric(alma$BestEMatr[1,]))),c(unique(as.numeric(alma$BestEMatr[1,]))),las=2)




