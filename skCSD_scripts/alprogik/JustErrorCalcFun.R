#source(paste0(ScriptLocation,"alprogik/JustErrorCalcFun.R"))
#JustErrorCalcFun(SimPath,locationsKernel,locationsData,SNR, lammbda)
JustErrorCalcFun<-function(SimPath,locationsKernel,locationsData,SNR, lammbda){
Reg<-TRUE

ErrorResults <- vector("list", length(locationsKernel))
SmoothR<- vector("list", length(locationsKernel))
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


BandwidtAll<-c(sqrt(SmoothPar[2]/SmoothPar[3]),2*sqrt(SmoothPar[2]/SmoothPar[3]), SmoothElWidth/2, SmoothElWidth/4, sqrt(SmoothElWidth)/2 , 10,30,50,70,90)

BandwidtAllName<-c(paste("Smoothed",c(1:length(BandwidtAll))))




EL2ErrorVarR<-array(0,c((length(BandwidtAll)+1)*length(R.all),length(M.all)))

#Reading in the LFP and adding noise
    LFP<-as.matrix(read.table(paste(outname,"/myLFP",sep="")))
    
if (SNR!=0){
      LFPvariance<-var(c(LFP))
      Noisevariance<-LFPvariance/SNR
      NoiseLFP<-rnorm(length(c(LFP)), mean = 0, sd = sqrt(Noisevariance))
      LFP<-LFP+NoiseLFP
    }
LFPOriginal<-LFP

for(m in 1:length(M.all)){
  for(r in 1:length(R.all)){
    LFP<-LFPOriginal
    M<-M.all[m]
    R<-R.all[r]
 el2ignorename<-paste(inname,outname1,"/El2Ignore_M",M,"_R",R,sep="")
    if(file.exists(el2ignorename)){
      El2Ignore<-c(as.matrix(read.table(el2ignorename)))
      LFP<-LFP[-El2Ignore,] }
      
      
     
 skCSD.all<-as.matrix(read.table(paste0(outname,outname1,"/skCSDall_M",M,"_R",R,"_SNR", SNR, "_lambda",lambda)))
    #Reading in tge currents
    
    #skCSD.Smoothed<-KernSmoothDistanceOther(SmoothElWidth/4,inname, outname, outname1,"skCSD",skCSD.all)
    #######################Calculate Currents
    ############################################
    #Calculate Errors
 
for(bw in 1:length(BandwidtAll)){

memb.currents.smoothed<-KernSmoothDistance(BandwidtAll[bw],inname, outname, outname1,BandwidtAllName[bw])
EL2ErrorVarR[r+(bw-1)*length(M.all),m]<-L2Error(skCSD.all,memb.currents.smoothed)[[1]]
#EL2ErrorVarR[r+(bw-1)*length(M.all),m]<-L2Error(skCSD.Smoothed,memb.currents.smoothed)[[1]]
#image(t(memb.currents.smoothed),col=rainbow(40))
}


######################CV error 
Kmatr<-as.matrix(read.table(paste(inname,outname1,"/K_M",M,"_R",R,sep="")))
    Ktildematr<-as.matrix(read.table(paste(inname,outname1,"/Ktilda_M",M,"_R",R,sep="")))
      
      Tmatr<-try(1/(4*pi*0.5)*t(Ktildematr)%*%ginv(Kmatr))
      
EL2ErrorVarR[r+bw*length(M.all),m]<-try(cross.valid(LFP,Tmatr,1:dim(skCSD.all)[1])[[1]],silent=TRUE)



   }
   }

    ErrorResults[[Dlength]] <-EL2ErrorVarR
SmoothR[[Dlength]]<-BandwidtAll
    }
return(list(ErrorResults=ErrorResults,SmoothR=SmoothR ))

}

alma<-JustErrorCalcFun(SimPath,locationsKernel,locationsData,SNR, lammbda)

break

Which2Plot<-4
matplot(alma[[1]][[Which2Plot]],t="l",ylim=c(0,2))
 #abline(v=c((1:3)*length(M.all)+0.5))
#axis(side=1, labels=rep(M.all,4),at=c(1:(length(M.all)*4)))
mtext(round(alma[[2]][[Which2Plot]],3), side=3, at=c(1:length(alma[[2]][[Which2Plot]])*length(M.all)-length(M.all)/2),line=-1)
  abline(v=c((1:9)*length(M.all)+0.5),col="ORANGE")
  #par(xpd=TRUE)
  #legend(4*length(R.all)+2,median(whichError2Plot),c(paste(M.all)),pch=20,col=1:6,title="M")

#BandwidtAll