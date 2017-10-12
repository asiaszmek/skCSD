



ImageAllParCurrRegBS<-function(SimPath,locationsKernel,locationsData,SNR,  ElPlot,NonLinPar){
  TSize<-1.3
  
  png(paste0(SumPlotPlace,"/",CellType,"skCSDall.png"),width=14, height=10,units="in" ,res=300)
  Errors<-array(0,c(2,length(locationsKernel)))
  for(loci in 1: length(locationsKernel)){
    inname<-paste0(SimPath, locationsKernel[loci],"/")
    outname<-paste0(SimPath, locationsData[loci])
    ETable<-as.matrix(read.table(paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,sep="")))
    CVTable<-as.matrix(read.table(paste(outname,outname1,"/ErrorCV_Noise","_SNR",SNR,sep="")))
    
    Errors[1,loci]<-min(ETable)
    Errors[2,loci]<-min(CVTable)
  }
  
  
  
  plot.new()
  par(mfrow=c(length(1)+2,3),oma=c(4,2,4,4),cex=1.1)
  #Elplot<-c(1,3,5)
  
  Dlength<-1 #ElPlot[valami]
  inname<-paste0(SimPath, locationsKernel[Dlength],"/")
  outname<-paste0(SimPath, locationsData[Dlength])
  #if(file.exists(paste(inname,outname1,"/params",sep=""))==FALSE) next
  #params<-as.matrix(read.table(paste(inname,outname1,"/params",sep="")))
  OrigCurr<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/membcurr",sep="")))
  SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
  seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
  SegNb <-seg.nb
  seg.length<-as.matrix(read.table(paste(SimPath,"/",locationsData[Dlength],"/","seglength",sep="")))
  funaramvonal<-function(x) x/seg.length
  memb.currents.line1<-array(0, c(dim(OrigCurr)))
  memb.currents.line1<-apply(OrigCurr,2,funaramvonal) 
  
  
  source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))
  
  CurrNew<-KernSmoothDistance(10,inname, outname, outname1,"30")
  Currnew1<-CurrNew
  #CurrNew[which(CurrNew<(0.5*min(CurrNew)))]<-0.5*min(CurrNew)
  
  
  #memb.currents.line1[which(memb.currents.line1<(0.5*min(memb.currents.line1)))]<-0.5*min(memb.currents.line1)
  ##########Change!!!!
  memb.currents.line1<- CurrNew[c(2,1,3:52),]
  
  
  #memb.currents.line1[which(memb.currents.line1<(0.7*max(memb.currents.line1)))]<-0.5*max(memb.currents.line1)
  #OrigCurrS<-as.matrix(read.table(paste(outname,outname1,"/membcurr_smoothed",sep="")))
 
  times2Plot<-as.matrix(read.table(paste0(inname,"time")))
  
  #Select just few Rs and Ms for plotting
  R.all<-2^(3:7)
  M<-512
  lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
  
  
  
  par(mfg=c(1,1))
  col2<-ColoursDori(memb.currents.line1)[[1]]
  ExtVal<-ColoursDori(memb.currents.line1)[[2]]
  BreaksIm<-tan(seq(-ExtVal, ExtVal,,201)/ExtVal*(pi/2-0.1)*NonLinPar)/tan((pi/2-0.1)*NonLinPar)*ExtVal
  
  par(mar=c(0.1,0.1,0.1,0.1),cex=1.1)
  SegNb<-dim(memb.currents.line1)[1]
  image(times2Plot,1:SegNb,t(memb.currents.line1),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", yaxt="n",breaks=BreaksIm)#,axis.args=list( at=c(-ExtVal, ExtVal), labels=c("Sink", "Source")) )
  mtext("Smoothed Ground Truth",3,cex=TSize)
  axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
  axis(4,round(seq(1,SegNb-10,,4)),las=2)
  mtext("Segment ID",4,cex=TSize,line=2.5)
  mtext("Time (ms)",1,cex=TSize,line=2.5)
  col2<-ColoursDori(memb.currents.line1)[[1]]
  ExtVal<-ColoursDori(memb.currents.line1)[[2]]
  BreaksIm<-tan(seq(-ExtVal, ExtVal,,201)/ExtVal*(pi/2-0.1)*NonLinPar)/tan((pi/2-0.1)*NonLinPar)*ExtVal
  
  par(mfg=c(1,2),mar=c(0.1,0.1,0.1,18)) 
  image.plot(times2Plot,1:SegNb,t(memb.currents.line1),col=col2,zlim=c(-ExtVal,ExtVal),breaks=BreaksIm,xaxt="n", yaxt="n",axis.args=list( at=c(-ExtVal, ExtVal), labels=c("Sink", "Source")), legend.only=TRUE,cex=TSize*1.2)
  
  
  par(mar=c(0.1,0.1,0.1,0.1))
   par(mar=c(0.1,0.1,0.1,0.1),oma=c(4,2,4,4),cex=1.1)
  par(mfg=c(1,3))
 
  
  matplot(t(Errors),pch=4, xlab="Number of Electrodes", ylab=" Error", cex=1, lwd=2,,t="b",xaxt="n")
  legend("topright",c("L1 Error", "CV Error"),bty="n", col=1:2, pch=20,cex=1)
  #plot(1:length(data.matrix(BestEMatr[3,])),data.matrix(BestEMatr[3,]),pch=4, 
  #xlab="Number of Electrodes", ylab="CV Error", cex=2, lwd=2,,t="b",xaxt="n")
  mtext("Reconstruction Error",3,cex=TSize)
  if (CellType=="BSOLD") axis(1, at=c(1:5),labels=c(8,16,32,64, 128),cex=1,las=2)
  if (CellType=="Yreg") axis(1, at=c(1:4),labels=c(8,16,32,64),las=2)
  mtext("Number of electrodes",1,cex=TSize,line=2.5)
  
  mtext("Error",2,cex=TSize,line=2.5)
  
  for(valami in 1:length(ElPlot)){
    
    Dlength<-valami #ElPlot[valami]
    inname<-paste0(SimPath, locationsKernel[Dlength],"/")
    outname<-paste0(SimPath, locationsData[Dlength])
    
    par(mfg=c(3,1))
    
    skCSD.all<-array(0,c(seg.nb, dim(memb.currents.line1)[2]))
    
    
    ETable<-as.matrix(read.table(paste(outname,outname1,"/ErrorCV_Noise","_SNR",SNR,sep="")))#as.matrix(read.table(paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,sep="")))
    
    SameplaceAll<-readLines(paste0(outname,outname1,"/sameplace.txt"))
    #   ETable<-read.table(paste0(outname,outname1,"/ErrorCV_Noise_SNR0"))
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
    skCSD.all<-skCSD.all[c(2,1,3:52),]
    col2<-ColoursDori(skCSD.all)[[1]]
    ExtVal<-ColoursDori(skCSD.all)[[2]]
    BreaksIm<-tan(seq(-ExtVal, ExtVal,,201)/ExtVal*(pi/2-0.1)*NonLinPar)/tan((pi/2-0.1)*NonLinPar)*ExtVal
    
    par(mfg=c(3,valami))
    par(mar=c(0.1,0.1,0.1,0.1))
    Els<-c(8, 32, 128)
    image(times2Plot,1:seg.nb,t(skCSD.all),col=col2,zlim=c(-ExtVal,ExtVal),breaks=BreaksIm,xaxt="n",yaxt="n")
    mtext(paste(Els[valami], "Electrodes"),3,cex=TSize)
    mtext("Time (ms)",1,cex=TSize,line=2.5)
    #if(m==1) mtext(paste0("R",R),side=3)
    #if(r==1) mtext(paste0("Lambda", Lambda),side=2)
    if(valami==3) {axis(4,at=round(seq(1,seg.nb-10,,4)),las=2)}
    axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
    mtext("Segment ID",4,line=2.5,cex=TSize)
    
  }
  dev.off()
  
}




ImageAllParCurrRegBS(SimPath,locationsKernel,locationsData, 0,c(8,32,128),0.7)

#
#c(2,1,3:52)

