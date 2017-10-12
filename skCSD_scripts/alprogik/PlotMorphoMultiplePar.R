#Plot at a certain time points the various reconstruction

source("alprogik/Colors_BlueRed.R")
library(data.table)
#image scale function for plotting color bars for image plots


#DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/Hex_12x10_inter50/"#BS_d50_el128/"
#DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/Y_el4x16_RotS_symm_d50_ver0/"
DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/Ganglion_d17.5_256Regular/"

#DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_el4x16_Rot_symm_d50_ver0/"
# DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/Sim_Yshaped/Y_el4x16_Elrand_symm_d50_ver4/"
#DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Hex_8x8_inter70_d30/"
plotFileName<-"Comparing"
dir.create(paste0(DataFolder,plotFileName))
CellType<-"gangNew" #"Y-rotated" #


SimTime<-as.matrix(read.table(paste0(DataFolder,'/time')))
TimeMax<-length(SimTime)
#TimeInstant<-which(SimTime==45) #45


if(CellType=="gangNew"){
  xA<-1
  yA<-2
  xAEl<-1 #1
  yAEl<-2
  PlotTitle<-"Ganglion Cell"
  outname1<-"/kernelOut_poli_CP3/"
  ToPlot<--c(1:10,70)
  TimeInstant<-266
}

if(CellType=="BS"){
  xA<-1
  yA<-3
  xAEl<-2 #1
  yAEl<-3
  ToPlot<-1:52
  PlotTitle<-"BS"
  outname1<-"/kernelOut_poli_CP3/"
}

if(CellType=="Y"){
  xA<-1
  yA<-3
  xAEl<-1
  yAEl<-3
  ToPlot<-1:86
  PlotTitle<-"Y-shaped"
  outname1<-"/kernelOut_poli_CP2/"
}


if(CellType=="Y-rotated"){
  xA<-1
  yA<-3
  xAEl<-2 #1
  yAEl<-3
  ToPlot<-1:86
  PlotTitle<-"Y-shaped"
  outname1<-"/kernelOut_poli_CP3/"
  TimeInstant<-which(SimTime==45) #45
  TimeInstant<-TimeInstant + 2
}


El2Ignore<-numeric()
segstart<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsstart_x_y_z'), colClasses='numeric')),ncol=3)
segend<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsend_x_y_z'), colClasses='numeric')),ncol=3)
segmid<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsmid_x_y_z'), colClasses='numeric')),ncol=3)
#Plotratio<-Xrange/max(segwidth)
segwidth1<-as.matrix(read.table(paste0(DataFolder,'segdiam_x_y_z'), colClasses='numeric'))
segwidth<-segwidth1/max(segwidth1)*30 # not realistic?
elec<-matrix(as.matrix(read.table(paste0(DataFolder,'elcoord_x_y_z'), colClasses='numeric')),ncol=3)
membCurr<-as.matrix(read.table(paste0(DataFolder,'/membcurr')))
funaramvonal<-function(x) x/seg.length
membcurr0<-apply(membCurr,2,funaramvonal)
membcurr<-membcurr0[ToPlot,] 

dimseg<-dim(segstart)[1]
Xrange<-range(segstart[,xA],elec[,xA])*1.05
Yrange<-range(segstart[,yA],elec[,yA])*1.05
XrangeEl<-range(segstart[,xAEl],elec[,xAEl])*1.05
YrangeEl<-range(segstart[,yAEl],elec[,yAEl])*1.05

seg.nb<-dim(segmid)[1]



PlotMultiPar(14)
PlotMultiPar(98*2)


PlotMultiPar<-function(TimeInstant){

########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# run a  for loop through the various reconstructions
VarOptionsR<-c(seq(5,130,25))
VarOptionsL<-c(1e-4,0.1)
png(paste0(DataFolder,"DiffParTime_",TimeInstant,".png"),width=1200, height=800)
par(mfrow=c(length(VarOptionsL)+1,length(VarOptionsR)))


#Plot Ground Truth
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main="Ground Truth",asp=1)


rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
# col1<-ColoursDori(membcurr[,TimeInstant])[[1]]
# ExtVal<-ColoursDori(membcurr[,TimeInstant])[[2]]
col1<-ColoursDori(membcurr[,TimeInstant])[[1]]
ExtVal<-ColoursDori(membcurr[,TimeInstant])[[2]]

#highlight the interesting part of the plot
symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(Xrange)/10* abs(membcurr[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
points(elec[,xA],elec[,yA], pch=8) #diff(Xrange)

#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))

szinskala2<-color.scale(c(membcurr[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
segments(segstart[c(1:dimseg)[ToPlot],xA], segstart[c(1:dimseg)[ToPlot],yA], segend[c(1:dimseg)[ToPlot],xA], segend[c(1:dimseg)[ToPlot],yA], col= szinskala2,lwd=segwidth[c(1:dimseg)[ToPlot]],lend=1)
#segments(segstart[1:dimseg,xA], segstart[1:dimseg,yA], segend[1:dimseg,xA], segend[1:dimseg,yA], col= szinskala2,lwd=segwidth[1:dimseg])
#par(mar=c(0.5,0.1,1,0.5)) 
mtext("mA/um ??",side=4,line=1,cex=0.7)

ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
# add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))


plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')


for(VarL in 1:length(VarOptionsL)){
for(VarRec in 1:length(VarOptionsR)){

  #skCSD<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSD_plain")))#paste0(outname,outname1,"/skCSDall_M",256,"_R",55,"_SNR0_lambda0")))
  C.calc<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSDall_M512_R", 
                                      VarOptionsR[VarRec], "lambda",VarOptionsL[VarL])))#30lambda1e-04")))#80lambda0.001")))
  SameplaceAll<-readLines(paste0(DataFolder,outname1,"/sameplace.txt"))
  skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
  for(i in 1: seg.nb){
    sameplace<-numeric() 
    sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
    #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
    if(length(sameplace)==1) skCSD.all[i,]<-C.calc[sameplace,]
    if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums(C.calc[c(sameplace),]),nrow=1)
  }
  skCSD<-skCSD.all[ToPlot,]
  ##############!!!!!!!!!!!!!!!!!!!!!
  
  
  
  
  
  coordsEnd<-as.matrix(read.table(paste0(DataFolder,"/coordsend_x_y_z")))
  
  
  
  #TimeInstant<-TimeInstant+1
  plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main=paste("R",VarOptionsR[VarRec],
                                                                                 "L", VarOptionsL[VarL] ),asp=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  # col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
  # ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]
  col1<-ColoursDori(skCSD)[[1]]
  ExtVal<-ColoursDori(skCSD)[[2]]
  symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(Xrange)/10* abs(skCSD[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
  
  
  #col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))
  
  szinskala2<-color.scale(c(skCSD[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
  points(elec[,xA],elec[,yA], pch=8)
  segments(segstart[c(1:dimseg)[ToPlot],xA], segstart[c(1:dimseg)[ToPlot],yA], 
           segend[c(1:dimseg)[ToPlot],xA], segend[c(1:dimseg)[ToPlot],yA], 
           col= szinskala2,lwd=segwidth[c(1:dimseg)[ToPlot]],lend=1)
  #par(mar=c(0.5,0.1,1,0.5)) 
  
  #mtext("mA/um ??",side=4,line=1,cex=0.7)
  
  ColMatr<-as.matrix(c(1:256))
  col2<-ColoursDori(ColMatr)[[1]]
  #add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
  # add.image(t(ColMatr),adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
  # #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
  # # axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
  # axis(4,labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
  add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, 
            adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
  #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
  axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c("sink","source"))
  
  
}
}
dev.off()

}

#####################################x
######################Plot Effect of Noise

for(Dlength in 1:length(locationsKernel)){
  png(paste0(SumPlotPlace,'/NoiseCompare_',locationsKernel[Dlength],'_R30_L0.01.png'), height=800, width=600)
  NoisePar<-c(0,100,50,20,10,5)
  par(mfrow=c(length(NoisePar),3),mar=c(1,1,1,1),oma=c(2,2,2,2))
  for (NS in 1: length(NoisePar)){
    SNR<-NoisePar[NS]
  inname<-paste0(SimPath, locationsKernel[Dlength],"/")
  outname1<-"/kernelOut_poli_CP3"
  MinIndex<-data.matrix(read.table(paste(inname,outname1,"/MinIndex_Noise","_SNR",SNR,sep="")))
  lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
  R.all<-seq(5,130,25)
  Rmin<-R.all[MinIndex[1]]
  Lmin<-lambda.all[MinIndex[2]]
  RChosen<-30
  LChosen<-0.01 #1e-04
  M<-512
LFP<-data.matrix(read.table( paste(inname,outname1,"/LFP_Noise","_SNR",SNR,sep="")))
skCSD.all.part<-data.matrix(read.table( paste(inname,outname1,"/ksCSD_Matr",M,"_R",Rmin,"lambda",
                               Lmin, "SNR",SNR,sep="")))
skCSD.all.part.Chosen<-data.matrix(read.table( paste(inname,outname1,"/ksCSD_Matr",M,"_R",RChosen,"lambda",
                                             LChosen, "SNR",SNR,sep="")))


SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

skCSD<-array(0,c(seg.nb, dim(LFP)[2]))
skCSDChosen<-array(0,c(seg.nb, dim(LFP)[2]))
#Reading in tge currents

for(i in 1: seg.nb){
  sameplace<-numeric() 
  sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
  #Plot the current related ti the different segments
  #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
  if(length(sameplace)==1) skCSD[i,]<- skCSD.all.part[sameplace,]
  if(length(sameplace)>1) skCSD[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
  
  if(length(sameplace)==1) skCSDChosen[i,]<- skCSD.all.part.Chosen[sameplace,]
  if(length(sameplace)>1) skCSDChosen[i,]<-as.matrix(colSums( skCSD.all.part.Chosen[c(sameplace),]),nrow=1)
  
}





#LimLFP<-max(abs(LFP))


#LimCSD<-max(abs(ksCSD))
#LimCSDChosen<-max(abs(skCSDChosen))

col1<-ColoursDoriLFP(LFP)[[1]]
LimLFP<-ColoursDoriLFP(LFP)[[2]]
image(t(LFP),col=col1,zlim=c(-LimLFP,LimLFP),xaxt='n',main=paste('LFP',SNR))

col1<-ColoursDori(skCSD)[[1]]
LimCSD<-ColoursDori(skCSD)[[2]]
image(t(skCSD),col=col1,main=paste(Rmin, Lmin), zlim=c(-LimCSD,LimCSD),xaxt='n',yaxt='n')
col1<-ColoursDori(skCSDChosen)[[1]]
LimCSDChosen<-ColoursDori(skCSDChosen)[[2]]
image(t(skCSDChosen),col=col1,main=paste(RChosen, LChosen), zlim=c(-LimCSDChosen,LimCSDChosen),xaxt='n',yaxt='n')


  }
  
  dev.off()
}

#Plot curves of CV
#x axis = SNR, y error, lines= number of electrodes
NoisePar<-c(0,50,20,10,5)
MinErrorArray<-array(0,c(length(locationsKernel),length(NoisePar)))
for(Dlength in 1:length(locationsKernel)){
 # png(paste0(SumPlotPlace,'/NoiseCompare_',locationsKernel[Dlength],'.png'), height=800, width=600)
 
  for (NS in 1: length(NoisePar)){
    SNR<-NoisePar[NS]
    inname<-paste0(SimPath, locationsKernel[Dlength],"/")
    outname1<-"/kernelOut_4"
    MinIndex<-data.matrix(read.table(paste(inname,outname1,"/MinIndex_Noise","_SNR",SNR,sep="")))
    ErrorCV<-data.matrix(read.table(paste(inname,outname1,"/ErrorCV_Noise","_SNR",SNR,sep="")))
    MinError<-ErrorCV[MinIndex[1]]
    MinErrorArray[Dlength,NS]<-MinError
    lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
    R.all<-2^(3:7)

    
  }
  
}

png(paste0(SumPlotPlace,'/NoiseCompare_',CellType,'_Symm.png'))
matplot(c(8,16,32,64),MinErrorArray,t='l', xlab='Number of Electrodes', ylab='CV Error', main="Effect of Noise", lwd=2)
legend("topright",c('No noise',paste(NoisePar[-1])), pch=20, col=c(1:length(NoisePar)),title='SNR')
dev.off()


write.table(ErrorCV, paste(inname,outname1,"/ErrorCV_Noise","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)

MinIndex<-which(ErrorCV==min(ErrorCV),arr.ind=TRUE)
write.table(MinIndex, paste(inname,outname1,"/MinIndex_Noise","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)

