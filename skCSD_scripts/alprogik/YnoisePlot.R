#source("/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/alprogik/CurrentsonMorpho2.R")
#source("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/alprogik/Colors_BlueRed.R")
source("alprogik/Colors_BlueRed.R")
library(data.table)
library(RColorBrewer)
#image scale function for plotting color bars for image plots


#DataFolder<-#"/home/csdori/ksCSD_2014/trunk/simulation/Ganglion_d17.5_128Regular/"#BS_d50_el128/"
DataFolder<-"/home/csdori/ksCSD_2014/trunk/Sim_Yshaped/Y_el4x8_Elrand_symm_d50_ver0/"

#DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/Y_el4x16_RotS_symm_d50_ver0/"
#DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_el4x16_Rot_symm_d50_ver0/"
 #DataFolder<-"/home/csdori/ksCSD_2014/trunk/simulation/GanglionOsc_128Regular/"
 #ffmpeg  -framerate 10 -pattern_type glob -i 'CSD_Morpho*.png' CSD.mp4
#DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Hex_8x8_inter70_d30/"
plotFileName<-"YNOise"
dir.create(paste0(DataFolder,plotFileName))
CellType<-"Y"#"Y" # "gangNew" #"Y-rotated" #




if(CellType=="Y"){
xA<-1
yA<-3
xAEl<-1
yAEl<-3
ToPlot<-1:86
PlotTitle<-"Y-shaped"
  outname1<-"/kernelOut_4/"
}




segstart<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsstart_x_y_z'), colClasses='numeric')),ncol=3)
segend<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsend_x_y_z'), colClasses='numeric')),ncol=3)
segmid<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsmid_x_y_z'), colClasses='numeric')),ncol=3)
elec<-matrix(as.matrix(read.table(paste0(DataFolder,'elcoord_x_y_z'), colClasses='numeric')),ncol=3)
seg.nb<-dim(segmid)[1]

seglength<-as.matrix(read.table(paste0(DataFolder,'/seglength')))
funaramvonal<-function(x) x/seglength

membCurr<-as.matrix(read.table(paste0(DataFolder,'/membcurr')))
membcurr<-apply(membCurr,2,funaramvonal)

SimTime<-as.matrix(read.table(paste0(DataFolder,'/time')))
TimeMax<-length(SimTime)
TimeInstant<-which(SimTime==5) #45
TimeInstant<-TimeInstant+2
########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#skCSD<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSD_plain")))#paste0(outname,outname1,"/skCSDall_M",256,"_R",55,"_SNR0_lambda0")))


 # Error2Plot<-as.matrix(read.table(paste0(DataFolder,"/ErrorTable")))
#32-es símításos
#CurrMin<- min(Error2Plot[(2*length(R.all)+1):(3*length(R.all)),])
#MinIndex<-which(Error2Plot[(2*length(R.all)+1):(3*length(R.all)),]==CurrMin,arr.ind=TRUE)
  
  
  
#MinIndex<-data.matrix(read.table(paste(DataFolder,outname1,"/MinIndex_Noise","_SNR",SNR,sep="")))
M<-512
lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
R.all<-2^(3:7)





Rread<-32#R.all[MinIndex[1]]
Lread<-0.01 #lambda.all[MinIndex[2]]
        
#if(CellType=="Y"  | CellType=="Y-rotated" | CellType=="BS" ){            
#png(plotname, height=500, width=1200, pointsize=20)
#TimeInstant<-500
#par(mfrow=c(1,6))
#}

Plotwidth<-10
membcurr<-SmoothedCurr

myPng <- function(..., width=8, height=8, res=300) {
  png(..., width=width*res/50, height=height*res/50, res=res,units="in")
}


png(paste0('/home/csdori/ksCSD_2014/trunk/plots/summary/','Y_NoisePlotErrorGray.png'), width=9, height=9,units='in',res=300)
Plotwidth<-3
#par(cex=1.3)
par(mfrow=c(2,3),cex=1.2)
#par(mar=c(3,3,3,3))
#Plot Ground truth
dimseg<-dim(segstart)[1]
Xrange<-range(segstart[,xA],elec[,xA])*1.05
Yrange<-range(segstart[,yA],elec[,yA])*1.05
XrangeEl<-range(segstart[,xAEl],elec[,xAEl])*1.05
YrangeEl<-range(segstart[,yAEl],elec[,yAEl])*1.05

plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main="Smoothed GT",asp=1)


rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
 col1<-ColoursDori(membcurr[,TimeInstant])[[1]]
 ExtVal<-ColoursDori(membcurr[,TimeInstant])[[2]]
#col1<-ColoursDori(membcurr)[[1]]
#ExtVal<-ColoursDori(membcurr)[[2]]

#highlight the interesting part of the plot
#symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(membcurr[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
 symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(membcurr[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
 points(elec[,xA],elec[,yA], pch=8) #diff(Xrange)

#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))

szinskala2<-color.scale(c(membcurr[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
segments(segstart[c(1:dimseg)[ToPlot],xA], segstart[c(1:dimseg)[ToPlot],yA], segend[c(1:dimseg)[ToPlot],xA], segend[c(1:dimseg)[ToPlot],yA], col= szinskala2,lend=1,lwd=5)#segwidth[c(1:dimseg)[ToPlot]])
#segments(segstart[1:dimseg,xA], segstart[1:dimseg,yA], segend[1:dimseg,xA], segend[1:dimseg,yA], col= szinskala2,lwd=segwidth[1:dimseg])
#par(mar=c(0.5,0.1,1,0.5)) 
#mtext("mA/um ",side=4,line=1,cex=0.7)

ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
# add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)
axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c(round(-ExtVal,3),round(ExtVal,3)))




par(mfg=c(1,3))
ErrorTableNoise<-as.matrix(read.table( paste("/home/csdori/ksCSD_2014/trunk/Sim_Yshaped","/EL1Error_Noise",sep="")))

matplot(t(ErrorTableNoise),t='b',lwd=2,pch=20,xaxt='n',xlab='SNR', ylab='L1 Error', main='L1 Error')
axis(1,at=1:5,c("No Noise", "64","16","4","1"),las=2)
legend("topleft", c("2x4","4x4","4x8","4x16"),col=1:4,pch=20)#title="Electrodes"






for(snr in 1:4){
SNR<-c(0,16,4,1)[snr]#c(0,100,50,20,10,5)[snr]


if(file.exists("/home/csdori/ksCSD_2014/trunk/Sim_Yshaped/EL1Error_RTableNoise")){
BestR<-as.matrix(read.table("/home/csdori/ksCSD_2014/trunk/Sim_Yshaped/EL1Error_RTableNoise"))[,-2]
BestL<-as.matrix(read.table("/home/csdori/ksCSD_2014/trunk/Sim_Yshaped/EL1Error_LTableNoise"))[,-2]
Lread<-BestL[3,snr]
Rread<-BestR[3,snr]

}




if (snr==1){C.calc<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSDall_M512_R",Rread,"lambda",Lread))) #4*pi*0.5*
} 
if (snr!=1){
C.calc<-as.matrix(read.table(paste0(DataFolder,outname1,"/ksCSD_Matr512_R",Rread,"lambda",Lread,"SNR",SNR)))
}
SameplaceAll<-readLines(paste0(DataFolder,outname1,"/sameplace.txt"))
skCSD.all<-array(0,c(seg.nb, 561))
    for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<-C.calc[sameplace,]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums(C.calc[c(sameplace),]),nrow=1)
    }
skCSD<-skCSD.all










#screen(5)
#TimeInstant<-TimeInstant+1
if (snr==1) {
par(mfg=c(1,2))
plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main=paste("No Noise"),asp=1)}

if (snr!=1) {
par(mfg=c(2,snr-1))
plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main=paste("SNR =",SNR),asp=1)}
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
 col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
 ExtVal<-ColoursDori(skCSD[,TimeInstant])[[2]]
#col1<-ColoursDori(skCSD)[[1]]
#ExtVal<-ColoursDori(skCSD)[[2]] #ColoursDori(skCSD)[[2]]
#symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(skCSD[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)


#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))
ExtVal<-0.018
szinskala2<-color.scale(c(skCSD[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
symbols(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], circles=c(diff(range(segstart[,xA]))/Plotwidth* abs(skCSD[,TimeInstant])/ExtVal), add=TRUE, bg=colors()[182], fg=colors()[182], inches=FALSE)
points(elec[,xA],elec[,yA], pch=8)
segments(segstart[c(1:dimseg)[ToPlot],xA], segstart[c(1:dimseg)[ToPlot],yA], segend[c(1:dimseg)[ToPlot],xA], segend[c(1:dimseg)[ToPlot],yA], col= szinskala2,lend=1,lwd=5)#segwidth[c(1:dimseg
#par(mar=c(0.5,0.1,1,0.5)) lend=1,lwd=5)#segwidth[c(1:dimseg
#mtext("mA/um ",side=4,line=1,cex=0.7)

ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
#add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# add.image(t(ColMatr),adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col2)
# #axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
# # axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
# axis(4,labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))


#add.image(max(XrangeEl-10),mean(YrangeEl), t(ColMatr), adj.x = 0.0, adj.y = 0.5,  image.width = 1, image.height = 0.05, col = col2)



#axis(4,at=c(min(YrangeEl),max(YrangeEl)),labels=c(round(-ExtVal,3),round(ExtVal,3)))
}
#po
dev.off()
####################################
################### Error calculation
##################################################x
#goal: plot for error for best reconstruction in case of noise: 4x2-4x16 elecrodes in case of various noise level

DataFolder<-"/home/csdori/ksCSD_2014/trunk/Sim_Yshaped/Y_el4x8_Elrand_symm_d50_ver0/"
ScriptLocation<-"/home/csdori/ksCSD_2014/trunk/"
source(paste0(ScriptLocation,"alprogik/errorCalcFun.R"))
ErrorTableNoise<-array(0,c(4,5)) #Compare to smoothed?
LTableNoise<-array(0,c(4,5)) #Compare to smoothed?
RTableNoise<-array(0,c(4,5)) #Compare to smoothed?
M<-512
lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" ) #c(0.1,0.01,0.001)#
R.all<-2^(3:7) #c(8,16,128)#
segmid<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsmid_x_y_z'), colClasses='numeric')),ncol=3)

seg.nb<-dim(segmid)[1]
seglength<-as.matrix(read.table(paste0(DataFolder,'/seglength')))
funaramvonal<-function(x) x/seglength

membCurr<-as.matrix(read.table(paste0(DataFolder,'/membcurr')))
membcurr<-apply(membCurr,2,funaramvonal)
morpho<-as.matrix(read.table(paste0(DataFolder,'/segcoordinates.txt')))
DistMorpho<-as.matrix(dist(morpho))
#image(as.matrix(DistMorpho))
SmoothedCurr<-array(0,dim(membcurr))


SmoothingMatrix<-array(0,c(length(seglength),length(seglength)))
BandWidth<-32
for(i in 1:length(seglength)) SmoothingMatrix[i,]<-exp(-DistMorpho[i,]^2/2/BandWidth^2)
SmoothingMatrix<-SmoothingMatrix/colSums(SmoothingMatrix)
SmoothedCurr<-SmoothingMatrix%*%membcurr





for (dind in 1:4){

DataFolder<-paste0(ScriptLocation, "/Sim_Yshaped/",c("Y_el4x2_Elrand_symm_d50_ver0/","Y_el4x4_Elrand_symm_d50_ver0/","Y_el4x8_Elrand_symm_d50_ver0/", "Y_el4x16_Elrand_symm_d50_ver0/"))[dind]


Indeces<-numeric()






for(snr in 1:5){
SNR<-c(0,64,16,4,1)[snr] #c(0,100,50,20,10,5)[snr]
ErrmatrSNR<-array(0,c(length(R.all),length(lambda.all)))
for(Rl in 1:length(R.all)){
  for(Ll in 1:length(lambda.all)){
  Rread<-R.all[Rl]
  Lread<-lambda.all[Ll]


if (snr==1){C.calc<-as.matrix(read.table(paste0(DataFolder,outname1,"/skCSDall_M512_R",Rread,"lambda",Lread))) #4*pi*0.5*
} 
if (snr!=1){
C.calc<-as.matrix(read.table(paste0(DataFolder,outname1,"/ksCSD_Matr512_R",Rread,"lambda",Lread,"SNR",SNR)))
}
SameplaceAll<-readLines(paste0(DataFolder,outname1,"/sameplace.txt"))
skCSD.all<-array(0,c(seg.nb, 561))
    for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<-C.calc[sameplace,]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums(C.calc[c(sameplace),]),nrow=1)
    }
#skCSD<-skCSD.all[ToPlot,]

ErrmatrSNR[Rl,Ll]<-L1Error(skCSD.all,SmoothedCurr)[[1]]


  }
}
write.table(ErrmatrSNR, paste(DataFolder,"/EL1Error","_SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)

ErrorTableNoise[dind,snr]<-min(ErrmatrSNR)
Indeces<-which(ErrmatrSNR==min(ErrmatrSNR),arr.ind=TRUE)
RTableNoise[dind,snr]<-R.all[Indeces[1]]
LTableNoise[dind,snr]<-lambda.all[Indeces[2]]

cat(snr)
} #snr

}

write.table(ErrorTableNoise, paste("/home/csdori/ksCSD_2014/trunk/Sim_Yshaped","/EL1Error_Noise",sep=""),row.names = FALSE,col.names = FALSE)
write.table(RTableNoise, paste("/home/csdori/ksCSD_2014/trunk/Sim_Yshaped","/EL1Error_RTableNoise",sep=""),row.names = FALSE,col.names = FALSE)
write.table(LTableNoise, paste("/home/csdori/ksCSD_2014/trunk/Sim_Yshaped","/EL1Error_LTableNoise",sep=""),row.names = FALSE,col.names = FALSE)

png(paste0('/home/csdori/ksCSD_2014/trunk/plots/summary/','Y_EffectN2.png'), width=600, height=600)
par(cex=1.3)
matplot(t(ErrorTableNoise),t='b',lwd=2,pch=20,xaxt='n',xlab='SNR', ylab='L1 Error', main='Influence of Noise on the skCSD reconstruction')
axis(1,at=1:5,c("No Noise", "64","16","4","1"))
legend("bottomright", c("2x4","4x4","4x8","4x16"),title="Electrodes",col=1:4,pch=20)
dev.off()


###############################
########xx How do we get from Fig 5 to Fig 6
##################################################x
png(paste0('/home/csdori/ksCSD_2014/trunk/plots/summary/','Representation2.png'), ,width=6, height=7, res=400, unit="in" )
 layout(matrix(c(1,1,2), 1, 3, byrow = TRUE))
par(cex=1)
plot(1,type='n',xlim=c(-150,150),ylim=c(-50,550),xlab='x (um)',ylab='y (um)',main='Branching Morphology \n Representation')

 #plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main=paste("SNR=",SNR),asp=1)
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4])#,col = "gray")
ToPlot<-seq(1,86,3) 


text(segmid[c(1:dimseg)[ToPlot],xA],segmid[c(1:dimseg)[ToPlot],yA], label=c(1:86)[ToPlot])

plot(1,type='n',xlim=c(-150,150),ylim=c(-50,550),xlab='Time',ylab='Segment ID',main='Interval \n Representation',xaxt='n',yaxt='n')
text(0,seq(-50,550,,86)[ToPlot], label=c(1:86)[ToPlot])
dev.off()




