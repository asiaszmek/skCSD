source("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/alprogik/Colors_BlueRed.R")
source("/home/csdori/ksCSD_2014/trunk/alprogik/Colors_BlueRed.R")

#image scale function for plotting color bars for image plots





xA<-1
yA<-2

Dlength<-1
Title<-'Ganglion'
#Title<-'Ball&Stick'
#where2save<-"/media/BA0ED4600ED416EB/agy/ksCSD_SVN/trunk/plots/morphologies.png"
#where2saveName<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/plots/Y4x2_ver1.png"

#par(mfcol=c(1,3),cex=1.2,lwd=2)
#setwd(paste(mainDir ,Morpholocations[j],sep=""))
#setwd("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/2015Y_8x8Eliminate/")
#setwd("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/2015Y_8x8Eliminate_random/")
# setwd("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_el4x8_Elrand_symm_d50_ver1/")
#setwd("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/BS_10_d50/")
#setwd("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/Sim_Yshaped/Y_el4x2_Elrand_random_d50_ver0/")
#setwd("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/GangNew_11x11/")
setwd("/home/csdori/ksCSD_2014/trunk/simulation/")
CellType<-"gangNew"
if(CellType=="gangNew"){
  SimPath<-"/home/csdori/ksCSD_2014/trunk/simulation/"
  #ElNumb<-c(5,7,10,20)
  #ElNumb<-c(10,20,50,100)
  locationsKernel<-c("cell_rand50_gangl_eld50","cell_rand100_gangl_eld50")#,"_m100_600")#_d",100)
  ElNumb<-c(50,100)
  #locationsKernel<-c(locationsKernel,paste0("Y_el",ElNumb,"_m200_600"))
  locationsData<-locationsKernel
  outname1<-"/kernelOut_poli_CP2"
}

Dlength<-1 #length(locationsKernel)
outname<-paste0(SimPath, locationsData[Dlength])
setwd(outname)

El2Ignore<-numeric()
segstart<-matrix(as.matrix(read.table('coordsstart_x_y_z', colClasses='numeric')),ncol=3)
segend<-matrix(as.matrix(read.table('coordsend_x_y_z', colClasses='numeric')),ncol=3)
elec<-matrix(as.matrix(read.table('elcoord_x_y_z', colClasses='numeric')),ncol=3)
SimTime<-as.matrix(read.table(paste0(outname,'/time')))
TimeInstant<-which(SimTime==45) #45
TimeInstant<-TimeInstant+2


inname<-paste0(SimPath, locationsKernel[Dlength],"/")
outname<-paste0(SimPath, locationsData[Dlength])

SimTime<-as.matrix(read.table(paste0(outname,'/time')))

#if((file.exists(paste0(outname,outname1,"/skCSDall_M",BestEMatr[2,Dlength],"_R",
 #                      BestEMatr[1,Dlength],"_SNR0_lambda0"))==FALSE )) 
params<-as.matrix(read.table(paste(inname,outname1,"/params",sep="")))
#skCSDall_M256_R5_SNR0_lambda0
#skCSD<-as.matrix(read.table(paste0(outname,outname1,"/skCSDall_M",256,"_R",55,"_SNR0_lambda0")))
skCSD<-as.matrix(read.table(paste0(outname,outname1,"/skCSD_plain")))

coordsEnd<-as.matrix(read.table(paste0(outname,"/coordsend_x_y_z")))


seglength<-as.matrix(read.table(paste0(outname,'/seglength')))
funaramvonal<-function(x) x/seglength

membCurr<-as.matrix(read.table(paste0(outname,'/membcurr')))
membcurr<-apply(membCurr,2,funaramvonal)




xA<-1
yA<-2
#Plotratio<-Xrange/max(segwidth)
segwidth1<-as.matrix(read.table('segdiam_x_y_z', colClasses='numeric'))
segwidth<-segwidth1/max(segwidth1)*30 # not realistic?
#membcurr<-membcurrO
Not2Plot<-c(1:10,70)
#membcurrO<-membcurr
membcurr<-membcurr[-Not2Plot,] 
for(Tcounter in 1:250){
  TimeInstant<-1+Tcounter*2
plotname<-paste(paste("CSD_MorphoCSD", "-", 
            formatC(Tcounter, digits = 3, flag = "0"), sep = ""), "png", sep = ".")

png(plotname, height=500, width=800, pointsize=20)

par(mfrow=c(1,2))
dimseg<-dim(segstart)[1]
Xrange<-range(segstart[,xA],elec[,xA])
Yrange<-range(segstart[,yA],elec[,yA])


#Grand Truth
#TimeInstant<-TimeInstant+1
#######################x
plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main="CSD")

rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
col1<-ColoursDori(membcurr[,TimeInstant])[[1]]
ExtVal<-ColoursDori(membcurr)[[2]]#[,TimeInstant])[[2]]
#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))

szinskala2<-color.scale(c(membcurr[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
segments(segstart[c(1:dimseg)[-Not2Plot],xA], segstart[c(1:dimseg)[-Not2Plot],yA], segend[c(1:dimseg)[-Not2Plot],xA], segend[c(1:dimseg)[-Not2Plot],yA], col= szinskala2,lwd=segwidth[c(1:dimseg)[-Not2Plot]])
#segments(segstart[1:dimseg,xA], segstart[1:dimseg,yA], segend[1:dimseg,xA], segend[1:dimseg,yA], col= szinskala2,lwd=segwidth[1:dimseg])
#par(mar=c(0.5,0.1,1,0.5)) 
ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col1)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))







#Reconstruction


#TimeInstant<-TimeInstant+1
plot(1,type='n',xlim=Xrange,ylim=Yrange,xlab='x (um)',ylab='y (um)',main="ksCSD")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
col1<-ColoursDori(skCSD[,TimeInstant])[[1]]
ExtVal<-ColoursDori(skCSD)[[2]]#[,TimeInstant])[[2]]
#col.limits<-max(abs(range(skCSD, membrane.curr.line)))#,memb.currents.vonal[,min.indexes])))

szinskala2<-color.scale(c(skCSD[,TimeInstant]),col=col1,zlim=c(-ExtVal*1.01,ExtVal*1.01))
segments(segstart[1:dimseg,xA], segstart[1:dimseg,yA], segend[1:dimseg,xA], segend[1:dimseg,yA], col= szinskala2,lwd=segwidth[1:dimseg])
#par(mar=c(0.5,0.1,1,0.5)) 
ColMatr<-as.matrix(c(1:256))
col2<-ColoursDori(ColMatr)[[1]]
add.image(max(Xrange),min(Yrange)-25, t(ColMatr), adj.x = 0.0, adj.y = 0.0,  image.width = 1, image.height = 0.05, col = col1)
#axis(4,at=c(min(Yrange),max(Yrange)),labels=c("sink","source"))
axis(4,at=c(min(Yrange),max(Yrange)),labels=c(round(-ExtVal*1.01,4),round(ExtVal*1.01,4)))
#points(elec[,xA],elec[,yA])
#if(file.exists(elignorename))  points(elec[ El2Ignore$x,xA],elec[ El2Ignore$x,yA],pch=4,col="RED")

dev.off()
}
