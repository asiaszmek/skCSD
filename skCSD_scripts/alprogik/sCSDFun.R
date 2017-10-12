#sCSDFun, ver ysimilar to sCSDCalc, but it doesn't run over several folders....
#Elcord and Cellcord should be arrays with 3 columns with resect to the coordinates
library(MASS)
#CurrDir<-getwd()
#DataFold<-paste0(SimPath,locationsData)
#for(t1 in 1: length(DataFold)){
#setwd(DataFold[t1])
sCSDTCalc<-function(Elcord,CellCord,const){
  if(dim(Elcord)[1]!=dim(CellCord)[1]) {cat( "Elcord and CellCord have different dimensions")
  break}
  DimMatr<-dim(Elcord)[1]
  TransMatr<-array(0, c(DimMatr,DimMatr))
  for(i in 1:DimMatr){
    for(j in 1:DimMatr){
      TransMatr[j,i]<-1/const*1/sqrt(sum((Elcord[j,]-CellCord[i,])^2))
    }
  }
  return(TransMatr)
}


Pstart<-matrix(as.numeric(as.matrix(read.table('coordsstart_x_y_z'))),ncol=3)
Pend<-matrix(as.numeric(as.matrix(read.table('coordsend_x_y_z'))),ncol=3)
Pmid<-matrix(as.numeric(as.matrix(read.table('coordsmid_x_y_z'))),ncol=3)

 elec<-matrix(as.matrix(read.table('elcoord_x_y_z', colClasses='numeric')),ncol=3)
membcurr<-as.matrix(read.table("membcurr"))
seg.length<-as.matrix(read.table("seglength"))
funaramvonal<-function(x) x/seg.length
memb.currents.line<-array(0, c(dim(membcurr)))
memb.currents.line<-apply(membcurr,2,funaramvonal) 
LFP<-as.matrix(read.table("myLFP"))
#SegCord<-as.matrix(read.table("segcoordinates.txt"))
SegNumbers<-1:length(seg.length)

#SegsChoosen<-sort(sample(SegNumbers, min(dim(elec)[1],dim(SegCord)[1]),replace=FALSE))
#SegsChoosen<-seq(range(Pmid[,3])[1],range(Pmid[,3])[2],length.out=dim(elec)[1])
#SegCordChoosen<-SegCord[SegsChoosen,]
BallsStick<-"yes"
if(BallsStick=="yes"){
SegCordChoosen<-array(0,c(dim(elec)[1],3))
SegCordChoosen[,3]<-seq(range(Pmid[,3])[1],range(Pmid[,3])[2],length.out=dim(elec)[1])
} else {
  #I know that in case of the Y shaped, it is standing towards ööö...
  SegCordChoosen<-array(0,c(dim(elec)[1],3))
  SegCordChoosen[,3]<-seq(range(Pmid[,3])[1],range(Pmid[,3])[2],length.out=dim(elec)[1])
  
}
write.table(SegCordChoosen,"sCSDCoord")
TransMatrix<-sCSDTCalc(elec,SegCordChoosen,4*pi*0.5)
#sCSD<-solve(TransMatrix)%*%LFP
sCSDginv<-ginv(TransMatrix)%*%LFP
write.table(sCSDginv,"sCSDginv")
#}

#setwd(CurrDir)
#par(mfrow=c(3,1))
#image(t(membcurr))
#image(t(sCSD))
#image(t(sCSDginv))


