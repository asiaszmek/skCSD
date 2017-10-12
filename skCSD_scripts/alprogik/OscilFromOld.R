###################
#bas
library(data.table)
BaseFile<-"/home/csdori/ksCSD_2014/trunk/simulation/GanglionOsc_256Regular/"
TargetFile<-"/home/csdori/ksCSD_2014/trunk/simulation/GanglionOsc_128Regular/"#Berd40_21x21/"
membcurr<-as.matrix(fread(paste0(BaseFile,"membcurr")))
examplelfp<-as.matrix(fread(paste0(BaseFile,"myLFP")))
segcoord<-as.matrix(read.table(paste0(BaseFile,"segcoordinates.txt")))

elcoord<-matrix(as.matrix(read.table(paste0(TargetFile,"elcoord_x_y_z"))),ncol=3)



sCSDTCalc<-function(Elcord,CellCord,const){
  
  d1<-dim(Elcord)[1]
d2<-dim(CellCord)[1]
  TransMatr<-array(0, c(d1,d2))
  for(i in 1:d2){
    for(j in 1:d1){
      TransMatr[j,i]<-1/const*1/sqrt(sum((Elcord[j,]-CellCord[i,])^2))
    }
  }
  return(TransMatr)
}



TransMatrix<-sCSDTCalc(elcoord,segcoord,4*pi*0.5)
#sCSD<-solve(TransMatrix)%*%LFP
NewLFP<-TransMatrix%*%membcurr

write.table(NewLFP,paste0(TargetFile,"myLFP"),row.names = FALSE, col.names = FALSE)






file.copy(paste0(BaseFile,"membcurr"), paste0(TargetFile,"membcurr"), overwrite = TRUE)
file.copy(paste0(BaseFile,"somav.txt"), paste0(TargetFile,"somav.txt"), overwrite = TRUE)
file.copy(paste0(BaseFile,"time"), paste0(TargetFile,"time"), overwrite = TRUE)


