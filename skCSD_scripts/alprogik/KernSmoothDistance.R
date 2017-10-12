#source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))

KernSmoothDistance<-function(BandWidth,inname, outname, outname1, Name){
membcurr<-as.matrix(read.table(paste(outname,"/membcurr",sep="")))
seg.length<-as.matrix(read.table(paste(inname,"seglength",sep="")))
funaramvonal<-function(x) x/seg.length
memb.currents.line<-array(0, c(dim(membcurr)))
memb.currents.line<-apply(membcurr,2,funaramvonal)
morpho<-as.matrix(read.table(paste0(inname,'segcoordinates.txt')))
DistMorpho<-as.matrix(dist(morpho))
#image(as.matrix(DistMorpho))
SmoothedCurr<-array(0,dim(memb.currents.line))


SmoothingMatrix<-array(0,c(length(seg.length),length(seg.length)))

for(i in 1:length(seg.length)) SmoothingMatrix[i,]<-exp(-DistMorpho[i,]^2/2/BandWidth^2)
SmoothingMatrix<-SmoothingMatrix/colSums(SmoothingMatrix)
SmoothedCurr<-SmoothingMatrix%*%memb.currents.line

write.table(SmoothedCurr, paste0(inname,outname1,"/",Name))
return(SmoothedCurr)
}
#############################

KernSmoothDistanceOther<-function(BandWidth,inname, outname, outname1, Name,skCSD.all){



morpho<-as.matrix(read.table(paste0(inname,'segcoordinates.txt')))
DistMorpho<-as.matrix(dist(morpho))
#image(as.matrix(DistMorpho))


SmoothedCurr<-array(0,dim(skCSD.all))


SmoothingMatrix<-array(0,c(dim(SmoothedCurr)[1],dim(SmoothedCurr)[1]))

for(i in 1:length(dim(SmoothedCurr)[1])) SmoothingMatrix[i,]<-exp(-DistMorpho[i,]^2/2/BandWidth^2)
SmoothingMatrix<-SmoothingMatrix/colSums(SmoothingMatrix)
SmoothedCurr<-SmoothingMatrix%*%skCSD.all

write.table(SmoothedCurr, paste0(inname,outname1,"/",Name))
return(SmoothedCurr)
}


#####################
#For spike triggered data


KernSmoothDistanceSpTrig<-function(BandWidth,inname, outname, outname1, Name){
membcurr<-as.matrix(fread(paste(outname,"/membSpTrig",sep="")))
seg.length<-as.matrix(read.table(paste(inname,"seglength",sep="")))
funaramvonal<-function(x) x/seg.length
memb.currents.line<-array(0, c(dim(membcurr)))
memb.currents.line<-apply(membcurr,2,funaramvonal)
morpho<-as.matrix(read.table(paste0(inname,'segcoordinates.txt')))
DistMorpho<-as.matrix(dist(morpho))
#image(as.matrix(DistMorpho))
SmoothedCurr<-array(0,dim(memb.currents.line))


SmoothingMatrix<-array(0,c(length(seg.length),length(seg.length)))

for(i in 1:length(seg.length)) SmoothingMatrix[i,]<-exp(-DistMorpho[i,]^2/2/BandWidth^2)
SmoothingMatrix<-SmoothingMatrix/colSums(SmoothingMatrix)
SmoothedCurr<-SmoothingMatrix%*%memb.currents.line

write.table(SmoothedCurr, paste0(inname,outname1,"/",Name))
return(SmoothedCurr)
}

