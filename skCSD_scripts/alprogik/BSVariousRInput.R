#Generate simple source width different  widths for ball&stick

R<-c()



DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/BS_d50_el128/"

for(R in seq(10,200,20)){
NewFolder<-paste0("/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/BS_d50_el128_Gauss",R,"/")
dir.create(NewFolder)
FilesInFolder<-list.files(DataFolder)

file.copy(paste0(DataFolder,FilesInFolder),NewFolder)



segmid<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsmid_x_y_z'), colClasses='numeric')),ncol=3)
seg.nb<-dim(segmid)[1]
elec<-matrix(as.matrix(read.table(paste0(DataFolder,'elcoord_x_y_z'), colClasses='numeric')),ncol=3)

MembCurr<-1/(R*sqrt(2*pi))*exp(-(segmid[,3]-segmid[27,3])^2/R^2)

plot(segmid[,3],MembCurr)


Dmatr<-proxy::dist(segmid,elec)
Tmatr<-1/(4*pi*0.5)*1/Dmatr
LFP<-t(Tmatr)%*%MembCurr



myLFP<-cbind(LFP,LFP,LFP,LFP,LFP)
myMembCurr<-cbind(MembCurr,MembCurr,MembCurr,MembCurr,MembCurr)

write.table(myLFP, paste0(NewFolder,"myLFP"), row.names = FALSE, col.names = FALSE )
write.table(myMembCurr, paste0(NewFolder,"membcurr"), row.names = FALSE, col.names = FALSE )
}

#####################################################x
###################### for Y-shaped
############################################################



R<-c()



DataFolder<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/BS_d50_el128/"


  NewFolder<-paste0("/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/BS_d50_el128_Gauss",R,"/")
  dir.create(NewFolder)
  FilesInFolder<-list.files(DataFolder)
  
  file.copy(paste0(DataFolder,FilesInFolder),NewFolder)
  
  
  
  segmid<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsmid_x_y_z'), colClasses='numeric')),ncol=3)
  seg.nb<-dim(segmid)[1]
  elec<-matrix(as.matrix(read.table(paste0(DataFolder,'elcoord_x_y_z'), colClasses='numeric')),ncol=3)

    
  MembCurr<-1/(R*sqrt(2*pi))*exp(-(segmid[,3]-segmid[27,3])^2/R^2)
  
  plot(segmid[,3],MembCurr)
  
  
  Dmatr<-proxy::dist(segmid,elec)
  Tmatr<-1/(4*pi*0.5)*1/Dmatr
  LFP<-t(Tmatr)%*%MembCurr
  
  
  
  myLFP<-cbind(LFP,LFP,LFP,LFP,LFP)
  myMembCurr<-cbind(MembCurr,MembCurr,MembCurr,MembCurr,MembCurr)
  
  write.table(myLFP, paste0(NewFolder,"myLFP"), row.names = FALSE, col.names = FALSE )
  write.table(myMembCurr, paste0(NewFolder,"membcurr"), row.names = FALSE, col.names = FALSE )
}



#####################################
#Generate simple source width different Domi

R<-c()



DataFolder<-"/home/zoe/agy/ksCSD_2014/trunk/simulation/Domi14/"
setwd(DataFolder)
#morpho<-as.matrix(read.table('segcoordinates.txt'))
morpho<-matrix(as.numeric(as.matrix(read.table('segcoordinates.txt'))),ncol=3)
plot(morpho[,1], morpho[,2], col="WHITE")#,xlim=c(-500,-200))
text(morpho[,1], morpho[,2],paste(1:length(morpho[,1])),cex=0.7)
text(morpho[149:180,1], morpho[149:180,2],paste(149:180),cex=0.7,col="RED")
#  membcurr<-matrix(as.numeric(as.matrix(read.table('simulation/cell_1/membcurr'))))
connections<-as.matrix(read.table('connections.txt')  )

  NewFolder<-paste0("/home/zoe/agy/ksCSD_2014/trunk/simulation/cell_Cos3/")
  dir.create(NewFolder)
  FilesInFolder<-list.files(DataFolder)
  
  file.copy(paste0(DataFolder,FilesInFolder),NewFolder)
  
  
  
  segmid<-matrix(as.matrix(read.table(paste0(DataFolder,'coordsmid_x_y_z'))),ncol=3)
  seglength<-as.matrix(read.table(paste0(DataFolder,'seglength')))
  SpFrNb<-5
  MembMatr<-array(0,c(length(seglength),SpFrNb))
  #fr<-1
  
  #Csak adott reszre nem 0 aramokat rakni
  #for(fr  in 1:50){
  selStim<-149:180
  for(fr  in 1:SpFrNb){
    
  plot(cumsum(seglength[selStim]),cos(cumsum(seglength[selStim])/max(cumsum(seglength[selStim]))*2*pi*fr))
  
  MembCline<-cos(cumsum(seglength[selStim])/max(cumsum(seglength[selStim]))*pi*fr+pi/2) #*2
  
  Membcurr<-MembCline#/seglength[selStim]
  MembMatr[selStim,fr]<-Membcurr
  }
  image(t(MembMatr),col = rainbow(100)) 
  matplot(MembMatr,t="l",xlim=c(149,180))
  
  seg.nb<-dim(segmid)[1]
  elec<-matrix(as.matrix(read.table(paste0(DataFolder,'elcoord_x_y_z'), colClasses='numeric')),ncol=3)
  
  
  Dmatr<-array(0,  c(length(seglength),dim(elec)[1]))
  for(i1 in 1:length(seglength)){
    for(i2 in 1:dim(elec)[1]){
      Dmatr[i1,i2]<-sqrt(sum((segmid[i1,]-elec[i2,])^2))
    }
    
  }
  
  
  Tmatr<-1/(4*pi*0.3)*1/Dmatr
  LFP<-t(Tmatr)%*%MembMatr
  
  
  
  
  write.table(LFP, paste0(NewFolder,"myLFP"), row.names = FALSE, col.names = FALSE )
  write.table(MembMatr, paste0(NewFolder,"membcurr"), row.names = FALSE, col.names = FALSE )
  write.table(1:50, paste0(NewFolder,"time"), row.names = FALSE, col.names = FALSE )
  
# 

