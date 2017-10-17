
library(data.table)
library(RColorBrewer)
library(fields)




M<-512
lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" ) #c(0.1,0.01,0.001)#
R.all<-2^(3:7) #c(8,16,128)#

dir = getwd()

ScriptLocation<-"skCSD_scripts/alprogik"
cat(paste0(ScriptLocation,"/Colors_BlueRed.R","\n"))
source(paste0(ScriptLocation,"/Colors_BlueRed.R"))



SimPath<-"simulation/" #paste0(ScriptLocation,"simulation/")
ElNumb<-c(8,16, 128) #16???
#ElNumb<-c(10,20,50,100)
locationsData<-c("BS_d50el8_CosChanging/",  "BS_d50el16_CosChanging/","BS_d50el128_CosChanging/")#paste0("BS_d50el",ElNumb) 
#locationsKernel<-c(locationsKernel,paste0("Y_el",ElNumb,"_m200_600"))
locationsKernel<-list()
outname1<-"skCSDreconstruct"
i <- 1
for(loc in locationsData){
    locationsKernel[[i]] <- paste(loc,outname1,sep="")
    i = i+1
}
#outname1<-"/kernelOut_poli_CP3"
#lapply(paste0(SimPath,locationsData, outname1), dir.create)


source(paste0(ScriptLocation,"/errorCalcFun.R"))



figfname <-"Fig_1.png"

png(paste0(figfname),width=10, height=9, res=400, unit="in" )
cat(paste('Figure saved as ',getwd(),"/",figfname,"\n",sep=''))


plot.new()
par(mfrow=c(5,1),oma=c(4,5,4,10))#,mar=c(1,1,1,1))
##par(mfrow=c(length(lambda.all)+2,length(R.all)),oma=c(2,3,4,2))#,mar=c(1,1,1,1))

##layout(rbind(c(1,2,3,4,5)),widths=c(0.5,0.2,1,3,0.2),respect=FALSE) 
##mtext(c("oma[SOUTH<-1]", "oma[WEST<-2]", "oma[NORTH<-3]", "oma[EAST<-4]"),     c(SOUTH<-1, WEST<-2, NORTH<-3, EAST<-4),      line=0.4, cex=1.2, col="green", outer=TRUE)
##title(main="skCSD Reconstruction", outer=TRUE, cex=3)




L1errors<-numeric()


MinError <- numeric()
MinErrorR<-numeric()
MinErrorL<-numeric()




for(Dlength in 1:length(locationsKernel)){
    inname<-paste0(SimPath, locationsKernel[Dlength],"/")
    outname<-paste0(SimPath, locationsData[Dlength])
    OrigCurr<-as.matrix(read.table(paste(outname,"membcurr",sep="")))
    SameplaceAll<-readLines(paste0(inname,"sameplace.txt"))
    coord<-as.matrix(read.table((paste0(outname,"coordsmid_x_y_z"))))[106:156]
    seg.nb<-length(as.vector(read.table(paste(outname,'segdiam_x_y_z',sep=''))))
    cat(coord)
    seg.length<-as.matrix(read.table(paste(outname,"seglength",sep="")))
    funaramvonal<-function(x) x/seg.length
    memb.currents.line1<-array(0, c(dim(OrigCurr)))
    memb.currents.line1<-apply(OrigCurr,2,funaramvonal) 
    ##memb.currents.line1[which(memb.currents.line1<(0.5*min(memb.currents.line1)))]<-0.5*min(memb.currents.line1)
    times2Plot<-seq(1,25,by=3)#1:10 #as.matrix(read.table(paste0(inname,"time")))
  
    coordsMid<-matrix(as.matrix(read.table(paste0(outname,"coordsmid_x_y_z"))),ncol=3)
  
    ## BesteEmatr<-as.matrix(read.table(paste(inname,"EL1ErrorVarR_SNR0",sep="")))
    ## CurrMin<- min(BesteEmatr[26:30,])
    ## IndecesMin<-which(BesteEmatr==CurrMin,arr.ind=TRUE)
    ## MinErrorR<-c(MinErrorR,R.all[IndecesMin[1]-25])
    ## MinErrorL<-c(MinErrorL,lambda.all[IndecesMin[2]])
  
  
    ## R2Plot<- MinErrorR[Dlength]
    ## lambda2Plot<- MinErrorL[Dlength]
  
  
  
    ## ##Select just few Rs and Ms for plotting
    R.allplot<-seq(5,130,25)
    M<-512
    lambda.allplot<-c("1e-05","1e-04","0.001","0.01","0.1" )
  
  
  
    if(Dlength==1){
        par(mfg=c(1,1))
        par(mar=c(1,0.1,1,3))
        
        SegNb<-dim(memb.currents.line1)[1]
        
        ##matplot(sort(coord),memb.currents.line1[c(2,1,3:52),times2Plot],t="l",ylim=c(-0.13,0.13),xaxt="n",main="Ground Truth",ylab="CSD (nA/um)")
    
        col2<-ColoursDori(memb.currents.line1*1.1)[[1]]
        ExtVal<-range(memb.currents.line1*1.1)[2] #ColoursDori(memb.currents.line1)[[2]]
        par(mar=c(1,0.1,0.1,0.1))
        SegNb<-dim(memb.currents.line1)[1]
        cat(length(times2Plot/2))
        cat(length(coord))
        image(times2Plot/2, sort(coord),t(memb.currents.line1[c(2,1,4:52),times2Plot]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", ylab="z (um)")
        mtext("Ground Truth",3)
        mtext("z (um)",side=2,line=2.5)
  #  axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
    
    
    
    }
    ##par(mar=c(0,8,0,8)) 
    ##if(R2Plot<4) par(mar=c(0,8,0,8)) 
    ##if(R2Plot>3) par(mar=c(0,4,0,4)) 
  
  
    par(mfg=c(Dlength+1,1))
  
    skCSD.all<-array(0,c(seg.nb, length(times2Plot)))
    ##for(m in 1:length(lambda.all)){
    ## for(r in 1:length(R.all)){
    m<-1
    r<-1
    Lambda<-lambda.allplot[m]
    R<-R.allplot[r]
    ##Reading in the LFP and adding noise
    ##Ktilda_M512_R55lambda1e-04
    currName<-paste0(outname,outname1,"/skCSDall_M",M,"_R",R,"lambda",Lambda)
    cat(paste(currName,'\n'))
    if(file.exists(currName)) {
        skCSD.all.part<-as.matrix(read.table(currName)) 
    
    
        for(i in 1: seg.nb){
            sameplace<-numeric() 
            sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
            ##Plot the current related ti the different segments
            ##  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
            if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,times2Plot]
            if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),times2Plot]),nrow=1)
        }      
        
    
    
    }else next
  
  
      
    L1errors<-c(L1errors,L1Error(skCSD.all, memb.currents.line1[,times2Plot])$lerrorInTime)
    cat(paste(L1errors))
    par(mar=c(1,0.1,1,0.1))
    if(Dlength!=3) #matplot(coordsMid[c(2,1,3:52),3],skCSD.all[c(2,1,3:52),times2Plot],t="l",ylim=c(-0.13,0.13),xaxt="n",main=paste(ElNumb[Dlength], "Electrodes"),ylab="CSD (nA/um)")
    {
        ##col2<-ColoursDori(memb.currents.line1)[[1]]
        ##ExtVal<-0.11 #ColoursDori(memb.currents.line1)[[2]]
        par(mar=c(1,0.1,0.1,0.1))
        SegNb<-dim(memb.currents.line1)[1]
        image(times2Plot, sort(coord),t(skCSD.all[c(2,1,3:52),]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", ylab="z (um)")
        mtext(paste(ElNumb[Dlength], "Electrodes"),3)
        mtext("z (um)",side=2,line=2.5)
        ## axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
    
    
    }
    if(Dlength==3) #matplot(coordsMid[c(2,1,3:52),3],skCSD.all[c(2,1,3:52),times2Plot],t="l",ylim=c(-0.13,0.13),xlab="x (um)",main=paste(ElNumb[Dlength], "Electrodes"),ylab="CSD (nA/um)")
    {
      
        image(times2Plot, sort(coord),t(skCSD.all[c(2,1,3:52),]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", ylab="z (um)")
        mtext(paste(ElNumb[Dlength], "Electrodes"),3)
        mtext("z (um)",side=2,line=2.5)
        image.plot(times2Plot, sort(coord),t(skCSD.all[c(2,1,3:52),]),col=col2,zlim=c(-ExtVal,ExtVal),xaxt="n", ylab="z (um)",                   legen.only=TRUE,add = TRUE, legend.line = -2   )
    
    
        mtext("z (um)",side=2,line=2.5)
        ##axis(1,round(seq(0,max(times2Plot)-10,,5)),las=2)
        
    }
    
}

## L1errors<-matrix(L1errors,ncol=Dlength)

## matplot(times2Plot/2, L1errors,t="l",ylab="L1 Error",xlab="Spatial frequency",lwd=2)
## legend("topleft",paste(ElNumb),col=1:3, lty=1:3)
## ## #dev.off()
## ## #}


dev.off()















