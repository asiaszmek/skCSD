#source("alprogik/ksCSD_compare_error.R")

library('scatterplot3d')
library('foreach')
library('doMC') 
library('fields')
library('MASS')

##############################
#SimPath<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/simulation/"
#SumPlotPlace<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/plots/summary/"
#SimPath<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation"/

#SimPath<-"/media/zoe/PUBLIC/DataCopy/"


#ScriptLocation<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/"
ScriptLocation<-"/home/csdori/ksCSD_2014/trunk/"
SumPlotPlace<-paste0(ScriptLocation,"plots/summary/")
source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))
	dir.create(paste0(ScriptLocation,"plots/"))
	dir.create(SumPlotPlace)

#outname<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_dense/"

locationsKernel<-numeric()
locationsData<-numeric()

CellType<-"Domi" #"gangNew"#"GangBN"#"YVer0"#"Y"# "gangNew" #Yreg"#"BS" #"Yreg" #"BSOLD"
SNR<-0
lambda<-0


#Domi
if(CellType=="Domi"){
  SimPath<-paste0(ScriptLocation,"simulation/")
  ElNumb<-14#15#c(128, 256)
  #locationsKernel<-c("Ganglion_128Regular","Ganglion_256Regular") 
  locationsKernel<-"Domi14"#c("Ganglion_128Close","Ganglion_256Close","Ganglion_128Regular","Ganglion_256Regular","Ganglion_128Random","Ganglion_256Random") 
  
  locationsData<-"Domi14"# "cell_Cos3"#"Domi14" #"cell_Cos2" #locationsKernel
  outname1<-"/kernelOut_4" ##_kCSD"

}


#BS
if(CellType=="BSOLD"){
  SimPath<-paste0(ScriptLocation,"simulation/")
ElNumb<-c(8,16,32,64, 128)
#ElNumb<-c(10,20,50,100)
locationsKernel<-c("BS_d50_el8", "BS_d50_el16","BS_d50_el32","BS_d50_el64", "BS_d50_el128")#paste0("BS_d50_el",ElNumb) 
#locationsKernel<-c(locationsKernel,paste0("Y_el",ElNumb,"_m200_600"))
locationsData<-locationsKernel
outname1<-"/kernelOut_4"
#outname1<-"/kernelOut_poli_CP3"
#lapply(paste0(SimPath,locationsData, outname1), dir.create)

}


if(CellType=="BS"){
  SimPath<-paste0(ScriptLocation,"simulation/")
ElNumb<-c(8,16,32,64, 128)
#ElNumb<-c(10,20,50,100)
locationsKernel<-c("BS_d50_el8_GaussChanging190", "BS_d50_el16_GaussChanging190","BS_d50_el32_GaussChanging190","BS_d50_el64_GaussChanging190", "BS_d50_el128_GaussChanging190")#paste0("BS_d50_el",ElNumb) 
#locationsKernel<-c(locationsKernel,paste0("Y_el",ElNumb,"_m200_600"))
locationsData<-c("BS_d50_el8_GaussChanging190", "BS_d50_el16_GaussChanging190","BS_d50_el32_GaussChanging190","BS_d50_el64_GaussChanging190", "BS_d50_el128_GaussChanging190")#paste0("BS_d50_el",ElNumb,"_Gauss10")  #locationsKernel
outname1<-"/kernelOut_4"
#outname1<-"/kernelOut_poli_CP3"
#lapply(paste0(SimPath,locationsData, outname1), dir.create)

}

if(CellType=="Y"){
  #SimPath<-paste0(ScriptLocation,"simulation/")
  #SimPath<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/Sim_Yshaped/"
   SimPath<-paste0(ScriptLocation,"Sim_Yshaped/")

    #SimPath<-paste0("/home/zoe/Dropbox/skCSD/data/")
  
  Versions<-c(0:5)
Yelectrodes<-rep(c(2,4,8,16),each=length(Versions))
#Yelectrodes<-16 #rep(c(16),each=length(Versions))
ElNumb<-4*Yelectrodes
#locationsKernel<-paste0("Y_el4x",Yelectrodes,"_Rot_symm_d50_ver",Versions) 
locationsKernel<-paste0("Y_el4x",Yelectrodes,"_Elrand_symm_d50_ver",Versions)
#locationsKernel<-"Y_el4x16_Elrand_symm_d50_ver0_kiserlet" #"Mainen3D_ElGrid2"
#locationsData<-paste0("Y_el4x",Yelectrodes,"_Elrand_symm_d50_ver",Versions)
locationsData<-locationsKernel
outname1<-"/kernelOut_4"
}



if(CellType=="YVer0"){
  #SimPath<-paste0(ScriptLocation,"simulation/")
  #SimPath<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/Sim_Yshaped/"
   SimPath<-paste0(ScriptLocation,"Sim_Yshaped/")

    #SimPath<-paste0("/home/zoe/Dropbox/skCSD/data/")
  
  Versions<-c(0)
Yelectrodes<-rep(c(2,4,8,16),each=length(Versions))
#Yelectrodes<-16 #rep(c(16),each=length(Versions))
ElNumb<-4*Yelectrodes
#locationsKernel<-paste0("Y_el4x",Yelectrodes,"_Rot_symm_d50_ver",Versions) 
locationsKernel<-paste0("Y_el4x",Yelectrodes,"_Elrand_symm_d50_ver",Versions)
#locationsKernel<-"Y_el4x16_Elrand_symm_d50_ver0_kiserlet" #"Mainen3D_ElGrid2"
#locationsData<-paste0("Y_el4x",Yelectrodes,"_Elrand_symm_d50_ver",Versions)
locationsData<-locationsKernel
outname1<-"/kernelOut_4"
}


if(CellType=="Yreg"){
#Y_el4x4_RotS_symm_d50_ver0
   SimPath<-SimPath<-paste0(ScriptLocation,"simulation/") 
  Versions<-c(0)
Yelectrodes<-rep(c(2,4,8,16),each=length(Versions))
ElNumb<-4*Yelectrodes
locationsKernel<-paste0("Y_el4x",Yelectrodes,"_RotS_symm_d50_ver",Versions)
locationsData<- paste0("Y_el4x",Yelectrodes,"_RotS_symm_d50_ver",Versions)#locationsKernel
outname1<-"/kernelOut_4"
}

if(CellType=="gang"){
 SimPath<-paste0(ScriptLocation,"simulation/")
  #ElNumb<-c(5,7,10,20)
  #ElNumb<-c(10,20,50,100)
  locationsKernel<-c("gang1","gang2","gang3")#,"_m100_600")#_d",100)
  ElNumb<-c(100,400,48)
  #locationsKernel<-c(locationsKernel,paste0("Y_el",ElNumb,"_m200_600"))
  locationsData<-locationsKernel
  outname1<-"/kernelOut_poli_kCSD"
}

if(CellType=="gangNew"){
  SimPath<-paste0(ScriptLocation,"simulation/")
  ElNumb<-128#c(128, 256)
  #locationsKernel<-c("Ganglion_128Regular","Ganglion_256Regular") 
  locationsKernel<-"GanglionOsc_128Regular"#c("Ganglion_128Close","Ganglion_256Close","Ganglion_128Regular","Ganglion_256Regular","Ganglion_128Random","Ganglion_256Random") 
  
  #ElNumb<-c(5,7,10,20)
  #ElNumb<-c(10,20,50,100)
  #locationsKernel<-c("GangNew_11x11","GangNew_16x16","GangNew_21x21")#,"_m100_600")#_d",100)
  #locationsKernel<- c("Hex_8x8_inter70_d30","Hex_6x5_inter100_d30","Hex_12x10_inter50_d15","Hex_12x10_inter50")[4]
  #ElNumb<-c(64,30,120,120)[4] #, 120, 144)
  #locationsKernel<-c(locationsKernel,paste0("Y_el",ElNumb,"_m200_600"))
  locationsData<-locationsKernel
   outname1<-"/kernelOut_4" ##_kCSD"
}

if(CellType=="GangBN"){
  SimPath<-paste0(ScriptLocation,"simulation/")
  ElNumb<-c(25,25,25,49,81,81,25)
  #locationsKernel<-c("Ganglion_128Regular","Ganglion_256Regular") 
  locationsKernel<-c("gang_5x5_50","gang_5x5_100","gang_5x5_200", "gang_7x7_200","gang_9x9_200","gang_9x9_100", "gang_5x5_400")#"Berd40_21x21"#c("gang_9x9_200","gang_9x9_50","gang_9x9_100","gang_5x5_25","gang_5x5_50","gang_5x5_100","gang_5x5_200","gang_5x5_400")#
  #locationsKernel<-c("gang_10x10_50","gang_10x10_100","gang_20x20_25","gang_20x20_50")#
  locationsData<-locationsKernel
   outname1<-"/kernelOut_4" ##_kCSD"
}


if(CellType=="Berdondini"){
  SimPath<-paste0(ScriptLocation,"simulation/")
  ElNumb<-441
  #locationsKernel<-c("Ganglion_128Regular","Ganglion_256Regular") 
  locationsKernel<-"Berd40_21x21"#c("gang_9x9_200","gang_9x9_50","gang_9x9_100","gang_5x5_25","gang_5x5_50","gang_5x5_100","gang_5x5_200","gang_5x5_400")#

  locationsData<-locationsKernel
  outname1<-"/kernelOut_4" ##_kCSD"
}

####Cross validation function
#ksCSD cross-validation
#cross.valid function gives back the value of the crossvalidation error and also the crossvalidated potential
#source("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/alprogik/errorCalcFun.R")
source(paste0(ScriptLocation,"alprogik/Colors_BlueRed.R"))
source(paste0(ScriptLocation,"alprogik/errorCalcFun.R"))
#source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))
 source(paste0(ScriptLocation,"alprogik/JustErrorCalcFun.R"))
#source(paste0(ScriptLocation,"alprogik/PlotFunctions.R"))
source(paste0(ScriptLocation,"alprogik/PlotFunctionsReg.R"))
source(paste0(ScriptLocation,"alprogik/ImageScale.R"))
source(paste0(ScriptLocation,"alprogik/sCSDFun.R"))#Calc.R"))
#source("/home/dcserpan/Documents/skCSD/skCSDnew/trunk/alprogik/errorCalcFun.R")
#source("/home/dcserpan/Documents/skCSD/skCSDnew/trunk/alprogik/PlotFunctions.R")


##########################
#This is how the currents are estimated from the potentials.. current at the location on the middle point od the segments as on the line, so sometimes the occur more.. For multiple occurances look at the sampleace.txt 
#C.calc<-1/const*t(K.tilda)%*%ginv(K)%*%LFP
#t(K.tilda)%*%ginv(K) is also saved as transfermatrix
#So it would be interesting to see, that for example on a asimple line, how look the currnets on "each side" of the line, whoch are otherwise summed up

#Don't forget, that some of the recordings are might not taken itno accout, look for these in El2Ignore
#Tikhonov regularization

#CalcErrorAndPlot<-function(SimPath,locationsKernel,locationsData,SNR, lammbda, withsoma)

for(l1 in 1:0){
# FUNCTION TO CALCULATE THE CURRENTS, eRRORS
#SNR.All<-c(0,10,100,1000,10000)
lammbda<-l1
#for (snr in 1:length(SNR.All)){
CalcErrorAndPlot(SimPath,locationsKernel,locationsData,SNR, 0)#lammbda)
}
##########################################################################
#####################################################################
#FOR PLOTTING THE RESULTS OF THE PREVIOUS FUNVTION... fOR R AN M THIS PLOTS THE RECONSTRUCTED skCSD
#ImageAllParCurr(SimPath,locationsKernel,locationsData,SNR, lammbda,6,6)
#ImageAllParCurrReg(SimPath,locationsKernel,locationsData, 0,5,5)
ImageAllParCurrReg(SimPath,locationsKernel,locationsData, 0,c(8,32,128),c(0.1,1e-5))
}
#only plot the best reconstruction in time


#break
##################################\###########################
###############################
#PLOTTING THE BEST ERROR IN CASE OF EACH ELECRODE
#ElErrorOut<-ElNum_ErrorPlot(SimPath,locationsKernel,locationsData,SNR, lammbda)
alma<-ElNum_ErrorPlot(SimPath,locationsKernel,locationsData,SNR)
BestEMatr<-alma$BestEMatr

matplot(t(matrix(data.matrix(BestEMatr[3,]),ncol=4)),pch=20)
ElNumb<-c(8,16,32,64)
#lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" )
#this plots the best error in case of different versions
png(paste0(SumPlotPlace,"BSelNbvsVar.png"), width=500, height=500)
par(mfrow=c(1,1), cex=1.2)
 plot(1:5,data.matrix(BestEMatr[3,]),pch=4, 
xlab="Number of Electrodes", ylab="CV Error", cex=2, lwd=2,,t="b",xaxt="n")
axis(1, at=c(1:5),labels=c(8,16,32,64, 128))
#legend("topright",c("Grid ",paste("Rand",1:5)),pch=20,col=1:6,title="El Geom")
dev.off()



png(paste0(SumPlotPlace,"YelRotNbvsVar.png"), width=900, height=500)
par(mfrow=c(1,3), cex=1.2)
 matplot(1:4,t(matrix(data.matrix(BestEMatr[3,]),ncol=4)),pch=4, 
xlab="Number of Electrodes", ylab="CV Error", cex=2, lwd=2,,t="b",xaxt="n")
axis(1, at=c(1:4),labels=c(8,16,32,64))
legend("topright",c("Grid ",paste("Rand",1:5)),pch=20,col=1:6,title="El Geom")


matplot(1:4,t(matrix(data.matrix(BestEMatr[1,]),ncol=4)),pch=4, 
xlab="Number of Electrodes", ylab="Best R", cex=2, t="b",lwd=2,xaxt="n")
axis(1, at=c(1:4),labels=c(8,16,32,64))
#legend("topright",c("Grid ",paste("Rand",1:5)),pch=20,col=1:6,title="El Geom")

matplot(1:4,t(matrix(data.matrix(BestEMatr[2,]),ncol=4)),pch=4, 
xlab="Number of Electrodes", ylab="Best lambda", cex=2, t="b",lwd=2, log="y",xaxt="n")
axis(1, at=c(1:4),labels=c(8,16,32,64))
dev.off()



#legennd("topright",lambda.all, col=1:)
BestEMatr<-alma$BestEMatr



BestEMatr<-ElErrorOut$BestEMatr
BestEMatrFix<-ElErrorOut$BestEMatrFix
BestEMatrCV<-ElErrorOut$BesterrorCV
#plot several things at the same time
#Error for same inputs, different SNR
#####################################


#locationsKernel<-numeric()
#locationsData<-numeric()


ElNumb<-c(8,16,32,64)
#locationsKernel<-paste0("Y_el",ElNumb,"_m200_600_d",50)
#locationsKernel<-c(locationsKernel,paste0("Y_el",ElNumb,"_m200_600_d",100))
#locationsData<-paste0("Y_el",ElNumb,"_symm_d",50)
#locationsData<-c(locationsData,paste0("Y_el",ElNumb,"_symm_d",100))[4]
locationsData<-locationsKernel
#outname1<-"/kernelOut"
SNR.All<-c(0,10,100,1000,10000)
SNR.Errors<-numeric()
for (snr in 1:length(SNR.All)){
  SNR<-SNR.All[snr]
  lammbda<-0
  ElErrorOut<-ElNum_ErrorPlot(SimPath,locationsKernel,locationsData,SNR, lammbda)
  BestEMatr<-ElErrorOut$BestEMatr
  BestEMatrFix<-ElErrorOut$BestEMatrFix
  SNR.Errors<-cbind(SNR.Errors,BestEMatr[3,])
}




png(paste0(SumPlotPlace,"ElError_Y_symmetric_d50.png"))
par(mar=c(4,4,4,6))
matplot(ElNumb,SNR.Errors, t="l", main="Y symmetric d=50",xlab="Number of Electrodes",ylab="fix smoothed error")
par(xpd=TRUE)
legend(max(ElNumb)+4,mean(SNR.Errors),c("Clean ",paste(SNR.All[-1])),pch=20,col=1:6,title="SNR")
par(xpd=FALSE)
dev.off()



#######################
#Plotting Kappa


#locationsKernel<-numeric()
#locationsData<-numeric()


ElNumb2<-c(8,16,32,64)
#locationsKernel<-paste0("Y_el",ElNumb,"_m200_600_d",100)
#locationsKernel<-c(locationsKernel,paste0("Y_el",ElNumb,"_m200_600_d",100))
#locationsData<-paste0("Y_el",ElNumb,"_symm_d",100)
#locationsData<-c(locationsData,paste0("Y_el",ElNumb,"_symm_d",100))[4]
#locationsData<-locationsKernel
RangeKappa<-numeric()
png(paste0(SumPlotPlace,"Kappas_Y_d100.png"),height=500, width=1200)
par(mfrow=c(1,length(locationsKernel)))
for(Dlength in 1:length(locationsKernel)){
  
  inname<-paste0(SimPath, locationsKernel[Dlength],"/")
  outname<-paste0(SimPath, locationsData[Dlength])
  params<-as.matrix(read.table(paste(inname,outname1,"/params",sep="")))
  R.all<-unique(params[1,])#[3]
  M.all<-unique(params[2,])#[5]
  KappaMatr<-read.table(paste(inname,outname1,"/kappaK","_SNR", SNR, "_lambda",lambda,sep=""))
  RangeKappa<-c(RangeKappa,range(KappaMatr))
  
  matplot(R.all,KappaMatr,log="y",t="l", xlab="R", ylab="Kappa", 
          main="Kappa value of K",sub=paste("Number of electrodes:", ElNumb2[Dlength-(Dlength%/%(length(ElNumb2)+1))*length(ElNumb2)] ))
  legend("bottomright", c(paste(M.all)),col=1:6, pch=20, title="M" )
}
dev.off()

RangeKappaMatr<-matrix(RangeKappa,nrow=2)
matplot(t(RangeKappaMatr), log="y",t="l")


#####################################
















library('scatterplot3d')
library('foreach')
library('doMC') 
#library('fields')
library('MASS')

##############################
#SimPath<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/simulation/"
#SimPath<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/"

#SumPlotPlace<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/plots/summary/"
#outname<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_dense/"

locationsKernel<-numeric()
locationsData<-numeric()

BestEAll<-numeric()
BestEFixAll<-numeric()


Versions<-c(0:5)
Yelectrodes<-rep(c(2,4,8,16),each=length(Versions))
#ElNumb<-4*unique(Yelectrodes)
#ElNumb<-4*Yelectrodes
locationsKernel<-paste0("Y_el4x",Yelectrodes,"_Elrand_random_d50_ver",Versions)
#locationsKernel<-c(locationsKernel,paste0("Y_el",ElNumb,"_m200_600_d",100))
locationsData<-paste0("Y_el4x",Yelectrodes,"_Elrand_symm_d50_ver",Versions)
#locationsData<-locationsKernel
outname1<-"/kernelOut_poli"

ElErrorOut<-ElNum_ErrorPlot(SimPath,locationsKernel,locationsData,SNR, lammbda)
BestEMatr<-ElErrorOut$BestEMatr # this gives back the results of 
BestEMatrFix<-ElErrorOut$BestEMatrFix



ErrorElectrodes<-matrix(c(BestEMatr[3,]),nrow=6)
ErrorFixElectrodes<-matrix(c(BestEMatrFix[3,]),nrow=6)

matplot(ElNumb2,t(ErrorElectrodes),xlab="El number",ylab="L2Error",
        main="Effect of increasing the number of Electrodes",t="l", lwd=2,ylim=c(0.8,1))
matplot(ElNumb,t(ErrorFixElectrodes),xlab="El number",ylab="L2Error Smoothed",
        main="Effect of increasing the number of Electrodes",t="l", lwd=2)

#SumPlotPlace<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/plots/summary/"
png(paste0(SumPlotPlace,"Y_symm_randomElectrodes_Error.png"))
par(mar=c(4,4,4,6),cex=1.3)
matplot(ElNumb2,t(ErrorElectrodes),xlab="El number",ylab="L2Error",
        main="Effect of Increasing the Number of Electrodes",t="l", lwd=2)
par(xpd=TRUE)
legend(max(ElNumb2)+4,mean(ErrorElectrodes),c("Grid ",paste("Rand",1:5)),pch=20,col=1:6,title="El Geom")
par(xpd=FALSE)
dev.off()
png(paste0(SumPlotPlace,"Y_symm_randomElectrodes_SmoothedError.png"))
par(mar=c(4,4,4,6),cex=1.3)
matplot(ElNumb,t(ErrorFixElectrodes),xlab="El number",ylab="L2Error Smoothed",
        main="Effect of increasing the number of Electrodes",t="l", lwd=2)
par(xpd=TRUE)
legend(max(ElNumb)+4,mean(ErrorFixElectrodes),c("Grid ",paste("Rand",1:5)),pch=20,col=1:6,title="El Geom")
par(xpd=FALSE)
dev.off()

####################################################################
########################################################


library('scatterplot3d')
library('foreach')
library('doMC') 
#library('fields')
library('MASS')

##############################
#SimPath<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/simulation/"
#SimPath<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/"

#SumPlotPlace<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/plots/summary/"
#outname<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_dense/"


#legend(max(ElNumb)+4,mean(ErrorFixElectrodes),c("Grid ",paste("Rand",1:5)),pch=20,col=1:6,title="El Geom")





#Imageplot at different time segments to compare the results
skCSDTime<-numeric()

for(Dlength in 1:length(locationsKernel)){
  # Dlength<-2
  inname<-paste0(SimPath, locationsKernel[Dlength],"/")
  outname<-paste0(SimPath, locationsData[Dlength])
  
  SimTime<-as.matrix(read.table(paste0(outname,'/time')))
  TimeInstant<-which(SimTime==5) #45
  TimeInstant<-TimeInstant+2
  #BestEMatr<-BestEMatrFix
  
if(file.exists(paste0(outname,outname1,"/skCSDall_M",512,"_R",BestEMatr[1,Dlength], "lambda",BestEMatr[2,Dlength]))==FALSE) next
      skCSD.all<-array(0,c(86, 561))
 skCSD.all.part<-as.matrix(read.table(paste0(outname,outname1,"/skCSDall_M",512,"_R",BestEMatr[1,Dlength], "lambda",BestEMatr[2,Dlength])))
    #Reading in tge currents
  cat(Dlength)
 
 
 SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))

 for(i in 1: seg.nb){
      sameplace<-numeric() 
      sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
      #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
      if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
      if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
    }
 
 
skCSDTime<-c(skCSDTime,c(skCSD.all[,TimeInstant])) #/max(c(skCSD[,TimeInstant]
  coordsEnd<-as.matrix(read.table(paste0(outname,"/coordsend_x_y_z")))
  
}

seglength<-as.matrix(read.table(paste0(outname,'/seglength')))
funaramvonal<-function(x) x/seglength

membCurr<-as.matrix(read.table(paste0(outname,'/membcurr')))
membCurr<-apply(membCurr,2,funaramvonal)[,TimeInstant]
skCSDTimeAll<-matrix(skCSDTime,ncol=length(locationsKernel))


png(paste0(SumPlotPlace,"Y_NewCool_time5.png"), width=800, height=800, pointsize=32)
    # clear objects  
#graphics.off()    
par(oma=c(2,2,2,2))  
#par(mfrow=c(1,2))
layout(rbind(c(1,2,3,4,5)),widths=c(0.5,0.2,1.5,3,0.2),respect=FALSE) 

#par(mfrow=c(1,1))
col1<-ColoursDori(membCurr)[[1]]
ExtVal<-ColoursDori(membCurr)[[2]]
par(mar=c(0.5,1,1,0.5)) 
image(0:1,1:86,t(as.matrix(membCurr)),col=col1, zlim=c(-ExtVal,ExtVal),main="CSD",xaxt="n",frame.plot=TRUE)
abline(h=0.5)
abline(v=1)
lines(c(membCurr/ExtVal/2)+0.5,1:86,col="BLACK")
par(mar=c(0.5,0.1,1,0.5)) 
image.scale(t(as.matrix(membCurr)),col=col1,axis.pos=4,zlim=c(-ExtVal,ExtVal))
#layout.show(2)
col2<-ColoursDori(skCSDTimeAll)[[1]]
ExtVal<-ColoursDori(skCSDTimeAll)[[2]]
par(mar=c(0.5,3,1,0.5)) 
image(1:4,1:dim(skCSDTimeAll)[1],t(skCSDTimeAll[,seq(1,24,by=6)]),col=col2, yaxt='n',zlim=c(-ExtVal,ExtVal),xaxt="n",main="ksCSD",frame.plot=TRUE)
abline(h=0.5)
abline(v=4.5)
axis(1, at=c(1:4),labels=c(8,16,32,64),las=2)
for(i in 1:4){
lines(t(skCSDTimeAll[,seq(1,24,by=6)])[i,]/(ExtVal*2)+i,1:86)
}

col2<-ColoursDori(skCSDTimeAll)[[1]]
par(mar=c(0.5,0.5,1,1.5)) 
image(1:dim(t(skCSDTimeAll[,-seq(1,24,by=6)]))[1],1:dim(t(skCSDTimeAll[,-seq(1,24,by=6)]))[2],t(skCSDTimeAll[,-seq(1,24,by=6)]),col=col2,yaxt='n',zlim=c(-ExtVal,ExtVal),xaxt="n",main="ksCSD - Random ",frame.plot=TRUE)
abline(h=0.5)
axis(1, at=c(1:4)*5-2,labels=c(8,16,32,64),las=2)
#,legend.shrink=0.3,legend.width=0.8,legend.mar=3)
abline(v=c(seq(0,dim(skCSDTimeAll)[2],by=5)+0.5),col="BLACK", lty=4)

for(i in 1:4){
  lines(t(rowMeans(skCSDTimeAll[,((i-1)*6+1):(i*6)]))/(ExtVal*2)+5*(i-1)+3,1:86,col="BLACK")
}
par(mar=c(0.5,0.1,1,0.5)) 
image.scale(t(skCSDTimeAll),col=col2,axis.pos=4,zlim=c(-ExtVal,ExtVal))


dev.off()



############################################################
################################################
#Imageplot at different time segments to compare the results
skCSDTime<-numeric()

for(Dlength in 1:length(locationsKernel)){
  # Dlength<-2
  
  inname<-paste0(SimPath, locationsKernel[Dlength],"/")
  outname<-paste0(SimPath, locationsData[Dlength])
  
  SimTime<-as.matrix(read.table(paste0(outname,'/time')))
  TimeInstant<-which(SimTime==5) #45
  TimeInstant<-TimeInstant+2
  #BestEMatr<-BestEMatrFix
  if((file.exists(paste0(outname,outname1,"/skCSDall_M",BestEMatr[2,Dlength],"_R",
                         BestEMatr[1,Dlength],"_SNR0_lambda0"))==FALSE )) 
params<-as.matrix(read.table(paste(inname,outname1,"/params",sep="")))
  #skCSDall_M256_R5_SNR0_lambda0
  skCSD<-as.matrix(read.table(paste0(outname,outname1,"/skCSDall_M",BestEMatr[2,Dlength],"_R",
                                     BestEMatr[1,Dlength],"_SNR0_lambda0")))
  skCSDTime<-c(skCSDTime,c(skCSD[,TimeInstant])) #/max(c(skCSD[,TimeInstant]
  coordsEnd<-as.matrix(read.table(paste0(outname,"/coordsend_x_y_z")))
  
}

seglength<-as.matrix(read.table(paste0(outname,'/seglength')))
funaramvonal<-function(x) x/seglength

membCurr<-as.matrix(read.table(paste0(outname,'/membcurr')))
membCurr<-apply(membCurr,2,funaramvonal)[,TimeInstant]
skCSDTimeAll<-matrix(skCSDTime,ncol=length(locationsKernel))


png(paste0(SumPlotPlace,"Y_Cool_time5.png"), width=800, height=800, pointsize=32)
    # clear objects  
#graphics.off()    
par(oma=c(2,2,2,2))  
#par(mfrow=c(1,2))
layout(rbind(c(1,2,3,4,5)),widths=c(0.5,0.2,1.5,3,0.2),respect=FALSE) 

#par(mfrow=c(1,1))
col1<-ColoursDori(membCurr)[[1]]
ExtVal<-ColoursDori(membCurr)[[2]]
par(mar=c(0.5,1,1,0.5)) 
image(0:1,1:86,t(as.matrix(membCurr)),col=col1, zlim=c(-ExtVal,ExtVal),main="CSD",xaxt="n",frame.plot=TRUE)
abline(h=0.5)
abline(v=1)
lines(c(membCurr/ExtVal/2)+0.5,1:86,col="BLACK")
par(mar=c(0.5,0.1,1,0.5)) 
image.scale(t(as.matrix(membCurr)),col=col1,axis.pos=4,zlim=c(-ExtVal,ExtVal))
#layout.show(2)
col2<-ColoursDori(skCSDTimeAll)[[1]]
ExtVal<-ColoursDori(skCSDTimeAll)[[2]]
par(mar=c(0.5,3,1,0.5)) 
image(1:4,1:dim(skCSDTimeAll)[1],t(skCSDTimeAll[,seq(1,24,by=6)]),col=col2, yaxt='n',zlim=c(-ExtVal,ExtVal),xaxt="n",main="ksCSD",frame.plot=TRUE)
abline(h=0.5)
abline(v=4.5)
axis(1, at=c(1:4),labels=c(8,16,32,64),las=2)
for(i in 1:4){
lines(t(skCSDTimeAll[,seq(1,24,by=6)])[i,]/(ExtVal*2)+i,1:86)
}

col2<-ColoursDori(skCSDTimeAll)[[1]]
par(mar=c(0.5,0.5,1,1.5)) 
image(1:dim(t(skCSDTimeAll[,-seq(1,24,by=6)]))[1],1:dim(t(skCSDTimeAll[,-seq(1,24,by=6)]))[2],t(skCSDTimeAll[,-seq(1,24,by=6)]),col=col2,yaxt='n',zlim=c(-ExtVal,ExtVal),xaxt="n",main="ksCSD - Random ",frame.plot=TRUE)
abline(h=0.5)
axis(1, at=c(1:4)*5-2,labels=c(8,16,32,64),las=2)
#,legend.shrink=0.3,legend.width=0.8,legend.mar=3)
abline(v=c(seq(0,dim(skCSDTimeAll)[2],by=5)+0.5),col="BLACK", lty=4)

for(i in 1:4){
  lines(t(rowMeans(skCSDTimeAll[,((i-1)*6+1):(i*6)]))/(ExtVal*2)+5*(i-1)+3,1:86,col="BLACK")
}
par(mar=c(0.5,0.1,1,0.5)) 
image.scale(t(skCSDTimeAll),col=col2,axis.pos=4,zlim=c(-ExtVal,ExtVal))


dev.off()

plot(membCurr,1:length(membCurr),t="l")
image(1:dim(skCSDTimeAll)[2],1:dim(skCSDTimeAll)[1],t(skCSDTimeAll),col=rainbow(200))
abline(v=c(seq(0,dim(skCSDTimeAll)[2],by=6)+0.5))
image(t(cbind(membCurr,skCSDTimeAll)),col=rainbow(200),zlim=range(cbind(membCurr,skCSDTimeAll)))





#######################################################3
#############################################################
#would be nice to have several lines above eachother representing the various distances...for example...



#WorkDir<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/2015_Morpho/kernelOut"
#Plot the kernel funtions
# 
#setwd(WorkDir)
#Kmatr<-as.matrix(read.table("K_M120_R10"))
#Ktildematr<-as.matrix(read.table("Ktilda_M120_R10"))


matplot(t(Kmatr),t='l',main="K matrix")
#Plotting the K tilde
matplot(t(Ktildematr),t='l', main="K.tilde", xlab="Where to calculate", ylab="value")
image(t(Ktildematr),xlab="Where to calculate",ylab="Connection with the electrode position")

#########################################################
########################################################


#Creates a plot, for 1 point in time compair smoothed ground through with reconstruction for the points where the reconstruction is small and high


mainDir<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/Sim_Yshaped//"
#locations<-c("2015Y_symmetric3x3","2015Y_symmetric4x4","2015Y_symmetric5x5", "2015Y_symmetric3x3closer","2015Y_symmetric4x4closer","2015Y_symmetric5x5closer","2015Y_8x8Eliminate_random","2015Y_8x8Eliminate")
#locations<-c("2015Y_symmetric4x4")

outname<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/Sim_Yshaped/Y_el4x2_Elrand_symm_d50_ver0/"
outname1<-"kernelOut_poli"
setwd(paste0(outname, outname1))
par(mfrow=c(1,1))
#membCurrSmooth<-as.matrix(read.table('membcurr_smoothed',colClasses='numeric'))
#membCurr<-as.matrix(read.table('membcurr_line'),colClasses='numeric'))
seglength<-as.matrix(read.table(paste0(outname,'/seglength')))
funaramvonal<-function(x) x/seglength

membCurr<-as.matrix(read.table(paste0(outname,'/membcurr')))
membCurrL<-funaramvonal(membCurr)
skCSD<-as.matrix(read.table('skCSD_MC',colClasses='numeric'))
ErrorTime<-as.matrix(read.table('BestSetupErrorinTime'))

plot(ErrorTime)
BaseError<-300
MinError<-which.min(ErrorTime)
MaxError<-which.max(ErrorTime)


#cbind(membCurrSmooth[,BaseError],skCSD[,BaseError])
PlotErrorTimeCompare<-function(BaseError,membCurrSmooth, skCSD, membCurr){
  Y1<-range(-abs(c(membCurrSmooth[,BaseError],skCSD[,BaseError])),abs(c(membCurrSmooth[,BaseError],skCSD[,BaseError])))
  par(mar=c(5, 4, 4, 6) + 0.1)
  matplot(cbind(membCurrSmooth[,BaseError],skCSD[,BaseError]), xlab="Compartment number", ylab="Normalized RMSE",main="Comparison",t='l',ylim=Y1)
  legend("topright",c("Simulated Currents","Smoothed Currents","skCSD"),col=c(1,2, "green"), pch=20)
  #legend("top",c("BaseError"))
  #maybe add the original current on the plot
  
  par(new=TRUE)
  
  Y2<-range(-abs(c(membCurr[,BaseError])),abs(c(membCurr[,BaseError])))
  plot( membCurr[,BaseError],   xlab="", ylab="", , 
        axes=FALSE, type="l", col="green",ylim=Y2)
  ## a little farther out (line=4) to make room for labels
  mtext("Original currents",side=4,col="green",line=4) 
  axis(4, ylim=c(0,1) ,col="green",col.axis="green",las=1)
  abline(h=0,col="BLUE")
}
PlotErrorTimeCompare(BaseError,membCurrSmooth, skCSD, membCurr)
PlotErrorTimeCompare(MinError,membCurrSmooth, skCSD, membCurr)
PlotErrorTimeCompare(MaxError,membCurrSmooth, skCSD, membCurr)
## Draw the time axis
#axis(1,pretty(range(time),10))


##################################
#Plot basis functions
where.t<-as.matrix(read.table(paste(outname,outname1,"/where_t",sep="")))
length.fitted.curve<-as.matrix(read.table(paste(outname,outname1,"/SmoothinKernelWidth",sep="")))[2]
j<-1
R<-50
M<-100
source.t<-seq(0,length.fitted.curve-3,length.out=M)
lam<-length.fitted.curve/M

fun1<-function(x) exp(-((source.t[i]-x)^2)/(R^2))
i<-1
Xlim<-c(where.t[j]-3*R,where.t[j]+3*R)

plot(fun1,xlim=Xlim,ylim=c(0,1),main=paste("R:",R,"\n"," Distance btw basises:",round(lam,2 )),xlab="t",yaxt="n",ylab="")
for(i in 2:M){
  plot(fun1,xlim=Xlim,add=TRUE,ylim=c(0,1))
}
abline(v=c(where.t),col=2)




####Not to good version
fun1<-function(x) {
  
  if (x>length.fitted.curve){ 
    x<-x-length.fitted.curve
  }
  else if (x<0){ 
    x<-x+length.fitted.curve
  }
  return(exp(-((source.t[i]-x)^2)/(R^2)))
}
i<-1
Xlim<-c(-2*R,length.fitted.curve+2*R)
plot(Vectorize(fun1),xlim=Xlim,ylim=c(0,1),main=paste("R:",R,"\n"," Distance btw basises:",round(lam,2 )),xlab="t",yaxt="n",ylab="")
for(i in 2:M){
  plot(fun1,xlim=Xlim,add=TRUE,ylim=c(0,1))
}
abline(v=c(where.t),col=2)


#which segments should taken into account at compariso?
ToCompare<-1:3
