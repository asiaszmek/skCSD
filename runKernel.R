library('foreach')
library('doMC')
library('fields')
library('MASS')

kernel_calculate <- function(simulation.location.name,
                             basis.width.min=3, #R min
                             basis.width.max=7, #R max
                             basis.width.step=1, #R step
                             basis.reg.min=-5, #lambda min
                             basis.reg.max=-1, #lambda max
                             basis.reg.step=1, #lambda step
                             basis.numb.g = 512 ) #Basis number
{
    wherearewenow<-getwd()

    ##where to save the outputs
    where2save<-simulation.location.name
    dir.create(where2save)
    outnameSub<-"skCSDreconstruct/"
    
    outname<-paste(where2save,"/",outnameSub,sep="")
    dir.create(outname)
    ##Parameters
    ##The method calculates the solutions in case of different parameters (width and number of Gaussian functions) it is possible to set the minimum and maximum value of these parameter and the step sizes
    M<-basis.numb.g
        
    DirData<-where2save
    if (getwd() != DirData)
      setwd(DirData)

    ##reading in electrode coordinates
    
    elec.kord <- as.matrix(read.table("elcoord_x_y_z"))
    elec.kord<-matrix(elec.kord,ncol=3)
    
    ##number of electrodes
    ##LFP
    ##lfp<-read.table('myLFP',colClasses='numeric' )
    LFP<-as.matrix(read.table("myLFP" ))
    lfpmatr<-dim(LFP)
    LFP<-matrix(as.numeric(LFP),nrow=lfpmatr[1])
    
    ##Don't forget to write out also these parameters
    ##m(lfpmatr)
    ##morphology
    morpho <-read.table("segcoordinates.txt")
    ##connection information
    connections <- read.table("connections.txt")
    ##cat(connections)
    ##membrane currents
    memb.currents <- as.matrix(read.table("membcurr"))
    cmatr<-dim(memb.currents )
    memb.currents <-matrix(as.numeric(memb.currents),nrow=cmatr[1])
    ##rm(lfpmatr)
    ##length of segments
    seg.length <- as.numeric(as.matrix(read.table("seglength")))
    
    rangeX<-range(elec.kord[,1],morpho[,1])
    rangeY<-range(elec.kord[,2],morpho[,2])
    rangeZ<-range(elec.kord[,3],morpho[,3])
#################################
    setwd(wherearewenow)
    
    source('utils/kernel_basisfunctions_Regularizal.R',local=TRUE,verbose=TRUE)
    sCSD_currents<-ksCSD_all( basis.width.min, basis.width.max, basis.width.step, basis.reg.min, 
                             basis.reg.max, basis.reg.step, LFP, elec.kord,
                             memb.currents, seg.length, where2save, M)
}  
