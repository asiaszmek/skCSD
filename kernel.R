library('scatterplot3d')
library('foreach')
library('doMC')
require(gWidgets)
library('fields')
library('MASS')

##options("guiToolkit"="RGtk2")


kernel_calc <- function() {
    ##Get the location of the current working directory
    wherearewenow<-getwd()
    ##for calculating the kernel functions
    kernel_calculate <- function(h,...) {
        ##where to save the outputs
        where2save<-svalue(simulation.location.name)
        dir.create(where2save)
    
    
        outnameSub<-"skCSDreconstruct/"
    
        outname<-paste(where2save,"/",outnameSub,sep="")
        dir.create(outname)
        ##Parameters
        ##The method calculates the solutions in case of different parameters (width and number of Gaussian functions) it is possible to set the minimum and maximum value of these parameter and the step sizes
        basis.width.min<-svalue(basis.width.min.g)
        basis.width.max<-svalue(basis.width.max.g)
        basis.width.step<-svalue(basis.width.step.g)
        
        basis.reg.min<-svalue(basis.reg.min.g)
        basis.reg.max<-svalue(basis.reg.max.g)
        basis.reg.step<-svalue(basis.reg.step.g)
  
    
        M<-svalue(basis.numb.g)
        
        DirData<-where2save
        setwd(DirData)
    
        ##reading in electrode coordinates
    
        elec.kord <- as.matrix(read.table("elcoord_x_y_z"))
        elec.kord<-matrix(elec.kord,ncol=3)
    
        outname<-paste(where2save,"/",outnameSub,sep="")
        ##if(file.exists(outname)) file.remove(list.files("outname"))
        dir.create(outnameSub)
    
    
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
                                 basis.reg.max, basis.reg.step, LFP , elec.kord ,memb.currents,seg.length  ,where2save,M) #,R.V ,source.V.t)
    }  
  
    win <- gwindow("Calculation of Kernel Functions")
    gp <- ggroup(horizontal=FALSE, cont=win)
    
    ##
    ##PLotting
    ##ggraphics(ps=6, container=win)
    ##     hist(rnorm(100))
  
  
    ##Setting the parameters
    tmp <- gframe("Width of basis functions R=2^k", container=gp, expand=TRUE)
  
  
    glabel("R min:", container=tmp)
    basis.width.min.g<-gedit("3", container=tmp, coerce.with=as.numeric)
  
    glabel("R max:", container=tmp)
    basis.width.max.g<-gedit("7", container=tmp, coerce.with=as.numeric)
  
    glabel("R step:", container=tmp)
    basis.width.step.g<-gedit("1", container=tmp, coerce.with=as.numeric)
  
    ##Lambda regularization parameter
    ##M=2^k....
    tmp <- gframe("Lambda regularization parameter 10^k", container=gp, expand=TRUE)
    glabel("lambda min:", container=tmp)
    basis.reg.min.g<-gedit("-5", container=tmp, coerce.with=as.numeric)
    
    glabel("lambda max:", container=tmp)
    basis.reg.max.g<-gedit("-1", container=tmp, coerce.with=as.numeric)
    
    glabel("lambda step:", container=tmp)
    basis.reg.step.g<-gedit("1", container=tmp, coerce.with=as.numeric)


    tmp <- gframe("Basis number", container=gp, expand=TRUE)
    glabel("M:", container=tmp)
    basis.numb.g<-gedit("512", container=tmp, coerce.with=as.numeric)
  
  
  
###############################################
###########################################
    ##LOading files
    glabel("Where is the data for skCSD?", container=gp)
    simulation.location.name<-gfilebrowse(paste0(wherearewenow, '/simulation/gang_9x9_200/'), 
                                        type='selectdir', quote=FALSE, cont=gp)
  
  
  
  
  
    ##run simulation
    startSimulation<-gbutton("Start simulation", container=gp,
                           handler =kernel_calculate )
  
  
  
    f_exit <-gbutton("Close simulation", container=gp,
                   handler = function( h,...) dispose( win ))
  
    
    setwd(wherearewenow)
  
  
}
############
kernel_calc() 
