source('runLFP.R')
source('runKernel.R')
source('BS_cos_simulation.R')  

celltype <- 1 #Ball and stick
abspath <- 'simulation' 
lfpysim <- 6 #simulate oscillations
eldistribute <- 1 #grid
orientation <-2 #y-axis
cellelectrodedist <- 50
colnb <- 1 #one column

dir <- getwd()
scripts_dir <- paste(dir,'skCSD_scripts/alprogik',sep='/')
cat(scripts_dir)
for (rownb in c(8,16,128)){
    
    cellname <- paste("BS_d",cellelectrodedist,"el",rownb,"_CosChanging",sep="")
     #load functions
    LFPy_run(cellname,celltype,abspath,lfpysim,eldistribute,orientation,colnb,rownb,cellelectrodedist,xmin=-100,xmax=600,ymin=0,ymax=200,ssNb=123456,triside=19)
    setwd(dir)
    path = paste('simulation/',cellname,sep="")
    
    simulate(path)
    
    setwd(dir)
    kernel_calculate(path)
    SNR <- 0
    source(paste(scripts_dir,'JustErrorCalcFunReg.R',sep="/"))
    setwd(dir)
    JustErrorCalcFunReg(path,'skCSDreconstruct','',0,scripts_dir)
    setwd(dir)
}

source('Figures/Fig1_BallStickOsc/BSCosPlot.R')
