source('runLFP.R') #load functions
source('runKernel.R')

celltype <- 1 #Ball and stick
abspath <- 'simulation' 
lfpysim <- 6 #simulate oscillations
eldistribute <- 1 #grid
orientation <-2 #y-axis
cellelectrodedist <- 50
colnb <- 1 #one column

dir = getwd()
for (rownb in c(8,16,128)){
    
    cellname <- paste("BS_d",cellelectrodedist,"el",rownb,"_CosChanging",sep="")

    LFPy_run(cellname,celltype,abspath,lfpysim,eldistribute,orientation,colnb,rownb,cellelectrodedist,xmin=-100,xmax=600,ymin=0,ymax=200,ssNb=123456,triside=19)
    cat('Run LFP')
    setwd(dir)
    path = paste('simulation/',cellname,sep="")
    #source('simulation/BS_cos_simulation.R')
    cat('simulation')
    kernel_calculate(path)
    setwd(dir)
}

source('Figures/Fig1_BallStickOsc/BSCosPlot.R')
