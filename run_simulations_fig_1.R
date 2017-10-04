library(parallel)
runOneBS<-function(rownb){
    source('runLFP.R') #load functions
    source('runKernel.R')
    
    celltype <- 1 #Ball and stick
    path <- 'simulation' 
    lfpysim <- 4 #simulate oscillations
    eldistribute <- 1 #grid
    orientation <-2 #y-axis
    cellelectrodedist <- 50
    colnb <- 1 #one column
    
    dir = getwd()
    
    cellname <- paste("BS_d",cellelectrodedist,"el",rownb,"_CosChanging",sep="")

    LFPy_run(cellname,celltype,path,lfpysim,eldistribute,orientation,colnb,rownb,cellelectrodedist,xmin=-100,xmax=600,ymin=0,ymax=200,ssNb=123456,triside=19)
    path = paste('simulation/',cellname,sep="")
    kernel_calculate(path)
    setwd(dir)

}



# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)
ele = c(8,16,128)
parLapply(cl,ele,runOneBS)


