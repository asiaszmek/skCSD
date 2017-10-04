source('runLFP.R') #load functions
celltype <- 1 #Ball and stick
path <- 'simulation' 
lfpysim <- 4 #simulate oscillations
eldistribute <- 1 #grid
orientation <-2 #y-axis
cellelectrodedist <- 50
colnb <- 1 #one column

for (rownb in c(8,16,128)){
    
    cellname <- paste("BS_d",cellelectrodedist,"el",rownb,"_CosChanging")

    LFPy_run(cellname,celltype,path,lfpysim,eldistribute,orientation,colnb,rownb,cellelectrodedist,xmin=-100,xmax=600,ymin=0,ymax=200,ssNb=123456,triside=19)
}
