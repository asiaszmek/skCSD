
cellTypes <-c(Ballstick=1, Y_shaped=2,Morpho1=3, Agasbogas=4, Mainen=5,User_Defined=6, Gang_Simple=7, Domi=8)
##electrode orientation
##the chosen coordinate is parallel to the normal vector of the plane of the electrode grid, this will be the 'X', the other follows based on the rigth hand rule
electrodeOrientation <-c(x=1, y=2, z=3)
electrodeDistribute<-c(Grid=1, Random=2, Hexagonal=3, Domi=4)
##which python code to ren
LFPysim<-c(random=1, Y_symmetric=2, Mainen=3, Oscill = 4, Const=5 )

LFPy_setup <- function(cellname,celltype,path,lfpysim,eldistribute,orientation,colnb,rownb,cellelectrodedist,xmin,xmax,ymin,ymax,ssNb,triside,sigma) {
    
    
     
    defaultdir <- getwd()
    main.folder <- defaultdir
########################
    ##this is pre-set, but the value doesn't matter much
    ##defining the position of the grid
    TriSide<-triside
 
    SetSeedNb<-ssNb
    set.seed(SetSeedNb)
    ##Creating the grid
    eleccoords<-array(0,c(rownb*colnb,3))
    eleccoords[,2]<-rep(seq(xmin,xmax,length.out=rownb),colnb)
    eleccoords[,3]<-rep(seq(ymin,ymax,length.out=colnb),each=rownb)
    eleccoords[,1]<-rep(cellelectrodedist,rownb*colnb)
    
    if (eldistribute==1){
        if (orientation==1){
            eleccoords[,2]<-rep(seq(xmin,xmax,length.out=rownb),colnb)
            eleccoords[,3]<-rep(seq(ymin,ymax,length.out=colnb),each=rownb)
            eleccoords[,1]<-rep(cellelectrodedist,rownb*colnb)
        }
        if (orientation==2){
            eleccoords[,3]<-rep(seq(xmin,xmax,length.out=rownb),colnb)
            eleccoords[,1]<-rep(seq(ymin,ymax,length.out=colnb),each=rownb)
            eleccoords[,2]<-rep(cellelectrodedist,rownb*colnb)
        }
        if (electrodeOrientation[(orientation)]==3){
            eleccoords[,1]<-rep(seq(xmin,xmax,length.out=rownb),colnb)
            eleccoords[,2]<-rep(seq(ymin,ymax,length.out=colnb),each=rownb)
            eleccoords[,3]<-rep(cellelectrodedist,rownb*colnb)
        }
    }
    if (eldistribute==2){
        if (orientation==1){
            
            eleccoords[,2]<-runif(rownb*colnb,xmin,xmax)
            eleccoords[,3]<-runif(rownb*colnb,ymin,ymax)
            eleccoords[,1]<-rep(cellelectrodedist,rownb*colnb)
        }
        if (orientation==2){
            eleccoords[,3]<-runif(rownb*colnb,xmin,xmax)
            eleccoords[,1]<-runif(rownb*colnb,ymin,ymax)
            eleccoords[,2]<-rep(cellelectrodedist,rownb*colnb)
        }
        if (orientation==3){
            eleccoords[,1]<-runif(rownb*colnb,xmin,xmax)
            eleccoords[,2]<-runif(rownb*colnb,ymin,ymax)
            eleccoords[,3]<-rep(cellelectrodedist,rownb*colnb)
        }
    }
    
    if (eldistribute==3){
        
        Xmin<-xmin
        Ymin<-ymin
        Zdist<-cellelectrodedist #Cell to electrode Distance
        ColNb<-colnb#number of rows, Even numb
        RowNb<-rownb/2 #number of columns, Even numb
        TriHeight<-TriSide*cos(pi/6)
        TriX1<-Xmin + TriSide*(1:(ColNb))-TriSide
        TriX2<-Xmin -TriSide/2 + TriSide*(1:(ColNb))
        TriY1<-Ymin + 2*TriHeight*(1:(RowNb)) - TriHeight
        TriY2<-Ymin + 2*TriHeight*(1:(RowNb))
        
        grid1<-expand.grid(TriX1,TriY1)
        grid2<-expand.grid(TriX2,TriY2)
        
        Xcord<-c(grid1[[1]] ,grid2[[1]])
        Ycord<-c(grid1[[2]] ,grid2[[2]])
        ##plot(Xcord,Ycord,asp=1)
        ##elccord<-c(Xcord, Ycord, rep(0,length(Ycord)))
        

        if (orientation==1){
            eleccoords[,2]<-Xcord
            eleccoords[,3]<-Ycord
            eleccoords[,1]<-rep(cellelectrodedist,length(Xcord))
        }
        if (orientation==2){
            eleccoords[,3]<-Xcord
            eleccoords[,1]<-Ycord
            eleccoords[,2]<-rep(cellelectrodedist,length(Xcord))
            }
        if (orientation==3){
            eleccoords[,1]<-Xcord
            eleccoords[,2]<-Ycord
            eleccoords[,3]<-rep(cellelectrodedist,length(Xcord))
        }
    }
        
    
    if (eldistribute==4){
        eleccoords<-as.matrix(read.table(paste0(main.folder,'/simulation/ElcoordsDomi14.txt')))
    }
    if (celltype==1){
        morpho.location<-paste(main.folder,'/simulation/morphology/ballstick.swc',sep='')
    }
    if (celltype==2){
        morpho.location<-paste(main.folder,'/simulation/morphology/villa.swc',sep='')
    }
    if (celltype==3){
        morpho.location<-paste(main.folder,'/simulation/morphology/morpho1.swc',sep='')
    }
    
    if (celltype==4){
        morpho.location<-paste(main.folder,'/simulation/morphology/neuron_agasbogas.swc',sep='')
    }
    if (celltype==5){
        morpho.location<-paste(main.folder,'/simulation/morphology/Mainen_swcLike.swc',sep='')
    }
    ##You can choose an
    if (celltype==6){
        morpho.location<-paste(main.folder,'/simulation/morphology/retina_ganglion.swc',sep='')
        ##morpho.location<-file.choose() }#.swc files work only!!
    }
    
    if (celltype==7){
        morpho.location<-paste(main.folder,'/simulation/morphology/Badea2011Fig2Du.CNG.swc',sep='')
        
    }
    if (celltype==8){
        morpho.location<-paste(main.folder,'/simulation/morphology/DomiCell.swc',sep='')
    }
        
###############################
    if (lfpysim==1){
        simName<-'LFP_calc.py'#_osc_sine.py'
    }
    if (lfpysim==2){
        simName<-'LFP_Y_symmetric.py' #'LFP_calc_active.py' 
    }
    if (lfpysim==3){
        simName<-'LFPymod_example6.py'
    }
    if (lfpysim==4){
        simName<-'LFP_calc_sine.py'#_osc_sine.py'
    }
    if (lfpysim==5){
        simName<-'LFP_calc_constInj.py'#_osc_sine.py'
    }
#####################################
    ##active channel description
    active.location<-paste(main.folder,'/simulation/morphology/active.hoc',sep='')
##############################
    
    simulation.location<-paste(main.folder,"/",path,sep='')
    cell.name<-cellname
    cat(cell.name)
    dir.create(paste(simulation.location,'/',cell.name,sep=''))
##############x
    
    
    ##plot(alak[,3],alak[,4])
    alak<-as.matrix(read.table(morpho.location, comment.char="#"))
    ##cat(alak)
    limx<-range(alak[,3],eleccoords[,1])
    limy<-range(alak[,4],eleccoords[,2])
    limz<-range(alak[,5],eleccoords[,3])
    ##plot(alak[,3],alak[,4])
    cat(paste('X range:', c(range(alak[,3]),'\n')),fill=TRUE)
    cat(paste('Y range:', c(range(alak[,4]),'\n')))
    cat(paste('Z range:', c(range(alak[,5]),'\n')))
    setupplot.name<-paste(simulation.location,'/',cell.name,'/setupplot.png',sep='')
    ## png(setupplot.name)
    ## sc<-scatterplot3d(alak[,3],alak[,4],alak[,5],pch=20,color='RED', xlim=limx, ylim=limy, zlim=limz,main='The cell and the electrode',xlab='x',ylab='y',zlab='z', scale.y=1)
    ## sc$points3d(eleccoords[,1],eleccoords[,2],eleccoords[,3],col='BLACK',pch=15)
    ## dev.off()
    ## sc<-scatterplot3d(alak[,3],alak[,4],alak[,5],pch=20,color='RED', xlim=limx, ylim=limy, zlim=limz,main='The cell and the electrode',xlab='x',ylab='y',zlab='z', scale.y=1)
    ## sc$points3d(eleccoords[,1],eleccoords[,2],eleccoords[,3],col='BLACK',pch=15)
    
    ## plotrange<-range(c(limx, limy,limz))
    ## plot3d(alak[,3],alak[,4],alak[,5],
    ##        col='RED',xlim=plotrange,
    ##        ylim=plotrange, zlim=plotrange,
    ##        aspect=TRUE, xlab="x (um)",
    ##        ylab="y (um)", "z (um)",
    ##        main="Simulational Setup")
    ## plot3d(eleccoords[,1],eleccoords[,2],eleccoords[,3],col='BLACK',add=TRUE, size=10)
    ## aspect3d(1,1,1)
    
    
    elcoord.location.name<-paste(simulation.location,'/',cell.name,'/elcoord_x_y_z',sep='')
    
    write.table(c(eleccoords),elcoord.location.name,col.names=FALSE, row.names=FALSE)
    write.table(c(cellelectrodedist, sigma),paste(simulation.location,'/',cell.name,'/elprop',sep=''),col.names=FALSE, row.names=FALSE)
    write.table(morpho.location,paste(simulation.location,'/',cell.name,'/morphology.txt',sep=''),col.names=FALSE, row.names=FALSE,quote=FALSE)
    write.table(active.location,paste(simulation.location,'/',cell.name,'/active.txt',sep=''),col.names=FALSE, row.names=FALSE,quote=FALSE)
    write.table(electrodeOrientation[orientation],paste(simulation.location,'/',cell.name,'/ElecOrient.txt',sep=''),col.names=FALSE, row.names=FALSE,quote=FALSE)
###########
    write.table(cell.name,paste(simulation.location,'/cellname',sep=''),col.names=FALSE, row.names=FALSE,quote=FALSE)
    ##return(list(elec=eleccoords,alak))
    cat(simName)
    return(list(simName1=simName,simulation.location=simulation.location,cellname=cellname,celltype=celltype))
}


LFPy_run<-function(cellname,celltype,path,lfpysim,eldistribute,orientation,colnb,rownb,cellelectrodedist=50,xmin=-100,xmax=600,ymin=0,ymax=200,ssNb=123456,triside=19) {
    ##simulation.location<-paste(svalue(main.folder.n),svalue(simulation.location.name),sep="")
    ##running the LFP simulation
        
    ##where to save the simulation results
    
    
    ##checking runtime
    ##Rprof(paste(outputfilename,'/profile_LFPsim.out',sep=''))
    sigma <-0.3
    main.folder.n <- getwd()
    simName2Be <- LFPy_setup(cellname,celltype,path,lfpysim,eldistribute,orientation,colnb,rownb,cellelectrodedist,xmin,xmax,ymin,ymax,ssNb,triside,sigma)
    simulation.location<-simName2Be$simulation.location
    cell.name.name <-cellname
    outputfilename<-paste(simulation.location,'/',cell.name.name,sep='')
    cat(outputfilename)
    setwd(simulation.location)
    cat(paste("\n", getwd()))
    cat(simName2Be$simName1)
    ##system(paste0('cd ', svalue(main.folder.n), '/simulation'))
    ##cat(paste('cd', svalue(main.folder.n)))
    ##export PYTHONPATH=/usr/local/nrn/lib/python
    system(paste0('python ', simName2Be$simName1))
    ##cat('This simulation is testing the Y shaped neuron only!!! ')
    ##system('ipython LFP_Y_symmetric.py')
    
    cat('LFPy simulation ready! ')
    
    setwd(outputfilename)
    ##Writing out the KCSD details
    if(celltype == 1){
        source(paste0(main.folder.n,"/utils/sCSDFun.R"))
        ##system(paste0('ipython ', simulation.location,'/', 'dori1DkCSD.py'))
    }## else system(paste0('ipython ', simulation.location,'/', 'KCSD_Chat/kCSD_dori.py'))
    ## 
    ##Calculate the connections
    
    ##Reading in the segment information from LFPy
    seg.cord<- matrix(as.matrix(read.table('coordsmid_x_y_z')),ncol=3)
    ##seg.kord<-seg.cord
    seg.start<-round(matrix(as.numeric(as.matrix(read.table('coordsstart_x_y_z'))),ncol=3),3)
    seg.end<-round(matrix(as.numeric(as.matrix(read.table('coordsend_x_y_z'))),ncol=3),3)
    seg.diam<-as.numeric(as.matrix(read.table('segdiam_x_y_z')))
    seg.db<-length(seg.diam)
    connections.1<-numeric()
    connections.2<-numeric()
    for (i in 1:(seg.db)){
        ##connected<-which(seg.start[i,1]==seg.end[,1] & seg.start[i,3]==seg.end[,3] & seg.start[i,2]==seg.end[,2],arr.ind=TRUE)
        connected<-which((seg.start[i,1]==seg.start[,1] & seg.start[i,3]==seg.start[,3] & seg.start[i,2]==seg.start[,2]) | (seg.start[i,1]==seg.end[,1] & seg.start[i,3]==seg.end[,3] & seg.start[i,2]==seg.end[,2]))
        
        melyik<-rep(i,length(connected))
        mihez<-connected
        if(connected[1]!=melyik[1]){
            connections.1<-c(connections.1,melyik )
            connections.2<-c(connections.2,mihez )
            ##cat(paste( i, 'is connected to ', connected,'\n'))
        }
    }
    ##calculating the connections from the coordinates of the segments
    conn.matr<-array(0,c(length(connections.1),2))
    conn.matr[,1]<-connections.1
    conn.matr[,2]<-connections.2
    conn.matr<-conn.matr[-which(conn.matr[,1]==conn.matr[,2]),]
    ##writing out the connections between the segments
    write.table(conn.matr, file="connections.txt",append=FALSE)
    ##write the coordinates in a different format
    write.table(seg.cord, file="segcoordinates.txt",append=FALSE)
    cat('Connection estimations ready!')
    
    
    SummaryofSims<-paste0(simulation.location,"/SummaryofSimulations.txt")
    What2Write2File<-paste(cell.name.name, celltype, lfpysim, cellelectrodedist, orientation, xmin, xmax, ymin, ymax, colnb, rownb, sigma, eldistribute,ssNb, "\n")
    if(file.exists(SummaryofSims)==FALSE)
        file.create(SummaryofSims)
    cat(What2Write2File,file=SummaryofSims,append=TRUE)
    
    
    setwd(main.folder.n)
    cat('Setup done')
    ##Rprof(NULL)
}

## celltype <- 1
## path <- 'simulation'
## lfpysim <- 4
## eldistribute <- 1
## orientation <-2
## colnb <- 1
## rownb <-8


## cellname <- "BS_"
## startpath = getwd()
## LFPy_run(cellname,celltype,path,lfpysim,eldistribute,orientation,colnb,rownb,cellelectrodedist=50,xmin=-100,xmax=600,ymin=0,ymax=200,ssNb=123456,triside=19)
