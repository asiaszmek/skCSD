##This script is for the simulation of the extracellular potential of a single cell based on a morphology and some other defined parameter.

library(gWidgetsRGtk2)
options(guiToolkit = "RGtk2")
require(gWidgets)

##install.packages("gWidgetsRGtk2", dep=TRUE)
library('scatterplot3d')
library("rgl")
##library('rglwidget')
DefaultDirectory<-getwd()

runLFPy <- function() {
  
    sigma<-0.3 #0.5
    wherearewenow<-getwd()
    ##cell types, with some examples (1,2,3,4, or user definid)
    cellTypes <-c(Ballstick=1, Y_shaped=2,Morpho1=3, Agasbogas=4, Mainen=5,User_Defined=6, Gang_Simple=7, Domi=8)
    ##electrode orientation
    ##the chosen coordinate is parallel to the normal vector of the plane of the electrode grid, this will be the 'X', the other follows based on the rigth hand rule
    electrodeOrientation <-c(x=1, y=2, z=3)
    electrodeDistribute<-c(Grid=1, Random=2, Hexagonal=3, Domi=4)
    ##which python code to ren
    LFPysim<-c(random=1, Y_symmetric=2, Mainen=3, Oscill = 4, Const=5 )
  
#####################xx
    LFPy_setup <- function(h,...) {
        main.folder<-svalue(main.folder.n)
########################
        ##this is pre-set, but the value doesn't matter much
        ##defining the position of the grid
        cell.electrode.dist<-svalue(cellelectrodedist)
        TriSide<-2*svalue(triside)
        x.min<-svalue(xmin)
        x.max<-svalue(xmax)
        y.min<-svalue(ymin)
        y.max<-svalue(ymax)
        col.nb<-svalue(colnb) #how many column in the grid
        row.nb<-svalue(rownb) #how many rows
        SetSeedNb<-svalue(ssNb)
        set.seed(SetSeedNb)
        ##Creating the grid
        elec.coords<-array(0,c(row.nb*col.nb,3))
        elec.coords[,2]<-rep(seq(x.min,x.max,length.out=row.nb),col.nb)
        elec.coords[,3]<-rep(seq(y.min,y.max,length.out=col.nb),each=row.nb)
        elec.coords[,1]<-rep(cell.electrode.dist,row.nb*col.nb)
        
        if (electrodeDistribute[svalue(eldistribute)]==1){
            if (electrodeOrientation[svalue(orientation)]==1){
                elec.coords[,2]<-rep(seq(x.min,x.max,length.out=row.nb),col.nb)
                elec.coords[,3]<-rep(seq(y.min,y.max,length.out=col.nb),each=row.nb)
                elec.coords[,1]<-rep(cell.electrode.dist,row.nb*col.nb)
            }
            if (electrodeOrientation[svalue(orientation)]==2){
                elec.coords[,3]<-rep(seq(x.min,x.max,length.out=row.nb),col.nb)
                elec.coords[,1]<-rep(seq(y.min,y.max,length.out=col.nb),each=row.nb)
                elec.coords[,2]<-rep(cell.electrode.dist,row.nb*col.nb)
            }
            if (electrodeOrientation[svalue(orientation)]==3){
                elec.coords[,1]<-rep(seq(x.min,x.max,length.out=row.nb),col.nb)
                elec.coords[,2]<-rep(seq(y.min,y.max,length.out=col.nb),each=row.nb)
                elec.coords[,3]<-rep(cell.electrode.dist,row.nb*col.nb)
            }
        }
        if (electrodeDistribute[svalue(eldistribute)]==2){
            if (electrodeOrientation[svalue(orientation)]==1){
                
                elec.coords[,2]<-runif(row.nb*col.nb,x.min,x.max)
                elec.coords[,3]<-runif(row.nb*col.nb,y.min,y.max)
                elec.coords[,1]<-rep(cell.electrode.dist,row.nb*col.nb)
            }
            if (electrodeOrientation[svalue(orientation)]==2){
                elec.coords[,3]<-runif(row.nb*col.nb,x.min,x.max)
                elec.coords[,1]<-runif(row.nb*col.nb,y.min,y.max)
                elec.coords[,2]<-rep(cell.electrode.dist,row.nb*col.nb)
            }
            if (electrodeOrientation[svalue(orientation)]==3){
                elec.coords[,1]<-runif(row.nb*col.nb,x.min,x.max)
                elec.coords[,2]<-runif(row.nb*col.nb,y.min,y.max)
                elec.coords[,3]<-rep(cell.electrode.dist,row.nb*col.nb)
            }
        }
        
        if (electrodeDistribute[svalue(eldistribute)]==3){
            
            Xmin<-x.min
            Ymin<-y.min
            Zdist<-cell.electrode.dist #Cell to electrode Distance
            ColNb<-col.nb#number of rows, Even numb
            RowNb<-row.nb/2 #number of columns, Even numb
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
        

            if (electrodeOrientation[svalue(orientation)]==1){
                elec.coords[,2]<-Xcord
                elec.coords[,3]<-Ycord
                elec.coords[,1]<-rep(cell.electrode.dist,length(Xcord))
            }
            if (electrodeOrientation[svalue(orientation)]==2){
                elec.coords[,3]<-Xcord
                elec.coords[,1]<-Ycord
                elec.coords[,2]<-rep(cell.electrode.dist,length(Xcord))
            }
            if (electrodeOrientation[svalue(orientation)]==3){
                elec.coords[,1]<-Xcord
                elec.coords[,2]<-Ycord
                elec.coords[,3]<-rep(cell.electrode.dist,length(Xcord))
            }
        }
        
        
        if (electrodeDistribute[svalue(eldistribute)]==4){
            elec.coords<-as.matrix(read.table(paste0(main.folder,'/simulation/ElcoordsDomi14.txt')))
        }
        if (cellTypes[svalue(chooseMethods)]==1){
            morpho.location<-paste(main.folder,'/simulation/morphology/ballstick.hoc',sep='')
        }
        if (cellTypes[svalue(chooseMethods)]==2){
            morpho.location<-paste(main.folder,'/simulation/morphology/villa.swc',sep='')
        }
        if (cellTypes[svalue(chooseMethods)]==3){
            morpho.location<-paste(main.folder,'/simulation/morphology/morpho1.swc',sep='')
        }
        
        if (cellTypes[svalue(chooseMethods)]==4){
            morpho.location<-paste(main.folder,'/simulation/morphology/neuron_agasbogas.swc',sep='')
        }
        if (cellTypes[svalue(chooseMethods)]==5){
            morpho.location<-paste(main.folder,'/simulation/morphology/Mainen_swcLike.swc',sep='')
        }
        ##You can choose an
        if (cellTypes[svalue(chooseMethods)]==6){
            morpho.location<-paste(main.folder,'/simulation/morphology/retina_ganglion.swc',sep='')
            ##morpho.location<-file.choose() }#.swc files work only!!
        }
        
        if (cellTypes[svalue(chooseMethods)]==7){
            morpho.location<-paste(main.folder,'/simulation/morphology/Badea2011Fig2Du.CNG.swc',sep='')
            
        }
        if (cellTypes[svalue(chooseMethods)]==8){
            morpho.location<-paste(main.folder,'/simulation/morphology/DomiCell.swc',sep='')
        }
        
###############################
        if (LFPysim[svalue( lfpysim)]==1){
            simName<-'LFP_calc.py'#_osc_sine.py'
        }
        if (LFPysim[svalue( lfpysim)]==2){
            simName<-'LFP_Y_symmetric.py' #'LFP_calc_active.py' 
        }
        if (LFPysim[svalue( lfpysim)]==3){
            simName<-'LFPymod_example6.py'
        }
        if (LFPysim[svalue( lfpysim)]==4){
            simName<-'LFP_calc_sine.py'#_osc_sine.py'
        }
        if (LFPysim[svalue( lfpysim)]==5){
            simName<-'LFP_calc_constInj.py'#_osc_sine.py'
        }
#####################################
        ##active channel description
        active.location<-paste(main.folder,'/simulation/morphology/active.hoc',sep='')
##############################
        
        simulation.location<-paste(main.folder,"/",svalue(simulation.location.name),sep='')
        cell.name<-svalue(cell.name.name)
        cat(cell.name)
        dir.create(paste(simulation.location,'/',cell.name,sep=''))
##############x
    
    
        ##plot(alak[,3],alak[,4])
        alak<-as.matrix(read.table(morpho.location, comment.char="#"))
        ##cat(alak)
        limx<-range(alak[,3],elec.coords[,1])
        limy<-range(alak[,4],elec.coords[,2])
        limz<-range(alak[,5],elec.coords[,3])
        ##plot(alak[,3],alak[,4])
        cat(paste('X range:', c(range(alak[,3]),'\n')),fill=TRUE)
        cat(paste('Y range:', c(range(alak[,4]),'\n')))
        cat(paste('Z range:', c(range(alak[,5]),'\n')))
        setupplot.name<-paste(simulation.location,'/',cell.name,'/setupplot.png',sep='')
        png(setupplot.name)
        sc<-scatterplot3d(alak[,3],alak[,4],alak[,5],pch=20,color='RED', xlim=limx, ylim=limy, zlim=limz,main='The cell and the electrode',xlab='x',ylab='y',zlab='z', scale.y=1)
        sc$points3d(elec.coords[,1],elec.coords[,2],elec.coords[,3],col='BLACK',pch=15)
        dev.off()
        sc<-scatterplot3d(alak[,3],alak[,4],alak[,5],pch=20,color='RED', xlim=limx, ylim=limy, zlim=limz,main='The cell and the electrode',xlab='x',ylab='y',zlab='z', scale.y=1)
        sc$points3d(elec.coords[,1],elec.coords[,2],elec.coords[,3],col='BLACK',pch=15)
        
        plotrange<-range(c(limx, limy,limz))
        plot3d(alak[,3],alak[,4],alak[,5],
               col='RED',xlim=plotrange,
               ylim=plotrange, zlim=plotrange,
               aspect=TRUE, xlab="x (um)",
               ylab="y (um)", "z (um)",
               main="Simulational Setup")
        plot3d(elec.coords[,1],elec.coords[,2],elec.coords[,3],col='BLACK',add=TRUE, size=10)
        aspect3d(1,1,1)
    
    
        elcoord.location.name<-paste(simulation.location,'/',cell.name,'/elcoord_x_y_z',sep='')
    
        write.table(c(elec.coords),elcoord.location.name,col.names=FALSE, row.names=FALSE)
        write.table(c(cell.electrode.dist, sigma),paste(simulation.location,'/',cell.name,'/elprop',sep=''),col.names=FALSE, row.names=FALSE)
        write.table(morpho.location,paste(simulation.location,'/',cell.name,'/morphology.txt',sep=''),col.names=FALSE, row.names=FALSE,quote=FALSE)
        write.table(active.location,paste(simulation.location,'/',cell.name,'/active.txt',sep=''),col.names=FALSE, row.names=FALSE,quote=FALSE)
        write.table(electrodeOrientation[svalue(orientation)],paste(simulation.location,'/',cell.name,'/ElecOrient.txt',sep=''),col.names=FALSE, row.names=FALSE,quote=FALSE)
###########
        write.table(cell.name,paste(simulation.location,'/cellname',sep=''),col.names=FALSE, row.names=FALSE,quote=FALSE)
        ##return(list(elec=elec.coords,alak))
        cat(simName)
        return(list(simName1=simName,simulation.location=simulation.location))
    }
    ##cat(elec.coords)
####################################
    win <- gwindow("Simulation of EP", width=400, height=700)
    
    gp <- ggroup(horizontal=FALSE, cont=win)
    ##Main folder -all the data should be in its subfoldersglabel("Main folder:", container=gp)
    main.folder.n<-gedit(DefaultDirectory,
                         container=gp, handler=LFPy_setup)
    
    
    
    
    glabel("Give a name to this simulation:", container=gp)
    cell.name.name<-gedit("cell_1", container=gp)
    
    glabel("Where to save the results", container=gp)
    simulation.location.name<- gedit("simulation", container=gp)
    
    
    
    
    ##Choosing the cell types
    tmp <- gframe("Cell types", container=gp, expand=TRUE)
    chooseMethods <- gradio(names(cellTypes), horizontal=TRUE,
                            cont=tmp,
                            handler=LFPy_setup)
  
  
  
    tmp <- gframe("Choosing the electrode positioning", container=gp, expand=TRUE)
    orientation <- gradio(names(electrodeOrientation), horizontal=TRUE,
                          cont=tmp,
                          handler=LFPy_setup)
    
    tmp <- gframe("Electrode Distribution", container=gp, expand=TRUE)
    eldistribute <- gradio(names(electrodeDistribute), horizontal=TRUE,
                           cont=tmp,
                           handler=LFPy_setup)
    
    tmp <- gframe("Choosing the LFPY simulation", container=gp, expand=TRUE)
    lfpysim <- gradio(names(LFPysim), horizontal=TRUE,
                      cont=tmp,
                      handler=LFPy_setup)
  
  
    
##########################
###############################
    ##Defining coordinates of the electrode
    tmp <- gframe("Parameters of the electrode", container=gp, expand=TRUE
                 ,horizontal=TRUE)
    
    ##cell to electrode distance
    widthofglabel<-6
    glabel("Cell to electrode distance", container=tmp)
    cellelectrodedist<-gedit("50", container=tmp, coerce.with=as.numeric,  handler=LFPy_setup,,width = widthofglabel)
    
    
    glabel("number of rows (x)", container=tmp)
    rownb<-gedit("4", container=tmp, coerce.with=as.numeric, handler=LFPy_setup,width = widthofglabel)
    glabel("number of columns (x)", container=tmp)
    colnb<-gedit("4", container=tmp, coerce.with=as.numeric, handler=LFPy_setup,width = widthofglabel)
    
    glabel("Hex Grid Const", container=tmp)
    triside<-gedit("19", container=tmp, coerce.with=as.numeric, handler=LFPy_setup,width = widthofglabel)
    
    
    glabel("Set Sewd Number", container=tmp)
    ssNb<-gedit("123456", container=tmp, coerce.with=as.numeric, handler=LFPy_setup,width = widthofglabel)
  
    
    ##Defining coordinates of the electrode
    tmp <- gframe("Parameters of the electrode1", container=gp, expand=TRUE
                 ,horizontal=TRUE)
    widthofglabel<-10
    glabel("a min (rows):", container=tmp)
    xmin<-gedit("-200", container=tmp, coerce.with=as.numeric, handler=LFPy_setup,width = widthofglabel)
  
    glabel("a max (rows):", container=tmp)
    xmax<-gedit("600", container=tmp, coerce.with=as.numeric, handler=LFPy_setup,width = widthofglabel)
  
    glabel("b min (columns)", container=tmp)
    ymin<-gedit("-200", container=tmp, coerce.with=as.numeric, handler=LFPy_setup, width = widthofglabel)
    glabel("b max (columns)", container=tmp)
    ymax<-gedit("200", container=tmp, coerce.with=as.numeric, handler=LFPy_setup,width = widthofglabel)
  
####################################x
    LFPy_run<-function(h,...) {
        ##simulation.location<-paste(svalue(main.folder.n),svalue(simulation.location.name),sep="")
        ##running the LFP simulation
        
        ##where to save the simulation results
    
    
        ##checking runtime
        ##Rprof(paste(outputfilename,'/profile_LFPsim.out',sep=''))
        simName2Be<-LFPy_setup()
        
        simulation.location<-simName2Be$simulation.location
        outputfilename<-paste(simulation.location,'/',svalue(cell.name.name),sep='')
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
        if(cellTypes[svalue(chooseMethods)]==1){
            source(paste0(svalue(main.folder.n),"/utils/sCSDFun.R"))
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
        What2Write2File<-paste(svalue(cell.name.name),svalue(chooseMethods),svalue( lfpysim),svalue(cellelectrodedist),svalue(orientation),svalue(xmin),svalue(xmax), svalue(ymin), svalue(ymax), svalue(colnb), svalue(rownb),sigma,  electrodeDistribute[svalue(eldistribute)],svalue(ssNb), "\n")
        if(file.exists(SummaryofSims)==FALSE)
            file.create(SummaryofSims)
        cat(What2Write2File,file=SummaryofSims,append=TRUE)
    
    
        setwd(wherearewenow)
        cat('Setup done')
        ##Rprof(NULL)
    }
  
 
  
  
###########################
  
    ##gfilebrowse (text = "Simulation location...", type = "selectdir", quote = TRUE, 
    ##        container =gp) 
  
  
    startSimulation<-gbutton("Start simulation", container=gp,
                           handler =LFPy_run )
  
    f_exit <-gbutton("Close simulation", container=gp,
                   handler = function( h,...) dispose( win ))
  
  
  
}
############

runLFPy() 
