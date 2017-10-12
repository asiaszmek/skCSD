
DirDataName<-c("simulation/cell_DomiShifted","simulation/gang_9x9_50","simulation/gang_17x17_50","simulation/gang_9x9_100")#,"simulation/gang_5x5_400")#simulation/MCS",, "simulation/Ganglion_d17.5_128Regular", )
 # DirData<- DirDataName[1]
#source("test10/kernel_withoutGUI.R")


for(i in 1:4){
  DirData<- DirDataName[i]
source("test10/kernel_withoutGUI.R")
cat(i)
}
