#I took a random cell morphology from neuromorpho
#http://neuromorpho.org/neuroMorpho/neuron_info.jsp?neuron_name=03a_pyramidal9aFI
#This program:
#gets the coordinates
#finds the branching points
#calculates a function for the branching points

#melyik szegmenshez kapcsolódik,soma-apical-basal,x,y,z, d, melyik másik szegmenshez kapcsolódik

#alak:shape of the cell.
#alak<-as.matrix(read.table('morphology/pyramid.txt'))
#name of moprhology file

#setwd(workingdirectory)morphoname<-scan('morphology.txt',what='charachter')
#read in the shape of the cell
#setwd('..')
#alak<-as.matrix(read.table(morphoname, comment.char="#"))

#SWC: There are four default segment tags in SWC: 1 soma, 2 axon, 3 dendrite (basal when apical present), and 4 apical dendrite.
#The three dimensional structure of a neuron can be represented in a SWC format (Cannon et al., 1998). SWC is a simple Standardized format. Each line has 7 fields encoding data for a single neuronal compartment:
#an integer number as compartment identifier
#type of neuronal compartment 
#x coordinate of the compartment
#y coordinate of the compartment
#z coordinate of the compartment
#radius of the compartment
#parent compartment

comp.nb<-dim(alak)[1] #number of compartments
soma<-length(alak[which(alak[,2]==1)]) #how many points represent the soma
branching.points<-numeric() #finding the branching points
branches<-numeric() #which segment belongs to which branch
branch.nb<-numeric() #number of branches

cell.length<-numeric()
##########################
j<-0 #just some counters

branches<-c(rep(1,soma))
j<-j+1
for (i in (soma+1): (comp.nb)){
if(alak[i,1]!=(alak[i,7]+1) | alak[i,2]!=(alak[i-1,2])){ 
#branching.points<-c(branching.points,alak[which(alak[,1]==alak[i,7]),1]) #ez itt hibásan volt, de nem is értem miért
branching.points<-c(branching.points,i)
j<-j+1
}
branches<-c(branches,j)
}

remove(j,i)
####################################
branch.nb<-max(branches) #number of branches
branch.comp.nb<-numeric(branch.nb) #how many compments are there in a branch

for(k in 1: branch.nb){
branch.comp.nb[k]<-length(which(branches==k))
}
####################################
#notplot<-1
#if(notplot==0){
#colours<-color.scale(c(1:branch.nb),col=rainbow(branch.nb))
#limx<-range(alak[,3],elec.kord[,1])
#limy<-range(alak[,4],elec.kord[,2])
#limz<-range(alak[,5],elec.kord[,3])
#sc<-scatterplot3d(alak[,3],alak[,5],alak[,4],color=colours,pch=20, xlim=limx, ylim=limz, zlim=limy,main='The cell and the electrode')
#sc$points3d(alak[branching.points,3],alak[branching.points,5],alak[branching.points,4],col='BLACK',pch=20)
#sc$points3d(elec.kord[,1],elec.kord[,3],elec.kord[,2],col='BLACK',pch=15)
#legend('topright',as.character(c(1:branch.nb)),col=rainbow(branch.nb), lty=1,bg='WHITE')

#colours<-color.scale(branches,col=rainbow(branch.nb+1))
#connection.tree<-list()
#connection.tree[[1]]<-0
#for(i in 1:(branch.nb-1)){
#connection.tree[[i+1]]<-c(connection.tree[[branches[alak[branching.points[i],7]]]],branches[alak[branching.points[i],7]])
#}

#t.snap<-which(skCSD==max(skCSD),arr.ind=TRUE)[2] #which timesnap to plot
#maxdist<-numeric()
#comp.dist<-numeric()
#comp.dist<-list()
#comp.dist[[1]]<-0

#for(i in 1:branch.nb){
#connection.tree[[i]]<-connection.tree[[i]]+1
#comp.dist[[i+1]]<-comp.place[[i+1]][-1]+sum(branch.length[connection.tree[[i]]])
#maxdist<-max(maxdist,max(comp.dist[[i]]))
#}

#par(mfrow=c(2,1))

#memb.currents<-read.table('membcurr')
#seg.length<-as.matrix(read.table('seglength'))
#funaramvonal<-function(x) x/seg.length

#membcurr<-apply(memb.currents,2,funaramvonal) 
#currlimits<-range(skCSD[,t.snap],membcurr[,t.snap])


#plot(c(comp.dist[[1]]),skCSD[which(segbranches==1),t.snap],xlim=c(0,maxdist), ylim=range(skCSD[,t.snap]),t='l',col=colours[1])
#for(i in 2:branch.nb){
#lines(c(comp.dist[[i]]),skCSD[which(segbranches==i),t.snap],col=colours[i])
#}



#plot(c(comp.dist[[1]]),membcurr[which(segbranches==1),t.snap],xlim=c(0,maxdist), ylim=range(skCSD[,t.snap]),t='l',col=colours[1]) #range(membcurr)
#for(i in 2:branch.nb){
#lines(c(comp.dist[[i]]),membcurr[which(segbranches==i),t.snap],col=colours[i])
#}


#(seg.start+seg.end)/2
#}

