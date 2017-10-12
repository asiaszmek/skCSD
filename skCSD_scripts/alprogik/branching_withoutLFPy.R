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

#read in the shape of the cell


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

#############################x

branch.length<-numeric(branch.nb) #length of branches
comp.length<-list(branch.nb) #compartments in the branches
comp.place<-list(branch.nb) #the distance measured from the beginning of the branch

for(j in 1: branch.nb){
lengthcomp<-numeric()
compplace<-0
lengthcomp<-sqrt(apply((alak[which(branches==j)[-1],3:5]-alak[which(branches==j)[-1]-1,3:5])^2,1,sum))
comp.length[[j]]<-lengthcomp
for (i in 1: (branch.comp.nb[j]-1)){
length<-lengthcomp[i]
compplace<-c(compplace,compplace[i]+length)
}
comp.place[[j]]<-compplace
branch.length[j]<-sum(comp.length[[j]])
}#j
remove(j,i,lengthcomp,compplace,length)
cell.length<-sum(branch.length) #length of the cell


vonal<-list(branch.nb)
vonalfun<-list(branch.nb)
for(j in 1:branch.nb){
#soma
t.param<-comp.place[[j]] # a távolságoktól függ... így a görbén egyenletesen lehet elosztani a pontokat.

#line interpolation az approx-szal
coords<-alak[which(branches==j),3:5] #segbranch.coords[[j]]
ts<-seq( from = 0, branch.length[j], length=100 )
vonal[[j]]<-apply( coords, 2, function(u) approx( t.param, u, xout = ts ,method="linear")$y ) 
vonalfun[[j]]<-apply( coords, 2, function(u) approxfun( t.param, u,method="linear") ) 
 
#else {vonal[[j]]<-alak[which(branches==j),c(3,4,5)]
#vonalfun[[j]]<-NA
#}
}
#curve(vonalfun[[j]]$V3(x),0,max(t.param))






