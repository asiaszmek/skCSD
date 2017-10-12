#LFPy branches
setwd(data.location)
seg.start<-matrix(as.matrix(read.table('coordsstart_x_y_z')),ncol=3)
seg.end<-matrix(as.matrix(read.table('coordsend_x_y_z')),ncol=3)
seg.length<-as.matrix(read.table('seglength'))
seg.db<-length(seg.length)

segbranches<-1
br<-1
for (i in 1:(seg.db-1)){
if (seg.start[i+1,1]==seg.end[i,1] & seg.start[i+1,2]==seg.end[i,2] & seg.start[i+1,3]==seg.end[i,3]){
 segbranches<-c(segbranches,br)
}

if (seg.start[i+1,1]!=seg.end[i,1] | seg.start[i+1,2]!=seg.end[i,2] | seg.start[i+1,3]!=seg.end[i,3]) {
br<-br+1
segbranches<-c(segbranches,br)
}
}

#segbranch.num<-length(seg.diam)
segbranch.nb<-max(segbranches)

segbranch.comp.nb<-numeric(segbranch.nb)
for(k in 1: segbranch.nb){
segbranch.comp.nb[k]<-length(which(segbranches==k))
}



segbranch.coords<-list()
for (j in 1:segbranch.nb){

cordmatr<-array(0,c(segbranch.comp.nb[[j]]+1,3))
cordmatr[1:segbranch.comp.nb[j],]<-seg.start[which(segbranches==j),]
cordmatr[segbranch.comp.nb[j]+1,]<-seg.end[which(segbranches==j)[segbranch.comp.nb[j]],]
segbranch.coords[[j]]<-cordmatr
}


branch.length<-numeric(segbranch.nb) #length of branches
comp.length<-list(segbranch.nb) #compartments in the branches
comp.place<-list(segbranch.nb) #the distance measured from the beginning of the branch

for(j in 1: segbranch.nb){
lengthcomp<-numeric()
compplace<-0

for (i in 1: (segbranch.comp.nb[j])){
wh<-sum(segbranch.comp.nb[0:(j-1)])+i 
length<-sqrt(sum((seg.end[wh,]-seg.start[wh,])^2))

lengthcomp<-c(lengthcomp,length)
compplace<-c(compplace,compplace[i]+length)
}
comp.place[[j]]<-compplace
comp.length[[j]]<-lengthcomp
branch.length[j]<-sum(comp.length[[j]])
}#j
remove(j,i,lengthcomp,compplace,length)
cell.length<-sum(branch.length) #length of the cell



vonal<-list(segbranch.nb)
vonalfun<-list(segbranch.nb)
for(j in 1:segbranch.nb){
#soma


t.param<-comp.place[[j]] # a távolságoktól függ... így a görbén egyenletesen lehet elosztani a pontokat.

#line interpolation az approx-szal
coords<-segbranch.coords[[j]]
ts<-seq( from = 0, branch.length[j], length=100 )
vonal[[j]]<-apply( coords, 2, function(u) approx( t.param, u, xout = ts ,method="linear")$y ) 
vonalfun[[j]]<-apply( coords, 2, function(u) approxfun( t.param, u,method="linear") ) 
 
#else {vonal[[j]]<-alak[which(branches==j),c(3,4,5)]
#vonalfun[[j]]<-NA
#}
}


