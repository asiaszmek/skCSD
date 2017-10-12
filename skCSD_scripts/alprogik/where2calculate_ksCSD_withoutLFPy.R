########################################
############x Basis functions
#########################################
source.branch.db<-round(branch.length/cell.length*M)
branch.source<-rep(1:branch.nb,source.branch.db)
where.branch.db<-branch.comp.nb
# real number of the basis functions
M<-sum(source.branch.db) #?
source.cord<-array(0,c(M,3))

b.cordx<-numeric()
b.cordy<-numeric()
b.cordz<-numeric()
where.cordx<-numeric()
where.cordy<-numeric()
where.cordz<-numeric()
source.cord.t<-numeric()
where.cord.t<-numeric()
where.which<-numeric()#which "where" is on which branch
for (j in 1:branch.nb){
#t<-seq(0,max(comp.place[[j]]),length.out=source.branch.db[j])
t<-seq(3,max(comp.place[[j]])-3,length.out=source.branch.db[j])

fcx<-vonalfun[[j]][[1]]
fcy<-vonalfun[[j]][[2]]
fcz<-vonalfun[[j]][[3]]

where.t<-t
#where.t<-seq(1,max(comp.place[[j]])-1,length.out=where.branch.db[j])
#where.t<-(comp.place[[j]][1:branch.comp.nb[j]]+seg.length[which(branches==j)]/2) #in case there is no simulation, this can be different


b.cordx<-c(b.cordx,fcx(t))
b.cordy<-c(b.cordy,fcy(t))
b.cordz<-c(b.cordz,fcz(t))
where.cordx<-c(where.cordx,fcx(where.t))
where.cordy<-c(where.cordy,fcy(where.t))
where.cordz<-c(where.cordz,fcz(where.t))
source.cord.t<-c(source.cord.t,t)
where.cord.t<-c(where.cord.t,where.t)
where.which<-c(where.which,rep(j,where.branch.db[j]))
}
where.db<-length(where.cordx)
where.cord<-array(0,c(where.db,3))#array(0,c(sum(where.branch.db),3))
where.cord[,1]<-where.cordx
where.cord[,2]<-where.cordy
where.cord[,3]<-where.cordz
source.cord[,1]<-b.cordx
source.cord[,2]<-b.cordy
source.cord[,3]<-b.cordz
remove(t,j,b.cordx,b.cordy,b.cordz,where.cordx,where.cordy,where.cordz)
###############################################

