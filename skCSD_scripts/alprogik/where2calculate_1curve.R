###################################################################
########################################
############x Basis functions
#########################################


# real number of the basis functions

#t<-seq(0,max(comp.place[[j]]),length.out=source.branch.db[j])
source.t<-seq(0,length.fitted.curve-3,length.out=M)

fcx<-vonalfun[[1]]
fcy<-vonalfun[[2]]
fcz<-vonalfun[[3]]

#where.t<-seq(1,max(comp.place[[j]])-1,length.out=where.branch.db[j])
#if the currents are known in specific points
where.t<-t.param #or user defined

where.cord<-matrix(c(fcx(where.t),fcy(where.t),fcz(where.t)),ncol=3)
source.cord<-matrix(c(fcx(source.t),fcy(source.t),fcz(source.t)),ncol=3)
###############################################

where.db<-dim(where.cord)[1]

cat('Position of basis functions determined! \n')
