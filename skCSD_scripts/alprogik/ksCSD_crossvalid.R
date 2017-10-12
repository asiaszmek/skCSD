#ksCSD cross-validation
#cross.valid function gives back the value of the crossvalidation error and also the crossvalidated potential
cross.valid<-function(fi, transfermatrix){
N<-dim(fi)[1]
t<-dim(fi)[2]
M<-dim(transfermatrix)[1]
fi_star<-array(0,c(N,t))
for(electrode in 1:N){
cat(electrode)
skCSD_star<-array(0,c(M,t))
skCSD_star<-transfermatrix[,-electrode]%*%fi[-electrode,]
fi_star[electrode,]<-(ginv(transfermatrix)%*%skCSD_star)[electrode,]
}#N
cross.valid.error<-sum(abs(fi-fi_star))/sum(abs(fi))
return(list(cross.valid.error,fi_star))
}#function
#zelim<-range(LFP,fi_star)
#image(t(as.matrix(LFP)),zlim=zelim)
#X11()
#image(t(fi_star),zlim=zelim)

