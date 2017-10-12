########################################################
######################################################
###################################################
#Calculation of basis functions

b.tilda.i<-function(i,R,source.cord.t,cord.t,j)
{
out<-numeric()
if (base=='gauss') {out<-exp(-((source.cord.t[i]-cord.t[j])^2)/(R^2))}
return(out)
}


#b.i(i,R,source.cord,elec.kord[j,])
b.i<-function(i,R,source.cord.t,r,length.fitted.curve){

#x,y,z szerint integrálunk
fun<-function(t2){
#t1 0-tól 1ig megy
#a t-nek
#t<-0+(max(comp.place[[branchwhich]-0)]*t1
#t2<-max(comp.place[[branchwhich]])*t1

xl<-fcx(t2)
yl<-fcy(t2)
zl<-fcz(t2)

if (base=='gauss'){ ff<- (exp(-((source.cord.t[i]-t2)^2)/(R^2)))/sqrt((r[1]-xl)^2+(r[2]-yl)^2+(r[3]-zl)^2)}

return(ff)
}


#integralt<-integrate(Vectorize(fun), -Inf,Inf,rel.tol=0.001)
#it is possible,that the
#integralt<-integrate(Vectorize(fun),source.cord.t[i]-3*R,source.cord.t[i]+3*R)
integralt<-integrate(Vectorize(fun), 0,length.fitted.curve)
#integralt<-integrate(Vectorize(fun), (source.cord[i,3]-a[3]-3*R)/(b[3]-a[3]),(source.cord[i,3]-a[3]+3*R)/(b[3]-a[3]))
#integralt<-cuhre(1,1,fun,lower=as.vector((source.cord[i,3]-10*R)/(b[3]-a[3])),upper=as.vector( (source.cord[i,3]+10*R)/(b[3]-a[3])),flags=list(verbose=0))
#integralt<-cuhre(1,1,Vectorize(fun),lower=as.vector((-0.5)),upper=as.vector(1 ),flags=list(verbose=0))
#bi.value<-const*integralt$value
bi.value<-integralt$value
#cat(integralt$subdivision,"\n")
return(bi.value)
}


cat('Functions to calculate ksCSD: called! \n')

