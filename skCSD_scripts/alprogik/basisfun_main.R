####################### kCSD
#source('alprogik/basisfun_branching.R')

#i function,in r position, i source's coordinate: source.cord
#r_t ahol ki akarom számolni
#lets use parametrized version

#The calculations should just be done in the same branches...
#b.tilda.i(i,R,source.cord,where.cord[j,])
b.tilda.i<-function(i,R,source.cord.t,cord.t,j)
{
out<-numeric()
if(branch.source[i] == where.which[j]) {
if (base=='gauss') {out<-exp(-((source.cord.t[i]-cord.t[j])^2)/(R^2))}
#if (base=='step'){ 
#thresh1<-((source.cord[i,1]-xl)^2+(source.cord[i,2]-yl)^2+(source.cord[i,3]-zl)^2)
#if ( thresh1<= R^2)  out<-1 #itt eredetileg (R/2 volt)
#else out<-0
#} #step

if (base=='cos'){ 
distan1<-sqrt((source.cord.t[i]-cord.t[j])^2)
value<-pi/R*distan1
if (distan1<= R/2) out<-cos(value)
if(distan1 > R/2) out<-0}
}
if(branch.source[i] != where.which[j]) out<-0
return(out)
}


#b.i(i,R,source.cord,elec.kord[j,])
b.i<-function(i,R,source.cord.t,r){
branchwhich<-branch.source[i]
#x,y,z szerint integrálunk
fun<-function(t2){
#t1 0-tól 1ig megy
#a t-nek
#t<-0+(max(comp.place[[branchwhich]-0)]*t1
#t2<-max(comp.place[[branchwhich]])*t1
fcx<-vonalfun[[branchwhich]][[1]]
fcy<-vonalfun[[branchwhich]][[2]]
fcz<-vonalfun[[branchwhich]][[3]]

xl<-fcx(t2)
yl<-fcy(t2)
zl<-fcz(t2)

if (base=='gauss'){ ff<- (exp(-((source.cord.t[i]-t2)^2)/(R^2)))/sqrt((r[1]-xl)^2+(r[2]-yl)^2+(r[3]-zl)^2)}

#if (base=='step'){
#thresh1<-((source.cord[i,1]-xl)^2+(source.cord[i,2]-yl)^2+(source.cord[i,3]-zl)^2)
#if ( thresh1<= R^2) {ff<- 1/sqrt((r[1]-xl)^2+(r[2]-yl)^2+(r[3]-zl)^2)*sqrt((b[1]-a[1])^2+(b[2]-a[2])^2+(b[3]-a[3])^2) }
#else ff<-0
#}

if (base=='cos'){ 
distan1<-sqrt((source.cord.t[i]-t2)^2)
value<-pi/R*distan1
if (distan1<= R/2) ff<-cos(value)/sqrt((r[1]-xl)^2+(r[2]-yl)^2+(r[3]-zl)^2)
if(distan1 > R/2) ff<-0}

return(ff)
}


#integralt<-integrate(Vectorize(fun), -Inf,Inf,rel.tol=0.001)
#it is possible,that the
#integralt<-integrate(Vectorize(fun),source.cord.t[i]-3*R,source.cord.t[i]+3*R)
integralt<-integrate(Vectorize(fun), 0,max(comp.place[[branchwhich]]),subdivisions=200)
#integralt<-integrate(Vectorize(fun), (source.cord[i,3]-a[3]-3*R)/(b[3]-a[3]),(source.cord[i,3]-a[3]+3*R)/(b[3]-a[3]))
#integralt<-cuhre(1,1,fun,lower=as.vector((source.cord[i,3]-10*R)/(b[3]-a[3])),upper=as.vector( (source.cord[i,3]+10*R)/(b[3]-a[3])),flags=list(verbose=0))
#integralt<-cuhre(1,1,Vectorize(fun),lower=as.vector((-0.5)),upper=as.vector(1 ),flags=list(verbose=0))

#bi.value<-const*integralt$value
bi.value<-integralt$value
#cat(integralt$subdivision,"\n")
return(bi.value)
}


