library('scatterplot3d')
#morpho<-read.table('/media/BA0ED4600ED416EB/agy/kCSD/progik/bs_futtat/branching_simple2014/network_Hela/morpho_pyr_dori.swc') #properties of segments
#connections<-read.table('/media/BA0ED4600ED416EB/agy/kCSD/progik/bs_futtat/branching_simple2014/network_Hela/suppyrRS_template_dori.hoc') #connection information

connection.function<-function(morpho,connections){
seg.nb<-max(connections[,2:3])
cat(dim(connections))

connection.matrix<-array(0,c(seg.nb,seg.nb))
for(i in 1: dim(connections)[1]) {connection.matrix[connections[i,3],connections[i,2]]<-1}
connection.matrix<-connection.matrix+t(connection.matrix)
connection.matrix.history<-connection.matrix
rowsum.connenction.matrix<-rowSums(connection.matrix)
connection.line<-numeric()
current<-1
previous.point<-0
#while (sum(connection.matrix)>1){
while (previous.point!=current){

rowsum.connection.matrix<-rowSums(connection.matrix)
rowsum.connection.matrix.history<-rowSums(connection.matrix.history)
if (rowsum.connection.matrix[current] >1) {
#if there is a point we haven't touched that, go to that direction
	if (rowsum.connection.matrix.history[current] >0){
	which.connected<-which(connection.matrix.history[current,]==1)
	}
	if (rowsum.connection.matrix.history[current] ==0){
	which.connected<-which(connection.matrix[current,]==1)
	}
#point.now<-which.connected[which(which.connected!=previous.point)][1]

point.now<-which.connected[length(which.connected)]
}
#it can only vonnect to the previous one
if (rowsum.connection.matrix[current] ==1) {point.now<-which(connection.matrix[current,]==1)[1]
}

connection.matrix[current,point.now]<-0
connection.matrix.history[current,point.now]<-0
connection.matrix.history[point.now,current]<-0
connection.line<-c(connection.line,point.now)
previous.point<-current
current<-point.now
}
connection.line[length(connection.line)]<-connection.line[1]


cat('Calculation of connection information: ready! \n')


#line interpolation az approx-szal
coords<-morpho[connection.line,3:5]
length<-0
comp.place<-0
for(i in 2:dim(coords)[1]){
length<-sqrt(sum((coords[i,]-coords[i-1,])^2))
comp.place<-c(comp.place,comp.place[i-1]+length)
}
t.param<-comp.place
length.fitted.curve<-max(comp.place)

ts<-seq( from = 0, length.fitted.curve, length=1000 )
vonal<-apply( coords, 2, function(u) approx( t.param, u, xout = ts ,method="linear")$y ) 
vonalfun<-apply( coords, 2, function(u) approxfun( t.param, u,method="linear") )

cat('ALMA') 
# vonalfun<-apply(coords,2,function(u) spline(t.param,u,xout=ts)$y)
cell.plot<-scatterplot3d(vonal, type="p", lwd=1,xlab='x (um)',ylab='y (um)', zlab='z (um)',
main='Fitted curve on cell morphology',angle=40,highlight.3d=TRUE,pch=20)
cell.plot<-scatterplot3d(vonal, type="l", lwd=3,xlab='x (um)',ylab='y (um)', zlab='z (um)',
main='Fitted curve on cell morphology',angle=40)
cell.plot$points3d(coords,type='h')

lines(cell.plot$xyz.convert(coords))
return(list(connection.line=connection.line,t.param=t.param, length.fitted.curve=length.fitted.curve,vonalfun=vonalfun))
}
