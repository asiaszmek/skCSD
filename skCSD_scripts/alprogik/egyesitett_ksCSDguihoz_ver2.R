library('scatterplot3d')
library('foreach')
library('doMC') 
library('fields')
library('MASS')

###############################
# TO DO
#kipróbálni, a többi morfológiára
#vizsgálni, hogy mennyit módosít a helyzeten, ha pl a hálózatosnál az alsó elektródát nem vesszük figyelembe.
#toconnect the fitted line, that the calculation-integration may go for a longer time


####Cross validation function
#ksCSD cross-validation
#cross.valid function gives back the value of the crossvalidation error and also the crossvalidated potential
cross.valid<-function(fi, transfermatrix){
  N<-dim(fi)[1]
  t<-dim(fi)[2]
  M<-dim(transfermatrix)[1]
  fi_star<-array(0,c(N,t))
  for(electrode in 1:N){
    
    skCSD_star<-array(0,c(M,t))
    skCSD_star<-transfermatrix[,-electrode]%*%fi[-electrode,]
    fi_star[electrode,]<-(ginv(transfermatrix)%*%skCSD_star)[electrode,]
  }#N
  cross.valid.error<-sqrt(sum((fi-fi_star)^2))/sqrt(sum((fi)^2))
  cross.valid.error.time<-sqrt(colSums((fi-fi_star)^2))/sqrt(colSums((fi)^2))
  return(list(cross.valid.error,fi_star,cross.valid.error.time))
}


##########################################
b.tilda.i<-function(i,R,source.cord.t,cord.t,j)
{
  out<-numeric()
  out<-exp(-((source.cord.t[i]-cord.t[j])^2)/(R^2))
  return(out)
}
#b.i(i,R,source.cord,elec.kord[j,])
b.i<-function(i,R,source.cord.t,r,length.fitted.curve){
  #x,y,z szerint integrálunk
  fun<-function(t2){
    if(0<t2 && t2<length.fitted.curve){
    xl<-fcx(t2)
    yl<-fcy(t2)
    zl<-fcz(t2)
    }
    if (t2>length.fitted.curve){ 
    xl<-fcx(t2-length.fitted.curve)
    yl<-fcy(t2-length.fitted.curve)
    zl<-fcz(t2-length.fitted.curve)
    }
    if (t2<0){ 
    xl<-fcx(t2+length.fitted.curve)
    yl<-fcy(t2+length.fitted.curve)
    zl<-fcz(t2+length.fitted.curve)
    }


    ff<- (exp(-((source.cord.t[i]-t2)^2)/(R^2)))/sqrt((r[1]-xl)^2+(r[2]-yl)^2+(r[3]-zl)^2)
    return(ff)
  }
  #integralt<-integrate(Vectorize(fun), lower=0,upper=length.fitted.curve,stop.on.error = FALSE)
  integralt<-integrate(Vectorize(fun), lower=-2*R,upper=length.fitted.curve+2*R,stop.on.error = FALSE)
  bi.value<-integralt$value
  #cat(integralt$subdivision,"\n")
  return(bi.value)
}
#########################################

cat('Functions to calculate ksCSD: called! \n')
##############################################
############################################
#############################################



#calculation of the curve fitten onto the morphology
library('scatterplot3d')


connection.function<-function(morpho,connections){
#morpho<-read.table('segcoordinates.txt')
#connections<-read.table('connections.txt')  
#skCSD<-read.table


seg.nb<-max(connections)


connection.matrix<-array(0,c(seg.nb,seg.nb))
for(i in 1: dim(connections)[1]) {connection.matrix[connections[i,2],connections[i,1]]<-1}
connection.matrix<-connection.matrix+t(connection.matrix)
connection.matrix[which(connection.matrix==2)]<-1

connection.matrix.history<-connection.matrix
rowsum.connenction.matrix<-rowSums(connection.matrix)
connection.line<-numeric()
connection.touched.once<-numeric()
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
  if (any(connection.line==point.now)==FALSE) connection.touched.once<-c(connection.touched.once,point.now)
  connection.line<-c(connection.line,point.now)
  
  previous.point<-current
  current<-point.now
}
connection.line[length(connection.line)]<-connection.line[1]


cat('Calculation of connection information: ready! \n')


#line interpolation az approx-szal
#coords.l<-morpho[connection.line,]
coords<-morpho[connection.line,]
length<-0
length.l<-0
comp.place<-0
comp.place.l<-0
for(i in 2:dim(coords)[1]){
  length<-sqrt(sum((coords[i,]-coords[i-1,])^2))
  comp.place<-c(comp.place,comp.place[i-1]+length)

}


t.param<-comp.place
length.fitted.curve<-max(comp.place)


ts<-seq( from = 0, length.fitted.curve, length=1000 )
vonal<-apply( coords, 2, function(u) approx( t.param, u, xout = ts ,method="linear")$y ) 
vonalfun<-apply( coords, 2, function(u) approxfun( t.param, u,method="linear") )
fcx<-vonalfun[[1]]
fcy<-vonalfun[[2]]
fcz<-vonalfun[[3]]


cat('ALMA') 
# vonalfun<-apply(coords,2,function(u) spline(t.param,u,xout=ts)$y)
cell.plot<-scatterplot3d(vonal, type="p", lwd=1,xlab='x (um)',ylab='y (um)', zlab='z (um)',
                         main='Fitted curve on cell morphology',angle=40,highlight.3d=TRUE,pch=20)
cell.plot<-scatterplot3d(vonal, type="l", lwd=3,xlab='x (um)',ylab='y (um)', zlab='z (um)',
                         main='Fitted curve on cell morphology',angle=40)
cell.plot$points3d(coords,type='h')

lines(cell.plot$xyz.convert(coords))
#return(list(connection.line=connection.line,t.param=t.param, length.fitted.curve=length.fitted.curve,vonalfun=vonalfun))
return(list(connection.line=connection.line,t.param=t.param, length.fitted.curve=length.fitted.curve,
            vonalfun=vonalfun, connection.touched.once=connection.touched.once, seg.nb=seg.nb))
}

################################################x
###################################################x
################################################################x

branching.result<-connection.function(morpho, connections)

connection.line<-branching.result$connection.line
t.param<-branching.result$t.param
seg.nb<-branching.result$seg.nb

length.fitted.curve<-branching.result$length.fitted.curve
vonalfun<-branching.result$vonalfun
fcx<-vonalfun[[1]]
fcy<-vonalfun[[2]]
fcz<-vonalfun[[3]]

###########################################################
######################################################
############################################

######################################x #memb.currents is optional
ksCSD_all<-function( basis.width.min, basis.width.max, basis.width.step, basis.number.min, basis.number.max, basis.number.step, LFP, elec.kord, memb.currents, seg.length, where2save){
  
  elec.kord<-matrix(elec.kord,ncol=3)
  el.nb<-dim(elec.kord)[1]
  
  
  
  
  ######
  
  
  
  cat('Entering the loop')
  
  ####
  ###
  
  #Loop for finding the best parameterss
  sigma<-0.5
  const<-1/(4*pi*sigma)
  base<-"gauss" #basis function type
  #to remember the best parameters
  M.best<-numeric()
  R.best<-numeric()
  params<-numeric()
  cv.error.best<-Inf #starting value of the best error
  mc.error.best<-Inf
  m.nb<-1
  
  
  R.all<-seq(basis.width.min,basis.width.max,basis.width.step)
  M.all<-seq(basis.number.min,basis.number.max,basis.number.step)
  errormatrix<-array(0,c(length(M.all),length(R.all)))
  proctimematrix<-array(0,c(length(M.all),length(R.all)))
  for(m in 1:length(M.all)){
    for(r in 1:length(R.all)){
      
      #profiling for each R and M
      proctime.m.r<-proc.time()
      M<-M.all[m]
      R<-R.all[r]
      
      
      base<-"gauss"
      #calculation of the location of basis functions
      ###################################################################
      ########################################
      ############x Basis functions
      #########################################
      
      
      # real number of the basis functions
      
      #t<-seq(0,max(comp.place[[j]]),length.out=source.branch.db[j])
      source.t<-seq(0,length.fitted.curve-3,length.out=M)
      
      #where.t<-seq(1,max(comp.place[[j]])-1,length.out=where.branch.db[j])
      #if the currents are known in specific points
      where.t<-t.param #or user defined
      #vonalfun
      #cat(vonalfun[[1]])
      #fcx<-vonalfun[[1]]
      #fcy<-vonalfun[[2]]
      #fcz<-vonalfun[[3]]
      
      
      
      
      where.cord<-matrix(c(fcx(where.t),fcy(where.t),fcz(where.t)),ncol=3)
      source.cord<-matrix(c(fcx(source.t),fcy(source.t),fcz(source.t)),ncol=3)
      ###############################################
      
      where.db<-dim(where.cord)[1]
      
      cat('Position of basis functions determined! \n')
      
      
      
      ######################x
      registerDoMC(cores=4)
      ################# 
      #Számoljuk ki a B illetve B.tilda mátrixot
      #egy sor egy adott i-hez tartozó fgv, oszlopokban azonos helyekhez tartozó
      #[i,j] : az i.függvény a j-dik helyen
      
      cat(paste('M:',M, 'R:',R, '\n'))
      source.nb<-M
      B.tilda<-array(0,c(source.nb,where.db))
      B<-array(0,c(source.nb,el.nb))
      for(i in 1:source.nb){
        Bj.result<-numeric(el.nb)
        #cat(paste(i,R,source.t,length.fitted.curve))
        Bj.result<-foreach(j=1:el.nb,.combine=c) %dopar% {
          b.i(i,R,source.t,elec.kord[j,],length.fitted.curve)
        }
        
        
        B[i,]<-Bj.result
        j<-0
        B.t.j.result<-numeric(where.db)
        B.t.j.result<-foreach(j=1:where.db,.combine=c) %dopar% {
          b.tilda.i(i,R, source.t,where.t,j) 
        }
        B.tilda[i,]<-B.t.j.result
        
        
      } #i
      
      
      K<-array(0,c(el.nb,el.nb))
      K.tilda<-array(0,c(source.nb,el.nb))
      
      K<-t(B)%*%B
      K.tilda<-t(B)%*%B.tilda
      
      #cat(dim(K))
      #cat(dim(K.tilda))
      
      #Let's use oseudoinverse
      C.calc<-1/const*t(K.tilda)%*%ginv(K)%*%LFP 
      
      cat('\n')
      image(C.calc)
      C.connection.touched.once<-C.calc
      #Lets put the current in the same order as they were in the simulation and add
      #up the current regarding the same spot
    
      skCSD.all<-array(0,c(seg.nb,dim(LFP)[2]))
      for(i in 1: seg.nb){
        sameplace<-numeric()
        sameplace<-c(which(branching.result$connection.line==i))
        for(j in 1:length(sameplace)){
          skCSD.all[i,]<-skCSD.all[i,]+C.calc[sameplace[j],]
        }      
        #skCSD.all[i,]<-skCSD.all[i,]/length(sameplace)
      }
      C.calc<-skCSD.all
    
      #C.calc<-C.calc[c(branching.result$connection.touched.once),]
      par(mfrow=c(1,2))
      #image(C.connection.touched.once,col=rainbow(200))
      image(C.calc,col=rainbow(200))
      ######################################
      #Is this result better than the previous ones
      cross.valid.version<-cross.valid(LFP,t(K.tilda)%*%ginv(K))
      cross.valid.error<-cross.valid.version[[1]]
      cross.valid.error.time<-cross.valid.version[[3]]
      fi_star<-cross.valid.version[[2]] #crossvalidated potential
      
      #in case we doesnt know the original membrane currents
      whicherror<-cross.valid.error
      #if we know
      
      errormatrix[m,r]<-cross.valid.error  
        
      params<-c(params,R,M,cross.valid.error)
      #ha a mostani jobb az elozo szimulacional
      if( cv.error.best >cross.valid.error){
        M.best<-M
        R.best<-R
        T.matrix<-t(K.tilda)%*%ginv(K)
        cv.error.best<-whicherror
	cv.error.best.time<-cross.valid.error.time
        C.best<-C.calc
        
        #writing it to file
      }
      ####################################
      ######## if we know the membrane currents
      #######################################x
      #f(exists("memb.currents") | exists("seg.length")){
      #szamol<-0
      #if(szamol==1){
      outname<-where2save
      
      
      if(memb.currents[1,1]!="NA" | seg.length!="NA"){
        #it is possible to compare the original and the  calculated currents
        #Estimation error of the currents
        
        #memb.currents<-as.matrix(read.table('membcurr'))
        #cat(paste('Segments lengths:',c(seg.length)))
        #cat(paste('dim of memb currents:', c(dim(memb.currents))))
        #cat(paste('dim of C:', c(dim(C.calc))))
        funaramvonal<-function(x) x/seg.length
        memb.currents.vonal<-apply(memb.currents,2,funaramvonal) 
        current.error.abs<-sqrt(sum((C.calc-memb.currents.vonal)^2))/sqrt(sum((memb.currents.vonal)^2))
        
        write.table(memb.currents.vonal,paste(outname,'/membcurr_line',sep=''),col.names=FALSE, row.names=FALSE)
        
    
        
        
        #Smoothed membrane currents
        
        lambda<-R
        #how far are from eachother the neighbouring segments... smoothing based on this
        length<-numeric()
        comp.place.curr<-0
        for(i in 2:dim(morpho)[1]){
          length<-sqrt(sum((morpho[i,]-morpho[i-1,])^2))
          comp.place.curr<-c(comp.place.curr,comp.place.curr[i-1]+length)
          
        }
        memb.currents.smoothed<-array(0,dim(memb.currents.vonal))
        for (t in 1: dim(memb.currents.vonal)[2]){
        memb.currents.smoothed[,t]<-ksmooth(comp.place.curr,memb.currents.vonal[,t],"normal", bandwidth=lambda,x.points=comp.place.curr)$y
        }
        #error compared to the smoothed membrane currents
        current.error.smoothed.abs<-sqrt(sum((C.calc-memb.currents.smoothed)^2))/sqrt(sum((memb.currents.smoothed)^2))
	#error compared to the smoothed membrane currents in time
        current.error.smoothed.time<-sqrt(colSums((C.calc-memb.currents.smoothed)^2))/sqrt(colSums((memb.currents.smoothed)^2))
        params<-c(params,current.error.abs,current.error.smoothed.abs)
        
        if( mc.error.best > current.error.smoothed.abs){
          M.best.mc<-M
          R.best.mc<-R
          T.matrix.mc<-t(K.tilda)%*%ginv(K)
          mc.error.best<- current.error.smoothed.abs
          C.best.mc<-C.calc
	  current.error.smoothed.time.best<-current.error.smoothed.time
          
          #writing it to file
        }
        
        
        
      } #current known
      
      #printing out the profiling result and writing it to file
      proc.time()-proctime.m.r
      proctimematrix[m,r]<-(proc.time()-proctime.m.r)[1]
     
    }} #for m and r
  
  proc.name<-paste(outname,"/proctime",sep="")
  colnames(proctimematrix)<-R.all
  rownames(proctimematrix)<-M.all
  write.table(proctimematrix,proc.name)
  #writing things to file
  #different error values
  #the best parameters
  
  write.table(c(M.best, R.best,cv.error.best),paste(outname,'/bestparams_M_R_error',sep=''),col.names=FALSE, row.names=FALSE)
  write.table(cv.error.best.time,paste(outname,'/cvErrorBestTime',sep=''),col.names=FALSE, row.names=FALSE)
  #all the parameters
  params<-matrix(params,ncol=(length(M.all)*length(R.all)))
  write.table(params,paste(outname,'/params',sep=''),col.names=FALSE, row.names=FALSE)
  #the transfer matrix in the best situation
  write.table(T.matrix,paste(outname,'/transfermatrix',sep=''),col.names=FALSE, row.names=FALSE)
  write.table(ginv(K),paste(outname,'/ginv_K',sep=''),col.names=FALSE, row.names=FALSE)
  write.table(t(K.tilda),paste(outname,'/t_K_tilda',sep=''),col.names=FALSE, row.names=FALSE)
  #the current in the best case
  write.table(C.best,paste(outname,'/skCSD',sep=''),col.names=FALSE, row.names=FALSE)
  #the error in case of smoothed current comparison- best case
  write.table(current.error.smoothed.time.best,paste(outname,'/BestSetupErrorinTime',sep=''),col.names=FALSE, row.names=FALSE)
  #smoothed membrane current
  write.table(memb.currents.smoothed,paste(outname,'/membcurr_smoothed',sep=''),col.names=FALSE, row.names=FALSE)
  #write to file, if have the comparison with the membrane currents
  write.table(c(M.best.mc, R.best.mc, mc.error.best),paste(outname,'/bestparams_MC__M_R_error',sep=''),col.names=FALSE, row.names=FALSE)
  #the transfer matrix in the best situation
  write.table(T.matrix.mc,paste(outname,'/transfermatrix_MC',sep=''),col.names=FALSE, row.names=FALSE)
  #the current in the best case
  write.table(C.best.mc,paste(outname,'/skCSD_MC',sep=''),col.names=FALSE, row.names=FALSE)
  
  #where did we calculate the currents?
  write.table(where.cord, paste(outname,'/where_cord',sep=''),col.names=FALSE, row.names=FALSE)
  #how do the points follow eachother
  write.table(branching.result$connection.touched.once, paste(outname,'/connection_touched_once',sep=''),col.names=FALSE, row.names=FALSE)
  write.table(C.connection.touched.once, paste(outname,'/C_T1',sep=''),col.names=FALSE, row.names=FALSE)
  
  return(C.best)
  
} #big ksCSD function






