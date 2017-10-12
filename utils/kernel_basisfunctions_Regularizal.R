library('scatterplot3d')
library('foreach')
library('doMC') 
#library('fields')
library('MASS')

###############################

#if(file.exists(outname)) file.remove(list.files("outname"))

base <-"gauss" #"poli"
# TO DO
#kipróbálni, a többi morfológiára
#vizsgálni, hogy mennyit módosít a helyzeten, ha pl a hálózatosnál az alsó elektródát nem vesszük figyelembe.
#toconnect the fitted lines, that the calculation-integration may go for a longer time


####Cross validation function
#ksCSD cross-validation
#cross.valid function gives back the value of the crossvalidation error and also the crossvalidated potential
cross.valid<-function(fi,
                      transfermatrix){
    N<-dim(fi)[1]
    t<-dim(fi)[2]
    M<-dim(transfermatrix)[1]
    cross.valid.error<-numeric()
    cross.valid.error.time<-numeric()
    fi_star<-array(0,c(N,t))
    for(electrode in 1:N){
    
        ##skCSD_star<-array(0,c(M,t))
        skCSD_star<-transfermatrix[,-electrode]%*%fi[-electrode,]
        fi_star<-(ginv(transfermatrix)%*%skCSD_star)
	cross.valid.error<-c(cross.valid.error,sqrt(sum((fi-fi_star)^2))/sqrt(sum((fi)^2)))
  	cross.valid.error.time<-c(cross.valid.error.time,sqrt(colMeans((fi-fi_star)^2))/sqrt(colSums((fi)^2)))

    }#N
    cross.valid.error.all<-mean(cross.valid.error)
    cross.valid.error.time.all<-mean(cross.valid.error.time)
    return(list(cross.valid.error.all,fi_star,cross.valid.error.time.all))
}


##########################################
b.tilda.i<-function(i,
                    R,
                    source.cord.t,
                    cord.t,j){
    out<-numeric()
    out<-exp(-((source.cord.t[i]-cord.t[j])^2)/(R^2))
    return(out)
}
#b.i(i,R,source.cord,elec.kord[j,])
b.i<-function(i,
              R,
              source.cord.t,
              r,
              length.fitted.curve){
    ##x,y,z szerint integrálunk
    fun<-function(t2){
        ##periodic boundary condition
        ##if(0<t2 && t2<length.fitted.curve){
        ##xl<-fcx(t2)
        ##yl<-fcy(t2)
        ##zl<-fcz(t2)
        ##}
        ##if (t2>length.fitted.curve){ 
        ##xl<-fcx(t2-length.fitted.curve)
        ##yl<-fcy(t2-length.fitted.curve)
        ##zl<-fcz(t2-length.fitted.curve)
        ##}
        ##if (t2<0){ 
        ##xl<-fcx(t2+length.fitted.curve)
        ##yl<-fcy(t2+length.fitted.curve)
        ##zl<-fcz(t2+length.fitted.curve)
        ##}
  
        if (t2>length.fitted.curve){ 
            t2<-t2-length.fitted.curve
        }
        else if (t2<0){ 
            t2<-t2+length.fitted.curve
        }
        xl<-fcx.l(t2)
        yl<-fcy.l(t2)
        zl<-fcz.l(t2)
        
        if(base =="poli"){
            if(abs(source.cord.t[i]-t2)<2*R){
                ff<-(1-(abs(source.cord.t[i]-t2)/(2*R))^4)/sqrt((r[1]-xl)^2+(r[2]-yl)^2+(r[3]-zl)^2)}
            else {
                ff<-0}
        }
        if(base=="gauss" )
            ff<- (exp(-((source.cord.t[i]-t2)^2)/(R^2)))/sqrt((r[1]-xl)^2+(r[2]-yl)^2+(r[3]-zl)^2)
            #cat(paste(r[1],xl,r[2],yl,r[3],zl,"\n",sep=" "))
        return(ff)
    }
    ##integralt<-integrate(Vectorize(fun), lower=0,upper=length.fitted.curve,stop.on.error = FALSE)
    ##the integration
    integralt<-integrate(Vectorize(fun), lower=-2*R,upper=length.fitted.curve+2*R,stop.on.error = TRUE)
    ##integralt<-integrate(Vectorize(fun), lower=source.cord.t[i]-2*R,upper=source.cord.t[i]+2*R,stop.on.error = FALSE) #!!!!!!!!
    
    bi.value<-integralt$value #/count.touched[i]
    ##cat(integralt$subdivision,"\n")
    return(bi.value)
}
#########################################

cat('Functions to calculate ksCSD: called! \n')
##############################################
############################################
#############################################



#calculation of the curve fitten onto the morphology
library('scatterplot3d')

#connection.function(morpho,connections)
connection.function<-function(morpho,
                              connections){
    ##morpho<-as.matrix(read.table('segcoordinates.txt'))
    ##morpho<-matrix(as.numeric(as.matrix(read.table('segcoordinates.txt'))),ncol=3)
    ##  membcurr<-matrix(as.numeric(as.matrix(read.table('simulation/cell_1/membcurr'))))
    ##connections<-as.matrix(read.table('connections.txt')  )
    ##skCSD<-read.table


    seg.nb<-max(connections)


    connection.matrix<-array(0,c(seg.nb,seg.nb))
    for(i in 1: dim(connections)[1]) {
        connection.matrix[connections[i,2],connections[i,1]]<-1 }
    connection.matrix<-connection.matrix+t(connection.matrix)
    
    rowsum.connenction.matrix<-rowSums(connection.matrix)
    
    
    connection.line<-numeric()
    connection.touched.once<-numeric()
    current<-1 
    Touched<-rep(0,seg.nb)

    while (rowSums(connection.matrix)[current]==0)
        current<-current + 1 #in case the first segment is not connected
    connection.line<-current
    connection.touched.once<-current
    previous.point<-0
    which.connected<-1
    ##while (sum(connection.matrix)>1){
    while (length(which.connected)>0){
        ##while(any(Touched==0)){
        ##if(rowsum.connenction.matrix.t[current]==0) # 
        which.connected<-which(connection.matrix[current,]!=0)
        ##rowsum.connenction.matrix[which.connected]
        if (length(which.connected)>1) {
            Ordered<-order(Touched[which.connected])
            which.connected<-which.connected[Ordered]
            cat(which.connected)
        } 
        ## which.connected<-which.connected[which(which.connected!=previous.point)]
        point.now<-which.connected[1]
        if(is.na(point.now)) break
  
        Touched[point.now]<-Touched[point.now]+1
  
        if (any(connection.line==point.now)==FALSE)
            connection.touched.once<-c(connection.touched.once,point.now)
        connection.line<-c(connection.line,point.now)
        connection.matrix[current,point.now]<-0
  
        previous.point<-current
        current<-point.now
  
        rowsum.connection.matrix<-rowSums(connection.matrix)
    }

    ##if (any(Touched<1) ) cat("Fitting line on morphology ------ well done :)")


    cat('Calculation of connection information: ready! \n')


    ##line interpolation az approx-szal
    coords.l<-morpho[connection.line,] #multiple times the same coordinates, as the line touches all the points
    coords<-morpho[connection.touched.once,]

    ##coords<-morpho[connection.line,]
    length<-0
    length.l<-0
    comp.place<-0
    comp.place.l<-0
    for(i in 2:dim(coords)[1]){
        length<-sqrt(sum((coords[i,]-coords[i-1,])^2))
        comp.place<-c(comp.place,comp.place[i-1]+length)
    }
    
    for(i in 2:dim(coords.l)[1]){
        length.l<-sqrt(sum((coords.l[i,]-coords.l[i-1,])^2))
        comp.place.l<-c(comp.place.l,comp.place.l[i-1]+length.l)
    }

    t.param<-comp.place 
    t.param.l<-comp.place.l
    length.fitted.curve<-max(comp.place.l) 


    ts<-seq( from = 0, length.fitted.curve, length=1000 )
    vonal<-apply( coords, 2, function(u) approx( t.param, u, xout = ts ,method="linear")$y ) 
    vonal.l<-apply( coords.l, 2, function(u) approx( t.param.l, u, xout = ts ,method="linear")$y ) 
    vonalfun<-apply( coords, 2, function(u) approxfun( t.param, u,method="linear") )
    vonalfun.l<-apply( coords.l, 2, function(u) approxfun( t.param.l, u,method="linear") )
    fcx<-vonalfun[[1]]
    fcy<-vonalfun[[2]]
    fcz<-vonalfun[[3]]

    ##connections<-read.table('connections.txt')  
    cat('ALMA') 
    ##vonalfun<-apply(coords,2,function(u) spline(t.param,u,xout=ts)$y)
    cell.plot<-scatterplot3d(vonal, type="p", lwd=1,xlab='x (um)',ylab='y (um)', zlab='z (um)', main='Fitted curve on cell morphology',angle=40,highlight.3d=TRUE,pch=20)
    
    png(paste(outname,"/Fitted_curveMorpho.png",sep=""))
    cell.plot<-scatterplot3d(vonal.l, type="l", lwd=3,xlab='x (um)',ylab='y (um)', zlab='z (um)', main='Fitted curve on cell morphology',angle=40)
    dev.off()
    ##cell.plot$points3d(coords,type='h')

    ##lines(cell.plot$xyz.convert(coords))
    ##calculating how many times a segment was touched by the line
    count.touched<-numeric()
    for(ct in 1:max(connection.line)){
        count.touched<-c(count.touched,length(which(connection.line==ct) ) ) 
    }

    ##return(list(connection.line=connection.line,t.param=t.param, length.fitted.curve=length.fitted.curve,vonalfun=vonalfun))
    return(list(connection.line=connection.line,
                t.param=t.param,
                t.param.l=t.param.l,
                length.fitted.curve=length.fitted.curve,
                vonalfun=vonalfun,
                vonalfun.l=vonalfun.l,
                connection.touched.once=connection.touched.once,
                seg.nb=seg.nb,
                count.touched=count.touched))
}

################################################x
###################################################x
################################################################x

branching.result<-connection.function(morpho, connections)
count.touched<-branching.result$count.touched
cat("How many times each segment were touched")
cat(count.touched)
connection.line<-branching.result$connection.line
t.param<-branching.result$t.param
t.param.l<-branching.result$t.param.l
seg.nb<-branching.result$seg.nb

length.fitted.curve<-branching.result$length.fitted.curve
vonalfun<-branching.result$vonalfun #!!!
fcx<-vonalfun[[1]]
fcy<-vonalfun[[2]]
fcz<-vonalfun[[3]]

vonalfun.l<-branching.result$vonalfun.l #!!!
fcx.l<-vonalfun.l[[1]]
fcy.l<-vonalfun.l[[2]]
fcz.l<-vonalfun.l[[3]]


###########################################################
######################################################
############################################

######################################x #memb.currents is optional
ksCSD_all<-function(basis.width.min,
                    basis.width.max,
                    basis.width.step,
                    basis.number.min,
                    basis.number.max,
                    basis.number.step,
                    LFP,
                    elec.kord,
                    memb.currents,
                    seg.length,
                    where2save,M){
    ##LFP<-as.matrix(read.table("myLFP"))
    ##elec.kord<-as.matrix(read.table("elcoord_x_y_z"))
    ##memb.currents<-as.matrix(read.table("membcurr"))
    ##seg.length<-as.matrix(read.table("seglength"))
    ##M<-70
    ##R<-50
    elec.kord<-matrix(elec.kord,ncol=3)
    el.nb<-dim(elec.kord)[1]
    LFPoriginal<-LFP
    
    
    ##
  
  
  
    cat('Entering the loop')
  
    ##
    ##
    
    ##Loop for finding the best parameterss
    sigma<-0.5
    const<-1/(4*pi*sigma)
    ##base<-"gauss" #basis function type
  
    ##to remember the best parameters
    M.best<-numeric()
    R.best<-numeric()
    params<-numeric()
    determinants<-numeric()
    cv.error.best<-Inf #starting value of the best error
    mc.error.best<-Inf
  
    mr.nb<-0 #number of loops
  #####################
    ##Initializing for comparison with "plain" currents
    mc.abs.error.best<-Inf #the previous best solution
    current.error.abs <-numeric()#the current error
    M.best.abs.mc<-numeric()# M in case of the current best soultion
    R.best.abs.mc<-numeric()## R in case of the current best soultion
    T.matrix.abs.mc<-numeric() #best T in case of the current best soultion
    C.best.abs.mc<-numeric()  #C.calc in case of the current best soultion
    current.error.abs.time <-numeric()#error in time in case of best soultion
    
    
    R.all<-2^seq(basis.width.min,basis.width.max,basis.width.step)
    kpower<-seq(basis.number.min,basis.number.max,basis.number.step)
    M.all<-10^kpower ### M atvette a lambda helyet
    errormatrix<-array(0,c(length(M.all),length(R.all)))
    proctimematrix<-array(0,c(length(M.all),length(R.all)))
    for(m in 1:length(M.all)){
        for(r in 1:length(R.all)){
            LFP<-LFPoriginal
            mr.nb<-mr.nb+1
            ##profiling for each R and M
            proctime.m.r<-proc.time()
            lambda<-M.all[m]
            R<-R.all[r]
            
            if(missing(M)) M<-512
            ##base <-"poli"
            ##base<-"gauss"
            ##calculation of the location of basis functions
###################################################################
########################################
############x Basis functions
#########################################
      
      
            ## real number of the basis functions
      
            ##t<-seq(0,max(comp.place[[j]]),length.out=source.branch.db[j])
            source.t<-seq(0,length.fitted.curve-3,length.out=M)
      
            ##where.t<-seq(1,max(comp.place[[j]])-1,length.out=where.branch.db[j])
            ##if the currents are known in specific points
            where.t<-t.param.l #or user defined
            write.table(where.t, paste(outname,'/where_t',sep=''))
            ##vonalfun
            ##cat(vonalfun[[1]])
            ##fcx<-vonalfun[[1]]
            ##fcy<-vonalfun[[2]]
            ##fcz<-vonalfun[[3]]
      
      
      
            ##where.cord <-in whoch point do we want to calculate
            ##where.cord<-matrix(c(fcx(where.t),fcy(where.t),fcz(where.t)),ncol=3)
            where.cord<-matrix(c(fcx.l(where.t),fcy.l(where.t),fcz.l(where.t)),ncol=3)
            source.cord<-matrix(c(fcx.l(source.t),fcy.l(source.t),fcz.l(source.t)),ncol=3)
      ###############################################
            
            where.db<-dim(where.cord)[1]
      
            cat('Position of basis functions determined! \n')
      
      
      
######################x
            ##registerDoMC(cores=4)
################# 
            ##Számoljuk ki a B illetve B.tilda mátrixot
            ##egy sor egy adott i-hez tartozó fgv, oszlopokban azonos helyekhez tartozó
            ##[i,j] : az i.függvény a j-dik helyen
      
            cat(paste('M:',M, 'R:',R, '\n'))
            source.nb<-M
            B.tilda<-array(0,c(source.nb,where.db))
            B<-array(0,c(source.nb,el.nb))
            cat(source.nb, where.db)
            for(i in 1:source.nb){
                Bj.result<-numeric(el.nb)
                ##cat(paste(i,R,source.t,length.fitted.curve))
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
            cat(paste(' Dim of B: ',dim(B)))
            cat(paste(' Dim of B.tilda: ',dim(B.tilda)))
            K.tilda<-t(B)%*%B.tilda
            
            ##Investigating K
            RegMatr<-array(0,c(dim(K)))
            diag(RegMatr)<-lambda
            K<-K+RegMatr
            ##condition number:
            CondEpsilon<-.Machine$double.eps
            KappaFact<-1
            if(kappa(K)>(KappaFact*1/CondEpsilon)){
                El2Ignore<-numeric()
                Ktest<-K
                Electrodes<-1:el.nb
                while(kappa(Ktest)>(KappaFact*1/CondEpsilon)){
                    KcondNb<-numeric()
          
                    for(Which2LeaveChoose in 1: length(Electrodes)){
                        Which2Leave<-Electrodes[Which2LeaveChoose]
                        ##We check again in which case decreases the kappa most an leave out that electrode
                        KcondNb<-c(KcondNb,kappa(K[-c(El2Ignore,Which2Leave),-c(El2Ignore,Which2Leave)]))
            
                    }
          
                    Ktest<-Ktest[-which.min(KcondNb),-which.min(KcondNb)]
                    El2Ignore<-c(El2Ignore,Electrodes[which.min(KcondNb)])
                    Electrodes<-Electrodes[-which.min(KcondNb)]
          
                    ##cat("eliminating the ",which.min(KcondNb),"th Electrode")
                }
                K<-K[-El2Ignore,-El2Ignore]
                K.tilda<-K.tilda[-El2Ignore,]
                LFP<-LFP[-El2Ignore,]
                cat("dim K", dim(K),"\n")
                cat("dim K.tilda", dim(K.tilda),"\n")
                cat("dim LFP", dim(LFP),"\n")
                cat("El to remove:", El2Ignore, "\n")
                ##also remove the affected part of K.tilda
                write.table(El2Ignore,paste(outname,"/El2Ignore_M",M,"_R",R,sep=""))
                if(dim(LFP)[1]<3){
                    cat("Too many electrodes removed for calculation")
                    next	
                }
                
            }
     
      
            ##cat(dim(K))
            ##cat(dim(K.tilda))
            
            ##Let's use oseudoinverse
            
            ##cat("K determinana:", det(K))
            
            cat("kortefa")
            
            cat(is.numeric(K))
            if ( any(is.na(K))){# | round(det(K),100)==0){
                
                cat("Baj van!")
                next
            } 
            cat("Kappa of K:",kappa(K))
            determinants<-c(determinants,kappa(K))
            ##if (round(det(K),20)==0 ) next
      

      
            C.calc<-1/const*t(K.tilda)%*%ginv(K)%*%LFP 
            ##C.calc<-1/const*t(K.tilda)%*%solve(K)%*%LFP 
            write.table(t(K.tilda)%*%ginv(K),paste(outname,"/TransferMatr_M",M,"_R",R,"lambda",lambda,sep=""))
            write.table(K,paste(outname,"/K_M",M,"_R",R,"lambda",lambda,sep=""))
            write.table(K.tilda,paste(outname,"/Ktilda_M",M,"_R",R,"lambda",lambda,sep=""))
            write.table(B,paste(outname,"/B_M",M,"_R",R,"lambda",lambda,sep=""))
            write.table(B.tilda,paste(outname,"/Btilda_M",M,"_R",R,"lambda",lambda,sep=""))
            write.table(C.calc,paste(outname,"/skCSDall_M",M,"_R",R,"lambda",lambda,sep=""))

            cat('\n')
            C.connection.touched.once<-C.calc
            ##Lets put the current in the same order as they were in the simulation and add
            ##up the current regarding the same spot
            cat('The dimension of the matrix ')
            cat(dim(C.calc))
            cat(' the number of segments:')
            cat(seg.nb)
            ##cat ('/n')
            skCSD.all<-array(0,c(seg.nb,dim(LFP)[2]))
            ##!!!
            kimaradt<-numeric()
            if(file.exists(paste(outname,"/sameplace.txt",sep=""))!=TRUE ){
                file.create(paste(outname,"/sameplace.txt",sep=""))
                for(i in 1: seg.nb){
                    sameplace<-numeric()
                    sameplace<-which(branching.result$connection.line==i)
                    ##cat('\n')
	
                    ##f(mr.nb>1 ) {
	  
                    cat(sameplace,file=paste(outname,"/sameplace.txt",sep=""),append=TRUE)
                    cat("\n",file=paste(outname,"/sameplace.txt",sep=""),append=TRUE)
#}
                    ##reorder the results, that the first row will be the current on the first segment etc...
                    if(length(sameplace)==0){ cat(i, "\n")
                        kimaradt<-c(kimaradt,i)
                    }
                    if(length(sameplace)==1)
                        skCSD.all[i,]<-C.calc[sameplace,]
                    if(length(sameplace)>1)
                        skCSD.all[i,]<-as.matrix(colSums(C.calc[c(sameplace),]),nrow=1)
                }
            }       
            for(i in 1: seg.nb){
                sameplace<-numeric()
                sameplace<-which(branching.result$connection.line==i)

                ##reorder the results, that the first row will be the current on the first segment etc...
                if(length(sameplace)==0){ cat(i, "\n")
                    kimaradt<-c(kimaradt,i)
                }
                if(length(sameplace)==1)
                    skCSD.all[i,]<-C.calc[sameplace,]
                if(length(sameplace)>1)
                    skCSD.all[i,]<-as.matrix(colSums(C.calc[c(sameplace),]),nrow=1)
            }
    
            ##The value should be devided by how many time that segment was touched
            ## skCSD.all[i,]<-skCSD.all[i,]/length(sameplace) ##!!!
            ##}
            C.calc<-skCSD.all
            ##write.table(C.calc,paste(outname,"/Ccalc_M",M,"_R",R,sep="")) ##might be too big data

            ##C.connection.touched.once<-C.calc
            ##C.calc<-C.calc[c(branching.result$connection.touched.once),]
            images.name<-paste(outname,'/image_M',M,'_R',R,"lambda",lambda,'.png',sep='')
            png(images.name)
            ##par(mfrow=c(1,2))
            ##image(C.connection.touched.once,col=rainbow(200))
            rangemax.C.calc<-max(abs(C.calc))
            image(C.calc,col=rainbow(200),zlim=c(-rangemax.C.calc,rangemax.C.calc))
            dev.off()
######################################
            ##Is this result better than the previous ones
            
####################################
######## if we know the membrane currents
#######################################x
            ##if(exists("memb.currents") | exists("seg.length")){
            
            szamol<-1
            if(szamol==1){
                
                if(memb.currents[1,1]!="NA" | seg.length!="NA"){
                    ##it is possible to compare the original and the  calculated currents
                    ##Estimation error of the currents
                    
                    ##memb.currents<-as.matrix(read.table('membcurr'))
                    
                    funaramvonal<-function(x) x/seg.length
                    memb.currents.vonal<-apply(memb.currents,2,funaramvonal) 
                    

                    ##         if( mc.abs.error.best > current.error.abs){
                    ##           M.best.abs.mc<-M 
                    ##           R.best.abs.mc<-R
                    ##           T.matrix.abs.mc<-t(K.tilda)%*%ginv(K)
                    ##           mc.abs.error.best<- current.error.abs
                    ##           C.best.abs.mc<-C.calc
                    ##           current.error.abs.time.best<-current.error.abs.time
                    ##           
                    ##           ##writing it to file
                    ##         }
                    
        
        
        
                    write.table(memb.currents.vonal,paste(outname,'/membcurr_line',sep=''),col.names=FALSE, row.names=FALSE)
        
    
        
        
                    ##Smoothed membrane currents
        
                    ##Rsmooth<-R#1/(sqrt(2)*R)
                    Rsmooth<-length.fitted.curve/dim(LFP)[1]
                    write.table(c(length.fitted.curve/dim(LFP)[1],length.fitted.curve,  dim(LFP)[1]) ,paste(outname, '/SmoothinKernelWidth', sep=''), col.names=FALSE, row.names=FALSE)
        
        
                } }##current known
      
            ##printing out the profiling result and writing it to file
            proc.time()-proctime.m.r
            proctimematrix[m,r]<-(proc.time()-proctime.m.r)[1]
            
        }} #for m and r
  
    proc.name<-paste(outname,"/proctime",sep="")
    colnames(proctimematrix)<-R.all
    rownames(proctimematrix)<-M.all
    write.table(proctimematrix,proc.name)
    ##writing things to file
    ##different error values
    ##the best parameters in case of comparison with the smoothed currents
    
    ##where did we calculate the currents?
    write.table(where.cord, paste(outname,'/where_cord',sep=''),col.names=FALSE, row.names=FALSE)
  
    ##write to files the determinants of matrix K
  
    write.table(determinants, paste(outname,'/determinants',sep=''),col.names=FALSE, row.names=FALSE)
  
    ##how do the points follow eachother
    write.table(branching.result$connection.touched.once, paste(outname,'/connection_touched_once',sep=''),col.names=FALSE, row.names=FALSE)
  
    
} #big ksCSD function






