simulate <- function(DataFolder) {
  #DataFolder<-"simulation/BS_d50el8_CosChanging/" #where is a simulation of a BS cell to get coordinates ets.
  #  cat(DataFolder)
  setwd(DataFolder)
  #morpho<-as.matrix(read.table('segcoordinates.txt'))
  
  morpho <-
    matrix(as.numeric(as.matrix(read.table(
      'coordsmid_x_y_z'
    ))), ncol = 3)
  plot(morpho[, 1], morpho[, 2], col = "WHITE")#,xlim=c(-500,-200))
  text(morpho[, 1], morpho[, 2], paste(1:length(morpho[, 1])), cex = 0.7)
  
  #  membcurr<-matrix(as.numeric(as.matrix(read.table('simulation/cell_1/membcurr'))))
  connections <- as.matrix(read.table('connections.txt'))
  
  NewFolder <- paste(getwd(),"skCSDreconstruct/",sep="/") #where to save
  dir.create(NewFolder)
  FilesInFolder <- list.files(DataFolder)
  
  file.copy(paste(getwd(), FilesInFolder, sep = "/"), NewFolder)
  
  
  
  segmid <- matrix(as.matrix(read.table('coordsmid_x_y_z'), ncol = 3))
  seglength <- as.matrix(read.table('seglength'))
  SpFrNb <- 30
  MembMatr <- array(0, c(length(seglength), SpFrNb))
  #fr<-1
  
  #Csak adott reszre nem 0 aramokat rakni
  #for(fr  in 1:50){
  selStim <- 1:dim(seglength)[1]
  for (fr  in 1:SpFrNb) {
    plot(cumsum(seglength[selStim]), cos(cumsum(seglength[selStim]) / max(cumsum(seglength[selStim])) *
                                           2 * pi * fr))
    
    MembCline <-
      cos(cumsum(seglength[selStim]) / max(cumsum(seglength[selStim])) * pi *
            fr + pi / 2) #*2
    
    Membcurr <- MembCline#/seglength[selStim]
    MembMatr[selStim, fr] <- Membcurr
  }
  image(t(MembMatr), col = rainbow(100))
  matplot(MembMatr, t = "l", xlim = c(149, 180))
  
  seg.nb <- dim(segmid)[1]
  elec <-
    matrix(as.matrix(read.table('elcoord_x_y_z', colClasses = 'numeric')), ncol =
             3)
  
  
  Dmatr <- array(0,  c(length(seglength), dim(elec)[1]))
  for (i1 in 1:length(seglength)) {
    for (i2 in 1:dim(elec)[1]) {
      Dmatr[i1, i2] <- sqrt(sum((segmid[i1, ] - elec[i2, ]) ^ 2))
    }
    
  }
  
  
  Tmatr <- 1 / (4 * pi * 0.3) * 1 / Dmatr
  LFP <- t(Tmatr) %*% MembMatr
  
  
 
  
  write.table(LFP, "myLFP", row.names = FALSE, col.names = FALSE)
  write.table(MembMatr,
              "membcurr",
              row.names = FALSE,
              col.names = FALSE)
  write.table(1:50, "time", row.names = FALSE, col.names = FALSE)
}