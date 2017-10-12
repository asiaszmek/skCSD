#Triangular grid generator
#19-μm hexagonal center-to-center pitch


TriSide<-38 #unit, um whatever
Xmin<--100
Ymin<--100
Zdist<-50 #Cell to electrode Distance
ColNb<-12 #number of rows, Even numb
RowNb<-6 #number of columns, Even numb
TriHeight<-TriSide#*sin(60)
TriX1<-Xmin + TriSide*(1:(ColNb/2))-TriSide
TriX2<-Xmin-TriSide/2 + TriSide*(1:(ColNb/2))
TriY1<-Ymin + TriHeight*(1:(RowNb)) - TriHeight
TriY2<-Ymin - TriHeight/2 + TriHeight*(1:(RowNb))

grid1<-expand.grid(TriX1,TriY1)
grid2<-expand.grid(TriX2,TriY2)

Xcord<-c(grid1[[1]] ,grid2[[1]])
Ycord<-c(grid1[[2]] ,grid2[[2]])
plot(Xcord,Ycord,asp=1)
elccord<-c(Xcord, Ycord, rep(0,length(Ycord)))

write.table(elccord,"elcoord_x_y_z")
