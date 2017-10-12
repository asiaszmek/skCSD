#BLUE and RED colour scale
library(fields)
library("RColorBrewer")

ColoursDori<-function(Xmatr){
  ExtValue<-max(abs(Xmatr))
  #ColPalette<-designer.colors( 256, rainbow(5,start=0,end=4/6,s=0.7), x= c( -ExtValue, -ExtValue/2,0, ExtValue/2,ExtValue))
  #ColPalette<-designer.colors( 256, col=c("RED","WHITE","BLUE"), x= seq( -1, 1,,3),alpha=1)#alpha=0.7)
  #
  #Col2Use<-colorRampPalette(c("RED","WHITE","BLUE"))( 64 )
  #ColPalette<-designer.colors( 256, col=Col2Use, x= seq( 0, 1,,3))
  #ColPalette<-designer.colors( 256, col=c(colors()[556:552],"WHITE",colors()[562:566]), x= seq( 0, 1,,11))
 
  #ColPalette<-designer.colors( 256, col=c("RED","WHITE","BLUE"), x= seq( 0, 1,,3))
  #ColPalette<-colorRampPalette(c("red", "white", "blue"),   space = "Lab")(200)
  #two.colors(n=256, start="BLUE", end="red", middle="white", x= c( -ExtValue,0,ExtValue))
  Cols<-brewer.pal(11, "RdBu")
  ColPalette<-colorRampPalette(c(Cols), space = "Lab")(200)
  return(list(ColPalette,ExtValue))
}


ColoursDoriLFP<-function(Xmatr){
  ExtValue<-max(abs(Xmatr))
  Cols<-brewer.pal(11,"PiYG")
  #ColPalette<-designer.colors( 256, col=c("PURPLE","WHITE","GREEN"), x= seq( 0, 1,,3))
  #ColPalette<-colorRampPalette(c("purple", "white", "green"), space = "Lab")(200)
  ColPalette<-colorRampPalette(c(Cols), space = "Lab")(200)
#ColPalette<-designer.colors( 256, col=c(colors()[551:547], "WHITE",colors()[254:258]), x= seq( 0, 1,,11))
return(list(ColPalette,ExtValue))
}

# 
# 
# colorRampPalette
# 
# 
# 
# colorRamp(c("red", "green"))( (0:4)/4 ) ## (x) , x in [0,1]
#      colorRampPalette(c("blue", "red"))( 4 ) ## (n)
#      ## a ramp in opacity of blue values
#      colorRampPalette(c(rgb(0,0,1,1), rgb(0,0,1,0)), alpha = TRUE)(8)
#      
#      require(graphics)
#      
#      ## Here space="rgb" gives palettes that vary only in saturation,
#      ## as intended.
#      ## With space="Lab" the steps are more uniform, but the hues
#      ## are slightly purple.
#      filled.contour(volcano,
#                     color.palette =
#                         colorRampPalette(c("red", "white", "blue")),
#                     asp = 1)
#      filled.contour(volcano,
#                     color.palette =
#                         colorRampPalette(c("red", "white", "blue"),
#                                     space = "Lab"),
#                     asp = 1)
#      
#      ## Interpolating a 'sequential' ColorBrewer palette
#      YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
#      filled.contour(volcano,
#                     color.palette = colorRampPalette(YlOrBr, space = "Lab"),
#                     asp = 1)
#      filled.contour(volcano,
#                     color.palette = colorRampPalette(YlOrBr, space = "Lab",
#                                                      bias = 0.5),
#                     asp = 1)
#      
#      ## 'jet.colors' is "as in Matlab"
#      ## (and hurting the eyes by over-saturation)
#      jet.colors <-
#        colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
#                           "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#      filled.contour(volcano, color = jet.colors, asp = 1)
#      ## space="Lab" helps when colors don't form a natural sequence
#      m <- outer(1:20,1:20,function(x,y) sin(sqrt(x*y)/3))
#      rgb.palette <- colorRampPalette(c("red", "orange", "blue"),
#                                      space = "rgb")
#      Lab.palette <- colorRampPalette(c("red", "orange", "blue"),
#                                      space = "Lab")
#      filled.contour(m, col = rgb.palette(20))
#      filled.contour(m, col = Lab.palette(20))
# 
