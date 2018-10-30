
plot_N_map<-function(cur.t,N,Grid,highlight=NULL,leg.title="Abundance"){
  require(rgeos)
  require(ggplot2)
  library(RColorBrewer)
  myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
  Tmp=Grid[[1]]
  if(is.null(highlight)==FALSE){
    midpoints=data.frame(gCentroid(Tmp[highlight,],byid=TRUE))
    colnames(midpoints)=c("Easting","Northing")
  }
  Abundance=N[,cur.t]
  Cur.df=cbind(data.frame(gCentroid(Tmp,byid=TRUE)),Abundance)
  new.colnames=colnames(Cur.df)
  new.colnames[1:2]=c("Easting","Northing")
  colnames(Cur.df)=new.colnames
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  p1=ggplot(Cur.df)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_gradientn(colours=myPalette(100),name=leg.title)
  if(is.null(highlight)==FALSE){
    #p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067,xmax=Easting,ymin=Northing,ymax=Northing+25067))
    p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
  }
  p1
}

plot_N_map_xy<-function(N,XY,leg.title="Abundance"){
  require(ggplot2)
  library(RColorBrewer)
  myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
  Abundance=N
  Cur.df=cbind(data.frame(Easting=XY[,1],Northing=XY[,2],Abundance))
  colnames(Cur.df)=c("Easting","Northing","Abundance")
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  p1=ggplot(Cur.df)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_gradientn(colours=myPalette(100),name=leg.title)
  p1
}

  

