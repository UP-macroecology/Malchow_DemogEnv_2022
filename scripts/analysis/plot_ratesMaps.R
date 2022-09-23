# require packages
require(raster)
require(rasterVis)
require(ggplot2)
library(colorspace)
require(ggnewscale)

# function for multiple plots
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# function for generating one figure
generate.figure <- function (file) {
  print(file)
  specID <- as.integer(strsplit(strsplit(unlist(strsplit(file,"/"))[5],"_")[[1]][2],"spec")[[1]][2])
  #specID <- as.integer(strsplit(unlist(strsplit(file,"/"))[5],"_")[[1]][1])
  print(specID)
  load(file)
  load("data/hillshade/hill_ch.Rdata")
  hill_ch <- projectRaster(hill_ch,growthr_stack,method = 'ngb')
  hill_ch <- mask(hill_ch, growthr_stack$rInf)
  # hill_ch <- setExtent(hill_ch, extent(growthr_stack))
  growthr_stack <- addLayer(growthr_stack, hill_ch$layer)
  fec_max <- c(NA,2.0,2.0,1.0,1.0,NA,2.0,2.0,3.0,5.0)[specID]
  gr_max  <- 2.0
  
  #todo: adapt font sizes
  
  plot.fec=gplot(growthr_stack$fec)+ #note: gplot() from rasterVis
  	geom_raster(data=as.data.frame(growthr_stack$layer, xy=T),aes(x=x, y=y,fill=layer),alpha=1.0)+
  	scale_fill_continuous(guide="none",low="#00000059", high="#FFFFFF59", na.value="white")+
  	new_scale_fill()+
    geom_raster(aes(fill=ifelse(value < fec_max, value, fec_max)),alpha=0.75)+
    scale_fill_viridis_c(limits = c(0,fec_max), breaks=seq(0,round(fec_max),length.out=3), direction = 1, option = "viridis",na.value="white")+
    theme_void() +
    labs(x="",y="", fill="fec.", title="(b)")+
    coord_equal(expand=F)
  
  plot.jSv=gplot(growthr_stack$jSv)+ #note: gplot() from rasterVis
  	geom_raster(data=as.data.frame(growthr_stack$layer, xy=T),aes(x=x, y=y,fill=layer),alpha=1.0)+
  	scale_fill_continuous(guide="none",low="#00000059", high="#FFFFFF59", na.value="white")+
  	new_scale_fill()+
    geom_raster(aes(fill=value),alpha=0.75)+
    scale_fill_viridis_c(limits = c(0, 1.0),breaks=c(0,0.5,1.0), direction = 1, option = "viridis",na.value="white")+
    theme_void() +
    labs(x="",y="", fill="juv. surv.", title="(a)")+
    coord_equal(expand=F)
  
  plot.aSv=gplot(growthr_stack$aSv)+ #note: gplot() from rasterVis
  	geom_raster(data=as.data.frame(growthr_stack$layer, xy=T),aes(x=x, y=y,fill=layer),alpha=1.0)+
  	scale_fill_continuous(guide="none",low="#00000059", high="#FFFFFF59", na.value="white")+
  	new_scale_fill()+
    geom_raster(aes(fill=value),alpha=0.75)+
    scale_fill_viridis_c(limits = c(0, 1.0), breaks=c(0,0.5,1.0), direction = 1, option = "viridis",na.value="white")+
  	#scale_fill_gradientn(colours=rev(cols_sequential), breaks = c(0,0.5,1.0), labels = c("0","0.5","1.0"), na.value = "white") +
    theme_void() +
    labs(x="",y="", fill="ad. surv.", title="(c)")+
    coord_equal(expand=F)
  
  plot.rInf=gplot(growthr_stack$rInf)+ #note: gplot() from rasterVis
    geom_raster(data=as.data.frame(growthr_stack$layer, xy=T),aes(x=x, y=y,fill=layer),alpha=1.0)+
    scale_fill_continuous(guide="none",low="#00000059", high="#FFFFFF59", na.value="white")+
    new_scale_fill()+
    geom_raster(aes(fill=ifelse(value < gr_max, value, gr_max)),alpha=0.75)+
    scale_fill_continuous_divergingx(palette = 'RdBu', limits = c(0, gr_max), breaks=c(0,1,gr_max), mid = 1.0, na.value="white") +
    theme_void() +
    labs(x="",y="", fill="growth rate", title="(d)")+
    coord_equal(expand=F)
  
  # adapt plot size if needed
  pdf(paste("results/analysis/rates_maps/maps/pdf/", unlist(strsplit(unlist(strsplit(file,"/"))[5], ".Rdata")),".pdf",sep=""), width=12, height=8)
  #png(paste("results/analysis/rates_maps/maps/png/", unlist(strsplit(unlist(strsplit(file,"/"))[5], ".Rdata")),".png",sep=""), width=24, height=16, type = "cairo",units='cm', res = 200)
	multiplot(plot.jSv, plot.aSv, plot.fec, plot.rInf, cols=2)
	Sys.sleep(0.1)
  dev.off()
  
}

# loop over all rasterstacks
for (file in list.files(path="results/analysis/rates_maps/stacks", pattern="Rdata", full.names = T)[-1]){
  generate.figure(file)
}

# test
generate.figure(file = "results/analysis/rates_maps/stacks/stack_spec10_year2016.Rdata")
