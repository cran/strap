geoscalePhylo<-function(tree,ages,cex.age=0.3,cex.ts=0.3,cex.tip=0.3,width=1,ranges=FALSE,scale=boxes,units=c("Age","Epoch","Period"),ts.col=TRUE,boxes="Age",nscale,x.lim,vers="ICS2013",...){

  if(missing(ages) == TRUE && ranges == TRUE){
    cat("\n ages file is not provided, ranges is to FALSE")
      ranges <- FALSE
  }
  
  if(missing(ages) == FALSE){
    ages<-ages[tree$tip.label,]    
  }
  
  if(ranges == TRUE && missing(ages) == FALSE){
     missing.tip.names <- setdiff(tree$tip.label,row.names(ages))
      if(length(missing.tip.names) > 0){
                
        cat(paste("\n",missing.tip.names,"not present in ages file, ranges set to FALSE"))
        cat("\n ranges set to FALSE")
          ranges <- FALSE
        
      }    
  }
   
  timescales <- NULL
  data(timescales,envir=environment())
  
  timescale <- timescales[[vers]]
  
  tscale_data<-matrix(ncol=3,nrow=6)
    colnames(tscale_data) <-c("srt","Depth","size")
    rownames(tscale_data) <-c("Eon","Era","Period","Epoch","Age","Other")
      tscale_data[,1] <- c(0,0,0,90,90,90)
      tscale_data[,2] <- c(0.5,0.5,0.5,1,1,1)
      tscale_data[,3] <- c(1,1,1,0.8,0.7,0.5)

  ## GEOLOGICAL RANGES OF TAXA

  units<-rownames(tscale_data)[sort(match(units,rownames(tscale_data)),decreasing=T)] 
  	
  p.dist <- dist.nodes(tree)[, length(tree$tip.label) + 1]

   if(is.null(tree$root.time)){     
      return(cat("\n tree$root.time is missing, check tree is time scaled."))
   } else (root.age <- tree$root.time)
  
    if(missing(x.lim) == TRUE){x.lim <- NULL
    } else(x.lim=sort(root.age - x.lim))

    timescale<-timescale[order(timescale[,1],decreasing=T),]
    	timescale.rescaled <- timescale
  		timescale.rescaled[,c(1,2,4)] <- root.age - timescale[,c(1,2,4)]
	
    offset<-array(dim=length(tree$tip.label),data=1)

    if(ranges==TRUE && missing(ages) == FALSE){
      
      offset.correction <- diff(range(ages)) * 0.01 
      taxon.ranges <- root.age - ages[,c("FAD","LAD")]
       offset <- array(dim=length(ages[,"FAD"]),data=(ages[,"FAD"] - ages[,"LAD"])+offset.correction)
    }

  ### ADDING A TIMESCALE

  if(scale != "n" | scale != "no"){
    if(scale=="myr"){    
      scale1=1
      scale2=10
        ticks <- root.age - seq(0,4600,scale1)
          time <- seq(0,4600,scale2)
          time2 <- root.age - time
            lwd<-c(1,0.5,0.5,0.5,0.5,0.7,0.5,0.5,0.5,0.5)
            col<-c("black","grey","grey","grey","grey")}

    if(scale != "myr"){
      time<-subset(timescale,timescale[,"Type"] == scale & timescale[,"Source"] == "ICS")
        time<-sort(unique(c(time[,1],time[,2])))
        	time2<-ticks <- root.age - time
            lwd=1
              col="black"
    }
  }
  
  # PLOT 1 - TIMESCALE
    par(fig=c(0,1,0,0.2))
    par(mar=c(1,3,0,2))

	    plot.phylo(tree,plot=FALSE,no.margin=T,x.lim=x.lim)

        vals<-tscale_data[units,"Depth"]
          vals<-c(vals,0.8)
            vals<-cumsum(vals/sum(vals))
              vals<-c(par()$usr[4],par()$usr[4]-(vals*(par()$usr[4]-par()$usr[3])))
                val<-vals[length(vals)-1] - vals[length(vals)]

    if(scale != "n" && scale != "no"){
      text(time2,(vals[length(vals)]+val*0.3),time, cex=cex.age,srt=90)
  	 segments(ticks,(vals[length(vals)-1]),ticks,(vals[length(vals)]+val*0.75),lwd=lwd,col=col)
    }

    for(t in 1:length(units)){
 	
 	    if(units[t] == "Other"){
 	    	tscale_n<-nscale
 	    } else tscale_n<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == units[t])
 	
 	 	
 	    if(ts.col == TRUE & units[t] != "Other"){rect(tscale_n[,"Start"],vals[t],tscale_n[,"End"],vals[t+1],col=rgb(tscale_n[,"Col_R"],tscale_n[,"Col_G"],tscale_n[,"Col_B"],maxColorValue=255))
 		    } else rect(tscale_n[,"Start"],vals[t],tscale_n[,"End"],vals[t+1],col="white")
          text(tscale_n[,"Midpoint"],(vals[t] + vals[t+1])/2,tscale_n[,"Name"],cex=cex.ts*tscale_data[match(units[t],rownames(tscale_data)),"size"],srt=tscale_data[match(units[t],rownames(tscale_data)),"srt"])
    }
 
  ## PLOT 2: PHYLOGENY

    par(fig=c(0,1,0.2,1),new=T)
    par(mar=c(0,3,2,2))
  
      plot.phylo(tree,plot=FALSE,no.margin=T,x.lim=x.lim)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
      if (! missing(boxes)){
        if(boxes == "Other"){
          tscale_x<-nscale
        } else {tscale_x<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == boxes)} 
  	        rect(tscale_x[,"Start"],par()$usr[3],tscale_x[,"End"],par()$usr[4],col=c("grey90","white"),border=NA)}


    par(fig=c(0,1,0.2,1),new=T)
    par(mar=c(0,3,2,2))
  
      plot.phylo(tree,label.offset=offset,edge.width=width,no.margin=T,x.lim=x.lim,cex=cex.tip,...)
        if (ranges == TRUE){
          segments(taxon.ranges[,"FAD"],lastPP$yy[c(1:length(tree$tip.label))],taxon.ranges[,"LAD"],lastPP$yy[c(1:length(tree$tip.label))],col="black",lwd=width*2)
        }
  }
