geoscalePhylo<-function(tree, ages, upwards=FALSE, units=c("Period", "Epoch", "Age"), boxes="Age", tick.scale=boxes, user.scale, cex.age=0.3, cex.ts=0.3, cex.tip=0.3, width=1, ts.col=TRUE, ranges=FALSE, vers="ICS2013", x.lim, quat.rm=FALSE,...){
    
  if(missing(ages) == TRUE && ranges == TRUE){
    cat("\n ages file is not provided, ranges is to FALSE")
      ranges <- FALSE
  }
  
  if(missing(user.scale) & any(units == "User")){
    units <- units[units != "User"]
    cat("\n user.scale not provided, 'Other' removed from units.")
  }
  
  if(missing(ages) == FALSE){
    ages<-ages[tree$tip.label,]    
  }
  
  if(any(units == "User") & !missing(user.scale)){
    
    Midpoint <- matrix(ncol=1,nrow=length(user.scale[,1]))
      Midpoint[,1] <- (user.scale[,"Start"] + user.scale[,"End"])/2
        user.scale <- cbind(user.scale,Midpoint)  
  }
  
  if(length(units) == 1){
  	ts.width=0.15
  } else if(length(units) == 2){
  	ts.width=0.2
  } else if(length(units) >= 3){
  	ts.width=0.25
  }
  
  # Standardizing the names of temporal units
  
  units[units == "Eonothem"] <- "Eon"
  units[units == "Erathem"] <- "Era"
  units[units == "Series"] <- "Epoch"
  units[units == "System"] <- "Period"  
  units[units == "Stage"] <- "Age"
    units <- unique(units)
      
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
      if(quat.rm == TRUE){
       timescale[(timescale[,"Midpoint"] < 3),"Name"] <- NA
      }
  
  tscale.data<-matrix(ncol=3,nrow=6)
    colnames(tscale.data) <-c("srt","Depth","size")
    rownames(tscale.data) <-c("Eon","Era","Period","Epoch","Age","User")
      if(upwards == TRUE){
        tscale.data[,"srt"] <- c(90,90,90,0,0,0)
      } else tscale.data[,"srt"] <- c(0,0,0,90,90,90)
      
      tscale.data[,"Depth"] <- c(1,1,1,2,3.5,3.5)
      tscale.data[,"size"] <- c(1,1,1,0.8,0.8,0.8)

  ## GEOLOGICAL RANGES OF TAXA

  units<-rownames(tscale.data)[sort(match(units,rownames(tscale.data)),decreasing=T)] 
  
   if(is.null(tree$root.time)){     
      return(cat("\n tree$root.time is missing, check tree is time scaled."))
   } else (root.age <- tree$root.time)
  
  
  # if x.lim set by user -  
  
  if(!missing(x.lim)){
    x.lim <- sort(root.age - x.lim)
  } else if(ranges == TRUE && !missing(ages) && missing(x.lim)){
    x.lim <- (root.age - min(ages)) + diff(range(ages))*0.05
  } else {
    x.lim <- NULL
  }
  
# if(missing(x.lim) == TRUE){
 #     x.lim <- NULL
 #   } else(x.lim=sort(root.age - x.lim))
  
    timescale<-timescale[order(timescale[,1],decreasing=T),]
    	timescale.rescaled <- timescale
  		timescale.rescaled[,c("Start","End","Midpoint")] <- root.age - timescale[,c("Start","End","Midpoint")]
	
    offset<-array(dim=length(tree$tip.label),data=1)

    if(ranges==TRUE && missing(ages) == FALSE){
      
      offset.correction <- diff(range(ages)) * 0.01 
      taxon.ranges <- root.age - ages[,c("FAD","LAD")]
       offset <- array(dim=length(ages[,"FAD"]),data=(ages[,"FAD"] - ages[,"LAD"])+offset.correction)
    }

  ### ADDING A TIMESCALE

  if(tick.scale != "n" | tick.scale != "no"){
    if(tick.scale == "myr" | is.numeric(tick.scale)){    
      scale.ticks=1
      
      if(is.numeric(tick.scale)){
       scale.ages <- tick.scale
      } else {scale.ages=10}
      
        ticks <- root.age - seq(0,4600,scale.ticks)
          time <- seq(0,4600,scale.ages)
          time2 <- root.age - time
            lwd<-c(1,0.5,0.5,0.5,0.5,0.7,0.5,0.5,0.5,0.5)
            col<-c("black","grey","grey","grey","grey")
    }

    if(tick.scale != "myr" & is.numeric(tick.scale) == FALSE){
      time<-subset(timescale,timescale[,"Type"] == tick.scale & timescale[,"Source"] == "ICS")
        time<-sort(unique(c(time[,1],time[,2])))
        	time2<-ticks <- root.age - time
            lwd=1
              col="black"
    }
  }
  
  # Plotting the tree vertically
  if(upwards == TRUE){
    
    # ADDING THE TIME SCALE
    par(fig=c(0,ts.width,0,1))
     par(mar=c(3,1,2,5))
    
        plot.phylo(tree,plot=FALSE,no.margin=T,y.lim=x.lim,direction="upwards")
          
          timescale.rescaled.names <- timescale.rescaled
           timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,"End"] > par()$usr[3],]
            timescale.rescaled.names[timescale.rescaled.names[,"Start"] < par()$usr[3],"Start"] <- par()$usr[3]
    
          if(min(timescale.rescaled.names[,"End"]) < par()$usr[4]){
            timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,"Start"] < par()$usr[4],]
             timescale.rescaled.names[timescale.rescaled.names[,"End"] > par()$usr[4],"End"] <- par()$usr[4]
          }
                timescale.rescaled.names[,"Midpoint"] <-(timescale.rescaled.names[,"Start"] + timescale.rescaled.names[,"End"])/2
    
    
          unit.depths <- tscale.data[units,"Depth"]
            if(tick.scale == "n" | tick.scale == "no"){
              unit.depths <- c(unit.depths,0.5)          
            } else if(length(units) <= 3){
              unit.depths <- c(unit.depths,2)
            } else if(length(units) > 3) {
              unit.depths <- c(unit.depths,2)
            }

           unit.depths <- cumsum(unit.depths/sum(unit.depths))
            unit.depths<-c(par()$usr[2],par()$usr[2]-(unit.depths*(par()$usr[2]-par()$usr[1])))
             
         depth <- unit.depths[length(unit.depths) - 1] - unit.depths[length(unit.depths)]
    
        if(tick.scale != "n" && tick.scale != "no"){
          text((unit.depths[length(unit.depths)]+depth*0.3),time2,time,cex=cex.age,srt=0)           
           segments((unit.depths[length(unit.depths)-1]),ticks,(unit.depths[length(unit.depths)]+depth*0.75),ticks,lwd=lwd,col=col)
        }
    
    for(t in 1:length(units)){
      
      if(units[t] == "User"){
        tscale<-user.scale
         tscale[,c("Start","End","Midpoint")] <- root.age - tscale[,c("Start","End","Midpoint")]
         tscale.names <- tscale
      } else {
        tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == units[t])
        tscale.names<-subset(timescale.rescaled.names,timescale.rescaled.names[,"Type"] == units[t])
      }
          
      if(ts.col == TRUE & units[t] != "User"){
        rect(unit.depths[t],tscale[,"Start"],unit.depths[t+1],tscale[,"End"],col=rgb(tscale[,"Col_R"],tscale[,"Col_G"],tscale[,"Col_B"],maxColorValue=255))
      } else rect(unit.depths[t],tscale[,"Start"],unit.depths[t+1],tscale[,"End"],col="white")
      text((unit.depths[t] + unit.depths[t+1])/2,tscale.names[,"Midpoint"],tscale.names[,"Name"],cex=cex.ts*tscale.data[match(units[t],rownames(tscale.data)),"size"],srt=tscale.data[match(units[t],rownames(tscale.data)),"srt"])
    }
    
    # ADDING THE PHYLOGENY
    par(fig=c((ts.width)+0.001,1,0,1),new=T)
     par(mar=c(3,0,2,2))
    
      plot.phylo(tree,plot=FALSE,no.margin=T,y.lim=x.lim,direction="upwards")
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    
      if (! missing(boxes)){
       if(boxes == "User"){
         tscale<-user.scale
       } else {
         tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == boxes)} 
        rect(par()$usr[3],tscale[,"Start"],par()$usr[4],tscale[,"End"],col=c("grey90","white"),border=NA)}
        
      par(fig=c(ts.width,1,0,1),new=T)
        par(mar=c(3,0,2,2))
    
        plot.phylo(tree,label.offset=offset,edge.width=width,no.margin=T,y.lim=x.lim,cex=cex.tip,direction="upwards",...)
          if (ranges == TRUE){
            segments(lastPP$xx[c(1:length(tree$tip.label))],taxon.ranges[,"FAD"],lastPP$xx[c(1:length(tree$tip.label))],taxon.ranges[,"LAD"],col="black",lwd=width*2)
          }
    
  } else{
    # Plotting the tree horizontally
    
  # PLOT 1 - TIMESCALE
    par(fig=c(0,1,0,(ts.width)+0.001))
    par(mar=c(1,3,0,2))

	    plot.phylo(tree,plot=FALSE,no.margin=T,x.lim=x.lim,direction="rightwards")
    
        timescale.rescaled.names <- timescale.rescaled
         timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,"End"] > par()$usr[1],]
          timescale.rescaled.names[timescale.rescaled.names[,"Start"] < par()$usr[1],"Start"] <- par()$usr[1]
    
        if(min(timescale.rescaled.names[,"End"]) < par()$usr[2]){
          timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,"Start"] < par()$usr[2],]
           timescale.rescaled.names[timescale.rescaled.names[,"End"] > par()$usr[2],"End"] <- par()$usr[2]
        }
            timescale.rescaled.names[,"Midpoint"] <-(timescale.rescaled.names[,"Start"] + timescale.rescaled.names[,"End"])/2
    
       unit.depths <- tscale.data[units,"Depth"]
        if(tick.scale == "n" & tick.scale == "no"){
          unit.depths <- c(unit.depths,0.5)          
        } else if(length(units) <= 3){
          unit.depths <- c(unit.depths,2)
        } else if(length(units) > 3) {
          unit.depths <- c(unit.depths,2)
        }
     
         unit.depths <- cumsum(unit.depths/sum(unit.depths))
          unit.depths<-c(par()$usr[4],par()$usr[4]-(unit.depths*(par()$usr[4]-par()$usr[3])))
    
       depth <- unit.depths[length(unit.depths) - 1] - unit.depths[length(unit.depths)]
    
    if(tick.scale != "n" && tick.scale != "no"){
      text(time2,(unit.depths[length(unit.depths)]+depth*0.3),time, cex=cex.age,srt=90)
  	  segments(ticks,(unit.depths[length(unit.depths)-1]),ticks,(unit.depths[length(unit.depths)]+depth*0.6),lwd=lwd,col=col)
    }

    for(t in 1:length(units)){
 	
 	    if(units[t] == "User"){
 	    	tscale<-user.scale
  	    	tscale[,c("Start","End","Midpoint")] <- root.age - tscale[,c("Start","End","Midpoint")]
     	    	tscale.names <- tscale
      } else {
 	      tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == units[t])
 	      tscale.names<-subset(timescale.rescaled.names,timescale.rescaled.names[,"Type"] == units[t])
 	    }
 	 	
 	    if(ts.col == TRUE & units[t] != "User"){rect(tscale[,"Start"],unit.depths[t],tscale[,"End"],unit.depths[t+1],col=rgb(tscale[,"Col_R"],tscale[,"Col_G"],tscale[,"Col_B"],maxColorValue=255))
 		    } else rect(tscale[,"Start"],unit.depths[t],tscale[,"End"],unit.depths[t+1],col="white")
          text(tscale.names[,"Midpoint"],(unit.depths[t] + unit.depths[t+1])/2,tscale.names[,"Name"],cex=cex.ts*tscale.data[match(units[t],rownames(tscale.data)),"size"],srt=tscale.data[match(units[t],rownames(tscale.data)),"srt"])
    }
 
  ## PLOT 2: PHYLOGENY

    par(fig=c(0,1,(ts.width)+0.001,1),new=T)
    par(mar=c(0,3,2,2))
  
      plot.phylo(tree,plot=FALSE,no.margin=T,x.lim=x.lim)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
      if (!missing(boxes)){
        if(boxes == "User"){
          tscale<-user.scale
        } else {tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == boxes)} 
  	        rect(tscale[,"Start"],par()$usr[3],tscale[,"End"],par()$usr[4],col=c("grey90","white"),border=NA)}

    par(fig=c(0,1,(ts.width)+0.001,1),new=T)
    par(mar=c(0,3,2,2))
  
      plot.phylo(tree,label.offset=offset,edge.width=width,no.margin=T,x.lim=x.lim,cex=cex.tip,...)
        if (ranges == TRUE){
          segments(taxon.ranges[,"FAD"],lastPP$yy[c(1:length(tree$tip.label))],taxon.ranges[,"LAD"],lastPP$yy[c(1:length(tree$tip.label))],col="black",lwd=width*2)
        }
  }
}