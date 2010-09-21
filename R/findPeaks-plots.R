plot_peak <- function(object,id, raw, cutoff, sample=NULL,residuals = FALSE,island=FALSE)
{
  if(!is.null(sample)){
    peaks <- object@peaks[object@peaks[,'sample']==sample,]
  }else{
    peaks <- object@.Data
  }
  peak <- peaks[id,]
  t =c(peak['rtmin'],peak['rtmax'])
  range <- profRange(raw,m = peak["mz"],rtrange=t)
  prof <- raw@env$profile
  mz <- peak['mz']
  df <- as.data.frame(peak)
  y <- prof[range$massidx, range$scanidx]

  cutoff <- cutoff[range$massidx, range$scanidx]
  x <- raw@scantime[range$scanidx]
  fitted <- egh(x, peak["rt"], peak["maxf"], peak["sigma"], peak["tau"])
  ymax <- max(y,peak['maxf'])
  if ((residuals|island)&!(residuals&island))
    par(mfrow=c(2,1),mar=rep(3,4),mgp=c(2,1,0))
  if(residuals&island)
    par(mfrow=c(3,1),mar=rep(3,4),oma=c(5,5,4,1),mgp=c(2,1,0))
  plot(x, y, xlab = "time", ylab = "intensity",ylim=c(0,ymax))
  lines(x, fitted, col = "red", lwd = 2)
  lines(x, cutoff, col = "blue")
  points(peak['rt'],peak['maxf'],col="blue",pch=19,cex=1.5)
  if (residuals)
    {
      if (!any(is.na(fitted)))
        plot(x, y - fitted, xlab="time", ylab = "residuals")
      abline(h=0)
      lines(x, cutoff, col = "blue")
      
    }
  if(island)   cplotPeaks2(object,raw,mz)
}

## the infamous GGobi/R peak visualization
gg_peaks <- function(object, raw, cutoff, sample=NULL,residuals = TRUE, island=TRUE, gg = ggobi(),dev = dev.new())
{
  
  if(!is.null(sample))
    peaks <- object@peaks[object@peaks[,'sample']==sample,]
  else
    peaks <- object@.Data
  peaks_df <- as.data.frame(peaks)
  rownames(peaks_df) <- seq_len(nrow(peaks_df))
  gg["peaks"] <- peaks_df
  d <- gg["peaks"]
  
  gSignalConnect(gg, "identify-point", function(gg, plot, id, dataset) {
    if (id == -1)
      return()
    if (!(as.RGtkObject(d) == dataset))
      return()
    if (is.function(dev))
      dev.set(dev())
    else dev.set(dev)
    plot_peak(object,id=id+1,raw=raw, cutoff=cutoff, sample=sample,residuals=residuals,island=island)
  })

  gg
}

gg_peaks_box <- function(object, raw, cutoff, sample=NULL ,residuals = TRUE,island=TRUE, gg=ggobi(), ...)
{
  ## the GGobi peak visualization
  
  
  peak_da <- gtkDrawingArea()
  asCairoDevice(peak_da)

  ## lazily determine device ID after realization
  dev <- -1
  gSignalConnect(peak_da, "realize", function(wid) dev <<- dev.cur())
  dev_fun <- function() dev
  
### FIXME: how to pass the device to gg_peaks?
  gg <- gg_peaks(object,raw,cutoff,sample,residuals,island,gg,dev_fun)
  disp <- display(gg[1], embed = TRUE)
  variables(disp) <- list(X = "rt", Y = "mz")
  imode(disp) <- "Identify"
  vbox <- gtkVBox()
  peak_hbox <- gtkHBox()
  peak_hbox$add(disp)
  peak_hbox$add(peak_da)
  peak_hbox
}

## profile image, chromatogram with cutoff, and the GGobi visualization
setMethod("explore", c("cpPeaks", "ProtoFindPeaksGauss"),
          function(object, protocol, raw, sample=NULL, residuals = TRUE, island=TRUE, gg = ggobi(),
                   dev = dev.cur())
          {
            
            ## FIXME: Use qt later, the profile image
            ## p_prof <- image(raw)
            ## the chromatogram
            
            prof <- raw@env$profile
            quant <- quantile(prof, protocol@alpha, na.rm = TRUE)
            ## p_chrom_max <- cplotPeaks(object,raw,mz)
            peak_hbox <- gg_peaks_box(object, raw=raw,
                                      cutoff=matrix(quant, nrow(prof), ncol(prof)),
                                      sample=sample,
                                      residuals=residuals,island=island)
            
            


            stack <- stack_plots(paste("Peak Detection"),
                                 list(peak_hbox))
            ## build the visualization
            ## stack <- stack_plots(paste("Peak Detection  ","(MZ = ",mz,")  "),
            ##                                  list(p_prof, p_chrom_max, peak_hbox))

            stack
          })

setMethod("explore", c("cpPeaks", "ProtoFindPeaksParabola"),
          function(object, protocol, raw, sample=NULL, residuals = TRUE, island=TRUE, gg = ggobi(),
                   dev = dev.cur())
          {
            
            ## FIXME: Use qt later, the profile image
            ## p_prof <- image(raw)
            ## the chromatogram
            
            prof <- raw@env$profile
            quant <- quantile(prof, protocol@alpha, na.rm = TRUE)
            ## p_chrom_max <- cplotPeaks(object,raw,mz)
            peak_hbox <- gg_peaks_box(object, raw=raw,
                                      cutoff=matrix(quant, nrow(prof), ncol(prof)),
                                      sample=sample,
                                      residuals=residuals,island=island)
            
            


            stack <- stack_plots(paste("Peak Detection"),
                                 list(peak_hbox))
            ## build the visualization
            ## stack <- stack_plots(paste("Peak Detection  ","(MZ = ",mz,")  "),
            ##                                  list(p_prof, p_chrom_max, peak_hbox))

            stack
          })


######################################################################
## peaks: low level plot functions
######################################################################
cplotPeaks <- function(obj,raw,mz,ylim=NULL){
  prof <- raw@env$profile
  mzmin <- raw@mzrange[1]
  p_chrom_res <- qplot(raw@scantime,prof[mz-mzmin+1,],xlab="scantime",ylab="residuals")
  if(class(obj)=="cpPeaks")
    alpha <- obj@pipeline@.Data[[4]]@alpha
  else
    alpha <- obj@pipeline@.Data[[1]]@pipeline@.Data[[4]]@alpha
  if(!is.null(alpha)){
    quant <- quantile(prof,alpha,na.rm=TRUE)
  }
  peaks <- peaks[peaks[,"mz"]==mz,]
  p_chrom_max <- p_chrom_res+geom_point(size=4,colour="blue",data=data.frame(peaks),aes(x=rt,y=maxf))
  ifelse(is.null(alpha),p_chrom_max,p_chrom_max <- p_chrom_max+geom_hline(yintercept=quant,colour="red",size=1))
  if(!is.null(ylim)) p_chrom_max <- p_chrom_max+scale_y_continuous(limits=ylim)
  p_chrom_max$title <- NULL
  p_chrom_max
}

## use base R graphic device, faster than ggplot
cplotPeaks2 <- function(obj,raw,mz,rt_range=NULL){
  prof <- raw@env$profile
  mzmin <- raw@mzrange[1]
  if(class(obj)=="cpPeaks")
    {alpha <- obj@pipeline@.Data[[4]]@alpha
     peaks <- obj}
  else{
    alpha <- obj@pipeline@.Data[[1]]@pipeline@.Data[[4]]@alpha
    peaks <- obj@peaks}
  quant <- quantile(prof,alpha,na.rm=TRUE)
  peaks <- peaks[peaks[,"mz"]==mz,]
  df <- as.data.frame(peaks)
  if(!is.null(rt_range)) {
    df <- df[df$rt<max(rt_range)&df$rt>min(rt_range),]
    idx <- which(raw@scantime>min(rt_range)&raw@scantime<max(rt_range))
    x <- raw@scantime[idx]
    y <- prof[mz-mzmin+1,idx]
  }else{
    x <- raw@scantime
    y <- prof[mz-mzmin+1,]}
  plot(x,y,xlab="scantime",ylab="intensity",pch=19,cex=1,main=paste('mz= ',mz),ylim=c(0,max(df$maxf,y)))
  points(df$rt,df$maxf,col="blue",pch=19,cex=1.5)
  abline(h=quant,col="red",lwd=2)
}


## plot peaks max
setGeneric('cplotPeaksMax',function(object,protocal,...) standardGeneric('cplotPeaksMax'))
setMethod('cplotPeaksMax',c("cpPeaks"),
          function(object,protocal,sample=NULL,chain=FALSE,...){
            if(!is.null(sample)){
              peaks <- object@peaks[object@peaks[,'sample']==sample,]
            }else{
              peaks <- object@.Data
            }
            dd <- data.frame(peaks)
            p <- ggplot(data=dd,aes(x=rt,y=mz,color=log(maxf)))+scale_color_gradient(low='yellow',high='black')+geom_point()+opts(legend.position='none')
            p
          })
## for ProtoFindComps
## plot peaks max
setMethod('cplotPeaksMax',c("xcmsSet"),
          function(object,protocal,sample=1,chain=TRUE,...){
            if(!is.null(sample)){
              peaks <- object@peaks[object@peaks[,'sample']==sample,]
            }else{
              peaks <- object@.Data
            }
            dd <- data.frame(peaks)
            p <- ggplot(data=dd,aes(x=rt,y=mz,color=log(maxf)))+scale_color_gradient(low='yellow',high='black')+geom_point()+opts(legend.position='none')
            ifelse(chain,p <- p+geom_line(aes(group=comp),color='black',alpha=0.5),p <- p)
            p
          })


