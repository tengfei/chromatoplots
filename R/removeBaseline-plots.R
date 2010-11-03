setMethod("explore", c("cpSample", protocolClass("removeBaseline")),
          function(object, protocol, raw,geom='l')
          {
            gg_vbox <- gg_baseline_box(object,raw=raw,geom=geom)
            stack_plots("Baseline Subtraction",
              list(gg_vbox),
              )
         
          })

## gg_baseline_box to contain the ggobi device
gg_baseline_box <- function(object,raw,gg=ggobi(),geom=geom){
  base_da <- gtkDrawingArea()
  asCairoDevice(base_da)
  dev <- -1
  gSignalConnect(base_da,'realize',function(wid) dev<<-dev.cur())
  dev_fun <- function() dev
  gg <- gg_baseline(object,raw,gg=gg,dev=dev_fun,geom=geom)
  disp <- display(gg[1],embed=TRUE)
  variables(disp) <- list(X='mz')
  imode(disp) <- "Identify"
  vbox <- gtkVBox()
  vbox$add(disp)
  vbox$add(base_da)
  vbox
}


gg_baseline <- function(object,raw,gg=ggobi(),dev=dev.new(),geom=geom){
  mz_range <- object@mzrange
  mz_range <- mz_range[1]:mz_range[2]
  mz_range <- as.data.frame(mz_range)
 # rownames(mz_range) <- seq_len(nrow(mz_range))
  rownames(mz_range) <- seq_len(nrow(mz_range))-1+object@mzrange[1]
  colnames(mz_range) <- 'mz'
  gg['mz'] <- mz_range
  d <- gg['mz']
  gSignalConnect(gg,'identify-point',function(gg,plot,id,dataset){
    if (id == -1)
      return()
    if (!(as.RGtkObject(d) == dataset))
      return()
    if (is.function(dev))
      dev.set(dev())
    else dev.set(dev)
    if(id>-1) mz <- mz_range$mz[id+1]
    baseplot(object, raw=raw, mz=mz, subtract = TRUE,geom=geom)
  })
  gg
}


baseplot <- function(object,raw, mz, subtract = TRUE,geom=geom)
{
  ## FIXME: replace the code with faster qtpaint function
  time <- raw@scantime
  mzmin <- raw@mzrange[1]
  rawint <- raw@env$profile[mz-mzmin+1,]
  chrom <- data.frame(time = time, intensity = rawint)
  #p <- qplot(time,rawint,data=chrom,geom="point")
  #p$title=NULL
  corint <- object@env$profile[mz-mzmin+1,]
  if (subtract) {
    res <- corint
    fit <- rawint - corint
  } else {
    res <- rawint - corint
    fit <- corint
  }
  ylim=range(rawint)
  par(mfrow=c(2,1),mar=c(2,2,1,1))
  plot(time,rawint,ylab='intensity',pch=20,ylim=ylim,type=geom)
  lines(time,fit,col='gray50')
  plot(time,res,ylab='residuals',pch=20,ylim=ylim,type=geom)
}

## cplot by ggplot2 used for plotting single slice at certain mz
setMethod('cplot',c("cpSample", protocolClass("removeBaseline")),
          function(obj,protocol,raw,mz,ylim=NULL,geom=c('p','l')){
  baseplot(obj,raw=raw,mz=mz,geom=geom)
})






