## cplotViewBL2 by basic graphci, so could be easy to embed into animation
## cplotViewBL2 <- function(object,protocal,raw,mz,xlim=NULL,ylim=NULL,type='l',subtract=T){
##   time <- raw@scantime
##   mzmin <- raw@mzrange[1]
##   rawint <- raw@env$profile[mz-mzmin+1,]
##   corint <- object@env$profile[mz-mzmin+1,]
##   if (subtract) {
##     res <- corint
##     fit <- rawint - corint
##   } else {
##     res <- rawint - corint
##     fit <- corint
##   }
##   if(is.null(xlim)){xlim=range(time)}
##   if(is.null(ylim)){ylim=c(min(res),max(rawint))}
##   par(mfrow=c(2,1),mar=c(5,2,1,1))
##   plot(x=time,y=rawint,xlim=xlim,ylim=ylim,type=type,pch=20)
##   lines(time,fit,col='gray50',lwd=2)
##   plot(time,res,ylab='residuals',xlim=xlim,ylim=ylim,type=type,pch=20)
## }


## cplotAniBL <- function(raw,by=10,ylim=c(0,6000),size=c(800,300),mzrange=NULL,outdir=NULL){
##   require(animation)
##   if(is.null(mzrange)){
##     mzrange <- raw@mzrange
##   }
##   mz <- seq(mzrange[1],mzrange[2],by=by)
##   n <- length(mz)
##   if(is.null(outdir)){
##     ani.options(nmax=n+1,ani.width=size[1],ani.height=size[2])
##   }else{
##     ani.options(nmax=n+1,ani.width=size[1],ani.height=size[2],outdir=outdir)
##   }
##   ani.start()
##   for(i in mz){
##     cplotViewBL(raw,mz=i,ylim=ylim)
##   }
##   ani.stop()
## }

## cplotAniBL2 <- function(object,raw,by=10,ylim=c(0,6000),size=c(800,300),mzrange=NULL,outdir=NULL){
##   require(animation)
##   if(is.null(mzrange)){
##     mzrange <- raw@mzrange
##   }
##   mz <- seq(mzrange[1],mzrange[2],by=by)
##   n <- length(mz)
##   if(is.null(outdir)){
##     ani.options(nmax=n+1,ani.width=size[1],ani.height=size[2])
##   }else{
##     ani.options(nmax=n+1,ani.width=size[1],ani.height=size[2],outdir=outdir)
##   }
##   ani.start()
##   for(i in mz){
##     cplotViewBL2(object=object,raw=raw,mz=i,xlim=xlim,ylim=ylim)
##   }
##   ani.stop()
## }

## #cplotAniBL(raw_prof,ylim=c(0,2000),mzrange=c(50,80),by=1,outdir='~/Desktop/by1_50_70')

## cplotAniPK <- function(object,raw,by=10,size=c(800,300),xlim=NULL,mzrange=NULL,outdir=NULL){
##   require(animation)
##   if(is.null(xlim)) {xlim <- range(raw@scantime)}
##   if(is.null(mzrange)){
##     mzrange <- raw@mzrange
##   }
##   mz <- seq(mzrange[1],mzrange[2],by=by)
##   n <- length(mz)
##   if(is.null(outdir)){
##     ani.options(nmax=n+1,ani.width=size[1],ani.height=size[2])
##   }else{
##     ani.options(nmax=n+1,ani.width=size[1],ani.height=size[2],outdir=outdir)
##   }
##   ani.start()
##   for(i in mz){
##     cplotPeaks2(obj=object,raw=raw,mz=i,rt_range=xlim)
##   }
##   ani.stop()
## }

                                        #cplotAniPK(peaks,raw1,outdir='~/Desktop/pk/')
