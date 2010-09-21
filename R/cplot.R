######################################################################
## quantity build ##
######################################################################
input_quantity <- function(object){
  peaks <- as.data.frame(object@peaks)
  lst <- by(peaks,peaks[,"sample"],function(peaks){
    tapply(peaks[,"maxf"],peaks[,"comp"],sum)
  })
  quantity <- unlist(lst)
  object@comps <- cbind(object@comps,quantity=quantity)
  object
}

######################################################################
## TIC ##
######################################################################
cplotTIC <- function(object,xintercept=NULL,xscale=NULL,color=NULL){

  if(is.null(color)){
    color <- c("red","blue","green")
  }
  if(!("quantity" %in% names(object@comps))){
    object <- input_quantity(object)
  }
  comps <- as.data.frame(object@comps)
  if(is.null(xscale)){
    xscale=range(comps[,'rt'])
  }
  p <- ggplot(comps,aes(x=rt,y=quantity,xend=rt,yend=0,group=sample))
  p <- p+geom_segment()+facet_grid(sample~.)
  if(is.null(xintercept)){
    print(p)
  }else{
    for(i in 1:length(xintercept)){
      p <- p+geom_vline(xintercept=xintercept[[i]],colour=color[i],alpha=0.5)
      cat(color[i]," : ",xintercept[[i]],"\n")
    }
    
    print(p)
  }
}


######################################################################
## RT ##
######################################################################
cplotRT <- function(object,xscale=NA,geom=NA,sample=NA,log=F){
  path <- object@filepaths
  tic <- list()
  sample.names <- rownames(object@phenoData)
  if(is.na(sample)) {
    sample=unique(object@comps[,'sample'])
  }
  sample=sample.names[sample]
  for(i in 1:length(path)){
    raw <- loadSample(path[i])
    ##     raw <- genProfile(raw)
    ##     prof <- raw@env$profile
    ##     tic[[i]] <- apply(prof,2,sum)
    tic[[i]] <- raw@tic
  }
  rt <- object@rt
  if(is.na(xscale)){
    xscale <- range(unlist(rt))
  }
  df2 <- data.frame()
  for(i in 1:length(rt)){
    x <- unlist(rt[[i]])
    df <- data.frame(rt=x,tic=unlist(tic),sample=rep(sample.names,each=length(rt[[1]][[1]])),
                     type=names(rt)[i],level=rep(1:length(rt[[1]]),each=length(rt[[1]][[1]])))
    df2 <- rbind(df2,df)
  }
  df2 <- df2[df2$sample%in%sample,]
  df2$sample <- factor(df2$sample)
  if(is.na(geom)){
    if(log){
      p<- qplot(rt,log(tic+1),group=sample,colour=sample,data=df2,geom="path",facets=type~.)+
        scale_x_continuous(limits=xscale)+
          scale_y_continuous(limits=range(df2[which(df2$rt<xscale[2]&df2$rt>xscale[1]),'tic']))+
            opts(legend.position="none")
    }else{
      p<- qplot(rt,tic,group=sample,colour=sample,data=df2,geom="path",facets=type~.)+
        scale_x_continuous(limits=xscale)+
          scale_y_continuous(limits=range(df2[which(df2$rt<xscale[2]&df2$rt>xscale[1]),'tic']))+
            opts(legend.position="none")
    }
    
  }else{
    if(geom=="heatmap"){
      if(log==F){
        p <- qplot(x=rt,y=level,xend=rt,yend=level-0.9,data=df2,colour=tic,geom="segment",facets=type~.,ylab="sample")+
          scale_colour_gradient(low="yellow",high="red")+
            scale_x_continuous(limits=xscale)+
              scale_y_continuous(breaks=1:length(rt[[1]]))+
                opts(legend.position="none")
      }else{
        p <- qplot(x=rt,y=level,xend=rt,yend=level-0.9,data=df2,colour=log(tic+1),geom="segment",facets=type~.,ylab="sample")+
          scale_colour_gradient(low="yellow",high="red")+
            scale_x_continuous(limits=xscale)+
              scale_y_continuous(breaks=1:length(rt[[1]]))+
                opts(legend.position="none")
      }
    }
  }
  p
}



######################################################################
## plot spectrum
######################################################################

cplotSpec <- function(object,sample=NA,comp=NA,group=NA,...){
  comps <- as.data.frame(object@comps)
  peaks <- as.data.frame(object@peaks)
  comps[,'inter'] <- as.character(interaction(comps[,'sample'],comps[,'comp']))
  peaks[,'inter'] <- as.character(interaction(peaks[,'sample'],peaks[,'comp']))
  if(!is.na(group))inter <- comps[comps[,'group']==group,'inter']
  if(is.na(group)) inter <- comps[comps[,'sample']==sample&
                                  comps[,'comp']==comp,'inter']
  idx <- which(peaks[,'inter']%in%inter)
  peaks <- peaks[idx,]
  par(mfrow=c(length(inter),1),mar=rep(0.5,4),oma=c(5,5,4,1),mgp=c(1.9,0,0))
  inter=sort(inter)
  for(i in 1:length(inter)){
    interSpec(object,inter=inter[i])
  }
}





interSpec <- function(object,inter=NA){
  comps <- as.data.frame(object@comps)
  peaks <- as.data.frame(object@peaks)
  comps[,'inter'] <- as.character(interaction(comps[,'sample'],comps[,'comp']))
  peaks[,'inter'] <- as.character(interaction(peaks[,'sample'],peaks[,'comp']))
  if(is.na(inter)) stop("specify a inter")
  idx <- which(peaks[,'inter']==inter)
  peaks <- peaks[idx,]
  plot(x=peaks[,'mz'],y=2*peaks[,'maxf']*sqrt(pi/2)*peaks[,'sigma'],
       xlab="mz",ylab="intensity",type="h",xlim=range(object@peaks[,'mz']))
  mtext(inter,side=3,line=-3,cex=3,col="gray")
}


######################################################################
## cplotRtFit
######################################################################
cplotRtFit <- function(object,raw,xscale=NA,sample=NA){
  sample.names <- rownames(object@phenoData)
  if(is.na(sample)) sample=unique(object@comps[,'sample'])
  s <- object@comps[,'sample']%in%sample
  retcor_d <- data.frame(raw=raw@comps[s,'rt'],
                         cor=object@comps[s,'rt'],
                         sample=sample.names[object@comps[s,'sample']])
  retcor_d$sample <- factor(retcor_d$sample)
  if(is.na(xscale)) xscale=range(c(retcor_d$raw,retcor_d$cor))
  p_retcor_fit <- qplot(raw,cor-raw,data=retcor_d,facets=sample~.)+
    scale_x_continuous(limits=xscale)+geom_smooth()
  p_retcor_fit
}


######################################################################
## cplotResult ##
######################################################################
## FIX ME: HOW TO FIND STH IN COMMON MORE PRECISELY
cplotResult <- function(object,id.diff=NULL,id.common=NULL,group.name=F){
  final <- as(object,"ExpressionSet")
  mxt <- exprs(final)
  cl <- ncol(mxt)
  if(group.name==F){
    if(!is.null(id.diff)){
      mx <- mxt[rownames(mxt)%in%id.diff,]
      mx <- cbind(mx,diff=apply(mx,1,function(x){
        mean(x[1:cl/2])-mean(x[(cl/2+1):cl])
      })
                  )
      mx1 <- mx[mx[,'diff']>0,]
      mx2 <- mx[mx[,'diff']<0,]
      mx1 <- mx1[,-(ncol(mx))]
      df1 <- data.frame(level=as.numeric(mx1))
      df1$sample <- rep(colnames(mx1),each=nrow(mx1))
      df1$id <- rep(rownames(mx1),ncol(mx1))
      df1$group <- 'A.High'
      mx2 <- mx2[,-(ncol(mx))]
      df2 <- data.frame(level=as.numeric(mx2))
      df2$sample <- rep(colnames(mx2),each=nrow(mx2))
      df2$id <- rep(rownames(mx2),ncol(mx2))
      df2$group <- 'B.High'
      df <- rbind(df1,df2)
      if(!is.null(id.common)){
        mx <- mxt[rownames(mxt)%in%id.common,]
        df3 <- data.frame(level=as.numeric(mx))
        df3$sample <- rep(colnames(mx),each=nrow(mx))
        df3$id <- rep(rownames(mx),ncol(mx))
        df3$group <- rep('Same level',nrow(mx))
        df <- rbind(df,df3)
      }
      p <- qplot(sample,level,data=df,colour=group,geom="line",group=id)
      p
    }
  }
  if(group.name==T){
    id <- c(id.diff,id.common)
    id <- id[!is.null(id)]
    mx <- mxt[rownames(mxt)%in%id,]
    if(length(id)==1) mx <- t(mx)
    df <- data.frame(level=as.numeric(mx))
    df$sample <- rep(colnames(mx),each=nrow(mx))
    df$group <- rep(id,ncol(mx))
    if(length(id)>1) p <- qplot(sample,level,data=df,colour=factor(group),geom="line",group=group)
    if(length(id)==1) p <- qplot(sample,level,data=df,geom="line",group=group)
  }
  p
}


######################################################################
## view aninmation based on parameter ##
######################################################################

finalObj <- function(object,tic.cutoff,npeaks.cutoff,dist.cutoff,rt_window){
  require(chromatoplots)
  require(Biobase)
  xset_comps_filt <- findComps(object,"sigma_filt",tic.cutoff=tic.cutoff,npeaks.cutoff=npeaks.cutoff)
  xset_groups <- groupComps(xset_comps_filt, "angle",dist.cutoff=dist.cutoff,rt_window=rt_window)
  xset_sum <- summarize(xset_groups, "common")
  xset_norm <- normalize(xset_sum, "scale")
  xset_norm
}

getPvalue <- function(object){
  require(limma)
  final <- as(object,"ExpressionSet")
  design <- final@phenoData@data
  df <- exprs(final)
  design <- model.matrix(~class,data=design)
  fit <- lmFit(df,design)
  fit <- eBayes(fit)
  df <- cbind(df,pvalue=fit$p.value[,2])
  df <- as.data.frame(df)
  df
}


## dir <- "~/Desktop/data/"
## cplotPvalue(dir,xlim=c(1,20))


cplotPvalue <- function(dir,xlim=NA,output=NA){
  require(animation)
  if(is.na(output)){
    ani.options(nmax = length(list.files(dir))+1)
  }else{
    ani.options(nmax = length(list.files(dir))+1,outdir=output)
  }
  ani.start()
  for(nm in list.files(dir)){
    file <- file.path(dir,nm)
    df <- read.csv(file,header=T)
    df <- df[order(df$pvalue),]
    df2 <- apply(df[,1:8],1,sum)
    summ <- sort(df2,decreasing=T)
    if(is.na(xlim)) xlim=c(1,nrow(df))
    xrange=c(min(xlim),min(nrow(df),max(xlim)))
    par(mfrow=c(2,1),mar=c(4,2,1,1))
    plot(df$pvalue,type='l',col="blue",xlim=xlim,ylab="p-value",
         ylim=range(df$pvalue[min(xrange):max(xrange)]))
    abline(v=seq(min(xlim),max(xlim),by=round(diff(xlim)/10)),
           lty="dotted",col="lightgray")
    mtext(nm,side=3,line=-6,cex=3,col="gray")
    plot(summ,type="l",col="red",xlim=xlim,ylab="Sum for each row",
         ylim=range(summ[min(xrange):max(xrange)]))
    abline(v=seq(min(xlim),max(xlim),by=round(diff(xlim)/10)),lty="dotted",col="lightgray")
    mtext(nm,side=3,line=-6,cex=3,col="gray")
  }
  ani.stop()
}

cplotHist <- function(dir){
  require(ggplot2)
  t <- numeric()
  for(nm in list.files(dir)){
    file <- file.path(dir,nm)
    df <- read.csv(file,header=T)
    top <- nrow(subset(df,df$pvalue<0.0001))
    t <- c(t,top)
  }
  d <- data.frame(t=t)
  qplot(data=d,x=t,xlab="Number of metabolites[pvlaue<0.00001]",ylab="count",geom="histogram")
                                        #      geom_text(x=max(t),y=max(table(t),label=names(table(t))[table(t)==as.numeric(max(table(t)))]))
}

## plot sigma_d vs sigma_t
cplotDt <- function(object){
  require(ggplot2)
  group <- as.data.frame(object@groups)
  p <- ggplot(data=group,aes(x=sigma_t,y=sigma_d))+geom_point()
  p
}


getFeatures <- function(object){
  final <- as(object,"ExpressionSet")
  df <- getPvalue(object)
  design <- object@phenoData
  f <- factor(design[,'class'])
  p <- df$pvalue
  mx <- exprs(final)
  diff <- apply(mx,1,function(mx){
    sp <- split(mx,f)
    diff <- mean(sp[[1]])-mean(sp[[2]])
  }
                )
  feat <- final@featureData@data
  feat$diff <- diff
  feat$pvalue <- p
  feat
}

## find top different metabolite
## make it return group id
topDiff <- function(object,scale=NULL,pvalue=NULL){
  ob <- getFeatures(object)
  if(is.null(scale)) scale <- 1:nrow(ob)
  ob <- ob[order(ob$pvalue,decreasing=F),][scale,]
  if(!is.null(pvalue)) 
    ob <- ob[ob$pvalue<pvalue,]
  ob
}

## object from this could be used for dump2msp
## find top common metabolite
topCommon <- function(object,scale=NULL,method="diff"){
  final <- as(object,"ExpressionSet")
  design <- object@phenoData
  ob <- getFeatures(object)
  if(is.null(scale)) scale=1:nrow(ob)
  ob$absdiff <- abs(ob$diff)
  mx <- exprs(final)
  ob$sum <- apply(mx,1,sum)
  if(method=="diff") idx <- order(ob$absdiff,decreasing=F)
  if(method=="sum") idx <- order(ob$sum,decreasing=F)
  ob <- ob[idx,][scale,]
  ob
}

cplotTotalP <- function(dir,method="diff",xscale=NA){
  require(ggplot2)
  final <- data.frame()
  if(method=="diff"){
    for(nm in list.files(dir)){
      file <- file.path(dir,nm)
      df <- read.csv(file,header=T)
      p <- df$pvalue
      this <- data.frame(pvalue= sort(p),
                         id=seq_len(nrow(df)),
                         group = rep(nm,nrow(df)))
      final <- rbind(final,this)
    }
    p <- qplot(x=id,y=pvalue,data=final,geom="line",colour=group,group=group)+
      opts(legend.position="none")
    if(!is.na(xscale)) {
      sub <- subset(final,id<max(xscale))
      yscale <- range(sub$pvalue)
      p <- p+scale_x_continuous(limits=xscale)+scale_y_continuous(limits=yscale)
    }}

  if(method=="sum"){
    for(nm in list.files(dir)){
      file <- file.path(dir,nm)
      df <- read.csv(file,header=T)
      df <- df[,1:8]
      df$sum <- apply(df,1,sum)
      this <- data.frame(sum =sort(df$sum,decreasing=T),
                         id=seq_len(nrow(df)),
                         group=rep(nm,nrow(df)))
      final <- rbind(final,this)
    }
    
    p <- qplot(x=id,y=sum,data=final,geom="line",colour=group,group=group)+
      opts(legend.position="none")
    if(!is.na(xscale)) {
      sub <- subset(final,id<max(xscale))
      yscale <- range(sub$sum)
      p <- p+scale_x_continuous(limits=xscale)+scale_y_continuous(limits=yscale)
    }}
  print(p)
}


######################################################################
## retention time ##
######################################################################
                                        # function rtmatch
rtmatch=function(old,new,mode,noise,range.min=-0.1,range.max=0.1,plot.match=T){
  cat("Retention time unit: minute\n\n")
  old=sort(old)
  n=0
  x=numeric()
  y=numeric()
  for(i in 1:length(old)){
    z=0
    cat("old rt(",old[i],") matches: ")
    for(k in 1:length(new)){
      if(new[k]<=(old[i]+range.max)&new[k]>=(old[i]+range.min)){
        x=old[i]
        y=c(y,new[k])
        cat(new[k]," ")
        z=1
      }       
    }
    if(z==1)n=n+1
    cat("\n")  
  }
  cat("\n")
  percent=n/length(old)*100
  cat("Total match:",percent,"%(",n,"in",length(old),") of top",length(old),"\n")
  d.old=old
  d.new=new
  if(plot.match){
    par(mfrow=c(1,1),mar=c(4,2,1,1))
    ddd <- c(d.old,d.new)
    d.min <- min(ddd)
    d.max <- max(ddd)
    plot(c(d.min,d.max),c(1,2),xlim=c(d.min,d.max),xlab="RT",ylim=c(0.5,2.5),ylab="",
         type="n",axes=F,main=paste("match range: [",range.min,",",range.max,"]","Total match:",percent,"%(",n,"in",length(old),") of top",length(old)))
    rect(0,0,d.max+5,3,col="grey90")
    abline(v=seq(0,d.max+5,0.25),col="white")
    abline(v=seq(0,d.max+5,1),col="white",lwd=2)
    abline(h=1:2,col="white")
    points(d.old,pch=16,jitter(rep(1,length(d.old)),4))
    points(d.new,pch=16,jitter(rep(2,length(d.new)),4))
    axis(side=1,at=seq(0,d.max+5,5),labels=seq(0,d.max+5,5))
    axis(side=2,at=1:2,labels=c("old","new"))
    box(col="grey80")
    if(plot.match==T){
      abline(v=mode,col="red") 
      abline(v=noise,col="blue")
    }}
  return(percent)
}

## Directly read raw dir, and compare tic, used to compare RT
cplotTIC3 <- function(dir,xscale=NULL,facets=FALSE){
  lst <- list.files(dir,full.names=T)
  prof.lst <- list()
  rt.lst <- list()
  for(i in lst){
    x <- loadSample(i)
    x <- genProfile(x)
    prof.lst[[i]] <- x@env$profile
    rt.lst[[i]] <- x@scantime
  }
  lst <- lapply(prof.lst,function(x){
    apply(x,2,sum)
  })
  df <- cbind(rt=unlist(rt.lst),int=unlist(lst))
  df <- cbind(df,sample=rep(1:length(lst),each=length(rt.lst[[1]])))
  df <- as.data.frame(df)
  p <- qplot(rt,int,data=df,group=factor(sample),geom='line',colour=factor(sample))+ylab('intensity')
  if(facets) p <- p+facet_grid(sample~.)
  if(is.null(xscale)) {print(p)} else {
    p <- p+scale_x_continuous(limits=xscale)
    print(p)
  }
}


