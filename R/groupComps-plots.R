setMethod("explore", c("xcmsSet", protocolClass("groupComps")),
          function(object, protocol, gg = ggobi()) 
          {
            groups_df <- as.data.frame(object@groups)
            groups_df$ncomps <- factor(groups_df$ncomps)
            
            group_id <- levels(factor(object@comps[,"group"]))
            ##group_id <- rep(as.numeric(levels(group_f)), table(group_f))
            groups_df <- cbind(group = group_id, groups_df)
            gg["groups"] <- groups_df
            gg_groups <- display(gg["groups"], embed = TRUE)
            imode(gg_groups) <- "Identify"
            variables(gg_groups) <- list(X = "sigma_t", Y = "mu_d")

            gg_comp_chain <- gg_all_comps_disp(object@comps, gg)
            
### FIXME: need dataset with vars for m/z and rows for component
### FIXME: exclusion specification and plot variables will change to show groups
            
            comps <- object@comps
            comps <- comps[order(comps[,"comp"]),]
            comps <- comps[order(comps[,"sample"]),]
            spectra <- t(find_spectra(object@peaks))
            ##spectra <- log(spectra + 1)
            colnames(spectra) <- paste("mz", seq_len(ncol(spectra)), sep="")
            gg["spectra"] <- as.data.frame(spectra)
            spectra_d <- gg["spectra"]
            gg_comp_peaks <- display(spectra_d, "Parallel Coordinates Display",
                                     list(X=colnames(spectra)[1]), embed = TRUE)
            ##variables(gg_comp_peaks) <- list(X=colnames(spectra)[1])

            ##group <- 17 # group to show in bottom plot
            ##spectra_g <- spectra[comps[,"group"] == group,]
            ##non_empty <- which(colSums(spectra_g) > 0)
            ##spectra_g <- spectra_g[,non_empty]
            ##colnames(spectra_g) <- paste("mz", non_empty, sep="")
            ##gg["spectra"] <- as.data.frame(spectra_g)
            ##gg_comp_peaks <- display(gg["spectra"],
            ##"Parallel Coordinates Display", embed = TRUE)
            ##variables(gg_comp_peaks) <- list(X=c(colnames(spectra_g)))
            final_hbox <- gtkHBox()
            groups_hbox <- gtkHBox()
            groups_hbox$add(gg_groups)
            groups_hbox$add(gg_comp_chain)
            spec_da <- gtkDrawingArea()
            asCairoDevice(spec_da)
            final_vbox <- gtkVBox()
            final_vbox$add(groups_hbox)
            final_vbox$add(gg_comp_peaks)
            final_hbox$add(final_vbox)
            final_hbox$add(spec_da)
            groups_d <- gg["groups"]
            all_comps_d <- gg["all_comps"]

            show_comp_peaks <- function(group)
            {
              comps_g <- comps[,"group"] == group
              spectra_g <- spectra[comps_g,,drop=FALSE]
              spectra_g <- spectra_g[,which(colSums(spectra_g) > 0),drop=FALSE]
              mz_vis <- rep(TRUE, nrow(spectra))
              mz_vis[comps_g] <- FALSE
              excluded(spectra_d) <- mz_vis
              ## mz_order <- order(colMeans(spectra_g), decreasing = TRUE)
##               variables(gg_comp_peaks) <- list(X=colnames(spectra_g)[mz_order])
            }

            gSignalConnect(gg, "identify-point",
                           function(gg, plot, id, dataset) {
                             if (id == -1)
                               return()
                             if (!(as.RGtkObject(groups_d) == dataset))
                               return()
                             group <- id+1
                             glyph_color(all_comps_d) <- 1
                             group_match <- all_comps_d[,"group",drop=TRUE] %in%
                               groups_d[group,"group", drop=TRUE]
                             glyph_color(all_comps_d)[group_match] <- 9
                             shadow <- rep(TRUE, nrow(all_comps_d))
                             shadow[group_match] <- FALSE
                             shadowed(all_comps_d) <- shadow
                             ## bug in ggobi prevents this
                             show_comp_peaks(group)
                             cplotSpec(object,group=group)
                           })

            ## gSignalConnect(all_comps_d, "identify-point",
            ##                function(gg, plot, id, dataset)
            ##                {
            ##                  browser()
            ##                  if (id == -1)
            ##                    return()
            ##                  if (!(as.RGtkObject(all_comps_d) == dataset))
            ##                    return()
            ##                  comp <- all_comps_d[id+1,,drop=TRUE]
            ##                  comps_g <- which(comps[,"group"] == comp$group)
            ##                  glyph_color(spectra_d) <- 1
            ##                  m <- comps_g[comps[comps_g, "comp"] == comp$comp]
            ##                  glyph_color(spectra_d)[m] <- 9
                             
            ##                })

            ## don't forget to scale the parallel coordinate plot
           # stack_plots("Component Grouping", list(groups_hbox, gg_comp_peaks))
             stack_plots("Component Grouping", list(final_hbox))
          })

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
  if(length(inter)){
  par(mfrow=c(length(inter),1),mar=rep(0.5,4),oma=c(5,5,4,1),mgp=c(1.9,0,0))
  inter=sort(inter)
  for(i in 1:length(inter)){
    interSpec(object,inter=inter[i])
  }
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

