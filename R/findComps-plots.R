gg_comps <- function(object, sample, gg = ggobi(), dev = dev.new())
{
  peaks <- object@peaks[object@peaks[,"sample"] == sample,]
  comp_f <- factor(peaks[,"comp"], levels = unique(peaks[,"comp"]))
  comp <- rep(unique(peaks[,"comp"]), table(comp_f))
  df <- data.frame(comp = factor(comp))
  rownames(df) <- seq_len(nrow(df))
  gg["e_comps"] <- df
  ed <- gg["e_comps"]
  edges(ed) <- assignments_to_edges(peaks[,"comp"])
  comps_df <- as.data.frame(object@comps[object@comps[,"sample"] == sample,])
  comps_df$comp <- factor(comps_df$comp)
  gg["comps"] <- comps_df
  d <- gg["comps"]
  mz_range <- range(peaks[,"mz"])
  gSignalConnect(gg, "identify-point", function(gg, plot, id, dataset) {
    if (id == -1)
      return()
    if (!(as.RGtkObject(d) == dataset))
      return()
    if (is.function(dev))
      dev.set(dev())
    else dev.set(dev)
    comp_peaks <- peaks[peaks[,"comp"] %in% comps_df[id+1,"comp"],,drop=FALSE]
    plot_comp(comp_peaks, mz_range)
  })
  gg
}

gg_comps_box <- function(object, sample, gg = ggobi())
{
  comp_da <- gtkDrawingArea()
  asCairoDevice(comp_da)

  ## lazily resolve device ID
  dev <- 1
  gSignalConnect(comp_da, "realize", function(wid) dev <<- dev.cur())
  dev_fun <- function() dev
  
  gg <- gg_comps(object, sample, gg, dev_fun)
  gg_comps <- display(gg["comps"], embed = TRUE)
  variables(gg_comps) <- list(X = "sigma", Y = "sigma_t")
  imode(gg_comps) <- "Identify"
  
  comp_hbox <- gtkHBox()
  comp_hbox$add(gg_comps)
  comp_hbox$add(comp_da)

  comp_hbox
}

setMethod("explore", c("xcmsSet", "ProtoFindComps"),
          function(object, protocol,sample, residuals=FALSE,island=FALSE,gg = ggobi()) 
          {
            peaks <- object@peaks[object@peaks[,'sample']==sample,]
            ##  peak_proto <- findPeaksProto(processRawsProto(object)@pipeline)
            raw <- loadSample(xset_comps@filepaths[sample])
            raw <- genProfile(raw)
            peak_proto <- raw
            prof <- raw@env$profile
            alpha <- xset_comps@pipeline@.Data[[1]]@pipeline@.Data[[4]]@alpha
            quant <- quantile(prof, alpha, na.rm = TRUE)
            cutoff<- matrix(quant,nrow(prof),nrow(prof))
            peak_hbox <- gg_peaks_box(object, peak_proto, prof, sample=sample,residuals=residuals, island=island, gg)
            
            comp_hbox <- gg_comps_box(object, sample, gg)
            
            ## get the display for the peaks and active its edges
            gg_peaks <- peak_hbox$getChildren()[[1]]
            edges(gg_peaks) <- gg["e_comps"]
            
            ## link the peaks and component visualizations
            comps_d <- gg["comps"]
            peaks_d <- gg["peaks"]
            e_comps_d <- gg["e_comps"]
            comp_peak_link <- function(gg, plot, id, dataset)
              {
                if (id == -1)
                  return()
                if (!(as.RGtkObject(comps_d) == dataset))
                  return()
                glyph_color(peaks_d) <- 1
                matching_peaks <- peaks[,"comp"] %in% comps_d[id+1,"comp"]
                glyph_color(peaks_d)[matching_peaks] <- 9
                ## no way to color edges
                ##glyph_color(e_comps_d) <- 1
                ##matching_peaks <- e_comps_d[,"comp"] %in% comps_d[id+1,"comp"]
                ##glyph_color(e_comps_d)[matching_peaks] <- 9
              }
            gSignalConnect(gg, "identify-point", comp_peak_link)

            ##comp_outlier <- which.max(gg["comps"][,"sigma",drop=TRUE])
            
            stack_plots("Component Detection", list(comp_hbox, peak_hbox))
          })

gg_all_comps_disp <- function(comps, gg = ggobi())
{
  comp_order <- order(comps[,"sample"])
  assignments <- comps[comp_order,"group"]
  e <- as.character(assignments_to_edges(assignments))
  e_groups <- matrix(e, ncol=2)
  edges(gg) <- e_groups
  comps <- comps[comp_order,]
  rownames(comps) <- seq_len(nrow(comps))
  gg["all_comps"] <- comps
  gg_comp_chain <- display(gg["all_comps"], embed = TRUE)
  variables(gg_comp_chain) <- list(X = "rt", Y = "sample")
  edges(gg_comp_chain) <- gg["e_groups"]
  gg_comp_chain
}

plot_comp <- function(peaks, mz_range)
{
  df <- data.frame(x = peaks[,"mz"], xend = peaks[,"mz"], y = peaks[,"maxf"], yend = 0)
  print(ggplot(df) + geom_segment(aes(x = x, xend = xend, y = y, yend = yend)) +
        scale_y_continuous("relative intensity") + scale_x_continuous("mz", limits=mz_range))
}



