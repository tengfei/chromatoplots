# Hierarchical (by exp. design) retention time correction using robust loess
# FIXME: Should probably have 'span' and 'runs' parameters here

rtcor.rloess <- function(object) {
  corrected <- retcor_hier_robust(object@comps, object@phenoData)
  raw <- split(object@comps[,"rt"], object@comps[,"sample"])
  split(object@comps[,"rt"], object@comps[,"sample"]) <- corrected
  res <- list()
  cor <- function(sample)
    {
      lo <- loess(corrected[[sample]] - raw[[sample]] ~ raw[[sample]])
      res <<- c(res, list(residuals(lo)))
      object@rt$raw[[sample]] +
        na.flatfill(predict(lo, object@rt$raw[[sample]]))
    }
  object@rt$corrected <- lapply(seq_along(object@rt$corrected), cor)
  rtcor_res <- unsplit(res, object@comps[,"sample"])
  object@comps <- cbind(object@comps, rtcor_res = rtcor_res)
### FIXME: Probably need to correct peaks here too
  object
}


retcor_hier_robust <- function(comps, design, subset = seq_len(nrow(design)),
                               ...)
{
  if (is.logical(subset)) # ensure subset are indices
    subset <- which(subset)
  
  comp_sub <- comps[comps[,"sample"] %in% subset,]
  if (ncol(design) > 0) { # recurse along first factor
    retcor_subset <- function(ind)
      retcor_hier_robust(comps, design[,-1,drop=F], ind, ...)
    sample_ind <- split(seq_len(nrow(design))[subset], design[subset,1])
    children <- lapply(sample_ind, retcor_subset)
    #samples <- vector("list", length(subset))
    #samples[unlist(sample_ind)] <- unlist(children, F)
    rt <- unsplit(children, design[subset,1])
    weights <- lapply(rt, attr, "weights")
  } else {
    # correct replicates
    weights <- split(rep(1, nrow(comp_sub)), comp_sub[,"sample"])
    rt <- split(comp_sub[,"rt"], comp_sub[,"sample"])
  }
  
  # find medians for each group
  #rt_medians <- tapply(comp_sub[,"rt"], comp_sub[,"group"], median)
  
  # find sample with highest sum of group sizes and use it as reference
  group_sizes <- table(comp_sub[,"group"])
  groups <- split(comp_sub[,"group"], comp_sub[,"sample"])
  reference <- which.max(lapply(groups, function(group) 
    sum(group_sizes[as.character(group)])))
  
  # return models (per sample) for predicting corrected retention times
  
  lapply(seq_along(subset), function(sample) {
    # in each group, subtract median rt
    #dev <- rt[[sample]] - rt_medians[groups[[sample]]]
    group_matching <- match(groups[[sample]], groups[[reference]])
    dev <- rt[[sample]] - rt[[reference]][group_matching]
    r_lo <- robust_loess(abs(dev) ~ rt[[sample]], weights = weights[[sample]], 
      na.action = "na.exclude", b = 2, runs = 10)
    model <- loess(dev ~ rt[[sample]], weights = weights(r_lo), na.action = "na.exclude")
    # save some plots for first sample matched to other replicate
    if (FALSE && sample == 2 && ncol(design) == 0) {
      rt_sorted <- sort(rt[[sample]])
      d <- data.frame(time = rt[[sample]], dev = dev, adev = abs(dev), weight = weights(r_lo))
      p_abs <- qplot(time, adev, colour = weight, data = d, ylab = expression(abs(dev))) +
        geom_line(x = rt_sorted, y = predict(r_lo, rt_sorted), size = 2, colour = "gray50")
      ggsave(p_abs, file = "robust-loess-abs.pdf", width = 9, height = 3)
      p <- qplot(time, dev, colour = weight, data = d) +
        geom_line(x = rt_sorted, y = predict(model, rt_sorted), size = 2, colour = "gray50")
      ggsave(p, file = "robust-loess.pdf", width = 9, height = 3)
      p_zoom <- p + scale_y_continuous(limits=c(-2,2))
      ggsave(p_zoom, file = "robust-loess-zoom.pdf", width = 9, height = 3)
      print(p_zoom)
      #browser()
    }
    ord <- order(rt[[sample]]) # flat interpolate NA's on ends
    corrected <- numeric(length(ord))
    corrected[ord] <- rt[[sample]][ord] - na.flatfill(predict(model, rt[[sample]][ord]))
    old_weights <- weights[[sample]]
    new_weights <- weights(r_lo)
    # override missings (no matches with reference) with previous weights
    new_weights[is.na(new_weights)] <- old_weights[is.na(new_weights)]
    attr(corrected, "weights") <- new_weights
    corrected
  })
}

