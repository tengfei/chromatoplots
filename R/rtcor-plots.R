plot_tic_rtcor <- function(raws, time, label = "corrected") {
  init <- data.frame(time = time[[1]], tic = raws[[1]]@tic)
  p <- ggplot(init, aes(x = time, y = tic)) + geom_line() +
    scale_x_continuous(limits=c(360,500)) +
    scale_y_continuous(limits=c(0,5e6))
  for (i in tail(seq_along(raws), -1))
  {
    p <- p + geom_line(x = time[[i]], y = raws[[i]]@tic)
  }
  p$title <- paste("TIC Chromatograms (", label, ")", sep="")
  p
}

setMethod("explore", c("xcmsSet", protocolClass("rtcor")),
          function(object, protocol, raw, xscale=NA,geom=NA,sample=NA,log=F,gg = ggobi())
          {
            gg_comp_chain <- gg_all_comps_disp(object@comps, gg)

##              s <- object@comps[,"sample"] == 5
##              retcor_d <- data.frame(raw = raw@comps[s,"rt"],
##                                     cor = object@comps[s,"rt"])
##              p_retcor_fit <- qplot(raw, cor-raw, data = retcor_d)+
##                scale_x_continuous(limits=c(360,1000)) + geom_smooth()
            
##             rawpipeline <- pipeline(pipeline(processRawsProto(object)),
##                                     outtype = "xcmsRaw")
##             raws <- lapply(object@filepaths, xcmsRaw, pipeline = rawpipeline)
##             p_tic_cor <- plot_tic_rtcor(raws, object@rt$corrected)
##             p_tic_cor$title <- NULL
            
##             stack_plots("Retention Time Correction",
##                         list(p_tic_cor, p_retcor_fit, gg_comp_chain))
            p_retcor_fit <- cplotRtFit(object,raw,xscale,sample)
            p_rt <- cplotRT(object,xscale=xscale,sample=sample,geom=geom,log=log)
            stack_plots("Retention Time Correction",
                        list(p_rt,p_retcor_fit,gg_comp_chain))
          })
