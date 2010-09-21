setStage("rtcor", "Correct Retention Time", "cpExperiment")

setProtocol("smooth", "Smooth",
            representation(missing = "numeric", extra = "numeric",
                           model = "character", span = "numeric",
                           family = "character"),
            retcor, "rtcor")

setProtocol("rloess", "Hierarchical Robust Loess Smoothing",
            fun = rtcor.rloess, parent = "rtcor")
