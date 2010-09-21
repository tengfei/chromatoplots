.onLoad <- function(libname, pkgname) {
  bioc <- getOption("BioC")
  peakmethods <- bioc$xcms$findPeaks.methods
  bioc$xcms$findPeaks.methods <- c(peakmethods, "pipeline", "islands")
  options(BioC = bioc)

  ## set protocol defaults
  defaultMethod(genProfile = "intbin", loadSample = "xcms",
                loadExperiment = "xcms", fillPeaks = "extract",
                findComps = "sigma", findPeaks = "islands",
                groupComps = "angle", normalize = "scale",
                removeBaseline = "median", rtcor = "rloess",
                summarize = "common")    
}
