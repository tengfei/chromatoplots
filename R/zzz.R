.onLoad <- function(libname, pkgname) {
  bioc <- getOption("BioC")
  peakmethods <- bioc$xcms$findPeaks.methods
  bioc$xcms$findPeaks.methods <- c(peakmethods, "pipeline", "islands")
  bioc$commandr$data_mode<- "memory"
  options(BioC = bioc)
  ## set protocol defaults
  getData(genProfile = NA, loadSample = NA,
                loadExperiment = NA, fillPeaks = NA,
                findComps = NA, findPeaks = NA,
                groupComps = NA, normalize = NA,
                removeBaseline = NA, rtcor = NA,
                summarize = NA)
  defaultMethod(genProfile = "intbin", loadSample = "xcms",
                loadExperiment = "xcms", fillPeaks = "extract",
                findComps = "sigma", findPeaks = "gauss",
                groupComps = "angle", normalize = "scale",
                removeBaseline = "median", rtcor = "rloess",
                summarize = "common")
}
