### Constructors for canonical pipelines

LoadPeaksPipeline <- function()
  Pipeline(Protocol("loadSample"), Protocol("genProfile"),
           Protocol("findPeaks"))
