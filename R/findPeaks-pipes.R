### Peak finding

setStage("findPeaks", "Find Peaks", "cpSample", "cpPeaks")

setProtocol("gauss", "Gaussian Fitting",
            representation(alpha = "numeric", egh = "logical"),
            findPeaks.gauss, "findPeaks")

setProtocol("parabola", "Parabola Fitting",
            representation(alpha = "numeric", egh = "logical"),
            findPeaks.parabola, "findPeaks")

setProtocol("matchedFilter", "Matched Filter",
            representation(fwhm = "numeric", sigma = "numeric", max = "numeric",
                           snthresh = "numeric", step = "numeric",
                           steps = "numeric", mzdiff = "numeric",
                           index = "logical"),
            getMethodLocal("findPeaks.matchedFilter", "xcmsRaw"), "findPeaks")

setProtocol("centWave", "Centroid Wavelet",
            representation(scanrange="numeric", minEntries="numeric",
                           dev="numeric", snthresh="numeric",
                           noiserange="numeric", minPeakWidth="numeric",
                           scales="numeric", maxGaussOverlap = "numeric",
                           minPtsAboveBaseLine="numeric",
                           scRangeTol="numeric", maxDescOutlier="numeric",
                           mzdiff="numeric", rtdiff="numeric",
                           integrate="numeric", fitgauss = "logical"),
            getMethodLocal("findPeaks.centWave", "xcmsRaw"), "findPeaks")

## setProtocol("MS1", "Load MS1 Precursor Peaks", representation(),
##            getMethodLocal("findPeaks.MS1", "xcmsRaw"), "findPeaks")

## lets us inject a pipeline into xcmsSet() (should not be exposed too much)
findPeaks.pipeline <- function(data, pipeline = Protocol("findPeaks"))
{
  perform(pipeline, data)
}
setProtocol("pipeline", "Via Pipeline", representation(pipeline = "Pipeline"),
            findPeaks.pipeline, "findPeaks")
setMethod("findPeaks.pipeline", "xcmsRaw", function(object, ...) {
  findPeaks.pipeline(as(object, "cpSample"), ...)
})
