# Stuff that's not quite there
if (FALSE) {
  # Correct retention time by matching landmarks in TIC
  
  # find TIC peaks detected using 'alpha' as the quantile cutoff,
  # match and fit rbe to differences
  .retcor.landmarks <- function(object, alpha = 0.99) {
    # find peaks in TIC for each sample
    landmarks <- lapply(object@raw, function(raw) {
      tic <- raw@tic
      fits <- fit_series(tic, quantile(tic, .95), raw@scantime)
      fits_to_peaks(fits, raw@scantime)
    })
    # landmarks ==> comps
    # comps ==> peaks
    matching <- match_components(comps, peaks, object@phenoData, FALSE)
    # correct rt by matching
  }
  
  setGeneric("retcor.landmarks", function(object, ...) standardGeneric("retcor.landmarks"))
  setMethod("retcor.landmarks", "xcmsSet", .retcor.landmarks)
  
  setProtocolClass("ProtoRetcorLandmarks",
    representation(alpha = "numeric"),
    c(formals(.retcor.landmarks), dispname = "TIC Lardmarks"),
    "ProtoRetcor")
}
