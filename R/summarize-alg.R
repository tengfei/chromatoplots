
summarize.common <- function(object) {
  common <- sum_common_peaks(object@peaks, object@comps)
  comps <- object@comps[,colnames(object@comps) != "quantity"]
  object@comps <- cbind(comps, quantity = common)
  object
}

find_spectra <- function(peaks)
{
  comp_peaks <- split(seq_len(nrow(peaks)), interaction(peaks[,"comp"], peaks[,"sample"]))
  existing_ints <- sapply(comp_peaks, length) > 0
  comp_peaks <- comp_peaks[existing_ints]
  mass <- peaks[,"mz"]
  height <- peaks[,"maxf"]
  sigma <- peaks[,"sigma"]
  max_mass <- max(mass)
  spectrum <- function(p) {
    y <- numeric(max_mass)
    y[mass[p]] <- 2*height[p]*sqrt(pi/2)*sigma[p]
    y
  }
  sapply(comp_peaks, spectrum)
}

sum_common_peaks <- function(peaks, comps)
{
  # ensure component order is compatible with find_spectra()
  comps <- cbind(comps, id = seq_len(nrow(comps)))
  comps <- comps[order(comps[,"comp"]),]
  comps <- comps[order(comps[,"sample"]),]
  spectra <- t(find_spectra(peaks))
  sums <- unsplit(by(spectra, comps[,"group"], function(group) {
    masses <- apply(group, 2, function(column) all(column != 0))
    apply(group[,masses,drop=FALSE], 1, sum)
  }), comps[,"group"])
  result <- numeric(length(sums))
  result[comps[,"id"]] <- sums
  result
}

scale_log_quantities <- function(comps)
{
  quantity <- suppressWarnings(log(comps[,"quantity"]))
  quantity[is.infinite(quantity)] <- NA
  unsplit(tapply(quantity, comps[,"sample"], function(sample)
  {
    sample - mean(sample, na.rm = TRUE)
  }), comps[,"sample"]) + mean(quantity, na.rm = TRUE)
}
