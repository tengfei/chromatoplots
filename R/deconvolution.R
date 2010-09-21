find_crossed_peaks <- function(peaks, tolerance = 3) {
  find_crossed_at_mass <- function(m) {
    peak_ind <- which(peaks[,"mz"] == m)
    p <- peaks[peak_ind,,drop=F]
    spread <- tolerance * abs(p[,"sigma"])
    left <- p[,"rt"] - spread
    right <- p[,"rt"] + spread
    crossed <- outer(left, right, "-") < 0 & outer(right, left, "-") > 0
    crossed[!lower.tri(crossed)] <- FALSE
    crossed_ind <- lapply(as.data.frame(crossed), function(col) peak_ind[col])
    cbind(a = rep(peak_ind, sapply(crossed_ind, length)), 
          b = unlist(crossed_ind))
  }
  crossed_list <- lapply(unique(peaks[,"mz"]), find_crossed_at_mass)
  first_col <- rep(rep(c(TRUE, FALSE), length(crossed_list)), 
    rep(sapply(crossed_list, nrow), each=2))
  crossed_vec <- unlist(crossed_list)
  cbind(a = crossed_vec[first_col], b = crossed_vec[!first_col])
}

find_crossed_peaks2 <- function(peaks) {
  start <- peaks[,"rtmin"]
  end <- peaks[,"rtmax"]
  mz <- peaks[,"mz"]
  start_ord <- order(start)
  peak_ord <- start_ord[order(mz[start_ord])]
  crossed <- which(tail(start[peak_ord], -1) <= head(end[peak_ord], -1))
  mat <- cbind(a = peak_ord[crossed], b = peak_ord[crossed+1])
  mat[mz[mat[,1]] == mz[mat[,2]],]
}
