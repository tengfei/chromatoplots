## There are two peak detection/fitting algorithms: one uses an NLS gaussian
## fit, the other is an approximation through parabolas.

findPeaks.gauss <- function(object, alpha = 0.99, egh = TRUE) {
  prof <- object@env$profile
  quant <- quantile(prof, alpha, na.rm = TRUE)
  fits <- apply(prof, 1, fit_profile, object@scantime, quant, egh)
  fits_flat <- unlist(fits, F)
  mz <- rep(profMz(object), sapply(fits, length))
  pks <- fits_to_peaks(fits_flat, mz)
  if(nrow(pks)==0) cat('No peaks are found, please adjust your parameters\n')
  pks
}

descendMin <- function(y, istart = which.max(y)) {
  if (!is.double(y)) y <- as.double(y)
  unlist(.C("DescendMin",
            y,
            length(y),
            as.integer(istart-1),
            ilower = integer(1),
            iupper = integer(1),
            DUP = FALSE, PACKAGE = "xcms")[4:5]) + 1
}



SSfixed <- selfStart(~ h*exp(-(x-mu)^2/(2*2.3^2)), function(mCall, data, LHS) {
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  len <- dim(xy)[1]
  maxpos <- which.max(xy[,2])
  mu <- xy[maxpos,1]
  h <- xy[maxpos,2]
  value <- c(mu, h)
  names(value) <- mCall[c("mu", "h")]
  value

}, c("mu", "h"))

## Exponential-Gaussian hybrid

egh <- function(x, mu, h, sigma, t) {
  sapply(x, function(x)
         if ((2*sigma^2 + t*(x - mu)) > 0)
         h*exp(-(x-mu)^2/(2*sigma^2 + t*(x-mu)))
         else 0
         )
}

## ~ egh(mu, h, sigma, t)
SShybrid <- selfStart(egh, function(mCall, data, LHS) {
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  
  integ <- function(a, b)
    sum((xy[(a+1):b,2]+xy[a:(b-1),2])*(xy[(a+1):b,1]-xy[a:(b-1),1]))/2
  
  maxpos <- which.max(xy[,2])
  maxval <- xy[maxpos,2]
  mu <- xy[maxpos,1]
  h <- xy[maxpos,2]
  sigma <- integ(1, nrow(xy))/(h*sqrt(2*pi))
  t <- 0
  value <- c(mu, h, sigma, t)
  names(value) <- mCall[c("mu", "h", "sigma", "t")]
  value

}, c("mu", "h", "sigma", "t"))

fit_profile <- 
  function(profile, times, 
           quant = quantile(profile, .99, na.rm = TRUE), egh = TRUE, min_points = 3, 
           noise_factor = 3, min_fit_radius = min_points*2) 
{
  ## fit a specific region of the profile, fill results into 'env'
  fit_region <- function(region, plot = FALSE) {
    if (diff(region) <= 1) # ignore empties
      return(NULL)

    noise.sd <- sd(profile[profile < quant], na.rm=T)
    noise <- noise_factor*noise.sd

    ## try a different way to define noise ratio.
    ##range <- (diff(region) + 1) * 2
    ## try a smaller region
    range <- max(round((diff(region) + 1)/2), min_fit_radius)
    maxpos <- region[1] + which.max(profile[region[1]:region[2]]) - 1
    ## intervals are trimmed, based on crossing above the quantile cutoff again and
    ## whether an increase happens in the wrong direction while within regions
    left <- max(1, region[1] - range)
    left_crossed <- which(profile[left:(region[1]-1)] > quant)
    if (length(left_crossed))
      left <- left + left_crossed[length(left_crossed)]
    left_down <- diff(profile[region[1]:maxpos]) < 0
    local_max <- which(diff(c(FALSE, left_down)) > 0) - 1 + region[1]
    local_min <- as.numeric(sapply(local_max, function(m) 
                                   descendMin(profile[m:maxpos], 1)[2] + m - 1))

    ## collect convoluted data
    ##  if((profile[local_max] - profile[local_min])/profile[local_max] > noise){
    ##   sink('~/Desktop/local2.txt',append=T)
    ##   cat('local_max:',profile[local_max],'local_min:',profile[local_min],'diff:',
    ##        profile[local_max] - profile[local_min],'noise:',noise,
    ##        'result:',profile[local_max] - profile[local_min] > noise,'\n')
    ##   sink()
    ## }
    ## normalize local diff
    ## local_min <- local_min[(profile[local_max] - profile[local_min])/profile[local_max] > noise]
    local_min <- local_min[profile[local_max] - profile[local_min] > noise]
    
    if (length(local_min)) {
      local_min <- local_min[1]
      fit_region(c(region[1], local_min))#, profile, times, quant, env, plot)
      fit_region(c(local_min, region[2]))#, profile, times, quant, env, plot)
      return()
    }
    right <- min(length(profile), region[2] + range)
    right_crossed <- which(profile[(region[2]+1):right] > quant)
    if (length(right_crossed))
      right <- region[2] + right_crossed[1] - 1
    right_up <- diff(profile[maxpos:region[2]]) > 0
    local_max <- which(diff(c(right_up, FALSE)) < 0) + maxpos
    local_min <- as.numeric(sapply(local_max, function(m) 
                                   descendMin(profile[maxpos:m], m - maxpos + 1)[1] + maxpos - 1))
    local_min <- local_min[profile[local_max] - profile[local_min] > noise]
    if (length(local_min)) {
      local_min <- local_min[1]
      fit_region(c(region[1], local_min))#, profile, times, quant, env, plot)
      fit_region(c(local_min, region[2]))#, profile, times, quant, env, plot)
      return()
    }
    ##  right <- maxpos + right_up[1] - 1
    ##xcms <- try(nls(y ~ SSfixed(x, mu, h), data.frame(x = times[left:right], 
    ##  y = profile[left:right])), T)
    hybrid <- character(0)
    if (egh)
      hybrid <- try(nls(y ~ SShybrid(x, mu, h, sigma, t), 
                        data.frame(x = times[left:right], y = profile[left:right])), T)
    if (is.character(hybrid))
      gauss <- try(nls(y ~ SSgauss(x, mu, sigma, h), data.frame(x = times[left:right], 
                                                                y = profile[left:right])), T)
    if (plot) {
      plot(left:right, profile[left:right], xlab = "scan", ylab = "residuals")
      if (is.character(hybrid) && !is.character(gauss)) {
        points(left:right, fitted(gauss), type = "l", col = "red", lwd = 2)
        print(coefficients(gauss))
      }
      ##if (!is.character(xcms))
      ##  points(left:right, fitted(xcms), type = "l", col = "blue", lwd = 2)
      if (!is.character(hybrid)) {
        points(left:right, fitted(hybrid), type = "l", col = "green", lwd = 2)
        print(coefficients(hybrid))
      }
    }
    name <- paste(region[1], ":", region[2], sep="")
    if (is.character(hybrid))
      model <- gauss
    else model <- hybrid
    assign(name, model, env)
  }
  filt <- profile > quant & !is.na(profile)
  diffs <- diff(c(FALSE, filt, FALSE))
  regions <- cbind(which(diffs == 1), which(diffs == -1)-1)
  ## here we ignore regions with less than min_points points
  regions <- regions[regions[,2] - regions[,1] >= min_points - 1,,drop=F]
  env <- new.env()
  apply(regions, 1, fit_region)
  as.list(env)
}

fits_to_peaks <- function(fits, masses = NULL)
{
  good_fits <- !sapply(fits, is.character)
  fits <- fits[good_fits]
  names(fits) <- NULL 
  coeff_names <- c("mu", "h", "sigma", "t")
  integ <- function(xy, a, b)
    sum((xy[(a+1):b,2]+xy[a:(b-1),2])*(xy[(a+1):b,1]-xy[a:(b-1),1]))/2
  peaks <- t(sapply(fits, function(model) {
    e <- model$m$getEnv()
    span <- diff(range(e$x))
    coeffs <- coefficients(model)[coeff_names]
    coeffs[is.na(coeffs)] <- 0 # if couldn't fit egh, assume 't' is zero
    coeffs[["sigma"]] <- abs(coeffs[["sigma"]])
    names(coeffs)<-coeff_names
    c(rt = coeffs[["mu"]],
      rtmin = e$x[1],
      rtmax = tail(e$x, 1),
      into = integ(cbind(e$x, e$y), 1, length(e$x)),
      intf = sqrt(2*pi)*coeffs[["h"]]*coeffs[["sigma"]],
      maxo = max(e$y),
      maxf = coeffs[["h"]],
      tau = coeffs[["t"]],
      sigma = coeffs[["sigma"]],
      error = sqrt(deviance(model)/df.residual(model)), 
      span = span / 5, # is this necessary?
      delta_i = coeffs[["h"]] - max(e$y),
      delta_t = coeffs[["mu"]] - e$x[which.max(e$y)],
      peak_sd = sd(e$y) ,
      "sigma/h" = coeffs[["sigma"]]/coeffs[["h"]]
      )
  }))

  if (!is.null(masses)) {
    masses <- masses[good_fits]
    peaks <- cbind(mz = masses, mzmin = masses, mzmax = masses, peaks)
  }

  ## try to get rid of the "flat" peaks
  ##  flat <- peaks[,"span"] < peaks[,"sigma.sigma"] | peaks[,"ret.mu"] < object@scantime[peaks[,"retmin"]] |
  ##    peaks[,"ret.mu"] > object@scantime[peaks[,"retmax"]]
  ##  peaks[!flat,]
  ## get rid of the 'bad' fit peaks and 'flat' peaks
  ## bad <- peaks[,'rt']<peaks[,'rtmin']|peaks[,'rt']>peaks[,'rtmax']|abs(peaks[,'sigma'])<=sigma.cutoff
  ## cat(length(bad),'\n')
  ## cat(nrow(peaks),'\n')
  ## peaks <- peaks[!bad,]
  new("cpPeaks", peaks)
}


#############################################################
## Peak detection using parabola approximation
#############################################################

findPeaks.parabola <- function(object) {
  raw <- object
  ## the white noise is assumed to have a poisson distribution
  ## if the noise is significant, we can approximate the poisson with the normal
  ## we then estimate sigma from the MAD (sigma = 1.4826*MAD)
  ## (could become an estimateNoise protocol)
  prof <- raw@env$profile
  mad_prof <- xcms::medianFilter(abs(prof), mrad = 0, nrad = 100)
  ## Filter out maxima that are piggybacking on peaks. With small peaks
  ## or the flat tails of large peaks, noise can cause the detection of
  ## spurious maxima, as the values are similar and above the noise
  ## cutoff. In a sense, the peak is adding an additional baseline. One
  ## could attempt to fit that baseline. AMDIS locally fits a line,
  ## which is just an approximation of fitting the actual peak. This
  ## approach is a bit complicated though, and it's not clear if it is
  ## needed since we are already subtracting a global baseline. AMDIS
  ## also checks whether a maximum rises above an adjacent minimum by
  ## some number of noise units. All of these seem to be variants of
  ## smoothing. What about just running a loess?  For now, trying a
  ## median filter of radius 2. Spurious maxima are avoided, but maxima
  ## tend to be misestimated. It's hard to say where the maxima is, but
  ## when a minimum is chosen, the fit suffers. For now, when the
  ## filtered signal is flat, we pick the maximum.
  filt_prof <- xcms::medianFilter(prof, 0, 2)
  ##filt_prof <- prof

  ## now we filter out the noise

  ## hack: convert profile matrix into a single series, compress to Rle
  ## the Rle has significant overhead -- how much memory must we save?
  ## right now, the compression is about 5X, while a sparse-vector might
  ## give 6X.

  rleProf <- function(x) Rle(as.integer(t(cbind(x, NA))))

  series <- rleProf(raw_prof@env$profile)
  cor_series <- rleProf(prof)
  mad_series <- rleProf(mad_prof)
  filt_series <- rleProf(filt_prof)

  ## We take the AMDIS approach of computing a
  ## global (median) noise factor, ignoring regions with imputed values
  ## (zeros). Then just multiple the sqrt(raw signal) by that factor.

  ## tricky: make sure we don't violate the profile row boundaries
  imputed <- runValue(series) == min(series, na.rm=TRUE)
  imputed[is.na(imputed)] <- FALSE
  breaks_prof <- which(is.na(runValue(series)))
  ## start_prof <- c(1L, head(start(series)[breaks_prof] + 1L, -1))
  ## start_prof <- rep(start_prof, diff(c(0L, breaks_prof)))
  ## end_prof <- end(series)[breaks_prof] - 1L
  ## end_prof <- rep(end_prof, diff(c(0L, breaks_prof)))
  ## prof_r <- IRanges(start_prof, end_prof)
  breaks_r <- PartitioningByEnd(start(series)[breaks_prof])
  prof_r <- IRanges(start(breaks_r), width = width(breaks_r) - 1)
  prof_r <- prof_r[rep(seq_len(length(prof_r)), diff(c(0L, breaks_prof)))]
  imputed_ir <- IRanges(start(series)[imputed] - 50, end(series)[imputed] + 50)
  imputed_ir <- pintersect(imputed_ir, prof_r[imputed])
  primary_ir <- gaps(imputed_ir, 1L, length(series))
### FIXME: should not use the single point intensity here, instead
### take value from median fit in baseline subtraction.
  noise <- mad_series[primary_ir] / sqrt(series[primary_ir])
  nf <- median(noise, na.rm = TRUE)

  ## a reasonable peak cut-off is 3*sd, multiply this by sqrt(intensity)
  cutoff_series <- nf*3*1.4826*sqrt(series)

  ## peak fitting (per ion and then TIC)

  ## 1) find local maxima that pass a noise filter
  ##max_series <- which(diff(diff(filt_series) > 0) == -1) + 1
  diff_series <- diff(filt_series)
  dd_series <- diff(diff_series > 0)
  max_series <- which(dd_series == -1) + 1
  max_r <- findRange(max_series, filt_series)
  max_views <- Views(series, max_r[which(diff_series[end(max_r)] < 0)])
  max_series <- viewWhichMaxs(max_views)
  ## maximum must be above noise cutoff
  max_series <- max_series[which(cor_series[max_series] >
                                 cutoff_series[max_series])]
  
  ## filter out maxima without flanking scans above threshold could
  ## reduce the threshold (maybe 2*sd), since it is less likely for
  ## three points to be above the threshold (from the geometric), but we
  ## are not just trying to detect peaks. we are trying to characterize
  ## them, which is tough if there are not enough non-noisy points.
  max_series <- max_series[which(cor_series[max_series-1] >
                                 cutoff_series[max_series-1])]
  max_series <- max_series[which(cor_series[max_series+1] >
                                 cutoff_series[max_series+1])]
  max_runs <- findRun(max_series, series)

  ## 2) fit each maximum with a gaussian

  ## could just fit a parabola, which matches the top-half of the
  ## gaussian fairly well. then estimate sigma from sigma =
  ## FWHM/(2*sqrt(2*log(2))). We could take the 5 points centered on the
  ## maximum. For small peaks, the fit will likely be affected by noise,
  ## but the measurements lack accuracy at low levels anyway.
  
  ## do we want to down-weight points (the ends) that fall under the
  ## threshold? Maybe by 0.5?

  fit_max <- start(series)[max_runs]
  fit_r <- IRanges(fit_max - 2L, fit_max + 2L)
  fit_r <- pintersect(fit_r, prof_r[max_runs])
  fit_fw <- width(fit_r) == 5L ## some are on the end of time course
  max_runs <- max_runs[fit_fw]
  fit_r <- fit_r[fit_fw] 
  ## avoid overlap between fit regions, splitting at minimum
  ## fit_crossed <- which(head(end(fit_r), -1) >= tail(start(fit_r), -1))
  ## imin <- viewWhichMins(Views(cor_series, fit_max[fit_crossed],
  ##                             fit_max[fit_crossed+1]))
  ## end(fit_r)[fit_crossed] <- imin
  fit_max_off <- integer(length(fit_r))
  ## fit_max_off[fit_crossed+1] <- start(fit_r)[fit_crossed+1] - imin
  ## start(fit_r)[fit_crossed+1] <- imin
  fit_f <- rep(seq_len(length(fit_r)), width(fit_r))
  fit_int <- as.integer(cor_series[fit_r])
  fit_y <- split(fit_int, fit_f)
  fit_r_scan <- IRanges(start(fit_r) - start(prof_r)[max_runs] + 1,
                        w = width(fit_r))
  fit_t <- raw@scantime[as.integer(fit_r_scan)]
  ## to avoid numerical overflow
  fit_t_off <- raw@scantime[start(fit_r_scan)]
  fit_t <- fit_t - rep(fit_t_off, width(fit_r))
  fit_x <- cbind(a = 1L, b = fit_t, c = fit_t^2)
  fit_d <- split.data.frame(fit_x, fit_f)
  fits <- do.call("rbind", lapply(seq_len(length(fit_r)), function(r) {
    coefficients(lm.fit(fit_d[[r]], fit_y[[r]]))
  }))

  ## determine mu, sigma, and h
  a <- fits[,"a"]
  b <- fits[,"b"]
  c <- fits[,"c"]
  mu <- -b/(2*c)
  h <- a + b*mu + c*(mu^2)
  mu <- mu + fit_t_off # correct 'mu' after using it with 'a', 'b', or 'c'
  sigma <- suppressWarnings((-sqrt(b^2-4*(a-h/2)*c)/c)/2.35703) # simplify?
  fit_t_i <- head(cumsum(c(1, width(fit_r))), -1)
  fit_t_max <- fit_t_i + fit_max_off
  span <- width(fit_r)
  fit_t_abs <- fit_t + rep(fit_t_off, span)

  ## integrate points under each peak
  xy <- cbind(fit_t, fit_int)
  m <- 1
  n <- nrow(xy)
  t_int <- (xy[(m+1):n,2]+xy[m:(n-1),2])*(xy[(m+1):n,1]-xy[m:(n-1),1])
  if (length(fit_t_i > 1)) # get rid of intervening points
    t_int <- t_int[-(tail(fit_t_i, -1) - 1)]
  into <- colSums(matrix(t_int, span-1))/2
  ##into <- viewSums(Views(Rle(t_int), successiveIRanges(span-1)))/2

  ## calculate residuals and error
  resid <- fit_int -
    (rep(a, span) + rep(b, span)*fit_x[,2] + rep(c, span) * fit_x[,3])
  ##error<-sqrt(viewSums(Views(Rle(resid^2), successiveIRanges(span))) /
  ##            (span - 3))
  error <- sqrt(colSums(matrix(resid^2, span))/(span - 3))

  minmz <- raw@mzrange[1]
  mz <- findInterval(start(fit_r), which(is.na(series))) + minmz
  
  peaks <-
    cbind(mz = mz, mzmin = mz, mzmax = mz,
          rt = mu, rtmin = mu - 3*sigma, rtmax = mu + 3*sigma,
          into = into, intf = sqrt(2*pi)*h*sigma,
          maxo = fit_int[fit_t_max], maxf = h,
          tau = 0, sigma = sigma, error = error, span = span,
          rto = fit_t_abs[fit_t_max])

  ## scanmin <- pmax(1, findInterval(peaks[,"rtmin"], raw@scantime))
  ## scanmax <- pmin(findInterval(peaks[,"rtmax"], raw@scantime)+1, ncol(prof))
  ## scan_off <- cumsum(c(0,rep(ncol(prof),nrow(prof)-1)))[peaks[,"mz"]-minmz+1]
  ## peaks <- cbind(peaks, scanmin = scanmin + scan_off,
  ##                scanmax = scanmax + scan_off)

  ## peak filtering: there are 3 types of peaks: the good fits, the weak
  ## fits and the bad fits. The bad fits are those that we cannot trust
  ## at all, and so should be discarded. The weak fits are often due to
  ## convolution and may be resolved using the good fits in the
  ## deconvolution stage.

  ## filter out obviously bad fits
  ## "upside down" or 'mu' estimated outside of fitted domain
  ## peaks <- peaks[c < 0 & mu <= fit_t[fit_t_i + width(fit_r) - 1] + fit_t_off &
  ##                mu >= fit_t[fit_t_i] + fit_t_off,]
  new("cpPeaks", peaks)
}
