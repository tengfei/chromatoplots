removeBaseline.median <- function(object, mzrad = 0, scanrad = 0) {
  object@env <- copyEnv(object@env)
  prof <- object@env$profile
  profMedFilt(object, mzrad, scanrad)
  object@env$profile <- prof - object@env$profile
  object
}

removeBaseline.rbe <- function(object, span = 0.15, runs = 3, b = 3,
                               progress = NULL) 
{
  rbe(object, attr(object, "proftime"), progress, span = span,
      runs = runs, b = b)$fitted
}

## rbe method
rbe <- function(object, time =NULL, progress = NULL, ...)
{
  prof <- object@env$profile
  time <- seq_len(ncol(prof))
  print(class(time))
  ## returns residuals for all points from loess fit only fit to 'scans'
  fit_rbe_to_mass <- function(mass)
    {
      intensity <- prof[mass,]
      lo <- robust_loess(intensity~time, ...)
      if (!is.null(progress))
        progress(count / nrow(prof))
      residuals(lo)
    }
  if (!is.null(progress))
    progress(0, "Fitting RBE baselines")
  id <- seq_len(nrow(prof))
  t(sapply(id, fit_rbe_to_mass))
}

robust_loess <- function(formula, data, weights, subset, na.action, runs = 3,
                         b = 3, ...)
{
  stopifnot(runs >= 0)
  stopifnot(b >= 0)
  mc <- match.call(expand.dots = FALSE)
  mc$runs <- mc$b <- mc$... <- NULL
  mc[[1]] <- as.name("model.frame")
  mf <- eval(mc, parent.frame())
  weights <- model.weights(mf)
  if (is.null(weights))
    weights <- rep(1, nrow(mf))
  
  ## for calling 'loess'
  lo_call <- as.call(c(loess, list(...)))
  lo_call$trace.hat <- "approximate"
  if ("span" %in% names(lo_call))
    span <- lo_call$span
  else span <- 0.1 # 500/length(weights)
  
  for (i in seq_len(runs+1)) {
    mf$"(weights)" <- weights[!is.na(weights)]
    lo_call$formula <- mf
    ## there is likely a problem here using a fixed 'span' - if a large number
    ## of points are excluded (zero weight) then the neighborhoods may become
    ## too small. We need to scale the 'span' according to the zero count.
    lo_call$span <- span * length(weights) / sum(weights > 0, na.rm = TRUE)
    model <- eval(lo_call, parent.frame())
    res <- residuals(model)
    s_mav <- ifelse(sum(res < b*mav(res), na.rm=TRUE) != 0,
                    mav(res[res < b*mav(res)]), 0)
    weights <- ifelse(res < 0, 1, pmax(1 - (res / s_mav / b) ^ 2, 0))
    weights[is.nan(weights)]<-0
    if (sum(weights, na.rm = TRUE) == 0)
      break
  }
  model
}

