
normalize.scale <- function(object) {
  object@comps[,"quantity"] <- scale_log_quantities(object@comps)
  object
}

