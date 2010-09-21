## shoule change tolerance, set it as argument

findComps.sigma <- function(object,tolerance=1) {
  ## we assume that the peaks are sorted by sample
  components <- by(object@peaks, object@peaks[,"sample"],
                   find_components,tolerance=tolerance)
  comp_list <- lapply(components, "[[", "comps")
  comps <- do.call("rbind", comp_list)
  samples <- rep(seq_along(components), sapply(comp_list, nrow))
  comps <- cbind(comps, sample = samples)
  object@comps <- comps
  peak_comps <- unlist(lapply(components, "[[", "assignments"))
  object@peaks <- cbind(object@peaks, comp = peak_comps)
  object
}

find_components <- function(peaks, tolerance = 1)
{
  ##crossed_ind <- find_crossed_peaks2(peaks)
  #crossed <- logical(nrow(peaks))
  ##crossed[crossed_ind] <- TRUE
  peaks <- cbind(peaks, comp = NA)
  ##peaks[!crossed,"comp"] <- componentize_peaks(peaks[!crossed,], tolerance)
  peaks[,"comp"] <- componentize_peaks(peaks, tolerance)
  comps <- assignments_to_components2(peaks[,"comp"], peaks)
  assignments <- components_to_assignments2(comps, peaks)
  peaks[,"comp"] <- assignments
  comps <- assignments_to_components2(peaks[,"comp"], peaks)
  resolve_orphans(comps, peaks)
}

componentize_peaks <- function(peaks, tolerance = 1)
{

  # add indices to peaks in order to keep track of them
  peaks_ind <- cbind(peaks, ind = seq(length = nrow(peaks)))
##  uncomped <- peaks_ind[order(peaks_ind[,"maxf"], decreasing = T),]
  sigma_ord <- order(peaks_ind[,"sigma"])
  ##q <- quantile(seq_along(sigma_ord), c(0.25, 0.5), type = 1)
  ##q <- quantile(seq_along(sigma_ord), c(0.25, 0.75), type = 1)
  ##q[is.na(q)] <- 0
  ## order by sigma away from first quartile, until median, then in order
  ## ord <- c(rbind(rev(head(sigma_ord, q[1])),
  ##                head(tail(sigma_ord, -q[1]), q[1])),
  ##          tail(sigma_ord, -(q[2] - (q[2] - 2*q[1]))))
  ## order by sigma from first to third quartile, then first to start, then rest
  ## ord <- c(sigma_ord[seq(q[1], q[2])], rev(head(sigma_ord, q[1]-1)),
  ##          tail(sigma_ord, -q[2]-1))
  ord <- sigma_ord ## just order by sigma
  uncomped <- peaks_ind[ord,]
  ## consider crossed/convoluted peaks last
  ## crossed_ind <- find_crossed_peaks2(peaks)
  ## crossed <- logical(nrow(peaks))
  ## crossed[crossed_ind] <- TRUE
  ## peaks <- rbind(peaks[!crossed,], peaks[crossed,])
  assignments <- rep(NA, nrow(peaks))
  comp_ind <- 1
  while(nrow(uncomped)) {
    ##sigma <- abs(uncomped[1, "sigma"])*tolerance
    rt <- uncomped[,"rt"]
    ##overlap <- (rt < rt[1] + sigma) & (rt > rt[1] - sigma)
    sigma <- uncomped[,"sigma"]*tolerance
    overlap <- (rt + sigma > rt[1]) & (rt - sigma < rt[1])
    overlap[1] <- TRUE
    comped <- uncomped[overlap,,drop=F]
    # make sure no duplicates across mass
    dup_mass <- comped[,"mz"] %in% unique(comped[duplicated(comped[,"mz"]),"mz"])
    if (any(dup_mass)) {
      dup_ind <- by(comped[dup_mass,], comped[dup_mass,"mz"], function(m)
        m[-which.min(abs(m[,"rt"] - rt[1])),"ind"]) 
      dups <-logical(max(uncomped[,"ind"]))
      dups[unlist(dup_ind)] <- TRUE
      overlap <- overlap & !dups[uncomped[,"ind"]]
    }
    overlap_ind <- uncomped[overlap, "ind"][order(uncomped[overlap,"mz"])]
    assignments[overlap_ind] <- comp_ind
    if (length(overlap_ind))
      comp_ind <- comp_ind + 1
    uncomped <- uncomped[!overlap,,drop=F]
  }
  assignments
}

### Far slower than initial implementation
componentize_peaks2 <- function(peaks, tolerance = 1) {
  mu <- peaks[,"rt"]
  sigma <- peaks[,"sigma"]
  ir <- IRanges((mu - tolerance*sigma)*1000, (mu + tolerance*sigma)*1000)
  ol <- as.matrix(overlap(ir, as.integer(mu*1000)))
  ## if comp contains peaks from same mz, only keep closest peak
  mu_delta <- abs(mu[ol[,1]] - mu[ol[,2]])
  delta_ord <- order(mu_delta)
  ord <- delta_ord[order(peaks[ol[,2],"maxf"][delta_ord], decreasing=TRUE)]
  ol <- ol[!duplicated(data.frame(peaks[ol[,1],"mz"], ol[,2])[ord,]),]
  ## make sure each peak is only assigned to one component
  ol <- ol[!duplicated(ol[,1]),]
  comps <- integer(nrow(peaks))
  comps[ol[,1]] <- as.integer(factor(ol[,2]))
  comps
}

nearest_neighbor <- function(x, vec) {
  nearest <- integer(length(x))
  vec_ord <- NULL
  if (is.unsorted(vec)) {
    vec_ord <- order(vec)
    vec <- vec[vec_ord]
  }
  int <- findInterval(x, vec)
  left <- int == 0
  nearest[left] <- 1
  int[int == length(vec)] <- int[int == length(vec)] - 1
  x <- x[!left]
  int <- int[!left]
  nearest[!left] <- ifelse(x - vec[int] < x - vec[int+1], int, int + 1)
  if (!is.null(vec_ord))
    nearest <- vec_ord[nearest]
  nearest
}

## 50X faster than initial
components_to_assignments2 <- function(components, peaks, tolerance = 1) {
  mu <- components[,"rt"]
  mu_peak <- peaks[,"rt"]
  sigma <- peaks[,"sigma"]
##  sigma <- components[,"sigma"]
##  ir <- IRanges((mu - tolerance*sigma)*1000, (mu + tolerance*sigma)*1000)
##  ol <- as.matrix(overlap(ir, as.integer(mu_peak*1000)))
  ir <- IRanges((mu_peak - tolerance*sigma)*1000,
                (mu_peak + tolerance*sigma)*1000)
  ol <- as.matrix(overlap(ir, as.integer(mu*1000)))
  ## if comp contains peaks from same mz, only keep closest peak
  mu_delta <- abs(mu_peak[ol[,2]] - mu[ol[,1]])
  ol <- ol[order(mu_delta),]
  ol <- ol[!duplicated(data.frame(peaks[ol[,2],"mz"], ol[,1])),]
  ## make sure each peak is only assigned to one component
  ol <- ol[!duplicated(ol[,2]),]
  comps <- peaks[,"comp"] + max(ol[,1])
  comps[ol[,2]] <- components[ol[,1], "comp"]
  as.integer(factor(comps))
}

components_to_assignments <- function(components, peaks, tolerance = 1) {
  assignments <- rep(NA, nrow(peaks))
  by(cbind(peaks, ind = seq_len(nrow(peaks))), peaks[,"mz"], function(m) {
    d <- abs(outer(m[,"rt"], components[,"rt"], "-"))
    sigma <- matrix(components[,"sigma"],nrow=nrow(d),ncol=ncol(d),byrow=T)
    d[d > abs(sigma)*tolerance] <- NA
    d_col <- col(d)
    d_row <- row(d)
    replicate(nrow(d), {
      closest <- which.min(d)
      d[,d_col[closest]] <<- NA
      d[d_row[closest],] <<- NA
      assignments[m[d_row[closest],"ind"]] <<- d_col[closest]
    })
  })
  assignments
}

## 2X over initial version, but already fast
assignments_to_components2 <- function(assignments, peaks) {
  sigma_ord <- order(peaks[,"sigma"])
  ord <- sigma_ord[order(assignments[sigma_ord], na.last = NA)]
  assign_f <- factor(assignments[ord])
  assign_tab <- table(assign_f)
  breaks <- cumsum(c(1, assign_tab))
  ir <- IRanges(head(breaks, -1), tail(breaks, -1) - 1)
  peaks <- peaks[ord,]
  ## views <- Views(Rle(peaks[,"maxf"]), ir)
  ## refs <- viewWhichMaxs(views)
  sharp_ir <- IRanges(start(ir), w = assign_tab/4)
  ## reference is peak w/ first quartile sigma
  refs <- start(ir) + assign_tab/4 
  views <- Views(Rle(peaks[,"rt"]), ir)
  rt_mean <- rep(viewSums(views)/assign_tab, assign_tab)
  views <- Views(Rle((peaks[,"rt"] - rt_mean)^2), ir)
  sigma_t <- sqrt(viewSums(views) / (assign_tab-1))
  cbind(comp = seq_along(refs), rt = peaks[refs, "rt"],
        sigma = peaks[refs,"sigma"], sigma_t = sigma_t,
        npeaks = assign_tab)
}

assignments_to_components <- function(assignments, peaks) {
  comps <- t(sapply(unique(assignments[!is.na(assignments)]), function(assignment) {
    g <- which(assignments == assignment)
    # try just using the sigma/mu of highest peak as comp parameters
    highest <- which.max(peaks[g, "maxf"])
    c(comp = assignment, peaks[g[highest],"rt"], 
      abs(peaks[g[highest],"sigma"]),
      sigma_t = sd(peaks[g, "rt"]),
      npeaks = length(g))
    #comps <<- rbind(comps, cbind(comp = assignment, mu = mean(peaks[g,"mu"]), 
    #  sigma = mean(abs(peaks[g,"sigma"]))))
  }))
  comps <- comps[order(comps[,1]),]
  #comps <- data.frame(comp = as.character(comps[,1]), comps[,-1])
  return(comps)
}

resolve_orphans <- function(comps, peaks) {
  assignments <- peaks[,"comp"]
  orphans <- is.na(assignments)
  if (any(orphans)) {
    assignments[orphans] <- max(assignments, na.rm=TRUE) + seq_len(sum(orphans))
    comps <- rbind(comps, assignments_to_components(assignments[orphans], peaks[orphans,]))
  }
  list(comps = comps, assignments = assignments)
}

