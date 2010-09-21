groupComps.angle <- function(object, rt_window = -1,dist.cutoff=0.05)
{
  matching <- match_components(object@comps, object@peaks, 
                               object@phenoData, , rt_window, dist.cutoff)
  comps <- object@comps[matching[,"id"], # sort, drop 'group'
                        colnames(object@comps) != "group"]
  object@comps <- cbind(comps, group = matching[,"cluster"])
  ## need a 'group' dataset here with distances, time variance, etc
  object@groups <- matchings_to_groups(object@comps, object@peaks)
  object
}

match_components <- 
function(comps, peaks, design = matrix(), subset = seq_len(nrow(design)), 
  rt_window = -1, dist.cutoff=1, maxdist = 0.9)
{
  # remember indices (past subsetting)
  if (!("id" %in% colnames(comps)))
    comps <- cbind(comps, id = seq_len(nrow(comps)))
  
  if (is.logical(subset)) # ensure subset are indices
    subset <- which(subset)
  
  if (ncol(design) > 0) { # recurse along first factor to get matchings
    match_subset <- function(ind)
      match_components(comps, peaks, design[,-1,drop=F], ind, rt_window,dist.cutoff)
    sample_ind <- split(seq_len(nrow(design))[subset], design[subset,1])
    children <- lapply(sample_ind, match_subset)
    #browser()
    find_median_comp <- function(cluster) 
      cluster[which.median(comps[cluster,"rt"])]
    children_split <- lapply(children, function(m) split(m[,1], m[,"cluster"]))
    median_list <- lapply(children_split, function(child) 
      sapply(child, find_median_comp))
    medians <- unlist(median_list)
    # rearrange peak and component matrices for the median "metacomponents"
    # medians need to be sorted to index into peak splitting below
    comps <- comps[order(comps[,"comp"]),]
    comps <- comps[order(comps[,"sample"]),]
    medians <- match(medians, comps[,"id"])
    peak_ind <- split(seq_len(nrow(peaks)), 
      interaction(peaks[,"comp"], peaks[,"sample"]))
    existing_ints <- sapply(peak_ind, length) > 0
    peak_ind <- peak_ind[existing_ints][medians]
    peaks <- peaks[unlist(peak_ind),]
    comps <- comps[medians,]
    meta_samples <- rep(seq_along(children), sapply(median_list, length))
    comps[,"sample"] <- meta_samples
    comps[,"comp"] <- seq_along(medians)
    peaks[,"sample"] <- rep(comps[,"sample"], sapply(peak_ind, length))
    peaks[,"comp"] <- rep(comps[,"comp"], sapply(peak_ind, length))
  } else { # subset here
    peaks <- peaks[peaks[,"sample"] %in% subset,]
    #comp_id <- comp_id[comps[,"sample"] %in% subset]
    comps <- comps[comps[,"sample"] %in% subset,]
  }
  
  comps <- comps[order(comps[,"comp"]),]
  comps <- comps[order(comps[,"sample"]),]
  
  # generates mass spectrum for a component
  spectra <- find_spectra(peaks)
  
  #spectra <- log(spectra+1)
  dists <- hopach_discosangle(t(spectra))
  
  #print(quantile(dists, 0.05))
  #dists[dists > maxdist] <- NA
  # instead of maxdist, why not a quantile?
  #dists[dists > quantile(dists, 0.05)] <- NA
  dists[dists > quantile(dists,dist.cutoff)] <- NA

  if (rt_window >= 0) {
    group_sep <- abs(outer(comps[,"rt"], comps[,"rt"], "-"))
    dists[group_sep > rt_window] <- NA
  }
  gc()
  
  # cluster the components
  clusters <- cluster_between(dists, comps[,"sample"])
  
  #browser()
  
  # map ID's from subset to original
  matching <- cbind(id = comps[,"id"], cluster = clusters)
  if (ncol(design) > 0) {
    matching_split <- split.data.frame(matching, comps[,"sample"])
    matching <- do.call("rbind", lapply(seq_along(children), function(child)
    {
      matching_sub <- matching_split[[child]]
      children_sub <- children[[child]]
      id_match <- match(matching_sub[,"id"], children_sub[,"id"])
      child_clusters <- children_split[[child]][as.character(children_sub[id_match,"cluster"])]
      id <- unlist(child_clusters)
      cluster <- rep(matching_sub[,"cluster"], sapply(child_clusters,length))
      cbind(id = id, cluster = cluster)
    }))
  }

  matching
}

cluster_between <- function(distances, grouping)
{
  # first, block matchings within groups
  grouping <- as.factor(grouping)
  groups <- split(seq_len(nrow(distances)), grouping)
  for (group in groups)
    distances[group,group] <- NA
  distances[lower.tri(distances)] <- NA # one distance per matching, please
  dist_col <- col(distances)
  dist_row <- row(distances)
  dist_order <- order(distances)
  clusters <- seq_along(grouping)
  for (i in dist_order) {  
    if (!is.na(distances[i])) {
      a <- dist_col[i]
      b <- dist_row[i]
      if (clusters[a] != a) {
        cluster <- clusters[a]
        cand <- b
      } else {
        cluster <- clusters[b]
        cand <- a
      }
      if (!is.na(distances[cluster,cand]) || 
          !is.na(distances[cand,cluster])) {
        # do not allow any more from candidate's group into cluster
        cand_group <- groups[[grouping[cand]]]
        distances[cluster, cand_group] <- NA
        distances[cand_group, cluster] <- NA
        clusters[clusters==cand] <- cluster
      }
    }
  }
  # for now at least, return clusters
  clusters
}

# about the same as the spectral angle but a lot faster
hopach_discosangle<-function(X,na.rm=FALSE){
	X <- as.matrix(X)
	dX<-dim(X)
	p<-dX[1]
	n<-dX[2]
	if(na.rm){
		N<-rowSums(!is.na(X))
		N2<-(!is.na(X))%*%t(!is.na(X))
		X[is.na(X)]<-0
		N<-sqrt(N%*%t(N))/N2
	}
	else
		N<-1
	out<-rowSums(X^2)
	out<-1-abs(N*tcrossprod(X)/sqrt(tcrossprod(out)))
	diag(out)<-0
	#suppressWarnings(out<-sqrt(out))
	#out[out=="NaN"]<-0
	return(out)
}

spectral_angle <- function(a, b)
{
  acos(min(1,sum(a * b) / sqrt(sum(a^2) * sum(b^2))))
}

matchings_to_groups <- function(comps, peaks)
{
  # FIXME: need to make component id's globally unique or add groups to peaks
  peaks <- cbind(peaks, key = interaction(peaks[,"comp"], peaks[,"sample"]))
  comps <- cbind(comps, key = interaction(comps[,"comp"], comps[,"sample"]))
  do.call("rbind", tapply(comps[,"key"], comps[,"group"], function(group)
  {
    group_peaks <- peaks[peaks[,"key"] %in% group,,drop=FALSE]
    spectra <- find_spectra(group_peaks)
    #spectra <- log(spectra + 1)
    distance <- hopach_discosangle(t(spectra))
    diag(distance) <- NA
    sigma_d <- NA
    if (any(!is.na(distance)))
      sigma_d <- sd(as.vector(distance), na.rm = TRUE)
    c(ncomps = length(group), npeaks_sum = nrow(group_peaks), 
      npeaks_mean = nrow(group_peaks)/length(group), 
      sigma_t = sd(group_peaks[,"rt"]),
      sigma_d = sigma_d,
      mu_d = mean(distance, na.rm = TRUE))
      #group = group)
  }))
}
