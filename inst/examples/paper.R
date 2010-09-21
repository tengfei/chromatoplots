#source("load.R")
library(chromatoplots)

# load raw data
cdf <- ".*\\.CDF$"
genotypes <- c("Control", "Const")
files <- dir(file.path("cdf/LD", genotypes), pattern = cdf, full.names = TRUE)
raws <- lapply(files, xcmsRaw, profmethod="bin")
raw <- raws[[5]]

# retrieve/create commonly used data elements
prof <- raw@env$profile
time <- raw@scantime
intensity <- prof[2,]
chrom_51 <- data.frame(time = time, intensity = intensity)

# Profile matrix image

gen_prof_image <- function(raw, filename)
{
  profile <- raw@env$profile
  image_d <- data.frame(time = rep(time, each = nrow(profile)), 
    mz = rep(profMz(raw), ncol(profile)), intensity = as.vector(profile))
  image_d <- image_d[profile > 150,]
  p <- ggplot(image_d, aes(x = time, y = mz, fill = intensity)) + geom_tile() + 
    scale_fill_gradient(low = "yellow", high = "black", trans = "log")
  ggsave(p, filename, width=12, height=8)
  p
}

p_profmat_raw <- gen_prof_image(raw, "profile-image.pdf")

# Raw chromatogram at 51 m/z

p <- ggplot(chrom_51, aes(x = time, y = intensity)) + geom_point()
p$title <- NULL
#p$title <- "Extracted Ion Chromatogram at 51 m/z"
ggsave(p, "chrom-51-raw.pdf", width = 9, height = 4)

# Smoothers with residuals underneath

# RBE

p_rbe <- p + stat_smooth(method = robust_loess)
ggsave(p_rbe, "chrom-51-rbe-fit.pdf",  width = 9, height = 4)
model <- robust_loess(intensity ~ time)
residuals <- residuals(model)
ggsave(qplot(time, residuals), "chrom-51-rbe-res.pdf",  
  width = 9, height = 4)

# Median filter

raw <- removeProfileBaseline(raw, "median", scanrad = 100)
prof_cor <- raw@env$profile
cor <- prof_cor[2,]
p_chrom_med <- p + geom_line(x = time, y = intensity - cor, 
  colour = "gray50", size = 2)
ggsave(p_chrom_med, "chrom-51-median-fit.pdf",  width = 9, height = 4)
p_chrom_res <- qplot(time, cor, ylab = "residuals")
ggsave(p_chrom_res, "chrom-51-median-res.pdf",  width = 9, height = 4)

# Corrected Profile Matrix

p_profmat_cor <- gen_prof_image(raw, "profile-image-cor.pdf")

# peak fits

peaks <- findPeaks(raw, "islands", egh = FALSE)
peak <- peaks[2520,]
plot_peak(peak, raw, TRUE, quant)
dev.copy(pdf, file = "peak-fit-gauss.pdf", width = 4, height = 8)
dev.off()

peaks <- findPeaks(raw, "islands")
peak <- peaks[2524,]
plot_peak(peak, raw, TRUE, quant)
dev.copy(pdf, file = "peak-fit-egh.pdf", width = 4, height = 8)
dev.off()

# Quantile cutoffs

quant <- quantile(prof_cor, .99)

# at 51 m/z
p_quant_51 <- p_chrom_res + geom_hline(intercept = quant, colour = "red", size = 2)
p_quant_51$title <- "Baseline corrected chromatogram at 51 m/z"
peaks_51 <- peaks[peaks[,"mz"] == 51,]
p_max_51 <- p_quant_51 + 
  geom_point(size = 4, colour = "blue", data = data.frame(peaks_51), 
  aes(x = rt, y = maxf))
ggsave(p_max_51, "chrom-51-quant.pdf",  width = 9, height = 4)

# at 81 m/z
mz <- 32
p_quant_81 <- qplot(time, prof_cor[mz,], ylab = "residuals") + 
  geom_hline(intercept = quant, colour = "red", size = 2)
p_quant_81$title <- "Baseline corrected chromatogram at 81 m/z"
peaks_81 <- peaks[peaks[,"mz"] == mz+49,]
p_max_81 <- p_quant_81 + 
  geom_point(size = 4, colour = "blue", data = data.frame(peaks_81), 
  aes(x = rt, y = maxf))
ggsave(p_max_81, "chrom-81-quant.pdf",  width = 9, height = 4)

# for looking at non-detects (72/73 m/z)
print(p_max_81 + scale_x_continuous(limits=c(1000,1500)))

# Build our xset
# FIXME: need way to retrieve pipeline from peak result
pipeline <- new("xcmsPipeline", rawpipeline = pipeline(raw), 
  findpeaksproto = xcmsProtocol("findPeaks", "islands"))
# xset_peaks <- xcmsSet(files, pipeline = pipeline)
# colnames(phenoData(xset_peaks))[1] <- "genotype"
# save(xset_peaks, file = "xset_peaks.rda")
load("xset_peaks.rda")

# still focus on the first sample
peaks <- xset_peaks@peaks[xset_peaks@peaks[,"sample"] == 5,]

# peak image

peak_df <- as.data.frame(peaks)
p_peaks <- ggplot(peak_df) + 
  geom_point(data = peak_df, aes(x = rt, y = mz, colour = maxf)) + 
  scale_colour_gradient(low = "yellow", high = "black", trans = "log")
ggsave(p_peaks, "peak-image.pdf",  width=12, height=8)

# component image

#xset_comps <- findComps(xset_peaks, "sigma")
# save(xset_comps, file = "xset_comps.rda")
load("xset_comps.rda")

peaks <- xset_comps@peaks[xset_comps@peaks[,"sample"] == 5,]

e <- assignments_to_edges(peaks[,"comp"])
e_df <- data.frame(start = peaks[e[,1],], stop = peaks[e[,2],])
p_comps <- p_peaks + geom_segment(data = e_df, colour = alpha("gray50", 0.7), 
  aes(x = start.rt, y = start.mz, xend = stop.rt, yend = stop.mz))
ggsave(p_comps, file = "component-image.pdf",  width = 12, height = 8)

# plot the components by their sigma and sigma_t

p_comps_sigma <- qplot(sigma, sigma_t, data = as.data.frame(xset_comps@comps), 
  xlab = expression(sigma), ylab = expression(sigma[t]))
ggsave(p_comps_sigma, "components-sigma.pdf",  width = 6, height = 6)

# histograms of distances

spectra <- find_spectra(xset_comps@peaks[xset_comps@peaks[,"sample"] %in% 1:4,])
#spectra <- log(spectra+1)
dists <- hopach_discosangle(t(spectra)) # actual distances

library(MASS)
sigma <- matrix(0, nrow(spectra), nrow(spectra))
diag(sigma) <- 1
mv <- mvrnorm(nrow(dists), rep(0, nrow(spectra)), sigma)
# take the log?
#mv[mv > 0] <- log(mv[mv > 0])
#mv[mv < 0] <- -log(abs(mv[mv < 0]))
dists_mv <- hopach_discosangle(mv)

prep_dists <- function(dists)
{
  # ignore distances to self
  diag(dists) <- NA
  distance <- as.vector(dists)
  #distance <- log(log(distance + 2))
  #distance[!is.na(distance)]
  distance
}

distance <- prep_dists(dists)
distance_mv <- prep_dists(dists_mv)

qplot(distance, geom = "histogram", binwidth = .01)
#  ylim = c(0, 1.25e6), xlim = c(0,1))

qplot(distance_mv, geom = "histogram", binwidth = .01) 
#  ylim = c(0, 1.25e6), xlim = c(0,1)) 

# compare distance distributions with box plot

both_distances <- data.frame(distance = c(distance, distance_mv), 
  source = rep(c("actual", "random"), each = length(distance)))
qplot(source, distance, data = both_distances, geom = "boxplot")

# grouping

#xset_group <- group(xset_comps, "angle")
#save(xset_group, file = "xset_group.rda")
load("xset_group.rda")

groups_none <- as.data.frame(xset_group@groups)
groups_none$ncomps <- factor(groups_none$ncomps)

# relationship between distance measures
p <- qplot(sigma_t, mu_d, data = groups_none, ylab=expression(mu[d]), 
  xlab=expression(sigma[t]))
ggsave(p, "group-comparison.pdf",  width = 6, height = 6)

# some histograms (not that interesting)
qplot(mu_d, data = groups_none, geom="histogram", binwidth=0.01, xlim=c(0,1))
qplot(sigma_t, data = groups_none, geom="histogram", binwidth=10)

# robust loess retention time correction fits are output with the correction

#xset_retcor <- retcor(xset_group, "robust")
#save(xset_retcor, file = "xset_retcor.rda")
load("xset_retcor.rda")

# show correction for the raw data

tic_overlay <- function(time, label) {
  init <- data.frame(time = time[[1]], tic = raws[[1]]@tic)
  p <- ggplot(init, aes(x = time, y = tic)) + geom_line() +
    scale_x_continuous(limits=c(360,500)) +
    scale_y_continuous(limits=c(0,5e6))
  for (i in tail(seq_along(raws), -1))
  {
    p <- p + geom_line(x = time[[i]], y = raws[[i]]@tic)
  }
  p$title <- paste("TIC Chromatograms (", label, ")", sep="")
  ggsave(p, file = paste("tic-overlay-", label, ".pdf", sep=""), 
    width = 9, height = 4)
  p
}

raw_time <- xset_retcor@rt$raw
corrected_time <- xset_retcor@rt$corrected

p_tic_raw <- tic_overlay(raw_time, "raw")
p_tic_cor <- tic_overlay(corrected_time, "corrected")

# plot showing all corrections
init <- data.frame(time = range(raw_time[[1]]), 
  dev = range(unlist(corrected_time) - unlist(raw_time)))
p <- ggplot(init, aes(x = time, y = dev)) # + geom_line()
for (i in seq_along(raws))
{
  df <- data.frame(time = raw_time[[i]], dev = corrected_time[[i]] - raw_time[[i]])
  p <- p + geom_line(data = df, aes(x = time, y = dev))
}

# second matching

#xset_group_rt <- group(xset_retcor, "angle", rt_window = 1)
#save(xset_group_rt, file="xset_group_rt.rda")
load("xset_group_rt.rda")

groups_rt <- as.data.frame(xset_group_rt@groups)
groups_rt$ncomps <- factor(groups_rt$ncomps)

# component chain plot
comps <- xset_group_rt@comps
rownames(comps) <- NULL # to avoid duplicate rownames
e <- assignments_to_edges(comps[,"group"])
comps_df <- as.data.frame(comps)
e_df <- data.frame(start = comps_df[e[,1],], stop = comps_df[e[,2],])
p_groups <- qplot(rt, sample, data = comps_df) + 
  geom_segment(data = e_df, colour = alpha("gray50", 0.7), 
    aes(x = start.rt, y = start.sample, xend = stop.rt, yend = stop.sample))
ggsave(p_groups, file = "group-image.pdf",  width = 12, height = 8)

# need to compare groups to those generated without an rt window

old_new <- factor(rep(c("None", "Size 1"), c(nrow(groups_none), nrow(groups_rt))))
group_comp <- cbind(rbind(groups_none, groups_rt), window = old_new)
p <- qplot(ncomps, data = group_comp, geom = "bar", fill = window, position = "dodge")
ggsave(p, "groups-ncomps.pdf",  width = 9, height = 4)
#group_comp <- group_comp[!is.na(group_comp[,"mu_d"]),]
#qplot(window, mu_d, data = group_comp, geom = "boxplot", ylab = expression(mu[d]))

# histogram of quantities

# unnormalized

quant_hist <- function(xset, scale)
{
  quants <- split(xset@comps[,"quantity"], xset@comps[,"sample"])
  p <- ggplot(data.frame(quantity = quants[[1]]), aes(x = quantity))
  if (!missing(scale))
    p <- p + scale
  for (quant in quants)
    p <- p + geom_density(data = data.frame(quantity = quant), aes(x = quantity))
  p
}

#xset_sum <- summarize(xset_group_rt, "common")
#save(xset_sum, file="xset_sum.rda")
load("xset_sum.rda")

p_sum <- quant_hist(xset_sum, scale_x_log())
ggsave(p_sum, file = "quantities-raw.pdf",  width=6, height=6)

# normalized

#xset_norm <- normalize(xset_sum, "scale")
#save(xset_norm, file = "xset_norm.rda")
load("xset_norm.rda")

p_norm <- quant_hist(xset_norm)
ggsave(p_norm, file = "quantities-norm.pdf",  width=6, height=6)

#common_groups <- tapply(xset_norm@comps[,"sample"], xset_norm@comps[,"group"], 
#  function(samples) all(5:8 %in% samples))

xset_quant <- xset_sum
#xset_quant <- xset_norm

common_groups <- which(xset_quant@groups[,"ncomps"] == length(files))
comps <- xset_quant@comps[order(xset_quant@comps[,"group"]),]
common_comps <- comps[,"group"] %in% as.numeric(common_groups)
result <- do.call("cbind", split(comps[common_comps,"quantity"], 
  comps[common_comps, "sample"]))
#colnames(result) <- sapply(files, basename)
colnames(result) <- paste(rep(c("MUT", "WT"),each=4), 1:4, sep=".")
result <- result[,5:8]

pdf(file = "quantity-pairs.pdf")
pairs(log(result))
dev.off()


plotmatrix(as.data.frame(result))

# diagnostic visualizations (to be incorporated into chromatoplots)

stack_plots <- function(title, plots, weights = rep(1, length(plots)))
{
  win <- gtkWindow(show = FALSE)
  win$setTitle(paste("Chromatoplots [Stage: ", title, "]", sep=""))
  win$setDefaultSize(800, 600)
  win$show()
  vbox <- gtkVBox()
  win$add(vbox)
  sizes <- weights / sum(weights) * win$getDefaultSize()$height
  for (i in seq_along(plots)) {
    plot <- plots[[i]]
    if (inherits(plot, "ggplot"))
      wid <- gtkDrawingArea()
    else # ggobi plot
      wid <- plot
    wid$showAll()
    vbox$add(wid)
    wid$setSizeRequest(-1, sizes[i])
    if (inherits(plot, "ggplot")) {
      asCairoDevice(wid)
      print(plot)
    }
  }
}

# baseline subtraction

# draw a line to make it look like we selected it
# may not be possible with geom_tile()
#p_profmat_raw_sel <- p_profmat_raw + geom_hline(intercept = 2, colour = "red")

p_chrom_med$title <- NULL
stack_plots("Baseline Removal",
  list(p_profmat_raw_sel, p_chrom_med, p_chrom_res), weights = c(2,1,1))

# peak detection

peaks <- xset_comps@peaks[xset_comps@peaks[,"sample"] == 5,]
pipeline <- new("xcmsPipeline", rawpipeline = pipeline(raw), 
  findpeaksproto = xcmsProtocol("findPeaks", "islands"))
gg <- explore(findPeaksProto(pipeline), raw, peaks, FALSE, dev = 3)

# create peak fit display GUI using GGobi

peak_da <- gtkDrawingArea()
asCairoDevice(peak_da)
gg_peaks <- display(gg[1], embed = TRUE)
variables(gg_peaks) <- list(X = "rt", Y = "mz")
imode(gg_peaks) <- "Identify"

peak_hbox <- gtkHBox()
peak_hbox$add(gg_peaks)
peak_hbox$add(peak_da)

# clear title of chromatogram
p_max_51$title <- NULL

stack_plots("Peak Detection", list(p_profmat_cor, p_max_51, peak_hbox))

# component detection

comp_da <- gtkDrawingArea()
asCairoDevice(comp_da)
explore(xset_comps, 5, gg, 2)

comps_d <- gg["comps"]
peaks_d <- gg["peaks"]
e_comps_d <- gg["e_comps"]
gSignalConnect(gg, "identify-point", function(gg, plot, id, dataset) {
  if (id == -1)
    return()
  if (!(as.RGtkObject(comps_d) == dataset))
    return()
  glyph_color(peaks_d) <- 1
  glyph_color(peaks_d)[peaks[,"comp"] %in% comps_d[id+1,"comp"]] <- 9
  # no way to color edges
  #glyph_color(e_comps_d) <- 1
  #glyph_color(e_comps_d)[e_comps_d[,"comp"] %in% comps_d[id+1,"comp"]] <- 9
})

#parent <- peak_hbox$getParent()
#parent$remove(peak_hbox)

gg_comps <- display(gg["comps"], embed = TRUE)
variables(gg_comps) <- list(X = "sigma", Y = "sigma_t")
imode(gg_comps) <- "Identify"

comp_hbox <- gtkHBox()
comp_hbox$add(gg_comps)
comp_hbox$add(comp_da)

edges(gg_peaks) <- gg["e_comps"]

#comp_outlier <- which.max(gg["comps"][,"sigma",drop=TRUE])

stack_plots("Component Detection", list(comp_hbox, peak_hbox))

# grouping

group_id <- levels(factor(xset_group@comps[,"group"]))
#group_id <- rep(as.numeric(levels(group_f)), table(group_f))
groups_none <- cbind(group = group_id, groups_none)
gg["groups"] <- groups_none
gg_groups <- display(gg["groups"], embed = TRUE)
variables(gg_groups) <- list(X = "sigma_t", Y = "mu_d")

comp_order <- order(xset_group@comps[,"sample"])
assignments <- xset_group@comps[comp_order,"group"]
e_groups <- matrix(as.character(assignments_to_edges(assignments)), ncol=2)
edges(gg) <- e_groups
comps <- xset_group@comps[comp_order,]
rownames(comps) <- seq_len(nrow(comps))
gg["all_comps"] <- comps
gg_comp_chain <- display(gg["all_comps"], embed = TRUE)
variables(gg_comp_chain) <- list(X = "rt", Y = "sample")
edges(gg_comp_chain) <- gg["e_groups"]

# eventually need dataset with vars for m/z and rows for component
# exclusion specification and plot variables will change to show groups

comps <- xset_group@comps
comps <- comps[order(comps[,"comp"]),]
comps <- comps[order(comps[,"sample"]),]
spectra <- t(find_spectra(xset_comps@peaks))
#spectra <- log(spectra + 1)
colnames(spectra) <- paste("mz", seq_len(ncol(spectra)), sep="")
gg["spectra"] <- as.data.frame(spectra)
spectra_d <- gg["spectra"]
gg_comp_peaks <- display(spectra_d, "Parallel Coordinates Display", list(X=colnames(spectra)[1]), 
  embed = TRUE)
#variables(gg_comp_peaks) <- list(X=colnames(spectra)[1])

#group <- 17 # group to show in bottom plot
#spectra_g <- spectra[comps[,"group"] == group,]
#non_empty <- which(colSums(spectra_g) > 0)
#spectra_g <- spectra_g[,non_empty]
#colnames(spectra_g) <- paste("mz", non_empty, sep="")
#gg["spectra"] <- as.data.frame(spectra_g)
#gg_comp_peaks <- display(gg["spectra"], "Parallel Coordinates Display", embed = TRUE)
#variables(gg_comp_peaks) <- list(X=c(colnames(spectra_g)))

groups_hbox <- gtkHBox()
groups_hbox$add(gg_groups)
groups_hbox$add(gg_comp_chain)

groups_d <- gg["groups"]
all_comps_d <- gg["all_comps"]

show_comp_peaks <- function(group)
{
  comps_g <- comps[,"group"] == group
  spectra_g <- spectra[comps_g,,drop=FALSE]
  spectra_g <- spectra_g[,which(colSums(spectra_g) > 0),drop=FALSE]
  mz_vis <- rep(TRUE, nrow(spectra))
  mz_vis[comps_g] <- FALSE
  excluded(spectra_d) <- mz_vis
  mz_order <- order(colMeans(spectra_g), decreasing = TRUE)
  variables(gg_comp_peaks) <- list(X=colnames(spectra_g)[mz_order])
}

gSignalConnect(gg, "identify-point", function(gg, plot, id, dataset) {
  if (id == -1)
    return()
  if (!(as.RGtkObject(groups_d) == dataset))
    return()
  group <- id+1
  glyph_color(all_comps_d) <- 1
  group_match <- all_comps_d[,"group",drop=TRUE] %in% 
    groups_d[group,"group", drop=TRUE]
  glyph_color(all_comps_d)[group_match] <- 9
  shadow <- rep(TRUE, nrow(all_comps_d))
  shadow[group_match] <- FALSE
  shadowed(all_comps_d) <- shadow
  # bug in ggobi prevents this
  #show_comp_peaks(group)
})

gSignalConnect(gg, "identify-point", function(gg, plot, id, dataset) {
  if (id == -1)
    return()
  if (!(as.RGtkObject(all_comps_d) == dataset))
    return()
  comp <- all_comps_d[id+1,,drop=TRUE]
  comps_g <- which(comps[,"group"] == comp$group)
  glyph_color(spectra_d) <- 1
  glyph_color(spectra_d)[comps_g[comps[comps_g, "comp"] == comp$comp]] <- 9
})

# don't forget to scale the parallel coordinate plot
stack_plots("Component Grouping", list(groups_hbox, gg_comp_peaks))

# retention time correction

gg_comp_chain <- display(gg["all_comps"], embed = TRUE)
variables(gg_comp_chain) <- list(X = "rt", Y = "sample")
edges(gg_comp_chain) <- gg["e_groups"]

s <- xset_group@comps[,"sample"] == 5
retcor_d <- data.frame(raw = xset_group@comps[s,"rt"], cor = xset_retcor@comps[s,"rt"])
p_retcor_fit <- qplot(raw, cor-raw, data = retcor_d) + 
 + scale_x_continuous(limits=c(360,1000)) + geom_smooth()
 
p_tic_cor$title <- NULL

stack_plots("Retention Time Correction", list(p_tic_cor, p_retcor_fit, gg_comp_chain))

# comps after rt correction

comps <- xset_group_rt@comps
comps <- comps[order(comps[,"comp"]),]
comps <- comps[order(comps[,"sample"]),]

# summarization

result <- matrix(NA, nrow(xset_sum@groups), nrow(phenoData(xset_sum)))
comps_n <- xset_sum@comps[!is.na(xset_sum@comps[,"quantity"]),]
group_n <- as.numeric(factor(comps_n[,"group"]))
result[(comps_n[,"sample"]-1)*nrow(result) + group_n] <- log(comps_n[,"quantity"]+1)
bad_groups <- apply(result, 1, function(x) all(is.na(x)))
result_groups <- comps_n[!bad_groups,"group"]
result <- result[!bad_groups,]

norm_diag <- t(apply(result, 1, function(g) c(mean = mean(g, na.rm=TRUE), 
  max_diff = max(abs(outer(g, g, "-")), na.rm=TRUE), 
  ncomps = sum(!is.na(g)))))

colnames(result) <- sapply(files, basename)
norm_diag <- as.data.frame(norm_diag)
norm_diag$ncomps <- factor(norm_diag$ncomps)

gg["result"] <- result
gg["result_diag"] <- norm_diag
gg_result_diag <- display(gg["result_diag"], embed=TRUE)
gg_result_bar <- display(gg["result_diag"], "Barchart", 
  list(X="ncomps"), embed=TRUE)
variables(gg_result_bar) <- list(X="ncomps")
gg_result <- display(gg["result"],  "Parallel Coordinates Display", embed=TRUE)
gg_comp_peaks <- display(spectra_d, "Parallel Coordinates Display", 
  list(X=colnames(spectra)[1]), embed = TRUE)

result_hbox <- gtkHBox()
result_hbox$add(gg_result_diag)
result_hbox$add(gg_result_bar)

result_diag_d <- gg["result_diag"]
result_d <- gg["result"]
gSignalConnect(gg, "identify-point", function(gg, plot, id, dataset) {
  if (id == -1)
    return()
  if (!(as.RGtkObject(result_diag_d) == dataset))
    return()
  result_ex <- rep(TRUE, nrow(result_d))
  result_ex[id+1] <- FALSE
  excluded(result_d) <- result_ex
  #glyph_color(result_d) <- 1
  #glyph_color(result_d)[id+1] <- 9
})

stack_plots("Summarization", list(result_hbox, gg_result, gg_comp_peaks))

# normalization

# pff, same as above

