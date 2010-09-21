oneStep <- function(dir_design,edit=F,mode='CML',diagnose=F){
  file <- list.files(dir_design[1],full.names=T)[1]
  raw <- loadSample(file)
  raw_prof <- genProfile(raw)
  ## explore(raw_prof)
  ## prof <- raw_prof@env$profile
  ## if(is.na(missing)) missing <- min(prof[prof!=0])
  ## prof2 <- as.numeric(prof)
  ## prof2[prof2==0]=missing
  ## prof <- matrix(prof2,nrow(prof),ncol(prof))
  ## raw_prof@env$profile <- prof
  cor_prof <- removeBaseline(raw_prof, "median", scanrad = 100)
  raw <- cor_prof
  peaks <- findPeaks(raw, "gauss")
  s1_exp <- loadExperiment(dir_design[1], pipeline = pipeline(peaks))
   samples <- dir_design
  ## xset_comps <- findComps(s1_exp,"sigma_filt",tic.cutoff=0,npeaks.cutoff=0,tic=0,npeaks=10)
  xset_comps <- findComps(s1_exp,"sigma")
   s1_s2_xset<- perform(pipeline(xset_comps),samples)
  #s1_s2_xset_filt<- perform(pipeline(xset_comps_filt),samples)
  ## xset_groups <- groupComps(s1_s2_xset,"angle",rt_window=1,dist.cutoff=0.05)
  ## xset_sum <- summarize(xset_groups, "common")set_groups_filt<- groupComps(s1_s2_xset_filt,"angle",rt_window=2,dist.cutoff=0.05)
  ## xset_sum <- summarize(xset_groups_filt, "common")
 # s1_s2_xset<- compsFilter(s1_s2_xset,tic.cutoff=0,npeaks.cutoff=0,tic=0,npeaks=15)
   xset_groups <- groupComps(s1_s2_xset,"angle",rt_window=-1,dist.cutoff=0.05)
#  browser()
  xset_rtcor <- rtcor(xset_groups, "rloess")
  ## cplotRT(xset_rtcor0,xscale=c(3350,3400))
  ## explore(xset_rtcor,raw=xset_groups)
  xset_sum <- summarize(xset_rtcor, "common")
  xset_norm <- normalize(xset_sum, "scale")
}




oneStep2 <- function(dir_design){
  load('/home/tengfei/Desktop/back01/final.rda')
  samples <- dir_design
  xset_norm <- perform(pipeline(final),samples)
}




