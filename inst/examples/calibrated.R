library(chromatoplots)
cdfdir <- "~/Datas/mets/suhyeon/exp"
cdffiles <- list.files(cdfdir,recursive=TRUE,full.names=TRUE)
cdffiles
############################################################
## Start from one single cdf file
############################################################
## Step 1: Load Sample
raw1 <- loadSample(cdffiles[1])
## Step 2: Generate Profile Matrix
raw_prof1 <- genProfile(raw1)
## Step 3: Remove baseline
cor_prof1 <- removeBaseline(raw_prof1, "median", scanrad = 100)
## Step 4: Find Peaks
peaks1 <- findPeaks(cor_prof1, "gauss")
## Step 5: Load Experiment
setwd(cdfdir)
## ci_exp <- loadExperiment(c("s1", "s2"))
## save(ci_exp,file="~/Datas/rdas/ci_exp.rda")
library(chromatoplots)
load("~/Datas/rdas/ci_exp.rda")
## Step 6: Find Components
xset_comps <- findComps(ci_exp, "sigma")
comps <- xset_comps@comps
head(comps)
comps[comps[,"comp"]==23,]
pks <- xset_comps@peaks
head(pks)
by(pks,pks[,"sample"],function())

pks <- as.data.frame(pks)
inter <- interaction(pks[,"sample"],pks[,"comp"])
pks$inter <- inter
lst <- by(pks,pks$inter,function(obj){
  ## 50-800
  mz.name <- paste("mz",50:800,sep="")
  df <- t(data.frame(rep(0,length(50:800))))
  colnames(df) <-mz.name
  rownames(df) <- unique(obj$inter)
  nms <- paste("mz",obj$mz,sep="")
  df[,nms] <- obj$maxf
  df
})
df <- do.call("rbind",lst)
df <- as.data.frame(df)
df$id <- rownames(df)
comps <- as.data.frame(xset_comps@comps)
comps$id <- interaction(comps$sample,comps$comp)
head(comps)
mets <- merge(comps,df,by="id")
dim(mets)
save(mets,file="~/Desktop/mets.csv")
head(mets[,1:10])
idx <- pks[,"comp"]==140

head(pks)
## 
head(xset_comps@comps)
write.csv(xset_comps@peaks,file="~/Desktop/comps.csv")
## Step 7: Grouping Components
xset_groups <- groupComps(xset_comps, "angle")
## Step 8: Retentio Time Correction
xset_rtcor <- rtcor(xset_groups, "rloess")
## Step 9: regroup after rt correction
xset_groups2 <- groupComps(xset_rtcor, "angle", rt_window = 1)
## Step 10: Summarize by common peaks
xset_sum <- summarize(xset_groups2, "common")
## Step 11: Normalize (by simple scaling)
xset_norm <- normalize(xset_sum, "scale")
