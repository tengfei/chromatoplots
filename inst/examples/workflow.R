library(commandr)
library(chromatoplots)
library(xcms)
library(IRanges)
library(Biobase)
library(ggplot2)
library(RGtk2)
library(cairoDevice)
library(rggobi)
source("~/Prolang/svn/repos/trunk/chromatoplots/R/cplot.R")
# source("~/Prolang/svn/repos/trunk/chromatoplots/R/cpExperiment-class.R")
source("~/Prolang/svn/repos/trunk/chromatoplots/R/simul.R")
source("~/Prolang/svn/repos/trunk/chromatoplots/R/onestep.R")
#source("~/Prolang/svn/repos/trunk/chromatoplots/R/findPeaks-alg.R")
#source("~/Prolang/svn/repos/trunk/chromatoplots/R/findPeaks-plots.R")
#source("~/Prolang/svn/repos/trunk/chromatoplots/R/findPeaks-pipes.R")
#source("~/Prolang/svn/repos/trunk/chromatoplots/R/genProfile-plots.R")
#source("~/Prolang/svn/repos/trunk/chromatoplots/R/genProfile-pipes.R")
#source("~/Prolang/svn/repos/trunk/chromatoplots/R/graphics.R")
## source("~/Prolang/svn/repos/trunk/chromatoplots/R/findComps-alg.R")
## source("~/Prolang/svn/repos/trunk/chromatoplots/R/findComps-pipes.R")

file <- '/home/tengfei/Data/tengfei/sampleraw/baseline/a1.CDF'
raw <- loadSample(file)
raw_prof <- genProfile(raw)

setwd('~/Desktop')
png('scan100.png',800,600)
cplotViewBL2(cor_prof,raw=raw_prof,mz=51)
dev.off()
png('scan10.png',800,600)
cplotViewBL2(cor_prof10,raw=raw_prof,mz=51)
dev.off()
png('scan50.png',800,600)
cplotViewBL2(cor_prof50,raw=raw_prof,mz=51)
dev.off()

file <- '/home/tengfei/Data/tengfei/sampleraw/rtshift5/s1/a1.CDF'

setwd('/home/tengfei/test')
dir <- c('s1','s2')
final <- oneStep(dir)


cplotViewBL2(cor_prof10,raw=raw_prof,mz=51)
raw <- cor_prof

load('/home/tengfei/Data/tengfei/sampleraw/sample_2_low0.1_intsd0.1/space.rda')
file <- '/home/tengfei/Data/tengfei/sampleraw/sample_2_low0.1_intsd0.5/s1/a1.CDF'

setwd( '/home/tengfei/Data/tengfei/cal/raw')
load('space.rda')

setwd('~/Desktop')
load('d.rda')
save(final,file='d.rda')

dir_design <- c('s1','s2')
final <- oneStep(dir_design)
final <- fi

cplotSpec(final,group=944)

file <- '/home/tengfei/Data/tengfei/cal/raw/s1/M0101A.CDF'
raw <- loadSample(file)
raw_prof <- genProfile(raw)
cor_prof <- removeBaseline(raw_prof, "median", scanrad = 100)
raw <- cor_prof
peaks <- findPeaks(object=raw,'gauss')

explore(peaks, raw = raw, sample=NA,residuals=TRUE,island=TRUE)
cplotSpec(final,group=23)
str(final)

warnings()

file <- "/home/tengfei/Data/tengfei/cal/raw/s2/M0201A.CDF"
rep <- 4
missing <- NULL
low <- 0.3
back_sd <- 0.1

lst <- simulator(dir='~/Desktop',model=file,rt_shift_sd=5,rep=rep,missing=missing,low=low,back_sd=back_sd)



lst$rt

pdf('~/Desktop/rt.pdf',width=8,height=4)

png('~/Desktop/rt.png',800,400)
cplotTIC3('~/Desktop/test',xscale=c(2190,2230),facets=F)
dev.off()

setwd('~/Data/tengfei/sampleraw/')
newdir <- 'sample_2_low0.3_intsd0'
system(paste('mkdir',newdir,sep=' '))
setwd(paste('./',newdir,sep=''))
system('mkdir s1')
system('mkdir s2')
system('mv ~/Desktop/a*.CDF ./s1')
system('mv ~/Desktop/b*.CDF ./s2')
system('touch log.txt')
cat(paste('model=',file,'rep=',rep,'missing=',missing,'low=',low,'int_sd=',int_sd,sep='\n'),file='log.txt')
dir_design <- c('s1','s2')
final <- oneStep2(dir_design)



nrow(tdiff)
tdiff <- topDiff(final,pvalue=0.0001)
tdiff
save(final,file='ffnew.rda')
setwd('~/Desktop')
id <- rownames(tdiff)
png('diff.png',800,800)
cplotResult(final,id.diff=id,group.name=T)
args(cplotResult)
dev.off()
rt.ori <- lst$rt[lst$coding%in%c(1,-1,2,-2)]/60
rt.cur <- tdiff$rt
tcommon <- topCommon(final,scale=1:5)
id.com <- rownames(tcommon)
png('common.png',800,800)
cplotResult(final,id.common=id.com,group.name=T)
dev.off()
png('total.png',800,800)
cplotResult(final,id.diff=id,id.common=id.com,group.name=F)
dev.off()
rt.cur.com <- tcommon$rt
rt.ori.com <- lst$rt[lst$coding==0]/60
rtmatch(rt.ori,rt.cur,plot.match=F)
rtmatch(rt.ori.com,rt.cur.com,plot.match=F)
png('rt.png',800,200)
rtmatch(c(rt.ori,rt.ori.com),c(rt.cur,rt.cur.com),rt.ori,rt.ori.com,plot.match=T)
dev.off()
png('rtdiff.png',800,200)
rtmatch(rt.ori,rt.cur,rt.ori,plot.match=T)
dev.off()
#save(lst,final,file='space2.rda')
sink('log.txt',append=T)
topDiff(final)
topCommon(final)
sink()

exprs(as(final,'ExpressionSet'))
final@comps

rm(final)
rep <- 2
missing <- NULL
low <- 0.3
int_sd <- 0.001
lst <- simulator(dir='~/Desktop',model=file,rep=rep,missing=missing,low=low,int_sd=int_sd)
setwd('~/Data/tengfei/sampleraw/')
newdir <- 'sample_2_low0.3_intsd0.001'
system(paste('mkdir',newdir,sep=' '))
setwd(paste('./',newdir,sep=''))
system('mkdir s1')
system('mkdir s2')
system('mv ~/Desktop/a*.CDF ./s1')
system('mv ~/Desktop/b*.CDF ./s2')
system('touch log.txt')
cat(paste('model=',file,'rep=',rep,'missing=',missing,'low=',low,'int_sd=',int_sd,sep='\n'),file='log.txt')
dir_design <- c('s1','s2')
final <- oneStep2(dir_design)
tdiff <- topDiff(final,pvalue=0.0001)
id <- rownames(tdiff)
png('diff.png',800,800)
cplotResult(final,id.diff=id)
dev.off()
rt.ori <- lst$rt[lst$coding%in%c(1,-1,2,-2)]/60
rt.cur <- tdiff$rt
tcommon <- topCommon(final,scale=1:5)
id.com <- rownames(tcommon)
png('common.png',800,800)
cplotResult(final,id.common=id.com,group.name=T)
dev.off()
png('total.png',800,800)
cplotResult(final,id.diff=id,id.common=id.com,group.name=F)
dev.off()
rt.cur.com <- tcommon$rt
rt.ori.com <- lst$rt[lst$coding==0]/60
rtmatch(rt.ori,rt.cur,plot.match=F)
rtmatch(rt.ori.com,rt.cur.com,plot.match=F)
png('rt.png',800,200)
rtmatch(c(rt.ori,rt.ori.com),c(rt.cur,rt.cur.com),rt.ori,rt.ori.com,plot.match=T)
dev.off()
png('rtdiff.png',800,200)
rtmatch(rt.ori,rt.cur,rt.ori,plot.match=T)
dev.off()
save(lst,final,file='space.rda')
sink('log.txt',append=T)
topDiff(final)
topCommon(final)
sink()
rm(final)
rep <- 2
missing <- NULL
low <- 0.3
int_sd <- 0.01
lst <- simulator(dir='~/Desktop',model=file,rep=rep,missing=missing,low=low,int_sd=int_sd)
setwd('~/Data/tengfei/sampleraw/')
newdir <- 'sample_2_low0.3_intsd0.01'
system(paste('mkdir',newdir,sep=' '))
setwd(paste('./',newdir,sep=''))
system('mkdir s1')
system('mkdir s2')
system('mv ~/Desktop/a*.CDF ./s1')
system('mv ~/Desktop/b*.CDF ./s2')
system('touch log.txt')
cat(paste('model=',file,'rep=',rep,'missing=',missing,'low=',low,'int_sd=',int_sd,sep='\n'),file='log.txt')
dir_design <- c('s1','s2')
final <- oneStep2(dir_design)
tdiff <- topDiff(final,pvalue=0.0001)
id <- rownames(tdiff)
png('diff.png',800,800)
cplotResult(final,id.diff=id)
dev.off()
rt.ori <- lst$rt[lst$coding%in%c(1,-1,2,-2)]/60
rt.cur <- tdiff$rt
tcommon <- topCommon(final,scale=1:5)
id.com <- rownames(tcommon)
png('common.png',800,800)
cplotResult(final,id.common=id.com,group.name=T)
dev.off()
png('total.png',800,800)
cplotResult(final,id.diff=id,id.common=id.com,group.name=F)
dev.off()
rt.cur.com <- tcommon$rt
rt.ori.com <- lst$rt[lst$coding==0]/60
rtmatch(rt.ori,rt.cur,plot.match=F)
rtmatch(rt.ori.com,rt.cur.com,plot.match=F)
png('rt.png',800,200)
rtmatch(c(rt.ori,rt.ori.com),c(rt.cur,rt.cur.com),rt.ori,rt.ori.com,plot.match=T)
dev.off()
png('rtdiff.png',800,200)
rtmatch(rt.ori,rt.cur,rt.ori,plot.match=T)
dev.off()
save(lst,final,file='space.rda')
sink('log.txt',append=T)
topDiff(final)
topCommon(final)
sink()
rm(final)
rep <- 2
missing <- NULL
low <- 0.3
int_sd <- 0.1
lst <- simulator(dir='~/Desktop',model=file,rep=rep,missing=missing,low=low,int_sd=int_sd)
setwd('~/Data/tengfei/sampleraw/')
newdir <- 'sample_2_low0.3_intsd0.1'
system(paste('mkdir',newdir,sep=' '))
setwd(paste('./',newdir,sep=''))
system('mkdir s1')
system('mkdir s2')
system('mv ~/Desktop/a*.CDF ./s1')
system('mv ~/Desktop/b*.CDF ./s2')
system('touch log.txt')
cat(paste('model=',file,'rep=',rep,'missing=',missing,'low=',low,'int_sd=',int_sd,sep='\n'),file='log.txt')
dir_design <- c('s1','s2')
final <- oneStep2(dir_design)
tdiff <- topDiff(final,pvalue=0.0001)
id <- rownames(tdiff)
png('diff.png',800,800)
cplotResult(final,id.diff=id)
dev.off()
rt.ori <- lst$rt[lst$coding%in%c(1,-1,2,-2)]/60
rt.cur <- tdiff$rt
tcommon <- topCommon(final,scale=1:5)
id.com <- rownames(tcommon)
png('common.png',800,800)
cplotResult(final,id.common=id.com,group.name=T)
dev.off()
png('total.png',800,800)
cplotResult(final,id.diff=id,id.common=id.com,group.name=F)
dev.off()
rt.cur.com <- tcommon$rt
rt.ori.com <- lst$rt[lst$coding==0]/60
rtmatch(rt.ori,rt.cur,plot.match=F)
rtmatch(rt.ori.com,rt.cur.com,plot.match=F)
png('rt.png',800,200)
rtmatch(c(rt.ori,rt.ori.com),c(rt.cur,rt.cur.com),rt.ori,rt.ori.com,plot.match=T)
dev.off()
png('rtdiff.png',800,200)
rtmatch(rt.ori,rt.cur,rt.ori,plot.match=T)
dev.off()
save(lst,final,file='space.rda')
sink('log.txt',append=T)
topDiff(final)
topCommon(final)
sink()
rm(final)
rep <- 2
missing <- NULL
low <- 0.3
int_sd <- 0.5
lst <- simulator(dir='~/Desktop',model=file,rep=rep,missing=missing,low=low,int_sd=int_sd)
setwd('~/Data/tengfei/sampleraw/')
newdir <- 'sample_2_low0.3_intsd0.5'
system(paste('mkdir',newdir,sep=' '))
setwd(paste('./',newdir,sep=''))
system('mkdir s1')
system('mkdir s2')
system('mv ~/Desktop/a*.CDF ./s1')
system('mv ~/Desktop/b*.CDF ./s2')
system('touch log.txt')
cat(paste('model=',file,'rep=',rep,'missing=',missing,'low=',low,'int_sd=',int_sd,sep='\n'),file='log.txt')
dir_design <- c('s1','s2')
final <- oneStep2(dir_design)
tdiff <- topDiff(final,pvalue=0.11)
id <- rownames(tdiff)
png('diff.png',800,800)
cplotResult(final,id.diff=id)
dev.off()
rt.ori <- lst$rt[lst$coding%in%c(1,-1,2,-2)]/60
rt.cur <- tdiff$rt
tcommon <- topCommon(final,scale=1:5)
id.com <- rownames(tcommon)
png('common.png',800,800)
cplotResult(final,id.common=id.com,group.name=T)
dev.off()
png('total.png',800,800)
cplotResult(final,id.diff=id,id.common=id.com,group.name=F)
dev.off()
rt.cur.com <- tcommon$rt
rt.ori.com <- lst$rt[lst$coding==0]/60
rtmatch(rt.ori,rt.cur,plot.match=F)
rtmatch(rt.ori.com,rt.cur.com,plot.match=F)
png('rt.png',800,200)
rtmatch(c(rt.ori,rt.ori.com),c(rt.cur,rt.cur.com),rt.ori,rt.ori.com,plot.match=T)
dev.off()
png('rtdiff.png',800,200)
rtmatch(rt.ori,rt.cur,rt.ori,plot.match=T)
dev.off()
save(lst,final,file='space.rda')
sink('log.txt',append=T)
topDiff(final)
topCommon(final)
sink()
system('aplay -q /usr/share/skype/sounds/TransferComplete.wav ')




lst <- simulator(dir='~/Desktop',model=file,rep=2,missing=150,int_sd=0.01)

setwd('~/Desktop')

f <- '/home/tengfei/Desktop/a1.CDF'
f <- loadSample(f)
f <- genProfile(f)
setwd('~/Desktop')
png('test.png',600,600)
image(f)
dev.off()
f <- '/home/tengfei/Desktop/b1.CDF'
f <- loadSample(f)
f <- genProfile(f)
setwd('~/Desktop')
png('test2.png',600,600)
image(f)
dev.off()

setwd("/home/tengfei/Prolang/R/data/raw")
setwd('~/Desktop/sample3')
dir_design <- c('s1','s2')
final <- oneStep(dir_design)
tdiff <- topDiff(final,pvalue=0.0001)
tdiff
id <- rownames(tdiff)
rt.diff <- tdiff$rt
topCommon(final,scale=1:5)
rt.common <- topCommon(final,scale=1:5)$rt
xintercept <- list(rt.diff*60,rt.common*60)
png('result2.png',600,600)
cplotResult(final,id)
dev.off()
png('tic2.png',600,600)
cplotTIC(final,xintercept=xintercept)
dev.off()
save(final,lst,file='last.rda')

length(rt.diff)



lst3 <- simulator(dir='~/Desktop',model=file,rep=4,missing=150,rep_sd=0.1,int_sd=0.1,low=0.3,low_sd=0.1)
system('mv ~/Desktop/a*.CDF ~/Desktop/sample3/s1')
system('mv ~/Desktop/b*.CDF ~/Desktop/sample3/s2')
setwd('~/Desktop')
dir_design <- c('sample3/s1','sample3/s2')
final <- oneStep(dir_design)
tdiff <- topDiff(final,pvalue=0.05)
id <- rownames(tdiff)
rt.diff <- tdiff$rt
rt.common <- topCommon(final,scale=1:5)$rt
xintercept <- list(rt.diff*60,rt.common*60)
png('result3.png',600,600)
cplotResult(final,id)
dev.off()
png('tic3.png',600,600)
cplotTIC(final,xintercept=xintercept)
dev.off()

setwd('~/Desktop')
save(final,file='ff.rda')

setwd('~/Desktop')
png('tic.png',600,600)
cplotTIC(final,xintercept=xintercept)
dev.off()

setwd('~/Prolang/R/data/raw')
dir_design <- c('s1','s2')
final2 <- oneStep(dir_design)

setwd('~/Desktop')
save(final2,file='ff2.rda')



setwd(datapath)
load("~/peaks.rda")
load("./xset_comps.rda")
load('./s1_s2.rda')
load("./group.rda")
load("./sum_norm.rda")
load("./total_peaks.rda")


file <- "/home/tengfei/Data/tengfei/cal/raw/s2/M0201A.CDF"

file <- "/home/tengfei/Desktop/sample2/s1/a1.CDF"
file <- "/home/tengfei/Desktop/sample/s2/b1.CDF"
### load raw data, without a profile matrix
raw <- loadSample(file)
raw2 <- loadSample(file2)



## explore it
## too slow by using ggplot2(at least 1 min), please use qtpaint instead.
## and make it interactive, I remove the legend.
file="~/Desktop/M0201A.CDF"
explore(raw)
str(raw)## generate a profile matrix using simple binning, step size of 1
raw_prof <- genProfile(raw)

## image for model(profile matrix)
rt <- raw_prof@scantime
missing <- 150
mz_range <- c(50,800)
  img <- function(profile){
    require(ggplot2)
    df <- data.frame(time=rep(rt,each=nrow(profile)),
                     mz=rep(mz_range[1]:mz_range[2],ncol(profile)),
                     int=as.numeric(profile))
    df <- df[profile>missing,]
    p <- ggplot(df,aes(x=time,y=mz,fill=int))+geom_tile()
    p <- p+scale_fill_gradient(low='yellow',high='black',trans='log')
    p <- p+opts(legend.position='none')+scale_x_continuous(limits=c(0,3700))
    print(p)
  }

img(prof)


## explore it,slow(1 min),without legend,please change to qtpaint
png("~/Desktop/test2.png",800.800)
image(raw_prof2)
dev.off()
sum(prof2==min(prof2))
## baseline subtraction using sliding window median smoother
cor_prof <- removeBaseline(raw_prof, "median", scanrad = 100)

## explore the correction,add argument mz
explore(cor_prof, raw = raw_prof,mz=150)

## accept correction
raw <- cor_prof

test <- profMz(object)

## find some peaks using our island-fitting method
#peaks <- findPeaks(raw, "gauss")
peaks <- findPeaks(object=raw, "gauss")
pk <- peaks@.Data
summary(pk[,'maxf']/pk[,'sigma'])
prof <- raw@env$profile
sum(apply(prof,2,sum)>0)

save(raw,raw_prof,peaks, file="~/peaks.rda")
#load("peaks.rda")
head(peaks@.Data)
## explore them
explore(peaks, raw = raw, sample=NA,residuals=TRUE,island=TRUE)
#use island=F, give quick response for peak analyse



## read in another file (a replicate) and generate its profile matrix
rep_file <- "raw/s1/M0102A.CDF"
rep_file <- "/home/tengfei/Desktop/sample/s1/a2.CDF"
rep_raw <- perform(pipeline(raw), rep_file)

## explore it, make sure it's OK like the other rep
#explore(rep_raw)

head(peaks@.Data)

## find its peaks in the same way
rep_peaks <- perform(findPeaksProto(peaks), rep_raw)

## explore them
explore(rep_peaks, raw = rep_raw,residuals=T,island=F)

## load all CI/Control reps at once into one experiment
s1_exp <- loadExperiment("/home/tengfei/Desktop/sample/s1/", pipeline = pipeline(peaks))
#save(s1_exp, file="s1_exp.rda")
#load("ci_exp.rda")
samples <- c("raw/s1", "raw/s2")
samples <- dir_design
total_peaks <-  perform(pipeline(s1_exp), samples)


# explore it
explore(s1_exp)## need to be fixed

# find components using our 'sigma' method
xset_comps <- findComps(s1_exp, "sigma")
xset_comps_filt <- findComps(s1_exp,"sigma_filt",tic.cutoff=0,npeaks.cutoff=0)
# explore
save(s1_exp,xset_comps,xset_comps_filt,file="xset_comps.rda")

explore(xset_comps_filt,sample=1,residual=TRUE,island=FALSE)
# load peak set for both CI and SD plants and find components

samples <- c("raw/s1", "raw/s2")

s1_s2_xset <- perform(pipeline(xset_comps), samples)

s1_s2_xset_filt<- perform(pipeline(xset_comps_filt),samples)
#save(s1_s2_xset, file="s1_s2_xset.rda")
#load("ci_ld_xset.rda")
save(s1_s2_xset,s1_s2_xset_filt,file="./s1_s2.rda")


explore(s1_s2_xset,sample=2,residual=T,island=F)

# find groups across samples using spectral angle distance
xset_groups <- groupComps(s1_s2_xset, "angle")
xset_groups_filt<- groupComps(s1_s2_xset_filt,"angle")

# explore the groups
explore(xset_groups)
explore(xset_groups_filt)

# correct retention time using robust loess smoother

xset_rtcor <- rtcor(xset_groups, "rloess")
xset_rtcor_filt<- rtcor(xset_groups_filt,"rloess")
# explore correction
explore(xset_rtcor_filt,raw=xset_groups_filt,geom="heatmap")
explore(xset_rtcor,raw=xset_groups,sample=1:2,xscale=c(1000,1500))

# regroup after rt correction
xset_groups2 <- groupComps(xset_rtcor, "angle", rt_window = 5)
xset_groups2_filt<- groupComps(xset_rtcor_filt, "angle", rt_window = 5)

xset_rtcor2 <- rtcor(xset_groups2, "rloess")
xset_rtcor_filt2<- rtcor(xset_groups2_filt,"rloess")

explore(xset_rtcor_filt2,raw=xset_groups2_filt,geom="heatmap")

xset_groups3 <- groupComps(xset_rtcor, "angle", rt_window = 1)
xset_groups_filt3<- groupComps(xset_rtcor_filt, "angle", rt_window = 1)
xset_rtcor3 <- rtcor(xset_groups3, "rloess")
xset_rtcor_filt3<- rtcor(xset_groups_filt3,"rloess")

explore(xset_rtcor3,raw=xset_groups3,xscale=c(2300,2500))

# directly filt by rt_window and dist.cutoff
xset_groups3 <- groupComps(s1_s2_xset, "angle",rt_window=5)
xset_groups_filt3<- groupComps(s1_s2_xset_filt,"angle",rt_window=5)

save(xset_groups,xset_groups_filt,xset_groups2,xset_groups2_filt,
     xset_groups3,xset_groups_filt3,xset_rtcor,xset_rtcor_filt,file="./group.rda")


# summarize by common peaks
xset_sum2 <- summarize(xset_groups2, "common")
xset_sum1 <- summarize(xset_groups,"common")
xset_sum3 <- summarize(xset_groups3,"common")
xset_sum2_filt<- summarize(xset_groups2_filt, "common")
xset_sum1_filt<- summarize(xset_groups_filt,"common")
xset_sum3_filt<- summarize(xset_groups_filt3,"common")
# explore our results
explore(xset_sum1)

# possibly normalize (by simple scaling)
xset_norm1 <- normalize(xset_sum1, "scale")
xset_norm2 <- normalize(xset_sum2, "scale")
xset_norm3 <- normalize(xset_sum2, "scale")
xset_norm1_filt<- normalize(xset_sum1_filt, "scale")
xset_norm2_filt<- normalize(xset_sum2_filt, "scale")
xset_norm3_filt<- normalize(xset_sum3_filt, "scale")

save(xset_sum1,xset_sum2,xset_sum3,xset_sum1_filt,xset_sum2_filt,xset_sum3_filt,
     xset_norm1,xset_norm2,xset_norm3,xset_norm1_filt,xset_norm2_filt,xset_norm3_filt,
     file="./sum_norm.rda")

explore(xset_norm3)

# again, explore
#explore(xset_norm)

# extract expressionSet and pass to explorase
#library(explorase)
#explorase(as(xset_norm, "ExpressionSet"))
#out put
dump2msp(xset_norm1_filt,id=c(183,80),file="~/interest.msp")
library(explorase)
explorase(as(xset_norm1_filt,"ExpressionSet"))
str(test)

datapath <- "~/Prolang/R/data/"
setwd(datapath)
load("total_peaks.rda")

gp <- d@groups
head(gp)
dim(gp)
nrow(gp[is.na(gp[,'sigma_d']),])

######################################################################
dt1 <- c(17.92,22.79,24.69,12.32,11.49,17.28,27,40.5)
dt2 <- c(17.67,6.99,7.72,9.63,40.53,40.87,40.34,25.62,21.83,12.54,39.1)
df <- c(dt1,dt2)
length(dt2)
length(df)
seq(from=0.5,to=0.9,by=0.05)
load('~/Data/tengfei/rda/s1s2xset.rda')

for(i in seq(0.01,0.14,by=0.01)){
   s1_s2_xset2<- compsFilter(s1_s2_xset,tic.cutoff=0,npeaks.cutoff=0,tic=0,npeaks=10)
   xset_groups <- groupComps(s1_s2_xset2,"angle",rt_window=1,dist.cutoff=i)
   xset_sum <- summarize(xset_groups, "common")
   xset_norm <- normalize(xset_sum, "scale")
   tdiff <- topDiff(xset_norm,pvalue=0.0001)
   rt <- tdiff[,'rt']
   cat(rt,'\n')
   percent<- rtmatch(df,rt,range.min=-1/60,range.max=1/60,plot.match=F)
   cat('new rt number:',length(rt),' dist_cutoff:',i,' percent:',percent,'\n',file="~/Desktop/para_dist",append=T)
 }
   

for(i in seq(from=0.01,to=0.20,by=0.005)){
  for(j in seq(from=0.50,to=0.85,by=0.05)){
    for(k in seq(from=0.50,to=0.85,by=0.05)){
    cat(i,j,'\n')
    d <- finalObj(total_peaks,tic.cutoff=j,npeaks.cutoff=k,
        dist.cutoff=i,rt_window=1)
    ob <- topDiff(d,pvalue=0.001)
    rt <- ob[,'rt']
    cat(rt,'\n')
    percent<- rtmatch(df,rt,range.min=-1/60,range.max=1/60,plot.match=F)
    cat(i,j,k,percent,'\n',file="~/Desktop/para2",append=T)
 }}}


for(j in seq(from=0.005,to=0.5,by=0.005)){
    d <- finalObj(total_peaks,tic.cutoff=0,npeaks.cutoff=0,
        dist.cutoff=j,rt_window=1)
    ob <- topDiff(d,pvalue=0.001)
    rt <- ob[,'rt']
    percent<- rtmatch(df,rt,range.min=-0.03,range.max=0.03,plot.match=F)
    cat(j,percent,'\n',file="~/Desktop/para_distonly2",append=T)
 }

for(j in seq(from=0,to=1,by=0.02)){
    d <- finalObj(total_peaks,tic.cutoff=0,npeaks.cutoff=j,
        dist.cutoff=0.05,rt_window=1)
    ob <- topDiff(d,pvalue=0.001)
    rt <- ob[,'rt']
    percent<- rtmatch(df,rt,range.min=-0.03,range.max=0.03,plot.match=F)
    cat(j,percent,'\n',file="~/Desktop/para_dist_npeaks2",append=T)
 }

for(j in seq(from=0,to=1,by=0.02)){
    d <- finalObj(total_peaks,tic.cutoff=j,npeaks.cutoff=0,
        dist.cutoff=0.05,rt_window=1)
    ob <- topDiff(d,pvalue=0.001)
    rt <- ob[,'rt']
    percent<- rtmatch(df,rt,range.min=-0.03,range.max=0.03,plot.match=F)
    cat(j,percent,'\n',file="~/Desktop/para_dist_tic2",append=T)
 }

setwd("~/Desktop")
png("distonly.png",800,800)
z <- read.table("~/Desktop/para_distonly")
qplot(V1,V2,data=z,xlab="Cutoff Value",ylab="Match(%)" ,main="Distance Cutoff Test",geom="line")
dev.off()
png("dist0.05_npeaks.png",800,800)
z <- read.table("~/Desktop/para_dist_npeaks")
qplot(V1,V2,data=z,xlab="Cutoff Value",ylab="Match(%)" ,main="Peaks Number(npeaks) Cutoff Test",geom="line")
dev.off()
png("dist0.05_tic.png",800,800)
z <- read.table("~/Desktop/para_dist_tic")
qplot(V1,V2,data=z,xlab="Cutoff Value",ylab="Match(%)" ,main="Tic Cutoff Test",geom="line")
dev.off()

png("dist0.05_tic2.png",800,800)
z <- read.table("~/Desktop/para_dist_tic2")
qplot(V1,V2,data=z,xlab="Cutoff Value",ylab="Match(%)" ,main="Tic Cutoff Test",geom="line")
dev.off()

png("dist0.05_npeaks2.png",800,800)
z <- read.table("~/Desktop/para_dist_npeaks2")
qplot(V1,V2,data=z,xlab="Cutoff Value",ylab="Match(%)" ,main="Npeaks Cutoff Test",geom="line")
dev.off()

png("distonly2.png",800,800)
z <- read.table("~/Desktop/para_distonly2")
qplot(V1,V2,data=z,xlab="Cutoff Value",ylab="Match(%)" ,main="Distance Cutoff Test",geom="line")
dev.off()



head(z)
plot(z,type='l')
z[z$V2>8,] #  0.64 <npeaks<0.72     0.72<tic<0.74


d <- finalObj(total_peaks,tic.cutoff=0.70,npeaks.cutoff=0,
        dist.cutoff=0.05,rt_window=1)
d4 <- finalObj(total_peaks,tic.cutoff=0,npeaks.cutoff=0,
        dist.cutoff=0.05,rt_window=1)

d1 <- finalObj(total_peaks,tic.cutoff=0.5,npeaks.cutoff=0,
        dist.cutoff=0.05,rt_window=1)

d2 <- finalObj(total_peaks,tic.cutoff=0,npeaks.cutoff=0,
        dist.cutoff=0.05,rt_window=1)

d3 <- finalObj(total_peaks,tic.cutoff=0.7,npeaks.cutoff=0,
        dist.cutoff=0.05,rt_window=2)



d5 <- finalObj(total_peaks,tic.cutoff=0.64,npeaks.cutoff=0,
        dist.cutoff=0.05,rt_window=2)

cplotResult(d,id=id)
cplotResult

head(ob)
ob[ob[,'ncomps']==8&ob$pvalue<0.1,]
id <- rownames(ob)
dim(ob1)
ob<- topDiff(d,pvalue=0.001)
ob1<- topDiff(d1,pvalue=0.001)
ob4 <- topDiff(d4,pvalue=0.001)
ob2 <- topDiff(d2,pvalue=0.001)
ob3 <- topDiff(d3,pvalue=1)
dim(ob)
dim(ob3)

exprs(as(d,"ExpressionSet"))
ob <- ob1
dim(ob)
ob[order(ob$rt),-c(4,5,6)]
 rt <- ob[,'rt']
rtmatch(df,rt,range.min=-0.03,range.max=0.03,plot.match=F)
head(ob)
ob
1/60
cplotResult(d,208)
ob[order(ob$rt),-c(4,5,6,8)]
ob2[order(ob2$rt),-c(4,5,6,8)]

cplotSpec(d2,group=245)
cplotSpec(d,group=208)




cat(0.1,0.2,file='~/Desktop/para')   
exprs(as(d,"ExpressionSet"))

x <- read.table('~/Desktop/para')
x[order(x[,4],decreasing=T),]

cplotDt(d)
head(gp)
summary(gp)
gp <- data.frame(gp)
qplot(ncomps,npeaks_mean,data=gp)
library(explorase)
explorase(as(d,"ExpressionSet"))

file.path(dir,"nl")
cplotSpec(d,group=94)
setwd("~/Desktop/data/quantity")
length(list.files(dir))
dir <- "~/Desktop/data/dist0.02"
cplotPvalue(dir,xlim=c(1,50))

pa <- seq(0,1,by=0.02)

for(i in pa){
  name <- paste("dist",i,".csv",sep="")
  d <-  finalObj(total_peaks,tic.cutoff=i,npeaks.cutoff=0,
        dist.cutoff=0.02,rt_window=2)
  d <- getPvalue(d)
  write.csv(d,file=name,row.names=F)
}


  d <- finalObj(total_peaks,tic.cutoff=0,npeaks.cutoff=0,
        dist.cutoff=0.02,rt_window=2)

  
 
topDiff(d17,1:40)
topCommon(d17,1:40,method="sum")
cplotSpec(d17,group=289) #15.07 min
cplotSpec(d17,group=147) #12.17 min
cplotResult(d,id=c(12,230,337,248))
head(getFeatures(d17))


test
id <- topDiff(d6,1:10)
id
id <- rownames(id)
diff <- topDiff(d17,1:40)
diff[order(diff$rt),c(1,7,9)]
cplotResult(d17,c(337,12,230,248))
cplotSpec(d13,group=13)
nrow(diff[diff$pvalue<0.001,])
test <- topCommon(d10)
names(test)
test[order(test$rt),c(1,7,10,11)]

cplotResult(d10,id=238)
cplotResult(d6,id=c(137,412,197,395))
cplotSpec(d9,group=121)
library(explorase)
explorase(as(d6,"ExpressionSet"))

cplotSpec(d2,group=248)

 df <- getPvalue(d)

test <- read.csv("~/Desktop/data/dist0.02.csv",header=T)
test <- df
head(test[order(test$pvalue),],n=16)

d@groups[order(d@groups[,'sigma_t'],decreasing=T),]
cplotSpec(d,group=166)
group=414
d@groups
length(inter)



library(chromatoplots)
source("~/Prolang/svn/repos/trunk/chromatoplots/R/cplot.R")
source("~/Prolang/svn/repos/trunk/chromatoplots/R/cpExperiment-class.R")               
######################################################################
## pvalue test for parameter ##
######################################################################
dir <- "~/Desktop/data/npeaks"
dir <- "~/Desktop/data/distonly"
dir <- "~/Desktop/data/dist_rt2"
dir <- "~/Desktop/data/dist_rt5"
dir <- "~/Desktop/data/dist_rt10"
dir <- "~/Desktop/data/quantity"
dir <- "~/Desktop/data/dist_0_0.01"
dir <- "~/Desktop/data/dist0.5tic0.5"
dir <- "~/Desktop/data/dist0.02"
dir <- "~/Desktop/data/tic0.69npeak0.7"
dir <- "~/Desktop/data/tic0.69npeaks0.7dist0.01_0.02"

setwd("~/Desktop")

png("dist2.png",800,800)
cplotTotalP(dir,method="sum")
dev.off()
png("dist2_diff.png",800,800)
cplotTotalP(dir,method="diff")
dev.off()
png("dist2_sum_50.png",800,800)
cplotTotalP(dir,method="sum",xscale=c(0,20))
dev.off()
png("dist2_diff_50.png",800,800)
cplotTotalP(dir,method="diff",xscale=c(0,20))
dev.off()
png("dist2_hist.png",800,800)
cplotHist(dir)
dev.off()


cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/dist_rt2")

cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/distonly")
cplotPvalue(dir,xlim=c(1,50),output="~/Prolang/R/chrom/figures/para/distonly_50")
cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/quantity")
cplotPvalue(dir,xlim=c(1,50),output="~/Prolang/R/chrom/figures/para/quantity_50")
cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/npeaks")
cplotPvalue(dir,xlim =c(1,50),output="~/Prolang/R/chrom/figures/para/npeaks_50")
cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/dist_rt5")
cplotPvalue(dir,xlim=c(1,50),output="~/Prolang/R/chrom/figures/para/dist_rt5_50")
cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/dist_rt10")
cplotPvalue(dir,xlim=c(1,50),output="~/Prolang/R/chrom/figures/para/dist_rt10_50")
cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/dist_0_0.01")
cplotPvalue(dir,xlim=c(1,50),output="~/Prolang/R/chrom/figures/para/dist50_0_0.01")
cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/dist0.5tic0.5")
cplotPvalue(dir,xlim=c(1,50),output="~/Prolang/R/chrom/figures/para/dist0.5tic0.5_50")
cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/dist0.02tic0.5")
cplotPvalue(dir,xlim=c(1,50),output="~/Prolang/R/chrom/figures/para/dist0.02tic0.5_50")
cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/tic0.69npeaks0.7")
cplotPvalue(dir,xlim=c(1,50),output="~/Prolang/R/chrom/figures/para/tic0.69npeaks0.7_50")
cplotPvalue(dir,output="~/Prolang/R/chrom/figures/para/dist0.02tic0.69npeaks0.7")
cplotPvalue(dir,xlim=c(1,50),output="~/Prolang/R/chrom/figures/para/dist0.02tic0.69npeaks0.7_50")


cplotHist(dir)

d <- data.frame()
t=c(1:5,5:8)
d <- data.frame(t=t)
d
max(table(t))
table(t)==2
names(table(t))[table(t)==as.numeric(max(table(t)))]

as.string(3)
library(explorase)
explorase(as(d,"ExpressionSet"))


z <- read.table("~/Desktop/para2")
z[order(z$V4,decreasing=T),][1:100,]
best <- z[z$V4>89,]
best[best$V2==0.8&best$V3==0.75,]
dim(best)

table(best$V1)
table(best$V2)
table(best$V3)
hist(best$V3)

######################################################################
## statistcis about genuine data
######################################################################
setwd('~/Data/tengfei/rda')
ff <- load('ff.rda')
ss <- topDiff(final,pvalue=0.000001)
al <- topDiff(final,pvalue=1)
pk <- final@peaks

id <- rownames(ss)

cplotSpec(final,group=552)

peaks <- final@peaks
comps <- final@comps

major <- al[rownames(al)%in%id,]
minor <- al[!(rownames(al)%in%id),]


hist(major[,'npeaks_mean'])
summary(major[,'npeaks_mean'])
hist(minor[,'npeaks_mean'])
summary(minor[,'npeaks_mean'])

df_mean <- data.frame(npeaks=c(major[,'npeaks_mean'],minor[,'npeaks_mean']),
                      group=rep(c('major','noise'),c(nrow(major),nrow(minor))))

png('denpeak.png',800,800)
p <- ggplot(df_mean,aes(x=npeaks,..density..,color=group))
p+geom_freqpoly()+facet_grid(group~.)+theme_bw()
dev.off()

png('denpeak.png',800,800)
p <- ggplot(df_mean,aes(x=npeaks,..density..,color=group))
p+geom_freqpoly()+facet_grid(group~.)+theme_bw()+scale_y_continuous(limits=c(0,0.3))
dev.off()

png('boxden.png',800,800)
p <- ggplot(df_mean,aes(factor(group),npeaks))
p+geom_boxplot()+theme_bw()
dev.off()

qq <- input_quantity(final)
cp <- comps
head(cp)
cp <- cp[,-8]

cp.mj <- cp[cp[,'group']%in%id,]
cp.mr <- cp[!cp[,'group']%in%id,]
cp.mj <- as.data.frame(cp.mj)
cp.mr <- as.data.frame(cp.mr)
inter <- interaction(cp.mj[,'sample'],cp.mj[,'comp'])
inter <- as.character(inter)
cp.mj <- cbind(cp.mj,inter=inter)
inter <- interaction(cp.mr[,'sample'],cp.mr[,'comp'])
inter <- as.character(inter)
cp.mr <- cbind(cp.mr,inter=inter)

inter.mj <- cp.mj[,'inter']
inter.mr <- cp.mr[,'inter']


plot(density(cp.mj[,'quantity']))
cp.mj.int <- cp.mj[,'quantity']
cp.mr.int <- cp.mr[,'quantity']

peak.mj <- peaks[as.character(interaction(peaks[,'sample'],peaks[,'comp']))%in%inter.mj,]
peak.mr <- peaks[as.character(interaction(peaks[,'sample'],peaks[,'comp']))%in%inter.mr,]

pk <- final@peaks
head(pk)
hist(peak.mj[,'rtmax']-peak.mj[,'rtmin'])
hist(peak.mj[,'span'])
plot(pk[,'maxf'],pk[,'sigma'])

png('~/Desktop/sigma.png',800,800)
hist(peak.mj[,'sigma'],main="Sigma for major components",xlab='Sigma')
dev.off()

peak.mj <- as.data.frame(peak.mj)
peak.mr <- as.data.frame(peak.mr)

df <- rbind(peak.mj,peak.mr)
df <- cbind(df,group=rep(c('major','noise'),c(nrow(peak.mj),nrow(peak.mr))))
library(ggplot2)
head(df)

png('denforall2.png',800,800)
p <- ggplot(df2,aes(x=maxf,..density..,color=group))
p+geom_freqpoly()+facet_grid(group~.)+xlab('Intensity/maxf')+theme_bw()
dev.off()

png('denforall2.png',800,800)
p <- ggplot(df2,aes(x=maxf,..density..,color=group))
p+geom_freqpoly()+facet_grid(group~.)+xlab('Intensity/maxf')+theme_bw()+scale_y_continuous(limits=c(0,5.0e-09))
dev.off()

p <- ggplot(df2,aes(x=maxf,color=group))
p+geom_histogram()+facet_grid(group~.)+xlab('Intensity/maxf')+theme_bw()+scale_y_continuous(limits=c(0,2.0e-07))
p <- ggplot(df,aes(y=maxf,x=group))
p+geom_boxplot(aes(group=group))

1.0e-02

p <- ggplot(df,aes(x=tau,color=group))
p+geom_density()+facet_grid(group~.)

p <- ggplot(df,aes(x=span,..density..,color=group))
p+geom_freqpoly()+facet_grid(group~.)

summary(peak.mj[,'tau'])
summary(peak.mr[,'tau'])

summary(final@peaks[,'tau'])

head(peak.mj)
plot(density(peak.mr[,'maxf']),type='l',col='red',xlim=c(0,100000))
lines(density(peak.mj[,'maxf']),type='l')


summary(peak.mj[,'maxf'])
peak.mr[,'maxf']

head(peaks)

prof <- raw_prof@env$profile
prof <- as.numeric(prof)
pf <- prof[prof>150]

df2 <- df[,c(10,19)]
head(df2)
df2 <- rbind(df2,data.frame(maxf=pf,group=rep('profile matrix',length(pf))))


prof2 <- prof[prof>1000000]
length(prof)
length(prof2)
plot(density(prof[prof>150]))
plot(density(prof2[prof2<3000]))
hist(prof[prof<10000])

getwd()
png('denforprof')
plot(density(prof[prof>150]),main='Density for intensity(over 150)',xlab='Intensity',ylab='Density')
dev.off()
hist(prof2)

summary(prof2)

length(prof)
length(prof[prof==150])


x <- 1:8e06

plot(x,dexp(x,0.005))

######################################################################
## check the SD of different levels, rt
######################################################################

file <- '/home/tengfei/Data/tengfei/cal/raw/s1/M0101A.CDF'
raw <- loadSample(file)
raw_prof <- genProfile(raw)
sum(raw_prof@env$profile>150)
cor_prof <- removeBaseline(raw_prof,'median',scanrad=100)

png('~/Desktop/raw.png',800,800)
image(raw_prof)
dev.off()


png('~/Desktop/raw_cor.png',800,800)
image(cor_prof)
dev.off()

peaks.filt <- peaks.filt[peaks.filt[,'sigma']<2.5,]
boxplot(peaks.filt[,'sigma/h'])
pk.f <- peaks.filt
head(pk.f[order(pk.f[,'sigma/h'],decreasing=T),])
plot(peaks.filt[,'sigma/h'],peaks.filt[,'maxf'])

raw1 <- cor_prof
raw <- cor_prof
peaks <- findPeaks(object=raw1,'gauss')

explore(peaks, raw = raw1, sample=NA,residuals=TRUE,island=TRUE)
save(peaks,raw1,raw,raw_prof,file='~/Data/tengfei/rda/peaks214.rda')

load('~/Data/tengfei/rda/peaks214.rda')

plot(abs(peaks[,'delta_t']),peaks[,'sigma'])
idx <- !(peaks[,'rt']<peaks[,'rtmin']|peaks[,'rt']>peaks[,'rtmax'])
peaks.filt <- peaks[idx,]



## focus on big sigma maybe flat??
peaks.bad <- peaks.filt[peaks.filt[,'sigma']>3,]
peaks.filt2 <- peaks.filt[peaks.filt[,'sigma']<3,]
boxplot(peaks.bad[,'sigma/h'])
plot(x=peaks.bad[,'sigma'],y=peaks.bad[,'maxf'])

pdf('~/Desktop/sigmaraw.pdf',width=8,height=4)
par(mfrow=c(1,2))
hist(peaks.filt2[,'sigma'],xlab=expression(sigma),main='Real data')
plot(dd,dnorm(dd,1.3,0.3),type='l',xlab=expression(sigma),ylab=density,main='Simulated data')
dev.off()

pks <- peaks.filt2[abs(peaks.filt2[,'tau'])<3,]
dd <- seq(-3,3,length=100)

pdf('~/Desktop/tauraw.pdf',width=8,height=4)
par(mfrow=c(1,2))
hist(pks[,'tau'],xlab=expression(tau),main='Real Data')
plot(dd,dnorm(dd,0,0.5),xlab=expression(tau),main='Simulated Data',type='l')
dev.off()


head(peaks.filt2)
sum(peaks.filt2[,'tau']>1)

plot(peaks.filt2[,'tau'],peaks.filt2[,'maxf'])

sum(abs(peaks.filt2[,'tau'])>3)

dd <- seq(0,3,length=100)
plot(dd,dnorm(dd,1.3,0.3),type='l')

head(peaks.bad)

boxplot(peaks.filt2[,'maxf'])

pk.filt2 <- as.data.frame(peaks.filt2)

qplot(data=pk.filt,x=maxf,y=sigma)
## test 10 15

class(cor_prof)

pk <- peaks@.Data
pk <- as.data.frame(pk)
pk.filt <- as.data.frame(peaks.filt)
pk.filt$class[pk.filt$sigma>3] <- 'sigma>3'
pk.filt$class[pk.filt$sigma<=3] <- 'sigma<=3'

qplot(data=pk.filt,x=class,y=maxf,geom='boxplot')
qplot(data=pk.filt,x=maxf,y=..density..,facets=class~.,geom='histogram')

library(ggplot2)

qplot(data=pk,x=maxf,y=error)
qplot(data=pk,x=maxf,y=delta_t)
qplot(data=pk,x=maxf,y=delta_i)
#qplot(data=pk,x=abs(delta_i),y=abs(delta_t))
qplot(data=pk,x=abs(delta_t),y=abs(tau))
qplot(y=delta_i,data=pk,geom="boxplot")
boxplot(pk$delta_i)
boxplot(pk$delta_t)

prof1 <- raw_prof@env$profile
prof2 <- cor_prof@env$profile
prof3 <- prof2
prof3[prof3<0]=0


pk[pk$mz==80,]
head(pk)

plot_peak(peaks,761,raw=raw1,cutoff=matrix(q2, nrow(prof2), ncol(prof2)),island=F)

head(pk)

cplotPeaks2(peaks,raw=raw1,mz=80,rt_range=c(433,440))
cplotPeaks2(peaks,raw=raw1,mz=81)
str(peaks)


q1 <- quantile(prof1,0.99,na.rm=T)
q2 <- quantile(prof2,0.99,na.rm=T)
q3 <- quantile(prof3,0.99,na.rm=T)

q1
q2
q3

sd(prof1[prof1 < q1])
sd(prof2[prof2 < q2])
sd(prof3[prof3 < q3])

sd(prof1[prof1 > q1])
sd(prof2[prof2 > q2])
sd(prof3[prof3 > q3])

lst1 <- apply(prof1,1,function(x){
  c(sd(x[x<q1]),sd(x[x>q1]))
})

lst2 <- apply(prof2,1,function(x){
  quantile(x,0.99,na.rm=T)
})

lst3 <- apply(prof3,1,function(x){
  quantile(x,0.99,na.rm=T)
})

boxplot(lst3)
plot(1:length(lst3),lst3)
summary(lst3)

lst1 <- (t(lst1))
plot(1:nrow(lst1),lst1[,1])

sd(prof1[prof1 > q1])
sd(prof2[prof2 > q2])
sd(prof3[prof3 > q3])



## try env
myfun <- function(x){
  env <- new.env()
  a=x+1
  b=a+1
  as.list(env)
}


file <- '~/Data/tengfei/rda/ff2.rda'
load(file)
tt <- final2@peaks[,'delta_t']
length(tt[tt>1])

pks <- final2@peaks
pks[pks[,'delta_t']<10&pks[,'delta_t']>5,]
idx <- pks[,'rt']>pks[,'rtmax']|pks[,'rt']<pks[,'rtmin']
pks2 <- pks[idx,]
pks2[pks2[,'delta_t']<10&pks2[,'delta_t']>5,]

head(pks2)




plot_peak(final2,13,raw=raw1,cutoff=matrix(q2, nrow(prof2), ncol(prof2)),island=F)
cplotPeaks2(final2,raw=raw1,mz=50,rt_range=c(210,220))

plot_peak(pks,3744,raw=raw1,cutoff=matrix(q2, nrow(prof2), ncol(prof2)),island=F)
cplotPeaks2(final2,raw=raw1,mz=84,rt_range=c(1032,1040))

load('~/Desktop/finalnoise4.rda')
tdiff <- topDiff(final,pvalue=0.0000001)
dim(tdiff)

load('~/Desktop/peakmove.rda')
tdiff2 <- topDiff(final,pvalue=0.0000001)
dim(tdiff2)

