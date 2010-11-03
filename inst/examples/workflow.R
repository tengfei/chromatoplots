library(commandr)
library(chromatoplots)
cdfpath <- system.file('cdf',package='faahKO')
cdffiles <- list.files(cdfpath,recursive=TRUE,full.names=TRUE)

getOption('BioC')$commandr

## options(warn=0)
raw1 <- loadSample(cdffiles[1])
## rawMat(raw1)
## FIXME: too long, maybe export to a graphic file
## explore(raw1)
raw_prof <- genProfile(raw1)
## explore(raw_prof)
## baseline subtraction using sliding window median smoother
cor_prof <- removeBaseline(raw_prof, "median", scanrad = 100)
##  Slow and need to be fixed
## cor_prof2 <- removeBaseline(raw_prof, "rbe")
## other baseline removement method

## getOption('BioC')
## explore the correction

## explore(cor_prof, raw = raw_prof,geom='l')
## ## default is 'l', because line is faster
## cplot(cor_prof,raw=raw_prof,mz=201,geom='p')

## accept correction
raw <- cor_prof
## find some peaks using our gaussian-fitting method
## peaks <- findPeaks(raw, "gauss")
## find some peaks using our parabola-fitting method
## peaks <- findPeaks(raw, "parabola")

## xcms has method like
## $findPeaks.methods
## [1] "centWave"      "matchedFilter" "MS1"           "MSW"          
## [5] "pipeline"      "islands"      
peaks <- findPeaks(raw, "gauss")        #slow
pipeline(peaks)
## peaks2 <- findPeaks(raw, "parabola")
## got a bug for xcms method
## peaks3 <- findPeaks(raw, "centWave")
## options(error=recover)
## save(peaks, file="peaks.rda")
## load("peaks.rda")
## explore them
## explore(peaks, raw = raw)

getOption('BioC')$commandr$stageData

##options(error=recover)
rep_file <- cdffiles[[2]]
rep_raw <- perform(pipeline(raw), rep_file)
pipes0 <- Pipeline(Protocol('loadSample','xcms'),
                  Protocol('genProfile','intbin'))
## it's the same with
pipes1 <- Pipeline(Protocol('loadSample'),
                  Protocol('genProfile'))
identical(pipes0,pipes1)
rep_raw0 <- perform(pipes1,rep_file)

## explore it, make sure it's OK like the other rep
## explore(rep_raw,raw=rep_raw0)

## find its peaks in the same way
rep_peaks <- perform(findPeaksProto(peaks), rep_raw)

## FIXME: should automatically change it.
## explore them
## explore(rep_peaks, raw = rep_raw)

## load all CI/Control reps at once into one experiment
##ci_exp <- loadExperiment("CI/Control", pipeline = pipeline(peaks))
## or LD Control and Const
setwd(cdfpath)
## ci_exp <- loadExperiment(c("WT", "KO"),pipeline = pipeline(peaks))
## a shorter form, default is cononical one
ci_exp <- loadExperiment(c("WT", "KO"))
save(ci_exp, file="~/Data/rda/ci_exp.rda")
load("~/Data/rda/ci_exp.rda")

# explore it
##explore(ci_exp)
# find components using our 'sigma' method
xset_comps <- findComps(ci_exp, "sigma")

explore(xset_comps,sample=1)

# load peak set for both CI and SD plants and find components
samples <- c("WT", "KO")

ci_ld_xset <- perform(pipeline(xset_comps), samples)
## save(ci_ld_xset, file="ci_ld_xset.rda")
## load("ci_ld_xset.rda")

## # explore
## explore(ci_ld_xset)

# find groups across samples using spectral angle distance
xset_groups <- groupComps(xset_comps, "angle")
perform(pipeline(xset_comps),c('WT','KO'))

# explore the groups
explore(xset_groups)

# correct retention time using robust loess smoother
xset_rtcor <- rtcor(xset_groups, "rloess")
save(xset_groups,xset_rtcor, file="~/Data/rda/xset_rtcor.rda")
bioc <- getOption('BioC')
save(bioc,file="~/Data/rda/bioc.rda")

load("~/Data/rda/bioc.rda")
load("~/Data/rda/rtcomp.rda")
options('BioC'=bioc)

# explore correction
explore(xset_rtcor)



# correct retention time using robust loess smoother
xset_rtcor <- rtcor(xset_groups, "rloess")
str(xset_rtcor@pipeline@.Data[[2]])
# explore correction
explore(xset_rtcor)

# regroup after rt correction
xset_groups2 <- groupComps(xset_rtcor, "angle", rt_window = 1)

# summarize by common peaks
xset_sum <- summarize(xset_groups2, "common")

# explore our results
explore(xset_sum)

# possibly normalize (by simple scaling)
xset_norm <- normalize(xset_sum, "scale")

# again, explore
explore(xset_norm)

## diagnostic without explorase

# extract expressionSet and pass to explorase
library(explorase)
explorase(as(xset_norm, "ExpressionSet"))

## metablomics identification



