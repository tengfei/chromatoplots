######################################################################
## outpout for downstream analysis ##
######################################################################

######################################################################
## First to make an output from chromatoplots to MSP file ,which can be
## read by NIST library
######################################################################
dump2msp <- function(object,id=NULL,file=NULL){
 if(is.null(id)|is.null(file))
  stop("Please specify the metabolite of interest and file(with file path) you want to write in.")
comps <- object@comps
peaks <- object@peaks
groups <- object@groups
lst <- by(comps,comps[,'group'],function(x){
  np<- sort(x[,'npeaks'])
  np.main <- np[round(median(1:length(np)))]
  comp.main <- x[x[,'npeaks']==np.main,'comp'][1]
  sample.main <-x[x[,'npeaks']==np.main,'sample'][1]
  df <- data.frame(comp=comp.main,
                   sample=sample.main,
                   group=rep(x[1,'group'],length(comp.main)))
   return (df)
})
df <- do.call(rbind,lst)
df.sub <- subset(df,df[,'group']%in%id)
for(i in 1:length(id)){
sample <- df.sub[df.sub[,'group']==id[i],'sample']
comp <- df.sub[df.sub[,'group']==id[i],'comp']
t.sub <- subset(peaks,peaks[,'sample']==sample&peaks[,'comp']==comp)
rt <- comps[comps[,'group']==id[i]&
            comps[,'sample']==sample&
            comps[,'comp']==comp,'rt']/60
npeaks <- comps[comps[,'group']==id[i]&
            comps[,'sample']==sample&
            comps[,'comp']==comp,'npeaks']
cat("NAME:ID=",id[i],"\r\n",file=file,append=T)
cat("COMMENT:\r\n",file=file,append=T)
cat("RI:\r\n",file=file,append=T)
cat("CASNO:\r\n",file=file,append=T)
cat("RT:",rt,"\r\n",file=file,append=T)
cat("SOURCE:\r\n",file=file,append=T)
cat("NUM PEAKS:",npeaks,"\r\n",file=file,append=T)
n=0
for(i in 1:nrow(t.sub)){
 # cat("(",t.sub[i,'mz']," ",t.sub[i,'maxf'],")",file="./test.txt",sep="",append=T)
 text <- sprintf("(%3d%10d)",round(t.sub[i,'mz']),round(t.sub[i,'maxf']))
 cat(text,file=file,sep=" ",append=T)
  n=n+1
  if(n==5&i!=nrow(t.sub)){
    cat("\r\n",file=file,sep="",append=T)
    n=0
  }
}
cat("\r\n\r\n",file=file,append=T)
}}


