\name{explore}
\alias{explore}
\alias{explore,cpPeaks,ProtoFindPeaks-method}
\alias{explore,cpPeaks,ProtoFindPeaksParabola-method}
\alias{explore,cpSample,ProtoLoadSample-method}
\alias{explore,cpSample,ProtoRemoveBaseline-method}
\alias{explore,xcmsRaw,ProtoGenProfile-method}
\alias{explore,xcmsSet,ProtoFindComps-method}
\alias{explore,xcmsSet,ProtoGroupComps-method}
\alias{explore,xcmsSet,ProtoNormalize-method}
\alias{explore,xcmsSet,ProtoRtcor-method}
\alias{explore,xcmsSet,ProtoSummarize-method}
\title{explore pipeline with graphics.}
\description{
  
}
\usage{
\S4method{explore}{cpPeaks,ProtoFindPeaks}(object, protocol, raw,
sample=NULL, residuals = TRUE, island=TRUE, gg = ggobi(), dev = dev.cur())
\S4method{explore}{cpPeaks,ProtoFindPeaksParabola}(object, protocol,
         raw, sample=NULL, residuals = TRUE,
                   island=TRUE, gg = ggobi(), dev = dev.cur())
\S4method{explore}{cpSample,ProtoLoadSample}(object, protocol,...)
\S4method{explore}{cpSample,ProtoRemoveBaseline}(object, protocol, raw,geom='l')
\S4method{explore}{xcmsRaw,ProtoGenProfile}(object, protocol) 
\S4method{explore}{xcmsSet,ProtoFindComps}(object, protocol,sample, residuals=FALSE,
                   island=FALSE,gg = ggobi()) 
\S4method{explore}{xcmsSet,ProtoGroupComps}(object, protocol, gg = ggobi()) 
\S4method{explore}{xcmsSet,ProtoNormalize}(object, protocol, gg = ggobi()) 
\S4method{explore}{xcmsSet,ProtoRtcor}(object, protocol, raw=NULL,
                   xscale=NULL,geom=c("line","heatmap"),sample=NULL,log=FALSE,gg = ggobi())
\S4method{explore}{xcmsSet,ProtoSummarize}(object, protocol,gg=ggobi())
}
\arguments{
}
\value{
}
\examples{
}
\author{Michael Lawrence, Tengfei Yin}

