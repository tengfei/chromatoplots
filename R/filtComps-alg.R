######################################################################
## component filter ##
######################################################################
filtComps.cutoff <- function(object,tic.cutoff,npeaks.cutoff,tic,npeaks,...){
  test <- object
  test <- input_quantity(test)
  comps <- test@comps
  peaks <- test@peaks
  compList <- by(comps,comps[,'sample'],function(comps){
    comps <- as.data.frame(comps)
    comps[,'inter']<- as.character(interaction(comps[,'sample'],comps[,'comp']))
    if(tic.cutoff!=0){
      tic.cut <- quantile(comps[,'quantity'],tic.cutoff)
    }else{
      tic.cut <- tic
    }
    
    if(npeaks.cutoff!=0){
      npeaks.cut <- quantile(comps[,'npeaks'],npeaks.cutoff)
    }else{
      npeaks.cut <- npeaks
    }
    print(npeaks.cut)
    rowLogic <- comps[,'quantity']>tic.cut & comps[,'npeaks']>npeaks.cut
    comps<- comps[rowLogic,]
    idx <- order(comps[,'sample'],comps[,'comp'])
    comps <- comps[idx,]
    comps
  })
  comps <- do.call("rbind",compList)
  peaks <- as.data.frame(peaks)
  peaks[,'inter'] <-  as.character(interaction(peaks[,'sample'],peaks[,'comp']))
  idx <- order(peaks[,'sample'],peaks[,'comp'])
  peaks <- peaks[idx,]
  rowLogic <- peaks[,'inter'] %in% comps[,'inter']
  peaks <- peaks[rowLogic,]
  coln <- colnames(comps)%in%c('quantity','inter')
  comps <- comps[,!coln]
  coln <- colnames(peaks)%in%c('inter')
  peaks <- peaks[,!coln]
  test@peaks <- data.matrix(peaks)
  test@comps <- data.matrix(comps)
  test
}
  

