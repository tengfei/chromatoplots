setMethod("explore", c("xcmsSet", "ProtoSummarize"),
          function(object, protocol,gg=ggobi())
          {
            result <- matrix(NA, nrow(object@groups),
                             nrow(object@phenoData))
            comps_n <- object@comps[!is.na(object@comps[,"quantity"]),]
            group_n <- as.numeric(factor(comps_n[,"group"]))
            result[(comps_n[,"sample"]-1)*nrow(result) + group_n] <-
              log(comps_n[,"quantity"]+1)
            bad_groups <- apply(result, 1, function(x) all(is.na(x)))
            result_groups <- comps_n[!bad_groups,"group"]
            result <- result[!bad_groups,]

            build_diag <- function(g)
              c(mean = mean(g, na.rm=TRUE), 
                max_diff = max(abs(outer(g, g, "-")), na.rm=TRUE), 
                ncomps = sum(!is.na(g)))
            norm_diag <- t(apply(result, 1, build_diag))
            files <- object@filepaths
            colnames(result) <- sapply(files, basename)
            norm_diag <- as.data.frame(norm_diag)
            norm_diag$ncomps <- factor(norm_diag$ncomps)

            spectra <- t(find_spectra(object@peaks))
            colnames(spectra) <- paste("mz",seq_len(ncol(spectra)),sep="")
            gg["spectra"] <- as.data.frame(spectra)
            spectra_d <- gg["spectra"]

            gg["result"] <- result
            gg["result_diag"] <- norm_diag
            gg_result_diag <- display(gg["result_diag"], embed=TRUE)
            gg_result_bar <- display(gg["result_diag"], "Barchart", 
                                     list(X="ncomps"), embed=TRUE)
            variables(gg_result_bar) <- list(X="ncomps")
            gg_result <- display(gg["result"],  "Parallel Coordinates Display",
                                 embed=TRUE)
            gg_comp_peaks <- display(spectra_d, "Parallel Coordinates Display", 
                                     list(X=colnames(spectra)[1]), embed = TRUE)

            result_hbox <- gtkHBox()
            result_hbox$add(gg_result_diag)
            result_hbox$add(gg_result_bar)

            result_diag_d <- gg["result_diag"]
            result_d <- gg["result"]
            gSignalConnect(gg, "identify-point",
                           function(gg, plot, id, dataset) {
                             if (id == -1)
                               return()
                             if (!(as.RGtkObject(result_diag_d) == dataset))
                               return()
                            ##  result_ex <- rep(TRUE, nrow(result_d))
##                              result_ex[id+1] <- FALSE
##                              excluded(result_d) <- result_ex
                             glyph_color(result_d) <- 1
                             glyph_color(result_d)[id+1] <- 9
                             glyph_color(spectra_d) <- 1
                             glyph_color(spectra_d)[id+1] <- 9
                           })

            stack_plots("Summarization",
                        list(result_hbox, gg_result, gg_comp_peaks))
          })
