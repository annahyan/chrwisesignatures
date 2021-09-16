## library(here)
## library(grid)
## library(lattice)
## library(gridExtra)
## library(tidyverse)
## library(RColorBrewer)

## a = readRDS(here('RDatas/checkmate.plot.input.RDS'))

## smp.mat = a[[1]][[1]]
## sig.dims = dim(smp.mat)

#' Function plots repeated estimates of correlation coefficients between
#' mutational singnatures
#' 
#' @param estimate.list List of repeated correlation estimates between signatures
#'
#' @import dplyr
#' @import tidyverse
#' @import ggplot2
#' 
#' @export

plot_all_experiments = function(all.estimates, title) {

    ## setting invariants 

    myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))
    sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(-1, 1))

    ## learning some parameters

    estimate.list = all.estimates[[1]]
    
    sig.dims = dim( estimate.list[[1]] )
    
    sig.lengths = sig.dims[1]
    
    active.squares = sig.lengths * (sig.lengths - 1) / 2

    mini.square.size = ceiling(sqrt(length(estimate.list)))
    smp.mat = estimate.list[[1]]

    exp.count.row = ceiling(sqrt(length(all.estimates)))
    
    ##     text.grobs = lapply(colnames(smp.mat), grid::textGrob)
    
    ## ggs is gonna be the list of plots

    ggs = list()
    k = 0

    for (i in 1:(sig.lengths - 1)) {
        for (j in (i+1):sig.lengths) {

            pps = lapply(all.estimates, function (estimate.list) {
                smp.line = sapply(estimate.list, function(x) x[i,j])
                
                pp = smp.line %>%
                    matrix(ncol = mini.square.size) %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column() %>%
                    tidyr::pivot_longer(-c(rowname)) %>%
                    ggplot(aes(x = rowname, y = name)) +
                    geom_raster(aes(fill = value)) +
                    sc + theme_void() +
                    theme(legend.position = "none")
            } ) 

            pps.arranged = gridExtra::arrangeGrob(grobs = pps, nrow = exp.count.row)

            k = k+1
            ggs[[k]] = pps.arranged
        }
    }
    ##  return(pps.arranged)
    layout.mat = matrix(NA, ncol = sig.lengths, nrow = sig.lengths)

    layout.mat[lower.tri(layout.mat)] = 1:active.squares

    layout.mat[2:sig.lengths, 1:(sig.lengths - 1)] =
        t(layout.mat[2:sig.lengths, 1:(sig.lengths - 1)])

    ##     grid.arrange(grobs = ggs, layout_matrix = layout.mat)
    
                                        # other.list = ggs
    ##     ggs[seq(k + 1, length = sig.lengths) ] = text.grobs[1:sig.lengths]
    
    layout.mat[1, 1:(sig.lengths - 1)] = seq(active.squares + sig.lengths, length = sig.lengths - 1)
    layout.mat[2:sig.lengths, sig.lengths] = seq(active.squares + 1, length = sig.lengths - 1)

    for (i in 1:(sig.lengths-1)) {
        ggs[[(active.squares + i)]] = grid::textGrob(colnames(smp.mat)[i])
    }

    for (i in 2:sig.lengths) {
        ggs[[(active.squares + sig.lengths - 2  + i )]] = grid::textGrob(colnames(smp.mat)[i])
    }

    gout = gridExtra::grid.arrange(grobs = ggs, layout_matrix = layout.mat,
                                   top = grid::textGrob(title))

    return(gout)
}
