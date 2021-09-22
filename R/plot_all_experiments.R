## library(here)
## library(grid)
## library(lattice)
## library(gridExtra)
## library(tidyverse)
## library(RColorBrewer)

## a = readRDS(here('RDatas/checkmate.plot.input.RDS'))

## smp.mat = a[[1]][[1]]
## sig.dims = dim(smp.mat)

#' Plots repeated estimates of correlation coefficients between mutational
#' signatures
#' 
#' @param estimate.list List of repeated correlation estimates between signatures. Required.
#' @param title Plot title. Required.
#' @param rect.lwd Line width around rectangles. Default: 0.8.
#'
#' @import dplyr
#' @import tidyverse
#' @import ggplot2
#' 
#' @export

plot_all_experiments = function(all.estimates, title, rect.lwd, mc.cores) {

    ## setting invariants 

    if (missing(title)) {
        stop("title is missing.")
    }

    if (missing(rect.lwd)) {
        rect.lwd = 0.8
    }

    if (missing(mc.cores)) {
        mc.cores = 1
    }
    
    myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))
    sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(-1, 1))

    rect = grid::grid.rect(.5,.5,width=unit(.99,"npc"), height=unit(0.99,"npc"), 
                           gp=grid::gpar(lwd=rect.lwd, fill=NA, col="black"), draw = FALSE)
    

    ## Simplifying the experiments

    simplify_by_exp = lapply(all.estimates, simplify2array)
    all_simplified = simplify2array(simplify_by_exp)
    
    ## learning some parameters

    sig.dims = dim( all_simplified ) [1]
    
    sig.lengths = sig.dims[1]
    
    active.squares = sig.lengths * (sig.lengths - 1) / 2

    mini.square.size = ceiling( sqrt( dim(all_simplified)[3] ) )
    smp.mat = all.estimates[[1]][[1]]

    exp.count.row = ceiling(sqrt(length(all.estimates)))
    
    ##     text.grobs = lapply(colnames(smp.mat), grid::textGrob)
    
    ## ggs is gonna be the list of plots

    ggs = list()
    k = 0

    if (sig.lengths < 3 ) {
        next
    } 
    
    for (i in 1:(sig.lengths - 1)) {
        for (j in (i+1):sig.lengths) {
            cat(i, j, "\n")

            ptm = proc.time()

            pps = lapply( 1:dim(all_simplified)[4], function (exp_iter) {
                smp.line = all_simplified[i,j, ,exp_iter]
                
                pp = smp.line %>%
                    matrix(ncol = mini.square.size) %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column() %>%
                    tidyr::pivot_longer(-c(rowname)) %>%
                    ggplot(aes(x = rowname, y = name)) +
                    geom_raster(aes(fill = value)) +
                    sc + theme_void() +
                    theme(legend.position = "none",
                          panel.border = element_rect(colour = "gray90", fill = NA, size = 0.5))
            })

            print(proc.time() - ptm)
            
            pps.arranged = gridExtra::arrangeGrob(grobs = pps, nrow = exp.count.row)

            pps.arranged.rect = grid::gTree(children = grid::gList(pps.arranged, rect))
            k = k+1
            ggs[[k]] = pps.arranged.rect
            rm(pps)
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
                                   top = grid::textGrob(title), padding = unit(1, "line"))

    return(gout)
}
