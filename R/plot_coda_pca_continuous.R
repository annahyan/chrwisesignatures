#' Classical or robust PCA after CoDa transformation for continuous annotations
#'
#' @details 
#' CoDa transformation is applied on 96-channel or 6-channel mutationa data,
#' as well as signature composition data before classical or robust PCA with
#' continuous annotations. The matrices have to be positive.
#'
#' @param dt_list Either list of mutational profiles or a data.frame with
#' compositional parts(signature intensities, mutation counts) as columns and
#' samples as rows.
#'
#' @param continuous A vector of annotation values to be color coded. The
#' length of this parameter has to be equal to the length of dt_list if list is
#' provided, or number of rows if data.frame provided.
#'
#' @param palette_name The palette to be passed to scale_color_distiller function.
#' default: Spectral. 
#' @param point_size Point size
#'
#' @param arrow_length the length of biplot arrows
#'
#' @param method The method for PCA: classical or robust. default:robust
#'
#' @importFrom robCompositions pcaCoDa 
#' 
#' @export


plot_coda_pca_continuous  <-  function(dt_list, continuous, palette_name, point_size, arrow_length, method) {

    if (missing(point_size)) {
        point_size = 3
    }
    
    if (missing(arrow_length)) {
        arrow_length = 2
    }
    if (missing(method)) {
        method = "robust"
    }

    if (missing(palette_name)) {
        palette_name = "Spectral"
    }
    
    if (inherits(dt_list, "list")) { 
        plot_material = do.call( rbind,
                                lapply(names(dt_list),
                                       function(x) {
                                           df = dt_list[[x]];
                                           rownames(df) = paste(x, rownames(df),
                                                                sep = ":");
                                           return(df)
                                       })
                                )

        ##        sample_types = rep(sample_classes, vapply(dt_list, nrow, numeric(1)))
    } else {
        plot_material = t(dt_list)
        ##        sample_types = sample_classes
    }


    
    pca_list = pcaCoDa(plot_material, method = method)
    outliers = outCoDa(plot_material)
    print(outliers)


    ## Getting proporition of variance explained
    proportion_variance = (pca_list$eigenvalues ^ 2) / sum((pca_list$eigenvalues ^ 2))
    percent_variance = 100 * proportion_variance

    
    rob_pca_scores = data.frame(pca_list$scores)
    rob_pca_scores$col = continuous

    rownames(rob_pca_scores) = rownames(plot_material)


    p = ggplot(rob_pca_scores, aes(x = Comp.1, y = Comp.2) ) + 
        geom_point(aes(colour = col), size = point_size) +
        scale_color_distiller(palette = palette_name ) # 
#        scale_color_manual(values = sample(CEMM_COLORS_ALL, 9 ) ) # +
    ## ggtitle(paste0("Rank:", nmf.rank))


    pca_loadings = as.data.frame(pca_list$loadings)
    rownames(pca_loadings) = colnames(plot_material)

    p = p + geom_segment(data = pca_loadings,
                         aes(x = 0, y = 0, xend = Comp.1 * arrow_length,
                             yend = Comp.2 * arrow_length),
                         size =0.5, 
                         arrow = arrow(length = unit(1/2, "picas") ),
                         color = "gray50") + 
        annotate("text", x = pca_loadings$Comp.1 * arrow_length * 1.1,
                 y = pca_loadings$Comp.2 * arrow_length * 1.1, 
                 label = rownames(pca_loadings)) +
        theme_bw() + theme(legend.title = element_blank()) +
        xlab(paste0("PC1", " (",round(percent_variance[1], 2), "%)" ) )  +
        ylab(paste0("PC2", " (",round(percent_variance[2], 2), "%)") )
    
    return(p)
}
