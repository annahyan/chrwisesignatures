#' Getting a matrix of interactions and returning an interaction network
#' 
#' @param data_matrix input matrix for which to calculate the metrics and 
#' demonstrate with a network
#' @param metric_functions is a list with elements as functions which return the
#' calculated interactions across signatures. 
#'
#'
#' @examples
#' set.seed(10)
#' mat = matrix(sample(100, 300, replace = TRUE), 30, 10,
#' 		   dimnames = list(paste0("c", 1:30), letters[1:10]))
#' mat[, 2] = 2 * mat[,1] + mat[,2]
#' mat[runif(300) < 0.5] = 0
#' 
#' metric_functions = list(CoDa = cor_coda,
#'                         cooccurrence = cooccurrence,
#'                         MI = bcmi,
#'                         Spearman = function(x) {cor_sigs(x, method = "spearman")})
#' 
#' plot_interaction_networks(mat, metric_functions, 
#'                                         filter_list = NULL,
#'                                         layout = "stress")
#' @export
plot_interaction_networks = function(data_matrix, metric_functions, 
                                        filter_list = NULL,
                                        layout = "stress") {

    calc_metrics = lapply(metric_functions, function(x) x(data_matrix))
    
    
    if (! is.null(filter_list)) {
        calc_metrics = apply_filter_list(calc_metrics = calc_metrics, 
                                         filter_list = filter_list)
    }
    
    network_graph = tissue_multi_graph(calc_metrics)
    p = print_multi_graphs(network_graph,  layout = layout)
    return(p)
}



