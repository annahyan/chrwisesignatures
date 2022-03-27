#' Applies filtering threshold on the matrix of interaction metrics
#'
#' @param calc_metrics A list of calculated metrics obtained by applying
#' metric_functions onto a data matrix.
#' @param filter_list A list of threshold values to be applied on the calc_metrics.
#' All the values with abslolute value smaller than the threshold is set to 0.
#' Default: NULL
#'
#' @examples
#' mat = matrix(sample(100, 300, replace = TRUE), 30, 10,
#' 		   dimnames = list(paste0("c", 1:30), letters[1:10]))
#' mat[, 2] = 2 * mat[,1] + mat[,2]
#' mat[runif(300) < 0.5] = 0
#' 
#' 
#' metric_functions = list(CoDa = cor_coda,
#'                         cooccurrence = cooccurrence,
#'                         MI = bcmi,
#'                         Spearman = function(x) {cor_sigs(x, method = "spearman")})
#' 
#' filter_threshold = list(MI = 0.2)
#' 
#' #'
#' metrics_out = lapply(metric_functions, function(x) x(mat))
#' 
#' 
#' filtered_metrics_out = apply_filter_list(calc_metrics = metrics_out,
#'                                          filter_list = filter_threshold)
#' 
#' @export 
apply_filter_list = function(calc_metrics, filter_list = NULL) {
    
    filtered_metrics = calc_metrics
    if (! is.null(filter_list)) {
        
        filter_list_absent = setdiff(names(filter_list), names(calc_metrics))
        
        if (length( filter_list_absent > 0))  {
            stop(paste("Filter_list has variable names not present:", 
                       paste(filter_list_absent, collapse = " ")) )
        } 
        
        for (var in names(filter_list)) {
            varlim = filter_list[[var]]
            mat = calc_metrics[[var]]
            mat[ abs(mat) < varlim ] = 0
            mat = mat[rowSums(mat) > 0, colSums(mat) > 0]
            filtered_metrics[[var]] = mat
        }        
    }
    return(filtered_metrics)
}
