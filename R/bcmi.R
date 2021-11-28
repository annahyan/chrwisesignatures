#' Calculating bias-corrected mutual information for signature pairs
#'
#' @param x Data.frame or a matrix of signature intensities(or counts).
#' Signatures need to correspond to columns.
#'
#' @param p.val All the correlations with >=p.val are set to 0.
#' @param p.adjust If TRUE, p.value adjustment with BH correction is applied
#' before setting correlations with p.value >= p.val to 0.
#'
#' @return BCMI matrix
#'
#' @import mpmi
#' 
#' @export

bcmi = function(x,  p.val = 0.05, p.adjust = TRUE) {

    if (any(x[!is.na(x)] < 0)) {
        stop("all elements of x must be greater than or equal to 0")
    }
    
    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    
    if (ncol(x) <= 2) 
        stop("calculation of average symmetric coordinates not possible. Less than two columns.")


    cmi_out = mpmi::cmi(x)

    bcmi_out = cmi_out$bcmi
    zvals_out = cmi_out$zvalues

    pvals = pnorm(-abs(zvals_out) )

    if ( p.adjust) {
        pvals = p.adjust(pvals, method = "BH")
    }

    pvals = matrix(pvals, ncol = ncol(zvals_out))

    filtered_bcmi = bcmi_out
    filtered_bcmi[pvals > p.val] = 0

    colnames(filtered_bcmi) = colnames(x)
    rownames(filtered_bcmi) = colnames(x)
    
    return(filtered_bcmi)
}
