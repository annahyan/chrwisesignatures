#' Simple correlation of signature
#'
#' @details All pairwise correlations are calculated by cor.test.
#' Non-CoDa version of \code{\link{cor_coda}} function, provided for comparison.
#'
#' @param x Data.frame or a matrix of signature intensities(or counts).
#' Signatures need to correspond to columns.
#'
#' @param p.val All the correlations with >=p.val are set to 0.
#' @param p.adjust If TRUE, p.value adjustment with BH correction is applied
#' before setting correlations with p.value >= p.val to 0.
#' @param ... arguments passed to cor.test
#'
#' @return Signature correlation matrix
#'
#' @export

cor_sigs = function(x, p.val = 0.05, p.adjust = TRUE,  ...) {

    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    ## if (any(x[!is.na(x)] <= 0)) 
    ##     stop("all elements of x must be greater than 0")
    ## if (ncol(x) <= 2) 
    ##     stop("calculation of average symmetric coordinates not possible")

    ind <- c(1:ncol(x))
    corav <- matrix(NA, ncol(x), ncol(x))
    corPvals <- matrix(NA, ncol(x), ncol(x))
    for (i in 1:(ncol(x) - 1)) {
        for (j in (i + 1):ncol(x)) {
            two_cols = x[, c(i, j)]

            corout = cor.test(two_cols[,1], two_cols[, 2], ...)

            corPvals[i, j]  <-  corout$p.value

            corav[i, j] <- corout$estimate
        }
    }
    
    if(p.adjust) corPvals = p.adjust(corPvals, method = "BH") 

    corav[corPvals > p.val ] = 0
    corav[lower.tri(corav)] <- t(corav)[lower.tri(corav)]
    diag(corav) <- 1

    colnames(corav) = colnames(x)
    rownames(corav) = colnames(x)
    
    return(corav)
}
