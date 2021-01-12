#' Pairwise signature correlations from symmetric pivot coordinates
#'
#' @details For each pair of signatures correlation based on symmetric pivot
#' coordinates. The non-CoDa version of the function is \code{\link{cor_sigs}}.
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
#' @import mpmi
#' 
#' @export


cor_coda = function(x, p.val = 0.05, p.adjust = TRUE, mi,  ...) {


    if ( ! missing(mi) ) {
        mi = TRUE
    }
    
    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    
    if (any(x[!is.na(x)] <= 0)) 
        stop("all elements of x must be greater than 0")
    
    if (ncol(x) <= 2) 
        stop("calculation of average symmetric coordinates not possible")
    
    balZav <- function(x) {
        D <- ncol(x)
        Z.av <- matrix(NA, ncol = 2, nrow = nrow(x))
        p1 <- sqrt(D - 1 + sqrt(D * (D - 2)))/sqrt(2 * D)
        if (D == 3) {
            p2 <- x[, 3]
        }
        else {
            p2 <- apply(x[, 3:D], 1, prod)
        }
        p3 <- (sqrt(D - 2) + sqrt(D))/(sqrt(D - 2) * (D - 1 + 
            sqrt(D * (D - 2))))
        p4 <- 1/(D - 1 + sqrt(D * (D - 2)))
        Z.av[, 1] <- p1 * (log(x[, 1]/(x[, 2]^p4 * p2^p3)))
        Z.av[, 2] <- p1 * (log(x[, 2]/(x[, 1]^p4 * p2^p3)))

        return(Z.av)
    }


    ind <- c(1:ncol(x))
    corZav <- matrix(NA, ncol(x), ncol(x))
    corPvals <- matrix(NA, ncol(x), ncol(x))
    for (i in 1:(ncol(x) - 1)) {
        for (j in (i + 1):ncol(x)) {
            balZavout = balZav(x[, c(i, j, ind[-c(i, j)])])

            balZavout = balZavout[is.finite(rowSums(balZavout)), ]


            if (mi == TRUE) { ### check if the mpmi MI calculation should be used

                mi.out = mpmi::cmi(balZavout[,1], balZavout[, 2])

                corZav[i, j] <- mi.out$bcmi
                corPvals[i, j] <- 2 * ( -abs( mi.out$zscore ))
                
            } else { ### cor.test is used instead
                
                
                corout = cor.test(balZavout[,1], balZavout[, 2], ...)
                corPvals[i, j]  <-  corout$p.value

                corZav[i, j] <- corout$estimate
            }
        }
    }

    if(p.adjust) corPvals = p.adjust(corPvals, method = "BH")
    
    corZav[corPvals > p.val ] = 0
    corZav[lower.tri(corZav)] <- t(corZav)[lower.tri(corZav)]
    diag(corZav) <- 1
    return(corZav)
}
