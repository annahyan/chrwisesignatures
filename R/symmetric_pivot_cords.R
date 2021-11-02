#' A function returning the symmetric pivots for i-th and j-th cols
#'
#'
#' @param x Data.frame or a matrix of signature intensities(or counts).
#' Signatures need to correspond to columns.
#'
#' @param i Index of a signature for which the pivots should be calculated
#' @param j Same as i
#'
#' @return A matrix, where the first two cols are the
#' symmetric coordinates
#'
#' @export


symmetric_pivot_coords = function(x, i, j,  ...) {

    if (any(x[!is.na(x)] <= 0)) 
        stop("all elements of x must be greater than 0")
    
    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    
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
    balZavout = balZav(x[, c(i, j, ind[-c(i, j)])])

    return(symm.pivots = balZavout)

    ## ind <- c(1:ncol(x))
    ## corZav <- matrix(NA, ncol(x), ncol(x))
    ## corPvals <- matrix(NA, ncol(x), ncol(x))
    ## for (i in 1:(ncol(x) - 1)) {
    ##     for (j in (i + 1):ncol(x)) {
    ##         balZavout = balZav(x[, c(i, j, ind[-c(i, j)])])

    ##         balZavout = balZavout[is.finite(rowSums(balZavout)), ]


    ##         if (mi == TRUE) { ### check if the mpmi MI calculation should be used

    ##             mi.out = mpmi::cmi(balZavout)

    ##             corZav[i, j] <- mi.out$bcmi[1,2]
    ##             corPvals[i, j] <- 2 * ( -abs( min(abs(mi.out$zvalues[1,2]), 100) ))
                
    ##         } else { ### cor.test is used instead
                
                
    ##             corout = cor.test(balZavout[,1], balZavout[, 2], ...)
    ##             corPvals[i, j]  <-  corout$p.value
    ##             corZav[i, j] <- corout$estimate
    ##         }
    ##     }
    ## }

    ## if(p.adjust) corPvals = p.adjust(corPvals, method = "BH")
    
    ## corZav[corPvals > p.val ] = 0
    ## corZav[lower.tri(corZav)] <- t(corZav)[lower.tri(corZav)]
    ## diag(corZav) <- 1
    ## return(corZav)
}
