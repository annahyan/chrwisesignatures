#' Pairwise signature correlations from symmetric pivot coordinates
#'
#' @details For each pair of signatures correlation based on symmetric pivot
#' coordinates. The non-CoDa version of the function is \code{\link{cor_sigs}}.
#'
#' @param x Data.frame or a matrix of signature intensities(or counts).
#' Signatures need to correspond to columns.
#'
#' @param min.mut Minimum number of mutations required from the signatures.
#' Samples where both i-th and j-th signatures have less than this number of
#' mutations are excluded before cor.test calculation.
#'
#' @param p.val All the correlations with >=p.val are set to 0.
#' @param rand.add If rand.add = TRUE, then the matrix is strictly positive,
#' otherwise rand.add = FALSE(default). Check what you're  
#' @param p.adjust If TRUE, p.value adjustment with BH correction is applied
#' before setting correlations with p.value >= p.val to 0.
#' @param ... arguments passed to cor.test
#'
#' @return Signature correlation matrix
#'
#' @import mpmi
#' 
#' @export


cor_coda = function(x, min.mut = 0, p.val = 0.05, rand.add = FALSE,
                    p.adjust = TRUE, mi = FALSE,  ...) {

    if (any(x[!is.na(x)] <= 0)) 
        stop("all elements of x must be greater than 0")
    
    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    
    if (ncol(x) <= 2) 
        stop("calculation of average symmetric coordinates not possible")

    ind <- c(1:ncol(x))
    corZav <- matrix(NA, ncol(x), ncol(x))
    corPvals <- matrix(NA, ncol(x), ncol(x))
    
    for (i in 1:(ncol(x) - 1)) {
        for (j in (i + 1):ncol(x)) {

            permuted_x =  x[, c(i, j, ind[-c(i, j)])]
            
            ## Selecting only rows where at least one of the
            ## signatures is positive (min.mut)
            permuted_x = permuted_x[ permuted_x[,1] > min.mut |
                                     permuted_x[,2] > min.mut, ]

            if (nrow(permuted_x) < 3 ) {
                corPvals[i, j]  <-  1
                corZav[i, j] <- 0
                next
            }


            if (rand.add) { ## If random noise was added, then symm_coords are
                ## calculated for samples individually
                ## ATM if one of the pair is 0, both are returned as 0
                balZavout = apply( as.data.frame(permuted_x), MARGIN = 1,
                                  function(rowvec) {
                                      if (rowvec[1] == 0 | rowvec[2] == 0 ) return(c(0,0))
                                      rowvec = rowvec[ which (rowvec > 0) ]
                                      if (length(rowvec) == 2) return(c(0,0))
                                      return(balZavRow(rowvec))
                                  } ) %>% t()
                balZavout = balZavout[rowSums(abs(balZavout)) > 0, ]
                
            } else {
            
                balZavout = balZav(permuted_x)

                balZavout = balZavout[is.finite(rowSums(balZavout)), ]
            }

            if (mi == TRUE) { ### check if the mpmi MI calculation should be used

                mi.out = mpmi::cmi(balZavout)

                corZav[i, j] <- mi.out$bcmi[1,2]
                corPvals[i, j] <- 2 * ( -abs( min(abs(mi.out$zvalues[1,2]), 100) ))
                
            } else { ### cor.test is used instead
                
                
                corout = cor.test(balZavout[,1], balZavout[, 2], ...)
                
                corPvals[i, j]  <-  corout$p.value
                corZav[i, j] <- corout$estimate
            }
        }
    }

    if(p.adjust) {
        corPvals = p.adjust(corPvals, method = "BH")
    }
    
    corZav[corPvals > p.val ] = 0
    corZav[lower.tri(corZav)] <- t(corZav)[lower.tri(corZav)]
    diag(corZav) <- 1
    return(corZav)
}
