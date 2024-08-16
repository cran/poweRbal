#' Function for computing the power given the region of acceptance
#'
#' \code{getPowerMultTSS} - Computes the power of one or multiple TSS by
#' calculating the proportion of values outside the region of acceptance for
#' a single alternative model.
#'
#' @param accept_regions Numeric matrix (one row per TSS) with two or four
#' columns: The first two columns contain the interval limits of the region
#' of acceptance, i.e., we reject the null hypothesis for values strictly
#' outside of this interval. The third and fourth columns contain the
#' probabilities to reject the null hypothesis if values equal the lower or
#' upper bound, respectively. If the last two columns are missing they are
#' interpreted as zeroes. See return value of \code{getAccRegion()}.
#' @param alt_data Numeric matrix  (one row per TSS) with values under the
#' alternative model.
#' If there is only one TSS, then it can be a simple vector of values instead
#' (returns a single unnamed value).
#'
#' @return \code{getPowerMultTSS} A vector containing the power regarding the
#' given TSS (retains row names of \code{accept_regions}).
#'
#' @export
#' @rdname powerCompute
#'
#' @examples
#' # Example with small data (with/without third and fourth column):
#' getPowerMultTSS(accept_regions = c(2,3, 0,0), alt_data = c(1,2,4,5))
#' getPowerMultTSS(accept_regions = c(2,3, 0.5,1), alt_data = c(1,2,4,5))
#' # Example with multiple rows/TSS:
#' getPowerMultTSS(accept_regions = matrix(c(2,3,0,0,
#'                                           20,30,0.5,0.5),
#'                                    nrow = 2, byrow = TRUE,
#'                                    dimnames = list(c("TSS1", "TSS2"), NULL)),
#'                 alt_data = matrix(c(  1,2,3,4,
#'                                     10,20,30,40),
#'                                   nrow = 2, byrow = TRUE,
#'                                   dimnames = list(c("TSS1", "TSS2"), NULL)))
#' # Example with generated TSS data:
#' getPowerMultTSS(accept_regions = getAccRegion(tss = c("Colless","SNI"),
#'                                               n = 6L),
#'                 alt_data = getTSSdata(tss = c("Colless", "SNI"), n = 6L,
#'                                       Ntrees = 20L, tm = list("aldous", -1)))
getPowerMultTSS <- function(accept_regions, alt_data){
  if(is.null(alt_data) || length(alt_data)<1){
    stop("alt_data not usable (not available).")
  }
  if(is.null(dim(accept_regions))){
    accept_regions <- matrix(accept_regions, nrow = 1)
    alt_data <- matrix(alt_data, nrow = 1)
  }
  if(sum(accept_regions[,1] > accept_regions[,2])>0){
    stop("accept_regions not usable (lower bound > upper bound).")
  }
  if(dim(accept_regions)[2]!=4){
    if(dim(accept_regions)[2]==2){
      accept_regions <- cbind(accept_regions,
                              rep(0,nrow(accept_regions)),
                              rep(0,nrow(accept_regions)))
    } else {
      stop("accept_regions not usable (wrong number of columns).")
    }
  }
  powers <- rep(NA, nrow(accept_regions))
  names(powers) <- dimnames(accept_regions)[[1]]
  for(tss in 1:nrow(accept_regions)){
    powers[tss] <- sum(alt_data[tss,]<accept_regions[tss,1] |
                         alt_data[tss,]>accept_regions[tss,2]) +
                   sum(alt_data[tss,]==accept_regions[tss,1]) *
                                  accept_regions[tss,3] +
                   sum(alt_data[tss,]==accept_regions[tss,2]) *
                                  accept_regions[tss,4]
    powers[tss] <- powers[tss] / ncol(alt_data)
  }
  return(powers)
}
