#' Functions for computing the power
#'
#' \code{getPowerMultAM} - Computes the power of a single TSS by calculating
#' the proportion of values outside the region of acceptance for one or multiple
#' alternative models.
#'
#' @param accept_region Numeric vector of length two (increasing) setting the
#' lower and upper limit of the region of acceptance (limits included).
#' @param alts_data Numeric matrix with one column of values per alternative
#' model. If there is only one alternative model, then it can be a simple
#' vector of values instead (returns a single unnamed value).
#'
#' @return \code{getPower} A vector containing the power regarding the given
#' alternative models (retains column names of \code{alts_data}).
#'
#' @export
#' @rdname powerCompute
#'
#' @examples
#' getPowerMultAM(accept_region = c(10,20),
#'          alts_data = matrix(c(9,11,13,15,17,19,21,
#'                              5, 5, 5,15,25,25,25), ncol = 2, byrow = FALSE,
#'                            dimnames = list(NULL, c("AltTM1", "AltTM2"))))
#' getPowerMultAM(accept_region = c(10,20), alts_data = c(9,11,13,15,17,19,21))
getPowerMultAM <- function(accept_region, alts_data){
  if(accept_region[1] >= accept_region[2]){
    stop(paste0("accept_region not usable: c(",accept_region[1],",",
                accept_region[2],")"))
  }
  if(is.null(alts_data) || length(alts_data)<1){
    stop(paste("alts_data not usable."))
  }
  if(is.null(dim(alts_data))){
    alts_data <- matrix(alts_data, ncol = 1)
  }
  powers <- rep(NA,ncol(alts_data))
  names(powers) <- dimnames(alts_data)[[2]]
  for(alt in 1:ncol(alts_data)){
    powers[alt] <- sum(alts_data[,alt]<accept_region[1] |
                         alts_data[,alt]>accept_region[2])  / nrow(alts_data)
  }
  return(powers)
}
#' Functions for computing the power
#'
#' \code{getPowerMultTSS} - Computes the power of one or multiple TSS by
#' calculating the proportion of values outside the region of acceptance for
#' a single alternative model.
#'
#' @param accept_regions Numeric matrix with two columns. Each column represents
#' an interval, i.e., lower and upper limit of the region of acceptance (limits
#' included), of a different TSS.
#' @param alt_data Numeric matrix with one row of values under the
#' alternative model per TSS.
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
#' getPowerMultTSS(accept_regions = matrix(c(10,20,
#'                                    110,120), ncol = 2, byrow = TRUE,
#'                            dimnames = list(c("TSS1", "TSS2"),NULL)),
#'          alt_data = matrix(c(  9, 14, 19, 24,
#'                              109,114,119,124), nrow = 2, byrow = TRUE,
#'                            dimnames = list(c("TSS1", "TSS2"),NULL)))
#' getPowerMultTSS(accept_regions = c(10,20), alt_data = c(9, 14, 19, 24))
getPowerMultTSS <- function(accept_regions, alt_data){
  if(is.null(alt_data) || length(alt_data)<1){
    stop(paste("alt_data not usable."))
  }
  if(is.null(dim(accept_regions))){
    accept_regions <- matrix(accept_regions, nrow = 1)
    alt_data <- matrix(alt_data, nrow = 1)
  }
  if(sum(accept_regions[,1] >= accept_regions[,2])>0){
    stop(paste("accept_regions not usable."))
  }
  powers <- rep(NA, nrow(accept_regions))
  names(powers) <- dimnames(accept_regions)[[1]]
  for(tss in 1:nrow(accept_regions)){
    powers[tss] <- sum(alt_data[tss,]<accept_regions[tss,1] |
                         alt_data[tss,]>accept_regions[tss,2])  / ncol(alt_data)
  }
  return(powers)
}
