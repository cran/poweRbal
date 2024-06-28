#' Functions for computing the region of acceptance
#'
#' \code{getAccRegion} - Computes the region of acceptance based on quantiles
#' for a specified level of significance and method. The critical region is
#' everything strictly outside this interval, i.e., the interval limits still
#' belong to the acceptance region.
#'
#' @author Sophie Kersting
#'
#' @param tss Vector containing the names (as character) of the tree shape
#' statistics that should be compared. You may either use the short names
#' provided in \code{tssInfo} to use the already included TSS, or use the
#' name of a list object containing similar information as the entries in
#' \code{tssInfo}. Example:\cr
#' Use \code{"new_tss"} as the name for the list object
#' \code{new_tss} containing at least the function
#' \code{new_tss$func = function(tree){...}},
#' and optionally also the information \code{new_tss$short},
#' \code{new_tss$simple}, \code{new_tss$name}, \code{new_tss$type},
#' \code{new_tss$only_binary}, and \code{new_tss$safe_n}.
#' @param null_model The null model that is to be used to determine the power
#' of the tree shape statistics. In general, it must be a function that
#' produces rooted binary trees in \code{phylo} format. \cr
#' If the respective model is included in this
#' package, then specify the model and its parameters by using a character
#' or list. Available are all options listed under parameter \code{tm} in
#' the documentation of function \code{genTrees} (type \code{?genTrees}).\cr
#' If you want to include your own tree model, then use the
#' name of a list object containing the function (with the two input parameters
#' \code{n} and \code{Ntrees}). Example: \cr
#' Use \code{"new_tm"} for the list object
#' \code{new_tm <- list(func = function(n, Ntrees){...})}.
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.
#' @param distribs Determines how the distributions and quantiles are computed.
#' Available are: \cr
#' - "sampled" (default): \code{sample_size} many trees are sample under the
#' respective models to determine the quantiles for the null model and how many
#' trees under the alternative models lie outside.
#' @param sample_size Sample size (integer >=10) if distributions are sampled
#' (default= 1000L).
#' @param test_type Determines the method. Available are: \cr
#' - "two-tailed" (default)
#' @param sig_lvl Level of significance (default=0.05, must be >0 and <1).
#'
#' @return \code{getAccRegion} Numeric matrix (one row per TSS) with two
#' columns for the interval limits of the acceptance region.
#'
#' @export
#' @rdname powerRegAcc
#'
#' @examples
#' getAccRegion(tss = c("Sackin", "Colless", "B1I"), n = 8L, sample_size = 100L)
getAccRegion <- function(tss,
                         null_model = "yule", n, distribs = "sampled",
                         sample_size = 1000L,
                         test_type = "two-tailed", sig_lvl= 0.05){
  if(distribs == "sampled"){
    tss_data <- getTSSdata(tss = tss, n = n, Ntrees = sample_size,
                           tm = null_model)
    return(getAccRegion_data(tss_data = tss_data, test_type = test_type,
                             sig_lvl = sig_lvl))
  } else {
    stop("Unknown distribs method.")
  }
}
#' Functions for computing the region of acceptance
#'
#' \code{getAccRegion_data} - Computes the region of acceptance for the
#' given TSS values based on quantiles for a specified level of significance and
#' method. The critical region is everything strictly outside this interval,
#' i.e., the interval limits still belong to the acceptance region.
#'
#' @param tss_data Numeric matrix of TSS values (one row per TSS). The row names
#' are used as names for the TSS.
#'
#' @author Sophie Kersting
#'
#' @return \code{getAccRegion_data} Numeric matrix (one row per TSS) with two
#' columns for the interval limits of the acceptance region.
#'
#' @export
#' @rdname powerRegAcc
#'
#' @examples
#' getAccRegion_data(getTSSdata(tss = c("Colless", "Sackin"), n = 8L,
#'                              Ntrees = 20L, tm = "yule"))
getAccRegion_data <- function(tss_data, test_type = "two-tailed",
                              sig_lvl = 0.05){
  if(is.null(dim(tss_data))){
    # Convert a vector of values to a matrix with one row.
    tss_data <- matrix(tss_data, nrow =1, byrow = TRUE,
                       dimnames = list("TSS", NULL))
  }
  acc_regions <- matrix(NA, nrow = nrow(tss_data), ncol = 2,
                        dimnames = list(rownames(tss_data),
                                        c("lower_limit","upper_limit")))
  # Calculate Acceptance Intervals
  if(test_type == "two-tailed"){
    for(i in 1:nrow(tss_data)){
      acc_regions[i,] <- unname(stats::quantile(tss_data[i,],
                                         probs = c(0+sig_lvl/2, 1-sig_lvl/2),
                                         na.rm = TRUE, names = TRUE, type = 7))
    }
    return(acc_regions)
  } else {
    stop("Unknown test type.")
  }
}
