#' Comparison of the power of TSS under different models
#'
#' \code{powerComp} - Compare the power of a set of TSS to identify trees
#' generated under different alternative models given a null model.
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
#' name of a list object containing the function (dependent on one parameter
#' \code{n}). Example: \cr
#' Use \code{"new_tm"} for the list object \cr
#' \code{new_tm <- list(func = function(n, Ntrees){...})}.
#' @param alt_models List containing the alternative models that are to be
#' used to determine the power of the tree shape statistics. Functions that
#' produce rooted binary trees in \code{phylo} format. The information of each
#' single model must be in the format described for \code{null_model}.
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.
#' @param distribs Determines how the distributions (and with that the
#' bounds of the critical region) are computed. Available are: \cr
#' - "exact_if_possible" (default): Tries to compute the exact distribution
#' under the null model if possible. Currently, this is only implemented for
#' \code{null_model = "yule"}, \code{"pda"}, or \code{"etm"}, and
#' \code{n}<=20. In all other cases the distribution is approximated
#' by sampling \code{N_null} many trees under the null model as in the
#' option "sampled" below. \cr
#' - "sampled": \code{N_null} many trees are sampled under the
#' null model to approximate the distribution.
#' @param N_null Sample size (integer >=10) if distributions are sampled
#' (default = 10000L).
#' @param N_alt Sample size (integer >=10) for the alternative models to
#' estimate the power (default = 1000L).
#' @param test_type Determines the method. Available are: \cr
#' - "two-tailed" (default): The lower and upper bound of the region of
#' acceptance are determined based on the (empirical) distribution function
#' such that P(TSS < lower bound) <= \code{sig_lvl}/2 and
#' P(TSS > upper bound) <= \code{sig_lvl}/2. See parameter \code{correction}
#' for specifying how conservative the test should be: the null
#' hypothesis can either be rejected only if the values are strictly outside of
#' this region of acceptance (can be too conservative) or it can also be
#' rejected (with certain probabilities) if the value equals the lower or
#' upper bound.\cr
#' - "two-tailed-unbiased": Experimental - Use with caution!\cr
#' The region of acceptance is optimized to yield an unbiased test, i.e., a test
#' that identifies non-null models with a probability of at least
#' \code{sig_lvl}.
#' The region of acceptance is determined similar to the default method.
#' However, it need not be symmetrical, i.e., not necessarily
#' cutting off \code{sig_lvl}/2 on both sides. Also see parameter
#' \code{correction} for specifying how conservative the test should be.
#' @param correction Specifies the desired correction method.
#' Available are: \cr
#' - "small-sample" (default): This method tries to ensure that the critical
#' region, i.e. the range of values for which the null hypothesis is rejected,
#' is as close to \code{sig_lvl} as possible (compared with "none" below, which
#' can be too conservative). The idea is that the null hypothesis is also
#' rejected with certain probabilities if the value matches the value of a
#' quantile. \cr
#' - "none": No correction method is applied. With that the test might be
#' slightly too conservative as the null hypothesis is maintained if the values
#' is >= the lower and <= the upper quantile.
#' @param sig_lvl Level of significance (default = 0.05, must be >0 and <1).
#'
#' @return \code{powerComp}  Returns an object of class 'poweRbal_data' which
#' is a list containing the following objects:\cr
#' - power: Numeric matrix containing the power values (one row per TSS and one
#' column per alternative model).\cr
#' - accept_regions: Numeric matrix containing the information on the region of
#' acceptance (one row per TSS and four columns).\cr
#' - CIradius: Numeric matrix containing the confidence interval radii (one row
#' per TSS and one column per alternative model).\cr
#' - actual_sample_sizes: Numeric vector containing the actual sample sizes
#' under each alternative model as some models do not always successfully
#' generate trees. \cr
#' - other input data.
#'
#' @references
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic models.
#'  arXiv:2406.05185, 2024.
#'
#' @export
#' @rdname powerCompare
#'
#' @examples
#' powerComp(tss = c("Sackin", "Colless", "B1I"),
#'           alt_models = list(list("aldous",-1), "pda", "etm"), n = 10L,
#'           distribs = "sampled", N_null = 40L, N_alt = 20L)
powerComp <- function(tss,
                      null_model = "yule",
                      alt_models, n,
                      distribs = "exact_if_possible",
                      N_null = 10000L, N_alt = 1000L,
                      test_type = "two-tailed", correction = "small-sample",
                      sig_lvl = 0.05){
  if(sum(is.na(c(tss,null_model,alt_models)))>0){
    stop(paste("Comparison impossible with given tss, null_model and",
               "alt_models."))
  }
  if(!distribs %in% c("exact_if_possible", "sampled")){
    stop(paste("Unknown distribs method."))
  }
  if(N_null < 10 || N_null %%1!=0 || N_alt < 10 || N_alt %%1!=0){
    stop(paste("The sample sizes must be an integer >=10."))
  }
  if(!test_type %in% c("two-tailed", "two-tailed-unbiased")){
    stop(paste("Unknown test type."))
  }
  if(!correction %in% c("small-sample", "none")){
    stop(paste("Unknown test correction method."))
  }
  if(sig_lvl >= 1 || sig_lvl <= 0){
    stop(paste("The level of significance must be >0 and <1."))
  }
  # Compute quantiles for null_model. ------------------------------------------
  acc_regions <- getAccRegion(tss = tss, null_model = null_model, n = n,
                              distribs = distribs, N_null = N_null,
                              N_alt = N_alt,
                              test_type = test_type, correction = correction,
                              sig_lvl = sig_lvl)
  return(powerComp_RegAcc(tss = tss, accept_regions = acc_regions,
              null_model = null_model, alt_models = alt_models, n = n,
              distribs = distribs, N_null = N_null, N_alt = N_alt,
              test_type = test_type, correction = correction,
              sig_lvl = sig_lvl))
}
#' Comparison of the power of TSS under different models
#'
#' \code{powerComp_RegAcc} - Compare the power of a set of TSS to identify trees
#' generated under different alternative models given a the region(s) of
#' acceptance.
#'
#' @param accept_regions Numeric matrix (one row per TSS) with two or four
#' columns: The first two columns contain the interval limits of the region
#' of acceptance, i.e., we reject the null hypothesis for values strictly
#' outside of this interval. The third and fourth columns contain the
#' probabilities to reject the null hypothesis if values equal the lower or
#' upper bound, respectively. If the last two columns are missing they are
#' interpreted as zeroes. See return value of \code{getAccRegion()}.
#'
#' @return \code{powerComp_RegAcc} Returns an object of class 'poweRbal_data'
#' similar to \code{powerComp}.
#'
#' @export
#' @rdname powerCompare
#'
#' @examples
#' powerComp_RegAcc(tss = c("Sackin", "Colless", "B1I"),
#'           accept_regions = getAccRegion(tss = c("Sackin", "Colless", "B1I"),
#'                                         n = 6L, null_model = "etm",
#'                                         N_null = 20L, distribs = "sampled"),
#'           null_model = "etm", distribs = "sampled",
#'           alt_models = list(list("aldous",-1), "pda", "yule"), n = 6L,
#'           N_null = 20L, N_alt = 20L)
powerComp_RegAcc <- function(tss, accept_regions, null_model, alt_models, n,
                             distribs = "exact_if_possible",
                             N_null = 10000L, N_alt = 1000L,
                             test_type = "two-tailed",
                             correction = "small-sample", sig_lvl = 0.05){
  if(sum(is.na(c(tss,null_model,alt_models)))>0){
    stop(paste("Comparison impossible with given tss, null_model and",
               "alt_models."))
  }
  if(!distribs %in% c("exact_if_possible", "sampled")){
    stop(paste("Unknown distribs method."))
  }
  if(N_null < 10 || N_null %%1!=0 || N_alt < 10 || N_alt %%1!=0){
    stop(paste("The sample sizes must be an integer >=10."))
  }
  if(!test_type %in% c("two-tailed", "two-tailed-unbiased")){
    stop(paste("Unknown test type."))
  }
  if(!correction %in% c("small-sample", "none")){
    stop(paste("Unknown test correction method."))
  }
  if(sig_lvl >= 1 || sig_lvl <= 0){
    stop(paste("The level of significance must be >0 and <1."))
  }
  acc_regions <- accept_regions
  # Compute the power of the TSS -----------------------------------------------
  powers <- matrix(NA, nrow = length(tss), ncol = length(alt_models),
                   dimnames = list(tss, rep(NA,length(alt_models))))

  radii <- matrix(NA, nrow = length(tss), ncol = length(alt_models),
                  dimnames = list(tss, rep(NA,length(alt_models))))

  act_sample_sizes <- rep(NA, length(alt_models))
  # Do this for each alternative model and fill each column --------------------
  for(i in 1:length(alt_models)){
    alt_data <- getTSSdata(tss = tss, n = n, Ntrees = N_alt,
                           tm = alt_models[[i]])
    if(!is.null(alt_data)){
      colnames(powers)[i] <- paste(unlist(alt_models[[i]]), collapse=", ")
      powers[,i] <- getPowerMultTSS(accept_regions = acc_regions,
                                    alt_data = alt_data)
      act_sample_sizes[i] <- ncol(alt_data)
      colnames(radii)[i] <- paste(unlist(alt_models[[i]]), collapse=", ")
      how_many_rejected <- round(act_sample_sizes[i] * powers[,i])
      radii[,i] <- 1.96/ sqrt(act_sample_sizes[i]) *
        sqrt(how_many_rejected*(act_sample_sizes[i]-how_many_rejected)/
               act_sample_sizes[i]/(act_sample_sizes[i]-1))

    }
  }
  names(act_sample_sizes) <- dimnames(powers)[[2]]
  result <- list(power = powers, accept_regions = acc_regions, CIradius = radii,
                 actual_sample_sizes = act_sample_sizes, tss = tss,
                 null_model = null_model, alt_models = alt_models, n = n,
                 distribs = distribs, N_null = N_null, N_alt = N_alt,
                 test_type = test_type, correction = correction,
                 sig_lvl = sig_lvl)
  class(result) <- "poweRbal_data"
  return(result)
}
