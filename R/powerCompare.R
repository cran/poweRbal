#' Comparison of the power of TSS under different models
#'
#' \code{powerComp} - Compare the power of a set of TSS to identify trees
#' generated under different alternative models given a null model.
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
#' @param with_CIradius Determines if the radii for the 95%-confidence
#' intervals of the power values should also be computed
#' and returned (default = FALSE).
#'
#' @return \code{powerComp}  Numeric matrix, one row per TSS and one
#' column per alternative model if \code{with_CIradius} is FALSE (default).
#' Else a list with two such numeric matrices, the first for the power and the
#' second for the confidence interval radius.
#'
#' @export
#' @rdname powerCompare
#'
#' @examples
#' powerComp(tss = c("Sackin", "Colless", "B1I"),
#'           alt_models = list(list("aldous",-1), "pda", "etm"), n=10L,
#'           sample_size = 40L, with_CIradius = TRUE)
powerComp <- function(tss,
                      null_model = "yule",
                      alt_models,
                      n, distribs = "sampled", sample_size = 1000L,
                      test_type = "two-tailed", sig_lvl = 0.05,
                      with_CIradius = FALSE){
  if(sum(is.na(c(tss,null_model,alt_models)))>0){
    stop(paste("Comparison impossible with given tss, null_model and",
               "alt_models."))
  }
  if(!distribs %in% c("sampled")){
    stop(paste("Unknown distribs method."))
  }
  if(sample_size < 10 || sample_size %%1!=0){
    stop(paste("The sample size must be an integer >=10."))
  }
  if(!test_type %in% c("two-tailed")){
    stop(paste("Unknown test type."))
  }
  if(sig_lvl >= 1 || sig_lvl <= 0){
    stop(paste("The level of significance must be >0 and <1."))
  }
  # Compute quantiles for null_model. ------------------------------------------
  acc_regions <- getAccRegion(tss = tss, null_model = null_model, n = n,
                              distribs = distribs, sample_size = sample_size,
                              test_type = test_type, sig_lvl = sig_lvl)
  # Compute the power of the TSS -----------------------------------------------
  powers <- matrix(NA, nrow = length(tss), ncol = length(alt_models),
                   dimnames = list(tss, rep(NA,length(alt_models))))
  if(with_CIradius){
    radii <- matrix(NA, nrow = length(tss), ncol = length(alt_models),
                    dimnames = list(tss, rep(NA,length(alt_models))))
  }
  act_sample_sizes <- rep(NA, length(alt_models))
  # Do this for each alternative model and fill each column --------------------
  for(i in 1:length(alt_models)){
    alt_data <- getTSSdata(tss = tss, n = n, Ntrees = sample_size,
                            tm = alt_models[[i]])
    if(!is.null(alt_data)){
      colnames(powers)[i] <- alt_models[[i]][[1]]
      powers[,i] <- getPowerMultTSS(accept_regions = acc_regions,
                                    alt_data = alt_data)
      act_sample_sizes[i] <- ncol(alt_data)
      if(with_CIradius){
        colnames(radii)[i] <- alt_models[[i]][[1]]
        how_many_rejected <- round(act_sample_sizes[i] * powers[,i])
        radii[,i] <- 1.96/ sqrt(act_sample_sizes[i]) *
          sqrt(how_many_rejected*(act_sample_sizes[i]-how_many_rejected)/
                 act_sample_sizes[i]/(act_sample_sizes[i]-1))
      }
    }
  }
  if(with_CIradius){
    return(list(power = powers, CIradius = radii,
                actual_sample_sizes = act_sample_sizes))
  } else {
    return(powers)
  }
}
