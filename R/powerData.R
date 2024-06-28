#' Functions for generating the TSS data under a tree model
#'
#' \code{getTSSdata} - Compute the tree shape statistics of trees generated
#' under a tree model for each given TSS.
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
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.
#' @param Ntrees Integer value (default = 1) that specifies the desired number
#' of generated trees.
#' @param tm If the respective model is included in this package, then specify
#' the model and its parameters by using a character or list. Available are all
#' options listed under parameter \code{tm} in
#' the documentation of function \code{genTrees} (type \code{?genTrees}).\cr
#' If you want to include your own tree model, then use the
#' name of a list object containing the function (with the two input parameters
#' \code{n} and \code{Ntrees}). Example: \cr
#' Use \code{"new_tm"} for the list object
#' \code{new_tm <- list(func = function(n, Ntrees){...})}.
#'
#' @return \code{getTSSdata} Numeric matrix of TSS values (one row per TSS).
#'
#' @export
#' @rdname powerData
#'
#' @examples
#' # Example using tree models included in this package.
#' getTSSdata(tss = c("Colless", "Sackin"), n = 5L, Ntrees = 3L,
#'            tm = "yule")
#' # Example using a "new" TM provided by the user.
#' my_aldous <- list(func = function(n, Ntrees){
#'                         trees <- lapply(1:Ntrees,
#'                            function(x){genAldousBetaTree(n = n, BETA =5L)})
#'                         attr(trees, "class") <- "multiPhylo"
#'                         return(trees)})
#' getTSSdata(tss = c("Colless", "Sackin"), n = 5L, Ntrees = 3L,
#'            tm = "my_aldous")
getTSSdata <- function(tss, n, Ntrees = 1L, tm) {
  # Create the list of trees under the model
  treeList <- NULL
  # If the model is included in this package.
  if((is.list(tm) && tm[[1]] %in% c("aldous","ford",
                                    "alt-birth-death", "density",
                                    "BiSSE", "DCO_sym", "DCO_asym",
                                    "IF_sym", "IF_asym", "IF-diff",
                                    "biased", "ASB", "simpleBrown_sym",
                                    "simpleBrown_asym", "lin-Brown_sym",
                                    "lin-Brown_asym", "lin-Brown-bounded_sym",
                                    "lin-Brown-bounded_asym","log-Brown_sym",
                                    "log-Brown_asym")) ||
     (is.character(tm) && tm %in% c("yule","pda","etm"))){
    treeList <- tryCatch(
      genTrees(n, Ntrees, tm),
      error = function(cond) { return(NULL) }
    )
  } else if (is.character(tm)) {
    # Else, if the model information is provided by the user.
    if(exists(tm)[[1]]){
      if(!is.null(get(tm)$func)) {
        treeList <- get(tm)$func(n=n, Ntrees = Ntrees)
      } else {
        stop(paste0("No function provided for tree model ",
                    paste(unlist(tm), collapse = ", ")))
      }
    }
  } else {
    stop(paste0("Unable to find tree model: ",
                paste(unlist(tm), collapse = ", ")))
  }
  # Get TSS values for the trees in the list.
  return(getTSSdata_trees(tss = tss, treeList = treeList))
}
#' Functions for generating the TSS data under a tree model
#'
#' \code{getTSSdata_trees} - Compute the tree shape statistics for each given
#' TSS and all given trees.
#'
#' @param treeList List of trees of class \code{multiphylo}.
#'
#' @return \code{getTSSdata_trees} Numeric matrix of TSS values
#' (one row per TSS).
#'
#' @export
#' @rdname powerData
#'
#' @examples
#' # Example using TSS provided in tssInfo.
#' getTSSdata_trees(tss = c("Colless", "Sackin"),
#'            treeList = lapply(1:6L, function(x) genYuleTree(5)))
#' # Example using a "new" TSS provided by the user.
#' my_avd <- list(func = treebalance::avgVertDep, short = "My AVD")
#' getTSSdata_trees(tss = c("Colless", "my_avd"),
#'            treeList = lapply(1:6L, function(x) genYuleTree(5)))
getTSSdata_trees <- function(tss, treeList) {
  if(is.null(treeList) || length(treeList) == 0){
    return(NULL)
  }
  # Initialize matrix for the TSS values.
  tss_data <- matrix(NA, nrow = length(tss), ncol = length(treeList),
                     dimnames = list(tss, NULL))
  for(i in 1:length(tss)){
    # Compute the TSS for all trees.
    if(!is.null(tssInfo[[tss[i]]])){
      # The TSS is already contained in tssInfo.
      tss_data[i,] <- sapply(1:length(treeList),
                             function(x){tssInfo[[tss[i]]]$func(treeList[[x]])})

    } else if(exists(tss[i])) {
      # The TSS is provided by the user.
      if(!is.null(get(tss[i])$func)) {
        tss_data[i,] <- sapply(1:length(treeList),
                               function(x){get(tss[i])$func(treeList[[x]])})
      } else {
        stop(paste0("No function provided for TSS ",tss[i],"."))
      }
      if(!is.null(get(tss[i])$short)) {
        # If a short name is provided use it.
        dimnames(tss_data)[[1]][i] <- get(tss[i])$short
      }
    } else {
      stop(paste("Unknown TSS:",tss[i]))
    }
  }
  return(tss_data)
}
