#' Generation of rooted binary trees under the BiSSE model
#'
#' \code{genBiSSETree} - Generates a rooted binary tree in \code{phylo}
#' format with the given number of \code{n} leaves under the BiSSE model.
#' In the BiSSE model all species have a state, either A or B, and depending on
#' the state a speciation rate \code{BIRTHRATES}, an extinction rate
#' \code{DEATHRATES} as well as a transition rate to the other state
#' \code{TRANSRATES}. \cr
#' Extinct species are removed from the tree, i.e., the generated tree contains
#' only species living at the present.\cr
#' Due to the restrictions of the \code{phylo} or \code{multiphylo} format
#' the number of leaves must be at least 2 since there must be at
#' least one edge.
#'
#' @author Sophie Kersting
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.
#' @param BIRTHRATES Numeric vector (default = c(1,1)) which specifies
#' the speciation rates in state A and B (vector with 2 values >=0, one
#' value >0).
#' @param DEATHRATES Numeric vector (default = c(0,0)) which specifies
#' the extinction rates in state A and B (vector with 2 values >=0).
#' @param TRANSRATES Numeric vector which specifies the
#' transition rates from A to B and from B to A (vector with 2 values >0).
#' @param TRIES Integer value (default = 5) that specifies
#' the number of attempts to generate a tree with \code{n} leaves.
#' @param TIMEperTRY Numeric value (default = 0.1) that specifies the maximum
#' amount of time (in seconds) invested per try.
#'
#' @return \code{genBiSSETree} A single tree of class \code{phylo} is
#' returned.
#'
#' @references
#'  - This function uses the \code{tree.bisse} function of the
#'  \code{diversitree} package
#'  (R. G. FitzJohn. Diversitree: Comparative Phylogenetic Analyses of
#'  Diversification in R.  Methods in Ecology and
#'  Evolution, 3(6):1084-1092, 2012). \cr
#'  - W. P. Maddison, P. E. Midford, and S. P. Otto. Estimating a binary
#'  character’s effect on speciation and extinction. Systematic Biology,
#'  56(5):701–710, 2007.
#'
#' @export
#' @rdname tmBiSSE
#'
#' @examples
#' if (requireNamespace("diversitree", quietly = TRUE)) {
#' genBiSSETree(n = 5, BIRTHRATES = c(1,2), DEATHRATES = c(0,0),
#'              TRANSRATES = c(0.1,0.3))
#' }
genBiSSETree <- function(n, BIRTHRATES = c(1,1), DEATHRATES = c(0,0),
                         TRANSRATES,
                         TRIES = 5, TIMEperTRY = 0.1){
  if(n < 2 || n%%1!=0){
    stop(paste("A tree must have at least 2 leaves, i.e., n>=2 and n must be",
               "an integer."))
  }
  if(length(BIRTHRATES)!=2 || sum(BIRTHRATES<0)> 0 || sum(BIRTHRATES=0)> 1){
    stop(paste("One of the two speciation rates must be >0 (both >=0)."))
  }
  if(length(DEATHRATES)!=2 || sum(DEATHRATES<0)> 0){
    stop(paste("The two extinction rates must be >=0."))
  }
  if(length(TRANSRATES)!=2 || sum(TRANSRATES<0)> 0){
    stop(paste("The two transition rates must be >=0."))
  }
  phy <- NULL
  i <- 1
  while(i<=TRIES && is.null(phy)){
    i <- i+1
    tryCatch({
      phy <- R.utils::withTimeout({
        test <- diversitree::tree.bisse(pars = c(BIRTHRATES, DEATHRATES,
                                                 TRANSRATES),
                                max.taxa = n, include.extinct = FALSE)
      }, timeout = TIMEperTRY)
    }, TimeoutException = function(ex) {
    NULL})
  }
  if(is.null(phy)){
    comment("Not able to generate tree. Returning NULL instead.")
  }
  return(phy)
}
