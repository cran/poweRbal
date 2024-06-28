#' Generation of the comb or caterpillar tree
#'
#' \code{genCombTree} - Generates the rooted binary comb tree (also known as
#' caterpillar tree) in \code{phylo} format with the given number of \code{n}
#' leaves. \cr
#' Due to the restrictions of the \code{phylo} or \code{multiphylo} format
#' the number of leaves must be at least 2 since there must be at
#' least one edge.
#'
#' @author Sophie Kersting
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.
#'
#' @return \code{genCombTree} A single tree of class \code{phylo} is returned.
#'
#' @references
#'  - D. Aldous. Probability Distributions on Cladograms. In Random Discrete
#'  Structures, pages 1–18. Springer New York, 1996.
#'
#'
#' @export
#' @rdname tmComb
#'
#' @examples
#' genCombTree(n = 6)
genCombTree <- function(n){
  if(n < 2 || n%%1!=0){
    stop(paste("A tree must have at least 2 leaves, i.e., n>=2 and n must be",
               "an integer."))
  }
  caterpillar_n <- "( , )"
  for(i in 3:n){
    caterpillar_n <- paste0("(",caterpillar_n,", )")
  }
  caterpillar_n <- paste0(caterpillar_n,";")
  return(phytools::read.newick(text = caterpillar_n))
}
