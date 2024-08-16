#' Generation of rooted binary trees under the PDA model
#'
#' \code{genPDATree} - Generates a rooted binary tree in \code{phylo}
#' format with the given number of \code{n} leaves under the
#' proportional-to-distinguishable-arrangements model.
#' Given \code{n}, all phylogenies (trees with labeled leaves) with
#' \code{n} leaves are equiprobable under the PDA.
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.\cr
#' Due to the restrictions of the \code{phylo} or \code{multiphylo} format,
#' the number of leaves must be at least 2 since there must be at
#' least one edge.
#'
#' @return \code{genPDATree} A single tree of class \code{phylo} is returned.
#'
#' @references
#'  - This function uses the \code{rtopology(..., rooted = T)} function of
#'  the \code{ape} package
#'  (E. Paradis, K. Schliep. “ape 5.0: an environment for modern
#'  phylogenetics and evolutionary analyses in R.”
#'  Bioinformatics, 35, 526-528, 2019). \cr
#'  - D. E. Rosen. Vicariant patterns and historical explanation in
#'  biogeography. Systematic Zoology, 27(2):159, 1978.
#'
#' @export
#' @rdname tmPDA
#'
#' @examples
#' genPDATree(n = 5)
genPDATree <- function(n){
  if(n < 2 || n%%1!=0){
    stop(paste("A tree must have at least 2 leaves, i.e., n>=2 and n must be",
               "an integer."))
  }
  phy <- ape::rtopology(n = n, rooted = TRUE,
                        tip.label = paste("t", sample.int(n,n), sep = ""))
  return(phy)
}
