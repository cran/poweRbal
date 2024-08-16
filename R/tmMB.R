#' Generation of the maximally balanced tree
#'
#' \code{genMBTree} - Generates the rooted binary maximally balanced tree in
#' \code{phylo} format with the given number of \code{n} leaves.
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.\cr
#' Due to the restrictions of the \code{phylo} or \code{multiphylo} format,
#' the number of leaves must be at least 2 since there must be at
#' least one edge.
#'
#' @return \code{genMBTree} A single tree of class \code{phylo} is returned.
#'
#' @references
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic models.
#'  arXiv:2406.05185, 2024.
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic
#'  models: Supplementary material. https://tinyurl.com/278cwdh8, 2024.
#'
#' @export
#' @rdname tmMB
#'
#' @examples
#' genMBTree(n = 6)
genMBTree <- function(n){
  if(n < 2 || n%%1!=0){
    stop(paste("A tree must have at least 2 leaves, i.e., n>=2 and n must be",
               "an integer."))
  }
  maxbal_n <- paste0(.auxMB_n(n),";")
  return(phytools::read.newick(text=maxbal_n))
}

# Auxiliary function to produce the Newick format of the MB-tree.
.auxMB_n <- function(n){
  if(n==1){
    return(" ")
  }else if(n==2){
    return("( , )")
  }else{
    return(paste0("(",.auxMB_n(ceiling(n/2)),","
                     ,.auxMB_n(floor(n/2)),")"))
  }
}
