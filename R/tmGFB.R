#' Generation of the greedy from the bottom tree
#'
#' \code{genGFBTree} - Generates the rooted binary greedy from the bottom  tree
#' in \code{phylo} format with the given number of \code{n} leaves. \cr
#' Due to the restrictions of the \code{phylo} or \code{multiphylo} format
#' the number of leaves must be at least 2 since there must be at
#' least one edge.
#'
#' @author Sophie Kersting
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.
#'
#' @return \code{genGFBTree} A single tree of class \code{phylo} is returned.
#'
#' @references
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic models.
#'  arXiv:2406.05185, 2024.
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic
#'  models: Supplementary material. https://tinyurl.com/278cwdh8, 2024.
#'
#' @export
#' @rdname tmGFB
#'
#' @examples
#' genGFBTree(n = 6)
genGFBTree <- function(n){
  if(n < 2 || n%%1!=0){
    stop(paste("A tree must have at least 2 leaves, i.e., n>=2 and n must be",
               "an integer."))
  }
  maxbal_n <- paste0(.auxGFB_n(n),";")
  return(phytools::read.newick(text=maxbal_n))
}

# Auxiliary function to produce the Newick format of the GFB-tree.
.auxGFB_n <- function(n){
  if(n==1){
    return(" ")
  }else if(n==2){
    return("( , )")
  }else if(n==3){
    return("(( , ), )")
  }else{
    full_lvls <- floor(log(n, base = 2))
    if((n-2^full_lvls)<=2^(full_lvls-1)){
      n_r <- 2^(full_lvls-1)
      n_l <- n-n_r
    }else{
      n_l <- 2^full_lvls
      n_r <- n-n_l
    }
    return(paste0("(",.auxGFB_n(n_l),","
                  ,.auxGFB_n(n_r),")"))
  }
}
