#' Generation of rooted binary trees under Aldous' beta splitting model
#'
#' \code{genAldousBetaTree} - Generates a rooted binary tree in
#' \code{phylo} format with the given number of \code{n} leaves under the
#' Aldous beta model.
#' The Aldous beta model is not a rate-based incremental evolutionary (tree)
#' construction and thus cannot generate edge lengths, only a topology.
#' Instead, the Aldous beta model works as follows: The idea is to start with
#' the root and the set of its descendant leaves, i.e., all \code{n} leaves.
#' Then, this set is partitioned into two subsets according to a density
#' function dependent on the parameter \code{beta}.
#' The two resulting subsets contain the leaves of the two maximal pending
#' subtrees of the root, respectively. The same procedure is then applied to the
#' root's children and their respective subsets, and so forth. \cr
#' Due to the restrictions of the \code{phylo} or \code{multiphylo} format
#' the number of leaves must be at least 2 since there must be at
#' least one edge.
#'
#' @author Sophie Kersting
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.
#' @param BETA Numeric value >=-2 which specifies how the leaf sets
#' are partitioned. For certain choices of \code{BETA} the Aldous beta model
#' coincides with known models:\cr
#' - \code{BETA} = 0: Yule model \cr
#' - \code{BETA} = -3/2: PDA model (all phylogenies equally probable) \cr
#' - \code{BETA} = -2: Caterpillar with \code{n} leaves
#'
#' @return \code{genAldousBetaTree} A single tree of class \code{phylo} is
#' returned.
#'
#' @references
#'  - D. Aldous. Probability Distributions on Cladograms. In Random Discrete
#'  Structures, pages 1–18. Springer New York, 1996.
#'
#'
#' @export
#' @rdname tmAldous
#'
#' @examples
#' genAldousBetaTree(n = 5, BETA = 1)
genAldousBetaTree <- function(n, BETA){
  if(n < 2 || n%%1!=0){
    stop(paste("A tree must have at least 2 leaves, i.e., n>=2 and n must be",
               "an integer."))
  }
  if(BETA < (-2)){
    stop(paste("The parameter BETA must be >=-2."))
  }
  # Use the auxiliary recursion to generate the edge matrix.
  edges <- .aldousBetaRecursion(BETA = BETA, n = n, nextFreeEnum = 1)
  # Create an environment that can be passed by reference and build the tree.
  phy <- list(edge = edges,
              tip.label = paste("t", sample.int(n,n), sep = ""),
              Nnode = n-1)
  attr(phy, "class") <- "phylo"
  # Last, change enumeration of nodes to fit 'cladewise' enumeration
  # and return the tree.
  return(enum2cladewise(phy, root = 1))
}

# Function for computing the probability of splitting n leaves into subsets
# of size i and n-i using the gamma function with parameter BETA.
.q_n <- function(i,BETA,n) {
  return(gamma(BETA+i+1)*gamma(BETA+n-i+1)/(gamma(i+1)*gamma(n-i+1)))
}

# Recursion for bipartitioning the set of leaves and creating the corresponding
# edge matrix.
.aldousBetaRecursion <- function(BETA, n, nextFreeEnum = 1) {
  if(n<2){
    stop("n must be >=2 for this recursion.")
  }
  # Partition into i and n-i leaves.
  if(BETA == (-2)){
    i <- 1
  } else {
    beta_probs <- sapply(1:(n-1), function(x) {.q_n(x,BETA,n)})
    not_NA <- !is.na(beta_probs)
    if(sum(not_NA)<1){
      i <- sample.int(n-1 ,1)
    } else {
      i <- sample(x=(1:(n-1))[not_NA] , size=1L , prob = beta_probs[not_NA])
    }
  }
  n_left <- i
  n_right <- n-i
  nextFreeEnum_left <- nextFreeEnum + 1
  nextFreeEnum_right <- nextFreeEnum + (2*n_left - 1) + 1
  # Initialize edge matrix.
  edges <- matrix(c(nextFreeEnum, nextFreeEnum_left,
                    nextFreeEnum, nextFreeEnum_right), byrow = TRUE, ncol = 2)
  # Add edges for left and right subtree if necessary.
  if(n_left > 1L){
    edges <- rbind(edges,.aldousBetaRecursion(BETA = BETA, n = n_left,
                                             nextFreeEnum = nextFreeEnum_left))
  }
  if(n_right > 1){
    edges <- rbind(edges,.aldousBetaRecursion(BETA = BETA, n = n_right,
                                             nextFreeEnum = nextFreeEnum_right))
  }
  return(edges)
}
