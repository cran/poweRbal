#' Generation of rooted binary trees under Ford's alpha model
#'
#' \code{genFordsAlphaTree} - Generates a rooted binary tree in
#' \code{phylo} format with the given number of \code{n} leaves under Ford's
#' alpha model.
#' Ford's alpha model is not a rate-based evolutionary (tree)
#' construction and thus cannot generate edge lengths, only a topology.
#' Instead, it works as follows: The idea is to start with a cherry and
#' incrementally increase the size of the tree by adding a new leaf with
#' a leaf edge to any edge (inner or leaf edge), one at a time.
#' Given a tree with i leaves, then each of the i-1 innner edges (includes
#' an additional root edge) is chosen with probability
#' \code{ALPHA}/(i-\code{ALPHA}).
#' Each of the i leaf edges is chosen with probability
#' (1-\code{ALPHA})/(i-\code{ALPHA}). \cr
#' Due to the restrictions of the \code{phylo} or \code{multiphylo} format
#' the number of leaves must be at least 2 since there must be at
#' least one edge.
#'
#' @author Sophie Kersting
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.
#' @param ALPHA Numeric value >=0 and <=1 which specifies the probabilites of
#' picking an inner or a leaf edge. For certain choices of
#' \code{ALPHA} Ford's alpha model coincides with known models:\cr
#' - \code{ALPHA} = 0: Yule model \cr
#' - \code{ALPHA} = 1/2: PDA model (all phylogenies equally probable) \cr
#' - \code{ALPHA} = 1: Caterpillar with \code{n} leaves
#'
#' @return \code{genFordsAlphaTree} A single tree of class \code{phylo} is
#' returned.
#'
#' @references
#'  - D. J. Ford. Probabilities on cladograms: introduction to the alpha
#'  model, 2005. \cr
#'  - G. Kaur, K. P. Choi, and T. Wu. Distributions of cherries and pitchforks
#'  for the Ford model. Theoretical Population Biology, 149:27–38, 2023.
#'
#' @export
#' @rdname tmFord
#'
#' @examples
#' genFordsAlphaTree(n = 5, ALPHA = 0.3)
genFordsAlphaTree <- function(n, ALPHA){
  if(n < 2 || n%%1!=0){
    stop(paste("A tree must have at least 2 leaves, i.e., n>=2 and n must be",
               "an integer."))
  }

  # Create the edge matrix -----------------------------------------------------
  m <- matrix(rep(NA,(2*n-1)*2), nrow = 2*n-1, ncol = 2)
  # Initialize matrix with cherry and root edge (first row) --------------------
  m[1,] <- c(1,2)
  m[2,] <- c(2,3)
  m[3,] <- c(2,4)
  # Initialize matrix with cherry and root edge --------------------------------
  is_leaf_edge <- c(FALSE, TRUE, TRUE, rep(NA, 2*n-4))
  # Do (n-2) speciation steps --------------------------------------------------
  if(n > 2){
    for(i in 2:(n-1)) {
      # Choose edge for leaf insertion -----------------------------------------
      numb_edges <- 2*i -1
      curr_probs <- rep(NA, numb_edges)
      curr_probs[is_leaf_edge[1:numb_edges]] <- (1-ALPHA)/(i-ALPHA)
      curr_probs[!is_leaf_edge[1:numb_edges]] <- ALPHA/(i-ALPHA)
      edge_index <- sample(1:numb_edges, size = 1, replace = F,
                           prob = curr_probs)
      # Numbers for incident and new nodes -------------------------------------
      edge_parent <- m[edge_index,1]
      edge_child <- m[edge_index,2]
      new_inner <- 2*i+1
      new_leaf <- 2*i+2
      # Edit edge matrix and types of edges ------------------------------------
      m[edge_index,] <- c(edge_parent, new_inner)
      m[numb_edges+1,] <- c(new_inner, edge_child)
      m[numb_edges+2,] <- c(new_inner, new_leaf)
      is_leaf_edge[numb_edges+1] <- is_leaf_edge[edge_index]
      is_leaf_edge[numb_edges+2] <- TRUE
      is_leaf_edge[edge_index] <- FALSE
    }
  }
  # Delete root edge and save the root -----------------------------------------
  curr_root <- m[1,2]
  m <- m[-1,]
  # Create the phylo object and enumerate cladewise ----------------------------
  phy <- list(edge = m, tip.label = paste("t", sample.int(n,n), sep = ""),
              Nnode = as.integer(n-1))
  attr(phy, "class") <- "phylo"
  return(enum2cladewise(phy, root = curr_root))
}
