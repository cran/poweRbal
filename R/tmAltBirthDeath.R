#' Generation of rooted binary trees under the alternative birth-death model
#'
#' \code{genAltBirthDeathTree} - Generates a rooted binary tree in \code{phylo}
#' format with the given number of \code{n} leaves under the alternative
#' birth-death model.
#' In the alternative birth-death process all species have the same speciation
#' \code{BIRTHRATE} and extinction rates \code{DEATHRATE}. Extinct species
#' remain as fossils inside the tree with zero speciation and extinction
#' rates.
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.\cr
#' Due to the restrictions of the \code{phylo} or \code{multiphylo} format,
#' the number of leaves must be at least 2 since there must be at
#' least one edge.
#' @param BIRTHRATE Positive numeric value (default = 1) which specifies the
#' rate at which the speciation events occur.
#' @param DEATHRATE Positive numeric value (default = 0) which specifies the
#' rate at which the extinction events occur.
#' @param TRIES Integer value (default = 5) that specifies
#' the number of attempts to generate a tree with \code{n} leaves.
#'
#' @return \code{genAltBirthDeathTree} A single tree of class \code{phylo} is
#' returned.
#'
#' @references
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic models.
#'  arXiv:2406.05185, 2024.
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic
#'  models: Supplementary material. https://tinyurl.com/278cwdh8, 2024.
#'
#' @export
#' @rdname tmAltBirthDeath
#' @examples
#' genAltBirthDeathTree(n = 7, DEATHRATE = 1)
genAltBirthDeathTree <- function(n, BIRTHRATE = 1, DEATHRATE = 0,
                                 TRIES = 5){
  if(n < 2 || n%%1!=0){
    stop(paste("A tree must have at least 2 leaves, i.e., n>=2 and n must be",
               "an integer."))
  }
  if(BIRTHRATE<=0){
    stop(paste("The speciation rate must be >0."))
  }
  if(DEATHRATE<0){
    stop(paste("The extinction rate must be >=0."))
  }
  phy <- NULL
  i <- 1
  while(i<=TRIES && is.null(phy)){
    i <- i+1
    # Create the edge matrix
    m <- matrix(rep(NA,(2*n-2)*2), nrow = 2*n-2, ncol = 2)
    # Initialize vector for current leaves and their rates
    curr_leaves <- 1 # Vector of the current leaves (start: only one node)
    birth_rates <- c(BIRTHRATE, rep(NA, 2*n-2)) # Rates of all nodes
    death_rates <- c(DEATHRATE, rep(NA, 2*n-2))
    free_numbers <- c(FALSE,rep(TRUE, 2*n-2))
    free_rows <- rep(TRUE, 2*n-2)
    anc_edge <- c(0, rep(NA, 2*n-2)) # Which edge leads from parent to node?
    desc_edges <- matrix(rep(NA,(2*n-1)*2), nrow = 2*n-1, ncol = 2)
    # Do speciation and extinction steps as long as necessary
    while(length(curr_leaves) < n &&
          sum(birth_rates[curr_leaves]>0)>=1){
      # Determine event type and the affected leaf
      is_speciation_event <- sample(c(TRUE, FALSE), size = 1, replace = F,
                                    prob = c(sum(birth_rates[curr_leaves]),
                                             sum(death_rates[curr_leaves])))
      if(is_speciation_event){ # Speciation event
        leaf_index <- sample(1:length(curr_leaves), size = 1, replace = F,
                             prob = birth_rates[curr_leaves]) # Choose leaf
        # New numbers for children
        child_num <- which(free_numbers)[c(1,2)]
        free_numbers[child_num] <- FALSE
        # Fill out matrix row by row (edges: parent->child)
        row_num <- which(free_rows)[c(1,2)]
        free_rows[row_num] <- FALSE
        m[row_num[1],] <- c(curr_leaves[leaf_index], child_num[1])
        m[row_num[2],] <- c(curr_leaves[leaf_index], child_num[2])
        anc_edge[child_num] <- row_num
        desc_edges[curr_leaves[leaf_index],] <- row_num
        # Remove parent rate and insert children's rates.
        birth_rates[curr_leaves[leaf_index]] <- NA
        death_rates[curr_leaves[leaf_index]] <- NA
        birth_rates[child_num] <- BIRTHRATE
        death_rates[child_num] <- DEATHRATE
        # Remove parent from current leaves and insert children.
        curr_leaves <- c(curr_leaves[-leaf_index], child_num[1], child_num[2])
      } else { # Extinction event
        leaf_index <- sample(1:length(curr_leaves), size = 1, replace = F,
                             prob = death_rates[curr_leaves]) # Choose leaf
        birth_rates[curr_leaves[leaf_index]] <- 0
        death_rates[curr_leaves[leaf_index]] <- 0
      }
    }
    if(length(curr_leaves) == n){ # If successful
      # Create the phylo object and enumerate cladewise
      phy <- list(edge = m, tip.label = paste("t", sample.int(n,n), sep = ""),
                  Nnode = as.integer(n-1))
      attr(phy, "class") <- "phylo"
      phy <- enum2cladewise(phy, root = 1)
    }
  }
  return(phy)
}

