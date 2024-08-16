#' Generation of rooted binary trees under the density model
#'
#' \code{genDensityTree} - Generates a rooted binary tree in \code{phylo}
#' format with the given number of \code{n} leaves under the
#' density-dependent model.
#' In the density-dependent tree generation process all species have the same
#' speciation \code{BIRTHRATE}, but the extinction rates depend on the
#' number of species (it increases linearly with the number of co-existing
#' lineages until an equilibrium number is reached at which speciation and
#' extinction rates are equal).
#' Extinct species are removed from the tree, i.e., the generated tree contains
#' only species living at the present.
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.\cr
#' Due to the restrictions of the \code{phylo} or \code{multiphylo} format,
#' the number of leaves must be at least 2 since there must be at
#' least one edge.
#' @param EQUILIB Integer value that specifies the equilibrium number.
#' @param BIRTHRATE Positive numeric value (default = 1) which specifies the
#' rate at which the speciation events occur.
#' @param TRIES Integer value (default = 5) that specifies
#' the number of attempts to generate a tree with \code{n} leaves.
#' @param TIMEperTRY Numeric value (default = 0.01) that specifies the maximum
#' amount of time (in seconds) invested per try.
#'
#' @return \code{genDensityTree} A single tree of class \code{phylo} is
#' returned.
#'
#' @references
#'  - P. H. Harvey, R. M. May, and S. Nee. Phylogenies without fossils.
#'  Evolution, 48(3):523–529, 1994.
#'
#'
#' @export
#' @rdname tmDensity
#' @examples
#' genDensityTree(n = 5, EQUILIB = 6)
genDensityTree <- function(n, BIRTHRATE = 1, EQUILIB,
                           TRIES = 5, TIMEperTRY = 0.01){
  if(n < 2 || n%%1!=0){
    stop(paste("A tree must have at least 2 leaves, i.e., n>=2 and n must be",
               "an integer."))
  }
  if(EQUILIB < 1 || EQUILIB %%1!=0){
    stop(paste("The equilibrium number must be an integer >0."))
  }
  if(BIRTHRATE<=0){
    stop(paste("The speciation rate must be >0."))
  }
  phy <- NULL
  i <- 1
  while(i<=TRIES && is.null(phy)){
    starting_time <- Sys.time()
    i <- i+1
    # Create the edge matrix ---------------------------------------------------
    m <- matrix(rep(NA,(2*n-2)*2), nrow = 2*n-2, ncol = 2)
    # Initialize vector for current leaves and their rates ---------------------
    curr_leaves <- 1 # Vector of the current leaves (start: only one node)
    birth_rates <- rep(BIRTHRATE, 2*n-1) # Birth rates remain constant
    free_numbers <- c(FALSE,rep(TRUE, 2*n-2))
    free_rows <- rep(TRUE, 2*n-2)
    anc_edge <- c(0, rep(NA, 2*n-2)) # Which edge leads from parent to node?
    desc_edges <- matrix(rep(NA,(2*n-1)*2), nrow = 2*n-1, ncol = 2)
    # Do speciation and extinction steps as long as necessary
    while(length(curr_leaves) < n && length(curr_leaves)>0 &&
          Sys.time()-starting_time < TIMEperTRY){
      # Death rates are density dependent and have to be updated
      if(length(curr_leaves) < EQUILIB){
        death_rates <- rep(BIRTHRATE*(length(curr_leaves)/EQUILIB), 2*n-1)
      } else {
        death_rates <- rep(BIRTHRATE, 2*n-1)
      }
      # Determine event type and the affected leaf
      is_speciation_event <- sample(c(TRUE, FALSE), size = 1, replace = F,
                                    prob = c(sum(birth_rates[curr_leaves]),
                                             sum(death_rates[curr_leaves])))
      leaf_index <- sample.int(length(curr_leaves),1) # Choose leaf
      if(is_speciation_event){ # Speciation event ------------------------------
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
        # Remove parent from current leaves and insert children.
        curr_leaves <- c(curr_leaves[-leaf_index], child_num[1], child_num[2])
      } else { # Extinction event  ---------------------------------------------
        if(length(curr_leaves) == 1){
          curr_leaves <- numeric(0)
        } else {
          # Numbers of parent and grandparent (with edges)
          node_leaf <- curr_leaves[leaf_index]
          edge_to_parent <- anc_edge[node_leaf]
          node_parent <- m[edge_to_parent,1]
          edge_to_grandpar <- anc_edge[node_parent]
          edge_to_sis <- desc_edges[node_parent,][which(!desc_edges[node_parent,]==edge_to_parent)]
          node_sis <- m[edge_to_sis,2]
          # Remove leaf, leaf edge, sister edge and parent
          if(edge_to_grandpar == 0){ # If the parent is the root.
            m[edge_to_sis, ] <- NA;    free_rows[edge_to_sis] <- TRUE
            m[edge_to_parent, ] <- NA; free_rows[edge_to_parent] <- TRUE
            anc_edge[c(node_leaf, node_sis, node_parent)] <- c(NA, 0, NA)
            desc_edges[node_parent, ] <- NA
            free_numbers[node_leaf] <- TRUE; free_numbers[node_parent] <- TRUE
            # Remove parent from current leaves and insert children.
            curr_leaves <- curr_leaves[-leaf_index]
          } else {
            m[edge_to_grandpar, 2] <- node_sis
            m[edge_to_sis, ] <- NA; free_rows[edge_to_sis] <- TRUE
            m[edge_to_parent, ] <- NA; free_rows[edge_to_parent] <- TRUE
            anc_edge[c(node_leaf, node_sis, node_parent)] <- c(NA,
                                                               edge_to_grandpar,
                                                               NA)
            desc_edges[node_parent, ] <- NA
            free_numbers[node_leaf] <- TRUE; free_numbers[node_parent] <- TRUE
            # Remove parent from current leaves and insert children.
            curr_leaves <- curr_leaves[-leaf_index]
          }
        }
      }
    }
    if(length(curr_leaves) == n){ # If successful
      # Create the phylo object and enumerate cladewise
      phy <- list(edge = m, tip.label = paste("t", sample.int(n,n), sep = ""),
                  Nnode = as.integer(n-1))
      attr(phy, "class") <- "phylo"
      phy <- enum2cladewise(phy, root = which(anc_edge==0))
    }
  }
  if(is.null(phy)){
    comment("Not able to generate tree. Returning NULL instead.")
  }
  return(phy)
}

