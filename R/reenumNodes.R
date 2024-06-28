#' Function to modify the node enumeration in rooted binary trees
#'
#' \code{enum2cladewise} - Changes the node enumeration to cladewise
#' enumeration, i.e., starting from the root we follow the rule:\cr
#' Go to the left child; if that does not exist or was already visited go (up
#' again and) to the right child.\cr
#' The nodes in the rooted binary tree can be nearly arbitrarily enumerated
#' (distinct nodes should have distinct values and the
#' values should be positive, i.e., >0).
#'
#' @param phy A rooted binary tree of class \code{phylo}.
#' @param root Integer value (default = NULL) that should only be specified if
#' the root is known precisely (not necessary, but speeds up computation).
#'
#' @return \code{enum2cladewise}  A single tree of class \code{phylo} is
#' returned with cladewise node enumeration.
#'
#' @export
#' @rdname reenumNodes
#'
#' @examples
#' # Example with cladewise enumeration:
#' phy_alreadycladew <- list(edge = matrix(c(6,7, 7,8, 8,1, 8,2,
#'                                           7,9, 9,3, 9,4, 6,5),
#'                                  byrow = TRUE, ncol = 2),
#'                           tip.label = rep(" ",5), Nnode = 4)
#' attr(phy_alreadycladew, "class") <- "phylo"
#' enum2cladewise(phy_alreadycladew, root = 6)$edge
#' ape::plot.phylo(phy_alreadycladew)
#' # Example with other node enumeration:
#' phy_example <- list(edge = matrix(c(1,55, 55,12, 12,2, 12,10, 55,9,
#'                                     9,13, 9,60, 1,3),
#'                                   byrow = TRUE, ncol = 2),
#'                     tip.label = rep(" ",5), Nnode = 4,
#'                     edge.length = rep(1, 8))
#' attr(phy_example, "class") <- "phylo"
#' # The reenumeration works with and without specifying the root:
#' enum2cladewise(phy_example, root = 1)$edge
#' ape::plot.phylo(enum2cladewise(phy_example))
enum2cladewise <- function(phy, root = NULL){
  # 1.) Check requirements:
  n <- length(phy$tip.label)
  m <- phy$Nnode
  if ((m + n - 1) != nrow(phy$edge)) {
    stop(paste("The input must fulfill |V|-1=|E|=nrow(phy$edge)",
               "to be a tree.\n"))
  }
  node_labs <- sort(unique(as.vector(phy$edge)))
  if (length(node_labs)!= m+n) {
    stop("Not all |V|-many nodes have a unique enumeration.")
  }
  if (sum(node_labs %%1 !=0)>0 || sum(node_labs <=0)>0) {
    stop("Nodes must be labeled with integers >0.\n")
  }
  # and if necessary change node enumeration to 1,...,|V|
  if (!identical(seq(1,(m+n)), node_labs)) {
    label_change <- rep(NA, max(node_labs))
    label_change[node_labs] <- 1:(n+m)
    phy$edge <- matrix(label_change[phy$edge],
                               byrow = FALSE, ncol=2)
    if(!is.null(root)){ # if root is known
      root <- label_change[root]
    }
  }
  # 2.) Determine the starting point for the enumeration, i.e., the root.
  if(is.null(root)){ # if root is not known
    parent_nodes <- sort(unique(phy$edge[,1]))
    child_nodes <- sort(unique(phy$edge[,2]))
    root_candidate <- setdiff(parent_nodes, child_nodes)
    if(length(root_candidate)!=1){
      stop(paste("The root cannot be correctly determined. Root candidates:",
                 paste(root_candidate, collapse = ", ")))
    }
  } else { # if root already known
    root_candidate <- root
  }
  if(length(root_candidate)!=1 || root_candidate <=0 ||
     root_candidate %%1 !=0){
    stop(paste0("Unable to determine the root. ",
                "It must be an integer >0. Current value(s):",
               paste(root_candidate, collapse = ", ")))
  }
  # 3.) Create the cladewise order.
  phy$node.descs <- .getDescs(phy) # ---------------- remove afterwards
  cladew_nodes <- NULL
  cladew_edges <- NULL
  node_stack <- root_candidate
  edge_stack <- NA
  is_leaf <- rep(F, n+m)
  corresp_leafedge <- rep(NA, n+m)
  while(length(node_stack)>0){
    curr_node <- node_stack[1]
    cladew_nodes <- c(cladew_nodes, curr_node)
    cladew_edges <- c(cladew_edges, edge_stack[1])
    child_temp <- .getChildren(phy, curr_node, method = "alsoEdges")
    if(is.null(child_temp)){
      is_leaf[curr_node] <- T
      corresp_leafedge[curr_node] <- edge_stack[1]
    }
    node_stack <- c(child_temp[1,], node_stack[-1])
    edge_stack <- c(child_temp[2,], edge_stack[-1])
  }
  phy$node.descs <- NULL # --------------------------- removed here
  # 4.) Now, modify the node enumeration and edge order accordingly.
  # Create new enumeration: i-th position holds desired new number of i-th node
  new_enum <- rep(NA,n+m)
  internal_counter <- n+1
  leaf_counter <- 1
  for(i in cladew_nodes){
    if(is_leaf[i]){
      new_enum[i] <- leaf_counter
      leaf_counter <- leaf_counter +1
    } else {
      new_enum[i] <- internal_counter
      internal_counter <- internal_counter +1
    }
  }
  # Change node enumeration:
  phy$edge <- matrix(new_enum[phy$edge], byrow = FALSE, ncol=2)
  # Change edge order:
  phy$edge <- phy$edge[cladew_edges[-1],]
  phy$edge.length <- phy$edge.length[cladew_edges[-1]]
  # Change order-attribute:
  attr(phy, "order") <- "cladewise"
  return(phy)
}


# Method to access the children of a node (copied from package 'treeDbalance').
.getDescs <- function(tree){
  n <- length(tree$tip.label)
  m <- tree$Nnode
  descs <- rep(NA,m+n) #for the m+n-1 descendants (the root is no descendant)+NA
  edge_to_desc <- rep(NA,m+n) #for the m+n-1 edges to descendants (as above) +NA
  descs_index <- rep(1,m+n) #where the nodes' descendants start
  for(i in 1:nrow(tree$edge)){ # count number of descendants
    source <- tree$edge[i,1]
    descs_index <- descs_index + c(rep(0,source),rep(1,m+n-source))
  }
  cur_index <- rep(0,m+n)
  for(i in 1:nrow(tree$edge)){ # transfer every row (edge) to descs
    source <- tree$edge[i,1]
    descs[descs_index[source]+cur_index[source]] <- tree$edge[i,2]
    edge_to_desc[descs_index[source]+cur_index[source]] <- i
    cur_index[source] <- cur_index[source]+1
  }
  # Add one column to facilitate getChildren.
  return(cbind(rbind(descs,edge_to_desc,descs_index), c(NA,NA,m+n)))
}

# Method to access the children of a node (copied from package 'treeDbalance').
.getChildren <- function(tree, node, method = "onlyNodes"){
  if(!"node.descs" %in% attributes(tree)$names){
    comment(paste("This may take longer as the attribute 'node.descs'",
                  "does not exist and has to be calculated first."))
    tree$node.descs <- .getDescs(tree)
  }
  desc <- tree$node.descs

  if(method=="onlyNodes"){
    if(desc["descs_index",node]<desc["descs_index",node+1]){
      return(unname(desc["descs",
                    desc["descs_index",node]:(desc["descs_index",node+1]-1)]))
    } else {
      return(NULL)
    }

  }else if(method=="alsoEdges"){
    if(desc["descs_index",node]<desc["descs_index",node+1]){
      return(as.matrix(unname(desc[c("descs","edge_to_desc"),
                                   desc["descs_index",node]:(desc["descs_index",
                                                                  node+1]-1)])))
    } else {
      return(NULL)
    }

  }else{
    stop("Unknown method to getChildren of a node.")
  }
}
