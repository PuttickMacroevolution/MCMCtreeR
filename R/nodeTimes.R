#' nodeTimes 
#'
#' Finds internal and external of a tree (ultrmetric or non-ultrametric)
#' @param phy user tree in ape format
#' @export

nodeTimes <- function(phy) {
    phy <- reorder(phy)
    depBranches <- node.depth.edgelength(phy)
    all <- cbind(depBranches[phy$edge[, 1]], depBranches[phy$edge[ ,2]])
    return(max(depBranches) - all)
}