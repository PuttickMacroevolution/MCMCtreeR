# nodeTimes 
#
# Finds internal and external of a tree (ultrmetric or non-ultrametric)

nodeTimes <- function(phy) {
    phy <- stats::reorder(phy)
    depBranches <- node.depth.edgelength(phy)
    all <- cbind(depBranches[phy$edge[, 1]], depBranches[phy$edge[ ,2]])
    return(max(depBranches) - all)
}
