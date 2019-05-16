# nodeTimes 
#
# Finds internal and external of a tree (ultrametric or non-ultrametric)

nodeTimes <- function(phy) {
	if(methods::is(phy$edge.length)[1] == "NULL") stop("no tree edge lengths")
    phy <- stats::reorder(phy)
    depBranches <- node.depth.edgelength(phy)
    all <- cbind(depBranches[phy$edge[, 1]], depBranches[phy$edge[ ,2]])
    return(max(depBranches) - all)
}
