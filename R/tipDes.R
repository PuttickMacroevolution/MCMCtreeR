#' Find Descendent Tips From A Common Ancestor
#'
#' This function finds tip descendants from a common ancestor
#' @param  phy user tree in ape format
#' @param  node ancestral, crown node of a clade
#' @keywords node ancestor
#' @export
#' @examples
#' phy <- rcoal(10)
#' node <- 13
#' tipDes(phy, node)

tipDes <- function(phy, node) {
	
	keep <- which(phy$edge[,1] == node)
	newKeep <- phy$edge[keep,]
	tips <- which(newKeep[,2] <= Ntip(phy))
	if(sum(tips) > 0) node <- newKeep[-tips,2]
	if(sum(tips) == 0) node <- newKeep[,2]
	growth <- newKeep[,2]
	
	while(any(node > Ntip(phy)) == TRUE) {
		keep <- unlist(sapply(node, function(k) which(phy$edge[,1] == k)))
		newKeep <- phy$edge[keep,]
		growth <- c(growth, newKeep[,2])
		tips <- which(newKeep[,2] <= Ntip(phy))
		if(sum(tips) > 0) node <- newKeep[-tips,2]
		if(sum(tips) == 0) node <- newKeep[,2]
	}
	
	growth <- growth[which(growth <= Ntip(phy))]
	return(phy$tip.label[growth])
}
