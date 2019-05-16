#' Find Descendent Tips From A Common Ancestor
#'
#' This function finds tip descendants from a common ancestor
#' @param  phy user tree in ape format
#' @param  node one or more nodes from the ape format that designate the crown monophyletic group
#' @details If a single node number is supplied the function returns a vector of tip labels. When more than one node number is supplied the function returns a list with each element a vector of tip labels for that node.
#' @keywords node ancestor
#' @author Mark Puttick
#' @export
#' @examples
#' set.seed(1029)
#' # one node
#' phy <- rcoal(10)
#' node <- 13
#' tipDes(phy, node)
#' ## multiple nodes
#' node <- c(13,14,15)
#' tipDes(phy, node)


tipDes <- function(phy, node) {
	
	list.out <- lapply(node, function(x.node) {
		keep <- which(phy$edge[,1] == x.node)
		newKeep <- phy$edge[keep,]
		tips <- which(newKeep[,2] <= Ntip(phy))
		if(sum(tips) > 0) x.node <- newKeep[-tips,2]
		if(sum(tips) == 0) x.node <- newKeep[,2]
		growth <- newKeep[,2]
	
		while(any(x.node > Ntip(phy)) == TRUE) {
			keep <- unlist(sapply(x.node, function(k) which(phy$edge[,1] == k)))
			newKeep <- phy$edge[keep,]
			growth <- c(growth, newKeep[,2])
			tips <- which(newKeep[,2] <= Ntip(phy))
			if(sum(tips) > 0) x.node <- newKeep[-tips,2]
			if(sum(tips) == 0) x.node <- newKeep[,2]
		}
	
		growth <- growth[which(growth <= Ntip(phy))]
		phy$tip.label[growth]
		}
	)
	if(length(node) == 1) {
		list.out <- list.out[[1]]
		} else {
		names(list.out) <- node
		}
	
	return(list.out)
}
