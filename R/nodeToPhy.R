#' Add Node Constraints to Tree in MCMCtree Format
#'
#' Produce tree with node labels in MCMCtree format
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param nodeCon node distributions in MCMCtreeR format
#' @param returnPhy logical specifying whether to return phy to console or write MCMCtree for input (default = TRUE)
#' @return If returnPhy=TRUE phylogeny with node labels in ape format
#' @return If returnPhy=FALSE phylogeny with node labels in MCMCtree format
#' @author Mark Puttick

nodeToPhy <- function(phy, monoGroups, nodeCon, returnPhy=TRUE) {
	if(methods::is(monoGroups)[1] != "list") monoGroups <- list(monoGroups)
	for(x in 1:length(nodeCon))	{
		# if(length(length(nodeCon)) == 1) phy$node.labels <- rep("", Nnode(phy))
		nodeLocation <- which((Ntip(phy) + 1) : (Ntip(phy) + Nnode(phy)) == ape::getMRCA(phy, monoGroups[[x]]))
		phy$node.labels[nodeLocation] <- nodeCon[x]
		}
	if(returnPhy == T) return(phy)
	if(returnPhy == F) {
		s <- write.tree(phy)	
		s <- gsub(")NA", ")", s)
		s <- gsub("\\[", "(", s)
		s <- gsub("\\]", ")", s)
		s <- gsub("~", ",", s)
		names(s) <- NULL
		treeFile <- data.frame(cbind(paste(Ntip(phy), "1"), s, "//end of file"))
		colnames(treeFile) <- rep("", 3)
		return(treeFile)
		}
	}
