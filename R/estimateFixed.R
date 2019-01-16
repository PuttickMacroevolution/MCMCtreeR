#' Fixed Age for PAML input
#'
#' Produce fixed age trees for PAML
#' @param minAge vector of fixed age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param writeMCMCTree logical whether to write tree in format that is compatible with mcmcTree to file
#' @param mcmcTreeName mcmcTree output file name
#' @return list containing node estimates for each distribution
#' \itemize{
#'  \item{"parameters"}{ estimated parameters for each node}
#'  \item{"apePhy"}{ phylogeny in ape format with node labels showing node distributions}
#'  \item{"mcmctree"}{ phylogeny in MCMCTree format}
#'  \item{"nodeLabels"}{ node labels in MCMCTreeR format}
#' }
#' @return If writeMCMCTree=TRUE tree in MCMCTree format in file "mcmcTreeName" written to current working directory
#' @examples
#' apeTree <- read.tree(text="((((human, (chimpanzee, bonobo)), gorilla), (orangutan, sumatran)), gibbon);")
#' monophyleticGroups <- list()
#' monophyleticGroups[[1]] <- c("human", "chimpanzee", "bonobo", "gorilla", "sumatran", "orangutan", "gibbon")
#' monophyleticGroups[[2]] <- c("human", "chimpanzee", "bonobo", "gorilla")
#' monophyleticGroups[[3]] <- c("human", "chimpanzee", "bonobo")
#' monophyleticGroups[[4]] <- c("sumatran", "orangutan")
#' minimumTimes <- c("nodeOne"=15, "nodeTwo"=6, "nodeThree"=8, "nodeFour"=13) / 10
#' maximumTimes <- c("nodeOne"=30, "nodeTwo"=12, "nodeThree"=12, "nodeFour"=20) / 10
#' estimateFixed(minAge=minimumTimes[1], monoGroups=monophyleticGroups[[1]], phy=apeTree)
#' @export

estimateFixed <- function(minAge, phy, monoGroups,  writeMCMCTree=FALSE, mcmcTreeName="estimateFixed.tre")	 {

	nodeCon <- paste0("'=", minAge, "'")
	output <- c()
	output$parameters <- c("fixed age" = minAge)
	output$apePhy <- nodeToPhy(phy, monoGroups, nodeCon, T) 
	output$mcmctree <- nodeToPhy(phy, monoGroups, nodeCon, F) 
	if(writeMCMCTree == T) 	write.table(output$mcmctree, paste0(mcmcTreeName), quote=F, row.names=F, col.names=F)	
	output$nodeLabels <- nodeCon
	return(output)
}