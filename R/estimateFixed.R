#' Fixed Age for MCMCtree analysis input
#'
#' Produce fixed age trees for MCMCtree analysis
#' @param minAge vector of fixed age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param writeMCMCtree logical whether to write tree in format that is compatible with MCMCTree to file
#' @param MCMCtreeName MCMCtree.output file name
#' @return list containing node estimates for each distribution
#' \itemize{
#'  \item{"parameters"}{ estimated parameters for each node}
#'  \item{"apePhy"}{ phylogeny in ape format with node labels showing node distributions}
#'  \item{"MCMCtree"}{ phylogeny in MCMCtreeR format}
#'  \item{"nodeLabels"}{ node labels in MCMCtreeR format}
#' }
#' @return If writeMCMCtree=TRUE tree in MCMCtree format in file "MCMCtreeName" written to current working directory
#' @examples
#' data(apeData)
#' attach(apeData)
#' monophyleticGroups <- list()
#' monophyleticGroups[[1]] <- c("human", "chimpanzee",
#' "bonobo", "gorilla", "sumatran", "orangutan", "gibbon")
#' monophyleticGroups[[2]] <- c("human", "chimpanzee",
#' "bonobo", "gorilla")
#' monophyleticGroups[[3]] <- c("human", "chimpanzee", "bonobo")
#' monophyleticGroups[[4]] <- c("sumatran", "orangutan")
#' minimumTimes <- c("nodeOne"=15, "nodeTwo"=6,
#' "nodeThree"=8, "nodeFour"=13) / 10
#' maximumTimes <- c("nodeOne"=30, "nodeTwo"=12,
#' "nodeThree"=12, "nodeFour"=20) / 10
#' estimateFixed(minAge=minimumTimes[1],
#' monoGroups=monophyleticGroups[[1]], phy=apeTree)
#' @export

estimateFixed <- function(minAge, phy, monoGroups, writeMCMCtree=FALSE, MCMCtreeName="estimateFixed.tre") {
	nodeCon <- paste0("'=", minAge, "'")
	output <- c()
	output$parameters <- c("fixed age" = minAge)
	output$apePhy <- nodeToPhy(phy, monoGroups, nodeCon, TRUE) 
	output$MCMCtree <- nodeToPhy(phy, monoGroups, nodeCon, FALSE) 
	if(writeMCMCtree == TRUE) 	utils::write.table(output$MCMCtree, paste0(MCMCtreeName), quote=FALSE, row.names=FALSE, col.names=FALSE)	
	output$nodeLabels <- nodeCon
	return(output)
}
