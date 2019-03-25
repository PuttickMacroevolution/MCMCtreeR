#' Estimate Upper Age for MCMCtree analysis
#'
#' Estimate the upper age distribution and output trees for MCMCtree input
#' @param maxAge vector of maximum age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param rightTail probability of right tail (maximum bound default = 0.025) 
#' @param writeMCMCtree logical whether to write tree in format that is compatible with MCMCTree to file
#' @param MCMCtreeName MCMCtree.output file name
#' @return list containing node estimates for each distribution
#' \itemize{
#'  \item{"parameters"}{ estimated parameters for each node}
#'  \item{"apePhy"}{ phylogeny in \pkg{APE} format with node labels showing node distributions}
#'  \item{"MCMCtree"}{ phylogeny in MCMCtree format}
#'  \item{"nodeLabels"}{ node labels in MCMCtreeR format}
#' }
#' @return If plot=TRUE plot of distributions in file 'pdfOutput' written to current working directory
#' @return If writeMCMCtree=TRUE tree in MCMCtree format in file "MCMCtreeName" written to current working directory
#' @export
#' @examples
#' data(apeData)
#' attach(apeData)
#' monophyleticGroups <- list()
#' monophyleticGroups[[1]] <- c("human", "chimpanzee", "bonobo", 
#' "gorilla", "sumatran", "orangutan", "gibbon")
#' getMRCA(apeTree, c("human", "chimpanzee", "bonobo", "gorilla"))
#' monophyleticGroups[[2]] <-  tipDes(apeTree, 10)
#' monophyleticGroups[[3]] <- tipDes(apeTree, 11)
#' monophyleticGroups[[4]] <- c("sumatran", "orangutan")
#' minimumTimes <- c("nodeOne"=15, "nodeTwo"=6,
#' "nodeThree"=8, "nodeFour"=13) / 10
#' maximumTimes <- c("nodeOne" = 30, "nodeTwo" = 12,
#' "nodeThree"=12, "nodeFour" = 20) / 10
#' estimateUpper(maxAge=maximumTimes, monoGroups=monophyleticGroups,
#' rightTail=0.025, phy=apeTree)

estimateUpper <- function(maxAge, rightTail=0.025,  phy, monoGroups, writeMCMCtree=FALSE, MCMCtreeName="estimateUpper.tre") {
    lMax <- length(maxAge)
    if(length(rightTail) < lMax) { rightTail <- rep_len(rightTail, lMax) ; print("warning - maxProb parameter value recycled") }
   	if(any(rightTail < 1e-7)) rightTail[which(rightTail < 1e-7)] <- 1e-300
    nodeFun <- function(x) {
        nodeCon <- paste0("'U[", maxAge[x], "~", rightTail[x], "]'")
        parameters <- c(maxAge[x], rightTail[x])
        return(list(nodeCon, parameters))
	    }
    out <- sapply(1:lMax, nodeFun)
    output <- c()
    prm <- matrix(unlist(out[2,]), ncol=2, byrow=TRUE)
    rownames(prm) <- paste0("node_", 1:lMax)
    colnames(prm) <-  c("tU", "pR")
    output$apePhy <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), TRUE)
    output$MCMCtree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), FALSE)
    output$parameters <- prm
	if(writeMCMCtree == TRUE) 	utils::write.table(output$MCMCtree, paste0(MCMCtreeName), quote=FALSE, row.names=FALSE, col.names=FALSE)	
	output$nodeLabels <- unlist(out[1,])
	return(output)
}
