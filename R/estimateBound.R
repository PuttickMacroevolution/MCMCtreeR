#' Estimate a Uniform Distribution for MCMCtree
#'
#' Estimate the paramaters of a soft-bounded uniform distribution and output trees for MCMCtree input
#' @param minAge vector of minimum age bounds for nodes matching order in monoGroups
#' @param maxAge vector of maximum age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param minProb probability of left tail (minimum bound) - default to hard minimum (minProb=0)
#' @param rightTail probability of right tail (maximum bound default = 0.975) 
#' @param plot logical specifying whether to plot to PDF
#' @param pdfOutput pdf output file name
#' @param writeMCMCtree logical whether to write tree in format that is compatible with MCMCTree to file
#' @param MCMCtreeName MCMCtree.output file name
#' @return list containing node estimates for each distribution
#' \itemize{
#'  \item{"parameters"}{ estimated parameters for each node}
#'  \item{"apePhy"}{ phylogeny in ape format with node labels showing node distributions}
#'  \item{"MCMCtree"}{ phylogeny in MCMCtreeR format}
#'  \item{"nodeLabels"}{ node labels in MCMCtreeR format}
#' }
#' @return If plot=TRUE plot of distributions in file 'pdfOutput' written to current working directory
#' @return If writeMCMCtree=TRUE tree in MCMCtree format in file "MCMCtreeName" written to current working directory
#' @export
#' @author Mark Puttick
#' @examples
#' data(apeData)
#' attach(apeData)
#' ## extract taxon descending from calibrated nodes 8, 10, 11, 13
#' ## these nodes can be visualised using plot.phylo
#' ## and nodelabels from ape
#' monophyleticGroups <- tipDes(apeData$apeTree, c(8,10,11,13))
#' minimumTimes <- c("nodeOne"=15, "nodeTwo"=6,
#' "nodeThree"=8, "nodeFour"=13) / 10
#' maximumTimes <- c("nodeOne" = 30, "nodeTwo" = 12,
#' "nodeThree"=12, "nodeFour" = 20) / 10
#' estimateBound(minAge=minimumTimes, maxAge=maximumTimes, 
#' monoGroups=monophyleticGroups, phy=apeTree, plot=FALSE)$MCMCtree

estimateBound <- function(minAge, maxAge, minProb=0.025, rightTail=0.025, phy, monoGroups, writeMCMCtree=FALSE, plot=TRUE, MCMCtreeName="bound.tre", pdfOutput="uniformPlot.pdf")	 {

	lMin <- length(minAge)
	lMax <- length(maxAge)
	if(lMin != lMax) stop("length of ages do not match")
	if(length(minProb) < lMin) { minProb <- rep_len(minProb, lMin) ; print("warning - minProb parameter value recycled") }
	if(length(rightTail) < lMin) { rightTail <- rep_len(rightTail, lMin) ; print("warning - maxProb parameter value recycled") }
	if(any(rightTail < 1e-7)) rightTail[which(rightTail < 1e-7)] <- 1e-300
	if(any(minProb < 1e-7)) minProb[which(minProb < 1e-7)] <- 1e-300
	
	nodeFun <- function(x){
		nodeCon <- paste0("'B[", minAge[x], "~", maxAge[x], "~", minProb[x], "~", rightTail[x], "]'")
		parameters <- c(minAge[x], maxAge[x], minProb[x], rightTail[x])
		return(list(nodeCon, parameters))
		}
		
	out <- sapply(1:lMin, nodeFun)
	output <- c()
	prm <- matrix(unlist(out[2,]), ncol=4, byrow=T)
	rownames(prm) <- paste0("node_", 1:lMin)
	colnames(prm) <-  c("tL", "tU", "pL", "pU")
	output$apePhy <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), T) 
	output$MCMCtree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), F) 
    output$parameters <- prm

		
	if(writeMCMCtree == T) {
		outputTree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), returnPhy=F) 
		utils::write.table(outputTree, paste0(MCMCtreeName), quote=F, row.names=F, col.names=F)
		}
	if(plot == T) {
		if(length(list.files(pattern=paste0(pdfOutput))) != 0) {
			cat(paste0("warning - deleting and over-writing file ", pdfOutput))
			file.remove(paste0(pdfOutput))
			}
	 	grDevices::pdf(paste0(pdfOutput), family="Times")
		for(k in 1:dim(prm)[1]) {
			plotMCMCtree(prm[k,], method="bound",  paste0(rownames(prm)[k], " bound"), upperTime = maxAge[k] +1)
			}
		grDevices::dev.off()
		}	
	
	output$nodeLabels <- unlist(out[1,])	
	return(output)	
	
}
