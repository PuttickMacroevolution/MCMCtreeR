#' Estimate Cauchy Distribution for MCMCtree
#'
#' Estimate the offset and scale paramaters of a soft-tailed cauchy distribution and output trees for MCMCtree input
#' @param minAge vector of minimum age bounds for nodes matching order in monoGroups
#' @param maxAge vector of maximum age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param offset offset value for cauchy distribution (default = 50) (p in PAML manual page 49)
#' @param scale scale value for cauchy distribution (default = 1.5) (c in PAML manual page 49)
#' @param estimateScale logical specifying whether to estimate scale with a given shape value (default = TRUE)
#' @param minProb probability of left tail (minimum bound) - default to hard minimum (minProb=0)
#' @param maxProb probability of right tail (maximum bound. default = 0.975) 
#' @param plot logical specifying whether to plot to PDF
#' @param pdfOutput pdf output file name
#' @param writeMCMCtree logical whether to write tree in format that is compatible with MCMCTree to file
#' @param MCMCtreeName MCMCtree.output file name
#' @return list containing node estimates for each distribution
#' \itemize{
#'  \item{"parameters"}{ estimated parameters for each node}
#'  \item{"apePhy"}{ phylogeny in ape format with node labels showing node distributions}
#'  \item{"MCMCtree"}{ phylogeny in MCMCtree format}
#'  \item{"nodeLabels"}{ node labels in MCMCtreeR format}
#' }
#' @return If plot=TRUE plot of distributions in file 'pdfOutput' written to current working directory
#' @return If writeMCMCtree=TRUE tree in MCMCtree format in file "MCMCtreeName" written to current working directory
#' @export
#' @examples
#' @author Mark Puttick
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
#' estimateCauchy(minAge=minimumTimes, maxAge=maximumTimes, 
#' monoGroups=monophyleticGroups, offset=0.5, 
#' phy=apeTree, plot=FALSE)$MCMCtree

estimateCauchy <- function(minAge, maxAge, phy, monoGroups, scale=1.5, offset=50, estimateScale=TRUE,  minProb=0, maxProb=0.975, plot=FALSE, pdfOutput="cauchyPlot.pdf", writeMCMCtree=FALSE, MCMCtreeName="cauchyInput.tre") {
	
	lMin <- length(minAge)
	lMax <- length(maxAge)
	if(lMin != lMax) stop("length of ages do not match")
	if(length(minProb) < lMin) { minProb <- rep_len(minProb, lMin) ; print("warning - minProb parameter value recycled") }
	if(length(maxProb) < lMin) { maxProb <- rep_len(maxProb, lMin) ; print("warning - maxProb parameter value recycled") }
	if(length(offset) < lMin) { offset <- rep_len(offset, lMin) ; print("warning - offset parameter value recycled") }
	if(length(scale) < lMin) { scale <- rep_len(scale, lMin) ; print("warning - scale parameter value recycled") }
	if(length(estimateScale) < lMin) { estimateScale <- rep_len(estimateScale, lMin) ; print("warning - estimateScale argument recycled") }

	nodeFun <- function(x)	{
		scaleInt <- scale[x]
		locationInt <- minAge[x]
		offsetInt <- offset[x]
		estimateScaleInt <- estimateScale[x]
		maxProbInt <- maxProb[x]
		minProbInt <- minProb[x]
		if(estimateScaleInt == FALSE) {
			scaleInt <- scaleInt ; offsetInt <- offsetInt
		} else {
			p <- offsetInt
			cEsts <- c()
			cTest <- seq (0.001, 10, by=1e-3)
			for(u in 1:length(cTest)) cEsts[u] = 1 * (minAge[x] + p + cTest[u] * 1/tan((pi * (0.5 + 1/pi * atan(p/cTest[u])) * (1-maxProbInt)/(1 - minProbInt))))
			closest <- which(abs(cEsts-maxAge[x]) == min(abs(cEsts-maxAge[x])))
			upperEst <- cEsts[closest]
			scaleInt <- cTest[closest]
			} 
		if(minProbInt < 1e-7) minProbInt <- 1e-300
		nodeCon <- paste0("'L[", locationInt, "~", offsetInt, "~", scaleInt, "~",  minProbInt, "]'")	
		parameters <- c(locationInt, offsetInt, scaleInt, minProbInt)
		names(parameters) <- c("tL", "p", "c", "pL")
		return(list(nodeCon, parameters))	
		}
		
	out <- sapply(1:lMin, nodeFun)
	output <- c()
	prm <- matrix(unlist(out[2,]), ncol=4, byrow=T)
	rownames(prm) <- paste0("node_", 1:lMin)
	colnames(prm) <-  c("tL", "p", "c", "pL")
	output$parameters <- prm
	
	output$apePhy <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), T) 
	output$MCMCtree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), F) 

	if(writeMCMCtree == TRUE) {
		outputTree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), returnPhy=F) 
		utils::write.table(outputTree, paste0(MCMCtreeName), quote=F, row.names=F, col.names=F)
		}	
		
	if(plot == TRUE) {
		cat("warning - cauchy plots will be approximations!")
		if(length(list.files(pattern=paste0(pdfOutput))) != 0) {
			cat(paste0("warning - deleting and over-writing file ", pdfOutput))
			file.remove(paste0(pdfOutput))
			}
	 	grDevices::pdf(paste0(pdfOutput), family="Times")
		for(k in 1:dim(prm)[1]) {
			plotMCMCtree(prm[k,], method="cauchy",  paste0(rownames(prm)[k], " cauchy"), upperTime = max(maxAge)+1)
			}
		grDevices::dev.off()
		}	
	output$nodeLabels <- unlist(out[1,])	
	return(output)
}
