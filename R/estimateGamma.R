#' Estimate Gamma Distribution for MCMCtree analysis
#'
#' Estimate the shape and rate paramaters of Gamma distribution and output trees for MCMCtree input
#' @param minAge vector of minimum age bounds for nodes matching order in monoGroups
#' @param maxAge vector of maximum age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param alpha alpha value for gamma distribution (default = 188)
#' @param beta beta value for gamma distribution (default = 2690)
#' @param offset distance of mean value from minimum bound
#' @param estimateAlpha logical specifying whether to estimate alpha with a given beta value (default = TRUE)
#' @param estimateBeta logical specifying whether to estimate beta with a given alpha value (default = FALSE)
#' @param plot logical specifying whether to plot to PDF
#' @param pdfOutput pdf output file name
#' @param writeMCMCtree logical whether to write tree in format that is compatible with MCMCTree to file
#' @param MCMCtreeName MCMCtree.output file name
#' @return list containing node estimates for each distribution
#' \itemize{
#'  \item{"parameters"}{ estimated parameters for each node}
#'  \item{"apePhy"}{ phylogeny in \pkg{APE} format with node labels showing node distributions}
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
#' minimumTimes <- c("8"=15, "10"=6,
#' "11"=8, "13"=13) / 10
#' maximumTimes <- c("8" = 30, "10" = 12,
#' "11"=12, "13" = 20) / 10
#' gamma.nodes <- estimateGamma(minAge=minimumTimes, maxAge=maximumTimes, 
#' monoGroups=monophyleticGroups, alpha=188, beta=2690, 
#' offset=0.1, phy=apeTree, plot=FALSE)
#' gamma.nodes

estimateGamma <- function(minAge, maxAge, phy, monoGroups, alpha=188, beta=2690, offset=0.1, estimateAlpha=TRUE, estimateBeta=FALSE, plot=FALSE, pdfOutput="gammaPlot.pdf", writeMCMCtree=FALSE, MCMCtreeName="gammaInput.tre") {
	
	lMin <- length(minAge)
	lMax <- length(maxAge)
	if(lMin != lMax) stop("length of ages do not match")
	if(length(alpha) < lMin) { alpha <- rep_len(alpha, lMin) ; print("warning - alpha parameter value recycled") }
	if(length(beta) < lMin) { beta <- rep_len(beta, lMin) ; print("warning - beta parameter value recycled") }
	if(length(offset) < lMin) { offset <- rep_len(offset, lMin) ; print("warning - offset parameter value recycled") }
	if(length(estimateAlpha) < lMin) { estimateAlpha <- rep_len(estimateAlpha, lMin) ; print("warning - estimateAlpha argument recycled") }
	if(length(estimateBeta) < lMin) { estimateBeta <- rep_len(estimateBeta, lMin) ; print("warning - estimateBeta argument recycled") }
	
	nodeFun <- function(x) {
		betaInt <- beta[x]
		alphaInt <- alpha[x]
		offsetInt <- offset[x]
		estimateAlphaInt <- estimateAlpha[x]
		estimateBetaInt <- estimateBeta[x]
		if(estimateAlphaInt == FALSE && estimateBetaInt == FALSE) {
			alphaInt <- alphaInt ; betaInt <- betaInt
		} else {			
			if(estimateAlphaInt) {
				alphaInt <- betaInt * (minAge[x] + offsetInt)
				} 
			if(estimateBetaInt) {
				betaInt <- (minAge[x] + offsetInt) / alphaInt
				}
			}
		nodeCon <- paste0("'G[", alphaInt, "~", betaInt, "]'")
		parameters <- c(alphaInt, betaInt)
		return(list(nodeCon, parameters))
		}

	out <- sapply(1:lMin, nodeFun)
	output <- c()
	prm <- matrix(unlist(out[2,]), ncol=2, byrow=TRUE)
	rownames(prm) <- paste0("node_", 1:lMin)
	colnames(prm) <-  c("alpha", "beta")
	output$parameters <- prm
	output$apePhy <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), TRUE) 
	output$MCMCtree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), FALSE) 

	if(writeMCMCtree == TRUE) {
		outputTree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), returnPhy=FALSE) 
		utils::write.table(outputTree, paste0(MCMCtreeName), quote=FALSE, row.names=FALSE, col.names=FALSE)
		}
		
	if(plot == TRUE) {
		if(length(list.files(pattern=paste0(pdfOutput))) != 0) {
			cat(paste0("warning - deleting and over-writing file ", pdfOutput))
			file.remove(paste0(pdfOutput))
			}
	 	grDevices::pdf(paste0(pdfOutput), family="Times")
		for(k in 1:dim(prm)[1]) {
			plotMCMCtree(prm[k,], method="gamma",  paste0(rownames(prm)[k], " gamma"), upperTime = max(maxAge)+1)
			}
		grDevices::dev.off()
		}	
	output$nodeLabels <- unlist(out[1,])	
	return(output)
}
