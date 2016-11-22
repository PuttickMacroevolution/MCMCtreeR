#' Estimate Gamma Distribution for MCMCTree
#'
#' Estimate the offset and scale paramaters of a soft-tailed cauchy distribution and output trees for MCMCTree input
#' @param minAge vector of minimum age bounds for nodes matching order in monoGroups
#' @param maxAge vector of maximum age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param alpha alpha value for gamma distribution (default = 188) (p in PAML manual page 49)
#' @param beta beta value for gamma distribution (default = 2690) (c in PAML manual page 49)
#' @param offset distance of mean value from minimum bound#' @param minProb probability of left tail (minimum bound) - default to hard minimum (minProb=0)
#' @param estimateAlpha logical specifying whether to estimate alpha with a given beta value (default = TRUE)
#' @param estimateShape logical specifying whether to estimate beta with a given alpha value (default = FALSE)
#' @param plot logical specifying whether to plot to PDF
#' @param pdfOutput pdf output file name
#' @param writeMCMCTree logical whether to write tree in format that is compatible with mcmcTree to file
#' @param mcmcTreeName mcmcTree output file name
#' @keywords 
#' @return list containing node estimates for each distribution
#' \itemize{
#'  \item{"parameters"}{ estimated parameters for each node}
#'  \item{"apePhy"}{ phylogeny in ape format with node labels showing node distributions}
#'  \item{"mcmctree"}{ phylogeny in MCMCTree format}
#'  \item{"nodeLabels"}{ node labels in MCMCTreeR format}
#' }
#' @return If plot=TRUE plot of distributions in file 'pdfOutput' written to current working directory
#' @return If writeMCMCTree=TRUE tree in MCMCTree format in file "mcmcTreeName" written to current working directory
#' @export
#' @examples
#' apeTree <- read.tree(text="((((human, (chimpanzee, bonobo)), gorilla), (orangutan, sumatran)), gibbon);")
#' monophyleticGroups <- list()
#' monophyleticGroups[[1]] <- c("human", "chimpanzee", "bonobo", "gorilla", "sumatran", "orangutan", "gibbon")
#' monophyleticGroups[[2]] <- c("human", "chimpanzee", "bonobo", "gorilla")
#' monophyleticGroups[[3]] <- c("human", "chimpanzee", "bonobo")
#' monophyleticGroups[[4]] <- c("sumatran", "orangutan")
#' minimumTimes <- c("nodeOne"=15, "nodeTwo"=6, "nodeThree"=8, "nodeFour"=13) / 10
#' maximumTimes <- c("nodeOne"=30, "nodeTwo"=12, "nodeThree"=12, "nodeFour"=20) / 10
#' estimateGamma(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, alpha=188, beta=2690, offset=0.1, phy=apeTree, plot=F)


estimateGamma <- function(minAge, maxAge, phy, monoGroups, alpha=188, beta=2690, offset=0.1, estimateAlpha=TRUE, estimateBeta=F,  plot=FALSE, pdfOutput="gammaPlot.pdf", writeMCMCTree=FALSE, mcmcTreeName="gammaInput.tre") {
	
	lMin <- length(minAge)
	lMax <- length(maxAge)
	if(lMin != lMax) stop("length of ages do not match")
	if(length(alpha) < lMin) { alpha <- rep_len(alpha, lMin) ; print("warning - alpha parameter value recycled") }
	if(length(beta) < lMin) { beta <- rep_len(beta, lMin) ; print("warning - beta parameter value recycled") }
	if(length(offset) < lMin) { offset <- rep_len(offset, lMin) ; print("warning - offset parameter value recycled") }
	if(length(estimateAlpha) < lMin) { estimateAlpha <- rep_len(estimateAlpha, lMin) ; print("warning - estimateAlpha argument recycled") }
	if(length(estimateBeta) < lMin) { estimateBeta <- rep_len(estimateBeta, lMin) ; print("warning - estimateBeta argument recycled") }

		
	
	nodeFun <- function(x){
	
		betaInt <- beta[x]
		alphaInt <- alpha[x]
		offsetInt <- offset[x]
		estimateAlphaInt <- estimateAlpha[x]
		estimateBetaInt <- estimateBeta[x]

		if(estimateAlphaInt == F && estimateBetaInt == F) {
			
			alphaInt <- alphaInt ; betaInt <- betaInt
			
			} else
			{
			
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
	prm <- matrix(unlist(out[2,]), ncol=2, byrow=T)
	rownames(prm) <- paste0("node_", 1:lMin)
	colnames(prm) <-  c("alpha", "beta")
	output$parameters <- prm
	
	output$apePhy <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), T) 
	output$mcmctree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), F) 

	if(writeMCMCTree == T) {
		outputTree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), returnPhy=F) 
		write.table(outputTree, paste0(mcmcTreeName), quote=F, row.names=F, col.names=F)
		}
	if(plot == T) {
		if(length(list.files(pattern=paste0(pdfOutput))) != 0) {
			cat(paste0("warning - deleting and over-writing file ", pdfOutput))
			file.remove(paste0(pdfOutput))
			}
	 	pdf(paste0(pdfOutput), family="Times")
		for(k in 1:dim(prm)[1]) {
			plotMCMCTree(prm[k,], method="gamma",  paste0(rownames(prm)[k], " gamma"), upperTime = max(maxAge)+1)
			}
		dev.off()
		}	
	
	output$nodeLabels <- unlist(out[1,])	
	return(output)
}