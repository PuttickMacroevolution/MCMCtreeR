#' mcmcTreePhy
#'
#' Wrapper function to estimate node distributions and add them to tree and output MCMCTree format phylogeny file. If parameter values are in vectors shorter than the number of nodes they are recycled. 
#' @param minAges vector of minimum age bounds for nodes matching order in monoGroups
#' @param maxAges vector of maximum age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param offset offset value for distribution (default = 50)
#' @param shape shape value for skew-t distribution (default = 50)
#' @param df degrees of freedom for skew-t distribution (default = 1)
#' @param scale scale value for distribution (default = 1.5) 
#' @param minProb probability of left tail (minimum bound) - default to hard minimum (minProb=0)
#' @param maxProb probability of right tail (maximum bound. default = 0.975) 
#' @param addMode addition to the minimum age to give the location of the distribution
#' @param rightTail probability of right tail (maximum bound default = 0.025) 
#' @param alpha alpha value for gamma distribution (default = 188) 
#' @param beta beta value for gamma distribution (default = 2690) 
#' @param estimateScale logical specifying whether to estimate scale with a given shape value (default = TRUE)
#' @param estimateMode logical speciftying whether to estimate the scale that produces probabilities of each tail that corresponds roughly to the values given by minProb (lower tail) and maxProb (upper tail)
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
#' mcmcTreePhy()


mcmcTreePhy <- function(phy, minAges, maxAges, monoGroups, method=c("cauchy", "upper", "bound", "gamma", "skewNormal", "skewT", "fixed"), offset=0.1, df=1, shape=50, scale=1.5, minProb=1e-8, addMode=0, maxProb=0.975, rightTail=0.025, alphaInput=188, betaInput=100, estimateScale=T, estimateShape=F, estimateMode=F, estimateAlpha=TRUE, estimateBeta=FALSE, plot=FALSE, plotPDF="nodeDistributions.pdf", writeMCMCTree=TRUE, writeTreeName="output.tre")	{
	
	lMin <- length(minAges)
	lMax <- length(maxAges)
	if(lMin != lMax) stop("length of ages do not match")
	if(lMin != length(monoGroups)) stop("length of ages and groups do not match")
	if(lMin != length(method)) {
		print(paste("using", method[1], "distirbution for all nodes"))
		method <- rep(method, lMin)
		}
	
	lengthParameters <- sapply(list(offset, df, shape, minProb, maxProb, scale, estimateScale, estimateShape, alphaInput, betaInput, estimateAlpha, estimateBeta), length)
	if(any(lengthParameters < lMin)) {
		print("length of some parameters and nodes do not match - first parameter will be used for each node")
		if(length(offset) < lMin) offset <- rep(offset[1], lMin)
		if(length(df) < lMin) df <- rep(df[1], lMin)
		if(length(shape) < lMin) shape <- rep(shape[1], lMin)
		if(length(minProb) < lMin)  minProb <- rep(minProb[1], lMin)
		if(length(maxProb) < lMin) maxProb <- rep(maxProb[1], lMin)
		if(length(scale) < lMin) scale <- rep(scale[1], lMin)
		if(length(estimateScale) < lMin) estimateScale <- rep(estimateScale[1], lMin)
		if(length(estimateShape) < lMin) estimateShape <- rep(estimateShape[1], lMin)
		if(length(alphaInput) < lMin) alphaInput <- rep(alphaInput[1], lMin)
		if(length(betaInput) < lMin) betaInput <- rep(betaInput[1], lMin)
		if(length(estimateAlpha) < lMin) estimateAlpha <- rep(estimateAlpha[1], lMin)
		if(length(estimateBeta) < lMin) estimateBeta <- rep(estimateBeta[1], lMin)
		if(length(rightTail) < lMin) rightTail <- rep(rightTail[1], lMin)
		if(length(addMode) < lMin) addMode <- rep(addMode[1], lMin)
		if(length(estimateMode) < lMin) estimateMode <- rep(estimateMode[1], lMin)

	}
	
	paramsAll <- list()
	paramsNode <- c()
	for(x in 1:	lMin) {		
		nodeCurrent <- mcmcPhyOneNode(minAgeOne=minAges[x], maxAgeOne=maxAges[x], phyOne=phy, monoGroupOne=monoGroups[[x]], methodOne=method[x], shapeOne=shape[x], offsetOne=offset[x],  dfOne=df[x], minProbOne=minProb[x], maxProbOne=maxProb[x], scaleOne=scale[x], estimateScaleOne=estimateScale[x], estimateShapeOne=estimateShape[x], alphaInputOne=alphaInput[x], betaInputOne=betaInput[x], estimateAlphaOne=estimateAlpha[x], estimateBetaOne=estimateBeta[x], rightTailOne=rightTail[x], addModeOne=addMode[x], estimateModeOne=estimateMode[x])
	  	paramsAll[[x]] <- nodeCurrent$parameters
	  	if(x == lMin && writeMCMCTree == TRUE) write.table(nodeCurrent$mcmctree, paste0(writeTreeName), quote=F, col.names=F, row.names=F)
 }
		
if(plot == TRUE) {
	pdf(paste(plotPDF), family="serif")		
	for(x in 1:lMin) plotMCMCTree(paramsAll[[x]], method=method[x], title=paste("node", x), upperTime=maximum[x]/100)
	dev.off()
}		
		
	names(paramsAll) <- paste0("node_", 1:lMin)	
	return(paramsAll)

}