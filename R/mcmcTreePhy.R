#' mcmcTreePhy
#'
#' Wrapper function to estimate node distributions and add them to tree and output MCMCTree format phylogeny file. If parameter values are in vectors shorter than the number of nodes they are recycled. 
#' @param minAge vector of minimum age bounds for nodes matching order in monoGroups
#' @param maxAge vector of maximum age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param offset offset value for distribution (default = 50) (p in PAML manual page 49)
#' @param scale scale value for distribution (default = 1.5) (c in PAML manual page 49)
#' @param estimateScale logical specifying whether to estimate scale with a given shape value (default = TRUE)
#' @param minProb probability of left tail (minimum bound) - default to hard minimum (minProb=0)
#' @param maxProb probability of right tail (maximum bound. default = 0.975) 
#' @param alpha alpha value for gamma distribution (default = 188) (p in PAML manual page 49)
#' @param beta beta value for gamma distribution (default = 2690) (c in PAML manual page 49)
#' @param estimateAlpha logical specifying whether to estimate alpha with a given beta value (default = TRUE)
#' @param estimateShape logical specifying whether to estimate beta with a given alpha value (default = FALSE)
#' @param rightTail probability of right tail (maximum bound default = 0.025) 
#' @param df degrees of freedom for skew-t distribution (default = 1)
#' @param shape shape value for skew-t distribution (default = 50)
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
#' estimateCauchy(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, offset=0.5, phy=apeTree, plot=F)

mcmcTreePhy <- function(phy, minAges, maxAges, monoGroups, method=c("cauchy", "upper", "bound", "gamma", "skewNormal", "skewT", "fixed"), offset=0.1, df=1, shape=50, scale=1.5, minProb=1e-8, addMode=0, maxProb=0.975, rightTail=0.025, alphaInput=188, betaInput=100, estimateScale=T, estimateShape=F, estimateMode=F, estimateAlpha=TRUE, estimateBeta=FALSE)	{
	
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
		phy <- nodeCurrent$apePhy
	  	paramsAll[[x]] <- nodeCurrent[[1]]
		paramsNode[[x]] <-  nodeCurrent$nodeLabels
	  	if(x == 	lMin) write.table(nodeCurrent$mcmctree, "mcmcTree.tre", quote=F, col.names=F, row.names=F)
 }
		
	return(paramsAll)

}