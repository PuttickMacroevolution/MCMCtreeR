#' MCMCtreePhy
#'
#' Wrapper function to estimate node distributions and add them to tree and output MCMCtree format phylogeny file. If parameter values are in vectors shorter than the number of nodes they are recycled. 
#' @param minAges vector of minimum age bounds for nodes matching order in monoGroups
#' @param maxAges vector of maximum age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param method vector of the type of calibration distribution for each node
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
#' @param estimateShape logical specifying whether to estimate shape with a given scale value
#' @param estimateMode logical speciftying whether to estimate the scale that produces probabilities of each tail that corresponds roughly to the values given by minProb (lower tail) and maxProb (upper tail)
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
#' monophyleticGroups[[2]] <-#' tipDes(apeTree, 10)
#' monophyleticGroups[[3]] <- tipDes(apeTree, 11)
#' monophyleticGroups[[4]] <- c("sumatran", "orangutan")
#' minimumTimes <- c("nodeOne"=15, "nodeTwo"=6,
#' "nodeThree"=8, "nodeFour"=13) / 10
#' maximumTimes <- c("nodeOne" = 30, "nodeTwo" = 12,
#' "nodeThree"=12, "nodeFour" = 20) / 10
#' # Cauchy, upper age, bound, and gamma applied individually to each node
#' MCMCtreePhy(phy=apeTree, minAge=minimumTimes, maxAge=maximumTimes, 
#' monoGroups=monophyleticGroups, plot=FALSE, 
#' method=c("cauchy", "upper", "bound", "gamma"), writeMCMCtree=FALSE)


MCMCtreePhy <- function(phy, minAges, maxAges, monoGroups, method=c("cauchy", "upper", "bound", "gamma", "skewNormal", "skewT", "fixed"), offset=0.1, df=1, shape=50, scale=1.5, minProb=1e-8, addMode=0, maxProb=0.975, rightTail=0.025, alpha=188, beta=100, estimateScale=TRUE, estimateShape=FALSE, estimateMode=FALSE, estimateAlpha=TRUE, estimateBeta=FALSE, plot=FALSE, pdfOutput="nodeDistributions.pdf", writeMCMCtree=TRUE, MCMCtreeName="output.tre")	{
	
	corMethod <- match(method, c("cauchy", "upper", "bound", "gamma", "skewNormal", "skewT", "fixed"))
	if(any(is.na(corMethod))) stop("Method ", paste0("'", method[which(is.na(corMethod))], "'", collapse=" "), " not found - check case and spelling ")
	lMin <- length(minAges)
	lMax <- length(maxAges)
	if(lMin != lMax) stop("length of ages do not match")
	if(lMin != length(monoGroups)) stop("length of ages and groups do not match")
	if(lMin != length(method)) {
		print(paste("using", method[1], "distirbution for all nodes"))
		method <- rep(method, lMin)
		}
	lengthParameters <- sapply(list(offset, df, shape, minProb, maxProb, scale, estimateScale, estimateShape, alpha, beta, estimateAlpha, estimateBeta), length)
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
		if(length(alpha) < lMin) alpha <- rep(alpha[1], lMin)
		if(length(beta) < lMin) beta <- rep(beta[1], lMin)
		if(length(estimateAlpha) < lMin) estimateAlpha <- rep(estimateAlpha[1], lMin)
		if(length(estimateBeta) < lMin) estimateBeta <- rep(estimateBeta[1], lMin)
		if(length(rightTail) < lMin) rightTail <- rep(rightTail[1], lMin)
		if(length(addMode) < lMin) addMode <- rep(addMode[1], lMin)
		if(length(estimateMode) < lMin) estimateMode <- rep(estimateMode[1], lMin)
		}
	paramsAll <- list()
	paramsNode <- c()
	for(x in 1:	lMin) {		
		nodeCurrent <- MCMCPhyOneNode(minAgeOne=minAges[x], maxAgeOne=maxAges[x], phyOne=phy, monoGroupOne=monoGroups[[x]], methodOne=method[x], shapeOne=shape[x], offsetOne=offset[x],  dfOne=df[x], minProbOne=minProb[x], maxProbOne=maxProb[x], scaleOne=scale[x], estimateScaleOne=estimateScale[x], estimateShapeOne=estimateShape[x], alphaOne=alpha[x], betaOne=beta[x], estimateAlphaOne=estimateAlpha[x], estimateBetaOne=estimateBeta[x], rightTailOne=rightTail[x], addModeOne=addMode[x], estimateModeOne=estimateMode[x])
		phy <- nodeCurrent$apePhy
	  	paramsAll[[x]] <- nodeCurrent$parameters
	  	if(x == lMin && writeMCMCtree == TRUE) utils::write.table(nodeCurrent$MCMCtree, paste0(MCMCtreeName), quote=FALSE, col.names=FALSE, row.names=FALSE)
	  	}
	if(plot == TRUE) {
		grDevices::pdf(paste(pdfOutput), family="serif")		
		for(x in 1:lMin) plotMCMCtree(paramsAll[[x]], method=method[x], title=paste("node", x), upperTime=maxAges[x] + 1.5)
		grDevices::dev.off()
		}
	names(paramsAll) <- paste0("node_", 1:lMin)	
	return(paramsAll)
}
