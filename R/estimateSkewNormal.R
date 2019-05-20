#' Estimate Skew Normal for MCMCtree analysis
#'
#' Estimate the shape, scale, and location paramaters of a skew normal distribution and output trees for MCMCtree input
#' @param minAge vector of minimum age bounds for nodes matching order in monoGroups
#' @param maxAge vector of maximum age bounds for nodes matching order in monoGroups
#' @param monoGroups list  with each element containing species that define a node of interest
#' @param phy fully resolved phylogeny in ape format
#' @param shape shape value for skew normal distribution (default = 50)
#' @param scale scale value for skew normal distribution (default = 1.5)
#' @param addMode addition to the minimum age to give the location of the distribution
#' @param maxProb probability of right tail (maximum bound default = 0.975) 
#' @param minProb probability of left tail (maximum bound default = 0.003) 
#' @param estimateScale logical specifying whether to estimate scale with a given shape value (default = TRUE)
#' @param estimateShape logical specifying whether to estimate shape with a given scale value (default = FALSE)
#' @param estimateMode logical speciftying whether to estimate the scale that produces probabilities of each tail that corresponds roughly to the values given by minProb (lower tail) and maxProb (upper tail)
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
#' @seealso \code{\link[sn]{qst}}
#' @author Mark Puttick
#' @examples
#' data(apeData)
#' attach(apeData)
#' monophyleticGroups <- list()
#' ## extract taxon descending from calibrated nodes 8, 10, 11, 13
#' ## these nodes can be visualised using plot.phylo
#' ## and nodelabels from ape
#' monophyleticGroups <- tipDes(apeData$apeTree, c(8,10,11,13))
#' minimumTimes <- c("nodeOne"=15, "nodeTwo"=6,
#' "nodeThree"=8, "nodeFour"=13) / 10
#' maximumTimes <- c("nodeOne" = 30, "nodeTwo" = 12,
#' "nodeThree"=12, "nodeFour" = 20) / 10
#' estimateSkewNormal(minAge=minimumTimes, maxAge=maximumTimes, 
#' monoGroups=monophyleticGroups, phy=apeTree, plot=FALSE)

estimateSkewNormal <- function(minAge, maxAge, monoGroups, phy, shape=50, scale=1.5, addMode=0, maxProb=0.975, minProb=0.003, estimateScale=TRUE, estimateShape=FALSE, estimateMode=FALSE, plot=FALSE, pdfOutput="skewNormalPlot.pdf", writeMCMCtree=FALSE, MCMCtreeName="skewNormalInput.tre")	{
	# parameters from the dst function in the sn package
	# dst(x, xi = 0, omega = 1, alpha = 0, nu = Inf, dp = NULL, log = FALSE) 
	# xi = location
	# alpha = shape
	# omega = scale
	# paml gives these parameters in the following order:
	# paml - location, scale, shape

	lMin <- length(minAge)
	lMax <- length(maxAge)
	if(lMin != lMax) stop("length of ages do not match")
	if(length(minProb) < lMin) { minProb <- rep_len(minProb, lMin) ; print("warning - minProb parameter value recycled") }
	if(length(maxProb) < lMin) { maxProb <- rep_len(maxProb, lMin) ; print("warning - maxProb parameter value recycled") }
	if(length(estimateScale) < lMin) { estimateScale <- rep_len(estimateScale, lMin) ; print("warning - estimateScale argument recycled") }
	if(length(estimateShape) < lMin) { estimateShape <- rep_len(estimateShape, lMin) ; print("warning - estimateShape argument recycled") }
	if(length(estimateMode) < lMin) { estimateMode <- rep_len(estimateMode, lMin); print("warning - estimateMode argument recycled") }
	if(length(shape) < lMin) { 
		shape <- rep_len(shape, lMin)
		if(any(estimateMode) || any(estimateScale)) print("warning - shape parameter value recycled")
		}
	if(length(scale) < lMin) { 
		scale <- rep_len(scale, lMin)
		if(any(estimateShape)) print("warning - scale parameter value recycled")
		}
	if(length(addMode) < lMin) { 
		addMode <- rep_len(addMode, lMin)
		if(any(estimateMode) == FALSE) print("warning - addMode parameter value recycled")
		}

	nodeFun <- function(x)	{
		upperEsts <- lowerEsts <- c()
		shapeInt <- shape[x]
		scaleInt <- scale[x]
		addModeInt <- addMode[x]
		estimateModeInt <- estimateMode[x]
		estimateShapeInt <- estimateShape[x]
		estimateScaleInt <- estimateScale[x]
		locationInt <- minAge[x] + addMode[x]
		maxProbInt <- maxProb[x]
		minProbInt <- minProb[x]
		
		if(estimateModeInt == TRUE) {
			estimateScaleInt=FALSE ; estimateShapeInt=FALSE
			}
		if(estimateShapeInt == TRUE && estimateScaleInt == TRUE) stop("need to set scale or shape")
		if(estimateShapeInt == FALSE && estimateScaleInt == FALSE && estimateModeInt == FALSE) {
			scaleInt <- scaleInt ; shapeInt <- shapeInt
		} else {
			if(estimateModeInt) {
				scaleTest <- seq (0.01, 2, by=1e-2)
				locTest <- seq(0.01, 1, by=0.01)
				for(y in 1:length(locTest)) {
					locationInt <- minAge[x] + locTest[y]
					for(u in 1:length(scaleTest)) {
						upperEsts[u] = sn::qsn(maxProbInt, xi=locationInt, omega=scaleTest[u], alpha=shapeInt)
						lowerEsts[u] = sn::qsn(minProbInt, xi=locationInt, omega=scaleTest[u], alpha=shapeInt)
						}
					closest <- which(abs(upperEsts-maxAge[x])==min(abs(upperEsts-maxAge[x])))
					closest2 <- which(abs(lowerEsts-minAge[x])==min(abs(lowerEsts-minAge[x])))
					if(sum(match(closest, ((closest2-2):(closest2+2))), na.rm=T) != 0) break()
					}
				if(sum(match(closest, ((closest2-2):(closest2+2))), na.rm=T) == 0) stop("offset and scale combination not found for Skew t distrbution. Maybe change shape and try again")
				upperEst <- round(upperEsts[closest], 2)
				lowerEst <- round(lowerEsts[closest2], 2)
				scaleInt <- round(scaleTest[closest], 2)
				locationInt <- round(locTest[y] + minAge[x], 2)
				}
				
			if(estimateScaleInt) {
				scaleTest <- seq (0.01, 2, by=1e-2)
				for(u in 1:length(scaleTest)) upperEsts[u] = sn::qsn(maxProbInt, xi=locationInt, omega=scaleTest[u], alpha=shapeInt)
				closest <- which(abs(upperEsts-maxAge[x])==min(abs(upperEsts-maxAge[x])))
				upperEst <- upperEsts[closest]
				scaleInt <- scaleTest[closest]
				}		
		
			if(estimateShapeInt) {
				shapeTest <- seq (1, 50, by=50)
				for(u in 1:length(shapeTest)) upperEsts[u] = sn::qsn(maxProbInt, xi=locationInt, omega=scaleInt, alpha=shapeTest[u])
				closest <- which(abs(upperEsts-maxAge[x])==min(abs(upperEsts-maxAge[x])))
				upperEst <- upperEsts[closest]
				shapeInt <-shapeTest[closest]
				}	
			}
			
			nodeCon <- paste0("'SN[", locationInt, "~", scaleInt, "~", shapeInt,  "]'")
			parameters <- c(locationInt, scaleInt, shapeInt)
			names(parameters) <- c("location", "scale", "shape")
			return(list(nodeCon, parameters))
		}
	out <- sapply(1:lMin, nodeFun)
	output <- c()
	prm <- matrix(unlist(out[2,]), ncol=3, byrow=T)
	rownames(prm) <- paste0("node_", 1:lMin)
	colnames(prm) <-  c("location", "scale", "shape")
	output$parameters <- prm 
	output$apePhy <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), T) 
	output$MCMCtree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), F) 
	if(writeMCMCtree == TRUE) {
		outputTree <- nodeToPhy(phy, monoGroups, nodeCon=unlist(out[1,]), returnPhy=F) 
		utils::write.table(outputTree, paste0(MCMCtreeName), quote=F, row.names=F, col.names=F)
		}
	if(plot == TRUE) {
		if(length(list.files(pattern=paste0(pdfOutput))) != 0) {
			cat(paste0("warning - deleting and over-writing file ", pdfOutput))
			file.remove(paste0(pdfOutput))
			}
	 	grDevices::pdf(paste0(pdfOutput), family="Times")
		for(k in 1:dim(prm)[1]) {
			plotMCMCtree(prm[k,], method="skewNormal",  paste0(rownames(prm)[k], " skewNormal"), upperTime = max(maxAge)+1)
			}
		grDevices::dev.off()
		}	
	output$nodeLabels <- unlist(out[1,])	
	return(output)
}
