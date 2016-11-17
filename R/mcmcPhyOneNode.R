#' mcmcPhyOneNode - internal function
#'
#' Produce fixed age trees for PAML
#' @param
#' @keywords 
#' @export
#' @examples


mcmcPhyOneNode <- function(minAgeOne, maxAgeOne, phyOne, monoGroupOne, methodOne, offsetOne, dfOne, shapeOne, minProbOne, maxProbOne, scaleOne, estimateScaleOne, estimateShapeOne, alphaInputOne, betaInputOne, estimateAlphaOne, estimateBetaOne, rightTailOne, addModeOne, estimateModeOne) {
	
	if(methodOne == "cauchy") {
		nodeEsts <- estimateCauchy(minAge = minAgeOne, maxAge = maxAgeOne, phy=phyOne, scale=scaleOne, offset=offsetOne, minProb=minProbOne, maxProb=maxProbOne,  monoGroups=monoGroupOne, estimateScale=estimateScaleOne)
	}
	
	if(methodOne == "upper") {
		nodeEsts <- estimateUpper(maxAge=maxAgeOne, rightTail=rightTailOne, phy=phyOne, monoGroups=monoGroupOne)
	}
	
	if(methodOne == "fixed") {
		nodeEsts <- estimateFixed(minAge=minAgeOne, phy=phyOne, monoGroups=monoGroupOne)
	}
	
	if(methodOne == "bound") {
		nodeEsts <- estimateBound(minAge=minAgeOne, maxAge=maxAgeOne, minProb=minProbOne, rightTail=rightTailOne, phy=phyOne, monoGroups=monoGroupOne)
	}
	
		if(methodOne == "gamma") {
		nodeEsts <- estimateGamma(minAge=minAgeOne, maxAge=maxAgeOne, phy=phyOne,  alpha=alphaInputOne, beta=betaInputOne, offset=offsetOne, estimateAlpha=alphaInputOne, estimateBeta=estimateBetaOne,  monoGroups=monoGroupOne)
	}
	
			if(methodOne == "skewNormal") {
		nodeEsts <- estimateSkewNormal(minAge=minAgeOne, maxAge=maxAgeOne, phy=phyOne, shape=shapeOne, scale=scaleOne, estimateScale=estimateScaleOne, estimateShape=estimateShapeOne, monoGroups=monoGroupOne, addMode=addModeOne, maxProb=maxProbOne, minProb=minProbOne, estimateMode=estimateModeOne)
	}
	
		if(methodOne == "skewT") {
		nodeEsts <- estimateSkewt(minAge=minAgeOne, maxAge=maxAgeOne, phy=phyOne, shape=shapeOne, scale=scaleOne, df=dfOne, estimateScale=estimateScaleOne, estimateShape=estimateShapeOne,  monoGroups=monoGroupOne, addMode=addMeanOne, maxProb=maxProbOne, minProb=minProbOne, estimateMode=estimateModeOne)
	}
	
	return(nodeEsts)
	
}	

