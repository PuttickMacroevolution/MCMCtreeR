# MCMCPhyOneNode - internal function
MCMCPhyOneNode <- function(minAgeOne, maxAgeOne, phyOne, monoGroupOne, methodOne, offsetOne, dfOne, shapeOne, minProbOne, maxProbOne, scaleOne, estimateScaleOne, estimateShapeOne, alphaOne, betaOne, estimateAlphaOne, estimateBetaOne, rightTailOne, addModeOne, estimateModeOne) {
	
	if(methodOne == "cauchy") {
		nodeEsts <- estimateCauchy(minAge = minAgeOne, maxAge = maxAgeOne, phy=phyOne, scale=scaleOne, offset=offsetOne, minProb=minProbOne, maxProb=maxProbOne,  monoGroups=monoGroupOne, estimateScale=estimateScaleOne, plot=F)
	}
	
	if(methodOne == "upper") {
		nodeEsts <- estimateUpper(maxAge=maxAgeOne, rightTail=rightTailOne, phy=phyOne, monoGroups=monoGroupOne)
	}
	
	if(methodOne == "fixed") {
		nodeEsts <- estimateFixed(minAge=minAgeOne, phy=phyOne, monoGroups=monoGroupOne)
	}
	
	if(methodOne == "bound") {
		nodeEsts <- estimateBound(minAge=minAgeOne, maxAge=maxAgeOne, minProb=minProbOne, rightTail=rightTailOne, phy=phyOne, monoGroups=monoGroupOne, plot=F)
	}
	
		if(methodOne == "gamma") {
		nodeEsts <- estimateGamma(minAge=minAgeOne, maxAge=maxAgeOne, phy=phyOne,  alpha=alphaOne, beta=betaOne, offset=offsetOne, estimateAlpha=alphaOne, estimateBeta=estimateBetaOne,  monoGroups=monoGroupOne, plot=F)
	}
	
			if(methodOne == "skewNormal") {
		nodeEsts <- estimateSkewNormal(minAge=minAgeOne, maxAge=maxAgeOne, phy=phyOne, shape=shapeOne, scale=scaleOne, estimateScale=estimateScaleOne, estimateShape=estimateShapeOne, monoGroups=monoGroupOne, addMode=addModeOne, maxProb=maxProbOne, minProb=minProbOne, estimateMode=estimateModeOne, plot=F)
	}
	
		if(methodOne == "skewT") {
		nodeEsts <- estimateSkewT(minAge=minAgeOne, maxAge=maxAgeOne, phy=phyOne, shape=shapeOne, scale=scaleOne, df=dfOne, estimateScale=estimateScaleOne, estimateShape=estimateShapeOne,  monoGroups=monoGroupOne, addMode=addModeOne, maxProb=maxProbOne, minProb=minProbOne, estimateMode=estimateModeOne, plot=F)
	}	
	return(nodeEsts)
}
