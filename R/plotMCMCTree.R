#' Plot distribution from MCMCTree node estimations
#'
#' Estimate the offset and scale paramaters of a soft-tailed cauchy distribution and output trees for MCMCTree input
#' @param parameters output parameters from node estimation function
#' @param method one of skewT, skewNormal, cauchy, gamma, or bound
#' @param title title for the plot
#' @param upperTime maxmimum age for x-axis plot
#' @keywords 
#' @return plot of the specified prior applied for a node
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

plotMCMCTree <- function(parameters, method=c("skewT", "skewNormal", "cauchy", "gamma", "bound"), title, upperTime) {
		
	if(method == "skewT") {
		#timeScale <- c(0, qst(0.975, xi=parameters[1], omega=parameters[2], alpha=parameters[3], nu=parameters[4]))
		timeScale <- c(0, upperTime)
		curve(dst(x, xi=parameters[1], omega=parameters[2], alpha=parameters[3], nu=parameters[4]),lty=1, col="grey20", lwd=2, xlim=timeScale, xaxs="i", yaxs="i", bty="l", las=1, ylab="density", xlab="Ma", main=paste(title), cex=0.8)
		}
	
	if(method == "skewNormal") {	
		#timeScale <- c(0, qst(0.975, xi=parameters[1], omega=parameters[2], alpha=parameters[3]))
		timeScale <- c(0, upperTime)
		curve(dsn(x, xi=parameters[1], omega=parameters[2], alpha=parameters[3]),lty=1, col="grey20", xlim=timeScale, lwd=1, xpd=T, xaxs="i", yaxs="i", bty="l", las=1, ylab="density", xlab="Ma", main=paste(title), cex=0.8)		}
	
	if(method == "cauchy") {	
		timeScale <- c(0, upperTime)
		plot(cauchyMCMCTree(xRange=timeScale, tL=parameters[1], p=parameters[2], c=parameters[3], minProb=parameters[4]), type="l", lty=1, col="grey20", xlim=timeScale, lwd=1, xpd=T, xaxs="i", yaxs="i", bty="l", las=1, ylab="density", xlab="Ma", main=paste(title), cex=0.8)
		}	
	
	if(method == "gamma") {	
		#timeScale <- c(0, qgamma(0.995, shape=parameters[1],rate=parameters[2]))
		timeScale <- c(0, upperTime)
		curve(dgamma(x, shape=parameters[1],rate=parameters[2]), lty=1, col="grey20", xlim=timeScale, lwd=1, xpd=T, xaxs="i", yaxs="i", bty="l", las=1, ylab="density", xlab="Ma", main=paste(title), cex=0.8)
		}		
		
	if(method == "bound") {	
		timeScale <- c(0, upperTime)
		plot(uniformMCMCTree(xRange=timeScale, tL=parameters[1], tU=parameters[2], minProb=parameters[3], maxProb=parameters[4]), lty=1, col="grey20", type="l", xlim=timeScale, lwd=1, xpd=T, xaxs="i", yaxs="i", bty="l", las=1, ylab="density", xlab="Ma", main=paste(title), cex=0.8)
		}	
}