#' Plot distribution from MCMCtree node estimations
#'
#' Estimate the offset and scale paramaters of a soft-tailed cauchy distribution and output trees for MCMCtree input
#' @param parameters output parameters from node estimation function
#' @param method one of skewT, skewNormal, cauchy, gamma, or bound
#' @param title title for the plot
#' @param upperTime maxmimum age for x-axis plot
#' @param lowerTime maxmimum age for x-axis plot (default = 0)
#' @param plotMCMCtreeData If TRUE returns co-ordinates to plot distributions to allow greater flexibility (default = TRUE)
#' @return plot of the specified prior applied for a node
#' @return If plotMCMCtreeData=TRUE x and y coordinates of distributions from 0 to upperTime on x axis
#' @export
#' @author Mark Puttick
#' @examples
#' data(apeData)
#' attach(apeData)
#' # create monophyletic groups descending from nodes 8, 10, 11, and 13
#' monophyleticGroups <- tipDes(apeData$apeTree, c(8,10,11,13))
#' minimumTimes <- c("nodeOne"=15, "nodeTwo"=6,
#' "nodeThree"=8, "nodeFour"=13) / 10
#' maximumTimes <- c("nodeOne" = 30, "nodeTwo" = 12,
#' "nodeThree"=12, "nodeFour" = 20) / 10
#' cauchy <- estimateCauchy(minAge=minimumTimes[1], maxAge=maximumTimes[1],
#' monoGroups=monophyleticGroups[[1]], offset=0.5, phy=apeTree, plot=FALSE)
#' ## un-comment to run
#' plotMCMCtree(parameters=cauchy$parameters, method="cauchy",
#' title="cauchyPlot", upperTime=maximumTimes[1]+1)

plotMCMCtree <- function (parameters, method = c("skewT", "skewNormal", "cauchy", 
    "gamma", "bound"), title, upperTime, lowerTime = 0, plotMCMCtreeData = TRUE) 
{
	
	x <- NULL
    if (method == "skewT") {
        timeScale <- c(lowerTime, upperTime)
        if (plotMCMCtreeData == TRUE) {
        graphics::curve(dst(x=x, xi = parameters[1], omega = parameters[2], 
            alpha = parameters[3], nu = parameters[4]), lty = 1, 
            col = "grey20", lwd = 2, xlim = timeScale, xaxs = "i", 
            yaxs = "i", bty = "l", las = 1, ylab = "density", 
            xlab = "Ma", main = paste(title), cex = 0.8)
  			} else {
            time <- seq(0, timeScale[2], 0.01)
            density <- dst(time, xi = parameters[1], omega = parameters[2], 
                alpha = parameters[3], nu = parameters[4])
            return(cbind(time, density))
        }
    }
    if (method == "skewNormal") {
        timeScale <- c(lowerTime, upperTime)
        if (plotMCMCtreeData == TRUE) {
        graphics::curve(dsn(x, xi = parameters[1], omega = parameters[2], 
            alpha = parameters[3]), lty = 1, col = "grey20", 
            xlim = timeScale, lwd = 1, xpd = T, xaxs = "i", yaxs = "i", 
            bty = "l", las = 1, ylab = "density", xlab = "Ma", 
            main = paste(title), cex = 0.8)
  			} else {
            time <- seq(0, timeScale[2], 0.01)
            density <- dsn(time, xi = parameters[1], omega = parameters[2], 
                alpha = parameters[3])
            return(cbind(time, density))
        }
    }
    if (method == "cauchy") {
        timeScale <- c(lowerTime, upperTime)
        if (plotMCMCtreeData == TRUE) {
        plot(cauchyMCMCtree(xRange = timeScale, tL = parameters[1], 
            p = parameters[2], c = parameters[3], minProb = parameters[4]), 
            type = "l", lty = 1, col = "grey20", xlim = timeScale, 
            lwd = 1, xpd = T, xaxs = "i", yaxs = "i", bty = "l", 
            las = 1, ylab = "density", xlab = "Ma", main = paste(title), 
            cex = 0.8)
  			} else {
        out <- cauchyMCMCtree(xRange = timeScale, tL = parameters[1], 
            p = parameters[2], c = parameters[3], minProb = parameters[4])
        return(out)
        }
    }
    if (method == "gamma") {
        timeScale <- c(lowerTime, upperTime)
        if (plotMCMCtreeData == TRUE) {
        graphics::curve(stats::dgamma(x, shape = parameters[1], rate = parameters[2]), 
            lty = 1, col = "grey20", xlim = timeScale, lwd = 1, 
            xpd = T, xaxs = "i", yaxs = "i", bty = "l", las = 1, 
            ylab = "density", xlab = "Ma", main = paste(title), 
            cex = 0.8)
        	} else {
            time <- seq(0, timeScale[2], 0.01)
            density <- stats::dgamma(time, shape = parameters[1], rate = parameters[2])
            return(cbind(time, density))
        }
    }
    if (method == "bound") {
        timeScale <- c(lowerTime, upperTime)
        if (plotMCMCtreeData == TRUE) {
        plot(uniformMCMCtree(xRange = timeScale, tL = parameters[1], 
            tU = parameters[2], minProb = parameters[3], maxProb = parameters[4]), 
            lty = 1, col = "grey20", type = "l", xlim = timeScale, 
            lwd = 1, xpd = T, xaxs = "i", yaxs = "i", bty = "l", 
            las = 1, ylab = "density", xlab = "Ma", main = paste(title), 
            cex = 0.8)
            } else {
            out <- uniformMCMCtree(xRange = timeScale, tL = parameters[1], 
                tU = parameters[2], minProb = parameters[3], 
                maxProb = parameters[4])
            return(out)
        }
    }
}
