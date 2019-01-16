#' Plot distribution from MCMCTree node estimations
#'
#' Estimate the offset and scale paramaters of a soft-tailed cauchy distribution and output trees for MCMCTree input
#' @param parameters output parameters from node estimation function
#' @param method one of skewT, skewNormal, cauchy, gamma, or bound
#' @param title title for the plot
#' @param upperTime maxmimum age for x-axis plot
#' @param lowerTime maxmimum age for x-axis plot (default = 0)
#' @param plotMCMCTreeData If TRUE returns co-ordinates to plot distributions to allow greater flexibility (default = TRUE)
#' @return plot of the specified prior applied for a node
#' @return If plotMCMCTreeData=TRUE x and y coordinates of distributions from 0 to upperTime on x axis
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
#' cauchy <- estimateCauchy(minAge=minimumTimes[1], maxAge=maximumTimes[1], monoGroups=monophyleticGroups[[1]], offset=0.5, phy=apeTree, plot=F)
#' plotMCMCTree(parameters=cauchy$parameters, method="cauchy", title="cauchyPlot", upperTime=maximumTimes[1]+1)
#' bound <- estimateBound(minAge=minimumTimes[1], maxAge=maximumTimes[1], monoGroups=monophyleticGroups[[1]], minProb=0.1, rightTail=0.1, phy=apeTree, plot=F)
#' plotMCMCTree(parameters=bound$parameters, method="bound", title="uniformPlot", upperTime=maximumTimes[1]+1)

plotMCMCTree <- function (parameters, method = c("skewT", "skewNormal", "cauchy", 
    "gamma", "bound"), title, upperTime, lowerTime = 0, plotMCMCTreeData = TRUE) 
{
    if (method == "skewT") {
        timeScale <- c(lowerTime, upperTime)
        if (plotMCMCTreeData == TRUE) {
        curve(dst(x, xi = parameters[1], omega = parameters[2], 
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
        if (plotMCMCTreeData == TRUE) {
        curve(dsn(x, xi = parameters[1], omega = parameters[2], 
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
        if (plotMCMCTreeData == TRUE) {
        plot(cauchyMCMCTree(xRange = timeScale, tL = parameters[1], 
            p = parameters[2], c = parameters[3], minProb = parameters[4]), 
            type = "l", lty = 1, col = "grey20", xlim = timeScale, 
            lwd = 1, xpd = T, xaxs = "i", yaxs = "i", bty = "l", 
            las = 1, ylab = "density", xlab = "Ma", main = paste(title), 
            cex = 0.8)
  			} else {
        out <- cauchyMCMCTree(xRange = timeScale, tL = parameters[1], 
            p = parameters[2], c = parameters[3], minProb = parameters[4])
        return(out)
        }
    }
    if (method == "gamma") {
        timeScale <- c(lowerTime, upperTime)
        if (plotMCMCTreeData == TRUE) {
        curve(dgamma(x, shape = parameters[1], rate = parameters[2]), 
            lty = 1, col = "grey20", xlim = timeScale, lwd = 1, 
            xpd = T, xaxs = "i", yaxs = "i", bty = "l", las = 1, 
            ylab = "density", xlab = "Ma", main = paste(title), 
            cex = 0.8)
        	} else {
            time <- seq(0, timeScale[2], 0.01)
            density <- dgamma(time, shape = parameters[1], rate = parameters[2])
            return(cbind(time, density))
        }
    }
    if (method == "bound") {
        timeScale <- c(lowerTime, upperTime)
        if (plotMCMCTreeData == TRUE) {
        plot(uniformMCMCTree(xRange = timeScale, tL = parameters[1], 
            tU = parameters[2], minProb = parameters[3], maxProb = parameters[4]), 
            lty = 1, col = "grey20", type = "l", xlim = timeScale, 
            lwd = 1, xpd = T, xaxs = "i", yaxs = "i", bty = "l", 
            las = 1, ylab = "density", xlab = "Ma", main = paste(title), 
            cex = 0.8)
            } else {
            out <- uniformMCMCTree(xRange = timeScale, tL = parameters[1], 
                tU = parameters[2], minProb = parameters[3], 
                maxProb = parameters[4])
            return(out)
        }
    }
}