#' Find Ancestral Node of a Clade
#'
#' This function allows you to find the common ancestor of a common clade
#' @param phy user tree in ape format
#' @param tipNames monophyletic members of a clade
#' @keywords 
#' @export
#' @examples
#' phy <- rcoal(10)
#' tips <- phy$tip.label[4:7]
#' commonAnc(phy, tipNames)


commonAnc <- function(phy, tipNames) {
	numTips <- match(tipNames, phy$tip.label)
	rez <- list()
	for(i in 1:length(numTips)) {
		store <- c()
		tp <- phy$edge[,2]
		hereLocal <- match(numTips[i], tp)
		numTip <- phy$edge[hereLocal, 1]
		stopLoop <- is.na(numTip)
		while(stopLoop == FALSE) {
				store <- c(store, numTip)
				hereLocal <- match(numTip, tp)
				numTip <- phy$edge[hereLocal, 1]
				stopLoop <- is.na(numTip)
			}
		rez[[i]] <- store	
		}
		
	targ <- length(numTips)
	nodeNames <- as.numeric(names(table(unlist(rez)))[which(table(unlist(rez)) == targ)])
	ans <- match(nodeNames, rez[[1]])	
	commonAnc <- nodeNames[which(ans == min(ans, na.rm=T))]
	return(commonAnc)
}