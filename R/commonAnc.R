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
		here <- match(numTips[i], tp)
		numTip <- phy$edge[here, 1]
		store <- c(store, numTip)
		while(1 == 1) {
				here <- match(numTip, tp)
				numTip <- phy$edge[here, 1]
				if(is.na(numTip)) break()
				store <- c(store, numTip)
			}
		rez[[i]] <- store	
		}
		
			targ <- length(numTips)

	k <- as.numeric(names(table(unlist(rez)))[which(table(unlist(rez)) == targ	)])
	ans <- match(k, rez[[1]])	
	commonAnc <- k[which(ans == min(ans))]
	return(commonAnc)
}