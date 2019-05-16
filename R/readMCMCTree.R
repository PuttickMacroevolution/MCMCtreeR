#' Read MCMCtree output tree into R
#'
#' Read MCMCtree output tree into R to produce time-scaled tree in APE format, and a table of the mean and 95% Highest Posterior Density
#' @param inputPhy file directory of 'Figtree' output from MCMCtree
#' @param forceUltrametric alters branch lengths at tips so tree is fully ultrametric (default = TRUE)
#' @param from.file Logical. Read a tree from file or locally from within R?
#' @return apePhy time-scaled output tree from MCMCtree in APE format
#' @return nodeAges mean and 95% Equal-tail CI ages for each node on the tree
#' @author Mark Puttick
#' @export
#' @examples
#' data(MCMCtree.output)
#' attach(MCMCtree.output)
#' # tree and node ages
#' readMCMCtree(MCMCtree.phy, from.file=FALSE)

 readMCMCtree <- function (inputPhy, forceUltrametric = TRUE, from.file=TRUE) 
{
    if(from.file) {
    	tree <- scan(paste0(inputPhy), what = "", sep = "\t", quiet=TRUE)
    	} else {
    		tree <- inputPhy
    	}
    phys <- gsub("\\[.*?\\]", "", tree)
    phy <- read.tree(text = phys[4])
    phyInt <- phy
    phyInt$node.label <- 1:Nnode(phy)
    phyInt$edge.length <- NULL
    phyInt$root.edge <- NULL
    phyInt$tip.label <- rep(paste0(sample(letters, replace = TRUE, 
        4), collapse = ""), Ntip(phyInt))
    nodingOrder <- as.numeric(as.character(unlist(strsplit(write.tree(phyInt)[[1]], 
        "[^0-9]+"))[-1]))
    openB <- gregexpr("[[]", tree[4])[[1]]
    closeB <- gregexpr("[]]", tree[4])[[1]]
    stepOne <- sapply(1:length(openB), function(r) substr(tree[4], 
        openB[r], closeB[r]))
    CIs <- t(sapply(1:length(stepOne), function(k) {
        endPoint <- gregexpr("[}]", stepOne[k])[[1]] - 1
        startPoint <- gregexpr("[{]", stepOne[k])[[1]] + 1
        as.numeric(strsplit(substr(stepOne[k], startPoint, endPoint), 
            ",")[[1]])
    }))
    reOrderNodes <- match(nodingOrder + Nnode(phy), phy$edge[, 
        1])
    reOrderNodes <- reOrderNodes + 1   
    reOrderNodes[which(is.na(reOrderNodes))] <- 1  
    CIs <- CIs[order(reOrderNodes), ]
    mean <- branching.times(phy)
    allAges <- cbind(mean, CIs)
    output <- list()
    
    if(forceUltrametric == TRUE)  {
		outer <- phy$edge[,2]
		inner <- phy$edge[,1]
		totalPath <- c()
			for(i in which(outer<=Ntip(phy))) {
				start <- i
				end <- inner[start]
				edgeTimes <- phy$edge.length[start]
				while(end != inner[1]) {
					start <- which(outer == end)
					end <- inner[start]
					edgeTimes <- c(edgeTimes, phy$edge.length[start])
					}
			  totalPath <- c(totalPath, sum(edgeTimes))
			}
			addLength <- max(totalPath) - totalPath
			phy$edge.length[which(outer <= Ntip(phy))] <- phy$edge.length[which(outer <= Ntip(phy))] + addLength
    	}
    
    output$apePhy <- phy
    colnames(allAges) <- c("mean", "95%_lower", "95%_upper")
    output$nodeAges <- allAges
    class(output) <- "MCMCtreer"
    return(output)
}
