#' priorPosterior
#'
#' Analyse prior and posterior node distributions from MCMCTree analysis
#' @param mcmcPrior prior of mcmc file from MCMCTree analysis using data=0
#' @param mcmcPosterior posterior of mcmc file from MCMCTree analysis using data
#' @param inputTree phylogeny in MCMCTree format used in MCMCTree analysis
#' @keywords 
#' @return list containing node estimates for each distribution
#' \itemize{
#'  \item{"prior"}{ distribution of effective prior}
#'  \item{"posterior"}{ distribution of effective posterior}
#'  \item{"specifiedPrior"}{ distribution of specified prior}
#' }
#' @export
#' @examples
#' priorPosterior()


priorPosterior <- function(mcmcPrior, mcmcPosterior, inputTree) {

	s <- read.table(inputTree, row.names=1)[,2]
	s <- strsplit(as.character(s), "'")[[1]]
	nodePr <- c(grep("SN[(]", s), grep("ST[(]", s), grep("B[(]", s), grep("U[(]", s), grep("L[(]", s), grep("U[(]", s))
	s[nodePr] <- gsub("[(]", "~",s[nodePr])
	s[nodePr] <- gsub("[)]", "",s[nodePr])
	s[nodePr] <- gsub(",", "~",s[nodePr])
	pl <- paste0(s, collapse="")
	phylo <- ladderize(read.tree(text=pl))
	nodeInfo <- which(phylo$node.label != "")
	priors <- posteriors <- givenPriors <- list()

	for(tt in 1:length(nodeInfo)) {
	
		colNum <- as.numeric(sort(nodeInfo + Ntip(phylo))[tt])
		noders <- strsplit(phylo$node.label[nodeInfo[tt]], "~")[[1]]
		nodeType <- as.character(noders[1])
		if(nodeType == "B") {
			nums <- as.numeric(noders[2:3])
			pr <- density(mcmcPrior[, paste0("t_n", colNum)])	
			post <- density(mcmcPost[, paste0("t_n", colNum)])
			names(nums) <- c("tL", "tU", "pL", "pU")[1:length(nums)]
			}
		if(nodeType == "SN") {
			nums <- as.numeric(noders[-1])
			pr <- density(mcmcPrior[, paste0("t_n", colNum)])	
			post <- density(mcmcPost[,paste0("t_n", colNum)])
			names(nums) <- c("location", "scale", "shape")
			}
		if(nodeType == "ST") {
			nums <- as.numeric(noders[-1])
			pr <- density(mcmcPrior[, paste0("t_n", colNum)])	
			post <- density(mcmcPost[,paste0("t_n", colNum)])
			names(nums) <- c("location", "scale", "shape", "df")
			}
		if(nodeType == "G") {
			nums <- as.numeric(noders[2:3])
			pr <- density(mcmcPrior[, paste0("t_n", colNum)])	
			post <- density(mcmcPost[,paste0("t_n", colNum)])
			names(nums) <- c("alpha", "beta")
			}
		if(nodeType == "L") {
			nums <- as.numeric(noders[-1])
			pr <- density(mcmcPrior[, paste0("t_n", colNum)])	
			post <- density(mcmcPost[,paste0("t_n", colNum)])
			names(nums) <- c("tL", "p", "c", "pL")
			}
		
		priors[[tt]] <- pr ; names(priors)[tt] <- paste0(colNum, "_", nodeType)
		posteriors[[tt]] <- post ; names(posteriors)[tt] <- paste0(colNum, "_", nodeType)
		givenPriors[[tt]] <- nums ; names(givenPriors)[tt] <- paste0(colNum, "_", nodeType)
		}
	nodeValues <- list()
	nodeValues$prior <- priors ; nodeValues$posterior <- posteriors ; nodeValues$specifiedPriors <- givenPriors
	return(nodeValues)
}