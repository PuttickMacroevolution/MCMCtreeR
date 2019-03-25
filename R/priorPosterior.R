#' priorPosterior
#'
#' Analyse prior and posterior node distributions from MCMCtree analysis
#' @param MCMCPrior prior of MCMC file from MCMCtree analysis using data=0
#' @param MCMCPosterior posterior of MCMC file from MCMCtree analysis using data
#' @param inputTree phylogeny in MCMCtree format used in MCMCtree analysis
#' @return list containing node estimates for each distribution
#' \itemize{
#'  \item{"prior"}{ distribution of effective prior}
#'  \item{"posterior"}{ distribution of effective posterior}
#'  \item{"specifiedPrior"}{ distribution of specified prior}
#' }
#' @export
#' @examples
#' data(MCMCtree.output)
#' # priorPosterior(MCMCPrior, 
#' # MCMCPosterior=MCMCtree.output$MCMCtree.posterior, 
#' # path.to.input.tree)

priorPosterior <- function(MCMCPrior, MCMCPosterior, inputTree) {


	s <- utils::read.table(inputTree, row.names=1)[,2]
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
			nums <- as.numeric(noders[2:5])
			pr <- stats::density(MCMCPrior[, paste0("t_n", colNum)])	
			post <- stats::density(MCMCPosterior[, paste0("t_n", colNum)])
			names(nums) <- c("tL", "tU", "pL", "pU")[1:length(nums)]
			}
		if(nodeType == "SN") {
			nums <- as.numeric(noders[-1])
			pr <- stats::density(MCMCPrior[, paste0("t_n", colNum)])	
			post <- stats::density(MCMCPosterior[,paste0("t_n", colNum)])
			names(nums) <- c("location", "scale", "shape")
			}
		if(nodeType == "ST") {
			nums <- as.numeric(noders[-1])
			pr <- stats::density(MCMCPrior[, paste0("t_n", colNum)])	
			post <- stats::density(MCMCPosterior[,paste0("t_n", colNum)])
			names(nums) <- c("location", "scale", "shape", "df")
			}
		if(nodeType == "G") {
			nums <- as.numeric(noders[2:3])
			pr <- stats::density(MCMCPrior[, paste0("t_n", colNum)])	
			post <- stats::density(MCMCPosterior[,paste0("t_n", colNum)])
			names(nums) <- c("alpha", "beta")
			}
		if(nodeType == "L") {
			nums <- as.numeric(noders[-1])
			pr <- stats::density(MCMCPrior[, paste0("t_n", colNum)])	
			post <- stats::density(MCMCPosterior[,paste0("t_n", colNum)])
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
