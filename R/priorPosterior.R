#' priorPosterior
#'
#' Analyse prior and posterior node distributions from MCMCtree analysis
#' @param MCMCPrior prior of MCMC file from MCMCtree analysis using data=0
#' @param MCMCPosterior posterior of MCMC file from MCMCtree analysis using data
#' @param inputTree phylogeny in MCMCtree format used in MCMCtree analysis
#' @param return.density Logical. Whether to return the density or original values for the effective prior and posterior.
#' @param rootCalibration = NULL If NULL, then behaves as by default. Alternatively, if a user has specified a root prior in the MCMCtreeR control file it can be added here as a vector in the form it would appear tree.
#' @return list containing node estimates for each distribution
#' \itemize{
#'  \item{"prior"}{ distribution of effective prior}
#'  \item{"posterior"}{ distribution of posterior}
#'  \item{"specifiedPrior"}{ distribution of specified prior}
#' }
#' @export
#' @author Mark Puttick, Pascal Title
#' @examples
#' data(MCMCtree.output)
#' # priorPosterior(MCMCPrior, 
#' # MCMCPosterior=MCMCtree.output$MCMCtree.posterior, 
#' # path.to.input.tree)

priorPosterior <- function(MCMCPrior, MCMCPosterior = NULL, inputTree, return.density=FALSE, rootCalibration = NULL) {

	if(methods::is(MCMCPrior)[1] == "character") MCMCPrior <- utils::read.csv(MCMCPrior, sep="\t")
	if(methods::is(MCMCPosterior)[1] == "character") MCMCPosterior <- utils::read.csv(MCMCPosterior, sep="\t")
	
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
	if (!is.null(rootCalibration)) {
		nodeInfo <- c(1, nodeInfo) # add root
		phylo$node.label[1] <- rootCalibration
	}

	for(tt in 1:length(nodeInfo)) {
		colNum <- as.numeric(sort(nodeInfo + Ntip(phylo))[tt])
		noders <- strsplit(phylo$node.label[nodeInfo[tt]], "~")[[1]]
		nodeType <- as.character(noders[1])
		if(nodeType == "B") {
			nums <- as.numeric(noders[2:5])
			pr <- MCMCPrior[, paste0("t_n", colNum)]
			if (!is.null(MCMCPosterior)) {
				post <- MCMCPosterior[, paste0("t_n", colNum)]
				}	
			names(nums) <- c("tL", "tU", "pL", "pU")[1:length(nums)]
			}
		if(nodeType == "SN") {
			nums <- as.numeric(noders[-1])
			pr <- MCMCPrior[, paste0("t_n", colNum)]
			if (!is.null(MCMCPosterior)) {
				post <- MCMCPosterior[,paste0("t_n", colNum)]
			}
			names(nums) <- c("location", "scale", "shape")
			}
		if(nodeType == "ST") {
			nums <- as.numeric(noders[-1])
			pr <- MCMCPrior[, paste0("t_n", colNum)]
			if (!is.null(MCMCPosterior)) {
				post <- MCMCPosterior[,paste0("t_n", colNum)]
				}
			names(nums) <- c("location", "scale", "shape", "df")
			}
		if(nodeType == "G") {
			nums <- as.numeric(noders[2:3])
			pr <- MCMCPrior[, paste0("t_n", colNum)]
			if (!is.null(MCMCPosterior)) {
				post <- MCMCPosterior[,paste0("t_n", colNum)]
			}
			names(nums) <- c("alpha", "beta")
			}
		if(nodeType == "L") {
			nums <- as.numeric(noders[-1])
			pr <- MCMCPrior[, paste0("t_n", colNum)]
			if (!is.null(MCMCPosterior)) {
				post <- MCMCPosterior[,paste0("t_n", colNum)]
			}
			names(nums) <- c("tL", "p", "c", "pL")
			}
		if(return.density) {
			pr <- stats::density(pr)
			post <- stats::density(post)
			}
		priors[[tt]] <- pr
		names(priors)[tt] <- paste0(colNum, "_", nodeType)
		posteriors[[tt]] <- post
		names(posteriors)[tt] <- paste0(colNum, "_", nodeType)
		givenPriors[[tt]] <- nums
		names(givenPriors)[tt] <- paste0(colNum, "_", nodeType)
		}
	nodeValues <- list()
	nodeValues$prior <- priors
	nodeValues$posterior <- posteriors
	nodeValues$specifiedPriors <- givenPriors
	return(nodeValues)
}
