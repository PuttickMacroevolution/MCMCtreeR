#' Read MCMCTree output tree into R
#'
#' Read MCMCTree output tree into R to produce time-scaled tree in APE format, and a table of the mean and 95% HPD ages
#' @param inputPhy file directory of 'Figtree' output from MCMCTree
#' @keywords 
#' @return apePhy time-scaled output tree from MCMCTree in APE format
#' @return nodeAges mean and 95% HPD ages for each node on the tree
#' @export
#' @examples
#' readMCMCTree()

 readMCMCTree <- function (inputPhy) 
{
         tree <- scan(paste0(inputPhy), what = "", sep = "\t")
    phys <- gsub("\\[.*?\\]", "", tree)
    phy <- read.tree(text = phys[4])
    phyInt <- phy
    phyInt$node.label <- 1:Nnode(phy)
    phyInt$edge.length <- NULL
    phyInt$root.edge <- NULL
    phyInt$tip.label <- rep(paste0(sample(letters, replace = T, 
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
    output$apePhy <- phy
    colnames(allAges) <- c("mean", "95%_lower", "95%_upper")
    output$nodeAges <- allAges
    return(output)
}