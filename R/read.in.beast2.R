read.in.beast2 <- function(directory.beast2) {

	setwd(directory.beast2)
	tree.files <- list.files(pattern=".tree$")
	if(length(tree.files) == 0) stop(paste0("no '.tree' files in the directory ", getwd()))
	
	phy <- read.nexus(tree.files)
	scan.in <- scan(tree.files, what="", sep="\t", quiet=TRUE)
	tre.in <- scan.in [-which(scan.in == "")]
	tree.locale  <- which(regexpr("tree ", scan.in) != -1)
	tree.full <- scan.in[tree.locale]
	tre.in <- read.tree(text=tree.full)
	output <- gregexpr("height=", tree.full)[[1]]
	
	all.out <- c()
	all.len <- length(strsplit(tree.full, "")[[1]])
	
	for(i in 1:length(output)) {
		start <- substring(tree.full, output[i], all.len)
		start2 <- gregexpr("[}]", start)[[1]][1]
		all.out <- c(all.out, substring(start, 1, start2))
	}
	
	tre.in$edge.length <- NULL
	tre.in$node.label <- (Ntip(tre.in) + 1) : (Ntip(tre.in) +   Nnode(tre.in))
	tre.in$root.edge <- NULL
	tre2 <- write.tree(tre.in)
	out.1 <- strsplit(tre2, "[(]")[[1]]
	out.1  <- out.1 [- which(out.1 == "")]
	out.2 <- unlist(sapply(out.1, function(q) strsplit(q, ",")))
	out.3 <- unlist(sapply(out.2, function(q) strsplit(q, "[)]")))
	all.order <- as.numeric(gsub(";", "", out.3))
	ord.out <- order(all.order[which(all.order > Ntip(tre.in))])
	xx <- all.order[which(all.order > Ntip(tre.in))][ord.out]

	node.label <- all.out [which(all.order > Ntip(tre.in))]
	nd.1 <- sapply(node.label[ord.out], function(x) strsplit(x, "height_95%_HPD=[{]")[[1]][2])
	nd.2 <- sapply(nd.1, function(xx) strsplit(xx, ",")[[1]])
	lower.est <- as.numeric(nd.2[1,])	
	upper.est <- as.numeric(gsub("[}]", "", nd.2[2,]))
	
	HPD.out <- cbind(upper.est, lower.est)
	
	ages.out <- apply(HPD.out, 1, function(x) paste0(x[1], ",", x[2]))
	node.ages <- cbind(ages.out, (Ntip(tre.in) + 1) : (Ntip(tre.in) +  Nnode(tre.in)))
	colnames(node.ages) <- c("age_95%HPD", "node.name")
	output <- list(phy = phy, node.ages =  node.ages)
	return(output)
}