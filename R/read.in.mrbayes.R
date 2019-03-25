read.in.mrbayes <- function(directory.mb) {
	## read in mrbayes trees
	# read in tree ape
	mrbayes.tree <- read.nexus(directory.mb)
	
	# scan in all data
	scan.in <- scan(directory.mb, what="", sep="\t", quiet=TRUE)
	scan.in <- scan.in[-which(scan.in == "")]

	names.tree <- (which(unlist(gregexpr("begin trees", scan.in)) != -1) + 2):which(unlist(gregexpr("end;", scan.in)) != -1)[2]
	tree.locale <- utils::tail(names.tree, 2)[1]
	names.tree <- names.tree[-c(which(tree.locale==names.tree):length(names.tree))]

	tip.names <- scan.in[names.tree[seq(2, length(names.tree), 2)]]
	tip.names <- as.character(tip.names)
	tip.names <- gsub(",", "", tip.names)

	# find maj rule tree

	tree.full <- scan.in[tree.locale]
	tre.in <- read.tree(text=tree.full)
	output <- gregexpr("[[]", tree.full)[[1]] - 1
	output2 <- gregexpr("[]]", tree.full)[[1]] - 1

	start <- substring(tree.full, output2[1]+2, output[2])
	for(u in 2:(length(output2) -1)) start <- c(start, substring(tree.full, output2[u]+2, output[u+1]))
	start.orig <- start
	all.but <- paste0(start[seq(1, length(start)-2, 2)], ":", 1:length(tre.in$edge[,1]))
	all.in <- paste0(c(all.but, ");"), collapse="")

	tres <- read.tree(text=all.in)
	phy <- read.tree(text=paste0(paste0(start[seq(1, length(start), 2)], collapse=""), ";"))
	values <- as.numeric(gsub(":", "", start[seq(2, length(start)-1, 2)]))
	tre.in$tip.label <- tip.names[as.numeric(tre.in$tip.label)]

	start <- regexpr("[[]", tree.full)
	end <- regexpr("[]]", tree.full)
	tree.full <- substring(tree.full, end+1)
	start <- regexpr("[[]", tree.full)
	end <- regexpr("[]]", tree.full)
	tree.full.orig <- tree.full
	minor.string <- substring(tree.full, start, end)
	count <- 1
	labels.nodes <- labels.edges <- c()

	while(minor.string[1] != "") {
		start.node <- regexpr("&", minor.string) + 1
		end.node <- regexpr("[,]", minor.string) - 1
		int.names <- unlist(sapply(letters, function(x) gregexpr(paste0(",", x), minor.string)))
		int.names <- sort(int.names[-which(int.names == -1)])
		nan.number <- unlist(gregexpr(",nan", minor.string))
		inf.number <- unlist(gregexpr(",inf", minor.string))
		if(nan.number != -1 || inf.number != -1) {
			nan.match <- match(nan.number, int.names)
			if(!is.na(nan.match)) int.names <- int.names[-nan.match]
			inf.match <- match(inf.number, int.names)
			if(!is.na(inf.match)) int.names <- int.names[-inf.match]
		}
	
		start.node <- c(start.node, int.names + 1)
		end.node <- c(end.node, int.names[-1] - 1, regexpr("[]]", minor.string) - 1)
	
		if(regexpr("&prob", minor.string) == 2) {
			labels.nodes <- rbind(labels.nodes, substring(minor.string, start.node, end.node))
		} else {
			labels.edges <- rbind(labels.edges, substring(minor.string, start.node, end.node))
		}
	
		tree.full <- substring(tree.full, end+1)	
		start <- regexpr("[[]", tree.full)
		end <- regexpr("[]]", tree.full)
		minor.string <- substring(tree.full, start, end)
		count <- count + 1
	}

	node.num <- apply(labels.nodes, 2, function(xu) {
	sapply(xu, function(xx) {
		split.n <- strsplit(xx, "=")[[1]][2]
		split.n <- gsub("[{]", "", split.n)
		split.n <- gsub("[}]", "", split.n)
				}
			)
		}
	)
	colnames(node.num) <- sapply(labels.nodes[1,], function(x) strsplit(x, "=")[[1]][1])
	node.location <- c(which(regexpr("[)]", start.orig[seq(1, length(start.orig)-2, 2)]) != -1), dim(tres$edge)[1] + 1)
	
	node.local <- c(utils::tail(node.location, 1), node.location[-length(node.location)])
	
	edges.num <- apply(labels.edges, 2, function(xu) {
	sapply(xu, function(xx) {
		split.n <- strsplit(xx, "=")[[1]][2]
		split.n <- gsub("[{]", "", split.n)
		split.n <- gsub("[}]", "", split.n)
				}
			)
		}
	)
	colnames(edges.num) <- sapply(labels.edges[1,], function(x) strsplit(x, "=")[[1]][1])

	data.here <- match(c("prob", "age_median", "age_mean", "age_95%HPD"), colnames(node.num))
	
	node.out <- node.num[rev(node.location),data.here]
	
	node.which <- sapply(as.numeric(node.out[,"age_median"]), function(l) which.min(abs(l - nodeTimes(mrbayes.tree)[,1])))
	node.out <- cbind(node.out, mrbayes.tree$edge[,1][node.which])
	colnames(node.out)[5] <- "node.name"
		
	data.here.2 <- match(c("age_mean", "age_95%HPD"), colnames(node.num))
	tips.info <- node.num[tres$edge.length,data.here.2]
	tips.info.out <- tips.info[which(phy$edge[,2] <= Ntip(phy)), ]
	re.tips <- match(mrbayes.tree$tip.label, tre.in$tip.label)
	tips.info.out <- tips.info.out[re.tips,]
	return(list(tip.ages=tips.info.out, node.ages=node.out, phy=mrbayes.tree))
}
