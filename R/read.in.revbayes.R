read.in.revbayes <- function(directory.rb) {
	## read in revbayes trees
	# read in tree ape
	revbayes.tree <- read.nexus(directory.rb)
	
	# scan in all data
	scan.in <- scan(directory.rb, what="", sep="\t", quiet=TRUE)
	scan.in <- scan.in[-which(scan.in == "")]

	names.tree <- (which(unlist(gregexpr("Taxlabels", scan.in)) != -1) + 1):(which(unlist(gregexpr("End;", scan.in)) != -1)[1] - 2)
	tip.names <- scan.in[names.tree]

	# find maj rule tree

	tree.locale <- which(regexpr(" = [[]&R[]]", scan.in) != -1)
	tree.full <- scan.in[tree.locale]
	tre.in <- read.tree(text=tree.full)
	output <- gregexpr("[[]", tree.full)[[1]] - 1
	output2 <- gregexpr("[]]", tree.full)[[1]] - 1

	start <- substring(tree.full, output2[1]+2, output[2])
	for(u in 2:(length(output2) -1)) start <- c(start, substring(tree.full, output2[u]+2, output[u+1]))
	
	
	
	for(x in 1:length(start)) {
		colon.start <- gregexpr(":", start[x])
		brack.end <- regexpr(")", start[x])
		comma.end <- gregexpr(",", start[x])
		if(colon.start != -1) {
			starting.point <- colon.start
			if(brack.end != -1) ending.point <- brack.end
			if(comma.end != -1) ending.point <- comma.end
			start[x] <- substring(start[x], ending.point)
		}
	}
	
	start.orig <- start
	start <- paste0(start.orig, collapse="")
	
	all.bracks <- strsplit(start, "")[[1]]
	closing <- which(all.bracks == ")")
	all.bracks[closing] <- paste0(")", (1:length(closing) + Ntip(revbayes.tree)))
	all.in <- paste0(paste0(all.bracks, collapse=""), ";")
	
	tre.in <- read.tree(text=all.in)
		
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
		index <- strsplit(minor.string, ",")[[1]]
		number.all <- number <- which(regexpr("[{]", index) != -1)
		while(length(number) > 0) {
			index[number] <- paste0(index[number], ",", index[number + 1])
			index <- index[-(number + 1)]
			number <- which(regexpr("[{]", index) != -1)
			number.all <- c(number, number.all)
			time.up <- match(number.all, number)
			if(is.numeric(time.up)) number <- number[-time.up]
		}
		
		index <- gsub("[]]", "", 	index)
		index <- gsub("[[]", "", 	index)
		index <- gsub("&", "", 	index)
		
	
		if(length(which(regexpr("age_95%_HPD", minor.string) != -1)) == 1) {
			labels.nodes <- rbind(labels.nodes, index)
		} else {
			labels.edges <- rbind(labels.edges, index)
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
	
	node.num.orig <- node.num[,1]
	node.num <- node.num[match(tre.in$node.label, node.num[,1]),]
	node.out <- node.num
	node.out <- cbind(node.out, node.num.orig)
	colnames(node.out)[4] <- "node.name"
	tips.info <- edges.num
	return(list(tip.ages=tips.info, node.ages=node.out, phy=revbayes.tree))
}
