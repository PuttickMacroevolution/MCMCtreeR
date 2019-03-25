#' Plot time-scaled phylogenies
#'
#' Plot time-scaled phylogenies with node uncertainty and timescale
#' @param phy A timescaled phylogeny, unless analysis.type="MCMCtree" and build.tree=TRUE
#' @param analysis.type The method used to generate the time-scale tree, one of MCMCtree, MrBayes, RevBayes, or User.
#' @param MCMC.chain The full posterior of age estimates for all nodes (default NULL)
#' @param node.ages List of user-supplied node ages applicable for analysis.type user. Either all nodes or a selection of nodes. Each list element must be named with its corresponding node label from the APE format.
#' @param directory.files The directory for files to summarise for MrBayes and RevBayes analyses
#' @param plot.type The plotting method for the phylogram corresponding to the APE definition. Phylogram is available for all analysis types, but cladogram is only avilable for MCMCtree analyses at present. Type distributions plots a phylogram with stats::density distributions on each of the nodes.
#' @param build.tree Logical. Only applicable to MCMCtree analyses, whether to timescale the phylogeny based on the full MCMC chain
#' @param node.method For plot.type phylogram the method to dispay age uncertainty on each node, either bar, node.length, or full.length. If 'none' is supplied plotting node uncertainty is suppressed. 
#' @param all.nodes If NULL (default) node uncertainty is plotted on all nodes. If node numbers are supplied, only these nodes will be labelled with uncertainty.
#' @param add.time.scale Logical. Adds a timescale to the plotted phylogeny.
#' @param add.abs.time Logical. Adds an absolute timescale alongside the geological timescale. Only applicable if add.time.scale is TRUE
#' @param scale.res The geological age designation to add to the plot can be one or a combination of Eon, Period, Epoch, Age. The order of plot (from bottom to top) is the same as the supplied order. Subsequent arguments add.abs.time and grey.bars are based on the last supplied age designation.
#' @param label.timescale.names Logical. Add names (Eon, Period, Epoch, Age, Ma) to axis if timescale used
#' @param time.correction Number to place branch lengths and age estimates in absolute time, deafult to one.
#' @param col.age The colouring of the node.method bars to summarise node age uncertainty
#' @param tip.lengths Logical. If the tree contains non-present tip ages, the age uncertainty associated with them will be added to the plot
#' @param density.col Colour of the node distributions (if analysis.type=distributions)
#' @param density.border.col Colour of the node distributions (if analysis.type=distributions) borders
#' @param cex.tips Size of the phylogeny tip labels
#' @param show.tip.label Logical. Should the tree tip labels be displayed
#' @param col.tree Colour of the phylogeny edges
#' @param tip.color Colour of the phylogeny tip labels
#' @param lwd.bar Width of the bar to summarise age uncertainty, applicable only if plot.type is phylogram and node.method is bar
#' @param grey.bars Logical. Should grey bars be used to signify time bins, applicable only if add.time.scale is TRUE
#' @param cex.age Size of the labels for the absolute timescale
#' @param cex.labels Size of the labels for the geological timescale
#' @param cex.names Logical. Add label.timescale.names, if applicable.
#' @param relative.height the relative height of the timescale labels
#' @param n.runs For analysis.type mrbayes, the number of independent chains to summarise
#' @param tip.bar.col The colour of uncertainty around non-contemporary tips
#' @param burn.in The number of points in the chain to discard for MrBayes and RevBayes analyses
#' @param distribution.height The relative height of node distributions when plot.type is distributions measured as the relative height of the descendent node
#' @param abs.age.mgp mgp values for the absolute age axis, only applicable if add.abs.time is TRUE
#' @param abs.age.lwd.ticks lwd values for the absolute age axis ticks, only applicable if add.abs.time is TRUE
#' @param abs.age.lwd lwd values for the absolute age axis horizontal line, only applicable if add.abs.time is TRUE
#' @param tck.abs.age tck values for the absolute age axis tick height, only applicable if add.abs.time is TRUE
#' @param abs.age.line line correction for the absolute age axis tick, only applicable if add.abs.time is TRUE
#' @param pos.age position position of absolute age axis
#' @param ... further arguments to be used in \code{\link[ape]{plot.phylo}}
#' @return If plot=TRUE plot of distributions in file 'pdfOutput' written to current working directory
#' @import ape
#' @import coda
#' @import sn
#' @export
#' @examples
#' data(MCMCtree.output)
#' attach(MCMCtree.output)
#' MCMCtree.file <- readMCMCtree(MCMCtree.phy, from.file=FALSE)
#' MCMC.tree.plot(MCMCtree.file[[1]],  analysis.type="MCMCtree",
#' MCMC.chain=MCMCtree.posterior, plot.type="distributions", cex.tips=0.5)

MCMC.tree.plot <- function(phy=NULL, analysis.type="MCMCtree", MCMC.chain=NULL, node.ages=NULL, directory.files=NULL, plot.type="phylogram", build.tree=FALSE, node.method="bar", all.nodes=NULL, add.time.scale=TRUE, add.abs.time=TRUE, scale.res="Epoch", label.timescale.names=FALSE, time.correction=1, col.age="blue", tip.lengths=FALSE, density.col="#00000050", density.border.col="#00000080", cex.tips=1, show.tip.label=TRUE, col.tree="black", tip.color="black", lwd.bar=1, grey.bars=TRUE, cex.age=1, cex.labels=1, cex.names=1, relative.height=0.08, tip.bar.col="#ff000050", burn.in=0.25, distribution.height=0.8, abs.age.mgp=c(3, 0.35, 0), abs.age.lwd.ticks=0.7, abs.age.lwd=0, tck.abs.age=-0.01, abs.age.line=-0.4, pos.age=NULL, n.runs=2, ...) {
	
	if(is.na(match(tolower(analysis.type), c("mcmctree", "mrbayes", "revbayes", "user")))) {
		stop("analysis.type must be one of 'MCMCtree', 'mrbayes', 'revbayes', or 'user'")
	}
	
	if(any(is.na(match(scale.res, c("Eon", "Period", "Epoch", "Age"))))) {
		stop("scale.res must be one/combination of 'Eon', 'Period', 'Epoch', or 'Age'")
	}
	
	mb.tree <- TRUE
	if(is.numeric(match(tolower(analysis.type), "mrbayes"))) mb.tree <- FALSE
	if(mb.tree && plot.type == "distributions" && is.null(MCMC.chain)) stop("The distribution plot.type options needs the full MCMC.chain")
	
	if(analysis.type == "mcmctree") {
		phy.in <- phy
		phy <- phy.in[[1]]
		}
	
	if(analysis.type == "mrbayes") {
		phy.data <- read.in.mrbayes(directory.files)
		t.name <- gsub(".con.tre", ".run1.t", directory.files)
		phy <- phy.data$phy
		if(plot.type == "cladogram") stop("sorry plot type 'cladogram' only avilable for MCMCtree trees")
		}
		
	if(analysis.type == "revbayes") {
		phy.data <- read.in.revbayes(directory.files)
		phy <- phy.data$phy
		if(plot.type == "cladogram") stop("sorry plot type 'cladogram' only avilable for MCMCtree trees")
		}
	
	if(analysis.type == "user") {
		if(plot.type == "cladogram") stop("sorry plot type 'cladogram' only avilable for MCMCtree trees")
		}	
		
	if(add.time.scale) {
		highest.level <- length(scale.res)
		check.level <- any(is.na(match(scale.res, c("Eon", "Period", "Epoch", "Age"))))
		if(is.na(check.level)) stop("scale.res must be one of Eon, Period, Epoch, Age")
		each.height <- -(relative.height * Ntip(phy)) / highest.level
		t.depth <- (each.height * highest.level)
		t.height <- Ntip(phy) + (relative.height * Ntip(phy))
		heights <- seq(from=t.depth, to=0.5, length.out=highest.level + 1)
	} else {
		t.depth <- 0
		t.height <- Ntip(phy)
	}
	
	if(analysis.type == "mcmctree") {
		if(!is.null(MCMC.chain)) {
			internal.node.local <- which(regexpr("t_n", names(MCMC.chain)) != -1)
			node.estimates <- MCMC.chain[,internal.node.local] * time.correction
			phy.matrix <- matrix(NA, ncol=2, nrow=Ntip(phy) * 2 - 2)
			node.in.ape <- as.numeric(gsub("t_n", "", names(node.estimates)))
			hpd.ages <- apply(node.estimates, 2, function(x) HPDinterval(as.mcmc(x))) 
			mean.ages <- apply(node.estimates, 2, function(x) mean(x))
			max.ages <- apply(node.estimates, 2, function(x) max(x))
			end.age <- hpd.ages[2,1]
			if(build.tree) {
				for(ii in 1:length(mean.ages)) phy.matrix[which(node.in.ape[ii] == phy$edge[,1]), 1] <- mean.ages[ii]
				ext.node <- which(phy$edge[,2] > Ntip(phy))
				phy.matrix[ext.node, 2] <- mean.ages[-1]
				phy.matrix[-ext.node, 2] <- 0
				b.lengths <- phy.matrix[,1] - phy.matrix[,2]
				phy$edge.length <- b.lengths
				phy$root.edge <- end.age - phy.matrix[1,1]
			} else {
				phy$edge.length <- phy$edge.length * time.correction
				phy$root.edge <- end.age - nodeTimes(phy)[1,1]
				}
			tip.range <- matrix(0, nrow=2, ncol=Ntip(phy))
		} else {
			phy$edge.length <- phy$edge.length * time.correction
			phy.in$nodeAges <- phy.in$nodeAges * time.correction
			hpd.ages <- t(phy.in$nodeAges[,-1])
			mean.ages <- phy.in$nodeAges[,1]
			end.age <- hpd.ages[2,1]
			node.in.ape <- as.numeric(rownames(phy.in$nodeAges))
			ext.node <- which(phy$edge[,2] > Ntip(phy))
			phy$root.edge <- end.age - mean.ages[1]
			tip.range <- matrix(0, nrow=2, ncol=Ntip(phy))
		}
	}
		
	if(analysis.type == "mrbayes") {
		node.hpd <- phy.data
		hpd.ages <- sapply(phy.data$node.ages[,"age_95%HPD"], function(xx) as.numeric(strsplit(xx, ",")[[1]]))
		end.age <- hpd.ages[2,1]
		node.in.ape <- as.numeric(phy.data$node.ages[,"node.name"])
		phy$root.edge <- end.age - nodeTimes(phy)[1,1]
		tip.range <- sapply(phy.data$tip.ages[,"age_95%HPD"], function(xx) as.numeric(strsplit(xx, ",")[[1]]))
		}
		
	if(analysis.type == "revbayes") {
		node.hpd <- phy.data
		hpd.ages <- sapply(phy.data$node.ages[,"age_95%_HPD"], function(xx) as.numeric(strsplit(xx, ",")[[1]]))
		end.age <- hpd.ages[2,1]
		node.in.ape <- as.numeric(phy.data$node.ages[,"node.name"])
		phy$root.edge <- end.age - nodeTimes(phy)[1,1]
		if(colnames(phy.data$tip.ages) == "age_95%_HPD") {
			tip.range <- sapply(phy.data$tip.ages[,"age_95%_HPD"], function(xx) as.numeric(strsplit(xx, ",")[[1]]))
		} else {
			tip.range <- matrix(0, nrow=2, ncol=Ntip(phy))
		}
	}
		
	if(analysis.type == "user") {
		node.ages <- lapply(node.ages, function(x) x * time.correction)
		phy$edge.length <- phy$edge.length * time.correction
		hpd.ages <- sapply(node.ages, function(xx) HPDinterval(as.mcmc(xx)))
		end.age <- hpd.ages[2,1]
		node.in.ape <- as.numeric(names(node.ages))
		phy$root.edge <- end.age - nodeTimes(phy)[1,1]
		if(length(node.ages) > (Nnode(phy) + 1)) {
			tip.range <- hpd.ages[-c(1:(Nnode(phy) + 1))]
		} else {
			tip.range <- matrix(0, nrow=2, ncol=Ntip(phy))
		}
		node.estimates <- node.ages
	}
		
	##########################################
	############ plot functions ##############
	##########################################
	
	############    phylogram   ##############
		
	if(plot.type == "phylogram") {
		plot.phylo(ladderize(phy), type="phylogram", edge.color=col.tree, root.edge=TRUE, y.lim=c(t.depth, t.height), cex=cex.tips, show.tip.label=show.tip.label, tip.color=tip.color, ...)
		last.plot.coord <- get("last_plot.phylo", envir = .PlotPhyloEnv)
		int <- last.plot.coord$yy[last.plot.coord$edge[,1]]
		ext <- last.plot.coord$yy[last.plot.coord$edge[,2]]
		structure <- cbind(int, ext)					
		node.in.ape <- c(node.in.ape, 1:Ntip(phy))
		mean.ages <- cbind(end.age - hpd.ages, end.age - tip.range)
		match.nodes <- last.plot.coord$edge[1,1]
		location.here <- match(match.nodes, node.in.ape)
		xx <- mean.ages[,location.here]
		match.edges <- which(last.plot.coord$edge[,1] == Ntip(phy) + 1)
		
		if(node.method != "none") {
		
		if(is.null(all.nodes)) {
			ext.node <- which(last.plot.coord$edge[,2] > Ntip(phy))
		} else {
			ext.node <- match(all.nodes, last.plot.coord$edge[,2])
		}	
			
			if(length(col.age) < length(ext.node)) col.age <- rep(col.age, length(ext.node) + 1)
			if(node.method == "full.length") {
				graphics::polygon(c(xx[1], xx[2], xx[2], xx[1]), matrix(c(0, 0, Ntip(phy), Ntip(phy)), ncol=2), col=col.age[1], border=FALSE)
				}
			if(node.method == "node.length") {
				graphics::polygon(c(xx[1], xx[2], xx[2], xx[1]), rbind(structure[match.edges,2], structure[match.edges,2]), col=col.age[1], border=FALSE)
			}
			if(node.method == "bar") {
				graphics::lines(c(xx[1], xx[2]), rep(structure[1,1], 2), col=col.age[1], lwd=lwd.bar)
			}			
			
			if(length(col.age) == length(ext.node)) col.age <- rep(col.age, length(ext.node) + 1)
			
		
		counter <- 2
		
		for(uu in ext.node) {
			match.nodes <- last.plot.coord$edge[uu,2]
			location.here <- match(match.nodes, node.in.ape)
			xx <- mean.ages[,location.here]
			match.edges <- range(which(last.plot.coord$edge[,1] == last.plot.coord$edge[uu,2]))
			if(node.method == "full.length") {
				graphics::polygon(c(xx[1], xx[2], xx[2], xx[1]), matrix(c(0, 0, Ntip(phy), Ntip(phy)), ncol=2), col=col.age[counter], border=FALSE)
				}
			if(node.method == "node.length") {
				graphics::polygon(c(xx[1], xx[2], xx[2], xx[1]), rbind(structure[match.edges,2], structure[match.edges,2]), col=col.age[counter], border=FALSE)
				}
			if(node.method == "bar") {
				graphics::lines(xx, rep(structure[uu,2], 2), col=col.age[counter], lwd=lwd.bar)
				}	
			counter = counter + 1
			}
		}
			
		if(tip.lengths) {
			phy.to <- nodeTimes(phy)[,2]
			tip.here <- which(phy$edge[,2] <= Ntip(phy))
			tip.to.label <- tip.here[which(nodeTimes(phy)[tip.here,2] > 0.01)]
			for(uu in phy$edge[tip.to.label,2] ) {
				match.nodes <- Nnode(phy) + uu
				xx <- mean.ages[,match.nodes]
				match.edges <- last.plot.coord$edge[uu,1]
				yy.loc <- last.plot.coord$yy[uu]
				graphics::lines(xx, rep(yy.loc, 2), col=tip.bar.col, lwd=lwd.bar)
			}
		}
	}
	
	############    cladogram   ##############

	if(plot.type == "cladogram") {
		plot.phylo(ladderize(phy), type="cladogram", edge.color=col.tree, root.edge=TRUE, y.lim=c(t.depth, t.height), cex=cex.tips, show.tip.label=show.tip.label, tip.color=tip.color, ...)
		last.plot.coord <- get("last_plot.phylo", envir = .PlotPhyloEnv)
		int <- last.plot.coord$yy[last.plot.coord$edge[,1]]
		ext <- last.plot.coord$yy[last.plot.coord$edge[,2]]
		structure <- cbind(int, ext)
		node.in.ape <- c(node.in.ape, 1:Ntip(phy))
		mean.ages <- cbind(end.age - hpd.ages, end.age - tip.range)
		match.nodes <- last.plot.coord$edge[1,1]
		location.here <- match(match.nodes, node.in.ape)
		xx <- mean.ages[,location.here]
		match.edges <- which(last.plot.coord$edge[,1] == Ntip(phy) + 1)
		match.nodes <- last.plot.coord$edge[1, 1]
		location.here <- which(match.nodes == last.plot.coord$edge[,1])
		xx <- mean.ages[,1]
		
		ext.node <- which(last.plot.coord$edge[,2] > Ntip(phy))
		if(length(col.age) < dim(last.plot.coord$edge)[1]) {
			col.age <- rep(col.age, dim(last.plot.coord$edge)[1])
			graphics::lines(xx, rep(structure[1,1], 2), lwd=lwd.bar, col=col.age[1])
		}		
		
		counter <- 1
		for(uu in 1:dim(last.plot.coord$edge)[1]) {
			match.nodes <- last.plot.coord$edge[uu, ]
			location.here <- match(match.nodes, node.in.ape)
			xx <- c(mean.ages[,location.here])
			graphics::polygon(c(xx[1], xx[3], xx[4], xx[2]), c(structure[uu,], rev(structure[uu,])), col=col.age[counter], border=FALSE)
			counter <- counter + 1
			}
		}
	
	############  distributions   ##############
	
	if(plot.type == "distributions") {
		
		plot.phylo(ladderize(phy), type="p", edge.color=col.tree, root.edge=TRUE, y.lim=c(t.depth, t.height), cex=cex.tips, show.tip.label=show.tip.label, tip.color=tip.color, ...)
		last.plot.coord <- get("last_plot.phylo", envir = .PlotPhyloEnv)
		last.plot.coord <- get("last_plot.phylo", envir = .PlotPhyloEnv)
		
		if(analysis.type == "mrbayes") {
			internal.nodes <- unique(phy$edge[,1])
			node.times <- vector("list", length(internal.nodes))
			for(kk in 1:n.runs) {
				tre.in.ape <- read.nexus(gsub("run1", paste0("run", kk), t.name))
				tre.in <- scan(gsub("run1", paste0("run", kk), t.name), what="", sep="\n", quiet=TRUE)
				start.names <- which(regexpr("translate", tre.in) != -1) + 1
				end.names <- which(regexpr("tree gen", tre.in) != -1)[1] - 1
				names.on.tree <- sapply(tre.in[start.names:end.names], function(xx) {
					xx <- gsub(",", "", utils::tail(strsplit(xx, " ")[[1]], 1))
					xx <- gsub(";", "", xx)
					}
				)
				clock.rate.in <- unlist(gregexpr("&clockrate", tre.in))
				clock.rate.loc <- which(clock.rate.in != -1)
				ceiling.rm <- ceiling(length(tre.in.ape) * burn.in)
				tre.in.ape <- tre.in.ape[-c(1:ceiling.rm)]
				clock.rate.loc <- clock.rate.loc[-c(1:ceiling.rm)]
				for(y in 1:length(tre.in.ape)) {
					tree.here <- tre.in[clock.rate.loc[y]]
					new.string <- substring(tree.here, regexpr("&clockrate", tree.here) + 1, )	
					new.string <- substring(new.string, 1, regexpr("[]]", new.string) - 1)	
					clock.rate <- as.numeric(strsplit(new.string, "=")[[1]][2])
					tree.in <- read.tree(text=tree.here)
					tree.in$tip.label <- as.character(names.on.tree[as.numeric(tree.in$tip.label)])
					tree.edge.unique <- unique(tree.in$edge[,1])
					node.times.now <- nodeTimes(tree.in) / clock.rate
					parts <- prop.part(tree.in, phy)
					match.on.tree <- which(attributes(parts)$number == 2)
					
					for(ii in 1:length(match.on.tree)) {
						name.now <- parts[[match.on.tree[ii]]]
						five <- match(tree.edge.unique[match.on.tree[ii]], tree.in$edge[,1])
						node.time.add <- node.times.now[five, 1]
						names.now <- attributes(parts)$labels[name.now]
						mrca.now <- getMRCA(phy, names.now)
						mrca.store <- match(mrca.now, internal.nodes)
						node.times[[mrca.store]] <- c(node.times[[mrca.store]], node.time.add)
					}
				}
			}
			
			names(node.times) <- internal.nodes
			all.nodes <- (Ntip(phy) + 1 ):(Nnode(phy) + Ntip(phy))
			int <- last.plot.coord$yy[last.plot.coord$edge[,1]]
			ext <- last.plot.coord$yy[last.plot.coord$edge[,2]]
			len.all.nd <- length(all.nodes)
			if(length(density.col) < len.all.nd) density.col <- rep(density.col, len.all.nd)
			if(length(density.border.col) < len.all.nd) density.border.col <- rep(density.border.col, len.all.nd)
		
			for(k in 1:len.all.nd) {
				location.here <- match(as.numeric(names(node.times)[k]), node.in.ape)
				lower.hpd <- hpd.ages[1, location.here]
				upper.hpd <- hpd.ages[2, location.here]
				MCMC.hpd <- node.times[[k]][intersect(which(node.times[[k]] >= lower.hpd), which(node.times[[k]] <= upper.hpd))]
				MCMC.den <- stats::density(MCMC.hpd)
				xx <- last.plot.coord$xx[1] - MCMC.den$x
				lower.y <- last.plot.coord$yy[all.nodes[k]]
				upper.y <- max(ext[which(all.nodes[k] == last.plot.coord$edge[,1])])	
				scale.den <- (upper.y - lower.y) * distribution.height
				y.scale <- MCMC.den$y * (scale.den / max(MCMC.den$y))
				yy <- y.scale + lower.y
				graphics::polygon(xx, yy, col=density.col[k], border=density.border.col[k])
				}
			}
			
			
		if(analysis.type == "revbayes") {
			internal.nodes <- unique(phy$edge[,1])
			node.times <- vector("list", length(internal.nodes))
			slash <- max(gregexpr("/", directory.files)[[1]])
			path.to.files <- substring(directory.files, 1, slash)
			log.files <- list.files(path=path.to.files, pattern=".trees")
			n.runs <- length(log.files)
			
			for(kk in 1:n.runs) {
				tre.in.ape <- read.tree(paste0(path.to.files, "/", log.files[kk]))
				ceiling.rm <- ceiling(length(tre.in.ape) * burn.in)
				tre.in.ape <- tre.in.ape[-c(1:ceiling.rm)]
			
				for(y in 1:length(tre.in.ape)) {
					tree.in <- tre.in.ape[[y]]
					tree.edge.unique <- unique(tree.in$edge[,1])
					node.times.now <- nodeTimes(tree.in)
					parts <- prop.part(tree.in, phy)
					match.on.tree <- which(attributes(parts)$number == 2)
					
					for(ii in 1:length(match.on.tree)) {
						name.now <- parts[[match.on.tree[ii]]]
						five <- match(tree.edge.unique[match.on.tree[ii]], tree.in$edge[,1])
						node.time.add <- node.times.now[five, 1]
						names.now <- attributes(parts)$labels[name.now]
						mrca.now <- getMRCA(phy, names.now)
						mrca.store <- match(mrca.now, internal.nodes)
						node.times[[mrca.store]] <- c(node.times[[mrca.store]], node.time.add)
					}
				}
			}
			
			names(node.times) <- internal.nodes
			all.nodes <- (Ntip(phy) + 1 ):(Nnode(phy) + Ntip(phy))
			int <- last.plot.coord$yy[last.plot.coord$edge[,1]]
			ext <- last.plot.coord$yy[last.plot.coord$edge[,2]]	
			len.all.nd <- length(all.nodes)
			if(length(density.col) < len.all.nd) density.col <- rep(density.col, len.all.nd)
			if(length(density.border.col) < len.all.nd) density.border.col <- rep(density.border.col, len.all.nd)
		
			for(k in 1:len.all.nd) {
				location.here <- match(as.numeric(names(node.times)[k]), node.in.ape)
				lower.hpd <- hpd.ages[1, location.here]
				upper.hpd <- hpd.ages[2, location.here]
				MCMC.hpd <- node.times[[k]][intersect(which(node.times[[k]] >= lower.hpd), which(node.times[[k]] <= upper.hpd))]
				MCMC.den <- stats::density(MCMC.hpd)
				xx <- last.plot.coord$xx[1] - MCMC.den$x
				lower.y <- last.plot.coord$yy[all.nodes[k]]
				upper.y <- max(ext[which(all.nodes[k] == last.plot.coord$edge[,1])])	
				scale.den <- (upper.y - lower.y) * distribution.height
				y.scale <- MCMC.den$y * (scale.den / max(MCMC.den$y))
				yy <- y.scale + lower.y
				graphics::polygon(xx, yy, col=density.col[k], border=density.border.col[k])
				}
			}	

		if(analysis.type == "mcmctree") {
			all.nodes <- (Ntip(phy) + 1 ):(Nnode(phy) + Ntip(phy))
			int <- last.plot.coord$yy[last.plot.coord$edge[,1]]
			ext <- last.plot.coord$yy[last.plot.coord$edge[,2]]	
			len.all.nd <- length(all.nodes)
			if(length(density.col) < len.all.nd) density.col <- rep(density.col, len.all.nd)
			if(length(density.border.col) < len.all.nd) density.border.col <- rep(density.border.col, len.all.nd)
			
			for(k in 1:len.all.nd) {
				lower.hpd <- hpd.ages[1, k]
				upper.hpd <- hpd.ages[2, k]
				MCMC.hpd <- node.estimates[,k][intersect(which(node.estimates[,k] >= lower.hpd), which(node.estimates[,k] <= upper.hpd))]
				MCMC.den <- stats::density(MCMC.hpd)
				xx <- last.plot.coord$xx[1] - MCMC.den$x
				lower.y <- last.plot.coord$yy[all.nodes[k]]
				upper.y <- max(ext[which(all.nodes[k] == last.plot.coord$edge[,1])])	
				scale.den <- (upper.y - lower.y) * distribution.height
				y.scale <- MCMC.den$y * (scale.den / max(MCMC.den$y))
				yy <- y.scale + lower.y
				graphics::polygon(xx, yy, col=density.col[k], border=density.border.col[k])
				}
			}
			
		if(analysis.type == "user") {
			all.nodes <- as.numeric(names(node.ages))
			int <- last.plot.coord$yy[last.plot.coord$edge[,1]]
			ext <- last.plot.coord$yy[last.plot.coord$edge[,2]]
			len.all.nd <- length(all.nodes)
			if(length(density.col) < len.all.nd) density.col <- rep(density.col, len.all.nd)
			if(length(density.border.col) < len.all.nd) density.border.col <- rep(density.border.col, len.all.nd)
						
			for(k in 1:len.all.nd) {
				lower.hpd <- hpd.ages[1, k]
				upper.hpd <- hpd.ages[2, k]
				node.now <- node.estimates[[k]]
				MCMC.hpd <- node.now[intersect(which(node.now >= lower.hpd), which(node.now <= upper.hpd))]
				MCMC.den <- stats::density(MCMC.hpd)
				xx <- last.plot.coord$xx[1] - MCMC.den$x
				lower.y <- last.plot.coord$yy[all.nodes[k]]
				upper.y <- max(ext[which(last.plot.coord$edge[,1] == all.nodes[k])])
				scale.den <- (upper.y - lower.y) * distribution.height
				y.scale <- MCMC.den$y * (scale.den / max(MCMC.den$y))
				yy <- y.scale + lower.y
				graphics::polygon(xx, yy, col=density.col[k], border=density.border.col[k])
				}
			}
		}
		
	##########################################
	############ add time scale ##############
	##########################################

		if(add.time.scale) {
			xlimit <- last.plot.coord$xx
			total.length <- max.xlim <- max(xlimit)
			eon.bin <- .bincode(max.xlim, Eon)
			period.bin <- .bincode(max.xlim, Period)
			epoch.bin <- .bincode(max.xlim, Epoch)
			if(is.na(epoch.bin)) epoch.bin <- length(Epoch)
			age.bin <- .bincode(max.xlim, Age)
			if(is.na(age.bin)) age.bin <- length(Age)
			
			bin.eon <- max.xlim - c(Eon[1:eon.bin], max.xlim)
			bin.period <- max.xlim - c(Period[1:period.bin], max.xlim)
			bin.epoch <- max.xlim - c(Epoch[1:epoch.bin], max.xlim)
			bin.age <- max.xlim - c(Age[1:age.bin], max.xlim)
			
			plot.bar <- function(bin.data, data.colour, name.data, ht, ht.pl, grey.line=FALSE, cex.labels.int) {
				names.in <- name.data[1:(length(bin.data) - 1)]
				for(u in 2:length(bin.data)) {
					graphics::polygon(rep(bin.data[(u):(u-1)], each=2),c(heights[ht:ht.pl], rev(heights[ht:ht.pl])), xpd=TRUE, col=data.colour[c(u-1)], border=FALSE)
					bin.length <- abs(diff(bin.data[(u):(u-1)]))
					percent <- bin.length / total.length
					if(percent >= 0.4) {
						graphics::text(mean(bin.data[(u):(u-1)]), mean(heights[ht:ht.pl]), names.in[c(u-1)], cex=cex.labels.int)
					} 
					if(percent >= 0.1 && percent < 0.4) {
						name.in <- paste0(paste0(strsplit(names.in[c(u-1)], "")[[1]][1:2], collapse=""), ".")
						graphics::text(mean(bin.data[(u):(u-1)]), mean(heights[ht:ht.pl]), name.in, cex=cex.labels.int)
					}
					if(percent < 0.1 && percent >= 0.04) {
						name.in <- paste0(paste0(strsplit(names.in[c(u-1)], "")[[1]][1], collapse=""), ".")
						graphics::text(mean(bin.data[(u):(u-1)]), mean(heights[ht:ht.pl]), name.in, cex=cex.labels.int)
					}	
					if(grey.line) {
						if(u %% 2 == 0) graphics::polygon(rep(bin.data[(u):(u-1)], each=2),c(0.5, Ntip(phy), Ntip(phy), 0.5), xpd=TRUE, col="#00000020", border=FALSE)
					}	
				}
			}
			
			plot.time.scale <- function(level, start.one, end.one, time.check=FALSE, grey.now=FALSE, cex.age.int) {
				
				
				
				
				if(level == "Eon") {
					plot.bar(bin.eon, Eon.colour, names(Eon.colour),start.one,end.one, grey.line=grey.now, cex.labels.int=cex.labels)
					# if(label.timescale.names) graphics::text(total.length, mean(heights[start.one:end.one]), "Eon", pos=4, font=4, cex=cex.labels)
					if(label.timescale.names) graphics::mtext("Eon", 4, at=mean(heights[start.one:end.one]), line=0, font=4, cex=cex.names, las=1)
					if(time.check) {
						if(is.null(pos.age)) pos.age <- Ntip(phy) + (Ntip(phy) / 100)
						graphics::axis(3, at=bin.eon, labels=signif(total.length-bin.eon, 2), line=abs.age.line, lwd=abs.age.lwd, lwd.ticks=abs.age.lwd.ticks, mgp=abs.age.mgp, cex.axis=cex.age.int, pos=pos.age, tck=tck.abs.age)
						# if(label.timescale.names) graphics::text(total.length, Ntip(phy) + (Ntip(phy) / 35), "Ma", pos=4, font=4, cex=cex.labels, adj=c(0,1), xpd=TRUE)
						if(label.timescale.names) graphics::mtext("Ma",4, line=0, at=Ntip(phy) + (Ntip(phy) / 35), xpd=TRUE, cex=cex.names, font=4, las=1)
					}
				}
			
				if(level == "Period") {
					plot.bar(bin.period, Period.colour, names(Period.colour),start.one,end.one, grey.line=grey.now, cex.labels.int=cex.labels)
					# if(label.timescale.names) graphics::text(total.length, mean(heights[start.one:end.one]), "Period", pos=4, font=4, cex=cex.labels)
					if(label.timescale.names) graphics::mtext("Period", 4, at=mean(heights[start.one:end.one]), line=0, font=4, cex=cex.names, las=1)
					if(time.check) {
						if(is.null(pos.age)) pos.age <- Ntip(phy) + (Ntip(phy) / 100)
						graphics::axis(3, at=bin.period, labels=signif(total.length-bin.period, 2), line=abs.age.line, lwd=abs.age.lwd, lwd.ticks=abs.age.lwd.ticks, mgp=abs.age.mgp, cex.axis=cex.age.int, pos=pos.age, tck=tck.abs.age)
						# if(label.timescale.names) graphics::text(total.length, Ntip(phy) + (Ntip(phy) / 35), "Ma", pos=4, font=4, cex=cex.labels, adj=c(0,1), xpd=TRUE)
						if(label.timescale.names) graphics::mtext("Ma",4, line=0, at=Ntip(phy) + (Ntip(phy) / 35), xpd=TRUE, cex=cex.names, font=4, las=1)
					}
				}
			
				if(level == "Epoch") {
					plot.bar(bin.epoch, Epoch.colour, names(Epoch.colour),start.one,end.one, grey.line=grey.now, cex.labels.int=cex.labels)
					# if(label.timescale.names) graphics::text(total.length, mean(heights[start.one:end.one]), "Epoch", pos=4, font=4, cex=cex.labels)
					if(label.timescale.names) graphics::mtext("Epoch", 4, at=mean(heights[start.one:end.one]), line=0, font=4, cex=cex.names, las=1)
					if(time.check) {
						if(is.null(pos.age)) pos.age <- Ntip(phy) + (Ntip(phy) / 100)
						graphics::axis(3, at=bin.epoch, labels=signif(total.length-bin.epoch, 2), line=abs.age.line, lwd=abs.age.lwd, lwd.ticks=abs.age.lwd.ticks, mgp=abs.age.mgp, cex.axis=cex.age.int, pos=pos.age, tck=tck.abs.age, ...)
						# if(label.timescale.names) graphics::text(total.length, Ntip(phy) + (Ntip(phy) / 35), "Ma", pos=4, font=4, cex=cex.labels, adj=c(0,1), xpd=TRUE)
						if(label.timescale.names) graphics::mtext("Ma",4, line=0, at=Ntip(phy) + (Ntip(phy) / 35), xpd=TRUE, cex=cex.names, font=4, las=1)
					}
				}
			
				if(level == "Age") {
					plot.bar(bin.age, Age.colour, names(Age.colour),start.one,end.one, grey.line=grey.now, cex.labels.int=cex.labels)
					# if(label.timescale.names) graphics::text(total.length, mean(heights[start.one:end.one]), "Age", pos=4, font=4, cex=cex.labels)
					if(label.timescale.names) graphics::mtext("Age", 4, at=mean(heights[start.one:end.one]), line=0, font=4, cex=cex.names, las=1)
					if(time.check) {
						if(is.null(pos.age)) pos.age <- Ntip(phy) + (Ntip(phy) / 100)
						graphics::axis(3, at=bin.age, labels=signif(total.length-bin.age, 2), line=abs.age.line, lwd=abs.age.lwd, lwd.ticks=abs.age.lwd.ticks, mgp=abs.age.mgp, cex.axis=cex.age.int, pos=pos.age, tck=tck.abs.age)
						# if(label.timescale.names) graphics::text(total.length, Ntip(phy) + (Ntip(phy) / 35), "Ma", pos=4, font=4, cex=cex.labels, adj=c(0,1), xpd=TRUE)
						if(label.timescale.names) graphics::mtext("Ma",4, line=0, at=Ntip(phy) + (Ntip(phy) / 35), xpd=TRUE, cex=cex.names, font=4, las=1)
					}
				}
			
			}
						
			start.one <- 1
			end.one <- 2
			if(length(add.abs.time) == 1 && add.abs.time[1] == TRUE) {
				add.abs.time <- rep(FALSE, length(scale.res))
				add.abs.time[length(add.abs.time)] <- TRUE
				}
			
			if(length(add.abs.time) == 1 && add.abs.time[1] == FALSE) {
				add.abs.time <- rep(FALSE, length(scale.res))
				}
				
			if(grey.bars) {
				grey.now.int <- rep(FALSE, length(scale.res))
				if(!any(add.abs.time)) {
					add.time <- length(scale.res)
					} else {
					add.time <- which(add.abs.time)	
					}
				grey.now.int[add.time] <- TRUE
			}
			
			for(uu in 1:length(scale.res)) {
				plot.time.scale(scale.res[uu], time.check=add.abs.time[uu], start.one=start.one, end.one=end.one, grey.now=grey.now.int[uu], cex.age.int=cex.age)
				start.one <- start.one + 1
				end.one <- end.one + 1
				}
		}
	}
