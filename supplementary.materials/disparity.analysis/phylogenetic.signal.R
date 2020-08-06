## Puttick et al. 2020. "The complex effects of mass extinctions on
## morphological disparity", Evolution.

## Take simulated data and pruned trees of all fossils, then estimate
## phylogenetic signal of mass extinctions using the phylo D statistic and
## estimate PGLS models of trait value ~ extinction.

library(motmot)
library(caper)
library(parallel)

## Read in files.
	all.files <- "tips_50.mu_0.4.severity_0.5.random_0.01.R"
	setwd("../precooked.output/")
	load(all.files ,verbose=TRUE)
	in.file <- eval(parse(text = gsub(".R", "", all.files)))

## Extract extinction times.
	mass.ext.time <- sapply(in.file, function(x)
	  x$mass.ext.time)
	phy.d <- lm.model <- ols.model <- list()

## For 100 simulates, estimate phylo D statistic and PGLS models.
	for(u in 1:100) {
		phy <- in.file[[u]][[1]]
		y <- in.file[[u]]$trait.matrix[[1]]
		mass.ext.time <- max(nodeTimes(phy)) - in.file[[u]]$mass.ext.time
		phy.2 <- contemporaryPhy(phy, mass.ext.time + 0.0001, mass.ext.time,
		  reScale = 0, allTraits=y[,2], closest.min = TRUE, traits.from.tip = FALSE)
		nt <- nodeTimes(phy)
		dead.lin <- phy$tip.label[phy$edge[which(signif(nodeTimes(phy)[, 2], 6) ==
		  signif(mass.ext.time, 6)), 2]]
		all.dead <- match(dead.lin, phy.2[[1]]$tip.label)
		all.tips <- rep(0, Ntip(phy.2[[1]]))
		all.tips[all.dead] <- 1
		yt <- as.matrix(cbind(phy.2[[1]]$tip.label, all.tips))
		yt <- data.frame(yt)
		rownames(yt) <- phy.2[[1]]$tip.label
		colnames(yt) <- c("tips", "ext")
		comp.data <- comparative.data(phy.2[[1]], yt, tips)

		try <- FALSE
		while(try[1] == FALSE) {
		done <- TRUE
		try <- tryCatch(phylo.d(comp.data, binvar = ext), error = function(eek) {
			print("phylo D error")
			FALSE
				}
			)
		}
		phy.d[[u]] <- try
		yt <- data.frame(cbind(phy.2[[1]]$tip.label, all.tips, phy.2[[3]]))
		rownames(yt) <- phy.2[[1]]$tip.label
		colnames(yt) <- c("tips", "ext", "y.val")
		yt[, 3] <- as.numeric(as.character(yt[, 3]))
		comp.data.lm <- comparative.data(phy.2[[1]], yt, tips, vcv = TRUE)
		mat <- as.matrix(yt[, 3])
		rownames(mat) <- phy.2[[1]]$tip.label
		lm.val <- transformPhylo.ML(mat, phy.2[[1]], model = "lambda",
		  modelCIs = FALSE)$Lambda[[1]]
		mod1 <- pgls(y.val ~ ext, data = comp.data.lm, lambda = lm.val)
		mod2 <- lm(y.val ~ ext, data = yt)
		lm.model[[u]] <- summary(mod1)
		ols.model[[u]] <- summary(mod2)
	}

	signif.lm <- length(which(sapply(1:100, function(x)
	  lm.model[[x]]$coefficients[2, 4]) < 0.05))
	slope.lm <- median(sapply(1:100, function(x)
	  lm.model[[x]]$coefficients[2, 1]))
	range.lm <- range(sapply(1:100, function(x)
	  lm.model[[x]]$coefficients[2, 1]))
	lm.stats <- paste0(signif(slope.lm), " (",
	  signif(range.lm[1]), ", ", signif(range.lm[2]),")")
	all.lm <- sapply(1:100, function(x)
	  lm.model[[x]]$coefficients[2, 1])

	signif.ols <- length(which(sapply(1:100, function(x)
    ols.model[[x]]$coefficients[2, 4]) < 0.05))
	slope.ols <- median(sapply(1:100, function(x)
	  ols.model[[x]]$coefficients[2, 1]))
	range.ols <- range(sapply(1:100, function(x)
	  ols.model[[x]]$coefficients[2, 1]))
	ols.stats <- paste0(signif(slope.ols), " (", signif(range.ols[1]), ", ",
	  signif(range.ols[2]),")")
	all.ols <- sapply(1:100, function(x)
	  ols.model[[x]]$coefficients[2, 1])

	all.dstat <- sapply(phy.d, function(x)
	  x[[1]])
	d.stat <- median(unlist(all.dstat))
	d.stat.range <- range(unlist(all.dstat))
	d.stats <- paste0(signif(d.stat), " (", signif(d.stat.range[1]), ", ",
	  signif(d.stat.range[2]),")")

	signif.bm <- length(which(sapply(1:100, function(x)
	  phy.d[[x]][[3]]) < 0.05))
	signif.rnd <- length(which(sapply(1:100, function(x)
	 phy.d[[x]][[2]]) < 0.05))

	all.stats <- cbind(all.lm, all.ols, all.dstat)
	colnames(all.stats) <- c("pgls.slope", "ols.slope", "dstatistic")

	output <- list(c(lm.stats, signif.lm, ols.stats, signif.ols, signif.bm,
	  signif.rnd, d.stats), all.stats)
	signif.difference.out  <- matrix(c(lm.stats, signif.lm, ols.stats, signif.ols,
	  signif.bm, signif.rnd, d.stats))
	colnames(signif.difference.out) <- gsub(".R", "", all.files)
	rownames(signif.difference.out) <- c("lm.slope", "proportion.lm.signif",
	  "ols.slope", "proportion.ols.signif", "phylo.d.support.bm",
		"phylo.d.support.random", "phylo.d")
	write.csv(signif.difference.out,
	  "../disparity.analysis/outputs/lm.phylo.d.csv")
	# names(output) <- all.files
	# saveRDS(signif.difference,
	#  file = "../disparity.analysis/outputs/phy.d.lm.R")
