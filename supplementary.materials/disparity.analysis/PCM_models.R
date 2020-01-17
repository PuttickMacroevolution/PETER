## XXX et al.
## Take simulated data and pruned trees of all fossils, then estimate
## parameters Phylogenetic Comparative models: Brownian motion,
## Ornstein-Uhlenbeck, and Early Burst using the motmot package.

library(motmot)

## n cores for parallel processing on Mac and Linux machines
	n.cores <- 1

## Input file.
	setwd("../precooked.output/")
	all.files <- "tips_50.mu_0.4.severity_0.5.random_0.01.R"

## Run over all iterations and output PCM model parameters
## and log-likelihood values.
	for (count in 1:length(all.files)) {
		load(all.files[count], verbose = TRUE)
		in.file <- eval(parse(text = gsub(".R", "", all.files[count])))
		bm.aicc <- eb.aicc <- ou.aicc <- c()
		for(phy.x in 1:100) {
			phy <- in.file[[phy.x]][[1]]
			extant.phy <- drop.fossil(phy)
			drop.tips <- which(is.na(match(phy$tip.label, extant.phy$tip.label)))
			phy.tip <- which(phy$edge[,2] <= Ntip(phy))
			y.all <- in.file[[phy.x]]$trait.matrix[[1]][phy.tip , 2]
			if(length(drop.tips) > 0) {
				phy.y <- as.matrix(y.all[-drop.tips])
			} else {
				phy.y <- as.matrix(y.all)
			}
		rownames(phy.y) <- extant.phy$tip.label
		bm.aicc <- rbind(bm.aicc,
		  unlist(transformPhylo.ML(y = phy.y, phy = extant.phy, model = "bm")))
		eb.aicc <- rbind(eb.aicc,
		  unlist(suppressWarnings(
			transformPhylo.ML(y = phy.y, phy = extant.phy, model = "ACDC",
			upperBound = -1e-6, modelCI = FALSE))))
		ou.aicc <- rbind(ou.aicc, suppressWarnings(
		  transformPhylo.ML(y = phy.y, phy = extant.phy, model = "OU",
			modelCI = FALSE)))
		cat("\r", paste0(phy.x, " %"))
	}
	out.all <- list(bm.aicc, eb.aicc, ou.aicc)
	name.out <- gsub(".R", ".models", all.files[count])
	assign(name.out, out.all)
	save(list = name.out,
	  file = paste0("../disparity.analysis/outputs/",
		gsub(".R", ".models.R", all.files[count])))
}
