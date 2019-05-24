## Puttick et al. 
## Code to simulate diversity and trait data through time with a background extinction and/or
## speciation biased according to trait values in a directional selection procedure.

## At each speciation or extinction event, the algorithm ranks all contemporary trait values and then 
## raised them to a power value (of 0.1, 0.5, 2, or 5). After being raised to this power the probabilities of 
## speciation and/or extinction are not even among lineages, so at each time step the probability of 
## speciation and/or extinction is biased towards lineages with specific trait values. 

## After running these analyses, we trimmed phylogenies and data be ultrametric, ‘extant’ lineages 
## present at the end of the simulation, and estimated the relative fit of Brownian motion, OU, 
## and Early Burst models. 

library(RSOLE)
library(parallel)

## n cores for parallel processing on Mac and Linux machines
	n.cores <- 1

## Function to generate background extinction models with directional selection and run PCM models
	background.ext <- function(mu=0, extant.tips=50, ext.size.prob=1, spp.size.prob=0, both.sides=TRUE, n.core=1) {
		outputs <- mclapply(1:100, mc.cores=n.core, function(x) {
			est.out <- build.tree.tips(lambda=1, mu=mu, extant.tips=extant.tips, mass.ext.occ=FALSE, ext.size.prob=ext.size.prob, spp.size.prob=spp.size.prob, both.sides=both.sides)
			tip.values <- est.out$trait.matrix[[1]][which(est.out$phy$edge[,2] <= Ntip(est.out$phy)),2]
			if(mu > 0) {
				names(tip.values) <- est.out$phy$tip.label
				extant.tree <- drop.fossil(est.out$phy)
				ext.tip.val <- 
	tip.values[match(extant.tree$tip.label, names(tip.values))] 
				ext.tip.val <- as.matrix(ext.tip.val)
			} else {
				ext.tip.val <- as.matrix(tip.values)
				rownames(ext.tip.val) <- est.out$phy$tip.label
				extant.tree <- est.out$phy
			}
			bm <- transformPhylo.ML(y=ext.tip.val, phy=extant.tree, model="bm")
			ou <- transformPhylo.ML(y=ext.tip.val, phy=extant.tree, model="OU")
			eb <- transformPhylo.ML(y=ext.tip.val, phy=extant.tree, model="ACDC", upperBound=-1e-6)
			c(bm, ou, eb, est.out[1:3])
			}
		)
	}

## lineages with smaller trait values were more likely to speciate compared to lineages with larger trait values (equal extinction probabilities); 

	for(mu.int in c(0, 0.4, 0.8)) {
		for(tips.int in c(50, 100, 200)) {
			for(ext.int in c(0, 0.1, 0.5, 2, 5)) {
				int.save <- background.ext(mu=mu.int, extant.tips=tips.int, ext.size.prob=ext.int, spp.size.prob=0, both.sides=FALSE, n.core=n.cores)
				name.out <- paste0("tips_", tips.int, ".mu_", mu.int, ".ext.size_", ext.int, ".spp.size_0")
				assign(name.out, int.save)
				}
			}
		}

## lineages with larger trait values were more likely to become extinct compared to lineages with 
## smaller trait values (equal speciation probabilities); 
	for(mu.int in c(0, 0.4, 0.8)) {
		for(tips.int in c(50, 100, 200)) {
			for(spp.int in c(0.1, 0.5, 2, 5)) {
				int.save <- background.ext(mu=mu.int, extant.tips=tips.int, spp.size.prob=-spp.int, ext.size.prob=0, both.sides=FALSE, n.core=n.cores)
				name.out <- paste0("tips_", tips.int, ".mu_", mu.int, ".ext.size_0", ".spp.size_", spp.int)
				assign(name.out, int.save)
				}
			}
		}

## lineages with comparatively smaller trait values were more likely to speciate, and lineages with 
## comparatively larger trait values were more likely to go extinct.

	for(mu.int in c(0, 0.4, 0.8)) {
		for(tips.int in c(50, 100, 200)) {
			for(spp.int in c(0.1, 0.5, 2, 5)) {
				int.save <- background.ext(mu=mu.int, extant.tips=tips.int, spp.size.prob=-spp.int, ext.size.prob=spp.int, both.sides=FALSE, n.core=n.cores)
				name.out <- paste0("tips_", tips.int, ".mu_", mu.int, ".ext.size_", spp.int, ".spp.size_", spp.int)
				assign(name.out, int.save)
				}
			}
		}

save(list=ls(), file="background.one.side.R")