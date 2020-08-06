## Puttick et al. 2020. "The complex effects of mass extinctions on
## morphological disparity", Evolution.

## Code to simulate diversity and trait data through time with no
## mass extinction simulation, in simulations with one, two, and five traits
## generated under Brownian motion and trees generated under homogeneous
## birth-death processes.

library(PETER)
library(parallel)

## Set speciation rate for all simulations.
	lambda.val <- 1

## Set Brownian motion variance for all simulations.
	bm.val <- 1

## N cores for parallel processing on Mac and Linux machines.
	n.cores <- 1

## Background extinction level.
  background.mu <- c(0, 0.4, 0.8)

## Number of tips to simulate.
  n.taxa <- c(50, 100, 200)

# Simulations with a single trait generated under Brownian motion.

for (mu.val in background.mu) {
	for (tip.n in n.taxa) {
		disp.out <- mclapply(1:100, mc.cores = n.cores, function(u) {
		  est.out <- build.tree.tips(lambda = lambda.val, mu = mu.val,
		    extant.tips = tip.n, mass.ext.occ = FALSE, sigma.rates = bm.val)
		})
	name.out <- paste0("tips_", tip.n, ".mu_", mu.val, ".sigma_", bm.val)
	assign(name.out, disp.out)
	save(list = name.out, file = paste0(name.out, ".R"), compress = "xz")
}

# Simulations with a two traits generated under Brownian motion with either
# absence (0) or presence of co-variance (0.75) between the traits.

# for (co.var in c(0, 0.75)) {
#   bm.val <- matrix(c(1, co.var, co.var, 1), 2, 2)
#	  for (mu.val in background.mu) {
#		  for (tip.n in n.taxa) {
#			  disp.out <- mclapply(1:100, mc.cores = n.cores, function(u) {
#				  est.out <- build.tree.tips(lambda = lambda.val, mu = mu.val,
#           n.characters = 2, extant.tips = tip.n, mass.ext.occ = FALSE,
#           sigma.rates = bm.val)
#				})
#		  name.out <- paste0("tips_", tip.n, ".two_traits.mu_", mu.val, ".covar_",
#       co.var)
#		  assign(name.out, disp.out)
#		  save(list = name.out, file = paste0(name.out, ".R"), compress = "xz")
#		  }
#	  }
# }

# Simulations with a five traits generated under Brownian motion with either
# absence (0) or presence of co-variance (0.75) between the traits.


# for(co.var in c(0, 0.75)) {
#   bm.val <- matrix(co.var, ncol = 5, nrow = 5); diag(bm.val) <- 1
#	  for (mu.val in background.mu) {
#		  for (tip.n in n.taxa) {
#			  disp.out <- mclapply(1:100, mc.cores = n.cores, function(u) {
#				  est.out <- build.tree.tips(lambda = lambda.val, mu = mu.val,
#         n.characters = 5, extant.tips = tip.n, mass.ext.occ = FALSE,
#         sigma.rates = bm.val)
#				})
#		  name.out <- paste0("tips_", tip.n, ".five_traits.mu_", mu.val,
#       ".covar_", co.var)
#		  assign(name.out, disp.out)
#		  save(list=name.out, file=paste0(name.out, ".R"), compress="xz")
#		  }
#	  }
# }
