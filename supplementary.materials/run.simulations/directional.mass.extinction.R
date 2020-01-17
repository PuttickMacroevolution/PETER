## XXX et al.
## Code to simulate diversity and trait data through time with a directional
## mass extinction when 50, 100, or 200 taxa are present.

## The extinction removes a proportion (0.5, 0.75, or 0.9) of contemporary
## species. Taxa are removed according to a strictly selective extinction in
## which all taxa with traits values above a quantile value become extinct
## (become tips). In a strong extinction there is the same quantile, but a 0.01
## chance taxa with traits above the quantile value survive, and 0.01 chance
## those with a value below it are lost. The simulation continues until the pre-
## extinction diversity (50, 100, 200 tips) is recovered. Simulations are
## repeated with one, two, and five traits. For these models, selection acts on
## the first trait and runs are repeated with no co-variance between trait
## evolution and 0.75 trait co-variation. All trees generated under homogeneous
## birth-death processes.

library(PETER)
library(parallel)

## Speciation rate.
	lambda.val <- 1

## Brownian variance.
	bm.val <- 1

## N cores for parallel processing on Mac and Linux machines.
	n.cores <- 1

## Chance of random extinction for strict, strong, and random extinctions.
	random <- c(0, 0.01, 1)

## Number of tips in the simulations.
	tip.number <- c(50, 100, 200)

## Proportion of contemporary taxa that go extinct at an extinction boundary.
	severity <- c(0.5, 0.75, 0.9)

## Background extinction level.
  background.mu <- c(0, 0.4, 0.8, 1)

## Example analysis with 100 simulates generated with lambda = 1, mu = 1,
## 50 extant tips, 25 taxa lost at the mass extinction, with strong selectivity
## in which there is a 0.01 chance of random extinction. 

  set.seed(101)
  disp.out <- mclapply(1:100, mc.cores = n.cores, function(in.lo) {
    est.out <- suppressWarnings(build.tree.tips(lambda = 1, mu = 0.4,
      extant.tips = 50, mass.ext.occ = TRUE, sigma.rates = 1,
      mass.ext.tips = 50, mass.ext.quantile = 1 - 0.5, mass.ext.severity = 0.5,
      chance.of.random = 0.01, trait.axes = 1))
	})
  name.out <- paste0("tips_50.mu_0.4.severity_0.5.random_0.01")
  assign(name.out, disp.out)
  save(list = name.out, file = paste0("../precooked.output/", name.out, ".R"),
    compress = "xz")


## Simulation of trees with a single trait with different levels of background
## extinction, number of tips, severities of extinction, and selectivity of the
## mass extinction event.

# for (mu.val in background.mu) {
#	  for (y in tip.number) {
#		  for (u in severity) {
#			  for (i in random) {
#				  disp.out <- mclapply(1:100, mc.cores = n.cores, function(in.lo) {
#					  est.out <- suppressWarnings(build.tree.tips(lambda = 1, mu = mu.val,
#           extant.tips = y, mass.ext.occ = TRUE, sigma.rates = 1,
#           mass.ext.tips = y, mass.ext.quantile = 1 - u, mass.ext.severity = u,
#           chance.of.random = i, trait.axes = 1))
#					})
#				name.out <- paste0("tips_", y, ".mu_", mu.val,
#         ".severity_", u, ".random_", i)
#				assign(name.out, disp.out)
#				save(list = name.out, file = paste0(name.out, ".R"), compress = "xz")
#				}
#			}
#		}
#	}

## Simulation of trees with a two traits with different levels of background
## extinction, number of tips, severities of extinction, and selectivity of the
## mass extinction event. Also, different simulations are run with two levels
## of covariance in trait evolution (0, 0.75).

# for (co.var in c(0, 0.75)) {
#	  bm.val <- matrix(co.var, ncol = 2, nrow = 2); diag(bm.val) <- 1 ;
#	  for (mu.val in background.mu) {
#		  for (y in tip.number) {
#			  for (u in severity) {
#				  for (i in random) {
#					  disp.out.int <- mclapply(1:100, mc.cores = 8, function(in.lo) {
#						  est.out <- build.tree.tips(lambda = 1, mu = mu.val,
#               extant.tips = y, mass.ext.occ = TRUE, sigma.rates = bm.val,
#               mass.ext.tips = y, mass.ext.quantile = 1-u,
#               mass.ext.severity = u, chance.of.random = i, n.characters = 2,
#               trait.axes = 1)
#						  est.out
#						})
#				    name.out <- paste0("tips_", y, ".two_traits.mu_", mu.val,
#             ".severity_", u, ".random_", i, ".covar_", co.var)
#				    assign(name.out, disp.out)
#				    save(list = name.out, file = paste0(name.out, ".R"),
#             compress = "xz")
#				  }
#			  }
#		  }
#	  }
# }



## Simulation of trees with a five traits with different levels of background
## extinction, number of tips, severities of extinction, and selectivity of the
## mass extinction event. Also, different simulations are run with two levels
## of covariance in trait evolution (0, 0.75).

# for (co.var in c(0, 0.75)) {
#	  bm.val <- matrix(co.var, ncol = 5, nrow = 5); diag(bm.val) <- 1 ;
#	  for (mu.val in background.mu) {
#		  for (y in tip.number) {
#			  for (u in severity) {
#				  for (i in random) {
#					  disp.out.int <- mclapply(1:100, mc.cores = 8, function(in.lo) {
#						  est.out <- build.tree.tips(lambda = 1, mu = mu.val,
#               extant.tips = y, mass.ext.occ = TRUE, sigma.rates = bm.val,
#               mass.ext.tips = y, mass.ext.quantile = 1-u,
#               mass.ext.severity = u, chance.of.random = i, n.characters = 2,
#               trait.axes = 1)
#						  est.out
#						})
#				    name.out <- paste0("tips_", y, ".two_traits.mu_", mu.val,
#             ".severity_", u, ".random_", i, ".covar_", co.var)
#				    assign(name.out, disp.out)
#				    save(list = name.out, file = paste0(name.out, ".R"),
#             compress = "xz")
#				  }
#			  }
#		  }
#	  }
# }
