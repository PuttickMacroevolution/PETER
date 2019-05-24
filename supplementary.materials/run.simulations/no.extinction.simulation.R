## Puttick et al. Code to simulate diversity and trait data through time with no mass extinction simulation, in 
## simulations with one, two, and five traits generated under Brownian motion and trees generated under 
## homogeneous birth-death processes

library(RSOLE)
library(parallel)

## set speciation rate for all simulations
	lambda.val <- 1
	
## set Brownian motion rate for all simulations
	bm.val <- 1
	
## n cores for parallel processing on Mac and Linux machines
	n.cores <- 1

# single trait

for(mu.val in c(0, 0.4, 0.8)) {
	for(tip.n in c(200, 100, 50)) {
		disp.out <- mclapply(1:100, mc.cores=n.cores, function(u) {
		est.out <- build.tree.tips(lambda=lambda.val, mu=mu.val, extant.tips=tip.n, mass.ext.occ=FALSE, sigma.rates=bm.val)
			}
		)
	name.out <- paste0("tips_", tip.n, ".mu_", mu.val, ".sigma_", bm.val)		
	assign(name.out, disp.out)
	save(list=name.out, file=paste0(name.out, ".R"), compress="xz")			
}

## two traits

# for(co.var in c(0, 0.75)) {
 	# bm.val <- matrix(c(1,co.var,co.var,1), 2, 2)
	# for(mu.val in c(0, 0.4, 0.8)) {
		# for(tip.n in c(200, 100, 50)) {
			# disp.out <- mclapply(1:100, mc.cores=n.cores, function(u) {
				# est.out <- build.tree.tips(lambda=lambda.val, mu=mu.val, n.characters=2, extant.tips=tip.n, mass.ext.occ=FALSE, sigma.rates=bm.val)
				# }
			# )
		# name.out <- paste0("tips_", tip.n, ".two_traits.mu_", mu.val, ".covar_", co.var)
		# assign(name.out, disp.out)
		# save(list=name.out, file=paste0(name.out, ".R"), compress="xz")	
		# }
	# }
# }
	
## five traits

# for(co.var in c(0, 0.75)) {
 	# bm.val <- matrix(co.var, ncol=5, nrow=5); diag(bm.val) <- 1
	# for(mu.val in c(0, 0.4, 0.8)) {
		# for(tip.n in c(200, 100, 50)) {
			# disp.out <- mclapply(1:100, mc.cores=n.cores, function(u) {
				# est.out <- build.tree.tips(lambda=lambda.val, mu=mu.val, n.characters=5, extant.tips=tip.n, mass.ext.occ=FALSE, sigma.rates=bm.val)
				# }
			# )
		# name.out <- paste0("tips_", tip.n, ".five_traits.mu_", mu.val, ".covar_", co.var)
		# assign(name.out, disp.out)
		# save(list=name.out, file=paste0(name.out, ".R"), compress="xz")	
		# }
	# }
# }