## Puttick et al. 2020. "The complex effects of mass extinctions on
## morphological disparity", Evolution.

## Estimate disparity in a series of time bins with disparity estimated using
## median pairwise euclidean distance, median pairwise distance from the root,
## Sum Of Variances, and Sum Of Ranges using the dispRity package.

library(parallel)
library(PETER)
library(dispRity)

## n cores for parallel processing on Mac and Linux machines
  n.cores <- 1

## input file for parallel processing on Mac and Linux machines
	sigma.files <- "tips_50.mu_0.4.severity_0.5.random_0.01.R"

## estimate disparity for timebins and by sub-sampling data
  all.outs <- mclapply(1:length(sigma.files), mc.cores = n.cores, FUN = function(count) {
    setwd("../precooked.output/")
    load(sigma.files[count], verbose = TRUE)
    in.file <- eval(parse(text = gsub(".R", "", sigma.files[count])))
    list.full <- list()
    counter <- 1
    
    # four, six, ten, bins
    for (bin.number in c(2, 4, 5, 8)) {
      # complete, half, tenth, hundreth sampling
      difference.in.metric <- c()
      for (xx in 1:100) {
      	ltt.out.temp <- in.file[[xx]]$ltt
      	mass.time.temp <- in.file[[xx]]$mass.ext.time
      	end.seq <- max(ltt.out.temp[, 1])
      	mass.to.end <- end.seq - mass.time.temp
      	pre.mass <- mass.time.temp - (mass.to.end / bin.number)
      	end.seq <- mass.time.temp + (mass.to.end / bin.number)
      	breaks <- c(pre.mass, mass.time.temp, end.seq)
      	bin.dat <- .bincode(ltt.out.temp[, 1], breaks, include.lowest = TRUE)
      	x2 <- tapply(in.file[[xx]]$node.traits, bin.dat, function(x)
      	  unlist(x))
      	if (is(x2[[1]])[1] == "numeric") {
      	  max.len <- sapply(x2, length)
      	} else {
      	  max.len <- sapply(x2, function(kk)
      	    kk[1])
      	}
      	sample.number <- 1
      	thr <- ceiling(sapply(x2 , length) * sample.number)
      	root.state <- in.file[[xx]]$root.state
      	if (length(x2) == 2) {
      	  x2 <- lapply(1:2, function(u)
      	    x2[[u]][sample(1:(length(x2[[u]])), thr[u])])
      	  x.out <- sapply(x2, function(x) {
      	    if (length(x) > 10000)
      	      x <- x[sample(1:length(x), 10000)]
      	    x.in <- as.matrix(x)
      	    n.col <- ncol(x.in)
      	    if (nrow(x.in) > 2) {
      	      euclid.out <- unlist(dispRity(x.in, metric = c(median, pairwise.dist),
      	        method = "euclidean")$disparity)
      	      if (n.col == 1)
      	        x.in <- cbind(as.matrix(x), 0)
      	      sov.out <- unlist(dispRity(x.in, metric = c(sum, variances))$disparity)
      	      sor.out <- unlist(dispRity(x.in, metric = c(sum, ranges))$disparity)
      	      c(euclid.out, sov.out, sor.out)
      	    } else {
      	      matrix(NA, nrow = 2, ncol = 2)
      	    }
      	  })
      	  temp.count <- apply(x.out, 1, diff)
		  if (!all(is.na(temp.count))) {
		    difference.in.metric <- rbind(difference.in.metric,
		      c(apply(x.out, 1, diff), sapply(x2, length)))
		  } else {
		    difference.in.metric <- rbind(difference.in.metric, rep(NA, 5))
		  }
	    } else {
		  difference.in.metric <- rbind(difference.in.metric, rep(NA, 5))
		}
	  }
	  colnames(difference.in.metric) <- c("euclidean", "sov", "sor", "n_bin1", "n_bin2")
	  list.full[[counter]] <- difference.in.metric
	  counter <- counter + 1
	}
	names(list.full) <- c("four bins", "eight bins", "ten bins", "sixteen bins")
	list.full
	})
	
  names(all.outs) <- sigma.files
  saveRDS(all.outs, "../disparity.analysis/outputs/tips_50.mu_0.4.severity_0.5.random_0.01_bin")