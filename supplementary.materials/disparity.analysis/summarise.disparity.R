## XXX et al.
## Estimate disparity for generated data using median pairwise euclidean
## distance, median pairwise distance from the root, Sum Of Variances, and
## Sum Of Ranges using the dispRity package.

library(parallel)
library(PETER)
library(dispRity)

## files to summarise
setwd("../precooked.output/")
input.file <- list.files(pattern = "tips_")

## n cores for parallel processing on Mac and Linux machines.
	n.cores <- 1

## For multiple files uncomment the lines below.
# all.outs <- mclapply(1:length(sigma.files),
#   mc.cores = n.cores, FUN = function(count) {

# Single file example.
load(input.file, verbose = TRUE)
in.file <- eval(parse(text = gsub(".R", "", input.file)))
ltt.out <- list()

for (xx in 1:100) {
  ltt.out.temp <- in.file[[xx]]$ltt
  going <- c()
  root.state <- in.file[[xx]]$root.state
  x.out <- sapply(in.file[[xx]]$node.traits, function(x) {
    x.in <- as.matrix(x)
	n.col <- ncol(x.in)
	if (nrow(x.in) > 2) {
	  euclid.out <- unlist(dispRity(x.in,
	    metric = c(median, pairwise.dist), method = "euclidean")$disparity)
	  if (n.col == 1)
	    x.in <- cbind(as.matrix(x), 0)
	  sov.out <- unlist(dispRity(x.in, metric = c(sum, variances))$disparity)
	  sor.out <- unlist(dispRity(x.in, metric = c(sum, ranges))$disparity)
	  cbind(euclid.out, sov.out, sor.out)
    } else {
	  cbind(rep(NA, 3))
	}
  })
  x.out <- t(x.out)
  node.info <- cbind(ltt.out.temp, x.out)
  colnames(node.info)[-c(1:2)] <- c("median.euclidean", "sov", "sor")
  ltt.out[[xx]] <- node.info
}

ltt.write <- gsub(".R", ".dvltt", input.file)
assign(ltt.write, ltt.out)
save(list = ltt.write, file =
  paste0("../disparity.analysis/outputs/", gsub(".R", ".dvltt.R", input.file)),
	compress = "xz")
