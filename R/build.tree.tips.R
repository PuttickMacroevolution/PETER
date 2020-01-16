#' @title Simulate simultaneously phylogeny and trait evolution
#' @description Simulates a birth-death phylogeny with simultaneous evolution under Brownian motion for one or more traits. The simulations continue until a user-set number of tips are sampled. Additionaly, a mass extinction can occur with a user-set intensity of extinction when a certain number of species are present. Traits can be correlated or independent, can be linked to increased rates of speciation and/or extinction at background times, and can drive selectivity at mass extinction event.
#' @param lambda Speciation rate.
#' @param mu Extinction rate. If mu=0, then the simulation follows a Yule pure birth process.
#' @param extant.tips The number of extant lineages with which to end the simulation. Must be > 0.
#' @param mass.ext.tips The number of contemporary lineages the simulation reaches before a mass extinction. Not applicable if mass.ext.occ = FALSE.
#' @param mass.ext.severity The proportion of lineages that succumb to extinction, not applicable if mass.ext.occ=FALSE.
#' @param n.characters The number of trait characters to simulate under a Brownian motion process.
#' @param sigma.rates Either a vector or matrix specifying the rates of character evolution. If a vector then all traits evolve independently, and values are recycled if the length is less than n.characters. If a square matrix, this specifies the relative rates (diagonals) and correlation between traits (off-diagonals).
#' @param spp.size.prob Probability of lineages with smaller trait values splitting during 'background' times. If the value is 0 then all lineages have an equal probability of splitting; larger values mean lineages with smaller trait values are more likely to split, and values < 0 mean lineages with larger trait values are more likely to split.
#' @param ext.size.prob Probability of lineages with larger trait values going extinct during 'background' times. If the value is 0 then all lineages have an equal probability of extinction; larger values mean lineages with relatively bigger trait values are more likely to go extinct, and values < 0 mean lineages with smaller trait values are more likely to go extinct.
#' @param mass.ext.occ Logical specifying if a mass extinction occurs at the time given by time.limit.mass.ext (default = TRUE).
#' @param mass.ext.quantile The cut-off for the trait values selected for extinction during the mass extinction event. Values above mass.ext.quantile will be marked for loss during the mass extinction event.
#' @param chance.of.random Numeric value to specify the selectivity of extinction. If chance.of.random = 0, extinction is strict and only lineages with trait values above or below the limit specified by mass.ext.quantile will go extinct. Values > 0 and < 1 indicate the trait is selected for with decreasing strength, until chance.of.random = 1 at which point all lineages are equally likely to go extinct.
#' @param larger.than Logical. If TRUE, then mass extinction selection occurs on traits larger than the proportion indicated by 'mass.ext.quantile'. If FALSE, the smaller traits are selected for extinction.
#' @param both.sides Logical. If TRUE, then mass extinction selection occurs on at the extreme ends of trait distribution. NOTE - this over-rides the argument 'larger.than' so the default is FALSE.
#' @param trait.axes The trait axes on which selection will occur. The default (NULL) acts on all trait axes.
#' @param return.trait.matrix Logical. If TRUE, the values at internal and external nodes is returned in a matrix corresponding to the APE phylo edge format.
#' @param return.branch.value.list Logical. Returns a list with all trait values at every time point. Default is FALSE as TRUE can results in huge files.
#' @param return.time.list Logical. Returns a list with every time points at which the simulation ran. Default is FALSE as TRUE can results in huge files.
#' @param upperBound The upper bound allowed for trait values. Trait space can be constrained by specifying a value rather than then default (Inf).
#' @param lowerBound The lower bound allowed for trait values. Trait space can be constrained by specifying a value rather than then default (-Inf).
#' @return A list containing:
#' \itemize{
#' \item {phy} {The output phylogeny}
#' \item {root.state} {value(s) for trait(s) at the root}
#' \item {trait.matrix} {List of trait values at all internal nodes and tips on the phylogeny. If n.characters > 1, then there are list entries for each trait value}
#' \item {branch.value.list.out} {The value of traits at each time unit in the simulation. Each list element corresponds to the edge in the phy$edge matrix, and the time given by time.list.out}
#' \item {vtt} {The variance measured at each node of the phylogeny to give 'variance-through-time' data. Only returned if return.through.time = TRUE}
#' \item {ltt} {The number of lineages measured at each node of the phylogeny to give 'lineage-through-time' data. Only returned if return.through.time = TRUE}
#' \item {vtt.extant} {The variance measured at each node of the phylogeny to give 'variance-through-time' data based on lineages that survive to present only. Only returned if return.through.time.extant = TRUE}
#' \item {ltt.extant} {The number of lineages measured at each node of the phylogeny to give 'lineage-through-time' data based on lineages that survive to present only. Only returned if return.through.time.extant = TRUE}
#' }
#' @author Mark Puttick
#' @import ape
#' @import MASS
#' @import motmot
#' @export

build.tree.tips <- function(lambda = 1, mu = 0, extant.tips = 100,
  mass.ext.tips = 50, mass.ext.quantile = 0.5, n.characters = 1,
  sigma.rates = 1, spp.size.prob = 0, ext.size.prob = 0, larger.than = TRUE,
  both.sides = FALSE, mass.ext.occ = TRUE, mass.ext.severity = 0.5,
  trait.axes = NULL, chance.of.random = 0.1, upperBound = Inf,
  lowerBound = -Inf, return.trait.matrix = TRUE,
  return.branch.value.list = FALSE, return.time.list = FALSE)
  {

  tree.alive <- FALSE
  total.time <- 0
  open.edges <- 2
  mass.ext.occ.input <- mass.ext.occ
  if (lambda < 1e-08)
    stop("speciation rates must be greater than zero")
  if (!is.matrix(sigma.rates)) {
    bm.rates <- diag(sigma.rates, n.characters, n.characters)
  } else {
    bm.rates <- sigma.rates
  }
  while (any(c(length(open.edges) != extant.tips, mass.ext.occ))) {
  if (tree.alive == FALSE) {
    start.edge <- matrix(c(1, 1, NA, NA), 2, 2)
    edge.length.list <- branch.value.list <- time.list <- list()
    open.edges <- which(is.na(start.edge[, 2]))
    edge.lengths <- c(0, 0)
    n.tips <- total.time <- 0
    root.state <- mvrnorm(n = 1, mu = rep(0, n.characters), Sigma = bm.rates)
    withinBounds <- FALSE
    if (all(c(root.state <= upperBound, root.state >= lowerBound)))
      withinBounds <- TRUE
    while (withinBounds == FALSE) {
      root.state <- mvrnorm(n = 1, mu = rep(0, n.characters), Sigma = bm.rates)
      if (all(c(root.state <= upperBound, root.state >= lowerBound)))
        withinBounds <- TRUE
    }
    branch.value.list[[1]] <- branch.value.list[[2]] <- rbind(root.state)
    time.list[[1]] <- time.list[[2]] <- total.time
    tree.alive <- TRUE
    mass.ext.occ <- mass.ext.occ.input
  }
  N <- length(open.edges)
  waiting.time <- rexp(1, N * (lambda + mu))
  uniform.dist <- runif(1, 0, 1)
  prop.lambda <- lambda / (lambda + mu)
  total.time <- total.time + waiting.time
  spp.or.ext <- uniform.dist <= prop.lambda
  edge.lengths[open.edges] <- edge.lengths[open.edges] + waiting.time
  branch.value.list[open.edges] <- lapply(branch.value.list[open.edges],
    function(branch.vals) {
      last.time <- branch.vals[dim(branch.vals)[1], ]
    withinBounds <- FALSE
    while (withinBounds == FALSE) {
      branch.vals.temp <- last.time + mvrnorm(n = 1, mu = rep(0, n.characters),
        Sigma = waiting.time * bm.rates)
      if (all(c(branch.vals.temp <= upperBound,
        branch.vals.temp >= lowerBound)))
        withinBounds <- TRUE
    }
    rbind(branch.vals, branch.vals.temp)
  })
  time.list[open.edges] <- lapply(time.list[open.edges], function(x)
    c(x, total.time))
  n.tips <- length(which(start.edge[, 2] > 0))
  if (spp.or.ext) {
    free.edges.values <- sapply(branch.value.list[open.edges], function(x)
      x[dim(x)[1]])
    if (both.sides) {
      diff.mean.val <- abs(free.edges.values - mean(free.edges.values))
      spp.chance.rank <- rank(diff.mean.val)
      if (spp.size.prob < 0)
        spp.chance.rank <- (max(spp.chance.rank) + 1) - spp.chance.rank
      spp.chance <- spp.chance.rank ^ abs(spp.size.prob)
    } else {
      spp.value.order <- rank(free.edges.values)
      if (spp.size.prob < 0)
        spp.value.order <- (max(spp.value.order) + 1) - spp.value.order
      spp.chance <- spp.value.order ^ abs(spp.size.prob)
    }
    prob.of.spp <- spp.chance / sum(spp.chance)
    spp.edge <- sample(1:length(open.edges), 1, prob = prob.of.spp)
    new.label <- max(start.edge, na.rm = TRUE) + 1
    start.edge[open.edges[spp.edge], 2] <- new.label
    start.edge <- rbind(start.edge,
      matrix(c(new.label, new.label, NA, NA), 2, 2))
    edge.lengths <- c(edge.lengths, 0, 0)
    node.temp <- tail(branch.value.list[[open.edges[spp.edge]]], 1)
    branch.value.list[[length(branch.value.list) + 1]] <- node.temp
    branch.value.list[[length(branch.value.list) + 1]] <- node.temp
    time.list[[length(time.list) + 1]] <- total.time
    time.list[[length(time.list) + 1]] <- total.time
  } else {
    free.edges.values <- sapply(branch.value.list[open.edges], function(x)
      x[dim(x)[1]])
    if (both.sides) {
      diff.mean.val <- abs(free.edges.values - mean(free.edges.values))
      ext.chance <- rank(diff.mean.val)
      if (ext.size.prob > 0)
        ext.chance <- (max(ext.chance) + 1) - ext.chance
        ext.chance <- ext.chance ^ abs(ext.size.prob)
      } else {
        ext.risk.order <- rank(free.edges.values)
        if (ext.size.prob < 0)
          ext.risk.order <- (max(ext.risk.order) + 1) - ext.risk.order
        ext.chance <- ext.risk.order ^ abs(ext.size.prob)
      }
   	prob.of.ext <- ext.chance / sum(ext.chance)
    ext.edge <- sample(1:length(open.edges), 1, prob = prob.of.ext)
    start.edge[open.edges[ext.edge], 2] <- 0
  }
  tips <- length(which(start.edge[, 2] == 0))
  n.tips.t <- tips
  if ((n.tips.t * 2 - 2) == dim(start.edge)[1])
    tree.alive <- FALSE
  open.edges <- which(is.na(start.edge[, 2]))
  if (length(open.edges) == (mass.ext.tips) && mass.ext.occ) {
    N <- length(open.edges)
    waiting.time <- rexp(1, N * (lambda + mu))
    total.time <- total.time + waiting.time
    edge.lengths[open.edges] <- edge.lengths[open.edges] + waiting.time
    branch.value.list[open.edges] <- lapply(branch.value.list[open.edges],
      function(branch.vals) {
        last.time <- branch.vals[dim(branch.vals)[1], ]
        withinBounds <- FALSE
        while (withinBounds == FALSE) {
          branch.vals.temp <- last.time + mvrnorm(n = 1,
            mu = rep(0, n.characters), Sigma = waiting.time * bm.rates)
          if (all(c(branch.vals.temp <= upperBound,
            branch.vals.temp >= lowerBound)))
              withinBounds <- TRUE
        }
        rbind(branch.vals, branch.vals.temp)
    })
    time.list[open.edges] <- lapply(time.list[open.edges], function(x)
      c(x, total.time))
    free.edges.values <- t(sapply(branch.value.list[open.edges], function(x)
      x[dim(x)[1], ]))
    if (n.characters == 1)
      free.edges.values <- t(free.edges.values)
      if (!is.null(trait.axes)) {
        trait.sel <- as.matrix(free.edges.values[, trait.axes])
      } else {
        trait.sel <- as.matrix(free.edges.values)
      }
      if (!both.sides) {
        if (larger.than) {
          ext.casualty <- apply(trait.sel, 2, function(x) which(x >=
            quantile(x, mass.ext.quantile)))
          ext.casualty <- unique(c(ext.casualty))
        } else {
          ext.casualty <- apply(trait.sel, 2, function(x)
            which(x <= quantile(x, mass.ext.quantile)))
          ext.casualty <- unique(c(ext.casualty))
        }
      } else {
        mass.ext.quantile <- mass.ext.quantile / 2
        ext.casualty <- apply(trait.sel, 2, function(x) {
          intersect(which(x >= quantile(x, mass.ext.quantile)),
            which(x <= quantile(x, (1 - mass.ext.quantile))))
        })
      }
      length.strict <- length(open.edges[ext.casualty])
      extinct.check <- floor(mass.ext.severity * length(open.edges))
      if (length(ext.casualty) > extinct.check) {
        ext.casualty <- sample(ext.casualty,
          floor(mass.ext.severity * length(open.edges)))
      }
      length.not.strict <- length(open.edges[-ext.casualty])
      vector.of.probs <- c(rep(1, length.strict),
        rep(chance.of.random, length.not.strict))
      vector.of.probs <- vector.of.probs / sum(vector.of.probs)
      if (chance.of.random <= 1e-08) {
        rand.ext <- 1:length(vector.of.probs)
        non.casualty <- rand.ext[-ext.casualty]
        ext.casualty <- sample(c(ext.casualty, non.casualty), length.strict,
          prob = vector.of.probs)
      } else {
        rand.ext <- 1:length(vector.of.probs)
        non.casualty <- rand.ext[-ext.casualty]
        to.go <- floor(mass.ext.severity * length(open.edges))
        ext.casualty <- sample(c(ext.casualty, non.casualty), to.go,
          prob = vector.of.probs)
      }
      start.edge[open.edges[ext.casualty], 2] <- 0
      mass.ext.occ <- FALSE
      mass.ext.time <- total.time
      open.edges <- which(is.na(start.edge[, 2]))
      tips <- length(which(start.edge[, 2] == 0))
      n.tips.t <- tips
      if ((n.tips.t * 2 - 2) == dim(start.edge)[1])
        tree.alive <- FALSE
    }
  }
  N <- length(open.edges)
  waiting.time <- rexp(1, N * (lambda + mu))
  total.time <- total.time + waiting.time
  edge.lengths[open.edges] <- edge.lengths[open.edges] + waiting.time
  branch.value.list[open.edges] <- lapply(branch.value.list[open.edges],
    function(branch.vals) {
      last.time <- branch.vals[dim(branch.vals)[1], ]
      withinBounds <- FALSE
      while (withinBounds == FALSE) {
        branch.vals.temp <- last.time +
          mvrnorm(n = 1, mu = rep(0, n.characters),
          Sigma = waiting.time * bm.rates)
        if (all(c(branch.vals.temp <= upperBound,
          branch.vals.temp >= lowerBound)))
            withinBounds <- TRUE
      }
      rbind(branch.vals, branch.vals.temp)
    })
  time.list[open.edges] <- lapply(time.list[open.edges], function(x)
    c(x, total.time))
  trait.list <- lapply(1:n.characters, function(all.traits) {
    anc.state <- t(sapply(branch.value.list, function(x) {
      x.one <- x[1, all.traits]
      x.two <- dim(x)[1]
      c(x.one, x[x.two, all.traits])
    }))
    colnames(anc.state) <- c("internal", "external")
    anc.state
  })
  tips <- c(which(start.edge[, 2] == 0), which(is.na(start.edge[, 2])))
  n.tips <- length(tips)
  edge.dim <- dim(start.edge)[1]
  if (((n.tips * 2 - 2) != edge.dim))
    warning("this is not a tree")
  continuous.seq <- match(seq(1:max(start.edge, na.rm = T)), c(start.edge))
  cont.na <- which(is.na(c(continuous.seq)))
  if (length(cont.na) > 0) {
    edge.one <- which(start.edge[, 1] > cont.na)
    start.edge[edge.one, 1] <- start.edge[edge.one, 1] - 1
    edge.two <- which(start.edge[, 2] > cont.na)
    start.edge[edge.two, 2] <- start.edge[edge.two, 2] - 1
  }
  start.edge[, 1] <- start.edge[, 1] + n.tips
  start.edge[-tips, 2] <- start.edge[-tips, 2] + n.tips
  start.edge[tips, 2] <- 1:n.tips
  start.edge[which(start.edge[, 2] <= n.tips), 2] <- 1:n.tips
  orig.outer <- start.edge[which(start.edge[, 2] > n.tips), 2]
  start.edge[which(start.edge[, 2] > n.tips), 2] <- sort(orig.outer)
  start.edge.orig.2 <- start.edge[, 1]
  for (p in 1:length(sort(orig.outer)))
    start.edge[, 1][which(orig.outer[p] == start.edge.orig.2)] <-
      sort(orig.outer)[p]
  phy <- list(edge = start.edge, tip.label = paste0("t", 1:n.tips),
    edge.length = edge.lengths, Nnode = n.tips - 1)
  class(phy) <- "phylo"
  re.order <- reorder(phy, order = "cladewise", index.only = TRUE)
  tree.object <- list()
  tree.object$phy <- read.tree(text = write.tree(phy))
  tree.object$root.state <- root.state
  if (return.trait.matrix)
    tree.object$trait.matrix <- lapply(trait.list, function(x) x[re.order, ])
  if (return.branch.value.list)
    tree.object$branch.value.list.out <- branch.value.list[re.order]
  if (return.time.list)
    tree.object$time.list.out <- time.list[re.order]
  time.ltt <- unlist(lapply(time.list[re.order], function(xx) xx[-1]))
  tab.time.ltt <- table(signif(time.ltt, 10))
  ltt.times <- as.numeric(names(tab.time.ltt))
  tot.ln <- length(ltt.times)
  re.scale <- diff(c(1, min(as.numeric(tab.time.ltt))))
  ltt.measure <- cbind(node.times = c(0, ltt.times),
    lineages = c(0, as.numeric(tab.time.ltt)))
  tree.object$ltt <- ltt.measure
  unlist.time <- unlist(time.list[re.order])
  trait.ltt <- lapply(signif(ltt.measure[, 1]), function(t.l) {
    time.point.node <- which(signif(unlist.time) == t.l)
    traits.out <- c()
    for (q in 1:n.characters) {
      unlist.traits <- unlist(sapply(branch.value.list[re.order], function(x)
        x[, q]))
      traits.out <- cbind(traits.out, unique(unlist.traits[time.point.node]))
    }
    traits.out
  })
  tree.object$node.traits <- trait.ltt
  if (mass.ext.occ.input)
    tree.object$mass.ext.time <- mass.ext.time
  return(tree.object)
}
