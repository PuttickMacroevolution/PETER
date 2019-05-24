README
R code to run all simulations and analyse presented in Puttick et al. 

FILE

All example data here use the 'precooked' file tips_50.mu_0.4.severity_0.5.random_0.01.R that is found in the folder 'precooked.output'

summarise.disparity.R
	Summarise all traits under a disparity model with median pairwise Euclidean
	distance, median pairwise distance from the root, Sum Of Variances, and
	Sum Of Ranges

analyse.time.bins.R
	Analyse disparity in separate time bins

phylogenetic.signal.R
	Estimation of the phylogenetic signal (phylo D) and evidence of trait-biased
	extinction (PGLS) of lineages at the mass extinction boundary.

PCM_models.R
	Data pruned of fossils so they represent 'extant' lineages present at the end 
	of the simulation only. These data are analysed using BM, OU, and EB models


