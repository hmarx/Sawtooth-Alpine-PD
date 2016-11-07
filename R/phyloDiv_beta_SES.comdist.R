#http://esapubs.org/archive/ecos/C006/023/R_code.R

#Description

#The code contains the function to estimate the relative contributions of selection, drift,
#dispersal limitation, and homogenizing dispersal to the community assembly of the heath
#vegetation, as well as the function to evaluate the importance of abiotic variables and spatial
#eigenvectors to the abundance of species through a Redundancy Analysis.


#######################################################################


### Code gently shared by James C Stegen with no guarantee of general application
### Original full code upon which the present one is based is stored at https://github.com/stegen/Stegen_etal_ISME_2013
### Presupposes existence of a species abundance matrix and a corresponding phylogenetic tree with names in the exact same order

library(picante)
library(psych)
library(vegan)
library(packfor) # Function forward.sel

## read in species table

comun = read.csv(file.choose(), sep = ";", header=T,row.names=1)
dim(comun); # this gives the dimensions
comun[1:5,1:5]; # this gives a look at the first 5 rows and columns

## read in the phylogeny

phylo = read.tree("arvore.final.txt")
phylo; # a summary of the phylogeny
plot.phylo(phylo,typ="fan"); # a quick plot

## make sure the names on the phylogeny are ordered the same as the names in OTU table

match.phylo.comun = match.phylo.data(phylo, t(comun));
str(match.phylo.comun);

## calculate empirical betaMNTD

bbeta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.comun$data),cophenetic(match.phylo.comun$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,'clbi_betaMNTD_weighted.csv', quote=F)

identical(colnames(match.phylo.comun$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.comun$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# calculate randomized betaMNTD

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.comun$data),ncol(match.phylo.comun$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.comun$data),taxaShuffle(cophenetic(match.phylo.comun$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.comun$data),ncol=ncol(match.phylo.comun$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.comun$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.comun$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.comun$data);
colnames(weighted.bNTI) = colnames(match.phylo.comun$data);
weighted.bNTI;
write.csv(weighted.bNTI,"clbi_weighted_bNTI.csv",quote=F);

pdf("clbi_weighted_bNTI_Histogram.pdf")
hist(weighted.bNTI)
dev.off()


############################



# Community turnover null model employed species abundance matrix using the code provided by
# CHASE, J.M., KRAFT, N.J.B., SMITH, K.G., VELLEND, M. & INOUYE, B.D. 2011.
#  Using null models to disentangle variation in community dissimilarity from variation in a-diversity. Ecosphere 2:art24. 

# Fractions of community assembly processes were calculated in a
#  spreadsheet using the distance matrices the pieces of code written
#  and cited above produced



#############################

# Environmental Relationships

#############################


# Read in overall file with species, environmental, and PCNM data as columns
#  The spreadsheet I was using used ";" instead of "," to separate
#  values in .csv files, hence the "sep = ";" below

all = read.csv(file.choose(), sep = ";", header=T)

# just checking
all[1:5, 1:5]

attach(all)

# data preparation
texture = Siltegkg.1 + Claygkg.1
env.brute = all[,c(6:7,18:19,21:29,32:36)]

# Transformations

env.brute2 = cbind(env.brute, texture)
log.env.brute2 = log(env.brute2 + 1)
log.convex = log(convexity + 8.5)
canopy = asin(sqrt(canopy.openness))
env = cbind(log.env.brute2, log.convex, canopy)
env.pcnm = cbind(env, all[,40:46])
names(env.pcnm)

# Standardization

env.pcnm.stand = decostand(env.pcnm, method = "standardize")

# PCA

pca = principal(env.pcnm.stand, nfactors = 9,rotate = "none", covar = FALSE, scores = TRUE)
pca$values  # eigenvalues
scores = pca$scores  # axes scores for future use
pca

## RDA
# Test and slection of environmental variables 

sp = all[,47:78]

sp.hel = decostand(sp, "hellinger")


rda.env = rda(sp.hel, scores)
anova.cca(rda.env)  # to check the overall explanatory power of pca axes
r2.env = RsquareAdj(rda.env)
r2.env
# forward selection of pca scores using adjusted R2 of rda.env
# as significance criterion
env.fwd = forward.sel(sp.hel, as.matrix(scores), adjR2thresh = r2.env$adj.r.square)
env.fwd  # checking significant pca axes
env.sign = sort(env.fwd$order)  # order of most significant forward-selected pca axes
env.final = env[,c(env.sign)]