###########################################  Phylo Beta Diversity ########################################### 

############################ Turnover of species within Clades for all alpine (talus + meadow)
############ Plotting observed and expected community structure across branching times of a phylogeny
# generate a plot of observed and expected PIst on summit community phylogeny 
# PIst: Bst analogue for presence/absence data, expressing phylogenetic turnover (independently of species turnover/richness).
com.beta <- spacodi.calc(sp.plot = t(sawMiseqBeta$comm), phy = sawMiseqBeta$phy, pairwise = T)
sp.permut=spacodi.by.nodes(sp.plot=com.beta$sp.plot, phy=com.beta$sp.tree, sp.parm="PIst", n.rep=999, method = "1s") #shuffling of abundances across entire dataset
head(sp.permut$observed.PIst)
head(sp.permut$randomization.test)

pdf(file="output/08_PhyloDiversity/MiSeq/beta/nodal/all/PIst.branching.time.pdf")
spacodi.permutplot(sp.permut, bty="n", add.id = F, sig.plot=TRUE)
dev.off()
write.csv(sp.permut$observed.PIst, file="output/08_PhyloDiversity/MiSeq/beta/nodal/all/observed.PIst.csv")
write.csv(sp.permut$expected.PIst, file="output/08_PhyloDiversity/MiSeq/beta/nodal/all/expected.PIst.csv")
write.csv(sp.permut$randomization.test, file="output/08_PhyloDiversity/MiSeq/beta/nodal/all/randomization.test.csv")

#### plotting diversity turnover on trees
pdf(file = "output/08_PhyloDiversity/MiSeq/beta/nodal/all/PIst.phylo.nodes.pdf", 10, 10)
spacodi.treeplot(sp.permut, com.beta$sp.tree, sig.plot=TRUE, add.id=FALSE, type="fan", tip=1, pch=.5)
  dev.off()


#If an ultrametric phylogeny is supplied, Dkk is equivalent to the mean pairwise phylogenetic distance
#(distance to MRCA) between two individuals drawn from a community,
#Dkl is the mean pairwise phylogenetic distance between individuals drawn from each of two communities, 
#and H is Dkl standardized to account for within-community diversity.
raoD(comm = sawMiseqAlp$comm, phy = sawMiseqAlp$phy)

unifrac(comm = sawMiseqAlp$comm, tree = sawMiseqAlp$phy)
phylosor.rnd(comm = sawMiseqAlp$comm, tree = sawMiseqAlp$phy)

comdist(comm = sawMiseqAlp$comm, dis = cophenetic(sawMiseqAlp$phy))


#Dnn (Calculates inter-community mean nearest taxon distance)
comdistnt(comm = sawMiseqAlp$comm, dis = cophenetic(sawMiseqAlp$phy))


############################ Turnover of species within Clades for talus only
# PIst: Bst analogue for presence/absence data, expressing phylogenetic turnover (independently of species turnover/richness).
com.beta.alp <- spacodi.calc(sp.plot = t(sawMiseqBetaAlp$comm), phy = sawMiseqBetaAlp$phy, pairwise = T)
sp.permut.alp=spacodi.by.nodes(sp.plot=com.beta.alp$sp.plot, phy=com.beta.alp$sp.tree, sp.parm="PIst", n.rep=999, method = "1s") #shuffling of abundances across entire dataset
head(sp.permut.alp$observed.PIst)
head(sp.permut.alp$randomization.test)

pdf(file="output/08_PhyloDiversity/MiSeq/beta/nodal/talus/PIst.branching.time.pdf")
spacodi.permutplot(sp.permut, bty="n", add.id = F, sig.plot=TRUE)
dev.off()
write.csv(sp.permut$observed.PIst, file="output/08_PhyloDiversity/MiSeq/beta/nodal/talus/observed.PIst.csv")
write.csv(sp.permut$expected.PIst, file="output/08_PhyloDiversity/MiSeq/beta/nodal/talus/expected.PIst.csv")
write.csv(sp.permut$randomization.test, file="output/08_PhyloDiversity/MiSeq/beta/nodal/talus/randomization.test.csv")

#### plotting diversity turnover on trees
pdf(file = "output/08_PhyloDiversity/MiSeq/beta/nodal/talus/PIst.phylo.nodes.pdf", 10, 10)
spacodi.treeplot(sp.permut, com.beta$sp.tree, sig.plot=TRUE, add.id=FALSE, type="fan", tip=1, pch=.5)
dev.off()



