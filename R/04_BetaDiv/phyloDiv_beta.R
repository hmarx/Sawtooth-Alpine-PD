#################################################################################################################
######## Phylo Beta Diversity (nodal) ###########################################################################
######## Hannah E. Marx, 20 Jan 2019 ############################################################################
#################################################################################################################

############################ Turnover of species within Clades for all alpine (talus + meadow)
############ Plotting observed and expected community structure across branching times of a phylogeny
# generate a plot of observed and expected PIst on summit community phylogeny 
# PIst: Bst analogue for presence/absence data, expressing phylogenetic turnover (independently of species turnover/richness).
com.beta.total <- spacodi.calc(sp.plot = t(sawTotal$comm), phy = sawTotal$phy, pairwise = T)
sp.permut.total=spacodi.by.nodes(sp.plot=com.beta.total$sp.plot, phy=com.beta.total$sp.tree, sp.parm="PIst", n.rep=999, method = "1s") #shuffling of abundances across entire dataset
head(sp.permut.total$observed.PIst)
head(sp.permut.total$randomization.test)

pdf(file="output/10_PhyloDiversity/beta/nodal/all/TotalPIst.branching.time.pdf")
spacodi.permutplot(sp.permut.total, bty="n", add.id = F, sig.plot=TRUE)
dev.off()
write.csv(sp.permut.total$observed.PIst, file="output/10_PhyloDiversity/beta/nodal/all/Totalobserved.PIst.csv")
write.csv(sp.permut.total$expected.PIst, file="output/10_PhyloDiversity/beta/nodal/all/Totalexpected.PIst.csv")
write.csv(sp.permut.total$randomization.test, file="output/10_PhyloDiversity/beta/nodal/all/Totalrandomization.test.csv")

#### plotting diversity turnover on trees
pdf(file = "output/10_PhyloDiversity/beta/nodal/all/TotalPIst.phylo.nodes.pdf", 10, 10)
spacodi.treeplot(sp.permut.total, com.beta.total$sp.tree, sig.plot=TRUE, add.id=FALSE, type="fan", tip=1, pch=.5)
dev.off()


############################ Turnover of species within Clades for talus only
# PIst: Bst analogue for presence/absence data, expressing phylogenetic turnover (independently of species turnover/richness).
com.beta.alp.total <- spacodi.calc(sp.plot = t(sawTotal.alpine$comm), phy = sawTotal.alpine$phy, pairwise = T)
sp.permut.alp.total=spacodi.by.nodes(sp.plot=com.beta.alp.total$sp.plot, phy=com.beta.alp.total$sp.tree, sp.parm="PIst", n.rep=999, method = "1s") #shuffling of abundances across entire dataset
head(sp.permut.alp.total$observed.PIst)
head(sp.permut.alp.total$randomization.test)

pdf(file="output/10_PhyloDiversity/beta/nodal/talus/Total.talus.PIst.branching.time.pdf")
spacodi.permutplot(sp.permut.alp.total, bty="n", add.id = F, sig.plot=TRUE)
dev.off()
write.csv(sp.permut.alp.total$observed.PIst, file="output/10_PhyloDiversity/beta/nodal/talus/Total.talus.observed.PIst.csv")
write.csv(sp.permut.alp.total$expected.PIst, file="output/10_PhyloDiversity/beta/nodal/talus/Total.talus.expected.PIst.csv")
write.csv(sp.permut.alp.total$randomization.test, file="output/10_PhyloDiversity/beta/nodal/talus/Total.talus.randomization.test.csv")

#### plotting diversity turnover on trees
pdf(file = "output/10_PhyloDiversity/beta/nodal/talus/Total.talus.PIst.phylo.nodes.pdf", 10, 10)
spacodi.treeplot(sp.permut.alp.total, com.beta.alp.total$sp.tree, sig.plot=TRUE, add.id=FALSE, type="fan", tip=1, pch=.5)
dev.off()


