#####################################################################################################################
############# Phylogenetic diversity within alpine summits (alpha) ################################################## 
############# Static Null Model ##################################################################################### 
############# Hannah E. Marx, 20 Jan 2019 ########################################################################### 
#####################################################################################################################

source("R/functions/chooseClade.R")

## Random resample from phylogeny pool, equal probability random draw from phylogeny pool 
## Pool = All Alpine Summits

# Tracheophyta
saw.total.sesmpd.phylonull <- ses.mpd(sawTotal$comm, cophenetic.phylo(sawTotal$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
saw.total.sesmntd.phylonull <- ses.mntd(sawTotal$comm, cophenetic.phylo(sawTotal$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

# Pruned phylogeny and community matrix for each clade 
tax=read.csv(file="data/fleshed_genera.csv", as.is=TRUE, row=1) 

sawtooth.asterales <- pruneCladeTaxonomyLookup(tip.labels = sawTotal$phy$tip.label, tax, level = "order", taxonomy = "Asterales")
pruned.tree.saw.total.aster <-drop.tip(sawTotal$phy, sawTotal$phy$tip.label[!sawTotal$phy$tip.label %in% sawtooth.asterales])
#37
saw.total.aster <- treedata(pruned.tree.saw.total.aster, t(sawTotal$comm))
astertotal.sesmpd.phylonull <- ses.mpd(t(saw.total.aster$data), cophenetic.phylo(saw.total.aster$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
astertotal.sesmntd.phylonull <- ses.mntd(t(saw.total.aster$data), cophenetic.phylo(saw.total.aster$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = sawTotal$phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
pruned.tree.saw.total.caryo <-drop.tip(sawTotal$phy, sawTotal$phy$tip.label[!sawTotal$phy$tip.label %in% sawtooth.caryophyllales])
#19
saw.total.caryo <- treedata(pruned.tree.saw.total.caryo, t(sawTotal$comm))
caryototal.sesmpd.phylonull <- ses.mpd(t(saw.total.caryo$data), cophenetic.phylo(saw.total.caryo$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
caryototal.sesmntd.phylonull <- ses.mntd(t(saw.total.caryo$data), cophenetic.phylo(saw.total.caryo$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.poales <- pruneCladeTaxonomyLookup(tip.labels = sawTotal$phy$tip.label, tax, level = "order", taxonomy = "Poales")
pruned.tree.saw.total.poales <-drop.tip(sawTotal$phy, sawTotal$phy$tip.label[!sawTotal$phy$tip.label %in% sawtooth.poales])
#19
saw.total.poales <- treedata(pruned.tree.saw.total.poales, t(sawTotal$comm))
poalestotal.sesmpd.phylonull <- ses.mpd(t(saw.total.poales$data), cophenetic.phylo(saw.total.poales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
poalestotal.sesmntd.phylonull <- ses.mntd(t(saw.total.poales$data), cophenetic.phylo(saw.total.poales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.lamiales <- pruneCladeTaxonomyLookup(tip.labels = sawTotal$phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
pruned.tree.saw.total.lamiales <-drop.tip(sawTotal$phy, sawTotal$phy$tip.label[!sawTotal$phy$tip.label %in% sawtooth.lamiales])
#12
saw.total.lamiales <- treedata(pruned.tree.saw.total.lamiales, t(sawTotal$comm))
lamialestotal.sesmpd.phylonull <- ses.mpd(t(saw.total.lamiales$data), cophenetic.phylo(saw.total.lamiales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
lamialestotal.sesmntd.phylonull <- ses.mntd(t(saw.total.lamiales$data), cophenetic.phylo(saw.total.lamiales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.brassicales <- pruneCladeTaxonomyLookup(tip.labels = sawTotal$phy$tip.label, tax, level = "order", taxonomy = "Brassicales")
pruned.tree.saw.total.brassicales <-drop.tip(sawTotal$phy, sawTotal$phy$tip.label[!sawTotal$phy$tip.label %in% sawtooth.brassicales])
#9
saw.total.brassicales <- treedata(pruned.tree.saw.total.brassicales, t(sawTotal$comm))
brassicalestotal.sesmpd.phylonull <- ses.mpd(t(saw.total.brassicales$data), cophenetic.phylo(saw.total.brassicales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
brassicalestotal.sesmntd.phylonull <- ses.mntd(t(saw.total.brassicales$data), cophenetic.phylo(saw.total.brassicales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.ericales <- pruneCladeTaxonomyLookup(tip.labels = sawTotal$phy$tip.label, tax, level = "order", taxonomy = "Ericales")
pruned.tree.saw.total.ericales <-drop.tip(sawTotal$phy, sawTotal$phy$tip.label[!sawTotal$phy$tip.label %in% sawtooth.ericales])
#8
saw.total.ericales <- treedata(pruned.tree.saw.total.ericales, t(sawTotal$comm))
ericalestotal.sesmpd.phylonull <- ses.mpd(t(saw.total.ericales$data), cophenetic.phylo(saw.total.ericales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
ericalestotal.sesmntd.phylonull <- ses.mntd(t(saw.total.ericales$data), cophenetic.phylo(saw.total.ericales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

##### Pool = only species found on talus

# Tracheophyta
saw.total.sesmpd.phylonull.talus <- ses.mpd(sawTotal.alpine$comm, cophenetic.phylo(sawTotal.alpine$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
saw.total.sesmntd.phylonull.talus <- ses.mntd(sawTotal.alpine$comm, cophenetic.phylo(sawTotal.alpine$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.asterales.talus <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.alpine$phy$tip.label, tax, level = "order", taxonomy = "Asterales")
pruned.tree.saw.total.aster.talus <-drop.tip(sawTotal.alpine$phy, sawTotal.alpine$phy$tip.label[!sawTotal.alpine$phy$tip.label %in% sawtooth.asterales])
saw.total.aster.talus <- treedata(pruned.tree.saw.total.aster.talus, t(sawTotal.alpine$comm))
astertotal.sesmpd.phylonull.talus <- ses.mpd(t(saw.total.aster.talus$data), cophenetic.phylo(saw.total.aster.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
astertotal.sesmntd.phylonull.talus <- ses.mntd(t(saw.total.aster.talus$data), cophenetic.phylo(saw.total.aster.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.caryophyllales.talus <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.alpine$phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
pruned.tree.saw.total.caryo.talus <-drop.tip(sawTotal.alpine$phy, sawTotal.alpine$phy$tip.label[!sawTotal.alpine$phy$tip.label %in% sawtooth.caryophyllales])
saw.total.caryo.talus <- treedata(pruned.tree.saw.total.caryo.talus, t(sawTotal.alpine$comm))
caryototal.sesmpd.phylonull.talus <- ses.mpd(t(saw.total.caryo.talus$data), cophenetic.phylo(saw.total.caryo.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
caryototal.sesmntd.phylonull.talus <- ses.mntd(t(saw.total.caryo.talus$data), cophenetic.phylo(saw.total.caryo.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.poales.talus <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.alpine$phy$tip.label, tax, level = "order", taxonomy = "Poales")
pruned.tree.saw.total.poales.talus <-drop.tip(sawTotal.alpine$phy, sawTotal.alpine$phy$tip.label[!sawTotal.alpine$phy$tip.label %in% sawtooth.poales])
saw.total.poales.talus <- treedata(pruned.tree.saw.total.poales.talus, t(sawTotal.alpine$comm))
poalestotal.sesmpd.phylonull.talus <- ses.mpd(t(saw.total.poales.talus$data), cophenetic.phylo(saw.total.poales.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
poalestotal.sesmntd.phylonull.talus <- ses.mntd(t(saw.total.poales.talus$data), cophenetic.phylo(saw.total.poales.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.lamiales.talus <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.alpine$phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
pruned.tree.saw.total.lamiales.talus <-drop.tip(sawTotal.alpine$phy, sawTotal.alpine$phy$tip.label[!sawTotal.alpine$phy$tip.label %in% sawtooth.lamiales])
saw.total.lamiales.talus <- treedata(pruned.tree.saw.total.lamiales.talus, t(sawTotal.alpine$comm))
lamialestotal.sesmpd.phylonull.talus <- ses.mpd(t(saw.total.lamiales.talus$data), cophenetic.phylo(saw.total.lamiales.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
lamialestotal.sesmntd.phylonull.talus <- ses.mntd(t(saw.total.lamiales.talus$data), cophenetic.phylo(saw.total.lamiales.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.brassicales.talus <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.alpine$phy$tip.label, tax, level = "order", taxonomy = "Brassicales")
pruned.tree.saw.total.brassicales.talus <-drop.tip(sawTotal.alpine$phy, sawTotal.alpine$phy$tip.label[!sawTotal.alpine$phy$tip.label %in% sawtooth.brassicales])
saw.total.brassicales.talus <- treedata(pruned.tree.saw.total.brassicales.talus, t(sawTotal.alpine$comm))
brassicalestotal.sesmpd.phylonull.talus <- ses.mpd(t(saw.total.brassicales.talus$data), cophenetic.phylo(saw.total.brassicales.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
brassicalestotal.sesmntd.phylonull.talus <- ses.mntd(t(saw.total.brassicales.talus$data), cophenetic.phylo(saw.total.brassicales.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.ericales.talus <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.alpine$phy$tip.label, tax, level = "order", taxonomy = "Ericales")
pruned.tree.saw.total.ericales.talus <-drop.tip(sawTotal.alpine$phy, sawTotal.alpine$phy$tip.label[!sawTotal.alpine$phy$tip.label %in% sawtooth.ericales])
saw.total.ericales.talus <- treedata(pruned.tree.saw.total.ericales.talus, t(sawTotal.alpine$comm))
ericalestotal.sesmpd.phylonull.talus <- ses.mpd(t(saw.total.ericales.talus$data), cophenetic.phylo(saw.total.ericales.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
ericalestotal.sesmntd.phylonull.talus <- ses.mntd(t(saw.total.ericales.talus$data), cophenetic.phylo(saw.total.ericales.talus$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)


##### Pool = only species found on meadow

# Tracheophyta
saw.total.sesmpd.phylonull.meadow <- ses.mpd(sawTotal.meadow$comm, cophenetic.phylo(sawTotal.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
saw.total.sesmntd.phylonull.meadow <- ses.mntd(sawTotal.meadow$comm, cophenetic.phylo(sawTotal.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.asterales.meadow <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.meadow$phy$tip.label, tax, level = "order", taxonomy = "Asterales")
pruned.tree.saw.total.aster.meadow <-drop.tip(sawTotal.meadow$phy, sawTotal.meadow$phy$tip.label[!sawTotal.meadow$phy$tip.label %in% sawtooth.asterales])
saw.total.aster.meadow <- treedata(pruned.tree.saw.total.aster.meadow, t(sawTotal.meadow$comm))
astertotal.sesmpd.phylonull.meadow <- ses.mpd(t(saw.total.aster.meadow$data), cophenetic.phylo(saw.total.aster.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
astertotal.sesmntd.phylonull.meadow <- ses.mntd(t(saw.total.aster.meadow$data), cophenetic.phylo(saw.total.aster.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.caryophyllales.meadow <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.meadow$phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
pruned.tree.saw.total.caryo.meadow <-drop.tip(sawTotal.meadow$phy, sawTotal.meadow$phy$tip.label[!sawTotal.meadow$phy$tip.label %in% sawtooth.caryophyllales])
saw.total.caryo.meadow <- treedata(pruned.tree.saw.total.caryo.meadow, t(sawTotal.meadow$comm))
caryototal.sesmpd.phylonull.meadow <- ses.mpd(t(saw.total.caryo.meadow$data), cophenetic.phylo(saw.total.caryo.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
caryototal.sesmntd.phylonull.meadow <- ses.mntd(t(saw.total.caryo.meadow$data), cophenetic.phylo(saw.total.caryo.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.poales.meadow <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.meadow$phy$tip.label, tax, level = "order", taxonomy = "Poales")
pruned.tree.saw.total.poales.meadow <-drop.tip(sawTotal.meadow$phy, sawTotal.meadow$phy$tip.label[!sawTotal.meadow$phy$tip.label %in% sawtooth.poales])
saw.total.poales.meadow <- treedata(pruned.tree.saw.total.poales.meadow, t(sawTotal.meadow$comm))
poalestotal.sesmpd.phylonull.meadow <- ses.mpd(t(saw.total.poales.meadow$data), cophenetic.phylo(saw.total.poales.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
poalestotal.sesmntd.phylonull.meadow <- ses.mntd(t(saw.total.poales.meadow$data), cophenetic.phylo(saw.total.poales.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.lamiales.meadow <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.meadow$phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
pruned.tree.saw.total.lamiales.meadow <-drop.tip(sawTotal.meadow$phy, sawTotal.meadow$phy$tip.label[!sawTotal.meadow$phy$tip.label %in% sawtooth.lamiales])
saw.total.lamiales.meadow <- treedata(pruned.tree.saw.total.lamiales.meadow, t(sawTotal.meadow$comm))
lamialestotal.sesmpd.phylonull.meadow <- ses.mpd(t(saw.total.lamiales.meadow$data), cophenetic.phylo(saw.total.lamiales.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
lamialestotal.sesmntd.phylonull.meadow <- ses.mntd(t(saw.total.lamiales.meadow$data), cophenetic.phylo(saw.total.lamiales.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.brassicales.meadow <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.meadow$phy$tip.label, tax, level = "order", taxonomy = "Brassicales")
pruned.tree.saw.total.brassicales.meadow <-drop.tip(sawTotal.meadow$phy, sawTotal.meadow$phy$tip.label[!sawTotal.meadow$phy$tip.label %in% sawtooth.brassicales])
saw.total.brassicales.meadow <- treedata(pruned.tree.saw.total.brassicales.meadow, t(sawTotal.meadow$comm))
brassicalestotal.sesmpd.phylonull.meadow <- ses.mpd(t(saw.total.brassicales.meadow$data), cophenetic.phylo(saw.total.brassicales.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
brassicalestotal.sesmntd.phylonull.meadow <- ses.mntd(t(saw.total.brassicales.meadow$data), cophenetic.phylo(saw.total.brassicales.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.ericales.meadow <- pruneCladeTaxonomyLookup(tip.labels = sawTotal.meadow$phy$tip.label, tax, level = "order", taxonomy = "Ericales")
pruned.tree.saw.total.ericales.meadow <-drop.tip(sawTotal.meadow$phy, sawTotal.meadow$phy$tip.label[!sawTotal.meadow$phy$tip.label %in% sawtooth.ericales])
saw.total.ericales.meadow <- treedata(pruned.tree.saw.total.ericales.meadow, t(sawTotal.meadow$comm))
ericalestotal.sesmpd.phylonull.meadow <- ses.mpd(t(saw.total.ericales.meadow$data), cophenetic.phylo(saw.total.ericales.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
ericalestotal.sesmntd.phylonull.meadow <- ses.mntd(t(saw.total.ericales.meadow$data), cophenetic.phylo(saw.total.ericales.meadow$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)



#################### Combine all output of the Total Dataset

phylogeny.poolSESmntd.total <- rbind((saw.total.sesmntd.phylonull %>%  mutate(summits = rownames(saw.total.sesmntd.phylonull), metric = "mntd", clade = "Tracheophyta", sig=ifelse(saw.total.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | saw.total.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                     (astertotal.sesmntd.phylonull %>%  mutate(summits = rownames(astertotal.sesmntd.phylonull), metric = "mntd", clade = "Asterales", sig=ifelse(astertotal.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | astertotal.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))), 
                                     (caryototal.sesmntd.phylonull %>%  mutate(summits = rownames(caryototal.sesmntd.phylonull), metric = "mntd", clade = "Caryophyllales", sig=ifelse(caryototal.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | caryototal.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                     (poalestotal.sesmntd.phylonull %>%  mutate(summits = rownames(poalestotal.sesmntd.phylonull), metric = "mntd", clade = "Poales", sig=ifelse(poalestotal.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | poalestotal.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                     (lamialestotal.sesmntd.phylonull %>%  mutate(summits = rownames(lamialestotal.sesmntd.phylonull), metric = "mntd", clade = "Lamiales", sig=ifelse(lamialestotal.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | lamialestotal.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                     (brassicalestotal.sesmntd.phylonull %>%  mutate(summits = rownames(brassicalestotal.sesmntd.phylonull), metric = "mntd", clade = "Brassicales", sig=ifelse(brassicalestotal.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | brassicalestotal.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                     (ericalestotal.sesmntd.phylonull %>%  mutate(summits = rownames(ericalestotal.sesmntd.phylonull), metric = "mntd", clade = "Ericales", sig=ifelse(ericalestotal.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | ericalestotal.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))))
names(phylogeny.poolSESmntd.total)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

phylogeny.poolSESmpd.total <- rbind((saw.total.sesmpd.phylonull %>%  mutate(summits = rownames(saw.total.sesmpd.phylonull), metric = "mpd", clade = "Tracheophyta", sig=ifelse(saw.total.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | saw.total.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                    (astertotal.sesmpd.phylonull %>%  mutate(summits = rownames(astertotal.sesmpd.phylonull), metric = "mpd", clade = "Asterales", sig=ifelse(astertotal.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | astertotal.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))), 
                                    (caryototal.sesmpd.phylonull %>%  mutate(summits = rownames(caryototal.sesmpd.phylonull), metric = "mpd", clade = "Caryophyllales", sig=ifelse(caryototal.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | caryototal.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                    (poalestotal.sesmpd.phylonull %>%  mutate(summits = rownames(poalestotal.sesmpd.phylonull), metric = "mpd", clade = "Poales", sig=ifelse(poalestotal.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | poalestotal.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))), 
                                    (lamialestotal.sesmpd.phylonull %>%  mutate(summits = rownames(lamialestotal.sesmpd.phylonull), metric = "mpd", clade = "Lamiales", sig=ifelse(lamialestotal.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | lamialestotal.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                    (brassicalestotal.sesmpd.phylonull %>%  mutate(summits = rownames(brassicalestotal.sesmpd.phylonull), metric = "mpd", clade = "Brassicales", sig=ifelse(brassicalestotal.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | brassicalestotal.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                    (ericalestotal.sesmpd.phylonull %>%  mutate(summits = rownames(ericalestotal.sesmpd.phylonull), metric = "mpd", clade = "Ericales", sig=ifelse(ericalestotal.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | ericalestotal.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))))
names(phylogeny.poolSESmpd.total)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")
phylogeny.pool.all.SES.total <- rbind(phylogeny.poolSESmntd.total, phylogeny.poolSESmpd.total)  
phylogeny.pool.all.SES.total <- cbind(phylogeny.pool.all.SES.total, pool = rep(x = "All Alpine", times = nrow(phylogeny.pool.all.SES.total)))
phylogeny.pool.all.SES.total <- cbind(phylogeny.pool.all.SES.total, data = rep(x = "total", times = nrow(phylogeny.pool.all.SES.total)))

# Talus
phylogeny.pool.talus.SESmntd.total <- rbind((saw.total.sesmntd.phylonull.talus %>%  mutate(summits = rownames(saw.total.sesmntd.phylonull.talus), metric = "mntd", clade = "Tracheophyta", sig=ifelse(saw.total.sesmntd.phylonull.talus[,"mntd.obs.p"] <= 0.05 | saw.total.sesmntd.phylonull.talus[,"mntd.obs.p"] > 0.95, 1,0))),
                                            (astertotal.sesmntd.phylonull.talus %>%  mutate(summits = rownames(astertotal.sesmntd.phylonull.talus), metric = "mntd", clade = "Asterales", sig=ifelse(astertotal.sesmntd.phylonull.talus[,"mntd.obs.p"] <= 0.05 | astertotal.sesmntd.phylonull.talus[,"mntd.obs.p"] > 0.95, 1,0))), 
                                            (caryototal.sesmntd.phylonull.talus %>%  mutate(summits = rownames(caryototal.sesmntd.phylonull.talus), metric = "mntd", clade = "Caryophyllales", sig=ifelse(caryototal.sesmntd.phylonull.talus[,"mntd.obs.p"] <= 0.05 | caryototal.sesmntd.phylonull.talus[,"mntd.obs.p"] > 0.95, 1,0))),
                                            (poalestotal.sesmntd.phylonull.talus %>%  mutate(summits = rownames(poalestotal.sesmntd.phylonull.talus), metric = "mntd", clade = "Poales", sig=ifelse(poalestotal.sesmntd.phylonull.talus[,"mntd.obs.p"] <= 0.05 | poalestotal.sesmntd.phylonull.talus[,"mntd.obs.p"] > 0.95, 1,0))),
                                            (lamialestotal.sesmntd.phylonull.talus %>%  mutate(summits = rownames(lamialestotal.sesmntd.phylonull.talus), metric = "mntd", clade = "Lamiales", sig=ifelse(lamialestotal.sesmntd.phylonull.talus[,"mntd.obs.p"] <= 0.05 | lamialestotal.sesmntd.phylonull.talus[,"mntd.obs.p"] > 0.95, 1,0))),
                                            (brassicalestotal.sesmntd.phylonull.talus %>%  mutate(summits = rownames(brassicalestotal.sesmntd.phylonull.talus), metric = "mntd", clade = "Brassicales", sig=ifelse(brassicalestotal.sesmntd.phylonull.talus[,"mntd.obs.p"] <= 0.05 | brassicalestotal.sesmntd.phylonull.talus[,"mntd.obs.p"] > 0.95, 1,0))),
                                            (ericalestotal.sesmntd.phylonull.talus %>%  mutate(summits = rownames(ericalestotal.sesmntd.phylonull.talus), metric = "mntd", clade = "Ericales", sig=ifelse(ericalestotal.sesmntd.phylonull.talus[,"mntd.obs.p"] <= 0.05 | ericalestotal.sesmntd.phylonull.talus[,"mntd.obs.p"] > 0.95, 1,0))))
names(phylogeny.pool.talus.SESmntd.total)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

phylogeny.pool.talus.SESmpd.total <- rbind((saw.total.sesmpd.phylonull.talus %>%  mutate(summits = rownames(saw.total.sesmpd.phylonull.talus), metric = "mpd", clade = "Tracheophyta", sig=ifelse(saw.total.sesmpd.phylonull.talus[,"mpd.obs.p"] <= 0.05 | saw.total.sesmpd.phylonull.talus[,"mpd.obs.p"] > 0.95, 1,0))),
                                           (astertotal.sesmpd.phylonull.talus %>%  mutate(summits = rownames(astertotal.sesmpd.phylonull.talus), metric = "mpd", clade = "Asterales", sig=ifelse(astertotal.sesmpd.phylonull.talus[,"mpd.obs.p"] <= 0.05 | astertotal.sesmpd.phylonull.talus[,"mpd.obs.p"] > 0.95, 1,0))), 
                                           (caryototal.sesmpd.phylonull.talus %>%  mutate(summits = rownames(caryototal.sesmpd.phylonull.talus), metric = "mpd", clade = "Caryophyllales", sig=ifelse(caryototal.sesmpd.phylonull.talus[,"mpd.obs.p"] <= 0.05 | caryototal.sesmpd.phylonull.talus[,"mpd.obs.p"] > 0.95, 1,0))),
                                           (poalestotal.sesmpd.phylonull.talus %>%  mutate(summits = rownames(poalestotal.sesmpd.phylonull.talus), metric = "mpd", clade = "Poales", sig=ifelse(poalestotal.sesmpd.phylonull.talus[,"mpd.obs.p"] <= 0.05 | poalestotal.sesmpd.phylonull.talus[,"mpd.obs.p"] > 0.95, 1,0))),
                                           (lamialestotal.sesmpd.phylonull.talus %>%  mutate(summits = rownames(lamialestotal.sesmpd.phylonull.talus), metric = "mpd", clade = "Lamiales", sig=ifelse(lamialestotal.sesmpd.phylonull.talus[,"mpd.obs.p"] <= 0.05 | lamialestotal.sesmpd.phylonull.talus[,"mpd.obs.p"] > 0.95, 1,0))),
                                           (brassicalestotal.sesmpd.phylonull.talus %>%  mutate(summits = rownames(brassicalestotal.sesmpd.phylonull.talus), metric = "mpd", clade = "Brassicales", sig=ifelse(brassicalestotal.sesmpd.phylonull.talus[,"mpd.obs.p"] <= 0.05 | brassicalestotal.sesmpd.phylonull.talus[,"mpd.obs.p"] > 0.95, 1,0))),
                                           (ericalestotal.sesmpd.phylonull.talus %>%  mutate(summits = rownames(ericalestotal.sesmpd.phylonull.talus), metric = "mpd", clade = "Ericales", sig=ifelse(ericalestotal.sesmpd.phylonull.talus[,"mpd.obs.p"] <= 0.05 | ericalestotal.sesmpd.phylonull.talus[,"mpd.obs.p"] > 0.95, 1,0))))
names(phylogeny.pool.talus.SESmpd.total)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")
phylogeny.pool.talus.SES.total <- rbind(phylogeny.pool.talus.SESmntd.total, phylogeny.pool.talus.SESmpd.total)                               
phylogeny.pool.talus.SES.total <- cbind(phylogeny.pool.talus.SES.total, pool = rep(x = "Talus", times = nrow(phylogeny.pool.talus.SES.total)))
phylogeny.pool.talus.SES.total <- cbind(phylogeny.pool.talus.SES.total, data = rep(x = "total", times = nrow(phylogeny.pool.talus.SES.total)))

##### Meadow
phylogeny.pool.meadow.SESmntd.total <- rbind((saw.total.sesmntd.phylonull.meadow %>%  mutate(summits = rownames(saw.total.sesmntd.phylonull.meadow), metric = "mntd", clade = "Tracheophyta", sig=ifelse(saw.total.sesmntd.phylonull.meadow[,"mntd.obs.p"] <= 0.05 | saw.total.sesmntd.phylonull.meadow[,"mntd.obs.p"] > 0.95, 1,0))),
                                             (astertotal.sesmntd.phylonull.meadow %>%  mutate(summits = rownames(astertotal.sesmntd.phylonull.meadow), metric = "mntd", clade = "Asterales", sig=ifelse(astertotal.sesmntd.phylonull.meadow[,"mntd.obs.p"] <= 0.05 | astertotal.sesmntd.phylonull.meadow[,"mntd.obs.p"] > 0.95, 1,0))), 
                                             (caryototal.sesmntd.phylonull.meadow %>%  mutate(summits = rownames(caryototal.sesmntd.phylonull.meadow), metric = "mntd", clade = "Caryophyllales", sig=ifelse(caryototal.sesmntd.phylonull.meadow[,"mntd.obs.p"] <= 0.05 | caryototal.sesmntd.phylonull.meadow[,"mntd.obs.p"] > 0.95, 1,0))),
                                             (poalestotal.sesmntd.phylonull.meadow %>%  mutate(summits = rownames(poalestotal.sesmntd.phylonull.meadow), metric = "mntd", clade = "Poales", sig=ifelse(poalestotal.sesmntd.phylonull.meadow[,"mntd.obs.p"] <= 0.05 | poalestotal.sesmntd.phylonull.meadow[,"mntd.obs.p"] > 0.95, 1,0))),
                                             (lamialestotal.sesmntd.phylonull.meadow %>%  mutate(summits = rownames(lamialestotal.sesmntd.phylonull.meadow), metric = "mntd", clade = "Lamiales", sig=ifelse(lamialestotal.sesmntd.phylonull.meadow[,"mntd.obs.p"] <= 0.05 | lamialestotal.sesmntd.phylonull.meadow[,"mntd.obs.p"] > 0.95, 1,0))),
                                             (brassicalestotal.sesmntd.phylonull.meadow %>%  mutate(summits = rownames(brassicalestotal.sesmntd.phylonull.meadow), metric = "mntd", clade = "Brassicales", sig=ifelse(brassicalestotal.sesmntd.phylonull.meadow[,"mntd.obs.p"] <= 0.05 | brassicalestotal.sesmntd.phylonull.meadow[,"mntd.obs.p"] > 0.95, 1,0))),
                                             (ericalestotal.sesmntd.phylonull.meadow %>%  mutate(summits = rownames(ericalestotal.sesmntd.phylonull.meadow), metric = "mntd", clade = "Ericales", sig=ifelse(ericalestotal.sesmntd.phylonull.meadow[,"mntd.obs.p"] <= 0.05 | ericalestotal.sesmntd.phylonull.meadow[,"mntd.obs.p"] > 0.95, 1,0))))
names(phylogeny.pool.meadow.SESmntd.total)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

phylogeny.pool.meadow.SESmpd.total <- rbind((saw.total.sesmpd.phylonull.meadow %>%  mutate(summits = rownames(saw.total.sesmpd.phylonull.meadow), metric = "mpd", clade = "Tracheophyta", sig=ifelse(saw.total.sesmpd.phylonull.meadow[,"mpd.obs.p"] <= 0.05 | saw.total.sesmpd.phylonull.meadow[,"mpd.obs.p"] > 0.95, 1,0))),
                                            (astertotal.sesmpd.phylonull.meadow %>%  mutate(summits = rownames(astertotal.sesmpd.phylonull.meadow), metric = "mpd", clade = "Asterales", sig=ifelse(astertotal.sesmpd.phylonull.meadow[,"mpd.obs.p"] <= 0.05 | astertotal.sesmpd.phylonull.meadow[,"mpd.obs.p"] > 0.95, 1,0))), 
                                            (caryototal.sesmpd.phylonull.meadow %>%  mutate(summits = rownames(caryototal.sesmpd.phylonull.meadow), metric = "mpd", clade = "Caryophyllales", sig=ifelse(caryototal.sesmpd.phylonull.meadow[,"mpd.obs.p"] <= 0.05 | caryototal.sesmpd.phylonull.meadow[,"mpd.obs.p"] > 0.95, 1,0))),
                                            (poalestotal.sesmpd.phylonull.meadow %>%  mutate(summits = rownames(poalestotal.sesmpd.phylonull.meadow), metric = "mpd", clade = "Poales", sig=ifelse(poalestotal.sesmpd.phylonull.meadow[,"mpd.obs.p"] <= 0.05 | poalestotal.sesmpd.phylonull.meadow[,"mpd.obs.p"] > 0.95, 1,0))),
                                            (lamialestotal.sesmpd.phylonull.meadow %>%  mutate(summits = rownames(lamialestotal.sesmpd.phylonull.meadow), metric = "mpd", clade = "Lamiales", sig=ifelse(lamialestotal.sesmpd.phylonull.meadow[,"mpd.obs.p"] <= 0.05 | lamialestotal.sesmpd.phylonull.meadow[,"mpd.obs.p"] > 0.95, 1,0))),
                                            (brassicalestotal.sesmpd.phylonull.meadow %>%  mutate(summits = rownames(brassicalestotal.sesmpd.phylonull.meadow), metric = "mpd", clade = "Brassicales", sig=ifelse(brassicalestotal.sesmpd.phylonull.meadow[,"mpd.obs.p"] <= 0.05 | brassicalestotal.sesmpd.phylonull.meadow[,"mpd.obs.p"] > 0.95, 1,0))),
                                            (ericalestotal.sesmpd.phylonull.meadow %>%  mutate(summits = rownames(ericalestotal.sesmpd.phylonull.meadow), metric = "mpd", clade = "Ericales", sig=ifelse(ericalestotal.sesmpd.phylonull.meadow[,"mpd.obs.p"] <= 0.05 | ericalestotal.sesmpd.phylonull.meadow[,"mpd.obs.p"] > 0.95, 1,0))))
names(phylogeny.pool.meadow.SESmpd.total)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")
phylogeny.pool.meadow.SES.total <- rbind(phylogeny.pool.meadow.SESmntd.total, phylogeny.pool.meadow.SESmpd.total)                               
phylogeny.pool.meadow.SES.total <- cbind(phylogeny.pool.meadow.SES.total, pool = rep(x = "Meadow", times = nrow(phylogeny.pool.meadow.SES.total)))
phylogeny.pool.meadow.SES.total <- cbind(phylogeny.pool.meadow.SES.total, data = rep(x = "total", times = nrow(phylogeny.pool.meadow.SES.total)))


###################### Master alpha phylogenetic diversity 
alpha.ses <- rbind(phylogeny.pool.all.SES.total, phylogeny.pool.talus.SES.total, phylogeny.pool.meadow.SES.total)
head(alpha.ses)     
alpha.ses$clade <- factor(alpha.ses$clade, levels = c("Ericales", "Brassicales", "Lamiales", "Caryophyllales", "Poales","Asterales", "Tracheophyta"))
alpha.ses <- cbind(alpha.ses, model = rep(x = "RD", times = nrow(alpha.ses)))

write.csv(alpha.ses, file="output/10_PhyloDiversity/alpha/sawtooth.alpha.RD.SES.csv")

