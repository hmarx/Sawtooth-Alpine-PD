########################################### Alpha Phylo Diversity ########################################### 

source("R/chooseClade.R")

## Random resample from phylogeny pool (==Sawtooth NF), equal probability random draw from phylogeny pool
saw.sesmpd.phylonull.talus.talus <- ses.mpd(sawMiseqAlp$comm, cophenetic.phylo(sawMiseqAlp$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
write.csv(saw.sesmpd.phylonull.talus, file="output/08_PhyloDiversity/MiSeq/alpha/static/saw.sesmpd.phylonull.talus.talus.csv")

saw.sesmntd.phylonull.talus.talus <- ses.mntd(sawMiseqAlp$comm, cophenetic.phylo(sawMiseqAlp$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
write.csv(saw.sesmntd.phylonull.talus, file="output/08_PhyloDiversity/MiSeq/alpha/static/saw.sesmntd.phylonull.talus.talus.csv")


# Pruned phylogeny and community matrix for each clade 
sawtooth.asterales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqAlp$phy$tip.label, tax, level = "order", taxonomy = "Asterales")
pruned.tree.saw.aster <-drop.tip(sawMiseqAlp$phy, sawMiseqAlp$phy$tip.label[!sawMiseqAlp$phy$tip.label %in% sawtooth.asterales])
saw.aster <- treedata(pruned.tree.saw.aster, sawtooth.com.miseq.alpine)
aster.sesmpd.phylonull <- ses.mpd(t(saw.aster$data), cophenetic.phylo(saw.aster$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
aster.sesmntd.phylonull <- ses.mntd(t(saw.aster$data), cophenetic.phylo(saw.aster$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqAlp$phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
pruned.tree.saw.caryo <-drop.tip(sawMiseqAlp$phy, sawMiseqAlp$phy$tip.label[!sawMiseqAlp$phy$tip.label %in% sawtooth.caryophyllales])
saw.caryo <- treedata(pruned.tree.saw.caryo, sawtooth.com.miseq.alpine)
caryo.sesmpd.phylonull <- ses.mpd(t(saw.caryo$data), cophenetic.phylo(saw.caryo$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
caryo.sesmntd.phylonull <- ses.mntd(t(saw.caryo$data), cophenetic.phylo(saw.caryo$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.brassicales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqAlp$phy$tip.label, tax, level = "order", taxonomy = "Brassicales")
pruned.tree.saw.brassicales <-drop.tip(sawMiseqAlp$phy, sawMiseqAlp$phy$tip.label[!sawMiseqAlp$phy$tip.label %in% sawtooth.brassicales])
saw.brassicales <- treedata(pruned.tree.saw.brassicales, sawtooth.com.miseq.alpine)
brassicales.sesmpd.phylonull <- ses.mpd(t(saw.brassicales$data), cophenetic.phylo(saw.brassicales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
brassicales.sesmntd.phylonull <- ses.mntd(t(saw.brassicales$data), cophenetic.phylo(saw.brassicales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.poales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqAlp$phy$tip.label, tax, level = "order", taxonomy = "Poales")
pruned.tree.saw.poales <-drop.tip(sawMiseqAlp$phy, sawMiseqAlp$phy$tip.label[!sawMiseqAlp$phy$tip.label %in% sawtooth.poales])
saw.poales <- treedata(pruned.tree.saw.poales, sawtooth.com.miseq.alpine)
poales.sesmpd.phylonull <- ses.mpd(t(saw.poales$data), cophenetic.phylo(saw.poales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
poales.sesmntd.phylonull <- ses.mntd(t(saw.poales$data), cophenetic.phylo(saw.poales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.saxifragales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqAlp$phy$tip.label, tax, level = "order", taxonomy = "Saxifragales")
pruned.tree.saw.saxifragales <-drop.tip(sawMiseqAlp$phy, sawMiseqAlp$phy$tip.label[!sawMiseqAlp$phy$tip.label %in% sawtooth.saxifragales])
saw.saxifragales <- treedata(pruned.tree.saw.saxifragales, sawtooth.com.miseq.alpine)
saxifragales.sesmpd.phylonull <- ses.mpd(t(saw.saxifragales$data), cophenetic.phylo(saw.saxifragales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
saxifragales.sesmntd.phylonull <- ses.mntd(t(saw.saxifragales$data), cophenetic.phylo(saw.saxifragales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)

sawtooth.lamiales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqAlp$phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
pruned.tree.saw.lamiales <-drop.tip(sawMiseqAlp$phy, sawMiseqAlp$phy$tip.label[!sawMiseqAlp$phy$tip.label %in% sawtooth.lamiales])
saw.lamiales <- treedata(pruned.tree.saw.lamiales, sawtooth.com.miseq.alpine)
lamiales.sesmpd.phylonull <- ses.mpd(t(saw.lamiales$data), cophenetic.phylo(saw.lamiales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
lamiales.sesmntd.phylonull <- ses.mntd(t(saw.lamiales$data), cophenetic.phylo(saw.lamiales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)


sawtooth.ericales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqAlp$phy$tip.label, tax, level = "order", taxonomy = "Ericales")
pruned.tree.saw.ericales <-drop.tip(sawMiseqAlp$phy, sawMiseqAlp$phy$tip.label[!sawMiseqAlp$phy$tip.label %in% sawtooth.ericales])
saw.ericales <- treedata(pruned.tree.saw.ericales, sawtooth.com.miseq.alpine)
ericales.sesmpd.phylonull <- ses.mpd(t(saw.ericales$data), cophenetic.phylo(saw.ericales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)
ericales.sesmntd.phylonull <- ses.mntd(t(saw.ericales$data), cophenetic.phylo(saw.ericales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 10000)


##### Combine all output of randomization from phylogeny pool
phylogeny.poolSESmntd <- rbind((saw.sesmntd.phylonull.talus.talus %>%  mutate(summits = rownames(saw.sesmntd.phylonull.talus.talus), metric = "mntd", clade = "Spermatophyta", sig=ifelse(saw.sesmntd.phylonull.talus.talus[,"mntd.obs.p"] <= 0.05 | saw.sesmntd.phylonull.talus.talus[,"mntd.obs.p"] > 0.95, 1,0))),
                               (aster.sesmntd.phylonull %>%  mutate(summits = rownames(aster.sesmntd.phylonull), metric = "mntd", clade = "Asterales", sig=ifelse(aster.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | aster.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))), 
                               (caryo.sesmntd.phylonull %>%  mutate(summits = rownames(caryo.sesmntd.phylonull), metric = "mntd", clade = "Caryophyllales", sig=ifelse(caryo.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | caryo.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                               (brassicales.sesmntd.phylonull %>%  mutate(summits = rownames(brassicales.sesmntd.phylonull), metric = "mntd", clade = "Brassicales", sig=ifelse(brassicales.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | brassicales.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                               (poales.sesmntd.phylonull %>%  mutate(summits = rownames(poales.sesmntd.phylonull), metric = "mntd", clade = "Poales", sig=ifelse(poales.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | poales.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                               (saxifragales.sesmntd.phylonull %>%  mutate(summits = rownames(saxifragales.sesmntd.phylonull), metric = "mntd", clade = "Saxifragales", sig=ifelse(saxifragales.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | saxifragales.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))))
names(phylogeny.poolSESmntd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

phylogeny.poolSESmpd <- rbind((saw.sesmpd.phylonull.talus.talus %>%  mutate(summits = rownames(saw.sesmpd.phylonull.talus.talus), metric = "mpd", clade = "Spermatophyta", sig=ifelse(saw.sesmpd.phylonull.talus.talus[,"mpd.obs.p"] <= 0.05 | saw.sesmpd.phylonull.talus.talus[,"mpd.obs.p"] > 0.95, 1,0))),
                               (aster.sesmpd.phylonull %>%  mutate(summits = rownames(aster.sesmpd.phylonull), metric = "mpd", clade = "Asterales", sig=ifelse(aster.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | aster.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))), 
                               (caryo.sesmpd.phylonull %>%  mutate(summits = rownames(caryo.sesmpd.phylonull), metric = "mpd", clade = "Caryophyllales", sig=ifelse(caryo.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | caryo.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                               (brassicales.sesmpd.phylonull %>%  mutate(summits = rownames(brassicales.sesmpd.phylonull), metric = "mpd", clade = "Brassicales", sig=ifelse(brassicales.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | brassicales.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                               (poales.sesmpd.phylonull %>%  mutate(summits = rownames(poales.sesmpd.phylonull), metric = "mpd", clade = "Poales", sig=ifelse(poales.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | poales.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                               (saxifragales.sesmpd.phylonull %>%  mutate(summits = rownames(saxifragales.sesmpd.phylonull), metric = "mpd", clade = "Saxifragales", sig=ifelse(saxifragales.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | saxifragales.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))))
names(phylogeny.poolSESmpd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")
phylogeny.poolSES <- rbind(phylogeny.poolSESmntd, phylogeny.poolSESmpd)                               
head(phylogeny.poolSES)     
phylogeny.poolSES$clade <- factor(phylogeny.poolSES$clade, levels = c("Saxifragales", "Poales", "Brassicales", "Caryophyllales", "Asterales", "Spermatophyta"))
phylogeny.poolSES <- cbind(phylogeny.poolSES, pool = rep(x = "Sawtooth Alpine", times = nrow(phylogeny.poolSES)))
write.csv(phylogeny.poolSES, file="output/08_PhyloDiversity/MiSeq/alpha/static/alpine.phylogeny.pool.SES.talus.csv")

