
##############################
## Dynamic null model of communtiy assembly in the Sawtooth NF, Idaho, USA

args <- commandArgs(TRUE)

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)

source("chooseClade.cluster.R")

sawtooths.asterales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqBeta$phy$tip.label, tax, level = "order", taxonomy = "Asterales")
sawtooths.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqBeta$phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
sawtooths.brassicales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqBeta$phy$tip.label, tax, level = "order", taxonomy = "Brassicales")
sawtooths.ericales <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqBeta$phy$tip.label, tax, level = "order", taxonomy = "Ericales")


## Repeat optimization of parameters 10 times to see if they converge 
# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs

initparsopt <- sample(c(1:100), size=2, replace = T)

## Spermatophyta
alps.damocles <- treedata(sawMiseqBeta$phy, t(sawMiseqBeta$comm))

opt.params.spermatophyta <- DAMOCLES_ML(phy = alps.damocles$phy,
                     pa = as.matrix(alps.damocles$data)[,9],
                     initparsopt = initparsopt,
                     idparsopt = 1:2,
                     parsfix = 0,
                     idparsfix = 3,
                     pars2 = c(1E-3,1E-4,1E-5,1000),
                     pchoice = 0)

write.csv(opt.params.spermatophyta, file=paste("opt.params.spermatophyta.", args[1], ".csv", sep=""))

############ just on asterales
pruned.tree.alps.aster <-drop.tip(sawMiseqBeta$phy, sawMiseqBeta$phy$tip.label[!sawMiseqBeta$phy$tip.label %in% sawtooths.asterales])
alps.damocles.aster <- treedata(pruned.tree.alps.aster, t(sawMiseqBeta$comm))

opt.params.asterales <- DAMOCLES_ML(phy = alps.damocles.aster$phy,
                                        pa = alps.damocles.aster$data[,9],
                                        initparsopt = initparsopt,
                                        idparsopt = 1:2,
                                        parsfix = 0,
                                        idparsfix = 3,
                                        pars2 = c(1E-3,1E-4,1E-5,1000),
                                        pchoice = 0)

write.csv(opt.params.asterales, file=paste("opt.params.asterales.", args[1], ".csv", sep=""))


############ just on caryophyllales
pruned.tree.alps.caryophyllales <-drop.tip(sawMiseqBeta$phy, sawMiseqBeta$phy$tip.label[!sawMiseqBeta$phy$tip.label %in% sawtooths.caryophyllales])
alps.damocles.caryophyllales <- treedata(pruned.tree.alps.caryophyllales, t(sawMiseqBeta$comm))

opt.params.caryophyllales <- DAMOCLES_ML(phy = alps.damocles.caryophyllales$phy,
                                         pa = alps.damocles.caryophyllales$data[,9],
                                         initparsopt = initparsopt,
                                         idparsopt = 1:2,
                                         parsfix = 0,
                                         idparsfix = 3,
                                         pars2 = c(1E-3,1E-4,1E-5,1000),
                                         pchoice = 0)

write.csv(opt.params.caryophyllales, file=paste("opt.params.caryophyllales.", args[1], ".csv", sep=""))



############ just on brassicales
pruned.tree.alps.brassicales <-drop.tip(sawMiseqBeta$phy, sawMiseqBeta$phy$tip.label[!sawMiseqBeta$phy$tip.label %in% sawtooths.brassicales])
alps.damocles.brassicales <- treedata(pruned.tree.alps.brassicales, t(sawMiseqBeta$comm))

opt.params.brassicales <- DAMOCLES_ML(phy = alps.damocles.brassicales$phy,
                                    pa = alps.damocles.brassicales$data[,9],
                                    initparsopt = initparsopt,
                                    idparsopt = 1:2,
                                    parsfix = 0,
                                    idparsfix = 3,
                                    pars2 = c(1E-3,1E-4,1E-5,1000),
                                    pchoice = 0)

write.csv(opt.params.brassicales, file=paste("opt.params.brassicales.", args[1], ".csv", sep=""))


############ just on ericales
pruned.tree.alps.ericales <-drop.tip(sawMiseqBeta$phy, sawMiseqBeta$phy$tip.label[!sawMiseqBeta$phy$tip.label %in% sawtooths.ericales])
alps.damocles.ericales <- treedata(pruned.tree.alps.ericales, t(sawMiseqBeta$comm))

opt.params.ericales <- DAMOCLES_ML(phy = alps.damocles.ericales$phy,
                                    pa = alps.damocles.ericales$data[,9],
                                    initparsopt = initparsopt,
                                    idparsopt = 1:2,
                                    parsfix = 0,
                                    idparsfix = 3,
                                    pars2 = c(1E-3,1E-4,1E-5,1000),
                                    pchoice = 0)

write.csv(opt.params.ericales, file=paste("opt.params.ericales.", args[1], ".csv", sep=""))








