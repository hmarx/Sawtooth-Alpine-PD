
## Loop DAMOCLES over all summits

##############################
## Dynamic null model of communtiy assembly in the Ecrin NP, France
#args <- commandArgs(trailingOnly = TRUE)
#i = as.numeric(args[1])

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)

setwd("output/08_PhyloDiversity/MiSeq//dynamic//eachSummit")

source("chooseClade.cluster.R")

sawtooth.spermatophyta <- pruneCladeTaxonomyLookup(tip.labels = sawMiseqBeta$phy$tip.label, tax, level = "Spermatophyta", taxonomy = "Spermatophyta")

#set parameter values
idparsopt = 1:2 #optimize extinction and colonization 

# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs
initparsopt = c(0.1,0.1)

############ just on spermatophyta in persistent source pool
pruned.tree.saw.sperma <-drop.tip(sawMiseqBeta$phy, sawMiseqBeta$tip.label[!sawMiseqBeta$tip.label %in% sawtooth.spermatophyta])
saw.damocles.sperma <- treedata(pruned.tree.saw.sperma, cbind(rownames(t(sawMiseqBeta$comm)), t(sawMiseqBeta$comm)))


for (i in 2:ncol(saw.damocles.sperma$data)){
  damo.bs.sperma <- DAMOCLES_bootstrap(phy = saw.damocles.sperma$phy, pa = as.matrix(saw.damocles.sperma$data[,c(1,i)]), idparsopt = 1:length(initparsopt), parsfix = 0, idparsfix = (1:3)[-idparsopt], pars2 = c(1E-3,1E-4,1E-5,1000), pchoice = 0, runs = 1000, estimate_pars = FALSE, conf.int = 0.95)
  
  write.csv(damo.bs.sperma$null_community_data, file=paste("spermatophyta.sawtooth.null_community_data.", i, ".csv", sep=""))
  write.csv(damo.bs.sperma$summary_table, file=paste("spermatophyta.sawtooth.summary_table.", i, ".csv", sep=""))
  
  
} #2:10



