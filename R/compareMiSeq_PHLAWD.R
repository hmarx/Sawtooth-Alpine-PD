######### Assess monophyletic clades (accuracy of MiSeq vs. PHLAWD)

## Get a taxonomic lookup table: https://github.com/wcornwell/TaxonLookup
head(plant_lookup())
head(plant_lookup(include_counts = TRUE))
tx <- lookup_table(genetree$tip.label, by_species=TRUE)

#tx <- read.csv("~/Desktop/TaxonLookup/plant_lookup.csv")
head(tx)
dim(tx)

## Assess monophyly: package MonoPhy
dat <- as.data.frame(treedata(genetree, tx)$dat)
dat <- cbind(rownames(dat), dat[1:3])
head(dat)
tr <- (treedata(genetree, tx)$phy)

mono <- AssessMonophyly(tree = tr, taxonomy = dat)
head(mono[1:5])
mono$taxonomy1$result #genus
mono$taxonomy2$result #family
mono$taxonomy3$result #order

head(mono$taxonomy1$result)
plot(na.omit(mono$taxonomy1$result[, 4]))
