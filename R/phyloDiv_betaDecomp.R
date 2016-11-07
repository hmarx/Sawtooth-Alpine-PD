
############################ Turnover of species between summits

###### Beta diveristy of Sawtooth National Forest Summits 
## Source pool = Alpine Summits (talus + meadow)
## Clade = Spermatophyta

source("R/beta_pd_decomp.r")

decompo_beta <- beta.pd.decompo(com = sawMiseqBeta$comm, tree = sawMiseqBeta$phy, type="both", output.dist=T, random=999)

decompo_beta_df <- beta.pd.decompo(com = sawMiseqBeta$comm, tree = sawMiseqBeta$phy, type="both", output.dist=F, random=999)
write.csv(decompo_beta_df, "output/08_PhyloDiversity/MiSeq/beta/spatial/all/decompo_beta_df.csv")

#### Distance Matricies take forever to run
write.csv(data.matrix(decompo_beta$betadiv$PhyloSor), file= "output/08_PhyloDiversity/MiSeq/beta/spatial/all/dist_PhloSor.csv")

write.csv(data.matrix(decompo_beta$betadiv$PhyloSor_turn), file="output/08_PhyloDiversity/MiSeq/beta/spatial/all/dist_PhyloSor_turn.csv")

write.csv(data.matrix(decompo_beta$betadiv$PhyloSor_PD), file="output/08_PhyloDiversity/MiSeq/beta/spatial/all/dist_PhyloSor_PD.csv")

write.csv(data.matrix(decompo_beta$betadiv$UniFrac), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/all/dist_UniFrac.csv")

write.csv(data.matrix(decompo_beta$betadiv$UniFrac_turn), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/all/dist_UniFrac_turn.csv")

write.csv(data.matrix(decompo_beta$betadiv$UniFrac_PD), file= "output/08_PhyloDiversity/MiSeq/beta/spatial/all/dist_UniFrac_PD.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_PhyloSor.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor_turn), file= "output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_PhyloSor_turn.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor_PD), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_PhyloSor_PD.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_UniFrac.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac_turn), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_UniFrac_turn.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac_PD), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_UniFrac_PD.csv")

###############
###### Beta diveristy of Sawtooth National Forest Summits 
## Source pool = Alpine Summits (talus)
## Clade = Spermatophyta

source("R/beta_pd_decomp.r")

decompo_beta_talus <- beta.pd.decompo(com = sawMiseqBetaAlp$comm, tree = sawMiseqBetaAlp$phy, type="both", output.dist=T, random=999)

decompo_beta_talus_df <- beta.pd.decompo(com = sawMiseqBetaAlp$comm, tree = sawMiseqBetaAlp$phy, type="both", output.dist=F, random=999)
write.csv(decompo_beta_talus_df, "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/decompo_beta_talus_df.csv")

#### Distance Matricies take forever to run
write.csv(data.matrix(decompo_beta_talus$betadiv$PhyloSor), file= "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/dist_PhloSor.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$PhyloSor_turn), file="output/08_PhyloDiversity/MiSeq/beta/spatial/talus/dist_PhyloSor_turn.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$PhyloSor_PD), file="output/08_PhyloDiversity/MiSeq/beta/spatial/talus/dist_PhyloSor_PD.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$UniFrac), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/dist_UniFrac.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$UniFrac_turn), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/dist_UniFrac_turn.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$UniFrac_PD), file= "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/dist_UniFrac_PD.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$SES_PhyloSor), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_PhyloSor.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$SES_PhyloSor_turn), file= "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_PhyloSor_turn.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$SES_PhyloSor_PD), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_PhyloSor_PD.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$SES_UniFrac), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_UniFrac.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$SES_UniFrac_turn), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_UniFrac_turn.csv")

write.csv(data.matrix(decompo_beta_talus$betadiv$SES_UniFrac_PD), file = "output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_UniFrac_PD.csv")

