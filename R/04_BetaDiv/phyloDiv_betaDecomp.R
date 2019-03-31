#################################################################################################################
######## Phylo Beta Diversity (spacial) #########################################################################
######## Hannah E. Marx, 20 Jan 2019 ############################################################################
#################################################################################################################

source("R/functions/beta_pd_decomp.r")

############################ Turnover of species between summits
###### Beta diveristy of Sawtooth National Forest Summits 
## Clade = Spermatophyta

###########################################  Total dataset ########################################### 
################# Source pool = Alpine Summits (talus + meadow)
decompo_beta_total <- beta.pd.decompo(com = sawTotal$comm, tree = sawTotal$phy, type="both", output.dist=T, random=999)
decompo_beta_total_df <- beta.pd.decompo(com = sawTotal$comm, tree = sawTotal$phy, type="both", output.dist=F, random=999)
write.csv(decompo_beta_total_df, "output/10_PhyloDiversity/beta/spatial/all/Total_decompo_beta_total_df.csv")

#### Distance Matricies take forever to run
write.csv(data.matrix(decompo_beta_total$betadiv$UniFrac), file = "output/10_PhyloDiversity/beta/spatial/all/Total_dist_UniFrac.csv")
write.csv(data.matrix(decompo_beta_total$betadiv$UniFrac_turn), file = "output/10_PhyloDiversity/beta/spatial/all/Total_dist_UniFrac_turn.csv")
write.csv(data.matrix(decompo_beta_total$betadiv$UniFrac_PD), file= "output/10_PhyloDiversity/beta/spatial/all/Total_dist_UniFrac_PD.csv")
write.csv(data.matrix(decompo_beta_total$betadiv$SES_UniFrac), file = "output/10_PhyloDiversity/beta/spatial/all/Total_SES_UniFrac.csv")
write.csv(data.matrix(decompo_beta_total$betadiv$SES_UniFrac_turn), file = "output/10_PhyloDiversity/beta/spatial/all/Total_SES_UniFrac_turn.csv")
write.csv(data.matrix(decompo_beta_total$betadiv$SES_UniFrac_PD), file = "output/10_PhyloDiversity/beta/spatial/all/Total_SES_UniFrac_PD.csv")

################# Source pool = Alpine Summits (talus)
decompo_beta_total_talus <- beta.pd.decompo(com = sawTotal.alpine$comm, tree = sawTotal.alpine$phy, type="both", output.dist=T, random=999)
decompo_beta_total_talus_df <- beta.pd.decompo(com = sawTotal.alpine$comm, tree = sawTotal.alpine$phy, type="both", output.dist=F, random=999)
write.csv(decompo_beta_total_talus_df, "output/10_PhyloDiversity/beta/spatial/talus/Total_decompo_beta_total_talus_df.csv")

write.csv(data.matrix(decompo_beta_total_talus$betadiv$UniFrac), file = "output/10_PhyloDiversity/beta/spatial/talus/Total_dist_UniFrac.csv")
write.csv(data.matrix(decompo_beta_total_talus$betadiv$UniFrac_turn), file = "output/10_PhyloDiversity/beta/spatial/talus/Total_dist_UniFrac_turn.csv")
write.csv(data.matrix(decompo_beta_total_talus$betadiv$UniFrac_PD), file= "output/10_PhyloDiversity/beta/spatial/talus/Total_dist_UniFrac_PD.csv")
write.csv(data.matrix(decompo_beta_total_talus$betadiv$SES_UniFrac), file = "output/10_PhyloDiversity/beta/spatial/talus/Total_SES_UniFrac.csv")
write.csv(data.matrix(decompo_beta_total_talus$betadiv$SES_UniFrac_turn), file = "output/10_PhyloDiversity/beta/spatial/talus/Total_SES_UniFrac_turn.csv")
write.csv(data.matrix(decompo_beta_total_talus$betadiv$SES_UniFrac_PD), file = "output/10_PhyloDiversity/beta/spatial/talus/Total_SES_UniFrac_PD.csv")

