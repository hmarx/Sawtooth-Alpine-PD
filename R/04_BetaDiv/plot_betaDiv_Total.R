#################################################################################################################
######## Plots of Phylo Beta Diversity (spatial) ################################################################
######## Hannah E. Marx, 20 Jan 2019 ############################################################################
#################################################################################################################
####### Decompose beta diveristy (PhyloSor) to get 'true' turnover (independent of Species Richness) ############

################################# Total dataset ################################# 
betadiv.dist.total.UniFrac.tmp <- read.csv("output/10_PhyloDiversity/beta/spatial/all/Total_SES_UniFrac.csv")
betadiv.dist.total.UniFracTurn.tmp <- read.csv("output/10_PhyloDiversity/beta/spatial/all/Total_SES_UniFrac_turn.csv")
betadiv.dist.total.UniFracPD.tmp <- read.csv("output/10_PhyloDiversity/beta/spatial/all/Total_SES_UniFrac_PD.csv")

betadiv.dist.total.UniFrac.matrix <- betadiv.dist.total.UniFrac.tmp[,-1]
colnames(betadiv.dist.total.UniFrac.matrix) <- betadiv.dist.total.UniFrac.tmp$X
rownames(betadiv.dist.total.UniFrac.matrix) <- betadiv.dist.total.UniFrac.tmp$X
betadiv.dist.total.UniFrac.matrix

betadiv.dist.total.UniFracTurn.matrix <- betadiv.dist.total.UniFracTurn.tmp[,-1]
colnames(betadiv.dist.total.UniFracTurn.matrix) <- betadiv.dist.total.UniFracTurn.tmp$X
rownames(betadiv.dist.total.UniFracTurn.matrix) <- betadiv.dist.total.UniFracTurn.tmp$X
betadiv.dist.total.UniFracTurn.matrix

betadiv.dist.total.UniFracPD.matrix <- betadiv.dist.total.UniFracPD.tmp[,-1]
colnames(betadiv.dist.total.UniFracPD.matrix) <- betadiv.dist.total.UniFracPD.tmp$X
rownames(betadiv.dist.total.UniFracPD.matrix) <- betadiv.dist.total.UniFracPD.tmp$X
betadiv.dist.total.UniFracPD.matrix

## aestitic values
#scale = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(98)
#fontsize_row = 10 - nrow(betadiv.dist.total.PhyloSor.matrix) / 15
#bk2 = unique(seq(-3,3, length=100))
fontsize_row = 15

#create the breaks
bk2 = unique(c(seq(-6, -.0001, length=100), 0, seq(.0001, 6, length=100)))

#set different color vectors for each interval
col1 = colorRampPalette(rev(brewer.pal(n = 9, name ="Blues")))(99) #set the order of greys
col2 <- rep("white", 1)
col3 = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(99)
colors2 <- c(col1, col2, col3)

######## Undecomposed
pdf(file="figs/beta/Total_SES_UniFrac.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.total.UniFrac.matrix > 1.96 | betadiv.dist.total.UniFrac.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.total.UniFrac.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.total.UniFrac.matrix, col=colors2, main="SES UniFrac", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, breaks=bk2)
dev.off()
length(which(betadiv.dist.total.UniFrac.matrix > 1.96))/2 #2
length(which(betadiv.dist.total.UniFrac.matrix < -1.96))/2 #0

########## TURN
pdf(file="figs/beta/Total_SES_UniFrac_turn.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.total.UniFracTurn.matrix > 1.96 | betadiv.dist.total.UniFracTurn.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.total.UniFracTurn.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.total.UniFracTurn.matrix, main="SES UniFrac_turn", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.total.UniFracTurn.matrix > 1.96))/2 #4
length(which(betadiv.dist.total.UniFracTurn.matrix < -1.96))/2 #0

######## PD
pdf(file="figs/beta/Total_SES_UniFrac_PD.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.total.UniFracPD.matrix > 1.96 | betadiv.dist.total.UniFracPD.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.total.UniFracPD.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.total.UniFracPD.matrix, main="SES UniFrac_PD", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.total.UniFracPD.matrix > 1.96))/2 #1
length(which(betadiv.dist.total.UniFracPD.matrix < -1.96))/2 #2


################################# Total dataset: Talus only ################################# 
betadiv.dist.total.UniFrac.talus.tmp <- read.csv("output/10_PhyloDiversity/beta/spatial/talus/Total_SES_UniFrac.csv")
betadiv.dist.total.UniFracTurn.talus.tmp <- read.csv("output/10_PhyloDiversity/beta/spatial/talus/Total_SES_UniFrac_turn.csv")
betadiv.dist.total.UniFracPD.talus.tmp <- read.csv("output/10_PhyloDiversity/beta/spatial/talus/Total_SES_UniFrac_PD.csv")

betadiv.dist.total.UniFrac.talus.matrix <- betadiv.dist.total.UniFrac.talus.tmp[,-1]
colnames(betadiv.dist.total.UniFrac.talus.matrix) <- betadiv.dist.total.UniFrac.talus.tmp$X
rownames(betadiv.dist.total.UniFrac.talus.matrix) <- betadiv.dist.total.UniFrac.talus.tmp$X
betadiv.dist.total.UniFrac.talus.matrix

betadiv.dist.total.UniFracTurn.talus.martrix <- betadiv.dist.total.UniFracTurn.talus.tmp[,-1]
colnames(betadiv.dist.total.UniFracTurn.talus.martrix) <- betadiv.dist.total.UniFracTurn.talus.tmp$X
rownames(betadiv.dist.total.UniFracTurn.talus.martrix) <- betadiv.dist.total.UniFracTurn.talus.tmp$X
betadiv.dist.total.UniFracTurn.talus.martrix

betadiv.dist.total.UniFracPD.talus.matrix <- betadiv.dist.total.UniFracPD.talus.tmp[,-1]
colnames(betadiv.dist.total.UniFracPD.talus.matrix) <- betadiv.dist.total.UniFracPD.talus.tmp$X
rownames(betadiv.dist.total.UniFracPD.talus.matrix) <- betadiv.dist.total.UniFracPD.talus.tmp$X
betadiv.dist.total.UniFracPD.talus.matrix

## aestitic values
#scale = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(98)
#fontsize_row = 10 - nrow(betadiv.dist.total.PhyloSor.matrix) / 15
#bk2 = unique(seq(-3,3, length=100))
fontsize_row = 15

#create the breaks
bk2 = unique(c(seq(-6, -.0001, length=100), 0, seq(.0001, 6, length=100)))

#set different color vectors for each interval
col1 = colorRampPalette(rev(brewer.pal(n = 9, name ="Blues")))(99) #set the order of greys
col2 <- rep("white", 1)
col3 = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(99)
colors2 <- c(col1, col2, col3)


######## Undecomposed
pdf(file="figs/beta/Total_SES_UniFrac.talus.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.total.UniFrac.talus.matrix > 1.96 | betadiv.dist.total.UniFrac.talus.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.total.UniFrac.talus.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.total.UniFrac.talus.matrix, col=colors2, main="SES UniFrac", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, breaks=bk2)
dev.off()
length(which(betadiv.dist.total.UniFrac.talus.matrix > 1.96))/2 #4
length(which(betadiv.dist.total.UniFrac.talus.matrix < -1.96))/2 #0

########## TURN
pdf(file="figs/beta/Total_SES_UniFrac_turn.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.total.UniFracTurn.talus.martrix > 1.96 | betadiv.dist.total.UniFracTurn.talus.martrix < -1.96, "*", ""), 
              nrow(betadiv.dist.total.UniFracTurn.talus.martrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.total.UniFracTurn.talus.martrix, main="SES UniFrac_turn", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.total.UniFracTurn.talus.martrix > 1.96))/2 #0
length(which(betadiv.dist.total.UniFracTurn.talus.martrix < -1.96))/2 #0

######## PD
pdf(file="figs/beta/Total_SES_UniFrac_PD.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.total.UniFracPD.talus.matrix > 1.96 | betadiv.dist.total.UniFracPD.talus.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.total.UniFracPD.talus.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.total.UniFracPD.talus.matrix, main="SES UniFrac_PD", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.total.UniFracPD.talus.matrix > 1.96))/2 #3
length(which(betadiv.dist.total.UniFracPD.talus.matrix < -1.96))/2 #0

################## Symmetric matrix:


cor_matrix <- lowerUpper(lower = betadiv.dist.total.UniFracTurn.talus.martrix, upper = betadiv.dist.total.UniFracTurn.matrix, diff=FALSE)
mat <- matrix(ifelse(cor_matrix > 1.96 | cor_matrix < -1.96, "*", ""), 
              nrow(cor_matrix))
mat[is.na(mat)] <-  ""
pdf(file="figs/beta/Total_SES_UniFrac_All_Talus.pdf", onefile=FALSE)
pheatmap(cor_matrix, main="SES UniFrac", cluster_rows=F,  cluster_cols = F,
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()






