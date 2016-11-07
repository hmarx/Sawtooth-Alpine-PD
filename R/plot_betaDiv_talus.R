


####################### Decompose beta diveristy (PhyloSor) to get 'true' turnover (independent of Species Richness) ####################### 


betadiv.dist.alp.PhyloSor.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_PhyloSor.csv")
betadiv.dist.alp.UniFrac.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_UniFrac.csv")

betadiv.dist.alp.PhyloSorTurn.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_PhyloSor_turn.csv")
betadiv.dist.alp.UniFracTurn.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_UniFrac_turn.csv")

betadiv.dist.alp.PhyloSorPD.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_PhyloSor_PD.csv")
betadiv.dist.alp.UniFracPD.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_UniFrac_PD.csv")

betadiv.dist.alp.PhyloSor.matrix <- betadiv.dist.alp.PhyloSor.tmp[,-1]
colnames(betadiv.dist.alp.PhyloSor.matrix) <- betadiv.dist.alp.PhyloSor.tmp$X
rownames(betadiv.dist.alp.PhyloSor.matrix) <- betadiv.dist.alp.PhyloSor.tmp$X
betadiv.dist.alp.PhyloSor.matrix

betadiv.dist.alp.UniFrac.matrix <- betadiv.dist.alp.UniFrac.tmp[,-1]
colnames(betadiv.dist.alp.UniFrac.matrix) <- betadiv.dist.alp.UniFrac.tmp$X
rownames(betadiv.dist.alp.UniFrac.matrix) <- betadiv.dist.alp.UniFrac.tmp$X
betadiv.dist.alp.UniFrac.matrix

betadiv.dist.alp.PhyloSorTurn.matrix <- betadiv.dist.alp.PhyloSorTurn.tmp[,-1]
colnames(betadiv.dist.alp.PhyloSorTurn.matrix) <- betadiv.dist.alp.PhyloSorTurn.tmp$X
rownames(betadiv.dist.alp.PhyloSorTurn.matrix) <- betadiv.dist.alp.PhyloSorTurn.tmp$X
betadiv.dist.alp.PhyloSorTurn.matrix

betadiv.dist.alp.UniFracTurn.matrix <- betadiv.dist.alp.UniFracTurn.tmp[,-1]
colnames(betadiv.dist.alp.UniFracTurn.matrix) <- betadiv.dist.alp.UniFracTurn.tmp$X
rownames(betadiv.dist.alp.UniFracTurn.matrix) <- betadiv.dist.alp.UniFracTurn.tmp$X
betadiv.dist.alp.UniFracTurn.matrix

betadiv.dist.alp.PhyloSorPD.matrix <- betadiv.dist.alp.PhyloSorPD.tmp[,-1]
colnames(betadiv.dist.alp.PhyloSorPD.matrix) <- betadiv.dist.alp.PhyloSorPD.tmp$X
rownames(betadiv.dist.alp.PhyloSorPD.matrix) <- betadiv.dist.alp.PhyloSorPD.tmp$X
betadiv.dist.alp.PhyloSorPD.matrix

betadiv.dist.alp.UniFracPD.matrix <- betadiv.dist.alp.UniFracPD.tmp[,-1]
colnames(betadiv.dist.alp.UniFracPD.matrix) <- betadiv.dist.alp.UniFracPD.tmp$X
rownames(betadiv.dist.alp.UniFracPD.matrix) <- betadiv.dist.alp.UniFracPD.tmp$X
betadiv.dist.alp.UniFracPD.matrix

## aestitic values
#scale = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(98)
#fontsize_row = 9 - nrow(betadiv.dist.alp.PhyloSor.matrix) / 9
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
pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_PhyloSor.pdf", onefile=FALSE)
### Matrix with significant SES values
mat <- matrix(ifelse(betadiv.dist.alp.PhyloSor.matrix > 1.96 | betadiv.dist.alp.PhyloSor.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.alp.PhyloSor.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.alp.PhyloSor.matrix, col=colors2, main="SES PhyloSor", cluster_rows=T, fontsize = fontsize_row * 1.75, 
         fontsize_col = fontsize_row, fontsize_row=fontsize_row, border_color=NA, display_numbers = mat, breaks=bk2, scale="none")
dev.off()

### PhyloSor Significant
length(which(betadiv.dist.alp.PhyloSor.matrix > 1.96))/2 #1
length(which(betadiv.dist.alp.PhyloSor.matrix < -1.96))/2 #5
### TOTAL pairwise comparisons
(nrow(betadiv.dist.alp.PhyloSor.matrix)*(nrow(betadiv.dist.alp.PhyloSor.matrix) - 1))/ 2
## 36

pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_UniFrac.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.alp.UniFrac.matrix > 1.96 | betadiv.dist.alp.UniFrac.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.alp.UniFrac.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.alp.UniFrac.matrix, col=colors2, main="SES UniFrac", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, breaks=bk2)
dev.off()
length(which(betadiv.dist.alp.UniFrac.matrix > 1.96))/2 #0
length(which(betadiv.dist.alp.UniFrac.matrix < -1.96))/2 #24

########## TURN
pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_PhyloSor_turn.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.alp.PhyloSorTurn.matrix > 1.96 | betadiv.dist.alp.PhyloSorTurn.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.alp.PhyloSorTurn.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.alp.PhyloSorTurn.matrix, main="SES PhyloSor_turn", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.alp.PhyloSorTurn.matrix > 1.96))/2 #8
length(which(betadiv.dist.alp.PhyloSorTurn.matrix < -1.96))/2 #7


pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_UniFrac_turn.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.alp.UniFracTurn.matrix > 1.96 | betadiv.dist.alp.UniFracTurn.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.alp.UniFracTurn.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.alp.UniFracTurn.matrix, main="SES UniFrac_turn", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.alp.UniFracTurn.matrix > 1.96))/2 #7
length(which(betadiv.dist.alp.UniFracTurn.matrix < -1.96))/2 #8

######## PD
pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_PhyloSor_PD.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.alp.PhyloSorPD.matrix > 1.96 | betadiv.dist.alp.PhyloSorPD.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.alp.PhyloSorPD.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.alp.PhyloSorPD.matrix, main="SES PhyloSor_PD", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.alp.PhyloSorPD.matrix > 1.96))/2 #8
length(which(betadiv.dist.alp.PhyloSorPD.matrix < -1.96))/2 #16


pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/talus/SES_UniFrac_PD.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.alp.UniFracPD.matrix > 1.96 | betadiv.dist.alp.UniFracPD.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.alp.UniFracPD.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.alp.UniFracPD.matrix, main="SES UniFrac_PD", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.alp.UniFracPD.matrix > 1.96))/2 #10
length(which(betadiv.dist.alp.UniFracPD.matrix < -1.96))/2 #14


