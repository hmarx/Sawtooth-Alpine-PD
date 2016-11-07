


####################### Decompose beta diveristy (PhyloSor) to get 'true' turnover (independent of Species Richness) ####################### 


betadiv.dist.PhyloSor.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_PhyloSor.csv")
betadiv.dist.UniFrac.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_UniFrac.csv")

betadiv.dist.PhyloSorTurn.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_PhyloSor_turn.csv")
betadiv.dist.UniFracTurn.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_UniFrac_turn.csv")

betadiv.dist.PhyloSorPD.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_PhyloSor_PD.csv")
betadiv.dist.UniFracPD.tmp <- read.csv("output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_UniFrac_PD.csv")

betadiv.dist.PhyloSor.martix <- betadiv.dist.PhyloSor.tmp[,-1]
colnames(betadiv.dist.PhyloSor.martix) <- betadiv.dist.PhyloSor.tmp$X
rownames(betadiv.dist.PhyloSor.martix) <- betadiv.dist.PhyloSor.tmp$X
betadiv.dist.PhyloSor.martix

betadiv.dist.UniFrac.matrix <- betadiv.dist.UniFrac.tmp[,-1]
colnames(betadiv.dist.UniFrac.matrix) <- betadiv.dist.UniFrac.tmp$X
rownames(betadiv.dist.UniFrac.matrix) <- betadiv.dist.UniFrac.tmp$X
betadiv.dist.UniFrac.matrix

betadiv.dist.PhyloSorTurn.matrix <- betadiv.dist.PhyloSorTurn.tmp[,-1]
colnames(betadiv.dist.PhyloSorTurn.matrix) <- betadiv.dist.PhyloSorTurn.tmp$X
rownames(betadiv.dist.PhyloSorTurn.matrix) <- betadiv.dist.PhyloSorTurn.tmp$X
betadiv.dist.PhyloSorTurn.matrix

betadiv.dist.UniFracTurn.matrix <- betadiv.dist.UniFracTurn.tmp[,-1]
colnames(betadiv.dist.UniFracTurn.matrix) <- betadiv.dist.UniFracTurn.tmp$X
rownames(betadiv.dist.UniFracTurn.matrix) <- betadiv.dist.UniFracTurn.tmp$X
betadiv.dist.UniFracTurn.matrix

betadiv.dist.PhyloSorPD.matrix <- betadiv.dist.PhyloSorPD.tmp[,-1]
colnames(betadiv.dist.PhyloSorPD.matrix) <- betadiv.dist.PhyloSorPD.tmp$X
rownames(betadiv.dist.PhyloSorPD.matrix) <- betadiv.dist.PhyloSorPD.tmp$X
betadiv.dist.PhyloSorPD.matrix

betadiv.dist.UniFracPD.matrix <- betadiv.dist.UniFracPD.tmp[,-1]
colnames(betadiv.dist.UniFracPD.matrix) <- betadiv.dist.UniFracPD.tmp$X
rownames(betadiv.dist.UniFracPD.matrix) <- betadiv.dist.UniFracPD.tmp$X
betadiv.dist.UniFracPD.matrix

## aestitic values
#scale = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(98)
#fontsize_row = 10 - nrow(betadiv.dist.PhyloSor.matrix) / 15
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
pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_PhyloSor.pdf", onefile=FALSE)
### Matrix with significant SES values
mat <- matrix(ifelse(betadiv.dist.PhyloSor.matrix > 1.96 | betadiv.dist.PhyloSor.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.PhyloSor.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.PhyloSor.matrix, col=colors2, main="SES PhyloSor", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, breaks=bk2, scale="none")
dev.off()

### PhyloSor Significant
length(which(betadiv.dist.PhyloSor.matrix > 1.96))/2 #1
length(which(betadiv.dist.PhyloSor.matrix < -1.96))/2 #5
### TOTAL pairwise comparisons
(nrow(betadiv.dist.PhyloSor.matrix)*(nrow(betadiv.dist.PhyloSor.matrix) - 1))/ 2
## 36

pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_UniFrac.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.UniFrac.matrix > 1.96 | betadiv.dist.UniFrac.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.UniFrac.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.UniFrac.matrix, col=colors2, main="SES UniFrac", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, breaks=bk2)
dev.off()
length(which(betadiv.dist.UniFrac.matrix > 1.96))/2 #0
length(which(betadiv.dist.UniFrac.matrix < -1.96))/2 #24

########## TURN
pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_PhyloSor_turn.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.PhyloSorTurn.matrix > 1.96 | betadiv.dist.PhyloSorTurn.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.PhyloSorTurn.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.PhyloSorTurn.matrix, main="SES PhyloSor_turn", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.PhyloSorTurn.matrix > 1.96))/2 #8
length(which(betadiv.dist.PhyloSorTurn.matrix < -1.96))/2 #7


pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_UniFrac_turn.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.UniFracTurn.matrix > 1.96 | betadiv.dist.UniFracTurn.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.UniFracTurn.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.UniFracTurn.matrix, main="SES UniFrac_turn", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.UniFracTurn.matrix > 1.96))/2 #7
length(which(betadiv.dist.UniFracTurn.matrix < -1.96))/2 #8

######## PD
pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_PhyloSor_PD.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.PhyloSorPD.matrix > 1.96 | betadiv.dist.PhyloSorPD.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.PhyloSorPD.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.PhyloSorPD.matrix, main="SES PhyloSor_PD", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.PhyloSorPD.matrix > 1.96))/2 #8
length(which(betadiv.dist.PhyloSorPD.matrix < -1.96))/2 #16


pdf(file="output/08_PhyloDiversity/MiSeq/beta/spatial/all/SES_UniFrac_PD.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.dist.UniFracPD.matrix > 1.96 | betadiv.dist.UniFracPD.matrix < -1.96, "*", ""), 
              nrow(betadiv.dist.UniFracPD.matrix))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.dist.UniFracPD.matrix, main="SES UniFrac_PD", cluster_rows=T, 
         fontsize_row=fontsize_row,  fontsize_col = fontsize_row, fontsize=fontsize_row *1.75, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.dist.UniFracPD.matrix > 1.96))/2 #10
length(which(betadiv.dist.UniFracPD.matrix < -1.96))/2 #14


