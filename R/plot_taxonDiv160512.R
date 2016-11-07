
head(taxonomy.table)
head(sawtooth.com.tax)

##### Diveristy of Taxonomic groups on summits PHLAWD

saw.com.taxonomy <- merge(x=sawtooth.com.phlawd.tax, y=taxonomy.table[1:4], by.x ="genus", by.y="X", all.x= T, sort = T)
head(saw.com.taxonomy)

saw.com.taxonomy.df <- gather(data =saw.com.taxonomy, "summit", "presence", 2:10)
head(saw.com.taxonomy.df)

saw.com.taxonomy.df.pres <- saw.com.taxonomy.df[saw.com.taxonomy.df$presence !=0, ]
unique(saw.com.taxonomy.df.pres$summit)

## Non-vasular plants
saw.com.taxonomy.df.pres[is.na(saw.com.taxonomy.df.pres$Spermatophyta),]
saw.com.taxonomy.df.pres.vasc <- saw.com.taxonomy.df.pres[!is.na(saw.com.taxonomy.df.pres$Spermatophyta),]
saw.com.taxonomy.df.pres.vasc$summit <- gsub(saw.com.taxonomy.df.pres.vasc$summit, pattern = "[.]", replacement = " ")

newOrder <- rev(names(sort(summary(as.factor(saw.com.taxonomy.df.pres.vasc$summit)))))
saw.com.taxonomy.df.pres.vasc$summit <- factor(saw.com.taxonomy.df.pres.vasc$summit, levels = newOrder)

levels(saw.com.taxonomy.df.pres.vasc$order)[levels(saw.com.taxonomy.df.pres.vasc$order) == ""] <- "Not Assigned"

summit.orders <- ggplot(saw.com.taxonomy.df.pres.vasc, aes(x=summit, fill = order)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  coord_cartesian(ylim=c(0,80))
ggsave("figs/orders.PHLAWD.pdf", summit.orders)

master.tax <- rbind(cbind(saw.com.collect.taxonomy.df.pres[c(1,17:21)], method = rep("Collected")), cbind(saw.com.taxonomy.df.pres.vasc, method = rep("PHLAWD")),
cbind(saw.com.miseq.taxonomy.df.pres[c(1,17:21)], method = rep("MiSeq")))

master.tax.plot <- ggplot(master.tax, aes(x=summit, fill = order)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(~ method) +
  coord_cartesian(ylim=c(0,80))
ggsave("figs/orders.pdf", master.tax.plot, width = 11, height = 8)

##### Diveristy of Taxonomic groups on summits MiSeq Data

saw.com.miseq.taxonomy <- merge(x=sawtooth.com.miseq.tax, y=taxonomy.table[1:4], by.x ="genus", by.y="X", all.x= T, sort = T)
head(saw.com.miseq.taxonomy)
colnames(saw.com.miseq.taxonomy)
dim(saw.com.miseq.taxonomy)

saw.com.miseq.taxonomy.df <- gather(data =saw.com.miseq.taxonomy, "summit", "presence", 17:25)
saw.com.miseq.taxonomy.df.pres <- saw.com.miseq.taxonomy.df[saw.com.miseq.taxonomy.df$presence !=0, ]
head(saw.com.miseq.taxonomy.df.pres)
str(saw.com.miseq.taxonomy.df.pres)
saw.com.miseq.taxonomy.df.pres$summit <- gsub(saw.com.miseq.taxonomy.df.pres$summit, pattern = "[.]", replacement = " ")
saw.com.miseq.taxonomy.df.pres$summit <- gsub(saw.com.miseq.taxonomy.df.pres$summit, pattern = "_", replacement = " ")
saw.com.miseq.taxonomy.df.pres$summit <- as.factor(saw.com.miseq.taxonomy.df.pres$summit)
newOrder <- rev(names(sort(summary(as.factor(saw.com.miseq.taxonomy.df.pres$summit)))))
saw.com.miseq.taxonomy.df.pres$summit <- factor(saw.com.miseq.taxonomy.df.pres$summit, levels = newOrder)

saw.com.miseq.taxonomy.df.pres[is.na(saw.com.miseq.taxonomy.df.pres$Spermatophyta),]
saw.com.miseq.taxonomy.df.pres[saw.com.miseq.taxonomy.df.pres$genus == "Unknown",]

levels(saw.com.miseq.taxonomy.df.pres$order)[levels(saw.com.miseq.taxonomy.df.pres$order) == ""] <- "Not Assigned"
saw.com.miseq.taxonomy.df.pres$order[is.na(saw.com.miseq.taxonomy.df.pres$order)] <- "Not Assigned"
newOrderClade <- rev(names(sort(summary(as.factor(saw.com.miseq.taxonomy.df.pres$order)))))
saw.com.miseq.taxonomy.df.pres$order <- factor(saw.com.miseq.taxonomy.df.pres$order, levels = newOrderClade)

summit.orders <- ggplot(saw.com.miseq.taxonomy.df.pres, aes(x=summit, fill = order)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 
ggsave("figs/orders.miseq.pdf", summit.orders)

summary(saw.com.miseq.taxonomy.df.pres$order)


## Just Alpine species
saw.com.miseq.taxonomy.df.pres.alpine <- saw.com.miseq.taxonomy.df.pres[saw.com.miseq.taxonomy.df.pres$Meadow == 0, ]

newOrderClade <- rev(names(sort(summary(as.factor(saw.com.miseq.taxonomy.df.pres.alpine$order)))))
saw.com.miseq.taxonomy.df.pres.alpine$order <- factor(saw.com.miseq.taxonomy.df.pres.alpine$order, levels = newOrderClade)
newOrder <- rev(names(sort(summary(as.factor(saw.com.miseq.taxonomy.df.pres.alpine$summit)))))
saw.com.miseq.taxonomy.df.pres.alpine$summit <- factor(saw.com.miseq.taxonomy.df.pres.alpine$summit, levels = newOrder)

summit.orders.alpine <- ggplot(saw.com.miseq.taxonomy.df.pres.alpine, aes(x=summit, fill = order)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 
ggsave("figs/orders.miseq.talus.pdf", summit.orders.alpine)


##### Diveristy of Taxonomic groups on summits Collected
saw.com.collect.taxonomy <- merge(x=sawtooth.com.collect.tax, y=taxonomy.table[1:4], by.x ="genus", by.y="X", all.x= T, sort = T)
head(saw.com.collect.taxonomy)
colnames(saw.com.collect.taxonomy)


saw.com.collect.taxonomy.df <- gather(data =saw.com.collect.taxonomy, "summit", "presence", 17:25)
saw.com.collect.taxonomy.df.pres <- saw.com.collect.taxonomy.df[saw.com.collect.taxonomy.df$presence !=0, ]
head(saw.com.collect.taxonomy.df.pres)
str(saw.com.collect.taxonomy.df.pres)
saw.com.collect.taxonomy.df.pres$summit <- gsub(saw.com.collect.taxonomy.df.pres$summit, pattern = "[.]", replacement = " ")
saw.com.collect.taxonomy.df.pres$summit <- gsub(saw.com.collect.taxonomy.df.pres$summit, pattern = "_", replacement = " ")
saw.com.collect.taxonomy.df.pres$summit <- as.factor(saw.com.collect.taxonomy.df.pres$summit)
newOrder <- rev(names(sort(summary(as.factor(saw.com.collect.taxonomy.df.pres$summit)))))
saw.com.collect.taxonomy.df.pres$summit <- factor(saw.com.collect.taxonomy.df.pres$summit, levels = newOrder)

saw.com.collect.taxonomy.df.pres[is.na(saw.com.collect.taxonomy.df.pres$Spermatophyta),]
saw.com.collect.taxonomy.df.pres[saw.com.collect.taxonomy.df.pres$genus == "Unknown",]

levels(saw.com.collect.taxonomy.df.pres$order)[levels(saw.com.collect.taxonomy.df.pres$order) == ""] <- "Not Assigned"
saw.com.collect.taxonomy.df.pres$order[is.na(saw.com.collect.taxonomy.df.pres$order)] <- "Not Assigned"
newOrderClade <- rev(names(sort(summary(as.factor(saw.com.collect.taxonomy.df.pres$order)))))
saw.com.collect.taxonomy.df.pres$order <- factor(saw.com.collect.taxonomy.df.pres$order, levels = newOrderClade)

summit.orders <- ggplot(saw.com.collect.taxonomy.df.pres, aes(x=summit, fill = order)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 
ggsave("figs/orders.collect.pdf", summit.orders)

## Just Alpine species
saw.com.collect.taxonomy.df.pres.alpine <- saw.com.collect.taxonomy.df.pres[saw.com.collect.taxonomy.df.pres$Meadow == 0, ]

newOrderClade <- rev(names(sort(summary(as.factor(saw.com.collect.taxonomy.df.pres.alpine$order)))))
saw.com.collect.taxonomy.df.pres.alpine$order <- factor(saw.com.collect.taxonomy.df.pres.alpine$order, levels = newOrderClade)
newOrder <- rev(names(sort(summary(as.factor(saw.com.collect.taxonomy.df.pres.alpine$summit)))))
saw.com.collect.taxonomy.df.pres.alpine$summit <- factor(saw.com.collect.taxonomy.df.pres.alpine$summit, levels = newOrder)

summit.orders.alpine <- ggplot(saw.com.collect.taxonomy.df.pres.alpine, aes(x=summit, fill = order)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 
ggsave("figs/orders.collect.talus.pdf", summit.orders.alpine)

summary(saw.com.collect.taxonomy.df.pres$order)


##### Diveristy of Taxonomic groups on summits in beta community matrix
sawtooth.com.miseq.beta

split2 <- strsplit(as.character(rownames(sawtooth.com.miseq.beta)), split="_", fixed=TRUE) #split names
genus.name2 <- sapply(split2, "[", 1L)
sawtooth.com.beta.tax  <- cbind("genus"=genus.name2, sawtooth.com.miseq.beta)
head(sawtooth.com.beta.tax)
dim(sawtooth.com.beta.tax) #143

saw.com.beta.taxonomy <- merge(x=sawtooth.com.beta.tax, y=taxonomy.table[1:4], by.x ="genus", by.y="X", all.x= T, sort = T)
head(saw.com.beta.taxonomy)
colnames(saw.com.beta.taxonomy)
dim(saw.com.beta.taxonomy)

saw.com.beta.taxonomy.df <- gather(data =saw.com.beta.taxonomy, "summit", "presence", 2:10)
saw.com.beta.taxonomy.df.pres <- saw.com.beta.taxonomy.df[saw.com.beta.taxonomy.df$presence !=0, ]
head(saw.com.beta.taxonomy.df.pres)
str(saw.com.beta.taxonomy.df.pres)
saw.com.beta.taxonomy.df.pres$summit <- gsub(saw.com.beta.taxonomy.df.pres$summit, pattern = "_", replacement = " ")
saw.com.beta.taxonomy.df.pres$summit <- as.factor(saw.com.beta.taxonomy.df.pres$summit)
newOrder <- rev(names(sort(summary(as.factor(saw.com.beta.taxonomy.df.pres$summit)))))
saw.com.beta.taxonomy.df.pres$summit <- factor(saw.com.beta.taxonomy.df.pres$summit, levels = newOrder)

saw.com.beta.taxonomy.df.pres[is.na(saw.com.beta.taxonomy.df.pres$Spermatophyta),]
saw.com.beta.taxonomy.df.pres[saw.com.beta.taxonomy.df.pres$genus == "Unknown",]

levels(saw.com.beta.taxonomy.df.pres$order)[levels(saw.com.beta.taxonomy.df.pres$order) == ""] <- "Not Assigned"
saw.com.beta.taxonomy.df.pres$order[is.na(saw.com.beta.taxonomy.df.pres$order)] <- "Not Assigned"
newOrderClade <- rev(names(sort(summary(as.factor(saw.com.beta.taxonomy.df.pres$order)))))
saw.com.beta.taxonomy.df.pres$order <- factor(saw.com.beta.taxonomy.df.pres$order, levels = newOrderClade)

summit.orders <- ggplot(saw.com.beta.taxonomy.df.pres, aes(x=summit, fill = order)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 
ggsave("figs/orders.beta.pdf", summit.orders)

summary(saw.com.beta.taxonomy.df.pres$order)

## Just Alpine species
split2 <- strsplit(as.character(rownames(sawtooth.com.miseq.talus.beta)), split="_", fixed=TRUE) #split names
genus.name2 <- sapply(split2, "[", 1L)
sawtooth.com.beta.talus.tax  <- cbind("genus"=genus.name2, sawtooth.com.miseq.talus.beta)
head(sawtooth.com.beta.talus.tax)
dim(sawtooth.com.beta.talus.tax) #114

saw.com.beta.talus.taxonomy <- merge(x=sawtooth.com.beta.talus.tax, y=taxonomy.table[1:4], by.x ="genus", by.y="X", all.x= T, sort = T)
head(saw.com.beta.talus.taxonomy)
colnames(saw.com.beta.talus.taxonomy)
dim(saw.com.beta.talus.taxonomy)

saw.com.beta.talus.taxonomy.df <- gather(data =saw.com.beta.talus.taxonomy, "summit", "presence", 2:10)
saw.com.beta.talus.taxonomy.df.pres <- saw.com.beta.talus.taxonomy.df[saw.com.beta.talus.taxonomy.df$presence !=0, ]
head(saw.com.beta.talus.taxonomy.df.pres)
str(saw.com.beta.talus.taxonomy.df.pres)
saw.com.beta.talus.taxonomy.df.pres$summit <- gsub(saw.com.beta.talus.taxonomy.df.pres$summit, pattern = "_", replacement = " ")
saw.com.beta.talus.taxonomy.df.pres$summit <- as.factor(saw.com.beta.talus.taxonomy.df.pres$summit)
newOrder <- rev(names(sort(summary(as.factor(saw.com.beta.talus.taxonomy.df.pres$summit)))))
saw.com.beta.talus.taxonomy.df.pres$summit <- factor(saw.com.beta.talus.taxonomy.df.pres$summit, levels = newOrder)

saw.com.beta.talus.taxonomy.df.pres[is.na(saw.com.beta.talus.taxonomy.df.pres$Spermatophyta),]
saw.com.beta.talus.taxonomy.df.pres[saw.com.beta.talus.taxonomy.df.pres$genus == "Unknown",]

levels(saw.com.beta.talus.taxonomy.df.pres$order)[levels(saw.com.beta.talus.taxonomy.df.pres$order) == ""] <- "Not Assigned"
saw.com.beta.talus.taxonomy.df.pres$order[is.na(saw.com.beta.talus.taxonomy.df.pres$order)] <- "Not Assigned"
newOrderClade <- rev(names(sort(summary(as.factor(saw.com.beta.talus.taxonomy.df.pres$order)))))
saw.com.beta.talus.taxonomy.df.pres$order <- factor(saw.com.beta.talus.taxonomy.df.pres$order, levels = newOrderClade)

summit.orders <- ggplot(saw.com.beta.talus.taxonomy.df.pres, aes(x=summit, fill = order)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 
ggsave("figs/orders.beta.talus.pdf", summit.orders)

summary(saw.com.beta.talus.taxonomy.df.pres$order)



