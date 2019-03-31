
library(Biostrings) # load package

# This function will take the full alignment from the total dataset (PHLAWD + MiSeq)
# that has been cleaned. It will remove the NCBI ID, and keep only the longest unique sequences 
# if there are multiple hits for a single species
#fasta.file <- "output/04_Clean/totalData/atpB.total.align.clean"

parsePHLAWDtotaldata <- function(fasta.file, splist.raw){
  GBseqs <- readDNAStringSet(fasta.file) #read .aln.full
  namesGB <- names(GBseqs) #get the full NCBI names
  print(length(namesGB))
  split <- strsplit(namesGB, split="|", fixed=TRUE) #split names
  for (i in 1:length(split)){
    if (split[i][[1]][[1]] == "ncbi"){
      genus.sp.name <- sapply(split[i], "[", 2L) 
      ncbi <- sapply(split[i], "[", 4L) 
      combinedname <- paste(genus.sp.name, "ncbi", ncbi, sep="|") 
      split[i] <- combinedname
    } else if (split[i][[1]][[1]] != "ncbi"){
      split2 <- strsplit(split[i][[1]], split="_", fixed=TRUE) 
      genus.sp.name <- paste(head(split2[[1]], (length(split2[[1]]) -1)), collapse ="_")
      combinedname <- paste(genus.sp.name, "miseq", tail(split2[[1]],1), sep="|")
      split[i] <- combinedname
      
    }
  }
  split <- gsub(pattern = "__", replacement = "_", split)
  split <- gsub(pattern = "cf_", replacement = "", split)
  names(GBseqs) <- split
    sizes <- rowSums(alphabetFrequency(GBseqs)[,c("A","C","T","G")]) #get the nucleotide lenght of each sequence
    ord <- order(names(GBseqs), -sizes)
    seqs <- GBseqs[ord] #order by length of sequence, longest first
    namesGBord <- names(seqs) #get the full NCBI names in correct order
    combinedname <- sapply(strsplit(names(GBseqs)[ord], split="|", fixed=TRUE), "[", 1L)
    ID <- duplicated(combinedname) # identify duplicated combined names
    uniques <- seqs[!ID] #get only the unique sequences, choosing the longest since it is the first in the list
    uniquesnames <- combinedname[!ID]
    
    # Get only species with full names in list
    sawspec <- uniques[sapply(strsplit(names(uniques), split="|", fixed=TRUE), "[", 1L)  %in% splist.raw]
    
    print(length(sawspec))
    file.name <- strsplit(fasta.file, split=".", fixed=TRUE)[[1]][[1]]
    species_uniques <- sawspec 
    names(species_uniques) <- sapply(strsplit(names(sawspec), split="|", fixed=TRUE), "[", 1L) 
    writeXStringSet(species_uniques, file=paste(file.name, "total.unique", sep=".", format="fasta"))
    #names(uniques) <- namesGBord[!ID] #full NCBI names
    writeXStringSet(sawspec, file=paste(file.name, "total.unique.GB", sep=".", format="fasta"))
  }



