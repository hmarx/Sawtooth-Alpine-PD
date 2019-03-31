## Last updated August 28, 2013
##With help from Matt Settles / Matt Pennell; Version of NEWfunction.R
##Modified May15,2013 by HEM to correct NCBI names, and add "_"

# Will have to install this package 
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

library(Biostrings) # load package

#setwd("~/Dropbox/Hannah-Dave/SanJuans/HannahFINAL/2_SpeciesList/") # Navigate to the directory with PHLAWD output to be parsed

# This function will take the full alignment from the PHLAWD output and remove the NCBI ID, 
# and keep only the longest unique sequences if there are multiple hits for a single species
parsePHLAWD <- function(fasta.file){
  GBseqs <- readDNAStringSet(fasta.file) #read .aln.full
  namesGB <- names(GBseqs) #get the full NCBI names
  print(length(namesGB))
  split <- strsplit(namesGB, split="_", fixed=TRUE) #split names
  genus.name <- sapply(split, "[", 1L) 
  species.name2 <- sapply(split, "[", 2L) 
  combinedname <- paste(genus.name, species.name2, sep="_") #get just genus_species
  sizes <- rowSums(alphabetFrequency(GBseqs)[,c("A","C","T","G")]) #get the nucleotide lenght of each sequence
  ord <- order(combinedname, -sizes)
  seqs <- GBseqs[ord] #order by lenght of sequence, longest first
  namesGBord <- names(seqs) #get the full NCBI names in correct order
  combinedname <- combinedname[ord]
  ID <- duplicated(combinedname) # identify duplicated combined names
  uniques <- seqs[!ID] #get only the unique sequences, choosing the longest since it is the first in the list
  uniquesnames <- combinedname[!ID]
  print(length(uniques))
  file.name <- strsplit(fasta.file, split=".", fixed=TRUE)[[1]][[1]]
  species_uniques <- uniques 
  names(species_uniques) <- uniquesnames
  writeXStringSet(species_uniques, file=paste(file.name, "unique", sep=".", format="fasta"))

}  

## To execute, run the above funtion, then call the file that you would like to parse. See the example atpB.FINAL.aln.full below:
#parsePHLAWD("atpB.FINAL.aln.rn") 
## Output: *unique.fasta == the alignment trimed to just the longest sequences, i.e. the unique seqs
