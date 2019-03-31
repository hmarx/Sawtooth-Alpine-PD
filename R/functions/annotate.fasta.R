# Annotates .fasta alignment to updated nomenclature

fastaFile= "output/03_totalData/atpB/atpB.v3.aln.fasta"
#taxonomy = df <- c("updated", "tips")

annotate.fasta <- function(fastaFile, taxonomy){
  GBseqs <- readDNAStringSet(fastaFile) #read .aln.rn
  ## Make a tmp file so don't overwrite 
  tmp.fasta <- DNAStringSet(GBseqs)
  fasta.names <- names(tmp.fasta)
  
  for (i in 1:length(tmp.fasta)){
    #print(i)
    if (names(tmp.fasta[i]) %in% taxonomy[,2]){
    fasta.names[i] = taxonomy[which(names(tmp.fasta[i])==taxonomy[,2]), 1]
        
    }
    else {
      fasta.names[i] = fasta.names[i]
       
    }
      
  }
  names(tmp.fasta) <- fasta.names
  file.name <- strsplit(fastaFile, split=".", fixed=TRUE)[[1]][[1]]
  writeXStringSet(tmp.fasta, file=paste(file.name, "MiSeq", sep=".", format="fasta"))
  return(tmp.fasta)
}
