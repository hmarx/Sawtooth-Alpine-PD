# Annotates .fasta alignment to updated nomenclature

#fastaFile= "output/03_totalData/atpB/atpB.v3.aln.fasta"
#taxonomy = df <- c("updated", "tips")

annotate.fasta <- function(fastaFile, taxonomy){
  GBseqs <- readDNAStringSet(fastaFile) #read .aln.rn
  ## Make a tmp file so don't overwrite 
  tmp.fasta <- DNAStringSet(GBseqs)
  fasta.names <- names(tmp.fasta)
  
  for (i in 1:length(tmp.fasta)){
    #print(i)
    if (names(tmp.fasta[i]) %in% taxonomy[,2]){
    fasta.names[i] = as.character(taxonomy[which(names(tmp.fasta[i])==taxonomy[,2]),1])
        
    }
    else {
      fasta.names[i] = fasta.names[i]
       
    }
      
  }
  names(tmp.fasta) <- fasta.names
  #file.name <- strsplit(fastaFile, split=".", fixed=TRUE)[[1]][[1]]
  writeXStringSet(tmp.fasta, file=paste(fastaFile, "rename", sep=".", format="fasta"))
  return(tmp.fasta)
}


annotate.fasta.total <- function(fastaFile, taxonomy){
  GBseqs <- readDNAStringSet(fastaFile) #read .aln.rn
  
  ## Make a tmp file so don't overwrite 
  tmp.fasta <- DNAStringSet(GBseqs)
  fasta.names <- names(tmp.fasta)
  
  split <- strsplit(fasta.names, split="|", fixed=TRUE) #split names
  for (i in 1:length(split)){
    print(i)
    if (split[i][[1]][[1]] == "ncbi"){
      combinedname <- sapply(split[i], "[", 2L) 
      if (combinedname %in% taxonomy[,4]){ #Annotated_Tip
        # if name is ncbi
        tmpname = unique(as.character(taxonomy[which(combinedname==taxonomy[,4]),1]))[1] #PHLAWD
        fasta.names[i] =  as.character(paste(split[i][[1]][[1]], tmpname, split[i][[1]][[3]], split[i][[1]][[4]], split[i][[1]][[5]], split[i][[1]][[6]], sep="|"))
      }
    } else if (split[i][[1]][[1]] != "ncbi"){
      # if not ncbi
      #split2 <- strsplit(split[i][[1]], split="_", fixed=TRUE) 
      #combinedname <- paste(head(split2[[1]], (length(split2[[1]]) -1)), collapse ="_")
      if (names(tmp.fasta[i]) %in% taxonomy[,5]){ #MiSeq.Label.190117
        fasta.names[i] =  as.character(taxonomy[which(split[i]==taxonomy[,5]),3]) #MiSeq_Label_190402
        
      } 
    } else {
      fasta.names[i] =  fasta.names[i] 
    }
  }
  names(tmp.fasta) <- fasta.names
  writeXStringSet(tmp.fasta, file=paste(fastaFile, "rename", sep=".", format="fasta"))
  return(tmp.fasta)
  }
  