if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("devtools")
library(devtools)

devtools::install_github("pievos101/PopGenome")
library(PopGenome)
mysequence <- system.file("sequence.bluewhale.fasta","/Users/phoebeburns/Documents/GitHub/Bioinformatics/sequence.bluewhale1.fasta",package = "msa")
mysequence1 <- readDNAStringSet("/Users/phoebeburns/Documents/GitHub/Bioinformatics/sequence.bluewhale1.fasta",format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE, with.qualities = FALSE)
mysequence1

mysequence2 <- system.file("sequence.bluewhale2.fasta","/Users/phoebeburns/Documents/GitHub/Bioinformatics/sequence.bluewhale2.fasta", package = "msa")
mysequence2 <- readDNAStringSet("/Users/phoebeburns/Documents/GitHub/Bioinformatics/sequence.bluewhale2.fasta", format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE, with.qualities = FALSE)
mysequence2

mysequence3 <- system.file("sequence.bluehwale3.fasta","/Users/phoebeburns/Documents/GitHub/Bioinformatics/sequence.bluewhale3.fasta",package = "msa")
mysequence3 <-readDNAStringSet("/Users/phoebeburns/Documents/GitHub/Bioinformatics/sequence.bluewhale3.fasta", format = "fasta",nrec = -1L,skip = 0L, seek.first.rec = FALSE,use.names = TRUE, with.qualities = FALSE)
mysequence3

mysequences4 <- system.file("sequence.bluewhale4.fasta","/Users/phoebeburns/Documents/GitHub/Bioinformatics/sequence.bluewhale4.fasta", package = "msa")
mysequences4 <- readDNAStringSet("/Users/phoebeburns/Documents/GitHub/Bioinformatics/sequence.bluewhale4.fasta", format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE, with.qualities = FALSE)
mysequences4

mysequence5 <- system.file("sequence.bluewhale5.fasta","/Users/phoebeburns/Documents/GitHub/Bioinformatics/sequence.bluewhale5.fasta", package = "msa")
mysequence5 <- readDNAStringSet("/Users/phoebeburns/Documents/GitHub/Bioinformatics/sequence.bluewhale5.fasta",format = "fasta", nrec = -1L,skip = 0L, seek.first.rec = FALSE, use.names = TRUE, with.qualities = FALSE)
mysequence5

allsequences <- c(mysequence1,mysequence2,mysequence3,mysequences4,mysequence5)
allsequences

aline <- msa(allsequences,"ClustalW")
aline

print(aline, show="complete") 

nchar(allsequences)

letterFrequency(allsequences[[1]], letters = "CG", OR= 0)
letterFrequency(allsequences[[2]], letters = "CG", OR= 0)
letterFrequency(allsequences[[3]], letters = "CG", OR= 0)
letterFrequency(allsequences[[4]], letters = "CG", OR= 0)
letterFrequency(allsequences[[5]], letters = "CG", OR= 0)

aline2<-msaConvert(aline, type="seqinr::alignment")
d<-dist.alignment(aline2, matrix="identity")
d

library("Biostrings")
myseq <- readDNAStringSet("sequence.bluewhale1.fasta", format = "fasta")
head(myseq)


#mysequence5_AA <- Biostrings::translate(mysequence5, genetic.code=GENETIC_CODE, no.init.codon=FALSE,if.fuzzy.codon="error")

for(n in 1:15){translate(myseq, (genetic.code = GENETIC_CODE))}

Alignment_phyDat <- msaConvert(aline,type="phangorn::phyDat")
write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")

#Second accession
library(UniprotR)
library(protti)
library(r3dmol)
BiocManager::install("GenomicAlignments")


