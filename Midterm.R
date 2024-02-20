#necessary packages needed
library(msa)
library(Biostrings)
library(seqinr)
library(phangorn)
library(tidyr)
library(dplyr)

#this is me labeling and reading the file into R
mt <- readDNAStringSet("sequences.fasta")
mt

#this is me aligning the sequence file
mtaline <-msa(mt,"ClustalW")

#the alinement will now be printed in the console
print(mtaline,show="complete")

#Now I saved my work from Rstudio as a file to my Bioinofrmatics folder on my computure. 
#I uploded the file to BLAST 
#I ran the file on BLAST and discovered a patient is experiencing a hbb gene for beta-globin
#There is a single point mutation on CD8 of (AAG>AAC)

#I am now extraxting signle sequences from the 20 sequences given to us from the file.
#This translation is nesseczsry to sequence to a protein.

#for the first sequence:
class(mt)
homo1 <-mt$Homo_sapiens_1
print(homo1)

#for the second sequence:
homo2 <-mt$Homo_sapiens_2
print(homo2)

#for the third sequence:
homo3<- mt$Homo_sapiens_3
print(homo3)

#for fourth sequence:
homo4 <- mt$Homo_sapiens_4
print(homo4)

#for fith sequence:
homo5<- mt$Homo_sapiens_5
print(homo5)

#for sixth sequence:
homo6 <-mt$Homo_sapiens_6
print(homo6)

#for seventh sequence:
homo7 <-mt$Homo_sapiens_7
print(homo7)

#for eigth sequence:
homo8<- mt$Homo_sapiens_8
print(homo8)

#for ninth sequence:
homo9<-mt$Homo_sapiens_9
print(homo9)

#for thenth sequence:
homo10<-mt$Homo_sapiens_10
print(homo10)

#for eleventh sequence:
homo11<-mt$Homo_sapiens_11
print(homo11)

#for twelth sequence:
homo12<- mt$Homo_sapiens_12
print(homo12)

#I am now writing a code to translate my original alinement sequence of homo6 into an amino acid sequence 
actualAA<- Biostrings::translate(homo6, genetic.code = GENETIC_CODE, no.init.codon = FALSE, if.fuzzy.codon = "error")

#The following will print my amino acid sequence into my console
print(actualAA)

#I will now be saving my sequence for homo6 as a fasta file by writing the following code.
for (i in names(actualAA))
  +  writeXStringSet(myAAStringSet[[i]], filepath = paste("/Users/phoebeburns/Documents/GitHub/Bioinformatics",i,".fasta"), format = "fasta")
write.fasta(names ="actualAA", sequences=actualAA, file.out = "actualAA.fasta")  

#I ran my fasta file through UniProt and recived the accesion number of A0A0J9YWK4
#There we multiple diseases that coulkd be associated with this gene mutation. I belive the most likely diagnosis is Beta Thalassemia
#I will now be uploading my protien image and code into GitHub. 
