library(UniprotR)
library(protti)
library(GenomicAlignments)
library(r3dmol)

a1<- read.csv("/Users/phoebeburns/Documents/GitHub/Bioinformatics/Accession_1_P0A799.csv")
a1
ls(a1

class(a1)

list<-as.list(df)
list

class(list)
class(list[["Accessions"]])

PGI_obj2<-GetProteinG0Info(List_2[["Accessions"]],directorypath = NULL)
PlotGoInfo(PGI_obj2,directorypath = NULL)
PlotGOALL(G00bj = PGI_obj2, Top = 10, directorypath = getwd(), width = 8, height -
Patho2 <-GetPathology_Biotech(List[["Accessions"]] ,directorypath=NULL)
Get. diseases(Patho2, directory=NULL)       

uni_protInfo2<-fetch_uniprot("P08839")
pdb_2<- fetch_pdb （pdb_ids=c（"2HWG"））
head(pdb_2)
