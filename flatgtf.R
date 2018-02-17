# This script is used to take a gtf file and "flatten" it.
# Sometimes overlapping annotations are problematic, this merges all overlapping exons.  

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
library(org.Hs.eg.db)


txdb<-makeTxDbFromGFF("refseq.hg19.gtf", format="gtf")
exonsByGene <- exonsBy(txdb,'gene')
mergedGene <- reduce(exonsByGene)
mergedGene <-
## Remove a few hundred 'tricky cases' like:
##  * trans-spliced genes (?? are there any in hg19.knownGene ??)
##  * genes on multiple chromosomes (eg: alt. haplotypes "chr_ctg9_hap1" )
##    Arguably(?) these should not be considered the 'same gene'.
mergedGene[1==elementLengths(runValue(strand(mergedGene))) &
                    1==elementLengths(runValue(seqnames(mergedGene)))]


mapped_seqs <- mappedkeys(org.Hs.egREFSEQ2EG)
xx <- as.list(org.Hs.egREFSEQ2EG[mapped_seqs])

refseqIDs<-names(mergedGene)
refseqIDsEG<-xx[refseqIDs]
refseqIDsEGdf<-do.call(rbind.data.frame, refseqIDsEG)
colnames(refseqIDsEGdf)<-"Entrez"
refseqIDsEGdf$refseq<-rownames(refseqIDsEGdf)

mapped_seqs <- mappedkeys(org.Hs.egSYMBOL)
xx <- as.list(org.Hs.egSYMBOL[mapped_seqs])

symbols <-xx[as.character(refseqIDsEGdf$Entrez)]
symbolsdf<-do.call(rbind.data.frame, symbols)
colnames(symbolsdf)<-"Gene"
symbolsdf$Entrez<-rownames(symbolsdf)

export(mergedGene,'refseq.hg19.merged.gff3')
gff<-read.table("refseq.hg19.merged.gff3")
colnames(gff)[9:10]<-c("ID", "refseq")
gff$refseq<-gsub("Name=", "", gff$refseq)
gff$refseq<-gsub(";", "", gff$refseq)

temp1<-merge(gff, refseqIDsEGdf, by="refseq", all.x=T)
full_data<-merge(temp1, symbolsdf, all.x=T, by="Entrez")
full_data<-(full_data[,c(3:10,11,2,1,12)])
write.table(full_data, "refseq.hg19.merged.gff3",quote=F, sep="\t", row.names=F)







save.image()
savehistory()

