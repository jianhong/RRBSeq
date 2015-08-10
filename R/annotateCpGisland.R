setwd("~/scratch/RRBseq/Sep2012/bsmap/")
library(ChIPpeakAnno)
data(TSS.mouse.NCBIM37)
comparison = c("STC1-ACvsWT", "T0-ACvsWT", "T14-ACvsWT", "T21-ACvsWT")
chr.column = 43
start.column = 45
end.column =46
library(org.Mm.eg.db)
for (i in 1:length(comparison))
{
	r2 = read.table(file=paste(comparison[i], "Sig_min10CpGisland_adjPlt0.05.xls",sep="-"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
	r3 = unique(cbind(r2[,chr.column], r2[,start.column], r2[,end.column]))
	r3.RD = RangedData(IRanges(start=as.numeric(r3[,2]), end=as.numeric(r3[,3])), space=r3[,1])
	r3.ann = annotatePeakInBatch(r3.RD, AnnotationData = TSS.mouse.NCBIM37, maxgap=5000, select="all", output="both")

	colnames(r2)[start.column] = "start"
	r4 = addGeneIDs(r3.ann,"org.Mm.eg.db",c("symbol","genename"))

	temp = as.data.frame(r4)
	colnames(temp)[1] = "chr"
	temp[,1] = paste("chr", as.character(temp[,1]),sep="")
	colnames(temp)
	colnames(r2)
	colnames(temp)[7] = "gene_strand"
	### get rid of the name of individual sequenced C
	r2 = unique(r2[,-44])
	
	changed.CpG = merge(r2, temp)

	write.table(changed.CpG , file=paste(comparison[i],"changed_min10CpGisland_adjPlt0.05_gene.xls",sep="-"), sep="\t", row.names=FALSE)
}
