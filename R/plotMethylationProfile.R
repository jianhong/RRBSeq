### need to add chromosome location (maybe as piled) and gene start and end, and transcription start direction etc ########

plotMethylationProfile <- function(Methylation.summary, chr = "chr8", start = 74200195, end = 74214474, gene ="Jak3", upstream = 5000, downstream = 5000, pvalue.col, p.value.cutoff = 0.05, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =10, act.ratio.col = 13, height =8, width =11, cex=4, wt.color ="black", act.color = "red", significant.color = "purple", type=c("absolute", "relative", "lullipop"))
{
	example= Methylation.summary[Methylation.summary[,chr.col ] == chr & Methylation.summary[,CpG.pos.col ] <= (end + downstream) & Methylation.summary[,CpG.pos.col ] >= max(0,(start - upstream)),]
	example  = example[order(example[,CpG.pos.col ], decreasing = FALSE),]
	example = cbind(example, order(example[,CpG.pos.col ]))

	sig.CpGs  = example[example[,pvalue.col] < p.value.cutoff, CpG.pos.col]

	pdf(paste(gene,"pdf",sep="."),width=11, height=8)

	if( type == "relative")
	{
		plot(rbind(cbind(order(example[,CpG.pos.col]),100 * example[,wt.ratio.col]),cbind(order(example[,CpG.pos.col ]),100 * example[,act.ratio.col])) , col=wt.color ,xlab="CpGs", ylab="methylation level %", main=paste(gene, " Methylation Profile"), cex=cex )
		lines(order(example[,CpG.pos.col]), 100 * example[,act.ratio.col], col=act.color, type="o")
		lines(order(example[,CpG.pos.col]), 100 * example[,wt.ratio.col], col=wt.color, type="o")
		points(order(example[,CpG.pos.col]),100 * example[,act.ratio.col], col=act.color, cex=cex)
		points(example[example[,CpG.pos.col] %in% sig.CpGs, dim(example)[2]], 100 * example[example[,CpG.pos.col ] %in% sig.CpGs, act.ratio.col], col=significant.color, cex=cex, pch="*")
	}
	if (type =="absolute")
	{
		plot(rbind(cbind(example[,CpG.pos.col],100 * example[,wt.ratio.col]),cbind(example[,CpG.pos.col ],100 * example[,act.ratio.col])) , col=wt.color ,xlab="CpGs", ylab="methylation level %", main=paste(gene, " Methylation Profile"), cex=cex )
		lines(example[,CpG.pos.col], 100 * example[,act.ratio.col], col=act.color, type="o")
		lines(example[,CpG.pos.col], 100 * example[,wt.ratio.col], col=wt.color, type="o")
		points(example[,CpG.pos.col],100 * example[,act.ratio.col], col=act.color, cex=cex)
		points(sig.CpGs, 100 * example[example[,CpG.pos.col ] %in% sig.CpGs, act.ratio.col], col=significant.color, cex=cex, pch="*")
	}
	if (type =="lullipop")
	{		
		plot(rbind(cbind(order(example[,CpG.pos.col]),1),cbind(order(example[,CpG.pos.col]),2)) , ylim=c(0,2.5), col=wt.color ,xlab="CpGs", ylab="", main=paste(gene, " Methylation Profile"), cex=cex , axes=F)
		lines(cbind(order(example[,CpG.pos.col]), 2), col=act.color, type="o")
		lines(cbind(order(example[,CpG.pos.col]), 1), col=wt.color, type="o")
		n = dim(example)[1]
		points(order(example[,CpG.pos.col]),rep(2,n), col=act.color, cex=cex)
		points(example[example[,CpG.pos.col] %in% sig.CpGs, dim(example)[2]], rep(2,length(sig.CpGs)), col=significant.color, cex=cex, pch="*")
		text(order(example[,CpG.pos.col]),1.2, 100 * example[,wt.ratio.col])
		text(order(example[,CpG.pos.col]),2.2, 100 * example[,act.ratio.col])
		text(3, 0.7, paste(paste(example[1,chr.col],example[1,CpG.pos.col], sep=":"), example[n, CpG.pos.col], sep="-"), cex = 1)
	}
	dev.off()
}

plotMethylationProfile(r1, chr = "chr8", start = 74200195, end = 74214474, gene ="Jak3", upstream = 5000, downstream = 5000, pvalue.col= dim(r1)[2], p.value.cutoff = 0.05, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =10, act.ratio.col = 13, height =8, width =11, cex=4, wt.color ="black", act.color = "red", significant.color = "purple", type="relative")

plotMethylationProfile(r1, chr = "chr8", start = 74200195, end = 74214474, gene ="Jak3", upstream = 0, downstream = 0, pvalue.col= dim(r1)[2], p.value.cutoff = 0.05, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =10, act.ratio.col = 13, height =8, width =22, cex=4, wt.color ="black", act.color = "red", significant.color = "purple", type="absolute")

plotMethylationProfile(r1, chr = "chr8", start = 74200195, end = 74214474, gene ="Jak3", upstream = -1000, downstream = -2000, pvalue.col= dim(r1)[2], p.value.cutoff = 0.05, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =10, act.ratio.col = 13, height =8, width =26, cex=4, wt.color ="black", act.color = "red", significant.color = "purple", type="lullipop")

genes=read.table("/Users/zhuj/Documents/ConsultingActivities/SolexaSeq/Methylation_seq/bsmapRRBseq/LowerPriority.txt", sep="\t",header=TRUE)

genes=read.table("/Users/zhuj/Documents/ConsultingActivities/SolexaSeq/Methylation_seq/bsmapRRBseq/individualCpG/LowerPriority-individualCpG.txt", sep="\t", header=TRUE)

genes=read.table("/Users/zhuj/Documents/ConsultingActivities/SolexaSeq/Methylation_seq/bsmapRRBseq/individualCpG/HigherPriority-individualCpG.txt", sep="\t", header=TRUE)
genes = genes[genes[,3] =="Individual",]

setwd("/Users/zhuj/Documents/ConsultingActivities/SolexaSeq/Methylation_seq/bsmapRRBseq/individualCpG/")

##### r1 is for individual CpG

for (i in 1:dim(genes)[1])
{	
	p.value.cutoff = 0.05
	pvalue.col = dim(r1)[2]
	upstream = genes[i,7]
	gene = as.character(genes[i,1])
	plotMethylationProfile(r1, chr = paste("chr",genes[i,4],sep=""), start = genes[i,5], end = genes[i,6], gene = gene, upstream = genes[i,7], downstream = genes[i,8], pvalue.col= pvalue.col, p.value.cutoff = p.value.cutoff, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =10, act.ratio.col = 13, height =8, width =11, cex=4, wt.color ="black", act.color = "red", significant.color = "purple", type="relative")
}


##### CpGisland.details is for CpG island

{if(CpGisland)
	{
	
	example = CpGisland.details[!is.na(CpGisland.details$symbol) & CpGisland.details$symbol == gene ,]
	example  = example[order(example[,28], decreasing = FALSE),]
	if (dim(example)[1] >0)
	{
		pdf(paste(gene,"pdf",sep="."),width=11, height=8)
		plot(rbind(cbind(order(example[,28]),100 * example[,37]),cbind(order(example[,28]),100 * example[,40])) , col="black" ,xlab="CpGs", ylab="methylation level %", main=paste(gene, " Methylation Profile"), cex=4 )
		lines(order(example[,28]), 100 * example[,40], col="red", type="o")
		lines(order(example[,28]), 100 * example[,37], col="black", type="o")
		points(order(example[,28]),100 * example[,40], col="red", cex=4)
		dev.off()
	}
	else
	{
		print(paste("this gene not near any CpG island: ", gene))
	}
	}
}

