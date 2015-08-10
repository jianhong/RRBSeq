### need to add chromosome location (maybe as piled) and gene start and end, and transcription start direction etc ########

plotMethylationProfile <- function(Methylation.summary, chr = "chr8", start = 74200195, end = 74214474, gene ="Jak3", upstream = 5000, downstream = 5000, pvalue.col, p.value.cutoff = 0.05, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =10, act.ratio.col = 13, height =8, width =11, cex=4, wt.color ="black", act.color = "red", significant.color = "purple", type=c("absolute", "relative", "lullipop"))
{
	maxup = max(0, (start-upstream))
	maxdown = end + downstream
	example= Methylation.summary[as.character(Methylation.summary[,chr.col ]) == chr & Methylation.summary[,CpG.pos.col] <= maxdown & Methylation.summary[,CpG.pos.col] >= maxup,]
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

plotMethylationProfile.noSig <- function(Methylation.summary, chr = "chr8", start = 74200195, end = 74214474, gene ="gene1", upstream = 5000, downstream = 5000,  chr.col = 6, CpG.pos.col = 7, wt.ratio.col =10, act.ratio.col = 13, height =8, width =11, cex=4, wt.color ="black", act.color = "red", significant.color = "purple", type=c("absolute", "relative", "lullipop"))
{
	example= Methylation.summary[Methylation.summary[,chr.col ] == chr & Methylation.summary[,CpG.pos.col ] <= (end + downstream) & Methylation.summary[,CpG.pos.col ] >= max(0,(start - upstream)),]
	example  = example[order(example[,CpG.pos.col ], decreasing = FALSE),]
	example = cbind(example, order(example[,CpG.pos.col ]))

	sig.CpGs  = example

	pdf(paste(gene,"pdf",sep="."),width=11, height=8)

	if( type == "relative")
	{
		plot(rbind(cbind(order(example[,CpG.pos.col]),100 * example[,wt.ratio.col]),cbind(order(example[,CpG.pos.col ]),100 * example[,act.ratio.col])) , col=wt.color ,xlab="CpGs", ylab="methylation level %", main=paste(gene, " Methylation Profile"), cex=cex )
		lines(order(example[,CpG.pos.col]), 100 * example[,act.ratio.col], col=act.color, type="o")
		lines(order(example[,CpG.pos.col]), 100 * example[,wt.ratio.col], col=wt.color, type="o")
		points(order(example[,CpG.pos.col]),100 * example[,act.ratio.col], col=act.color, cex=cex)
		#points(example[example[,CpG.pos.col] %in% sig.CpGs, dim(example)[2]], 100 * example[example[,CpG.pos.col ] %in% sig.CpGs, act.ratio.col], col=significant.color, cex=cex, pch="*")
	}
	if (type =="absolute")
	{
		plot(rbind(cbind(example[,CpG.pos.col],100 * example[,wt.ratio.col]),cbind(example[,CpG.pos.col ],100 * example[,act.ratio.col])) , col=wt.color ,xlab="CpGs", ylab="methylation level %", main=paste(gene, " Methylation Profile"), cex=cex )
		lines(example[,CpG.pos.col], 100 * example[,act.ratio.col], col=act.color, type="o")
		lines(example[,CpG.pos.col], 100 * example[,wt.ratio.col], col=wt.color, type="o")
		points(example[,CpG.pos.col],100 * example[,act.ratio.col], col=act.color, cex=cex)
		#points(sig.CpGs, 100 * example[example[,CpG.pos.col ] %in% sig.CpGs, act.ratio.col], col=significant.color, cex=cex, pch="*")
	}
	if (type =="lullipop")
	{		
		plot(rbind(cbind(order(example[,CpG.pos.col]),1),cbind(order(example[,CpG.pos.col]),2)) , ylim=c(0,2.5), col=wt.color ,xlab="CpGs", ylab="", main=paste(gene, " Methylation Profile"), cex=cex , axes=F)
		lines(cbind(order(example[,CpG.pos.col]), 2), col=act.color, type="o")
		lines(cbind(order(example[,CpG.pos.col]), 1), col=wt.color, type="o")
		n = dim(example)[1]
		points(order(example[,CpG.pos.col]),rep(2,n), col=act.color, cex=cex)
		#points(example[example[,CpG.pos.col] %in% sig.CpGs, dim(example)[2]], rep(2,length(sig.CpGs)), col=significant.color, cex=cex, pch="*")
		text(order(example[,CpG.pos.col]),1.2, 100 * example[,wt.ratio.col])
		text(order(example[,CpG.pos.col]),2.2, 100 * example[,act.ratio.col])
		text(3, 0.7, paste(paste(example[1,chr.col],example[1,CpG.pos.col], sep=":"), example[n, CpG.pos.col], sep="-"), cex = 1)
	}
	dev.off()
}

 r1 = read.table("~/Documents/ConsultingActivities/SolexaSeq/Methylation_seq/RRBSchrX-forSanchita.xls",sep="\t",header=TRUE)

plotMethylationProfile.noSig(r1, chr="chrX", start=103460373, end=103483233, upstream=0, gene="X:103460373-103483233", chr.col = 1, CpG.pos.col = 2, wt.ratio.col =8, act.ratio.col = 5, height =8, width =11, cex=4, wt.color ="black", act.color = "red",, type="relative")


plotMethylationProfile.noSig(r1, chr="chrX", start=74026821, end=74085636, upstream=0, downstream=0, gene="X:74026821-74085636", chr.col = 1, CpG.pos.col = 2, wt.ratio.col =8, act.ratio.col = 5, height =8, width =11, cex=4, wt.color ="black", act.color = "red",, type="relative")


#plotMethylationProfile.noSig(r1, chr="chrX", start=106187124, end=106203699, upstream=0, downstream=0, gene="X:106187124-106203699", chr.col = 1, CpG.pos.col = 2, wt.ratio.col =8, act.ratio.col = 5, height =8, width =11, cex=4, wt.color ="black", act.color = "red",, type="relative")

#plotMethylationProfile.noSig(r1, chr="chrX", start=52988078, end=53021660, upstream=0, downstream=0, gene="X:52988078-53021660",  chr.col = 1, CpG.pos.col = 2, wt.ratio.col =8, act.ratio.col = 5, height =8, width =11, cex=4, wt.color ="black", act.color = "red",, type="relative")

plotMethylationProfile.noSig(r1, chr="chrX", start=103431517, end=103484957, upstream=0, downstream=0, gene="X:103431517-103484957",  chr.col = 1, CpG.pos.col = 2, wt.ratio.col =8, act.ratio.col = 5, height =8, width =11, cex=4, wt.color ="black", act.color = "red",, type="relative")


plotMethylationProfile.noSig(r1, chr="chrX", start=103957167, end=103981284, upstream=0, downstream=0, gene="X:103957167-103981284",chr.col = 1, CpG.pos.col = 2, wt.ratio.col =8, act.ratio.col = 5, height =8, width =11, cex=4, wt.color ="black", act.color = "red",, type="relative")


############# plot with significance labeled #####


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

 #load('~/Documents/ConsultingActivities/SolexaSeq/Methylation_seq/Sept2012/RRBSmergedDataT0.RData')

r1 = read.table("/Users/zhuj/Documents/ConsultingActivities/SolexaSeq/Methylation_seq/Sept2012/T21-ACvsT21-WT-Allmin10C_stats.xls",sep="\t",header=TRUE)
 genes=read.csv("/Users/zhuj/Documents/ConsultingActivities/SolexaSeq/Methylation_seq/Sept2012/genes_plot.csv",sep=",",header=TRUE)
 genes
    Chr     Start       End   Gene
1  chr2 133377727 133379222   Bmp2
2 chr19  56471009  56472131  Casp7
3  chr4  88940111  88940848 Cdkn2a
4  chr7 150644823 150647374 Cdkn1c
5 chr15  79762842  79763232   Cbx7
6  chr7  35902799  35905760  Cebpa
7  chr4 151712524 151712625   Chd5
8 chr12   8304714   8306587   Gdf7
9 chr16  33829392  33830741  Itgb5

genes = genes[genes[,3] =="Individual",]

setwd("/Users/zhuj/Documents/ConsultingActivities/SolexaSeq/Methylation_seq/Sept2012/figures")

source("/Users/zhuj/Documents/ConsultingActivities/SolexaSeq/Methylation_seq/plotMethylationProfileFunction.R")

i =3
plotMethylationProfile(r1, chr = as.character(genes[i,1]), start = genes[i,2], end = genes[i,3], gene = as.character(genes[i,4]), upstream = 0, downstream = 0, pvalue.col= dim(r1)[2], p.value.cutoff = 0.05, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =43, act.ratio.col = 40, height =8, width =11, cex=4, wt.color ="green", act.color = "red", significant.color = "purple", type="absolute")

plotMethylationProfile(r1, chr = "chr3", start = 135354553, end = 135354945, gene = "Nfkb1", upstream = 0, downstream = 0, pvalue.col= dim(r1)[2], p.value.cutoff = 0.05, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =43, act.ratio.col = 40, height =8, width =11, cex=4, wt.color ="green", act.color = "red", significant.color = "purple", type="absolute")


for (i in 1:dim(genes)[1])
{
	plotMethylationProfile(r1, chr = as.character(genes[i,1]), start = genes[i,2], end = genes[i,3], gene = as.character(genes[i,4]), upstream = 0, downstream = 0, pvalue.col= dim(r1)[2], p.value.cutoff = 0.05, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =43, act.ratio.col = 40, height =8, width =11, cex=4, wt.color ="green", act.color = "red", significant.color = "purple", type="relative")
}


for (i in 1:dim(genes)[1])
{
	plotMethylationProfile(r1, chr = as.character(genes[i,1]), start = genes[i,2], end = genes[i,3], gene = paste(as.character(genes[i,4]),"T0",sep=""), upstream = 0, downstream = 0, pvalue.col= dim(r1)[2], p.value.cutoff = 0.05, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =25, act.ratio.col = 20, height =8, width =11, cex=4, wt.color ="green", act.color = "red", significant.color = "purple", type="relative")
}

for (i in 1:dim(genes)[1])
{
	plotMethylationProfile(r1, chr = as.character(genes[i,1]), start = genes[i,2], end = genes[i,3], gene = paste(as.character(genes[i,4]),"T14",sep=""), upstream = 0, downstream = 0, pvalue.col= dim(r1)[2], p.value.cutoff = 0.05, chr.col = 6, CpG.pos.col = 7, wt.ratio.col =35, act.ratio.col = 30, height =8, width =11, cex=4, wt.color ="green", act.color = "red", significant.color = "purple", type="relative")
}
