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

