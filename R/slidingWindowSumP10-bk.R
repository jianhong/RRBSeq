library(GenomicRanges)
#methySummary = read.table('~/scratch/RRBseq/P10WTSingleCpGall.xls', sep="\t", header=TRUE)
methySummary = read.table('~/Documents/hiseq/MichaelGreen/IchiroOnoyama/RRBseq/singleCpG/P10WTSingleCpGall.xls', sep="\t", header=TRUE)
dim(methySummary)
colnames(methySummary)

sliding.Window.Sum <- function(Ginput, window.size =100, step=20)
{
	steps.per.window = window.size / step - 1
	binned = data.frame()
	for (i in 0:steps.per.window)
	{
		  window.start = i * step + 1
		binned = rbind(binned,do.call(rbind, lapply(unique(as.character(seqnames(Ginput))), function(chr) 
		{
			temp = subset(Ginput, as.character(seqnames(Ginput))==chr & start(Ginput) >= window.start)
			if (length(temp) >0)
			{
				last.breaks = ceiling((max(start(temp)) - window.start)/window.size) * window.size + window.start
				breaks = seq(from=window.start, to = last.breaks , by=window.size)
				#print(last.breaks)
				#print(breaks)
				#print(dim(temp))
				if ( length(breaks) > 1)
				{
					hist.methyC = hist(rep(start(temp), temp$methyC), breaks= breaks, plot=FALSE)
					hist.totalC = hist(rep(start(temp), temp$totalC), breaks= breaks, plot=FALSE)
					cbind(rep(chr, length(hist.methyC$mids)), hist.methyC$mids, hist.methyC$counts, hist.totalC$counts)
				}
				else
					c(chr,breaks, sum(temp$methyC), sum(temp$totalC))
			}
		})))
	} ### end for loop
	temp = binned[order(as.numeric(as.character(binned[,2]))),]
	gr1 = GRanges(IRanges(start=as.numeric(as.character(temp[,2])) - window.size/2, width= rep(window.size,dim(temp)[1])),    seqnames=temp[,1],  
methyC= as.numeric(as.character( temp[,3])),	
 totalC= as.numeric(as.character( temp[,4])))
	gr1 = gr1[mcols(gr1)$totalC >0,]
	mcols(gr1)$ratio = mcols(gr1)$methyC / mcols(gr1)$totalC
	gr1
}
	
####### Easier to directly get data from the RRBS individual output
######## read.methySeq(file, format="RRBS")

P10.T0.1 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 13])),
totalC =  as.numeric(as.character( methySummary[, 12])))

P10.T0.2 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 10])),
totalC =  as.numeric(as.character( methySummary[, 9])))

P10.T0.3 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 7])),
totalC =  as.numeric(as.character( methySummary[, 6])))

P10.T14.1 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 22])),
totalC =  as.numeric(as.character( methySummary[, 21])))

P10.T14.2 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 19])),
totalC =  as.numeric(as.character( methySummary[, 18])))

P10.T14.3 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 16])),
totalC =  as.numeric(as.character( methySummary[, 15])))


P10.T21.1 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 31])),
totalC =  as.numeric(as.character( methySummary[, 30])))

P10.T21.2 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 28])),
totalC =  as.numeric(as.character( methySummary[, 27])))

P10.T21.3 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 25])),
totalC =  as.numeric(as.character( methySummary[, 24])))

wt.T0.1 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 40])),
totalC =  as.numeric(as.character( methySummary[, 39])))

wt.T0.2 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 37])),
totalC =  as.numeric(as.character( methySummary[, 36])))

wt.T0.3 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 34])),
totalC =  as.numeric(as.character( methySummary[, 33])))

wt.T14.1 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 49])),
totalC =  as.numeric(as.character( methySummary[, 48])))

wt.T14.2 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 46])),
totalC =  as.numeric(as.character( methySummary[, 45])))

wt.T14.3 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 43])),
totalC =  as.numeric(as.character( methySummary[, 42])))


wt.T21.1 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 58])),
totalC =  as.numeric(as.character( methySummary[, 57])))

wt.T21.2 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 55])),
totalC =  as.numeric(as.character( methySummary[, 54])))

wt.T21.3 = GRanges(IRanges(start=as.numeric(as.character(methySummary[,2])), width=rep(1, dim(methySummary)[1])), seqnames=methySummary[,1],
methyC =  as.numeric(as.character( methySummary[, 52])),
totalC =  as.numeric(as.character( methySummary[, 51])))


window.size = 200
step = 50

w200.s50.P10.T21.1 = sliding.Window.Sum(P10.T21.1, window.size =window.size, step=step)

w200.s50.P10.T21.2 = sliding.Window.Sum(P10.T21.2, window.size =window.size, step=step)

w200.s50.P10.T21.3 = sliding.Window.Sum(P10.T21.3, window.size =window.size, step=step)

w200.s50.P10.T14.1 = sliding.Window.Sum(P10.T14.1, window.size =window.size, step=step)

w200.s50.P10.T14.2 = sliding.Window.Sum(P10.T14.2, window.size =window.size, step=step)

w200.s50.P10.T14.3 = sliding.Window.Sum(P10.T14.3, window.size =window.size, step=step)

w200.s50.P10.T0.1 = sliding.Window.Sum(P10.T0.1, window.size =window.size, step=step)

w200.s50.P10.T0.2 = sliding.Window.Sum(P10.T0.2, window.size =window.size, step=step)

w200.s50.P10.T0.3 = sliding.Window.Sum(P10.T0.3, window.size =window.size, step=step)

save.image("P10window200step50.RData")

w200.s50.wt.T21.1 = sliding.Window.Sum(wt.T21.1, window.size =window.size, step=step)

w200.s50.wt.T21.2 = sliding.Window.Sum(wt.T21.2, window.size =window.size, step=step)
####
w200.s50.wt.T21.3 = sliding.Window.Sum(wt.T21.3, window.size =window.size, step=step)

w200.s50.wt.T14.1 = sliding.Window.Sum(wt.T14.1, window.size =window.size, step=step)

w200.s50.wt.T14.2 = sliding.Window.Sum(wt.T14.2, window.size =window.size, step=step)

w200.s50.wt.T14.3 = sliding.Window.Sum(wt.T14.3, window.size =window.size, step=step)

w200.s50.wt.T0.1 = sliding.Window.Sum(wt.T0.1, window.size =window.size, step=step)

w200.s50.wt.T0.2 = sliding.Window.Sum(wt.T0.2, window.size =window.size, step=step)

w200.s50.wt.T0.3 = sliding.Window.Sum(wt.T0.3, window.size =window.size, step=step)

save.image("P10WTwindow200step50.RData")
