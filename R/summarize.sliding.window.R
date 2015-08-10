summarize.sliding.window <- function(gr, window.size =200L, step=50L, summarize=c("mean.ratio", "total.methyC.over.total.totalC"))
{
	cat(date(), "Validating input ...\n");
	 summarize = match.arg(summarize)
	if (missing(gr)) {
        stop("Missing required argument gr!")
    }
    if (class(gr) != "GRanges" ) {
        stop("No valid gr passed in. It needs to be GRanges object")
    }
    if ( class(window.size) != "integer" || window.size < 2) {
        stop("window.size needs to be an integer greater than 1")
    }
	if ( class(step) != "integer" || step > window.size) {
        stop("step needs to be an integer less or equal to window")
    }

	steps.per.window = window.size / step - 1
	binned = data.frame()
	
	n.samples = dim(mcols(gr)[, grep("totalC", colnames(mcols(gr)))])[2]
	
	if (mode(n.samples) == "NULL") { n.samples = 1}
	for (i in 0:steps.per.window)
	{
		  window.start = i * step + 1
		binned = rbind(binned,do.call(rbind, lapply(unique(as.character(seqnames(gr))), function(chr) 
		{
			cat(date(), "summarizing sliding window for chromosome", chr, " ...\n");
			temp = subset(gr, as.character(seqnames(gr))==chr & start(gr) >= window.start)
			if (length(temp) >0)
			{
				last.breaks = ceiling((max(start(temp)) - window.start)/window.size) * window.size + window.start
				breaks = seq(from=window.start, to = last.breaks , by=window.size)
				if  (length(breaks) > 1)
				{
					methyC = matrix()
					totalC = matrix()
					ratio = matrix()
					this.methyC = mcols(temp)[, grep("methyC", colnames(mcols(temp)))]
					this.totalC = mcols(temp)[, grep("totalC", colnames(mcols(temp)))]
					this.ratio = mcols(temp)[, grep("ratio", colnames(mcols(temp)))]
					
					for (k in 1:n.samples)
					{
						if (n.samples > 1)
						{
						hist.methyC = hist(rep(start(temp), this.methyC[,k]), breaks= breaks, plot=FALSE)
						hist.totalC = hist(rep(start(temp), this.totalC[,k]), breaks= breaks, plot=FALSE)
						}
						else
						{
						hist.methyC = hist(rep(start(temp), this.methyC), breaks= breaks, plot=FALSE)
						hist.totalC = hist(rep(start(temp), this.totalC), breaks= breaks, plot=FALSE)
						}
						if (k >1)
						{
							methyC = cbind(methyC, hist.methyC$counts)
							totalC = cbind(totalC, hist.totalC$counts)
						}
						else
						{
							methyC = hist.methyC$counts
							totalC = hist.totalC$counts
						}
						if (summarize == "mean.ratio")
						{
							if (n.samples > 1)
							{
							hist.ratio = hist(rep(start(temp), round(ceiling(this.ratio[,k] * 10000), 0)), breaks= breaks, plot=FALSE)
							}
							else
							{
							hist.ratio = hist(rep(start(temp), round(ceiling(this.ratio * 10000), 0)), breaks= breaks, plot=FALSE)
							}
							hist.nC = hist(start(temp), breaks= breaks, plot=FALSE)
							if (k >1)
								ratio = cbind(ratio, hist.ratio$counts/10000/hist.nC$counts)
							else
								ratio = hist.ratio$counts/10000/hist.nC$counts
						}
						else if (summarize == "total.methyC.over.total.totalC")
						{
							if (k >1)
								ratio = cbind(ratio, hist.methyC$counts / hist.totalC$counts)
							else
								ratio = hist.methyC$counts / hist.totalC$counts
						}
					}  ### end for samples
					cbind(rep(chr, length(hist.methyC$mids)), hist.methyC$mids, methyC, totalC, ratio)
				} #### end if breaks > 1
			} ### end if length(temp) >0
			})))
	} ### end for loop windows

	temp = binned[order(as.numeric(as.character(binned[,2]))),]
	last.methyC = n.samples + 2
	first.totalC = last.methyC + 1
	last.totalC = last.methyC + n.samples
	first.ratio = last.totalC + 1
	last.ratio = last.totalC + n.samples

	GRanges(IRanges(start=as.numeric(as.character(temp[,2])) - window.size/2, width= rep(window.size,dim(temp)[1])),    seqnames=temp[,1],  methyC= temp[,3:last.methyC], totalC= temp[,first.totalC:last.totalC], ratio= temp[,first.ratio:last.ratio])
}
