
normalizeData <- function(gr, method=c("none", "quantiles", "quantiles.robust", "mean", "median"),  ... )
{
	exprs = mcols(gr)[, grep("ratio", colnames(mcols(gr)))]
	methyC = mcols(gr)[, grep("methyC", colnames(mcols(gr)))]
	totalC = mcols(gr)[, grep("totalC", colnames(mcols(gr)))]
		
	rownames(methyC) = paste(as.character(seqnames(gr[[1]])), as.character(start(gr[[1]])), as.character(end(gr[[1]])), sep="-")
	rownames(totalC) = paste(as.character(seqnames(gr[[1]])), as.character(start(gr[[1]])), as.character(end(gr[[1]])), sep="-")
	if (method == "quantile")
	{
		exprs  = normalize.quantiles(exprs)
	}
	else if (method == "quantiles.robust")
	{
		exprs  = normalize.quantiles.robust(exprs, ...)
	}
	else if (method == "mean")
	{
		avgs = colSummarizeAvg(exprs)$Estimates
		scaling.factors = avgs /avgs[1]
		scaling.factors = matrix(rep(scaling.factors, dim(exprs)[1]), ncol = length(scaling.factors), byrow=TRUE)
		exprs = exprs /scaling.factors
	}
	else if (method == "median")
	{
		avgs= colSummarizeMedian(exprs)$Estimates
		scaling.factors = avgs /avgs[1]
		scaling.factors = matrix(rep(scaling.factors, dim(exprs)[1]), ncol = length(scaling.factors), byrow=TRUE)
		exprs = exprs /scaling.factors
	}

	rownames(exprs) = paste(as.character(seqnames(gr[[1]])), as.character(start(gr[[1]])), as.character(end(gr[[1]])), sep="-")
	list(exprs = exprs,  methyC=methyC, totalC = totalC)
}
