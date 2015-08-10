compareMethylation <-
function(gr, gr.names,Treatment, norm.method=c("none", "quantiles", "quantiles.robust", "mean", "median"),design, ...)
{
	if(missing(gr))
	{
		stop("gr is required!")
	}
	if (class(gr) != "GRanges" ) {
        stop("No valid gr passed in. It needs to be GRanges object")
    }

	if (length(design) ==0)
	{
		stop("design is required, please refer to limma package for create design ")
	}
	
	norm.d  = normalizeData(gr, method=norm.method)

	exprs = norm.d$exprs
	colnames(exprs)=paste("normalized.ratio", gr.names, sep=".")	
        	methyC = norm.d$methyC
	colnames(methyC)=paste("methyC", gr.names, sep=".")
       	totalC = norm.d$totalC
	colnames(totalC)=paste("totalC", gr.names, sep=".")
	fr = cbind(rownames(methyC), exprs, methyC, totalC)
	colnames(fr)[1]= "ID"

	eSet = ExpressionSet(assayData=exprs)	
	fit1 = lmFit(eSet, design)
	fit2 <- eBayes(fit1)
        res = topTable(fit2, num = dim(exprs)[1], coef="typewt")
        res = res[,c(1:2, 5:6)]
        colnames(res)[2:4] = paste(Treatment, colnames(res)[2:4],sep=".")
	colnames(res) = gsub("logFC", "difference", colnames(res))
        fr = merge(fr, res,by = "ID")
}
