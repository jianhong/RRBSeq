transformData <- function(gr, transformation="asin")
{
	samples = grep("ratio", colnames(mcols(gr)))
	totalC.ind = grep("totalC", colnames(mcols(gr)))
	#mcols(gr)$ratio[mcols(gr)$ratio==0] = 1 / (4 *mcols(gr)$totalC[mcols(gr)$ratio==0,)
	#mcols(gr)$ratio[mcols(gr)$ratio]==1] = 1 - 1 / (4 *mcols(gr)$totalC[mcols(gr)$ratio==1,])
	for (i in 1:length(samples))
	{
		mcols(gr)[mcols(gr)[,samples[i]]==0,samples[i]] = 1 / (4 *mcols(gr)[mcols(gr)[,samples[i]]==0,totalC.ind[i]])
		mcols(gr)[mcols(gr)[,samples[i]]==1, samples[i]] = 1 - 1 / (4 *mcols(gr)[mcols(gr)[,samples[i]]==1,totalC.ind[i]])
	}
	if (transformation == "asin")
	{
		#mcols(gr)$ratio = asin(sqrt(mcols(gr)$ratio))
		for (i in samples)
			mcols(gr)[,i] = asin(sqrt(mcols(gr)[,i]))
	}
	 if (transformation == "log2")
       	 {
		for (i in samples)
			mcols(gr)[,i] = log2(mcols(gr)[,i])
		#mcols(gr)$ratio = log2(mcols(gr)$ratio)

	}
	gr
}
