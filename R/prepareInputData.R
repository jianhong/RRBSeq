prepareInputData <- function(methyRatioFiles, 
                             ind.chr=1L, ind.pos=2L,  
                             ind.strand=3L, ind.totalC=6L, 
                             ind.methyC=7L, ind.ratio=5L, 
                             min.totalC=5L, ...)
{
    if (min.totalC < 1){
        warning("min.totalC needs to be at least 1, so automatically set to 1")
        min.totalC = 1
    }
    allInputData = list()
    for (i in 1:length(methyRatioFiles))
    {
        x = read.table(methyRatioFiles[i], ...)
        allInputData[[i]] = x[, c(ind.chr, ind.pos, ind.strand, ind.totalC, ind.methyC, ind.ratio)]
    }
    mergedData = allInputData[[1]]
    for (i in 2:length(methyRatioFiles))
    {
        mergedData = merge(mergedData, allInputData[[i]], ...)
    }
    ind.methyC = numeric()
    ind.totalC = numeric()
    ind.ratio = numeric()
    for (i in 1:length(methyRatioFiles))
    {
        ind1 = 3* i +  1
        ind2 = ind1 + 1
        ind3 = ind2 + 1
        ind.totalC = c(ind.totalC, ind1)
        ind.methyC = c(ind.methyC, ind2)
        ind.ratio = c(ind.ratio, ind3)
        mergedData = subset(mergedData, as.numeric(as.character(mergedData[[,ind1]])) >= min.totalC)
    }
    GR  = GRanges(IRanges(start=as.numeric(as.character(mergedData[ , ind.pos])), width= rep(1, dim(temp)[1])),  seqnames=mergedData[, ind.chr],  methyC= mergedData[, ind.methyC], totalC=  mergedData[ , ind.totalC], ratio = mergedData[, ind.ratio])
    GR
}
