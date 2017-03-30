CpGislandScore <- function(data, genome, ...){
    stopifnot(class(data)=="GRanges")
    stopifnot(length(data)>0)
    stopifnot(class(genome)=="BSgenome")
    len <- width(data)
    seq <- getSeq(genome, data)
    gcContent <- letterFrequency(seq, letters="CG", OR=0)
    expect <- gcContent[, "C"] * gcContent[, "G"] / len
    CpGs <- vcountPattern("CG", seq, max.mismatch=0, min.mismatch=0, fixed=TRUE)
    observed2expectedCpGratio <- CpGs/expect
    observed2expectedCpGratio
}