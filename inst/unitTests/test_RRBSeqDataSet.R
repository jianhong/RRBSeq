test_RRBSeqDataSet <- function(){
    path <- system.file("extdata", package="RRBSeq")
    methyRatioFiles <- dir(path, "^test", full.names=TRUE)
    rrb <- RRBSeqDataSet(methyRatioFiles, ind.totalC=5, ind.methyC=4, ind.ratio=6, header=TRUE)
    gr <- as(rrb, "GRanges")
    rrb.resampled <- summarize.sliding.window(rrb)
    checkEqualsNumeric(methyC(rrb.resampled), 
                       c(18, 17, 25, 25, 22, 9, 1, 1, 1, 1, 3,
                         18, 17, 25, 25, 22, 9, 1, 1, 1, 1, 3))
    checkEqualsNumeric(totalC(rrb.resampled), 
                       c(64, 70, 78, 78, 66, 16, 10, 10, 10, 10, 30,
                         64, 70, 78, 78, 66, 16, 10, 10, 10, 10, 30))
}