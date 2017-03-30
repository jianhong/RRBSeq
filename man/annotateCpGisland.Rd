\name{annotateCpGisland}
\alias{annotateCpGisland}
\title{
    Annotate CpGs islands
}
\description{
    Add CpGs island information to data.
}
\usage{
annotateCpGisland(data, CpGs, ...)
}

\arguments{
  \item{data}{
    A \link[GenomicRanges]{GRanges} object.
}
  \item{CpGs}{
    A \link[GenomicRanges]{GRanges} object.
}
  \item{\dots}{
    Not used.
}
}
\value{
A GRanges object.
}
\author{
Jianhong Ou
}
\examples{
if(interactive()){
extdata <- system.file("extdata", package="RRBSeq")
mrf <- dir(extdata, "_T.*\\\\.txt", full.names=TRUE)
rrb <- RRBSeqDataSet(mrf)

##  transform data
rrb <- transformData(rrb, transformation="asin")

## summarize.sliding.window before compare
rrb <- summarize.sliding.window(rrb, window.size=200L, step=50L, summarize="mean.ratio")

## normalize the data
rrb <- normalizeData(rrb, norm.method="quantiles")

# do comparison by default
g <- factor(gsub("inst.extdata..(.*?)_.*.txt", "\\\\1", mrf))
design <- model.matrix(~-1+g)
design(rrb) <- design
design(rrb)

const.mat <- makeContrasts(contrasts=c("gBraf-gwt"), levels=design)
contrasts.matrix(rrb) <- const.mat
contrasts.matrix(rrb)

res <- compareMethylation(rrb, coef=1, robust=FALSE)

res$CpGislandScore <- CpGislandScore(res, genome=Mmusculus)
res.s <- res[res$CpGislandScore>=.65 & res$P.Value<0.05]

CpGs <- queryCpGFromUCSC("Mus musculus")
mcols(CpGs)
mm9.CpGs <- CpGs[[1]]
mm10.CpGs <- liftOverCpGs(mm9.CpGs, from="mm9", to="mm10")
mm10.CpGs <- unlist(mm10.CpGs)

res.s <- annotateCpGisland(res.s, mm10.CpGs)

library(EnsDb.Mmusculus.v79)
annotation <- genes(EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- seqlevelsStyle(res.s)
res.s <- annotateGeneInfo(res.s, annotation)

## plot the results
symbol <- res.s$symbol[1]
plotMethy(res.s, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db, gene=symbol)
plotMethy(res.s, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db, gene=symbol, type="circle")
plotMethy(raw, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db, gene=symbol)
}
}
