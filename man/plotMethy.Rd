\name{plotMethy}
\alias{plotMethy}
\title{
    Plot methylation profile
}
\description{
    Lollipop plot of the methylation profile.
}
\usage{
plotMethy(data, txdb, orgdb, range, gene, 
                      meanCutoff=c(totalC=0, methyC=0, ratio=0.05),
                      ...)
}

\arguments{
  \item{data}{
    A \link[GenomicRanges]{GRanges} or \link{RRBSeqDataSet} object.
}
  \item{txdb}{
    A \link[GenomicFeatures]{TxDb} object
  }
  \item{orgdb}{
    A \link[AnnotationDbi]{OrgDb-class} object
  }
  \item{range}{
    A \link[GenomicRanges]{GRanges} object. Giving the ranges to plot.
  }
  \item{gene}{
    Gene symbol. If range is not giving, gene will determine the ranges to plot.
  }
  \item{meanCutoff}{
  cutoff value for totalC, methyC and ratio.
  }
  \item{...}{
    Parameters could be passed to \link[trackViewer]{lolliplot}.
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
