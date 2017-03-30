\name{annotateGeneInfo}
\alias{annotateGeneInfo}
\title{
    Annotate CpGs islands
}
\description{
    Add CpGs island information to data.
}
\usage{
annotateGeneInfo(data, annotation, 
                  bindingType="nearestBiDirectionalPromoters",
                  bindingRegion=c(-5000, 500),
                  ignore.data.strand=TRUE, ...)
}

\arguments{
  \item{data}{
    A \link[GenomicRanges]{GRanges} object.
}
  \item{annotation}{
    annotation data, \link[GenomicRanges]{GRanges} object
  }
  \item{bindingType}{
    Specifying the criteria to associate methylation sites with annotation. Here is how to use it together with the parameter bindingRegion
    \itemize{
        \item To obtain methylation sites within 5kb upstream and up to 3kb downstream of TSS within the gene body, set bindingType = "startSite" and bindingRegion = c(-5000, 3000)
        \item To obtain  methylation sites  up to 5kb upstream within the gene body and 3kb downstream of gene/Exon End, set bindingType = "endSite" and bindingRegion = c(-5000, 3000)
        \item To obtain methylation sites from 5kb upstream to 3kb downstream of genes/Exons , set  bindingType = "fullRange" and bindingRegion = c(-5000, 3000)
        \item To obtain methylation sites with nearest bi-directional promoters within 5kb upstream and 3kb downstream of TSS, set  bindingType = "nearestBiDirectionalPromoters" and bindingRegion = c(-5000, 3000)
    }
    \describe{
        \item{startSite}{start position of the feature (strand is considered)} 
        \item{endSite}{end position of the feature (strand is considered)} 
        \item{fullRange}{whole range of the feature} 
        \item{nearestBiDirectionalPromoters}{nearest promoters from both direction of the methylation sites (strand is considered). It will report bidirectional promoters if there are promoters in both directions in the given region (defined by bindingRegion). Otherwise, it will report the closest promoter in one direction.}
    }
  }
  \item{bindingRegion}{
    Annotation range used together with bindingType, which is a vector with two integer values, default to c (-5000, 5000). The first one must be no bigger than 0. And the sec ond one must be no less than 1.  For details, see bindingType.
  }
  \item{ignore.data.strand}{
    ignore the methylation sites strand or not.
  }
  \item{...}{
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
