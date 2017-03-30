\name{compareMethylation}
\alias{compareMethylation}
\title{
    Compare methylation between groups
}
\description{
    Compare methylation between different groups.
}
\usage{
compareMethylation(rrb, coef=1, robust=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rrb}{
    A \link{RRBSeqDataSet} object.
}
  \item{coef}{
    Column number or column name specifying which coefficient or contrast of the linear model is of interest. See \link[limma]{topTable}.
}
  \item{robust}{
    Logical. See \link[limma]{eBayes}.
}
  \item{\dots}{
    Not used.
}
}
\details{
    The comparison will be made depend on the info in the design and contrasts.matrix.
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
}
}
