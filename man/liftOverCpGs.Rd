\name{liftOverCpGs}
\alias{liftOverCpGs}
\title{
    lift over CpGs between genome builds
}
\description{
    Use rtracklayer to lift overthe intervals between genome builds
}
\usage{
liftOverCpGs(x, from, to, ah=AnnotationHub())
}

\arguments{
  \item{x}{
    A \link[GenomicRanges]{GRanges} object.
}
  \item{from, to}{
    The genome builds such as hg18, hg19.
  }
  \item{ah}{
    A \link[AnnotationHub]{AnnotationHub} object.
  }
}
\value{
A GRangesList object.
}
\author{
Jianhong Ou
}
\examples{
if(interactive()){
CpGs <- queryCpGFromUCSC("Mus musculus")
mcols(CpGs)
mm9.CpGs <- CpGs[[1]]
mm10.CpGs <- liftOverCpGs(mm9.CpGs, from="mm9", to="mm10")
mm10.CpGs <- unlist(mm10.CpGs)
}
}
