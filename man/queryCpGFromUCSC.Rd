\name{queryCpGFromUCSC}
\alias{queryCpGFromUCSC}
\title{
    Query CpGs island from UCSC
}
\description{
    Query CpGs island from UCSC server.
}
\usage{
queryCpGFromUCSC(keywords, ah=AnnotationHub())
}

\arguments{
  \item{keywords}{
    keywords used to search.
}
  \item{ah}{
    A \link[AnnotationHub]{AnnotationHub} object.
  }
}
\value{
A \link[AnnotationHub]{AnnotationHub} object.
}
\author{
Jianhong Ou
}
\examples{
if(interactive()){
CpGs <- queryCpGFromUCSC("Mus musculus")
mcols(CpGs)
}
}
