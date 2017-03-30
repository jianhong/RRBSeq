\name{RRBSeqDataSet-class}
\docType{class}

% Class
\alias{RRBSeqDataSet-class}
\alias{class:RRBSeqDataSet}

% Constructor
\alias{RRBSeqDataSet}

% Accessors
\alias{ratio}
\alias{methyC}
\alias{totalC}
\alias{ratio<-}
\alias{methyC<-}
\alias{totalC<-}
\alias{ratio,RRBSeqDataSet-method}
\alias{methyC,RRBSeqDataSet-method}
\alias{totalC,RRBSeqDataSet-method}
\alias{ratio<-,RRBSeqDataSet,matrix-method}
\alias{methyC<-,RRBSeqDataSet,matrix-method}
\alias{totalC<-,RRBSeqDataSet,matrix-method}

\alias{design}
\alias{design<-}
\alias{design,RRBSeqDataSet-method}
\alias{design<-,RRBSeqDataSet,matrix-method}
\alias{design<-,RRBSeqDataSet,NULL-method}

\alias{contrasts.matrix}
\alias{contrasts.matrix<-}
\alias{contrasts.matrix,RRBSeqDataSet-method}
\alias{contrasts.matrix<-,RRBSeqDataSet,matrix-method}
\alias{contrasts.matrix<-,RRBSeqDataSet,NULL-method}

\alias{sampleNames}
\alias{sampleNames<-}
\alias{sampleNames,RRBSeqDataSet-method}
\alias{sampleNames<-,RRBSeqDataSet,character-method}

% data transform methods
\alias{transformData}
\alias{summarize.sliding.window}
\alias{normalizeData}
\alias{transformData,RRBSeqDataSet-method}
\alias{summarize.sliding.window,RRBSeqDataSet-method}
\alias{normalizeData,RRBSeqDataSet-method}

% Coercion
\alias{coerce,RRBSeqDataSet,GRanges-method}

\title{
    RRBSeqDataSet class
}
\description{
    The RRBSeqDataSet class is a matrix-like container where \link[SummarizedExperiment]{rowRanges}
    represent ranges of interest and columns represent samples (with sample data summarized as a 
    \link{DataFrame}). A RRBSeqDataSet contains three assays, they are ratio 
    (ratio of methylated to total), methyC (methylated counts) and totalC (total counts). 
    Each represented by a matrix-like object of numeric mode.
    
    RRBSeqDataSet is a subclass of \link[SummarizedExperiment]{RangedSummarizedExperiment} and,
    thus, it inherites all the methods of RangedSummarizedExperiment class. The methods documented
    below are additional methods that are specific to RRBSeqDataSet object.
}
\usage{
## Constructor
RRBSeqDataSet(methyRatioFiles, 
              ind.chr=1L, ind.pos=2L,  
              ind.strand=3L, ind.totalC=6L, 
              ind.methyC=7L, ind.ratio=5L, 
              min.totalC=5L, ...,
              colData=DataFrame(methyRatioFiles=methyRatioFiles),
              metadata=list(),
              design=NULL,
              contrasts.matrix=NULL)

## Accessors
methyC(object, ...)
methyC(object, ...) <- value

totalC(object, ...)
totalC(object, ...) <- value

ratio(object, ...)
ratio(object, ...) <- value

design(object, ...)
design(object, ...) <- value

contrasts.matrix(object, ...)
contrasts.matrix(object, ...) <- value

sampleNames(object)
sampleNames(object) <- value

## methods
transformData(object, transformation="asin", ...)
summarize.sliding.window(object, window.size =200L, 
       step=50L, summarize=c("mean.ratio", "total.methyC.over.total.totalC"),
       ...)
normalizeData(object, method=c("none", "quantiles", 
       "quantiles.robust", "mean", "median"),  ... )
}
\arguments{
  \item{object}{
  A \code{RRBSeqDataSet} object.
}
  \item{methyRatioFiles}{
  A character vector describe the file paths of methylation ratio files.
}
  \item{ind.chr, ind.pos, ind.strand, ind.totalC, ind.methyC, ind.ratio}{
  Columns of chromosome, position, strand, total counts, methylated counts, ratios.
}
  \item{min.totalC}{
  The minimal count of total reads in each position.
  }
  \item{\dots}{
  In constructor RRBSeqDataSet, all the parameters could be passed to read.table.
  In method normalizeData, all the parameters could be passed to \link[preprocessCore]{normalize.quantiles.robust}.
}
  \item{colData}{
    An optional \link[S4Vectors]{DataFrame} describing the samples.
  }
  \item{metadata}{
    An optional \code{list} of arbitrary content describing the overall experiment.
  }
  \item{design}{
    A \code{matrix} which expresses how the counts for each position depend on the samples. 
    See \link[stats]{model.matrix}
  }
  \item{contrasts.matrix}{
    A \code{matrix} which columns corressponding to contrasts. See \link[limma]{makeContrasts}
  }
  \item{transformation}{
    How to transform ratio values. Could be asin or log2.
  }
  \item{window.size}{
    Window size.
  }
  \item{step}{
    Step size to move the window
  }
  \item{summarize}{
    How to summarize the ratio in a window. Could be mean.ratio or total.methyC.over.total.totalC.
  }
  \item{method}{
    Normalization method. Could be none, quantiles, quantiles.robust, mean or median.
  }
  \item{value}{
    The proper values.
  }
}
\value{
A RRBSeqDataSet object
}

\section{Coercion}{
    RRBSeqDataSet object can be coerced to \link[GenomicRanges]{GRanges} object by \code{as(object, to="GRanges")}.
}

\author{
Jianhong Ou
}

\examples{
    path <- system.file("extdata", package="RRBSeq")
    methyRatioFiles <- dir(path, "^test", full.names=TRUE)
    rrb <- RRBSeqDataSet(methyRatioFiles, ind.totalC=5, 
                         ind.methyC=4, ind.ratio=6, header=TRUE)
}