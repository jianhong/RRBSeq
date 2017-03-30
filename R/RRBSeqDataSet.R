## require SummarizedExperiment
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClass("RRBSeqDataSet", contains = "RangedSummarizedExperiment", 
         representation = representation(
             design = "matrixOrNULL",
             contrasts.matrix = "matrixOrNULL"
         )
    )

setValidity2("RRBSeqDataSet", function(object){
    if(class(rowRanges(object))!="GRanges"){
        return("rowRanges must be an object of GRanges, but not GRangesList")
    }
    if(!mode(assays(object)$totalC) %in% c("numeric", "integer")){
        return("totalC must be an numeric matrix")
    }
    if(!mode(assays(object)$methyC) %in% c("numeric", "integer")){
        return("methyC must be an numeric matrix")
    }
    if(!mode(assays(object)$ratio) %in% c("numeric", "integer")){
        return("ratio must be an numeric matrix")
    }
    if(length(colData(object)$sampleNames)!=ncol(assays(object)$totalC)){
        return("sampleNames is not correct.")
    }
    return(TRUE)
})


## new
RRBSeqDataSet <- function(methyRatioFiles, 
                          ind.chr=1L, ind.pos=2L,  
                          ind.strand=3L, ind.totalC=6L, 
                          ind.methyC=7L, ind.ratio=5L, 
                          min.totalC=5L, ...,
                          colData=DataFrame(methyRatioFiles=methyRatioFiles),
                          metadata=list(),
                          design=NULL,
                          contrasts.matrix=NULL){
    if (min.totalC < 1){
        warning("min.totalC needs to be at least 1, so automatically set to 1")
        min.totalC = 1
    }
    if(missing(methyRatioFiles)){
        stop("methyRatioFiles is required.")
    }
    stopifnot(class(colData)=="DataFrame")
    if(length(colData$sampleNames)!=length(methyRatioFiles)){
        colData$sampleNames <- gsub("\\.[^.]*$", "", 
                                    basename(methyRatioFiles))
    }
    allInputData <- lapply(methyRatioFiles, function(.ele){
        read.table(.ele, ...)[, c(ind.chr, ind.pos, ind.strand, 
                                  ind.totalC, ind.methyC, ind.ratio)]
    })
    sharedPos.l <- lapply(allInputData, function(.ele) 
        paste(.ele[, 1], .ele[, 2], .ele[, 3], sep="_"))
    sharedPos <- table(unlist(sharedPos.l))
    sharedPos <- sharedPos[sharedPos>=length(methyRatioFiles)/2]
    sharedPos <- names(sharedPos)
    subData <- mapply(function(a, b){
        b[match(sharedPos, a), 4:6]
    }, sharedPos.l, allInputData, SIMPLIFY = FALSE)
    getData <- function(j){
        do.call(cbind, lapply(subData, function(.ele){
            .ele <- as.numeric(as.character(.ele[, j]))
            .ele[is.na(.ele)] <- 0
            .ele
        }))
    }
    totalC <- getData(1)
    methyC <- getData(2)
    ratio <- getData(3)
    colnames(totalC) <- colnames(methyC) <- colnames(ratio) <- 
        make.names(colData$sampleNames)
    fil <- apply(totalC, 1, function(.ele) all(.ele >= min.totalC))
    
    sharedPos <- sharedPos[fil]
    
    if(length(sharedPos)<1){
        stop("After filter, no useful data left.")
    }
    
    sharedPos <- do.call(rbind, strsplit(sharedPos, "_"))
    totalC <- totalC[fil, , drop=FALSE]
    methyC <- methyC[fil, , drop=FALSE]
    ratio <- ratio[fil, , drop=FALSE]
    
    GR  <- GRanges(seqnames=as.character(sharedPos[, 1]),
                   IRanges(start=as.numeric(as.character(sharedPos[ , 2])), 
                           width=1),  
                   strand=as.character(sharedPos[, 3]))
    
    se <- SummarizedExperiment(assays=list(totalC=totalC,
                                           methyC=methyC,
                                           ratio=ratio),
                               rowRanges=GR,
                               colData=colData,
                               metadata=metadata)
    
    object <- new("RRBSeqDataSet", se, 
                  design=design,
                  contrasts.matrix=contrasts.matrix)
    
    # stash the package version
    #metadata(object)[["version"]] <- packageVersion("RRBSeq")
    
    return(object)
}

## 
setAs(from="RRBSeqDataSet", to="GRanges", function(from){
    gr <- rowRanges(from)
    eset <- ratio(from)
    colnames(eset) <- paste("normalized.ratio", sampleNames(from), sep=".")
    methyC <- methyC(from)
    colnames(methyC) <- paste("methyC", sampleNames(from), sep=".")
    totalC <- totalC(from)
    colnames(totalC) <- paste("totalC", sampleNames(from), sep=".")
    mcols(gr) <- DataFrame(eset, methyC, totalC)
    gr
})

## accessors
setGeneric("ratio", function(object, ...) standardGeneric("ratio"))
setGeneric("methyC", function(object, ...) standardGeneric("methyC"))
setGeneric("totalC", function(object, ...) standardGeneric("totalC"))
setMethod("ratio", "RRBSeqDataSet", function(object, ...){
    assays(object)$ratio
})
setMethod("methyC", "RRBSeqDataSet", function(object, ...){
    assays(object)$methyC
})
setMethod("totalC", "RRBSeqDataSet", function(object, ...){
    assays(object)$totalC
})
setGeneric("ratio<-", function(object, ..., value) standardGeneric("ratio<-"))
setGeneric("methyC<-", function(object, ..., value) standardGeneric("methyC<-"))
setGeneric("totalC<-", function(object, ..., value) standardGeneric("totalC<-"))
setReplaceMethod("ratio", c("RRBSeqDataSet", "matrix"), function(object, ..., value){
    assays(object)$ratio <- value
    object
})
setReplaceMethod("methyC", c("RRBSeqDataSet", "matrix"), function(object, ..., value){
    assays(object)$methyC <- value
    object
})
setReplaceMethod("totalC", c("RRBSeqDataSet", "matrix"), function(object, ..., value){
    assays(object)$totalC <- value
    object
})

#setGeneric("design", function(object, ...) standardGeneric("design"))
#setGeneric("design<-", function(object, ..., value) standardGeneric("design<-"))
# importMethodFrom(BiocGenerics)
setMethod("design", "RRBSeqDataSet", function(object, ...){
    object@design
})
setReplaceMethod("design", c("RRBSeqDataSet", "matrix"), function(object, ..., value){
    stopifnot(nrow(value)==ncol(assays(object)$ratio))
    object@design <- value
    object
})
setReplaceMethod("design", c("RRBSeqDataSet", "NULL"), function(object, ..., value){
    object@design <- value
    object
})
setGeneric("contrasts.matrix", function(object, ...) standardGeneric("contrasts.matrix"))
setGeneric("contrasts.matrix<-", function(object, ..., value) standardGeneric("contrasts.matrix<-"))
setMethod("contrasts.matrix", "RRBSeqDataSet", function(object, ...){
    object@contrasts.matrix
})
setReplaceMethod("contrasts.matrix", c("RRBSeqDataSet", "matrix"), function(object, ..., value){
    object@contrasts.matrix <- value
    object
})
setReplaceMethod("contrasts.matrix", c("RRBSeqDataSet", "NULL"), function(object, ..., value){
    object@contrasts.matrix <- value
    object
})

## importMethodFrom Biobase
setMethod("sampleNames", c("RRBSeqDataSet"), function(object){
    colData(object)$sampleNames
})
setReplaceMethod("sampleNames", c("RRBSeqDataSet", "character"), function(object, value){
    colData(object)$sampleNames <- value
})

## method: transform data

setGeneric("transformData", function(object, transformation="asin", ...) standardGeneric("transformData"))
.transformData <- function(object, transformation="asin", ...)
{
    stopifnot(class(object)=="RRBSeqDataSet")
    stopifnot(ncol(assays(object)$ratio)>1)
    transformation <- match.arg(transformation, c("asin", "log2"))
    assays(object)$ratio[assays(object)$ratio==0] <- 
        1/(4*assays(object)$totalC[assays(object)$ratio==0])
    assays(object)$ratio[assays(object)$ratio==1] <- 
        1 - 1/(4*assays(object)$totalC[assays(object)$ratio==1])
    assays(object)$ratio <- switch(transformation, 
                                   asin=asin(sqrt(assays(object)$ratio)),
                                   log2=log2(assays(object)$ratio))
    object
}
setMethod("transformData", "RRBSeqDataSet", .transformData)

## method: summarize.sliding.window
setGeneric("summarize.sliding.window", function(object, window.size =200L, step=50L, summarize=c("mean.ratio", "total.methyC.over.total.totalC"), ...) standardGeneric("summarize.sliding.window"))
.summarize.sliding.window <- function(object, window.size =200L, step=50L, summarize=c("mean.ratio", "total.methyC.over.total.totalC"), ...)
{
    summarize = match.arg(summarize)
    if (missing(object)) {
        stop("Missing required argument object!")
    }
    if (class(object) != "RRBSeqDataSet" ) {
        stop("No valid object passed in. It needs to be RRBSeqDataSet object")
    }
    if ( class(window.size) != "integer" || window.size < 2) {
        stop("window.size needs to be an integer greater than 1")
    }
    if ( class(step) != "integer" || step > window.size) {
        stop("step needs to be an integer less or equal to window")
    }
    
    gr <- rowRanges(object)
    strand(gr) <- "*"
    rg <- range(gr)
    wid <- 10^nchar(as.character(window.size))
    start(rg) <- floor(start(rg)/wid)*wid
    end(rg) <- ceiling(end(rg)/wid)*wid
    slidingWindows <- function(object, width, step=1L, ...){ ## break from GRanges::slidingWindows
        n <- ceiling(pmax(width(object) - width, 0L)/step) + 1L
        window.starts <- as.integer(IRanges(rep(0L, length(n)), 
                                            width = n)) * step + 1L
        windows <- restrict(IRanges(window.starts, width = width), 
                            end = rep(width(object), n))
        windows.abs <- shift(windows, rep(start(object), n))
        win <- relist(windows.abs, PartitioningByWidth(n))
        gr <- rep(granges(object), lengths(win))
        ranges(gr) <- unlist(win)
        gr
    }
    windows <- slidingWindows(rg, width=window.size, step=step)
    ol <- findOverlaps(windows, gr, maxgap = 0L, minoverlap = 1L)
    methyC <- assays(object)$methyC
    totalC <- assays(object)$totalC
    ratio <- assays(object)$ratio
    if(length(ol)<1){
        stop("Bugs found! overlaps less than 1")
    }
    ugroup <- sort(unique(queryHits(ol), na.last=TRUE, method="quick"))
    ugroup.hits.no <- table(queryHits(ol))[as.character(ugroup)] ##need to comfirm that no 1e10 ...
    this.methyC <- rowsum(methyC[subjectHits(ol), , drop=FALSE], 
                          group=queryHits(ol), na.rm = TRUE)
    this.totalC <- rowsum(totalC[subjectHits(ol), , drop=FALSE],
                          group=queryHits(ol), na.rm = TRUE)
    this.ratio <- switch(summarize,
                         mean.ratio={
                             sums <- rowsum(ratio[subjectHits(ol), , drop=FALSE],
                                            group=queryHits(ol), na.rm = TRUE)
                             stopifnot(identical(rownames(sums), names(ugroup.hits.no)))
                             sums/as.numeric(ugroup.hits.no)
                         },
                         total.methyC.over.total.totalC={
                             this.methyC/this.totalC
                         })
    
    GR  <- windows[ugroup]
    
    se <- SummarizedExperiment(assays=list(totalC=this.totalC,
                                           methyC=this.methyC,
                                           ratio=this.ratio),
                               rowRanges=GR,
                               colData=colData(object),
                               metadata=metadata(object))
    
    rrb <- new("RRBSeqDataSet", se, 
               design=design(object),
               contrasts.matrix=contrasts.matrix(object))
    
    # stash the package version
    #metadata(rrb)[["version"]] <- packageVersion("RRBSeq")
    
    return(rrb)
}
setMethod("summarize.sliding.window", "RRBSeqDataSet", .summarize.sliding.window)

## method: normalization
setGeneric("normalizeData", function(object, method=c("none", "quantiles", "quantiles.robust", "mean", "median"), ...) standardGeneric("normalizeData"))
.normalizeData <- function(object, method=c("none", "quantiles", "quantiles.robust", "mean", "median"),  ... )
{
    stopifnot(class(object)=="RRBSeqDataSet")
    method <- match.arg(method)
    exprs <- assays(object)$ratio
    methyC <- assays(object)$methyC
    totalC <- assays(object)$totalC
    
    exprs <- switch(method,
                    quantiles=normalize.quantiles(exprs),
                    quantiles.robust=normalize.quantiles.robust(exprs, ...),
                    mean={
                        avgs <- colSummarizeAvg(exprs)$Estimates
                        scaling.factors <- avgs /avgs[1]
                        scaling.factors <- matrix(rep(scaling.factors, nrow(exprs)), 
                                                  ncol=length(scaling.factors), 
                                                  byrow=TRUE)
                        exprs / scaling.factors
                    },
                    median={
                        avgs <- colSummarizeMedian(exprs)$Estimates
                        scaling.factors <- avgs /avgs[1]
                        scaling.factors <- matrix(rep(scaling.factors, nrow(exprs)), 
                                                  ncol=length(scaling.factors), 
                                                  byrow=TRUE)
                        exprs <- exprs /scaling.factors
                    },
                    exprs)
    
    assays(object)$ratio <- exprs
    object
}
setMethod("normalizeData", "RRBSeqDataSet", .normalizeData)
