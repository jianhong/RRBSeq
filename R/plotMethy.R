plotMethy <- function(data, txdb, orgdb, range, gene, 
                      meanCutoff=c(totalC=0, methyC=0, ratio=0.05),
                      ...){
    if(missing(txdb) || missing(orgdb) || missing(data)){
        stop("data, txdb and orgdb are required")
    }
    if(missing(range) && missing(gene)){
        stop("range and gene both are missing.")
    }
    if(any(!c("totalC", "methyC", "ratio") %in% names(meanCutoff)) ||
       !inherits(meanCutoff, c("numeric", "integer"))){
        stop("meanCutoff must be a numeric vector with names totalC, methyC, ratio")
    }
    if(missing(range)){
        org <- gsub(".db", "", deparse(substitute(orgdb)))
        eg <- mget(gene, get(paste0(org, "SYMBOL2EG")), ifnotfound = NA)
        if(is.na(eg[[1]][1])){
            eg <- mget(gene, get(paste0(org, "ALIAS2EG")), ifnotfound = NA)
        }
        if(is.na(eg[[1]][1])){
            stop("can not retrieve location info by ", gene)
        }
        eg <- eg[[1]][1]
        start <- mget(eg, get(paste0(org, "CHRLOC")))[[1]][1]
        end <- mget(eg, get(paste0(org, "CHRLOCEND")))[[1]][1]
        chr <- mget(eg, get(paste0(org, "CHR")))[[1]][1]
        strand <- ifelse(start>0, "+", "-")
        range <- GRanges(chr, IRanges(abs(start)-2000, abs(end)+1000), 
                         strand=strand)
        seqlevelsStyle(range) <- "UCSC"
    }
    stopifnot(inherits(data, c("GRanges", "RRBSeqDataSet")))
    if(class(data)=="GRanges"){
        mc <- mcols(data)
        total.C <- grepl("^totalC", colnames(mc))
        methy.C <- grepl("^methyC", colnames(mc))
        ratio <- grepl("ratio", colnames(mc))
        if(sum(total.C)<1 || sum(methy.C)<1 || sum(ratio)<1){
            stop("Input could be output of compareMethylation which will generate",
                 "a GRanges object with metadata which colnames contains",
                 "totalC, methyC and ratio.")
        }
        total.C <- as.data.frame(mc[, total.C, drop=FALSE])
        methy.C <- as.data.frame(mc[, methy.C, drop=FALSE])
        ratio <- as.data.frame(mc[, ratio, drop=FALSE])
        sampleNames <- gsub("^methyC.", "", colnames(methy.C))
    }else{
        ## class(data)=="RRBSeqDataSet"
        raw <- data
        data <- as(raw, "GRanges")
        total.C <- totalC(raw)
        methy.C <- methyC(raw)
        ratio <- ratio(raw)
        sampleNames <- sampleNames(raw)
    }
    
    keep <- rowMeans(methy.C) > meanCutoff['methyC'] &
        rowMeans(total.C) > meanCutoff['totalC'] &
        rowMeans(ratio) > meanCutoff['ratio']
    data <- reCenterPeaks(data[keep], width=1)
    total.C <- total.C[keep, , drop=FALSE]
    methy.C <- methy.C[keep, , drop=FALSE]
    ratio <- ratio[keep, , drop=FALSE]
    
    if(length(unique(sampleNames))!=ncol(methy.C)){
        stop("Can not get sample names")
    }
    args <- as.list(match.call(expand.dots=FALSE))$`...`
    if(length(args$type)!=1){
        if(length(sampleNames)==1){
            type <- "circle"
        }else{
            type <- "pie.stack"
        }
    }else{
        type <- args$type
    }
    if(length(args$legend)==0){
        legend <- NULL
    }
    data.gr <- rep(data, length(sampleNames))
    color.set <- as.list(as.data.frame(rbind(rainbow(length(sampleNames)), 
                                             "#FFFFFFFF"), 
                                       stringsAsFactors=FALSE))
    names(color.set) <- sampleNames
    if(type %in% c("pie", "pie.stack")){
        type <- "pie.stack"
        unMethyC <- total.C - methy.C
        if(any(unMethyC < 0)){
            stop("Some total counts are smaller than methylated counts")
        }
        colnames(unMethyC) <- colnames(methy.C)
        
        mc <- lapply(as.list(as.data.frame(rbind(methy.C, unMethyC))), 
                     function(.ele) matrix(.ele, ncol=2, byrow=FALSE))
        mcols(data.gr) <- do.call(rbind, mc)
        colnames(mcols(data.gr)) <- c("methylated", "unmethylated")
        data.gr$stack.factor <- rep(sampleNames, each=length(data))
        data.gr$color <- color.set[data.gr$stack.factor]
        legend <- list(labels=sampleNames, col="gray80", 
                       fill=sapply(color.set, `[`, 1))
    }else{
        data.gr <- split(data.gr, rep(sampleNames, each=length(data)))
        data.gr <- mapply(function(.ele, id){
            mcols(.ele) <- methy.C[, id, drop=FALSE]
            .ele$color <- color.set[[id]][1]
            .ele
        }, data.gr, 1:length(data.gr))
    }
    ## get transcripts from ranges
    method <- ifelse(length(args$method)>0, args$method, "BH")
    if(length(args$features)==0){
        suppressMessages(trs <- geneModelFromTxdb(txdb, orgdb, gr=range))
        features <- sapply(trs, function(.ele) .ele$dat)
        features <- unlist(GRangesList(features))
        features$featureLayerID <- paste(features$symbol, 
                                         features$transcript, 
                                         sep=":")
        features$fill <- as.numeric(factor(features$transcript))
        features$height <- ifelse(grepl("utr", features$feature), 0.02, 0.04)
    }
    args$SNP.gr <- data.gr
    args$features <- features
    args$ranges <- range
    args$type <- type
    args$legend <- legend
    if(length(args$ylab)==0 && class(data.gr)=="GRanges") args$ylab <- "methylation"
    args <- as.list(args)
    do.call(lolliplot, args = args)
}
