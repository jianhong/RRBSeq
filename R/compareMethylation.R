compareMethylation <- function(rrb, coef=1, robust=FALSE, ...){
    if(missing(rrb)){
        stop("rrb is required!")
    }
    if (class(rrb) != "RRBSeqDataSet" ) {
        stop("No valid rrb passed in. It needs to be RRBSeqDataSet object")
    }
    
    if (length(design(rrb)) == 0){
        stop("design is required, please refer to limma package for create design ")
    }
    
    gr <- rowRanges(rrb)
    eset <- ratio(rrb)
    colnames(eset) <- paste("normalized.ratio", sampleNames(rrb), sep=".")
    methyC <- methyC(rrb)
    colnames(methyC) <- paste("methyC", sampleNames(rrb), sep=".")
    totalC <- totalC(rrb)
    colnames(totalC) <- paste("totalC", sampleNames(rrb), sep=".")
    
    if(any(colSums(design(rrb))>1)){ ## with duplicates
        fit1 <- lmFit(eset, design(rrb))
        if(length(contrasts.matrix(rrb))>0){
            fit1 <- contrasts.fit(fit1, contrasts = contrasts.matrix(rrb))
        }
        fit2 <- eBayes(fit1, robust=robust)
        ans <- topTable(fit2, number = nrow(eset), coef=coef, sort.by="none")
        colnames(ans) <- gsub("logFC", "difference", colnames(ans))
        ans <- ans[, !colnames(ans) %in% c("t", "B"), drop=FALSE]
    }else{
        ## only for simple comparison
        if(ncol(design(rrb))>2){
            stop("design is complicated. Please try to simplify the design")
        }
        if(any(totalC!=floor(totalC)) || any(methyC!=floor(methyC))){
            stop("The counts for total or methy should be integer.")
        }
        unMethyC <- totalC - methyC
        if(any(unMethyC<0)){
            stop("There are totalC smaller than methyC.",
                 "Please double check your inputs.")
        }
        ans <- apply(cbind(methyC, unMethyC), 1, function(.ele){
            ft <- fisher.test(matrix(.ele, nrow=2))
            c(odds.ratio=unname(ft$estimate),
              P.Value=ft$p.value,
              adj.P.Val=NA)
        })
        if(any(colnames(ans)!=c("odds.ratio", "P.Value", "adj.P.Val"))){
            ans <- t(ans)
        }
        if(nrow(ans)!=length(gr)){
            stop("Bug found.")
        }
        args <- as.list(match.call(expand.dots=TRUE))
        method <- ifelse(length(args$method)>0, args$method, "BH")
        ans[, "adj.P.Val"] <- p.adjust(ans[, "P.Value"], method=method)
    }
    
    mcols(gr) <- DataFrame(eset, methyC, totalC, ans)
    gr
}
