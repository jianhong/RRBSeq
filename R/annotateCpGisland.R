annotateCpGisland <- function(data, CpGs, ...){
    stopifnot(class(data)=="GRanges")
    stopifnot(class(CpGs)=="GRanges")
    stopifnot(length(CpGs$name)==length(CpGs))
    stopifnot(length(CpGs)>0)
    stopifnot(length(data)>0)
    CpGs$name <- as.character(CpGs$name)
    ol <- findOverlaps(data, CpGs)
    anno <- CpGs[subjectHits(ol)]
    anno.group <- split(anno, queryHits(ol))
    anno.name <- sapply(anno.group,
                        function(.ele) paste(.ele$name,
                                             collapse=";"))
    anno.ranges <- lapply(anno.group, ranges)
    anno.id <- as.numeric(names(anno.group))
    data$CpGsName <- rep("", length(data))
    data[anno.id]$CpGsName <- anno.name
    data$CpGsRanges <- IRangesList(rep(IRanges(), length(data)))
    data[anno.id]$CpGsRanges <- IRangesList(anno.ranges)
    data
}