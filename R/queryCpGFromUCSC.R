queryCpGFromUCSC <- function(keywords, ah=AnnotationHub()){
    query(ah, c("UCSC", keywords, "CpG Islands"))
}

liftOverCpGs <- function(x, from, to, ah=AnnotationHub()){
    stopifnot(class(x)=="GRanges")
    chainfiles <- query(ah , c(paste(from, "to", to), "chainfile"))
    if(length(chainfiles)>1){
        message("chainfiles more than 1. Returning record of AnnotationHub")
        return(chainfiles)
    }
    if(length(chainfiles)==0){
        stop("No results!")
    }
    chain <- chainfiles[[1]]
    liftOver(x, chain)
}