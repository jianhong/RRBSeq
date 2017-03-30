annotateGeneInfo <- function(data, annotation, 
                             bindingType="nearestBiDirectionalPromoters",
                             bindingRegion=c(-5000, 500),
                             ignore.data.strand=TRUE, ...){
    annoPeaks(data, annotation, 
              bindingType, bindingRegion, 
              ignore.data.strand, ...)
}