library(GenomicRanges)
assignCpG2CpGisland <-
function(input.GR, CpGisland.GR,maxgap = 0L, minoverlap = 1L, n.groups=2)
{
	cat(date(), "Validating input ...\n");
	if (missing(input.GR)) {
        stop("Missing required argument input.GR!")
    }
    if (class(input.GR) != "RangedData" && class(input.GR) != "GRanges" ) {
        stop("No valid input.RD passed in. It needs to be GRanges or RangedData object")
    }
    if (missing(CpGisland.GR)) {
        stop("Missing requirement argument CpGisland.RD!")
		}
    if (class(CpGisland.GR) != "RangedData" && class(input.GR) != "GRanges" ) {
        stop("No valid CpGisland.RD passed in. It needs to be GRanges or RangedData object")
    }		
	
	ol = findOverlaps(as(input.GR,"GRanges"), as(CpGisland.GR,"GRanges"))
	ol.matrix = as.matrix(ol)
	colnames(ol.matrix) = c("CpGName", "CpGislandName")
	t1 = as.data.frame(input.RD)
	t1 = cbind(as.character(rownames(t1)), t1)
	colnames(t1)[1:3] = c("CpGName", "chr", "pos")
	t2 = as.data.frame(CpGisland.RD)
	t2 = cbind(as.character(rownames(t2)), t2)
	colnames(t2)[1:2] = c("CpGislandName", "chr")
	t3= merge(ol.matrix, t1[,1:3])
	t4 = merge(t3, t2)	
	CpG.in.island = dim(t3)[1]/dim(t1)[1] * 100
	list(CpG.in.island.percent = CpG.in.island, C.in.CpGisland = t4)
}

summarizeCpGisland <- function(CinCpGisland, totalC.columns=c(6,11,16,21,26,31,36,41), ratio.columns = c(5,10,15,20,25,30,35,40), methyC.columns =c(7,12,17,22,27,32,37,42),treatments=c("STC1-AC", "STC1-WT","T0-AC", "T0-WT", "T14-AC", "T14-WT", "T21-AC", "T21-WT"), CpGisland.column=45)
{
	methylation.summary = data.frame()
        for (i in 1:length(totalC.columns))
        {
                total.Ratio = rowsum(CinCpGisland[,ratio.columns[i]],group=CinCpGisland[,CpGisland.column])
                total = table(CinCpGisland[,CpGisland.column])
                meanRatio = as.numeric(total.Ratio) /as.numeric(total)
		totalC = rowsum(CinCpGisland[,totalC.columns[i]],group=CinCpGisland[,CpGisland.column]) 
		methyC = rowsum(CinCpGisland[,methyC.columns[i]],group=CinCpGisland[,CpGisland.column])
		meanRatioFromTotalC = as.numeric(methyC)/as.numeric(totalC)
                if (i ==1)
                {
                        methylation.summary = cbind(rownames(total.Ratio), total, meanRatio,totalC, methyC,meanRatioFromTotalC)
                }
                else
                {
                        methylation.summary = cbind(methylation.summary,total, meanRatio,totalC, methyC, meanRatioFromTotalC)
                }
                k = 5*i - 3
                colnames(methylation.summary)[k:(k+4)] = c(paste("totalUniqueC", treatments[i],sep="-"),paste("mean.Ratio",treatments[i],sep="-"),  paste("totalC",treatments[i],sep="-"), paste("methyC",treatments[i], sep="-"),paste("mean.Ratio.FromTotalC", treatments[i], sep="-"))
        }
        colnames(methylation.summary)[1] =colnames(CinCpGisland)[CpGisland.column] 	
	ind = dim(CinCpGisland)[2]
	if (ind < CpGisland.column)
	{
		ind = CpGisland.column
	}
	unique.CpGisland = unique(cbind(CinCpGisland[,1],CinCpGisland[,CpGisland.column:ind]))
	colnames(unique.CpGisland)[1]="chr"
	CpGisland.summary = merge(methylation.summary, unique.CpGisland)
 	CpGisland.summary
}
