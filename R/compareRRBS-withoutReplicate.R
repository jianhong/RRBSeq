compareRRBS <-
function(RRBSsummary,col.methycount1=7, col.methycount2=12, col.totalcount1 = 6, col.totalcount2 =11, multiAdj = TRUE, multiAdjMethod="BH", maxP=1, minTotalCount=10)

{
	if(missing(RRBSsummary))
	{
		stop("REBSsummary is required!")
	}
	RRBSsummary[,col.methycount1] =as.numeric(as.character(RRBSsummary[,col.methycount1]))
	RRBSsummary[,col.methycount2] = as.numeric(as.character(RRBSsummary[,col.methycount2]))
	RRBSsummary[,col.totalcount1] = as.numeric(as.character(RRBSsummary[,col.totalcount1]))
	RRBSsummary[,col.totalcount2]= as.numeric(as.character(RRBSsummary[,col.totalcount2]))
	group1 = colnames(RRBSsummary)[col.methycount1]
	group2 = colnames(RRBSsummary)[col.methycount2]
	total1 = colnames(RRBSsummary)[col.totalcount1]
	total2 = colnames(RRBSsummary)[col.totalcount2]
	RRBSsummary = RRBSsummary[RRBSsummary[,col.totalcount1] >= minTotalCount & RRBSsummary[,col.totalcount2] >= minTotalCount,]
	unique.pair = unique(cbind(RRBSsummary[,col.methycount1], RRBSsummary[,col.methycount2],RRBSsummary[,col.totalcount1], RRBSsummary[,col.totalcount2]))
	colnames(unique.pair) = c(group1,group2, total1,total2)
	pvalue =do.call(rbind,lapply(1:dim(unique.pair)[1], function(i) {
		temp = matrix(c(unique.pair[i,2], unique.pair[i,1], unique.pair[i,4] - unique.pair[i,2], unique.pair[i,3]-unique.pair[i,1]),
      	 		nrow = 2,  dimnames =  list(c(group2, group1),  c("methylated", "not methylated")))
		r1 = fisher.test(temp)
		c(unique.pair[i,],r1$p.value, r1$est)
	}))
	colnames(pvalue) = c(group1,group2,total1, total2,"p.value", "odds.ratio")
	
	RRBSsummary = merge(RRBSsummary, pvalue)	
	if (multiAdj & dim(RRBSsummary)[1] >1)
	{
				procs = c(multiAdjMethod)
				ptemp = RRBSsummary[,dim(RRBSsummary)[2]-1]
        			res <- mt.rawp2adjp(ptemp, procs)
        			adjp = unique(res$adjp)
        			colnames(adjp)[1] = "p.value"
        			colnames(adjp)[2] = paste(multiAdjMethod, "adjusted.p.value", sep = ".")
        			temp = merge(RRBSsummary, adjp, all.x = TRUE)
				temp[is.na(temp[, dim(temp)[2]]),dim(temp)[2]] = 1
				temp[temp[,dim(temp)[2]] <=maxP,]
	}
	else
	{
				RRBSsummary[RRBSsummary[,dim(RRBSsummary)[2]-1] <=maxP,]
	}
}
