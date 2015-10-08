#' Calculate deltaWs
#'
#' This function will calculate deltaWs from a \code{\link{GRanges}} object with read fragments.
#'
#' @param frags A \code{\link{GRanges}} with read fragments (see \code{\link{bam2GRanges}}).
#' @param reads.per.window Number of reads in each dynamic window.
#' @import GenomicRanges
#' @author Aaron Taudt
#' @export
deltaWCalculator <- function(frags, reads.per.window=10) {

	frags.split <- split(frags, seqnames(frags))
	reads.per.chrom <- sapply(frags.split, length)
	chroms2parse <- names(reads.per.chrom)[reads.per.chrom>2*reads.per.window]
	chroms2skip <- setdiff(names(reads.per.chrom),chroms2parse)
	if (length(chroms2skip)>0) {
		warning(paste0("Not parsing chromosomes ",paste(chroms2skip, collapse=',')," because they do not have enough reads."))
	}

	frags.new <- GRangesList()
	for (chrom in chroms2parse) {
		f <- frags.split[[chrom]]
		f <- f[order(start(f))]
		f$pcsum <- NA
		f$mcsum <- NA
		f$preads <- NA
		f$mreads <- NA
		for (i1 in 1:reads.per.window) {
			## Get index of window positions
			windex <- seq(from=i1, by=reads.per.window, to=length(f))
			## Make plus and minus strand cumulative sum
			f$pcsum <- cumsum(strand(f)=='+')
			f$mcsum <- cumsum(strand(f)=='-')
			if (i1>1) {
				f$preads[windex] <- diff(as.numeric(f$pcsum[c(1,windex)]))
				f$mreads[windex] <- diff(as.numeric(f$mcsum[c(1,windex)]))
			} else {
				f$preads[windex] <- c(0,diff(as.numeric(f$pcsum[c(windex)])))
				f$mreads[windex] <- c(0,diff(as.numeric(f$mcsum[c(windex)])))
			}
		}
		## Right - left window with padding of the ends
		pfill <- (f$preads[length(f)]:0)[1:reads.per.window]
		pfill[is.na(pfill)] <- 0
		mfill <- (f$mreads[length(f)]:0)[1:reads.per.window]
		mfill[is.na(mfill)] <- 0
		f$pdelta <- abs(c(f$preads[-c(1:reads.per.window)],pfill) - f$preads)
		f$mdelta <- abs(c(f$mreads[-c(1:reads.per.window)],mfill) - f$mreads)

		## Minus and plus merged
		f$deltaW <- pmin(f$pdelta,f$mdelta)

		frags.new[[chrom]] <- f
	}
	frags.new <- unlist(frags.new)
	names(frags.new) <- NULL

	return(frags.new)

}




# ### Ashley Sanders & David Porubsky
# ### May 12, 2015
# 
# ################	################	################	################  ################	################	################	################
# ################  this function will calculate deltaWs from fileFreqs                                    ################
# ################  it will bin the data for this calculation and compare LHS and RHS of bin               ################
# ################ it will generate a file with new start/end Positions and deltaW for region              ################
# ################  ################	################	################  ################	################	################	################
# 
# 
# deltaWCalculator<- function(fileFreqs, chr, window=10, gapfile=F) 
#   
#   # fileFreqs contains chr in  2nd col, startPos  in 3rd col, and strand  in 4th col
#   # window is the % of readDepth used to calculate the sliding window  
#   # if gapfile is provided (length>1) then will remove these regions from the deltaWs (overlaps set to 0)
# 
# {
#   start <- 1
#   end <- nrow(fileFreqs)
#   
#   win<- round(nrow(fileFreqs)/window) # window size assigned based on a percentage of total reads in file
#   if (win %% 2 != 0) {win <- win+1} # ensures even number for window
#   
#   ## using a sliding window, calculate the change in W reads (note: checked and deltaC yields the same result)
#   ranges<- successiveIRanges(rep(win, (nrow(fileFreqs)-win+1)), gapwidth=1-win)    
#   
#   startPos<- fileFreqs[end(ranges)-win*0.5,3] # pulls out the startPos for the deltaW ranges
#   endPos<- fileFreqs[start(ranges)+win*0.5,3] # pulls out the endPos for the deltaW ranges
#   ## calculate deltaWs    
#   firstWs <- sapply(start(ranges), function(x) table(fileFreqs[start(ranges)[x]:(end(ranges)-win*0.5)[x],4])[2]) ## number of W in the first 1/2 of the ranges (i.e. start(ranges) to end(ranges)-win*05))
#   secondWs <- sapply(start(ranges), function(x) table(fileFreqs[(start(ranges)+win*0.5)[x]:end(ranges)[x],4])[2])  ## number of W in the second 1/2 of the ranges (i.e. start(ranges)+win*0.5 to end(ranges))
#   dWs<- abs( firstWs-secondWs) # calculate difference in Ws    
#   
#   deltaWs <- data.frame(startPos=startPos, endPos=endPos, deltaW=dWs) 
#   
#   
#   if (length(gapfile) > 1)
#   {
#     gapsChr<- gapfile[which(gapfile[,2] == as.character(fileFreqs[1,2])),]
#     if(nrow(gapsChr) > 0) #'*'#
#     {
#       gapRle  <- GRanges(chr, IRanges(start=gapsChr[,3], end=gapsChr[,4]))
#       
#       dWs = GRanges(chr, IRanges(start=deltaWs[,1], end=deltaWs[,2]), score=deltaWs[,3]) 
#       hits<- findOverlaps(dWs, gapRle)
#       deltaWs[queryHits(hits),3] <-0 # turns any deltaWs in gaps to 0
#     }   
#   }
#   return(deltaWs)
# } 
# 

