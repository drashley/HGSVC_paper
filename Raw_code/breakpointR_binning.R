# aneufinder - An R-package for CNV detection in whole-genome single cell sequencing data
# Copyright (C) 2015  Aaron Taudt
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


bam2fragment <- function(bamfile, bamindex=bamfile, pairedEndReads=FALSE, chromosomes=NULL, min.mapq=10, remove.duplicate.reads=TRUE) {

	library(Rsamtools)
	library(GenomicAlignments)
	library(GenomicRanges)
	## Check if bamindex exists
	bamindex.raw <- sub('\\.bai$', '', bamindex)
	bamindex <- paste0(bamindex.raw,'.bai')
	if (!file.exists(bamindex)) {
		bamindex.own <- Rsamtools::indexBam(bamfile)
		warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
		bamindex <- bamindex.own
	}
	file.header <- Rsamtools::scanBamHeader(bamfile)[[1]]
	chrom.lengths <- file.header$targets
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('The specified chromosomes ', chrstring, ' do not exist in the data.')
	}
	gr <- GenomicRanges::GRanges(seqnames=Rle(chroms2use), ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
	if (!remove.duplicate.reads) {
		if (pairedEndReads) {
			data <- GenomicAlignments::readGAlignmentPairsFromBam(bamfile, index=bamindex, param=ScanBamParam(which=range(gr), what='mapq'))
			data <- GenomicAlignments::first(data)	# take only first mapping fragment of each pair
		} else {
			data <- GenomicAlignments::readGAlignmentsFromBam(bamfile, index=bamindex, param=ScanBamParam(which=range(gr), what='mapq'))
		}
	} else {
		if (pairedEndReads) {
			data <- GenomicAlignments::readGAlignmentPairsFromBam(bamfile, index=bamindex, param=ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
			data <- GenomicAlignments::first(data)	# take only first mapping fragment of each pair
		} else {
			data <- GenomicAlignments::readGAlignmentsFromBam(bamfile, index=bamindex, param=ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
		}
	}
	## Filter by mapping quality
	if (!is.null(min.mapq)) {
		data <- data[mcols(data)$mapq >= min.mapq]
	}
	## Remove duplicate reads
	if (remove.duplicate.reads) {
		sp <- start(data)[as.logical(strand(data)=='+')]
		sp1 <- c(sp[length(sp)], sp[-length(sp)])
		sm <- start(data)[as.logical(strand(data)=='-')]
		sm1 <- c(sm[length(sm)], sm[-length(sm)])
		data <- c(data[strand(data)=='+'][sp!=sp1], data[strand(data)=='-'][sm!=sm1])
	}
	## Store the read fragments as GRanges
	fragments <- GRanges(seqnames=seqnames(data), ranges=IRanges(start=start(data), end=end(data)), strand=strand(data))
		return(fragments)

}



getDynamicWindows <- function(frags, reads.per.window=10) {

# 	## TEST
# 	# Make artificial strandseq
# 	pfrags <- frags[strand(frags)=='+']
# 	mfrags <- frags[strand(frags)=='-']
# 	frags <- c(pfrags[1:100], mfrags[201:300], mfrags[c(23,25,54,60)], pfrags[c(254,256)])
# 	frags <- frags[order(start(frags))]

	frags.new <- GRangesList()
	for (chrom in unique(seqnames(frags))) {
		f <- frags[seqnames(frags)==chrom]
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
		## Right - left window
		pfill <- (f$preads[length(f)]:0)[1:reads.per.window]
		pfill[is.na(pfill)] <- 0
		mfill <- (f$mreads[length(f)]:0)[1:reads.per.window]
		mfill[is.na(mfill)] <- 0
		f$pdelta <- abs(c(f$preads[-c(1:reads.per.window)],pfill) - f$preads)
		f$mdelta <- abs(c(f$mreads[-c(1:reads.per.window)],mfill) - f$mreads)
		
		## Minus and plus merged
		f$delta <- apply(matrix(c(f$pdelta,f$mdelta),ncol=2),1,min)

		frags.new[[chrom]] <- f
	}
	frags.new <- unlist(frags.new)
	names(frags.new) <- NULL

	return(frags.new)

}


