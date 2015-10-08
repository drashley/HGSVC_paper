#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param keep.duplicate.reads A logical indicating whether or not duplicate reads should be kept.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairsFromBam readGAlignmentsFromBam first
#' @author Aaron Taudt
#' @export
bam2GRanges <- function(file, bamindex=file, chromosomes=NULL, pairedEndReads=FALSE, keep.duplicate.reads=TRUE) {

	## Check if bamindex exists
	bamindex.raw <- sub('\\.bai$', '', bamindex)
	bamindex <- paste0(bamindex.raw,'.bai')
	if (!file.exists(bamindex)) {
		bamindex.own <- Rsamtools::indexBam(file)
		warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
		bamindex <- bamindex.own
	}
	file.header <- Rsamtools::scanBamHeader(file)[[1]]
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
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
	}
	## Import the file into GRanges
	gr <- GenomicRanges::GRanges(seqnames=Rle(chroms2use), ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
	if (keep.duplicate.reads) {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairsFromBam(file, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq'))
		} else {
			data.raw <- GenomicAlignments::readGAlignmentsFromBam(file, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq'))
		}
	} else {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairsFromBam(file, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
		} else {
			data.raw <- GenomicAlignments::readGAlignmentsFromBam(file, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
		}
	}
	if (pairedEndReads) {
		data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
		data.last <- as(GenomicAlignments::last(data.raw), 'GRanges')
		strand(data.last) <- strand(data.first)
		data <- sort(c(data.first, data.last))
	} else {
		data <- as(data.raw, 'GRanges')
	}
	return(data)

}


