#' Wrapper function for breakpointR
#'
#' This script will move through .bam files in a folder and perform several steps (see Details).
#'
#' 1. calculate deltaWs chromosome-by-chromsome
#' 2. localize breaks that pass zlim above the threshold
#' 3. genotype both sides of breaks to confirm whether strand state changes
#' 4. write a file of _reads, _deltaWs and _breaks in a chr fold -> can upload on to UCSC Genome browser
#' 5. write a file for each index with all chromosomes included -> can upload on to UCSC Genome browser
#'
#' @param bamfile A file with aligned reads in BAM format.
#' @param dataDirectory Output directory. If non-existent it will be created.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param chromosomes If only a subset of the chromosomes should be processed, specify them here.
#' @param windowsize The window size used to calculate deltaWs, either an integer or 'scale'.
#' @param trim The amount of outliers in deltaWs removed to calculate the stdev (10 will remove top 10\% and bottom 10\% of deltaWs).
#' @param peakTh The treshold that the peak deltaWs must pass to be considered a breakpoint.
#' @param zlim The number of stdev that the deltaW must pass the peakTh (ensures only significantly higher peaks are considered).
#' @param bg The amount of background introduced into the genotype test.
#' @param minReads The minimum number of reads required for genotyping.
#' @param writeBed If \code{TRUE}, will generate a bed of reads and deltaWs and breaks for upload onto the UCSC genome browser.
#' @param verbose Verbose messages if \code{TRUE}.
#' @param depthTable If \code{TRUE}, will generate a table that also contains reads/Mb between breaks (important for SCE calls).
#' @author Ashley Sanders, David Porubsky, Aaron Taudt
#' @export

runBreakpointrALL <- function(datapath=NULL, dataDirectory='BreakPointR_analysis', pairedEndReads=FALSE, chromosomes=NULL, windowsize=100, scaleWindowSize=T, trim=10, peakTh=0.33, zlim=3.291, bg=0.02, minReads=10, writeBed=T, verbose=T, depthTable=T) {

	files <- list.files(datapath, pattern=".bam$", full=T)

	if (!file.exists(dataDirectory)) {
		dir.create(dataDirectory)
	}

	deltas.all.files <- GenomicRanges::GRangesList()
	breaks.all.files <- GenomicRanges::GRangesList()
	for (bamfile in files) {
		message("Working on file ",bamfile)
		
		deltas.breaks.obj <- runBreakpointr(bamfile=bamfile, dataDirectory=file.path(dataDirectory, basename(bamfile)), pairedEndReads=pairedEndReads, chromosomes=chromosomes, windowsize=windowsize, trim=trim, peakTh=peakTh, zlim=zlim, bg=bg, minReads=minReads, writeBed=writeBed, verbose=verbose, depthTable=depthTable)

		deltas <- unlist(deltas.breaks.obj$deltas)
		breaks <- unlist(deltas.breaks.obj$breaks)

		suppressWarnings( deltas.all.files[[bamfile]] <- deltas )
		if (length(breaks)) {
			suppressWarnings( breaks.all.files[[bamfile]] <- breaks )
		}
	}
	return(list(deltas=deltas.all.files, breaks=breaks.all.files)) 
}

