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

runBreakpointr <- function(bamfile, dataDirectory='./BreakPointR_analysis/', pairedEndReads=FALSE, chromosomes=NULL, windowsize=100, scaleWindowSize=T, trim=10, peakTh=0.33, zlim=3.291, bg=0.02, minReads=10, writeBed=T, verbose=T, depthTable=T) {

	fileDestination <- dataDirectory
	if (!file.exists(dataDirectory)) {
		dir.create(dataDirectory)
	}

	SummaryfileDestination <- file.path(dataDirectory, 'Summary_files/')
	if (!file.exists(SummaryfileDestination)) {
		dir.create(SummaryfileDestination)
	}

	fragments <- bam2GRanges(bamfile, pairedEndReads=pairedEndReads, chromosomes=chromosomes, keep.duplicate.reads=FALSE)
	if (scaleWindowSize==F) {
		reads.per.window <- windowsize
		message("Calculating deltaWs")
		dw <- deltaWCalculator(fragments, reads.per.window=reads.per.window)
	}
	
	deltas.all.chroms <- GenomicRanges::GRangesList()
	breaks.all.chroms <- GenomicRanges::GRangesList()
	for (chr in unique(seqnames(fragments))) {

		message("Working on chromosome ",chr)

		message("Calculating deltaWs")
		if (scaleWindowSize==T){
			reads.per.window <- round(windowsize/(seqlengths(fragments)[1]/seqlengths(fragments)[chr])) # scales the bin to chr size, anchored to chr1 (249250621 bp)
			dw <- suppressWarnings( deltaWCalculator(fragments[seqnames(fragments)==chr], reads.per.window=reads.per.window) )
		}
		deltaWs <- dw[seqnames(dw)==chr]

		message("Finding breakpoints")
		breaks<- breakSeekr(deltaWs, trim=trim, peakTh=peakTh, zlim=zlim)

		message("Genotyping")
		if (length(breaks) >0 ) {
			## genotype the file between the breaks
			newBreaks <- GenotypeBreaks(breaks, fragments, backG=bg, minReads=minReads)
		} else {
			newBreaks<-NULL # assigns something to newBreaks if no peaks found
		}

		### write breaks and deltas into GRanges
		suppressWarnings( deltas.all.chroms[[chr]] <- deltaWs )
		if (length(newBreaks)) {
			suppressWarnings( breaks.all.chroms[[chr]] <- newBreaks )
		}
		
		### Write to BED ###
		chrfileDestination <- file.path(fileDestination, chr)
		insertchr <- function(hmm.gr) {
			mask <- which(!grepl('chr', seqnames(hmm.gr)))
			mcols(hmm.gr)$chromosome <- as.character(seqnames(hmm.gr))
			mcols(hmm.gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(hmm.gr)$chromosome[mask])
			mcols(hmm.gr)$chromosome <- as.factor(mcols(hmm.gr)$chromosome)
			return(hmm.gr)
		}
		if (depthTable==T && !is.null(newBreaks) && length(newBreaks)>0) # write a table with the depth information for everything between the breaks
		{
			message("Writing depthTable")
			outputFile<- GenotypeBreaks(newBreaks, fragments, backG=bg, minReads=minReads)
			outputFile <- insertchr(outputFile)
			depths <- round(outputFile$readNo / (end(outputFile)-start(outputFile))*1000000, digits=2)      
			depthT<- cbind(as.data.frame(outputFile), depths) # depths = chr, start, end, readNo, Ws, Cs, genoT, pVal
			if( chr == 'chr1'){ # write column names
				write.table(depthT, file=file.path(fileDestination, paste0(basename(bamfile), '_depthTable.txt')), row.names=FALSE, col.names=T, quote=FALSE, append=F)
			} else{
				write.table(depthT, file=file.path(fileDestination, paste0(basename(bamfile), '_depthTable.txt')), row.names=FALSE, col.names=F, quote=FALSE, append=T)
			}
		}
		if (writeBed==T) {
			## WRITE ALL THE DATA INTO A SINGLE BED FILE:
			# write table of fileFreqs
			deltaWs <- insertchr(deltaWs)
			bedfile<- as.data.frame(deltaWs)[c('chromosome','start','end','strand','preads')]
			bedfile$name <- 0
			bedfile <- bedfile[c('chromosome','start','end','name','preads','strand')]
			names(bedfile) <- c('chromosome','start','end','name','score','strand')
			head<- paste('track name=', basename(bamfile), '_reads_', chr, ' visibility=1 colorByStrand="103,139,139 243,165,97"', sep="")
			write.table(head, file=paste(chrfileDestination, basename(bamfile), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')   
			write.table(bedfile, file=paste(chrfileDestination, basename(bamfile), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')
			
			## append the bedgraph of deltaWs:
			bedG<- as.data.frame(deltaWs)[c('chromosome','start','end','deltaW')]
			head<- paste('track type=bedGraph name=', basename(bamfile),'_dWs', '_bin', reads.per.window, '_', chr, ' description=BedGraph_of_deltaWs_',basename(bamfile), '_',  chr, ' visibility=full color=200,100,10', sep="")
			write.table(head, file=paste(chrfileDestination, basename(bamfile), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')   
			write.table(bedG, file=paste(chrfileDestination, basename(bamfile), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')      
			
			## append the breakpoints:
			if (is.null(newBreaks)) {
				bpG<- cbind(chr,0,1,'na')
			} else {
				newBreaks <- insertchr(newBreaks)
				bpG<-  as.data.frame(newBreaks)[c('chromosome','start','end','genoT')]
			}  
			head<- paste('track name=', basename(bamfile), '_BPs', '_bin', reads.per.window, '_', chr, ' description=BedGraph_of_breakpoints_',basename(bamfile), '_',  chr, ' visibility=pack color=75,125,180', sep="")
			write.table(head, file=paste(chrfileDestination, basename(bamfile), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')   
			write.table(bpG, file=paste(chrfileDestination, basename(bamfile), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')

		}
		
		### Write an INDEX file of deltaWs and breaks for ALL CHR
		if( chr== 'chr1'){ # write header
			head_reads<- paste('track name=', basename(bamfile), '_reads visibility=1 colorByStrand="103,139,139 243,165,97"', sep="")
			write.table(head_reads, file=paste(SummaryfileDestination, basename(bamfile), '_reads.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=F, sep='\t')   
			head_dW<- paste('track type=bedGraph name=', basename(bamfile),'_dWs description=BedGraph_of_deltaWs_',basename(bamfile), '_allChr visibility=full color=200,100,10', sep="")
			write.table(head_dW, file=paste(SummaryfileDestination, basename(bamfile), '_DeltaWs.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=F, sep='\t')   
			head_br<- paste('track name=', basename(bamfile), '_BPs description=BedGraph_of_breakpoints_',basename(bamfile),'_allChr visibility=dense color=75,125,180', sep="")
			write.table(head_br,  file=paste(SummaryfileDestination, basename(bamfile), '_breakPoints.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=F, sep='\t')
		} # write data
		write.table(bedfile,file=paste(SummaryfileDestination, basename(bamfile), '_reads.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=T, sep='\t')   
		write.table(bedG, file=paste(SummaryfileDestination, basename(bamfile), '_DeltaWs.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=T, sep='\t')
		write.table(bpG, file=paste(SummaryfileDestination, basename(bamfile), '_breakPoints.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=T, sep='\t')
		
	}
	return(list(deltas=deltas.all.chroms, breaks=breaks.all.chroms))
}

