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

runBreakpointr <- function(input.data, dataDirectory='./BreakPointR_analysis/', pairedEndReads=FALSE, chromosomes=NULL, windowsize=100, scaleWindowSize=T, trim=10, peakTh=0.33, zlim=3.291, bg=0.02, minReads=10, writeBed=T, verbose=T, depthTable=T) {

	fileDestination <- dataDirectory
	if (!file.exists(dataDirectory)) {
		dir.create(dataDirectory)
	}

	SummaryfileDestination <- file.path(dataDirectory, 'Summary_files/')
	if (!file.exists(SummaryfileDestination)) {
		dir.create(SummaryfileDestination)
	}

	## check the class of the input data
	if ( class(input.data) != "GRanges" ) {
		fragments <- bam2GRanges(input.data, pairedEndReads=pairedEndReads, chromosomes=chromosomes, keep.duplicate.reads=FALSE)
	} else {
		fragments <- input.data
		input.data <- 'CompositeFile'
	}

	filename <- basename(input.data)
	
	if (scaleWindowSize==F) {
		reads.per.window <- windowsize
		message("Calculating deltaWs")
		dw <- deltaWCalculator(fragments, reads.per.window=reads.per.window)
	}
	
	deltas.all.chroms <- GenomicRanges::GRangesList()
	breaks.all.chroms <- GenomicRanges::GRangesList()
	counts.all.chroms <- GenomicRanges::GRangesList()
	for (chr in unique(seqnames(fragments))) {

		message("Working on chromosome ",chr)

		message("Calculating deltaWs")
		if (scaleWindowSize==T){
			## normalize only for size of the chromosome 1
			reads.per.window <- round(windowsize/(seqlengths(fragments)[1]/seqlengths(fragments)[chr])) # scales the bin to chr size, anchored to chr1 (249250621 bp) 
			dw <- suppressWarnings( deltaWCalculator(fragments[seqnames(fragments)==chr], reads.per.window=reads.per.window) )
			
			## normalize for size of chromosome one and for read counts of each chromosome
			#num.reads <- length( fragments[seqnames(fragments)==chr] )
			#reads.per.window <- round( num.reads/( windowsize*(seqlengths(fragments)[chr]/seqlengths(fragments)[1]) ) )
			#reads.per.window <- round(reads.per.window/2)
			#dw <- suppressWarnings( deltaWCalculator(fragments[seqnames(fragments)==chr], reads.per.window=reads.per.window) )
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

		### count reads between the breaks and write into GRanges
		if (is.null(newBreaks)) {
			Ws <- length(deltaWs[strand(deltaWs) == '-'])
			Cs <- length(deltaWs[strand(deltaWs) == '+'])
			start.pos <- start(deltaWs)[1]
			end.pos <- end(deltaWs)[length(end(deltaWs))]
			chrRange <- GenomicRanges::GRanges(seqnames=chr, ranges=IRanges(start=start.pos, end=end.pos))	
			counts <- cbind(Ws,Cs)
			mcols(chrRange) <- counts
			seqlengths(chrRange) <- seqlengths(deltaWs)[chr]

			## genotype entire chromosome
			WC.ratio <- (chrRange$Ws-chrRange$Cs)/sum(c(chrRange$Ws,chrRange$Cs))
			if (WC.ratio > 0.8) {
				state <- 'ww'
			} else if (WC.ratio < -0.8) {
				state <- 'cc'
			} else if (WC.ratio < 0.2 & WC.ratio > -0.2) {
				state <- 'wc'
			} else {
				state <- '?'
			}					

			chrRange$states <- state
			suppressWarnings( counts.all.chroms[[chr]] <- chrRange )
	
		} else {
			
			breaks.strand <- newBreaks
			strand(breaks.strand) <- '*'
			breakrange <- gaps(breaks.strand)
			breakrange <- breakrange[strand(breakrange)=='*']

			## pull out reads of each line, and genotype in the fragments
			strand(breakrange) <- '-'
			Ws <- GenomicRanges::countOverlaps(breakrange, fragments)
			strand(breakrange) <- '+'
			Cs <- GenomicRanges::countOverlaps(breakrange, fragments)
			strand(breakrange) <- '*'

			## pull out states for each region between the breakpoints
			concat.states <- paste(newBreaks$genoT, collapse="-")
			split.states <- unlist(strsplit(concat.states, "-"))
			states.idx <- seq(from=1,to=length(split.states), by=2)
			states.idx <- c(states.idx, length(split.states))
			states <- split.states[states.idx]

			counts <- cbind(Ws,Cs)
			mcols(breakrange) <- counts
			breakrange$states <- states
			suppressWarnings( counts.all.chroms[[chr]] <- breakrange )
		}

		### write breaks and deltas into GRanges
		suppressWarnings( deltas.all.chroms[[chr]] <- deltaWs[,'deltaW'] )  #select only deltaW metadata column to store
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
			outputFile <- insertchr(newBreaks)
			depths <- round(outputFile$readNo / (end(outputFile)-start(outputFile))*1000000, digits=2)      
			depthT<- cbind(as.data.frame(outputFile), depths) # depths = chr, start, end, readNo, Ws, Cs, genoT, pVal
			if( chr == 'chr1'){ # write column names
				write.table(depthT, file=file.path(fileDestination, paste0(basename(input.data), '_depthTable.txt')), row.names=FALSE, col.names=T, quote=FALSE, append=F)
			} else{
				write.table(depthT, file=file.path(fileDestination, paste0(basename(input.data), '_depthTable.txt')), row.names=FALSE, col.names=F, quote=FALSE, append=T)
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
			head<- paste('track name=', basename(input.data), '_reads_', chr, ' visibility=1 colorByStrand="103,139,139 243,165,97"', sep="")
			write.table(head, file=paste(chrfileDestination, basename(input.data), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')   
			write.table(bedfile, file=paste(chrfileDestination, basename(input.data), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')
			
			## append the bedgraph of deltaWs:
			bedG<- as.data.frame(deltaWs)[c('chromosome','start','end','deltaW')]
			head<- paste('track type=bedGraph name=', basename(input.data),'_dWs', '_bin', reads.per.window, '_', chr, ' description=BedGraph_of_deltaWs_',basename(input.data), '_',  chr, ' visibility=full color=200,100,10', sep="")
			write.table(head, file=paste(chrfileDestination, basename(input.data), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')   
			write.table(bedG, file=paste(chrfileDestination, basename(input.data), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')      
			
			## append the breakpoints:
			if (is.null(newBreaks)) {
				bpG<- cbind(chr,0,1,'na')
			} else {
				newBreaks <- insertchr(newBreaks)
				bpG<-  as.data.frame(newBreaks)[c('chromosome','start','end','genoT')]
			}  
			head<- paste('track name=', basename(input.data), '_BPs', '_bin', reads.per.window, '_', chr, ' description=BedGraph_of_breakpoints_',basename(input.data), '_',  chr, ' visibility=pack color=75,125,180', sep="")
			write.table(head, file=paste(chrfileDestination, basename(input.data), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')   
			write.table(bpG, file=paste(chrfileDestination, basename(input.data), '_', chr, '_bin', reads.per.window, '_breakpointR.bed', sep=""), row.names=FALSE, col.names=F, quote=F, append=T, sep='\t')

		}
		
		### Write an INDEX file of deltaWs and breaks for ALL CHR
		if( chr== 'chr1'){ # write header
			head_reads<- paste('track name=', basename(input.data), '_reads visibility=1 colorByStrand="103,139,139 243,165,97"', sep="")
			write.table(head_reads, file=paste(SummaryfileDestination, basename(input.data), '_reads.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=F, sep='\t')   
			head_dW<- paste('track type=bedGraph name=', basename(input.data),'_dWs description=BedGraph_of_deltaWs_',basename(input.data), '_allChr visibility=full color=200,100,10', sep="")
			write.table(head_dW, file=paste(SummaryfileDestination, basename(input.data), '_DeltaWs.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=F, sep='\t')   
			head_br<- paste('track name=', basename(input.data), '_BPs description=BedGraph_of_breakpoints_',basename(input.data),'_allChr visibility=dense color=75,125,180', sep="")
			write.table(head_br,  file=paste(SummaryfileDestination, basename(input.data), '_breakPoints.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=F, sep='\t')
		} # write data
		write.table(bedfile,file=paste(SummaryfileDestination, basename(input.data), '_reads.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=T, sep='\t')   
		write.table(bedG, file=paste(SummaryfileDestination, basename(input.data), '_DeltaWs.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=T, sep='\t')
		write.table(bpG, file=paste(SummaryfileDestination, basename(input.data), '_breakPoints.txt', sep=""), row.names=FALSE, col.names=F, quote=FALSE, append=T, sep='\t')
		
	}
	## creating list of list where filename is first level list ID and deltas, breaks and counts are second list IDs
	deltas.all.chroms <- unlist(deltas.all.chroms)
	breaks.all.chroms <- unlist(breaks.all.chroms)
	counts.all.chroms <- unlist(counts.all.chroms)
	names(deltas.all.chroms) <- NULL
	names(breaks.all.chroms) <- NULL
	names(counts.all.chroms) <- NULL
	deltas.list <- GenomicRanges::GRangesList()
	breaks.list <- GenomicRanges::GRangesList()
	counts.list <- GenomicRanges::GRangesList()
	deltas.list[[filename]] <- deltas.all.chroms
	breaks.list[[filename]] <- breaks.all.chroms
	counts.list[[filename]] <- counts.all.chroms
	
	## save set parameters for future reference
	parameters <- c(windowsize=windowsize, scaleWindowSize=scaleWindowSize, trim=trim, peakTh=peakTh, zlim=zlim, bg=bg, minReads=minReads)

	return(list(deltas=deltas.list, breaks=breaks.list, counts=counts.list, params=parameters))
}

