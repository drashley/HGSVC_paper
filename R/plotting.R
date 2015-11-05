#' Plotting function for binned read counts

#' @param file Bamfile with aligned reads.
#' @param bamindex Bam-index file with or without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.

#' @import ggplot2
#' @import reshape2
#' @import grid
#' @importFrom Aneufinder indexBam scanBamHeader ScanBamParam scanBamFlag
