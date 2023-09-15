#!/usr/bin/env Rscript

# bwa16S -- a reference based 16S analysis package
# By Jonathan Badger
# Version 0.95 September 11, 2023

library(optparse)
library(stringr)

# read FASTA file and return data frame
read_fasta <- function(file) {
  # Read the file line by line
  fasta <- readLines(file)
  # Identify header lines
  ind <- grep(">", fasta)
  # Identify the sequence lines
  s <- data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
  # Process sequence lines
  seqs <- rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse="")
  }
  # Create a data frame 
  DF <- data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)
  # Return the data frame as a result object from the function
  return(DF)
}

# run FLASh to merge forward and reverse reads
merge_reads <- function(input, threads) {
  flash <- Sys.which("flash")
  if (flash == "") {
    stop("No executable for FLASh found.")
  }
  r1 <- Sys.glob(str_c(input,"/*R1*fastq*"))
  r2 <- Sys.glob(str_c(input,"/*R2*fastq*"))
  if (length(r1) != length(r2) || length(r1)==0) {
    stop(str_c("Problem with fastqs in ", input))
  }
  dir.create("merged", showWarnings = FALSE)
  for(i in 1:length(r1)) {
    sample = strsplit(basename(r1[i]),"_")[[1]][1]
    merged <- str_c("merged/", sample, ".fastq.gz")
    if (!file.exists(merged)) {
      cmd <- str_c(flash, " ", r1, " ", r2, " -t ", threads, " -M 250 -z -c > ", merged)
      system(cmd)
    }
  }
}

# map merged reads against reference using bwa
map_reads <- function(input, reference, threads) {
  dir.create("mapped", showWarnings = FALSE)
  bwa <- Sys.which("bwa")
  if (bwa == "") {
    stop("No executable for bwa found.")
  }
  for(fq in Sys.glob(str_c(input, "/*.fastq.gz"))) {
    sample <- strsplit(basename(fq),".fastq")[[1]][1]
    sam <- str_c("mapped/", sample, ".sam")
    if (!file.exists(sam)) {
      cmd <- str_c(bwa, " mem -t ", threads, " -o ", sam, 
                   " ", reference, " ", fq)
      system(cmd)
    }
  }
}

# make table of counts from sam files
make_counts <- function(nerr, reference) {
  seqs <- read_fasta(normalizePath(reference))
  counts <- data.frame(row.names = seqs$name)
  samtools <- Sys.which("samtools")
  if (samtools == "") {
    stop("No executable for samtools found.")
  }
  for(sam in Sys.glob("./mapped/*.sam")) {
    sample <- strsplit(basename(sam),".sam")[[1]][1]
    cmd <- str_c(samtools, " view -e '[NM] <= ", nerr, "' ", sam, 
                 " | cut -f 3 | sort | uniq -c")
    ct <- system(cmd, intern = TRUE)
    for(c in ct) {
      fields <- strsplit(str_trim(c), " ")[[1]]
      counts[fields[2], sample] <- as.integer(fields[1])
    }
  }
  counts <- cbind(Taxon=row.names(counts), counts)
  counts[is.na(counts)] <- 0
  write.table(counts, file = "unit_table.tsv", sep = "\t", 
              quote = FALSE, row.names = FALSE)
}

# main function for parsing arguments and running pipeline
main <- function(args) {
    version <- "0.95 14 Sept 2023"
    if (length(args)==0) {
        args <- c("-h")
    }
    option_list <- list(
        make_option(c("-i", "--input"), action="store", default=NULL, 
                    type="character",
                    help="path of directory where fastq files are"),
	      make_option(c("-d", "--database"), default=NULL, action="store",
                    type="character",
                    help="database bwa index file with taxonomy"),
        make_option(c("-t", "--threads"), action="store", default=8, 
                    type="numeric",
                    help ="number of threads to use (default 8"),
        make_option(c("-n", "--nerr"), action="store",
                    type="numeric", default=2,
                    help="maximum number of errors allowed to reference (default 2)")
    )
    opt <- parse_args(OptionParser(option_list=option_list), args)
    
    if (is.null(opt$input)) {
        stop("option -i/--input is required\n")
    }
    if (is.null(opt$database)) {
        stop("option -d/--database is required\n")
    }
	
		merge_reads(opt$input, opt$threads)
		map_reads("merged/", opt$database, opt$threads)
		make_counts(opt$nerr, opt$database)
}

# run main function only if run as script
if (!interactive()) {
    main(commandArgs(trailingOnly = TRUE))
}


