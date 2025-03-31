# Install necessary packages if you haven't already

packages <- c("BiocManager", "Biostrings", "rtracklayer", "GenomicRanges", "derfinder")

# Check if each package is installed, install if missing
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "BiocManager") {
      install.packages("BiocManager")
    } else {
      BiocManager::install(pkg)
    }
  }
}

# Load libraries
library(Biostrings)
library(rtracklayer)
library(GenomicRanges)
library(derfinder)

# Function to calculate GC content in a sequence for a given bin size
calculate_gc_content_bins <- function(sequence, bin_size) {
  # Convert the DNA sequence to uppercase and ensure it's in a character vector form
  sequence <- toupper(as.character(sequence))
  
  # Get the total length of the sequence
  seq_length <- nchar(sequence)
  
  # Initialize an empty vector to store GC content values
  gc_contents <- NULL
  
  # Loop through the sequence in windows of size 'bin_size'
  for (start_pos in seq(1, seq_length, by = bin_size)) {
    end_pos <- min(start_pos + bin_size - 1, seq_length)
    subseq <- substr(sequence, start_pos, end_pos)
    
        # Count occurrences of "GC" dinucleotide
    gc_count <- gregexpr("GC", subseq, fixed = TRUE)
    gc_matches <- 0
    gc_matches <- sum(gc_count[[1]] > 0)  # Count valid matches

    gc_contents <- c(gc_contents, gc_matches)
  }
  
  return(gc_contents)
}

# Function to create a BigWig file from GC content for multiple chromosomes/contigs
create_bigwig <- function(fasta_file, bin_size, output_file) {
  # Read the FASTA file
  seqs <- readDNAStringSet(fasta_file)
  
  # Initialize an empty GRanges list to store all GC content data
  gr_sum <- list()
  
  # Loop through each sequence (chromosome/contig) in the FASTA file
  for (chrom_name in names(seqs)) {
    # Calculate GC content in bins for the current sequence
    gc_content_bins <- calculate_gc_content_bins(seqs[[chrom_name]], bin_size)
    
    # Create a GRanges object for the current sequence
    gr <- GRanges(
      seqnames = rep(chrom_name, length(gc_content_bins)),
      ranges = IRanges(start = seq(1, length(gc_content_bins)) * bin_size - bin_size + 1,
                       end = seq(1, length(gc_content_bins)) * bin_size),
      score = gc_content_bins
    )
    
    # Append the GRanges object to the list
    gr_sum[[ chrom_name ]] = gr
  }

  print(gr_sum )

  #gr_total = derfinder::coerceGR( gr_sum )
  #gr_total <- do.call(c, gr_sum)
  gr_total <- GRanges(
    seqnames = names(unlist(lapply(gr_sum, seqnames))),
    ranges = IRanges(
      start = unlist(lapply(gr_sum, start)),
      end = unlist(lapply(gr_sum, end))
    ),
    score = unlist(lapply(gr_sum, function(gr) mcols(gr)$score))
  )
  max_lengths <- tapply(end(gr_total), seqnames(gr_total), max)
  seqinfo(gr_total) <- Seqinfo(seqnames = names(max_lengths), seqlengths = as.vector(max_lengths))

  browser()
  # Export to BigWig
  export(gr_total, con = output_file, format = "BigWig")
  
  cat("BigWig file created at:", output_file, "\n")
}

# Main script
process_fasta_to_bigwig <- function(fasta_file, bin_size, output_file) {
  # Create the BigWig file from GC content for all sequences in the FASTA file
  create_bigwig(fasta_file, bin_size, output_file)
}


library(optparse)

# Set up command-line options
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", help = "Input FASTA file path"),
  make_option(c("-o", "--outfile"), type = "character", help = "Output BigWig file path"),
  make_option(c("-b", "--bin_size"), type = "integer", default = 50, help = "Bin size for GC content calculation [default 50]")
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opts$infile) || is.null(opts$outfile)) {
  print_help(opt_parser)
  stop("Both input file and output file must be provided")
}

# Process the FASTA file and generate the BigWig file
process_fasta_to_bigwig(opts$infile, opts$bin_size, opts$outfile)

## to test this
# source ('tests/data/Same_in_R.R')
# process_fasta_to_bigwig( "tests/data/chrM.fa.gz", 50,  "tests/data/out.R.bw")