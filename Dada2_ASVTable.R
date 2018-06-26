################################################################################
# Title: Dada2 ASVTable
# Author: Nicholas Bense
# Date: 6/26/18
# Description: R script to create Amplicon Sequence Variant Tables for 16s and ITS 
# microbial diversity surveys. These can be used for taxonomic alignment downstream.
################################################################################

# Import packages
library(dada2)
library(argparse)

# Initialize argument parser
parser <- ArgumentParser()

# Read arguments/parameters
parser$add_argument("-f","--fwd",type="character",help="Name/path to fastq or fastq.gz files for single end sequencing or forward reads if paired end sequencing")
parser$add_argument("-r","--rev",type="character",default=NULL,help="Name/path to fastq or fastq.gz files for reverse reads if paired end sequencing")
parser$add_argument("-o","--out",type="logical",default="Dada2_ASVTable.tsv",help="Name/path to ASV table file output")
parser$add_argument("--pyr",type="logical",default=FALSE,help="Specify whether or not the data is from 454 pyrosequencing")

# Parse arguments
args <- parser$parse_args()

# Check to see if paired end read was passed to -r or --rev argument
if (!is.null(args$rev)){
  
  # Extract sample names from filepath/name
  sample.names<-sapply(strsplit(basename(args$fwd), "[.]"), `[`, 1)
  sample.names<-sapply(sample.names,function(x) gsub("_R1.",".",as.character(x)))
  
  # Learn error model for forward and reverse reads
  err_fwd <- learnErrors(args$fwd, multithread=TRUE)
  err_rvs <- learnErrors(args$rev, multithread=TRUE)
  
  # Dereplicate forward and reverse reads
  derep_Fs <- derepFastq(args$fwd, verbose=TRUE)
  derep_Rs <- derepFastq(args$rev, verbose=TRUE)
  
  # Assign sample names to dereplicated reads
  names(derep_Fs) <- sample.names
  names(derep_Rs) <- sample.names
  
  # Perform sample inference from reads
  dada_Fs <- dada(derep_Fs, err=err_Fs, multithread=TRUE)
  dada_Rs <- dada(derep_Rs, err=err_Rs, multithread=TRUE)
  
  # Merge paired end reads
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  
  # Create amplicon sequence variant table
  seqtab <- makeSequenceTable(mergers)
  
} else {
  
  # Extract sample names from filepath/name
  sample.names<-sapply(strsplit(basename(args$fwd), "[.]"), `[`, 1)

  # Learn error model for single end reads
  err_fns <- learnErrors(args$fwd, multithread=FALSE)

  # Dereplicate single end reads
  derep_fns <- derepFastq(args$fwd, verbose=TRUE)

  # Assign sample names to dereplicated reads
  names(derep_fns) <- sample.names

  # Check if --pyr tag was set to TRUE to represent pyrosequencing data
  if (args$pyr==TRUE)
    {
    # Perform sample inference from reads with pyrosequencing parameters
    dada_se <- dada(derep_se, err=err_se, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
  
    }
  else{
    
    # Perform sample inference from reads
    dada_se <- dada(derep_se, err=err_se, multithread=TRUE)

    # Create amplicon sequence variant table
    seqtab <- makeSequenceTable(dada_se)
  }
}

# Remove chimeras from asv tables
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Write ASV table to user specified output file/path
write.table(seqtab.nochim, file = args$o, sep = "\t", row.names=TRUE, col.names=TRUE)
