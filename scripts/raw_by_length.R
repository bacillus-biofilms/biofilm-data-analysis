# read necessary libraries. Please make sure they are all previously installed.

library(data.table)
library(plyr)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("file_raw", help="Path to input file containing raw counts obtained from Deseq.")
parser$add_argument("file_pid_lt_gl", help="Path to input file containing protein_id, locus_tag and gene length obtained from feature table file from GenBank.")
parser$add_argument("file_pid_ps", help="Path to input file containing protein_id and ps obtained from phylostratigraphy.")
parser$add_argument("file_out", help="Path to output file containing locus_tag and raw counts divided by gene length.")
args <- parser$parse_args()
file_raw <- args$file_raw
file_pid_lt_gl <- args$file_pid_lt_gl
file_pid_ps <- args$file_pid_ps
file_out <- args$file_out

if( file.access(file_raw) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", file_raw))
}

if( file.access(file_pid_lt_gl) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", file_ps))
}

if( file.access(file_pid_ps) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", file_len))
}

# read raw counts obtained from Deseq

biofilm_raw_counts <- read.table(file=file_raw, header=T, sep = '\t')


# make row names as first column

setDT(biofilm_raw_counts, keep.rownames = TRUE)[]


# give names to columns

colnames(biofilm_raw_counts) <- c("locus_tag", "6h_1", "2d_1", "2d_2", "2d_3", "3d_1", "3d_2", "3d_3", "5d_1", "5d_2", "5d_3", "7d_1", "6h_2", "7d_2", "7d_3", "14d_1", "14d_2", "14d_3", "1m_1", "1m_2", "1m_3", "2m", "6h_3", "LC_1", "LC_2", "LC_3", "12h_1", "12h_2", "12h_3", "1d_1", "1d_2", "1d_3")


# change the order of columns

biofilm_raw_counts <- biofilm_raw_counts[,c(1,24:26,2,13,23,27:29,30:32,3:12,14:22)]


# remove genes with all zeros

biofilm_raw_counts <- biofilm_raw_counts[rowSums(biofilm_raw_counts[, -1])>0, ]


# read table containing protein_id, locus_tag and gene_length obtained from feature table file from GenBank

pid_lt_gl <- read.table(file=file_pid_lt_gl, header = T, sep = '\t')


# join raw counts with gene length by locus_tag

biofilm_raw_counts_gl <- join (pid_lt_gl, biofilm_raw_counts, by="locus_tag", type="inner")


# read table containing protein_id and ps obtained from phylostratigraphy

pid_ps <- read.table(file=file_pid_ps, header=T, sep = '\t')


# join raw counts and gene lenght information with ps

biofilm_raw_counts_gl_ps <- join (biofilm_raw_counts_gl, pid_ps, by="protein_id", type="inner")


# divide raw counts by gene length

biofilm_raw_counts_gl_ps_by_gene_length <- sweep(biofilm_raw_counts_gl_ps[,4:34], biofilm_raw_counts_gl_ps[,"gene_length"], MARGIN=1,"/")


# merge values with gene ids

biofilm_raw_counts_gl_ps_by_gene_length_id <- cbind (biofilm_raw_counts_gl_ps_by_gene_length, biofilm_raw_counts_gl_ps)


# remove unnecessary columns

biofilm_raw_counts_gl_ps_by_gene_length_id <- biofilm_raw_counts_gl_ps_by_gene_length_id[,c(33,1:31)]


# save output file containing raw locus_tag and raw counts divided by gene_length

write.table(biofilm_raw_counts_gl_ps_by_gene_length_id, file_out, sep="\t", quote=F, row.names=F)
