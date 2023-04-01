#!/usr/bin/env Rscript

# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
# Load libraries
library(biomaRt)
library(dplyr)
library(optparse)

# Set up command-line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input file path"),
  make_option(c("-o", "--output"), type="character", help="Output file path")
)

# Parse command-line options
parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)
combined_graph_nodes <- read.csv(args$input, header = F)$V1
combined_graph_nodes <- read.csv('./Temporary_files/combined_graph_nodes.lst', header = F)$V1


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

paralogs <- getBM(
  attributes = c("ensembl_gene_id",
                 "external_gene_name",
                 "hsapiens_paralog_ensembl_gene",
                 "hsapiens_paralog_associated_gene_name"),
  filters = "external_gene_name",
  values = combined_graph_nodes,
  mart = ensembl
)
paralogs %>%
  filter(hsapiens_paralog_associated_gene_name != '') %>%
  select(external_gene_name, hsapiens_paralog_associated_gene_name) %>%
  write.table(args$output, sep='\t', quote = F, row.names = F)

tmp <- paralogs %>%
  filter(hsapiens_paralog_associated_gene_name != '') %>%
  select(external_gene_name, hsapiens_paralog_associated_gene_name)

n_tmp <- tmp %>%
    group_by(external_gene_name) %>%
    summarise(n=n()) %>%
    pull(n)

hist(n_tmp)