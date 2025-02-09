# Load libraries
library(Biostrings)

# Snakemake variables
dna_fa <- snakemake@input[["dna"]]
rna_fa <- snakemake@output[["rna"]]

# Convert DNA fasta to RNA fasta
dna <- readDNAStringSet(dna_fa)
rna <- RNAStringSet(chartr("T", "U", as.character(dna)))

# Modify headers
names(rna) <- gsub("^(\\S+).*", "\\1", names(rna))

writeXStringSet(rna,rna_fa)