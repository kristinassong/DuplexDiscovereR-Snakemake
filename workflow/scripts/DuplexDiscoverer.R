# Load libraries
library(BiocManager)
library(DuplexDiscovereR)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

# Snakemake variables
sj_file <- snakemake@params[["sj"]]
gtf_file <- snakemake@input[["gtf"]]
chim_file <- snakemake@params[["chim"]]
gc_file <- snakemake@params[["gc"]]
fa_file <- snakemake@input[["fa"]]
dg_file <- snakemake@output[["dg"]]
sample <- snakemake@params[["sample"]]
lib <- snakemake@params[["lib"]]
table <- snakemake@params[["table"]]

# RNA fasta
genome_fasta <- c(fa_file)

# Gene counts
gene_counts <- c(gc_file)
GeneCounts <- read.delim(gene_counts, header = FALSE, stringsAsFactors = FALSE, skip=4)
GeneCounts <- GeneCounts[, c("V1", "V2")] 

# STAR Chimeric.out.junction
chim_out_junc <- c(chim_file)
ChimOutJunc <- read.delim(chim_out_junc, header = TRUE, stringsAsFactors = FALSE)

# Annotation gtf
genome_gtf <- c(gtf_file)
gtf = import.gff(genome_gtf)
GeneAnnoGR = gtf[gtf$type=='gene']
colnames(mcols(GeneAnnoGR))[colnames(mcols(GeneAnnoGR)) == "gene_biotype"] <- "gene_type"

# Splice junction
sj_out_tab <- c(sj_file)
SpliceJnc <- read.delim(sj_out_tab, header = FALSE, stringsAsFactors = FALSE)
colnames(SpliceJnc) <- c('chrom','start','end','strand','intron_motif','annot','uniq_mapped_reads','multi_mapped_reads','max_overhang')
SpliceJnc$strand[SpliceJnc$strand == 0] <- "*"
SpliceJnc$strand[SpliceJnc$strand == 1] <- "+"
SpliceJnc$strand[SpliceJnc$strand == 2] <- "-"
SpliceJncGR <- GRanges(seqnames = SpliceJnc$chrom, ranges = IRanges(start = SpliceJnc$start,
                        end = SpliceJnc$end), strand = SpliceJnc$strand)

# Run workflow
result <- DuplexDiscovereR::runDuplexDiscoverer(
  data = ChimOutJunc,
  junctions_gr = SpliceJncGR,
  anno_gr = GeneAnnoGR,
  df_counts = GeneCounts,
  sample_name = sample,
  lib_type = lib,
  table_type = table,
  fafile = genome_fasta,
)

# Save result
dg <- dd_get_duplex_groups(result)
dg_df <- as.data.frame(dg)
write.table(dg_df, file = dg_file, row.names = FALSE, sep = "\t")
