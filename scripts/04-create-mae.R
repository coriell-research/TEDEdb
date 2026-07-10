#!/usr/bin/env Rscript
#
# Create a MultiAssayExperiment object a user supplied annotation file
# and the quant.sf files for each sample processed by the REdiscoverTE pipeline.
#
# NOTE: There are hard-coded paths to annotation files
# ------------------------------------------------------------------------------
suppressMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(MultiAssayExperiment, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tximport, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(jsonlite, warn.conflicts = FALSE, quietly = TRUE))

# get commandline arguments
option_list <- list(
  make_option(c("-d", "--quants_dir"),
              type = "character",
              default = NULL,
              help = "Path to quants/ directory containing sub-directories for each SAMPLE",
              metavar = "quants_dir"
  ),
  make_option(c("-f", "--annotation"),
              type = "character",
              default = NULL,
              help = "Path to annotation file. Must at least contain a column called 'Run' and 'group'. Other metadata is optional",
              metavar = "annotation"
  ),
  make_option(c("-o", "--out_dir"),
              type = "character",
              default = ".",
              help = "Location to save the final MultiAssayExperiment object RDS file",
              metavar = "out_dir"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

message("Reading in annotation files for processing steps...")
tx2gene <- readRDS("resources/REdiscoverTE_hg38/tx2gene_REdiscoverTE.rds")

message("Creating count matrices from quant files...")
quant_files <- list.files(
  path = opt$quants_dir,
  pattern = "quant.sf*",
  recursive = TRUE,
  full.names = TRUE
)
names(quant_files) <- regmatches(quant_files, regexpr("SRR[0-9]+", quant_files))

txi <- tximport(
  files = quant_files,
  type = "salmon",
  countsFromAbundance = "lengthScaledTPM",
  tx2gene = tx2gene,
  importer = function(x) data.table::fread(x, showProgress = FALSE)
)

# Extract individual matrices
intron_mat <- txi$counts[grepl("__intron$", rownames(txi$counts)), ]
intergenic_mat <- txi$counts[grepl("__intergenic", rownames(txi$counts)), ]
exon_mat <- txi$counts[grepl("__exon", rownames(txi$counts)), ]
all_re <- Reduce(union, list(rownames(intron_mat), rownames(intergenic_mat), rownames(exon_mat)))
gene_mat <- txi$counts[!rownames(txi$counts) %chin% all_re, ]

# Combine the intronic and intergenic counts into a single RE matrix
intron_dt <- as.data.table(intron_mat, keep.rownames = "feature_id")
intergenic_dt <- as.data.table(intergenic_mat, keep.rownames = "feature_id")
re_dt <- rbind(intron_dt, intergenic_dt)
re_dt[, feature_id := gsub("__intron$|__intergenic$", "", feature_id)]
re_dt.m <- melt(re_dt, id.vars = "feature_id", variable.name = "Run", value.name = "count")
by_re <- re_dt.m[, .(count = sum(count)), by = .(Run, feature_id)]
re_mat <- as.matrix(dcast(by_re, feature_id ~ Run, value.var = "count", fill = 0.0), rownames = "feature_id")

# Read in the sample metadata
message("Reading in sample annotation information...")
meta_df <- fread(opt$annotation)
stopifnot("Run" %in% colnames(meta_df))
stopifnot("group" %in% colnames(meta_df))

# Collect meta_info from Salmon for each run
message("Collecting mapping information from Salmon meta_info.json files...")
meta_files <- list.files(
  path = opt$quants_dir,
  pattern = "meta_info.json",
  recursive = TRUE,
  full.names = TRUE
)
names(meta_files) <- regmatches(meta_files, regexpr("SRR[0-9]+", meta_files))

# Helper function for extracting JSON metadata into a data.table
getJsonData <- function(x) {
  jdata <- jsonlite::read_json(x)
  d <- lapply(jdata, \(x) paste(unlist(x), collapse = "; "))
  as.data.table(d)
}

# Collect into a single data.table
meta_info <- rbindlist(lapply(meta_files, getJsonData), idcol="Run")
cols <- names(meta_info)[names(meta_info) %like% "num_|frag"]
meta_info[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
meta_info[, percent_mapped := as.numeric(percent_mapped)]

# Join sample annotations and meta_info
meta_df <- merge(meta_df, meta_info, by="Run", all.x=TRUE)

# Convert to data.frame for colData of MAE
setDF(meta_df, rownames = meta_df$Run)
meta_df <- subset(meta_df, select = -Run)
stopifnot(all(names(quant_files) %in% rownames(meta_df)))

# Reorder the metadata to match colnames of matrices
stopifnot(all(colnames(gene_mat) == colnames(re_mat)))
meta_df <- meta_df[colnames(gene_mat), ]

# Create MultiAssayExperiment object from each assay and metadata
exp_list <- list(
  "gene counts" = gene_mat,
  "RE counts" = re_mat,
  "exon RE counts" = exon_mat,
  "intron RE counts" = intron_mat,
  "intergenic RE counts" = intergenic_mat
)

message("Creating MultiAssayExperiment object...")
mae <- MultiAssayExperiment(experiments = exp_list, colData = meta_df)

# write mae .rds file to out_dir
message(paste("Writing MultiAssayExperiment object to", file.path(opt$out_dir, "MultiAssayExperiment.rds")))
saveRDS(mae, file = file.path(opt$out_dir, "MultiAssayExperiment.rds"))
message("Done.")
