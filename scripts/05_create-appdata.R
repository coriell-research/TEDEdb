#!/usr/bin/env Rscript
# Create combined metadata SummarizedExperiment objects
#
# This script will import all SE files from each experiment and extract the DE
# results into a single SummarizedExperiment object used as the backend database
# for the Shiny app. Metadata information, annotated manually in Excel, is read
# in and aligned with columns of the imported DE data to create a finalized SE 
# object.
#
# ----------------------------------------------------------------------------
message("Loading libraries...")
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(HDF5Array))
source(here("scripts", "helpers.R"))
N_CORES <- 32


message("Getting all processed SummarizedExperiment files...")
se_files <- list.files(
  path = here("bioprojects"),
  pattern = "se.rds",
  recursive = TRUE,
  full.names = TRUE
)
names(se_files) <- str_extract(se_files, "PRJNA[0-9]+")
message("Found ", length(se_files), " results.")


# Extract DE results from SE objects and bind into single DT -----------------


message("Extracting DE data from SummarizedExperiments...")
de <- rbindlist(
  parallel::mclapply(se_files, extract_de, mc.cores = N_CORES), 
  idcol = "BioProject"
  )

# Add new column as unique identifier "id"
de[, id := str_c(BioProject, contrast, sep = ".")]


# Read in contrast-level metadata files --------------------------------------


message("Reading in contrast-level metadata...")
metafile <- here("appdata", "contrast-metadata - metadata.tsv")
drugfile <- here("appdata", "contrast-metadata - drugs.tsv")
cellfile <- here("appdata", "contrast-metadata - cells.tsv")

if (any(!file.exists(c(metafile, cellfile, drugfile)))) {
  message("One of: ", metafile, drugfile, " and ", cellfile, "do not exist!")
  stop()
}
metadata <- fread(metafile)
drugs <- fread(drugfile)
cells <- fread(cellfile)

# Which IDs are annotated in the metadata?
analyzed_ids <- de[, unique(id)]
annotated_ids <- metadata[, unique(id)]
missing_annotation <- setdiff(analyzed_ids, annotated_ids)
missing_data <- setdiff(annotated_ids, analyzed_ids)

if (length(missing_annotation) > 0) {
  message("There are ", length(missing_annotation), " IDs without an annotation!")
  message("Writing these IDs to: appdata/missing-annotations.tsv")
  fwrite(data.table(id = missing_annotation), here("appdata", "missing-annotations.tsv"), sep = "\t")
}

if (length(missing_data) > 0) {
  message("There are ", length(missing_data), " IDs annotated without data!")
  message("Writing these IDs to: appdata/missing-data.tsv")
  fwrite(data.table(id = missing_data), here("appdata", "missing-data.tsv"), sep = "\t")
}

message("Inner joining metadata onto differential expression data...")
metadata <- metadata[drugs, on = "drug", nomatch = 0L]
metadata <- metadata[cells, on = "cell_line", nomatch = 0L]

message("Flagging contrasts with potential outliers...")
outliers <- getOutlierContrasts(se_files, cores = N_CORES)
metadata[, `:=`(outlier_avg_frag_length_mean = id %chin% outliers$outlier_avg_frag_length_mean,
                outlier_avg_frag_length_sd = id %chin% outliers$outlier_avg_frag_length_sd,
                outlier_avg_num_eq_classes = id %chin% outliers$outlier_avg_num_eq_classes,
                outlier_avg_num_processed = id %chin% outliers$outlier_avg_num_processed,
                outlier_avg_num_mapped = id %chin% outliers$outlier_avg_num_mapped,
                outlier_avg_num_decoy_fragments = id %chin% outliers$outlier_avg_num_decoy_fragments,
                outlier_avg_num_dovetail_fragments = id %chin% outliers$outlier_avg_num_dovetail_fragments,
                outlier_avg_num_fragments_filtered_vm = id %chin% outliers$outlier_avg_num_fragments_filtered_vm,
                outlier_avg_num_alignments_below_threshold_vm = id %chin% outliers$outlier_avg_num_alignments_below_threshold_vm,
                outlier_avg_percent_mapped = id %chin% outliers$outlier_avg_percent_mapped)]

# Create a single catch-all flag
metadata[, outlier_flags := getOutlierFlags(.SD), .SDcols = names(metadata) %like% "outlier"]
metadata[outlier_flags == "", outlier_flags := "None"]

# Remove the intermediate columns
metadata[, `:=`(outlier_avg_frag_length_mean=NULL,
                outlier_avg_frag_length_sd=NULL,
                outlier_avg_num_eq_classes=NULL,
                outlier_avg_num_processed=NULL,
                outlier_avg_num_mapped=NULL,
                outlier_avg_num_decoy_fragments=NULL,
                outlier_avg_num_dovetail_fragments=NULL,
                outlier_avg_num_fragments_filtered_vm=NULL,
                outlier_avg_num_alignments_below_threshold_vm=NULL,
                outlier_avg_percent_mapped=NULL
              )]


# Shape into matrices --------------------------------------------------------


message("Creating assay data from differential expression results...")
assay_cols <- c("logFC", "AveExpr", "z", "P.Value", "adj.P.Val", "SE")
assays <- vector("list", length(assay_cols))

names(assays) <- assay_cols
for (i in seq_along(assay_cols)) {
  assays[[i]] <- col2assay(de, rows = "feature_id", cols = "id", vals = assay_cols[[i]])
}

# Impute missing values for assays used in dim reduction
assays[["logFC"]][is.na(assays[["logFC"]])] <- 0
assays[["z"]][is.na(assays[["z"]])] <- 0
assays[["P.Value"]][is.na(assays[["P.Value"]])] <- 1


# Create SummarizedExperiment object of DE results ---------------------------


setDF(metadata, rownames = metadata$id)
keep <- intersect(colnames(assays[[1]]), rownames(metadata))
metadata <- metadata[keep, ]
assays <- lapply(assays, \(x) x[, keep])
stopifnot("rownames of metadata do not match colnames of matrices!" = all(colnames(assays[[1]]) == rownames(metadata)))

message("Creating final SummarizedExperiment object from differential expression results...")
se <- SummarizedExperiment(assays = assays, colData = metadata)


# Create rowData -------------------------------------------------------------


# TODO: This hardcoded path needs a permanent home in a top-level directory
#  possibly add the analysis/ directory to the github repo as well
gene_lengths <- fread("/home/gennaro/data/cancer-rnaseq-database/analysis/results/data-files/01/v26-gene-length-annotation.tsv")
gene_lengths <- gene_lengths[gene_type == "protein_coding"]

rd <- data.table(feature_name = rownames(se))
rd <- merge(x = rd, y = gene_lengths, by.x = "feature_name", by.y="gene_name", all.x=TRUE)

# Some gene names are replicated and need to be collapsed (mostly double annotated on PAR Y)
rd <- rd[!gene_id %like% "PAR_Y"]

# Select the longest of the gene_ids to keep
multi <- rd[, .N, by=feature_name][N > 1][, feature_name]
all_ids <- rd[feature_name %chin% multi, unique(gene_id)]
keep_ids <- rd[feature_name %chin% multi][
               order(feature_name, -total_gene_length)][, 
               head(.SD, n=1), by=feature_name][, gene_id]
drop_ids <- setdiff(all_ids, keep_ids)
rd <- rd[!gene_id %chin% drop_ids]

# Add feature type labels and separate TE information
rd[, feature_type := fifelse(feature_name %like% ".*\\..*\\..*", "TE", "Gene")]
rd[feature_type == "TE", c("class", "family", "subfamily") := tstrsplit(feature_name, "\\.")]

setDF(rd, rownames = rd$feature_name)
rd$feature_name <- NULL

rd <- rd[rownames(se), ]
stopifnot("rownames(rowData) != rownames(se)" = all(rownames(rd) == rownames(se)))
rowData(se) <- cbind(rowData(se), rd)


# Save finalized object ------------------------------------------------------


message("Saving HDF5-backed SE for database backend...")
saveHDF5SummarizedExperiment(se, dir=here("appdata", "se_hdf5"), replace=TRUE)


# Select input data ----------------------------------------------------------


# Create data used by the shiny app for populating select inputs
select_inputs <- lapply(metadata, \(x) sort(unique(x)))
saveRDS(select_inputs, here("appdata", "select-inputs.rds"))
message("Appdata creation complete.")
message("Done.")
