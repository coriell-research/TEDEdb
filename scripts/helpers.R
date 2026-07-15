extract_de <- function(fpath) {
  se <- readRDS(fpath)
  metadata(se)[["de"]]
}


col2assay <- function(df, rows, cols, vals) {
  m <- data.table::dcast(df, get(rows) ~ get(cols), value.var = vals, fill = NA)
  as.matrix(m, rownames = "rows")
}


getSampleMetadata <- function(x) {
  se <- readRDS(x)
  result <- data.table::as.data.table(data.frame(colData(se)))

  return(result)
}


getContrasts <- function(x) {
  se <- readRDS(x)
  cm <- S4Vectors::metadata(se)$fit$contrasts
  cm <- data.table::as.data.table(cm, keep.rownames = "group")

  cm.m <- data.table::melt(
    cm,
    id.vars = "group",
    variable.name = "contrast",
    value.name = "is_member",
    variable.factor = FALSE
  )

  return(cm.m)
}


getOutlierContrasts <- function(sefiles, cores) {
  # Sample-level metadata with Salmon meta_info
  message("Reading in sample-level metadata...")
  metadata <- data.table::rbindlist(
    parallel::mclapply(sefiles, getSampleMetadata, mc.cores = cores),
    fill = TRUE,
    idcol = "BioProject"
  )
  metadata[, group_id := paste(BioProject, group, sep = ".")]
  keep <- c(
    "BioProject",
    "BioSample",
    "group",
    "group_id",
    "salmon_version",
    "library_types",
    "avg_frag_length_mean",
    "avg_frag_length_sd",
    "avg_num_eq_classes",
    "avg_num_processed",
    "avg_num_mapped",
    "avg_num_decoy_fragments",
    "avg_num_dovetail_fragments",
    "avg_num_fragments_filtered_vm",
    "avg_num_alignments_below_threshold_vm",
    "avg_percent_mapped"
  )
  metadata <- metadata[, ..keep]
  metadata[,
    library_layout := data.table::fifelse(
      library_types %like% "^[IOM]",
      "PAIRED",
      "SINGLE"
    )
  ]

  # Compute sample-level failures based on global information
  message("Computing outliers based on global distributions...")
  # fmt: skip
  metadata[, `:=`(is_outlier_avg_frag_length_mean = coriell::outliers_by_mad(avg_frag_length_mean),                                
                  is_outlier_avg_frag_length_sd = coriell::outliers_by_mad(avg_frag_length_sd),                                 
                  is_outlier_avg_num_eq_classes = coriell::outliers_by_mad(log10(avg_num_eq_classes+1), direction="low"),                               
                  is_outlier_avg_num_processed = coriell::outliers_by_mad(log10(avg_num_processed+1), direction="low"),                                
                  is_outlier_avg_num_mapped = coriell::outliers_by_mad(log10(avg_num_mapped+1), direction="low"),                                      
                  is_outlier_avg_num_decoy_fragments = coriell::outliers_by_mad(log10(avg_num_decoy_fragments+1), direction="high"),                     
                  is_outlier_avg_num_dovetail_fragments = coriell::outliers_by_mad(log10(avg_num_dovetail_fragments+1), direction="high"),                  
                  is_outlier_avg_num_fragments_filtered_vm = coriell::outliers_by_mad(log10(avg_num_fragments_filtered_vm+1), direction="high"),             
                  is_outlier_avg_num_alignments_below_threshold_vm = coriell::outliers_by_mad(log10(avg_num_alignments_below_threshold_vm+1), direction="high"),
                  is_outlier_avg_percent_mapped = coriell::outliers_by_mad(avg_percent_mapped, direction="low")
                  ), by = library_layout]

  # Get contrast to group mapping
  message("Reading in group to contrast mapping...")
  cm_dt <- data.table::rbindlist(
    parallel::mclapply(sefiles, getContrasts, mc.cores = cores),
    idcol = "BioProject"
  )
  cm_dt <- cm_dt[is_member != 0]
  cm_dt[, id := paste(BioProject, contrast, sep = ".")]
  cm_dt[, group_id := paste(BioProject, group, sep = ".")]

  # Collect groups with samples that failed QC for some reason
  message("Determining groups with outlier samples...")
  has_outlier_avg_frag_length_mean <- metadata[
    is_outlier_avg_frag_length_mean == TRUE,
    unique(group_id)
  ]
  has_outlier_avg_frag_length_sd <- metadata[
    is_outlier_avg_frag_length_sd == TRUE,
    unique(group_id)
  ]
  has_outlier_avg_num_eq_classes <- metadata[
    is_outlier_avg_num_eq_classes == TRUE,
    unique(group_id)
  ]
  has_outlier_avg_num_processed <- metadata[
    is_outlier_avg_num_processed == TRUE,
    unique(group_id)
  ]
  has_outlier_avg_num_mapped <- metadata[
    is_outlier_avg_num_mapped == TRUE,
    unique(group_id)
  ]
  has_outlier_avg_num_decoy_fragments <- metadata[
    is_outlier_avg_num_decoy_fragments == TRUE,
    unique(group_id)
  ]
  has_outlier_avg_num_dovetail_fragments <- metadata[
    is_outlier_avg_num_dovetail_fragments == TRUE,
    unique(group_id)
  ]
  has_outlier_avg_num_fragments_filtered_vm <- metadata[
    is_outlier_avg_num_fragments_filtered_vm == TRUE,
    unique(group_id)
  ]
  has_outlier_avg_num_alignments_below_threshold_vm <- metadata[
    is_outlier_avg_num_alignments_below_threshold_vm == TRUE,
    unique(group_id)
  ]
  has_outlier_avg_percent_mapped <- metadata[
    is_outlier_avg_percent_mapped == TRUE,
    unique(group_id)
  ]

  # Add indicators for each group with sample failures
  # fmt: skip
  cm_dt[, `:=`(outlier_avg_frag_length_mean = group_id %in% has_outlier_avg_frag_length_mean,
               outlier_avg_frag_length_sd = group_id %in% has_outlier_avg_frag_length_sd,
               outlier_avg_num_eq_classes = group_id %in% has_outlier_avg_num_eq_classes,
               outlier_avg_num_processed = group_id %in% has_outlier_avg_num_processed,
               outlier_avg_num_mapped = group_id %in% has_outlier_avg_num_mapped,
               outlier_avg_num_decoy_fragments = group_id %in% has_outlier_avg_num_decoy_fragments,
               outlier_avg_num_dovetail_fragments = group_id %in% has_outlier_avg_num_dovetail_fragments,
               outlier_avg_num_fragments_filtered_vm = group_id %in% has_outlier_avg_num_fragments_filtered_vm,
               outlier_avg_num_alignments_below_threshold_vm = group_id %in% has_outlier_avg_num_alignments_below_threshold_vm,
               outlier_avg_percent_mapped = group_id %in% has_outlier_avg_percent_mapped)]

  # Now collect the contrasts with the potential outliers
  message("Determining contrasts with outlier groups...")
  con_outlier_avg_frag_length_mean <- cm_dt[
    outlier_avg_frag_length_mean == TRUE,
    unique(id)
  ]
  con_outlier_avg_frag_length_sd <- cm_dt[
    outlier_avg_frag_length_sd == TRUE,
    unique(id)
  ]
  con_outlier_avg_num_eq_classes <- cm_dt[
    outlier_avg_num_eq_classes == TRUE,
    unique(id)
  ]
  con_outlier_avg_num_processed <- cm_dt[
    outlier_avg_num_processed == TRUE,
    unique(id)
  ]
  con_outlier_avg_num_mapped <- cm_dt[
    outlier_avg_num_mapped == TRUE,
    unique(id)
  ]
  con_outlier_avg_num_decoy_fragments <- cm_dt[
    outlier_avg_num_decoy_fragments == TRUE,
    unique(id)
  ]
  con_outlier_avg_num_dovetail_fragments <- cm_dt[
    outlier_avg_num_dovetail_fragments == TRUE,
    unique(id)
  ]
  con_outlier_avg_num_fragments_filtered_vm <- cm_dt[
    outlier_avg_num_fragments_filtered_vm == TRUE,
    unique(id)
  ]
  con_outlier_avg_num_alignments_below_threshold_vm <- cm_dt[
    outlier_avg_num_alignments_below_threshold_vm == TRUE,
    unique(id)
  ]
  con_outlier_avg_percent_mapped <- cm_dt[
    outlier_avg_percent_mapped == TRUE,
    unique(id)
  ]

  # Return a list of these IDs for use in the main annotation script
  result <- list(
    outlier_avg_frag_length_mean = con_outlier_avg_frag_length_mean,
    outlier_avg_frag_length_sd = con_outlier_avg_frag_length_sd,
    outlier_avg_num_eq_classes = con_outlier_avg_num_eq_classes,
    outlier_avg_num_processed = con_outlier_avg_num_processed,
    outlier_avg_num_mapped = con_outlier_avg_num_mapped,
    outlier_avg_num_decoy_fragments = con_outlier_avg_num_decoy_fragments,
    outlier_avg_num_dovetail_fragments = con_outlier_avg_num_dovetail_fragments,
    outlier_avg_num_fragments_filtered_vm = con_outlier_avg_num_fragments_filtered_vm,
    outlier_avg_num_alignments_below_threshold_vm = con_outlier_avg_num_alignments_below_threshold_vm,
    outlier_avg_percent_mapped = con_outlier_avg_percent_mapped
  )

  return(result)
}


getOutlierFlags <- function(l, sep = ",") {
  mat <- do.call(cbind, l)
  result <- apply(mat, 1, function(row) {
    cols <- names(l)[which(row)]
    paste(cols, collapse = sep)
  })

  return(result)
}
