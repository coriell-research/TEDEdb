## Differential Expression Analysis
##
## NOTE: in order to successfully run the analysis script a contrast matrix must
## be defined for each experiment individually. It is also necessary to change
## the BIOPROJECT variable below to reflect the desired bioproject.
## 
suppressPackageStartupMessages(library(coriell))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(PCAtools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(qsmooth))
suppressPackageStartupMessages(library(quantro))
suppressPackageStartupMessages(library(MultiAssayExperiment))

# Set up global defaults
BIOPROJECT <- "PRJNA578935"
WD <- here("bioprojects", BIOPROJECT)
CORES <- 12
QTEST_ALPHA <- 0.01
FC <- 1.2

# Create the results directory structure if it does not exist
if (!dir.exists(here(WD, "results"))) {
  message("Creating results directories...")
  dir.create(here(WD, "results"))
  dir.create(here(WD, "results", "rds-files"))
  dir.create(here(WD, "results", "figures"))
  dir.create(here(WD, "results", "data-files"))
}


# Filter counts -----------------------------------------------------------


# Read in the MultiAssayExperiment file which contains gene and TE counts
message("Reading in MultiAssayExperiment data...")
mae <- readRDS(here(WD, "MultiAssayExperiment.rds"))
stopifnot("'group' column not found in MAE object metadata" = "group" %in% colnames(colData(mae)))
stopifnot("'BioSample' column not found in MAE object metadata" = "BioSample" %in% colnames(colData(mae)))

# Import the GTF .rds file used in REdiscoverTE gene annotations to get
message("Reading in GTF annotation...")
gtf <- "/mnt/data/gdata/human/REdiscoverTE_hg38/gencode.v26.annotation.gtf.rds"
gtf_gr <- readRDS(gtf) 
coding <- gtf_gr[gtf_gr$type == "gene" & gtf_gr$gene_type == "protein_coding", ]
coding <- unique(coding$gene_name)

# Extract the count matrices from the MAE object
genes <- assay(mae, "gene counts")
repeats <- assay(mae, "RE counts")

# Subset for only the features of interest -- protein coding genes and
#  TEs from the 'main' families and combine to single counts matrix
message("Filtering counts...")
repeats <- repeats[rownames(repeats) %like% "^LINE\\.|^SINE\\.|^LTR\\.|.*\\.SVA\\..*|^DNA\\.", ]
genes <- genes[rownames(genes) %chin% coding, ]
counts <- rbind(genes, repeats)

# Sum technical replicates if present -- has the effect of renaming all samples
#  by their BioSample ID
message("Summarizing counts to the BioSample level...")
counts <- sumTechReps(counts, ID = mae$BioSample)
mae$bases <- as.integer(mae$bases)
metadata <- as.data.table(colData(mae))
metadata <- unique(metadata[, .(BioSample, group)])
setDF(metadata, rownames = metadata$BioSample)
metadata <- metadata[colnames(counts), ]
stopifnot("All rownames of metadata do not match colnames of counts" = all(rownames(metadata) == colnames(counts)))


# Test for global expression differences ----------------------------------


# Test for global expression differences between groups
message("Testing for violations of global scaling normalization...")
doParallel::registerDoParallel(cores = CORES)
logcounts <- cpm(counts, log = TRUE)
qtest <- quantro(logcounts, groupFactor = metadata$group, B = 1e3)
qaov <- anova(qtest)
perm_pval <- quantroPvalPerm(qtest)
qaov_pval <- qaov[["Pr(>F)"]][1]


# Differential expression testing -----------------------------------------


# Create experimental design
message("Differential expression testing pipeline...")
design <- model.matrix(~0 + group, data = metadata)
colnames(design) <- gsub(pattern = "^group", replacement = "", x = colnames(design))

# DEFINE CONTRAST MATRIX : THIS MUST BE MODIFIED FOR EACH EXPERIMENT
cm <- makeContrasts(
  H1373.MRTX849_QDX1_24hr_vs_H1373.Vehicle_QDX1_24hr = H1373.MRTX849_QDX1_24hr - H1373.Vehicle_QDX1_24hr,
  H1373.MRTX849_QDX5_24hr_vs_H1373.Vehicle_QDX5_24hr = H1373.MRTX849_QDX5_24hr - H1373.Vehicle_QDX5_24hr,
  H2030.MRTX849_QDX1_24hr_vs_H2030.Vehicle_QDX1_24hr = H2030.MRTX849_QDX1_24hr - H2030.Vehicle_QDX1_24hr,
  H2030.MRTX849_QDX7_24hr_vs_H2030.Vehicle_QDX7_24hr = H2030.MRTX849_QDX7_24hr - H2030.Vehicle_QDX7_24hr,
  H2122.MRTX849_Palbociclib_QDX1_24hr_vs_H2122.Vehicle_QDX1_24hr = H2122.MRTX849_Palbociclib_QDX1_24hr - H2122.Vehicle_QDX1_24hr,
  H2122.MRTX849_Palbociclib_QDX7_24hr_vs_H2122.Vehicle_QDX7_24hr = H2122.MRTX849_Palbociclib_QDX7_24hr - H2122.Vehicle_QDX7_24hr,
  H2122.MRTX849_QDx1_24h_vs_H2122.Vehicle_QDX1_24hr = H2122.MRTX849_QDx1_24h - H2122.Vehicle_QDX1_24hr,
  H2122.MRTX849_QDx7_24h_vs_H2122.Vehicle_QDX7_24hr = H2122.MRTX849_QDx7_24h - H2122.Vehicle_QDX7_24hr,
  H2122.MRTX849_QDX1_24hr_vs_H2122.Vehicle_QDX1_24hr = H2122.MRTX849_QDX1_24hr - H2122.Vehicle_QDX1_24hr,
  H2122.MRTX849_QDX7_24hr_vs_H2122.Vehicle_QDX7_24hr = H2122.MRTX849_QDX7_24hr - H2122.Vehicle_QDX7_24hr,
  H2122.Palbociclib_QDX1_24hr_vs_H2122.Vehicle_QDX1_24hr = H2122.Palbociclib_QDX1_24hr - H2122.Vehicle_QDX1_24hr,
  H2122.Palbociclib_QDX7_24hr_vs_H2122.Vehicle_QDX7_24hr = H2122.Palbociclib_QDX7_24hr - H2122.Vehicle_QDX7_24hr,
  MiaPaCa2.MRTX849_QDx1_24hr_vs_MiaPaCa2.Vehicle_QDx1_24hr = MiaPaCa2.MRTX849_QDx1_24hr - MiaPaCa2.Vehicle_QDx1_24hr,
  MiaPaCa2.MRTX849_QDx5_24hr_vs_MiaPaCa2.Vehicle_QDx5_24hr = MiaPaCa2.MRTX849_QDx5_24hr - MiaPaCa2.Vehicle_QDx5_24hr,
  levels = design
)

y <- DGEList(counts = counts, samples = metadata)
keep <- rowSums(y$counts > 10) >= ncol(y) * 0.7
y <- y[keep,, keep.lib.sizes = FALSE]

# If global normalization assumptions are violated then perform qsmooth and set
#  an offset value for each gene. If not, use TMM scaling factors
if (perm_pval < QTEST_ALPHA || qaov_pval < QTEST_ALPHA) {
  message("Global scaling normalizations violated. Performing QSmooth normalization...")
  qs <- qsmooth(y$counts, group_factor = y$samples$group)
  qsd <- qsmoothData(qs)
  offsets <- log(y$counts + 0.1) - log(qsd + 0.1)
  y <- scaleOffset(y, offsets)
} else {
  message("Performing TMM normalization...")
  y <- calcNormFactors(y, method = "TMM")
}
# Fit linear models with limma
message("Fitting models...")
lcpm <- cpm(y, offset = y$offset, log = TRUE)
fit <- lmFit(lcpm, design, robust = TRUE)
fit2 <- contrasts.fit(fit, cm)
fit2 <- treat(fit2, fc = FC, trend = TRUE, robust = TRUE)

# Collect all results into a single data.table
message("Performing differential expression against a logFC cutoff with `limma::treat`")
extractResults <- function(x) {
  res <- topTreat(fit2, coef = x, number = Inf, confint = TRUE)
  res$z <- zscoreT(res$t, df = fit2$df.total)
  d <- as.data.table(res, keep.rownames = "feature_id")
}
results <- lapply(colnames(cm), extractResults)
names(results) <- colnames(cm)
dt <- rbindlist(results, idcol = "contrast")

# Add standard error column
dt[, SE := (CI.R - CI.L) / 3.92]


# Save results ------------------------------------------------------------

message("Saving results as a SummarizedExperiment object...")
se <- SummarizedExperiment(
  assays = list("counts" = y$counts, "logcounts" = lcpm),
  colData = metadata,
  metadata = list("fit" = fit2, "de" = dt)
)
saveRDS(se, here(WD, "results", "rds-files", "se.rds"))


# Diagnostic plots --------------------------------------------------------

message("Creating diagnostic plots...")
png(here(WD, "results", "figures", "sa-plot.png"), width = 8, height = 6, units = "in", res = 150)
plotSA(fit2)
dev.off()

png(here(WD, "results", "figures", "logCPM-boxplots.png"), width = 11, height = 7, units = "in", res = 150)
print(plot_boxplot(lcpm, y$samples, fillBy = "group", rle = TRUE) + labs(y = "RLE") + theme_coriell())
dev.off()

png(here(WD, "results", "figures", "logCPM-density.png"), width = 11, height = 7, units = "in", res = 150)
print(plot_density(lcpm, y$samples, colBy = "group") + labs(x = "logCPM") + theme_coriell())
dev.off()

pca.res <- pca(lcpm, y$sample)
png(here(WD, "results", "figures", "pca-biplot.png"), width = 8, height = 6, units = "in", res = 150)
print(biplot(pca.res, colby = "group", legendPosition = "bottom", lab = NULL))
dev.off()

png(here(WD, "results", "figures", "quantro-plot.png"), width = 8, height = 6, units = "in", res = 150)
quantroPlot(qtest)
dev.off()

quickmap(
  lcpm, 
  removeVar = 0.9, 
  annotation_col = y$samples[, "group", drop = FALSE],
  main = "10% Most Variable Features", 
  filename = here(WD, "results", "figures", "heatmap.png")
)

message("Done.")
