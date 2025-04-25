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
BIOPROJECT <- "PRJNA379584"
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

# Drop non-estimable batches
keep <- rownames(subset(colData(mae), !batch %in% c("12B", "12C")))
mae <- mae[, keep]

# Drop outlier samples
outliers <- c("SAMN06611403", "SAMN06611688", "SAMN06611416", "SAMN06611395",
              "SAMN06611394", "SAMN06611862", "SAMN06611860")
mae <- mae[, !mae$BioSample %in% outliers]

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
metadata <- unique(metadata[, .(BioSample, group, batch)])
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
design <- model.matrix(~0 + group + batch, data = metadata)
colnames(design) <- gsub(pattern = "^group", replacement = "", x = colnames(design))

# DEFINE CONTRAST MATRIX : THIS MUST BE MODIFIED FOR EACH EXPERIMENT
cm <- makeContrasts(
  # Batch 1
  AZD8055_vs_DMSO_6hr.1 = AZD8055.6 - DMSO.6,
  AZD8055_vs_DMSO_24hr.1 = (AZD8055.24 - AZD8055.6) - (DMSO.24 - DMSO.6),
  Belinostat_vs_DMSO_6hr.1 = Belinostat.6 - DMSO.6,
  Belinostat_vs_DMSO_24hr.1 = (Belinostat.24 - Belinostat.6) - (DMSO.24 - DMSO.6),
  Entinostat_vs_DMSO_6hr.1 = Entinostat.6 - DMSO.6,
  Entinostat_vs_DMSO_24hr.1 = (Entinostat.24 - Entinostat.6) - (DMSO.24 - DMSO.6),
  Imatinib_vs_DMSO_6hr.1 = Imatinib.6 - DMSO.6,
  Imatinib_vs_DMSO_24hr.1 = (Imatinib.24 - Imatinib.6) - (DMSO.24 - DMSO.6),
  MK_2206_vs_DMSO_6hr.1 = MK_2206.6 - DMSO.6,
  MK_2206_vs_DMSO_24hr.1 = (MK_2206.24 - MK_2206.6) - (DMSO.24 - DMSO.6),
  Tivantinib_vs_DMSO_6hr.1 = Tivantinib.6 - DMSO.6,
  Tivantinib_vs_DMSO_24hr.1 = (Tivantinib.24 - Tivantinib.6) - (DMSO.24 - DMSO.6),
  Topotecan_vs_DMSO_6hr.1 = Topotecan.6 - DMSO.6,
  Topotecan_vs_DMSO_24hr.1 = (Topotecan.24 - Topotecan.6) - (DMSO.24 - DMSO.6),
  YK_4_279_vs_DMSO_6hr.1 = YK_4_279.6 - DMSO.6,
  YK_4_279_vs_DMSO_24hr.1 = (YK_4_279.24 - YK_4_279.6) - (DMSO.24 - DMSO.6),
  
  # Batch 2
  Alisertib_vs_DMSO_6hr.2 = Alisertib.6 - DMSO.6,
  Methylstat_vs_DMSO_24hr.2 = Methylstat.6 - DMSO.6,
  NVP_231_vs_DMSO_6hr.2 = NVP_231.6 - DMSO.6, 
  Obatoclax_vs_DMSO_6hr.2 = Obatoclax.6 - DMSO.6,
  SNX_2112_vs_DMSO_6hr.2 = SNX_2112.6 - DMSO.6,
  Sunitinib_vs_DMSO_6hr.2 = Sunitinib.6 - DMSO.6,
  YM_155_vs_DMSO_6hr.2 = YM_155.6 - DMSO.6,
  
  # Batch 3
  Alisertib_vs_DMSO_24hr.3 = Alisertib.24 - DMSO.24,
  Methylstat_vs_DMSO_24hr.3 = Methylstat.24 - DMSO.24,
  NVP_231_vs_DMSO_24hr.3 = NVP_231.24 - DMSO.24,
  Obatoclax_vs_DMSO_24hr.3 = Obatoclax.24 - DMSO.24,
  SNX_2112_vs_DMSO_24hr.3 = SNX_2112.24 - DMSO.24,
  Sunitinib_vs_DMSO_24hr.3 = Sunitinib.24 - DMSO.24,
  YM_155_vs_DMSO_24hr.3 = YM_155.24 - DMSO.24,
  
  # Batch 4
  At13387_vs_DMSO_6hr.4 = At13387.6 - DMSO.6,
  At13387_vs_DMSO_24hr.4 = (At13387.24 - At13387.6) - (DMSO.24 - DMSO.6),
  BI_2536_vs_DMSO_6hr.4 = BI_2536.6 - DMSO.6,
  BI_2536_vs_DMSO_24hr.4 = (BI_2536.24 - BI_2536.6) - (DMSO.24 - DMSO.6),
  BMS_754807_vs_DMSO_6hr.4 = BMS_754807.6 - DMSO.6,
  BMS_754807_vs_DMSO_24hr.4  = (BMS_754807.24 - BMS_754807.6) - (DMSO.24 - DMSO.6),
  Bortezomib_vs_DMSO_6hr.4 = Bortezomib.6 - DMSO.6,
  Brefeldin_A_vs_DMSO_6hr.4 = Brefeldin_A.6 - DMSO.6,
  Brefeldin_A_vs_DMSO_24hr.4 = (Brefeldin_A.24 - Brefeldin_A.6) - (DMSO.24 - DMSO.6),
  CD_437_vs_DMSO_6hr.4 = CD_437.6 - DMSO.6,
  CD_437_vs_DMSO_24hr.4 = (CD_437.24 - CD_437.6) - (DMSO.24 - DMSO.6),
  DBeQ_vs_DMSO_6hr.4 = DBeQ.6 - DMSO.6,
  DBeQ_vs_DMSO_24hr.4 = (DBeQ.24 - DBeQ.6) - (DMSO.24 - DMSO.6),
  Fingolimod_vs_DMSO_6hr.4 = Fingolimod.6 - DMSO.6,
  Fingolimod_vs_DMSO_24hr.4 = (Fingolimod.24 - Fingolimod.6) - (DMSO.24 - DMSO.6),
  Flavopiridol_vs_DMSO_6hr.4 = Flavopiridol.6 - DMSO.6,
  Flavopiridol_vs_DMSO_24hr.4 = (Flavopiridol.24 - Flavopiridol.6) - (DMSO.24 - DMSO.6),
  Fluvastatin_vs_DMSO_6hr.4 = Fluvastatin.6 - DMSO.6,
  Fluvastatin_vs_DMSO_24hr.4 = (Fluvastatin.24 - Fluvastatin.6) - (DMSO.24 - DMSO.6),
  Parthenolide_vs_DMSO_6hr.4 = Parthenolide.6 - DMSO.6,
  Parthenolide_vs_DMSO_24hr.4 = (Parthenolide.24 - Parthenolide.6) - (DMSO.24 - DMSO.6),
  
  # Batch 5
  Foretinib_vs_DMSO_6hr.5 = Foretinib.6 - DMSO.6,
  Foretinib_vs_DMSO_24hr.5  = (Foretinib.24 - Foretinib.6) - (DMSO.24 - DMSO.6),
  GDC_0941_vs_DMSO_6hr.5 = GDC_0941.6 - DMSO.6,
  GDC_0941_vs_DMSO_24hr.5  = (GDC_0941.24 - GDC_0941.6) - (DMSO.24 - DMSO.6),
  GSK1210151A_vs_DMSO_6hr.5 = GSK1210151A.6 - DMSO.6,
  GSK1210151A_vs_DMSO_24hr.5  = (GSK1210151A.24 - GSK1210151A.6) - (DMSO.24 - DMSO.6),
  GSK461364_vs_DMSO_6hr.5 = GSK461364.6 - DMSO.6,
  GSK461364_vs_DMSO_24hr.5 = (GSK461364.24 - GSK461364.6) - (DMSO.24 - DMSO.6),
  JNJ_26854165_vs_DMSO_6hr.5 = JNJ_26854165.6 - DMSO.6,
  JNJ_26854165_vs_DMSO_24hr.5 = (JNJ_26854165.24 - JNJ_26854165.6) - (DMSO.24 - DMSO.6),
  KHS101_vs_DMSO_6hr.5 = KHS101.6 - DMSO.6,
  KHS101_vs_DMSO_24hr.5  = (KHS101.24 - KHS101.6) - (DMSO.24 - DMSO.6),
  KW_2449_vs_DMSO_6hr.5 = KW_2449.6 - DMSO.6,
  KW_2449_vs_DMSO_24hr.5  = (KW_2449.24 - KW_2449.6) - (DMSO.24 - DMSO.6),
  MST_312_vs_DMSO_6hr.5 = MST_312.6 - DMSO.6,
  MST_312_vs_DMSO_24hr.5  = (MST_312.24 - MST_312.6) - (DMSO.24 - DMSO.6),
  Panobinostat_vs_DMSO_6hr.5 = Panobinostat.6 - DMSO.6,
  Panobinostat_vs_DMSO_24hr.5 = (Panobinostat.24 - Panobinostat.6) - (DMSO.24 - DMSO.6),
  Piperlongumine_vs_DMSO_6hr.5 = Piperlongumine.6 - DMSO.6,
  Piperlongumine_vs_DMSO_24hr.5  = (Piperlongumine.24 - Piperlongumine.6) - (DMSO.24 - DMSO.6),
  
  # Batch 6
  AZD6482_vs_DMSO_6hr.6 = AZD6482.6 - DMSO.6,
  AZD6482_vs_DMSO_24hr.6  = (AZD6482.24 - AZD6482.6) - (DMSO.24 - DMSO.6),
  AZD7762_vs_DMSO_6hr.6 = AZD7762.6 - DMSO.6,
  AZD7762_vs_DMSO_24hr.6  = (AZD7762.24 - AZD7762.6) - (DMSO.24 - DMSO.6),
  Afatinib_vs_DMSO_6hr.6 = Afatinib.6 - DMSO.6,
  Afatinib_vs_DMSO_24hr.6 = (Afatinib.24 - Afatinib.6) - (DMSO.24 - DMSO.6),
  Bafetinib_vs_DMSO_6hr.6 = Bafetinib.6 - DMSO.6,
  Bafetinib_vs_DMSO_24hr.6 = (Bafetinib.24 - Bafetinib.6) - (DMSO.24 - DMSO.6),
  Bardoxolone_methyl_vs_DMSO_6hr.6 = Bardoxolone_methyl.6 - DMSO.6,
  Bardoxolone_methyl_vs_DMSO_24hr.6 = (Bardoxolone_methyl.24 - Bardoxolone_methyl.6) - (DMSO.24 - DMSO.6),
  Bexarotene_vs_DMSO_6hr.6 = Bexarotene.6 - DMSO.6,
  Bexarotene_vs_DMSO_24hr.6 = (Bexarotene.24 - Bexarotene.6) - (DMSO.24 - DMSO.6),
  Brivanib_vs_DMSO_6hr.6 = Brivanib.6 - DMSO.6,
  Brivanib_vs_DMSO_24hr.6 = (Brivanib.24 - Brivanib.6) - (DMSO.24 - DMSO.6),
  CHM_1_vs_DMSO_6hr.6 = CHM_1.6 - DMSO.6,
  CHM_1_vs_DMSO_24hr.6  = (CHM_1.24 - CHM_1.6) - (DMSO.24 - DMSO.6),
  Cerulenin_vs_DMSO_6hr.6 = Cerulenin.6 - DMSO.6,
  Cerulenin_vs_DMSO_24hr.6  = (Cerulenin.24 - Cerulenin.6) - (DMSO.24 - DMSO.6),
  Ciclopirox_vs_DMSO_6hr.6 = Ciclopirox.6 - DMSO.6,
  Ciclopirox_vs_DMSO_24hr.6 = (Ciclopirox.24 - Ciclopirox.6) - (DMSO.24 - DMSO.6),
  x5Zx_7_Oxozeanol_vs_DMSO_6hr.6 = x5Zx_7_Oxozeanol.6 - DMSO.6,
  x5Zx_7_Oxozeanol_vs_DMSO_24hr.6 = (x5Zx_7_Oxozeanol.24 - x5Zx_7_Oxozeanol.6) - (DMSO.24 - DMSO.6),
  
  # Batch 7
  Crizotinib_vs_DMSO_6hr.7 = Crizotinib.6 - DMSO.6,
  Crizotinib_vs_DMSO_24hr.7  = (Crizotinib.24 - Crizotinib.6) - (DMSO.24 - DMSO.6),
  Dasatinib_vs_DMSO_6hr.7 = Dasatinib.6 - DMSO.6,
  Dasatinib_vs_DMSO_24hr.7  = (Dasatinib.24 - Dasatinib.6) - (DMSO.24 - DMSO.6),
  Dinaciclib_vs_DMSO_6hr.7 = Dinaciclib.6 - DMSO.6,
  Dinaciclib_vs_DMSO_24hr.7  = (Dinaciclib.24 - Dinaciclib.6) - (DMSO.24 - DMSO.6),
  Docetaxel_vs_DMSO_6hr.7 = Docetaxel.6 - DMSO.6,
  Docetaxel_vs_DMSO_24hr.7 = (Docetaxel.24 - Docetaxel.6) - (DMSO.24 - DMSO.6),
  Doxorubicin_vs_DMSO_6hr.7 = Doxorubicin.6 - DMSO.6,
  Doxorubicin_vs_DMSO_24hr.7 = (Doxorubicin.24 - Doxorubicin.6) - (DMSO.24 - DMSO.6),
  Erastin_vs_DMSO_6hr.7 = Erastin.6 - DMSO.6, 
  Erastin_vs_DMSO_24hr.7 = (Erastin.24 - Erastin.6) - (DMSO.24 - DMSO.6),
  GW_405833_vs_DMSO_6hr.7 = GW_405833.6 - DMSO.6, 
  GW_405833_vs_DMSO_24hr.7  = (GW_405833.24 - GW_405833.6) - (DMSO.24 - DMSO.6),
  Gemcitabine_vs_DMSO_6hr.7 = Gemcitabine.6 - DMSO.6,
  Gemcitabine_vs_DMSO_24hr.7 = (Gemcitabine.24 - Gemcitabine.6) - (DMSO.24 - DMSO.6),
  Gossypol_vs_DMSO_6hr.7 = Gossypol.6 - DMSO.6,
  Gossypol_vs_DMSO_24hr.7 = (Gossypol.24 - Gossypol.6) - (DMSO.24 - DMSO.6),
  HMN_214_vs_DMSO_6hr.7 = HMN_214.6 - DMSO.6,
  HMN_214_vs_DMSO_24hr.7  = (HMN_214.24 - HMN_214.6) - (DMSO.24 - DMSO.6),
  Ibrutinib_vs_DMSO_6hr.7 = Ibrutinib.6 - DMSO.6,
  Ibrutinib_vs_DMSO_24hr.7  = (Ibrutinib.24 - Ibrutinib.6) - (DMSO.24 - DMSO.6),
  
  # Batch 8
  KX2_391_vs_DMSO_6hr.8 = KX2_391.6 - DMSO.6,
  KX2_391_vs_DMSO_24hr.8  = (KX2_391.24 - KX2_391.6) - (DMSO.24 - DMSO.6),
  Linifanib_vs_DMSO_6hr.8 = Linifanib.6 - DMSO.6,
  Linifanib_vs_DMSO_24hr.8  = (Linifanib.24 - Linifanib.6) - (DMSO.24 - DMSO.6),
  Manumycin_A_vs_DMSO_6hr.8 = Manumycin_A.6 - DMSO.6,
  Manumycin_A_vs_DMSO_24hr.8  = (Manumycin_A.24 - Manumycin_A.6) - (DMSO.24 - DMSO.6),
  MG_132_vs_DMSO_6hr.8 = MG_132.6 - DMSO.6,
  MG_132_vs_DMSO_24hr.8  = (MG_132.24 - MG_132.6) - (DMSO.24 - DMSO.6),
  Mithramycin_A_vs_DMSO_6hr.8 = Mithramycin_A.6 - DMSO.6,
  Mithramycin_A_vs_DMSO_24hr.8  = (Mithramycin_A.24 - Mithramycin_A.6) - (DMSO.24 - DMSO.6),
  Mitomycin_C_vs_DMSO_6hr.8 = Mitomycin_C.6 - DMSO.6,
  Mitomycin_C_vs_DMSO_24hr.8  = (Mitomycin_C.24 - Mitomycin_C.6) - (DMSO.24 - DMSO.6),
  MK_1775_vs_DMSO_6hr.8 = MK_1775.6 - DMSO.6,
  MK_1775_vs_DMSO_24hr.8  = (MK_1775.24 - MK_1775.6) - (DMSO.24 - DMSO.6),
  ML210_vs_DMSO_6hr.8 = ML210.6 - DMSO.6,
  ML210_vs_DMSO_24hr.8  = (ML210.24 - ML210.6) - (DMSO.24 - DMSO.6),
  Necrosulfonamide_vs_DMSO_6hr.8 = Necrosulfonamide.6 - DMSO.6, 
  Necrosulfonamide_vs_DMSO_24hr.8  = (Necrosulfonamide.24 - Necrosulfonamide.6) - (DMSO.24 - DMSO.6),
  Neratinib_vs_DMSO_6hr.8 = Neratinib.6 - DMSO.6, 
  Neratinib_vs_DMSO_24hr.8  = (Neratinib.24 - Neratinib.6) - (DMSO.24 - DMSO.6),
  Nutlin_3_vs_DMSO_6hr.8 = Nutlin_3.6 - DMSO.6,
  Nutlin_3_vs_DMSO_24hr.8  = (Nutlin_3.24 - Nutlin_3.6) - (DMSO.24 - DMSO.6),
  
  # Batch 9
  FK866_vs_DMSO_6hr.9 = FK866.6 - DMSO.6,
  FK866_vs_DMSO_24hr.9  = (FK866.24 - FK866.6) - (DMSO.24 - DMSO.6),
  Istradefylline_vs_DMSO_6hr.9 = Istradefylline.6 - DMSO.6,
  Istradefylline_vs_DMSO_24hr.9  = (Istradefylline.24 - Istradefylline.6) - (DMSO.24 - DMSO.6),
  PAC_1_vs_DMSO_6hr.9 = PAC_1.6 - DMSO.6,
  PAC_1_vs_DMSO_24hr.9  = (PAC_1.24 - PAC_1.6) - (DMSO.24 - DMSO.6),
  Paclitaxel_vs_DMSO_6hr.9 = Paclitaxel.6 - DMSO.6,
  Paclitaxel_vs_DMSO_24hr.9 = (Paclitaxel.24 - Paclitaxel.6) - (DMSO.24 - DMSO.6),
  PF_3758309_vs_DMSO_6hr.9 = PF_3758309.6 - DMSO.6,
  PF_3758309_vs_DMSO_24hr.9  = (PF_3758309.24 - PF_3758309.6) - (DMSO.24 - DMSO.6),
  PHA_665752_vs_DMSO_6hr.9 = PHA_665752.6 - DMSO.6, 
  PHA_665752_vs_DMSO_24hr.9  = (PHA_665752.24 - PHA_665752.6) - (DMSO.24 -  DMSO.6),
  Rigosertib_vs_DMSO_6hr.9 =  Rigosertib.6 - DMSO.6,
  Rigosertib_vs_DMSO_24hr.9  = (Rigosertib.24 - Rigosertib.6) - (DMSO.24 - DMSO.6),
  SB_225002_vs_DMSO_6hr.9 = SB_225002.6 - DMSO.6, 
  SB_225002_vs_DMSO_24hr.9  = (SB_225002.24 - SB_225002.6) - (DMSO.24 - DMSO.6),
  SB_743921_vs_DMSO_6hr.9 = SB_743921.6 - DMSO.6,
  SB_743921_vs_DMSO_24hr.9  = (SB_743921.24 - SB_743921.6) - (DMSO.24 - DMSO.6),
  SCH_529074_vs_DMSO_6hr.9 = SCH_529074.6 - DMSO.6,
  SCH_529074_vs_DMSO_24hr.9  = (SCH_529074.24 - SCH_529074.6) - (DMSO.24 - DMSO.6),
  SCH_79797_vs_DMSO_6hr.9 = SCH_79797.6 - DMSO.6,
  SCH_79797_vs_DMSO_24hr.9  = (SCH_79797.24 - SCH_79797.6) - (DMSO.24 - DMSO.6),
  
  # Batch 10a
  Bortezomib_vs_DMSO_6hr.10A = Bortezomib.6 - DMSO.6,
  Bortezomib_vs_DMSO_24hr.10A  = (Bortezomib.24 - Bortezomib.6) - (DMSO.24 - DMSO.6),
  Erastin_vs_DMSO_6hr.10A = Erastin.6 - DMSO.6, 
  Erastin_vs_DMSO_24hr.10A  = (Erastin.24 - Erastin.6) - (DMSO.24 - DMSO.6),
  Gemcitabine_vs_DMSO_6hr.10A = Gemcitabine.6 -  DMSO.6,
  Gemcitabine_vs_DMSO_24hr.10A  = (Gemcitabine.24 - Gemcitabine.6) - (DMSO.24 - DMSO.6),
  GSK1210151A_vs_DMSO_6hr.10A = GSK1210151A.6 - DMSO.6, 
  GSK1210151A_vs_DMSO_24hr.10A  = (GSK1210151A.24 - GSK1210151A.6) - (DMSO.24 - DMSO.6),
  
  # Batch 10b
  Bortezomib_vs_DMSO_6hr.10B = Bortezomib.6 - DMSO.6,
  Bortezomib_vs_DMSO_24hr.10B  = (Bortezomib.24 - Bortezomib.6) - (DMSO.24 - DMSO.6),
  GSK461364_vs_DMSO_6hr.10B = GSK461364.6 - DMSO.6, 
  GSK461364_vs_DMSO_24hr.10B  = (GSK461364.24 - GSK461364.6) - (DMSO.24 - DMSO.6),
  Pluripotin_vs_DMSO_6hr.10B = Pluripotin.6 - DMSO.6, 
  Pluripotin_vs_DMSO_24hr.10B = (Pluripotin.24 - Pluripotin.6) - (DMSO.24 - DMSO.6),
  Sorafenib_vs_DMSO_6hr.10B = Sorafenib.6 - DMSO.6,
  Sorafenib_vs_DMSO_24hr.10B  = (Sorafenib.24 - Sorafenib.6) - (DMSO.24 - DMSO.6),
  TW_37_vs_DMSO_6hr.10B = TW_37.6 - DMSO.6,
  TW_37_vs_DMSO_24hr.10B  = (TW_37.24 - TW_37.6) - (DMSO.24 - DMSO.6),
  Vorinostat_vs_DMSO_6hr.10B = Vorinostat.6 - DMSO.6, 
  Vorinostat_vs_DMSO_24hr.10B = (Vorinostat.24 - Vorinostat.6) - (DMSO.24 - DMSO.6),
  
  # Batch 11a
  Arachidonyl_vs_DMSO_6hr.11A = Arachidonyl.6 - DMSO.6, 
  Arachidonyl_vs_DMSO_24hr.11A = (Arachidonyl.24 - Arachidonyl.6) - (DMSO.24 - DMSO.6), 
  Epigallocatechin_gallate_vs_DMSO_6hr.11A = Epigallocatechin_gallate.6 - DMSO.6, 
  Epigallocatechin_gallate_vs_DMSO_24hr.11A = (Epigallocatechin_gallate.24 - Epigallocatechin_gallate.6) - (DMSO.24 - DMSO.6),
  Parbendazole_vs_DMSO_6hr.11A = Parbendazole.6 - DMSO.6, 
  Parbendazole_vs_DMSO_24hr.11A  = (Parbendazole.24 - Parbendazole.6) - (DMSO.24 - DMSO.6),
  Rosiglitazone_vs_DMSO_6hr.11A = Rosiglitazone.6 - DMSO.6, 
  Rosiglitazone_vs_DMSO_24hr.11A  = (Rosiglitazone.24 - Rosiglitazone.6) - (DMSO.24 - DMSO.6),
  Vincristine_vs_DMSO_6hr.11A = Vincristine.6 -  DMSO.6, 
  Vincristine_vs_DMSO_24hr.11A  = (Vincristine.24 - Vincristine.6) - (DMSO.24 - DMSO.6), 
  Y_27632_vs_DMSO_6hr.11A = Y_27632.6 - DMSO.6, 
  Y_27632_vs_DMSO_24hr.11A = (Y_27632.24 - Y_27632.6) - (DMSO.24 - DMSO.6),
  
  # Batch 11b
  BIO_vs_DMSO_6hr.11B = BIO.6 - DMSO.6, 
  KU_0063794_vs_DMSO_6hr.11B = KU_0063794.6 - DMSO.6, 
  LY_2183240_vs_DMSO_6hr.11B = LY_2183240.6 - DMSO.6, 
  PF_573228_vs_DMSO_6hr.11B = PF_573228.6 - DMSO.6,
  
  # Batch 12a
  AUY922_vs_DMSO_6hr.12A = AUY922.6 - DMSO.6, 
  AUY922_vs_DMSO_24hr.12A  = (AUY922.24 - AUY922.6) - (DMSO.24 - DMSO.6), 
  BGJ398_vs_DMSO_6hr.12A = BGJ398.6 - DMSO.6, 
  BGJ398_vs_DMSO_24hr.12A  = (BGJ398.24 - BGJ398.6) - (DMSO.24 - DMSO.6), 
  BIO_vs_DMSO_6hr.12A = BIO.6 - DMSO.6, 
  BIO_vs_DMSO_24hr.12A  = (BIO.24 - BIO.6) - (DMSO.24 - DMSO.6), 
  BKM120_vs_DMSO_6hr.12A = BKM120.6 - DMSO.6, 
  BKM120_vs_DMSO_24hr.12A  = (BKM120.24 - BKM120.6) - (DMSO.24 - DMSO.6), 
  BYL719_vs_DMSO_6hr.12A = BYL719.6 - DMSO.6, 
  BYL719_vs_DMSO_24hr.12A  = (BYL719.24 - BYL719.6) - (DMSO.24 - DMSO.6),
  Ciclopirox_vs_DMSO_24hr.12A = Ciclopirox.24 - DMSO.24, 
  INC280_vs_DMSO_6hr.12A = INC280.6 - DMSO.6,
  INC280_vs_DMSO_24hr.12A  = (INC280.24 - INC280.6) - (DMSO.24 - DMSO.6),
  KU_0063794_vs_DMSO_24hr.12A = KU_0063794.24 - DMSO.24,
  LDK378_vs_DMSO_6hr.12A = LDK378.6 - DMSO.6, 
  LDK378_vs_DMSO_24hr.12A  = (LDK378.24 - LDK378.6) - (DMSO.24 - DMSO.6),
  LEE011_vs_DMSO_6hr.12A = LEE011.6 - DMSO.6, 
  LEE011_vs_DMSO_24hr.12A  = (LEE011.24 - LEE011.6) - (DMSO.24 - DMSO.6),
  LY_2183240_vs_DMSO_24hr.12A = LY_2183240.24 - DMSO.24, 
  Manumycin_A_vs_DMSO_6hr.12A = Manumycin_A.6 - DMSO.6, 
  NIlotinib_vs_DMSO_6hr.12A = NIlotinib.6 - DMSO.6, 
  NIlotinib_vs_DMSO_24hr.12A  = (NIlotinib.24 - NIlotinib.6) - (DMSO.24 - DMSO.6),
  PF_573228_vs_DMSO_24hr.12A  = PF_573228.24 - DMSO.24, 
  Piperlongumine_vs_DMSO_6hr.12A = Piperlongumine.6 - DMSO.6, 
  PKC412_vs_DMSO_6hr.12A = PKC412.6 - DMSO.6, 
  PKC412_vs_DMSO_24hr.12A  = (PKC412.24 - PKC412.6) - (DMSO.24 - DMSO.6),
  Pluripotin_vs_DMSO_24hr.12A = Pluripotin.24 - DMSO.24,
  RAD001_vs_DMSO_6hr.12A = RAD001.6 - DMSO.6, 
  RAD001_vs_DMSO_24hr.12A  = (RAD001.24 - RAD001.6) - (DMSO.24 - DMSO.6), 
  SCH_79797_vs_DMSO_24hr.12A = SCH_79797.24 - DMSO.24,
  SU11274_vs_DMSO_6hr.12A = SU11274.6 - DMSO.6, 
  SU11274_vs_DMSO_24hr.12A  = (SU11274.24 - SU11274.6) - (DMSO.24 - DMSO.6),
  Y_27632_vs_DMSO_6hr.12A = Y_27632.6 - DMSO.6,
  
  # Batch 12b -- not estimable
  #CAY10618_in_Methanol_vs_Methanol_6hr.12B = CAY10618_in_Methanol.6 - Methanol.6, 
  #CAY10618_in_Methanol_vs_Methanol_24hr.12B = (CAY10618_in_Methanol.24 - CAY10618_in_Methanol.6) - (Methanol.24 - Methanol.6),
  
  # Batch 12c -- not estimable
  #Leptomycin_B_in_Ethanol_vs_Ethanol_6hr.12C = Leptomycin_B_in_Ethanol.6 - Ethanol.6, 
  #Leptomycin_B_in_Ethanol_vs_Ethanol_24hr.12C = (Leptomycin_B_in_Ethanol.24 - Leptomycin_B_in_Ethanol.6) - (Ethanol.24 - Ethanol.6),
  
  levels = design
)

y <- DGEList(counts = counts, samples = metadata)
L <- min(y$samples$lib.size) / 1e6
before <- aveLogCPM(y)
keep <- rowSums(cpm(y) > L) >= (ncol(y) * 0.7)
y <- y[keep,, keep.lib.sizes = FALSE]
after <- aveLogCPM(y)

.pardefault <- par()
png(here(WD, "results", "figures", "abundance-hist.png"), width = 12, height = 6, units = "in", res = 150)
par(mfrow=c(1, 2))
hist(before, breaks=100, main="Before Filtering", xlab="AveLogCPM")
hist(after, breaks=100, main="After Filtering", xlab="AveLogCPM")
dev.off()
par(.pardefault)

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
print(plot_boxplot(lcpm, y$samples, fillBy = "group", rle = TRUE) + 
        labs(y = "RLE") + 
        theme_coriell() +
        theme(legend.position="none")
      )
dev.off()

png(here(WD, "results", "figures", "logCPM-density.png"), width = 11, height = 7, units = "in", res = 150)
print(plot_density(lcpm, y$samples, colBy = "group") + 
        labs(x = "logCPM") + 
        theme_coriell() +
        theme(legend.position="none")
      )
dev.off()

pca.res <- pca(lcpm, y$samples)
png(here(WD, "results", "figures", "pca-biplot.png"), width = 8, height = 6, units = "in", res = 150)
print(biplot(pca.res, colby = "group", lab = NULL))
dev.off()

png(here(WD, "results", "figures", "quantro-plot.png"), width = 8, height = 6, units = "in", res = 150)
quantroPlot(qtest)
dev.off()

# quickmap(
#   lcpm, 
#   removeVar = 0.9, 
#   annotation_col = y$samples[, "group", drop = FALSE],
#   main = "10% Most Variable Features", 
#   filename = here(WD, "results", "figures", "heatmap.png")
# )

message("Done.")
