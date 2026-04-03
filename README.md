# TEDEdb

Transposable element (TE) expression in cancer is a well defined characteristic of the cancer 
transcriptome ([source](https://www.nature.com/articles/nrc.2017.35)) however, TE expression is 
traditionally overlooked in RNA-seq analyses due to difficulties involving mappability, sequence 
similarity, and quality of annotations ([source](https://www.nature.com/articles/s41576-020-0251-y)). 

This repository contains the scripts needed for uniform processing of publicly available RNA-seq 
data from GEO with the explicit purpose of accurately quantifying gene and transposable element 
expression in cancer cell lines treated with epigenetic and cytotoxic drugs. This project has three 
main arms:

1. A pipeline for downloading and processing (quantifying) gene and TE expression from publicly 
available data hosted on NCBI GEO. This custom pipeline implements a modification of the
[REDiscoverTE](https://www.nature.com/articles/s41467-019-13035-2) method which has been shown to
accurately quantify transcript expression as well a TE expression at the loci level. One of the 
important modifications to this existing pipeline is the inclusion of decoy sequences to avoid 
spurious annotations as well a implementing proper offsets to transcript level counts when 
summarizing to the gene-level ([source](https://f1000research.com/articles/4-1521/v1)).
2. The second arm of this repository is a uniform differential expression analysis which assesses 
and attempts to control for violations of the assumptions made by global scaling normalization 
methods. An important step in comparing expression values between conditions during differential 
expression analysis is normalization. Violations of these assumptions (i.e. global expression 
changes) can lead to spurious DE results. Here, I have implemented an automated detection for 
violations of these global normalization assumptions using the [quantro](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4495646/) 
R package. If violations of these assumptions are detected then normalization is performed using
smooth quantile normalization ([qsmooth](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5862355/))
in order to preserve these global expression differences. Appropriately normalized data are then 
tested for differential expression against a fold-change threshold using the [limma::treat](https://pubmed.ncbi.nlm.nih.gov/25605792/)
method. These differential expression results, along with the fitted models, raw and normalized 
counts, are then exported to `SummarizedExperiment` objects for each individual bioproject.
3. The last arm of this project involves collecting each uniformly processed bioproject's 
differential expression results into a single database back-end and exposing the data with a [Shiny](https://www.rstudio.com/products/shiny/) application. The Shiny application implements 
several novel meta-analysis and visualization functions for the exploration of large-scale patterns 
of differential expression across datasets. It is also worth mentioning that considerable effort was
taken to standardize the cell type and drug annotations for all samples in this database. Cell types 
are given uniform labels from [Cellosaurus](https://www.cellosaurus.org/) and the [Cancer Cell Line Encyclopedia](https://sites.broadinstitute.org/ccle) and drugs are assigned functional labels 
from [PubChem](https://pubchem.ncbi.nlm.nih.gov/) and The [Broad Drug Repurposing Hub](https://repo-hub.broadinstitute.org/repurposing#home) when available. The code for the actual Shiny 
application implementation is located in a separate github [repository](https://github.com/coriell-research/TEDEdb-app).
The meta-analysis application has been useful for hypothesis generation. For instance, observing 
drug treatments that result in expression patterns which deviate from other samples in their drug 
class has led to further studies regarding CDK12 inhibition. 

## Directory structure

- appdata/: contains all data used to create and used by the meta-app. scripts/appdata.R details the
creation of the SummarizedExperiment object used as the backend for the database application. Since
this data is fairly large it is not present in this repository.
- bioprojects/: contains subdirectories named by their NCBI bioproject ID from GEO. Each subdirectory 
contains all raw and processed data, analysis scripts, and analysis output for that particular 
bioproject. Currently, only the analysis scripts and annotation information are present in this 
repository. Like the appdata above, raw data is housed on the bioinformatics servers at Coriell.
- scripts/: contains all processing scripts used to download data from SRA and
process for use in the meta-app. The scripts generally describe the processing steps for a single 
bioproject.
- analysis/: contains all analysis scripts for generating figures and data used in the publication
