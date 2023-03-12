
if (!require("BiocManager", quietly = TRUE)) {
    install.packages(
        "BiocManager",
        quiet = TRUE
    )
}
packages <- c(
  # "data.table",
  # "dplyr",
  "openxlsx",
  "factoextra",
  # "methylSig",
  "DSS",
  "bsseq",
  "ChIPpeakAnno",
  "irlba",
  # "ChIPseeker",
  # "methyAnalysis",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  # "EnhancedVolcano",
  # "data.table",
  "ggpubr",
  "ggrepel",
  # "ComplexHeatmap",
  "tidyr"
)
BiocManager::install(
    packages,
    update = FALSE,
    quiet = TRUE,
    Ncpus = parallel::detectCores()
)
invisible(
    lapply(
      packages,
      library,
      character.only = TRUE
    )
)

sample_metadata <- read.csv(
  "src/sample_sheets/sample_sheet.csv"
)
# Here I temporarily subset the sample sheets to already available samples
suffix <- ".markdup.sorted_CpG.bedGraph"
input_dir <- "out/nextflow/nfcore_methylseq/methyldackel/"
files <- list.files(
    path = input_dir,
    pattern = suffix
)
samples <- sub(
    suffix,
    "",
    files
)
sample_metadata <- sample_metadata[sample_metadata$sample %in% samples,]
sample_metadata

# This will be used to load files with read.bismarck()
sample_metadata$bedGraph <- paste0(
    'out/nextflow/nfcore_methylseq/methyldackel/',
    sample_metadata$sample,
    '.markdup.sorted_CpG.bedGraph'
)

outdir <- "out/r/compute_pca/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bsseq_stranded = bsseq::read.bismark(
    files = sample_metadata$bedGraph,
    colData = sample_metadata,
    rmZeroCov = FALSE,
    strandCollapse = FALSE,
    BACKEND = "HDF5Array",
    dir = "out/r/compute_pca/hdf5a",
    replace = TRUE
)
bsseq_stranded

mean_cov <- rowMeans(getCoverage(bsseq_stranded))
quantile(mean_cov)

cg_sd <- rowSds(
    getMeth(
        bsseq_stranded,
        type="raw"
    ),
)
# Replace NA by 0 to remove them later with low sd sites
cg_sd[is.na(cg_sd)] <- 0
quantile(cg_sd)

length(mean_cov)
length(cg_sd)

# Subset for the sole purpose to have 
bsseq_subset <- bsseq_stranded[
    mean_cov > 5 & cg_sd > 0.2,
]
bsseq_subset

pca_result <- prcomp_irlba(
    x = getMeth(
        bsseq_subset,
        type="raw"
    ),
    n = 3,
    scale = FALSE
)

# Extract the first two principal components
pc1 <- pca_result$v[, 1]
pc2 <- pca_result$v[, 2]

df <- data.frame(
    data.frame(pca_result$rotation),
    colData(bsseq_subset)
)

# Plot the first two principal components
p <- ggplot(
    df,
    aes(
        x = PC1,
        y = PC2,
        shape = tissue,
        label = sample,
        color = patient_number
    )
)
p <- p + geom_point()
p <- p + geom_text_repel()
# Strangely does not do what is needed with irlba input.
# p <- fviz_pca_ind(pca_result)

ggsave(
    p,
    file = "out/r/compute_pca/pca_ind.pdf",
    width = 10,
    height = 10
)
