library("bsseq")

message("Reading metadata...")
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

# row.names should be declared to be used as row.names of the BSseq object
# because DMLtest() uses these row.namesto select samples
row.names(sample_metadata) <- sample_metadata$sample
sample_metadata

message("Reading bedGraph files...")
bsseq_stranded = bsseq::read.bismark(
    files = sample_metadata$bedGraph,
    colData = sample_metadata,
    rmZeroCov = FALSE,
    strandCollapse = FALSE,
    BACKEND = "HDF5Array",
    dir = "out/r/prepare_smoothed_bsseq/hdf5a",
    replace = TRUE
)
bsseq_stranded
message("Smoothing...")
bsseq_smoothed <- BSmooth(
    BSseq = bsseq_stranded, 
    verbose = TRUE
)
bsseq_smoothed
saveRDS(
    bsseq_smoothed,
    file = "out/r/prepare_smoothed_bsseq/bsseq_smoothed.rds"
)