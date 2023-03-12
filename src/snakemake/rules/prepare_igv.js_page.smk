# Define variable containing all files in the directory
SAMPLES, = glob_wildcards("out/nextflow/nfcore_methylseq/methyldackel/{sample}.markdup.sorted_CpG.bedGraph")

rule prepare_igvjs_page:
    input:
        expand(
            "out/ucsc/bedGraphToBigWig_chrominfo-GRCh37/cut/_-f1-4/sort/bed/nextflow/nfcore_methylseq/methyldackel/{sample}.markdup.sorted_CpG.bw",
            sample = SAMPLES
        )

