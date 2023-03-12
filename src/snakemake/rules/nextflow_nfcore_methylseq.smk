rule nextflow_nfcore_methylseq:
    input:
        csv = "src/sample_sheet.csv"
    log:
        log = "out/nextflow/nfcore_methylseq/log"
    shell:
        """
        ~/nextflow run nf-core/methylseq \
            --input {input.csv} \
            --outdir `dirname {log}` \
            --genome GRCh37 \
            -profile docker > {log}
        """