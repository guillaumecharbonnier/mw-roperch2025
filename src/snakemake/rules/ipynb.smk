rule ipynb_compute_differential_methylation:
    output:
        xlsx = "out/ipynb/compute_differential_methylation/{design}.xlsx",
    log:
        rds = "out/ipynb/compute_differential_methylation/{design}.rds",
        # optional path to the processed notebook
        notebook="out/ipynb/compute_differential_methylation/{design}.r.ipynb"
    conda:
        "../envs/jupyter_rkernel_diff_methyl.yaml"
    threads: 32 # Huge memory usage, instances of this rule should not run in parallel
    notebook:
        "../../ipynb/compute_differential_methylation.r.ipynb"

rule ipynb_index:
    input:
        expand(
            "out/ipynb/compute_differential_methylation/{design}.xlsx",
            design = [
        #         "ODG_080",
        #         "ODG_082",
        #         # "ODG_084",
        #         # "ODG_086",
        #         # "ODG_088",
        #         # "ODG_090",
        #         # "ODG_092",
        #         # "ODG_094",
        #         # "ODG_096",
        #         "ODG_098",
        #         "ODG_100",
        #         "ODG_102",
        #         "ODG_104",
        #         "ODG_106",
        #         "ODG_108",
        #         # Currently bugged since the addition of top_dml_dmr_merge
		# #"dev_paired",
        #         "dev_unpaired"
                "all_unpaired"
            ]
        )
    # output:
    #     "out/ipynb/index.py.html"
    log:
        notebook="out/ipynb/index.py.ipynb"
    conda:
        "../envs/jupyter_rkernel_diff_methyl.yaml"
    notebook:
        "../../ipynb/index.py.ipynb"


rule nbconvert_ipynb_index:
    input:
        ipynb = "out/ipynb/index.py.ipynb"
    output:
        "out/ipynb/index.py.html"
    conda:
        "../envs/jupyter_rkernel_diff_methyl.yaml"
    shell:
        "jupyter nbconvert {input.ipynb} --to html --no-prompt --TagRemovePreprocessor.remove_input_tags='remove-input'"
