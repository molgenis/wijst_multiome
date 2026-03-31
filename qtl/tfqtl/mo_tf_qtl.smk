######################################################
# Authors: Roy Oelen
# Name: mo_tf_qtl.smk
# Function: Snakemake pipeline for TF-QTL hybrid analysis
######################################################

import os
import re
import glob

############################################
# load default config
############################################

configfile: "config.yaml"

############################################
# configuration
############################################

RESULTS_BASE     = config["results_base"]
GENO_TMPL        = config["genotype_prefix_template"]
RSCRIPT_LOC      = config["script_loc"]
MERGE_SCRIPT_LOC = config["merge_script_loc"]
R_COMMAND        = config["r_command"]
CHUNK_INPUT_DIR  = config["chunk_input_dir"]

############################################
# helpers
############################################

def chrom_from_chunk(chunk: str) -> str:
    return chunk.replace("chr", "")

############################################
# define chunks based on genotype files
############################################

CHUNKS = []
geno_dir = os.path.dirname(GENO_TMPL.format(chrom=""))
geno_pat = re.compile(r".*_chr([^/]+)\.bed$")

for bed_file in glob.glob(os.path.join(geno_dir, "*.bed")):
    m = geno_pat.match(bed_file)
    if m:
        CHUNKS.append(f"chr{m.group(1)}")

CHUNKS = sorted(set(CHUNKS))

if not CHUNKS:
    print(f"[WARN] No genotype-based chunks found in: {geno_dir}")

############################################
# default target
############################################

rule all:
    input:
        expand(f"{RESULTS_BASE}/{{chunk}}/result.tsv.gz", chunk=CHUNKS),
        f"{RESULTS_BASE}/merged/results_fdr.tsv.gz"

############################################
# run each chunk
############################################

rule run_interaction:
    input:
        confinement = config["confinement"],
        smf = config["smf"],
        covariates = config["covariates"],
        tf_file = config["tf_file"],
        genotype_bed = lambda wc: GENO_TMPL.format(chrom=chrom_from_chunk(wc.chunk)) + ".bed",
        genotype_bim = lambda wc: GENO_TMPL.format(chrom=chrom_from_chunk(wc.chunk)) + ".bim",
        genotype_fam = lambda wc: GENO_TMPL.format(chrom=chrom_from_chunk(wc.chunk)) + ".fam"
    output:
        tsv = f"{RESULTS_BASE}/{{chunk}}/result.tsv.gz"
    params:
        in_dir   = CHUNK_INPUT_DIR,
        outdir   = lambda wc: f"{RESULTS_BASE}/{wc.chunk}",
        fixed_effects  = config["fixed_effects"],
        random_effects = config["random_effects"],
        barcode_column = config["barcode_column"],
        genotype_prefix = lambda wc: GENO_TMPL.format(chrom=chrom_from_chunk(wc.chunk)),
        tf_gausnorm = "--tf_gausnorm" if bool(config.get("tf_gausnorm", False)) else "",
        rscript = RSCRIPT_LOC,
        rcmd    = R_COMMAND
    threads: 2
    shell:
        r"""
        mkdir -p "{params.outdir}"

        {params.rcmd} "{params.rscript}" \
            --in "{params.in_dir}" \
            --out "{params.outdir}" \
            --confinement "{input.confinement}" \
            --smf_loc "{input.smf}" \
            --tf_file "{input.tf_file}" \
            --covariates_file "{input.covariates}" \
            --fixed_effects "{params.fixed_effects}" \
            --random_effects "{params.random_effects}" \
            --barcode_column "{params.barcode_column}" \
            --genotype_loc "{params.genotype_prefix}" \
            {params.tf_gausnorm}

        test -s "{output.tsv}"
        """

############################################
# merge results
############################################

rule merge_results:
    input:
        expand(f"{RESULTS_BASE}/{{chunk}}/result.tsv.gz", chunk=CHUNKS)
    output:
        merged = f"{RESULTS_BASE}/merged/results_fdr.tsv.gz"
    params:
        in_dir = RESULTS_BASE,
        rscript = MERGE_SCRIPT_LOC,
        rcmd    = R_COMMAND
    shell:
        r"""
        mkdir -p "$(dirname {output.merged})"

        {params.rcmd} "{params.rscript}" \
            --in "{params.in_dir}" \
            --out "{output.merged}"
        """
