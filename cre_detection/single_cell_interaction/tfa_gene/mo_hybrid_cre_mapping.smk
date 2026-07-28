######################################################
# Authors: Roy Oelen
# Name: mo_hybrid_cre_mapping.smk
# Function: Snakemake pipeline for hybrid interaction analysis
######################################################

# for paths
import os
# for pattern matching
import re


############################################
# load default config
############################################

# default config file
configfile: "config.yaml"


############################################
# configuration and path setup
############################################

# grab variables from the configuration
CHUNK_BASE   = config["chunk_base"]
RESULTS_BASE = config["results_base"]
# also where R and the script are
RSCRIPT_LOC  = config["script_loc"]
MERGE_SCRIPT_LOC = config["merge_script_loc"]
R_COMMAND    = config["r_command"]


############################################
# helpers
############################################

def resolve_rel_or_abs(base_dir, maybe_rel_path):
    """Return absolute path unchanged; join relative path to base_dir."""
    if os.path.isabs(maybe_rel_path):
        return maybe_rel_path
    return os.path.join(base_dir, maybe_rel_path)

def chunk_dir_from_chunk(chunk: str) -> str:
    return os.path.join(CHUNK_BASE, chunk)

def chunk_outdir_from_chunk(chunk: str) -> str:
    return os.path.join(RESULTS_BASE, chunk)

def chrom_from_chunk(chunk: str) -> str:
    # chunk like "chr9-99361653-100377592" -> "9"
    return chunk.split("-")[0].replace("chr", "")

def is_set(val) -> bool:
    """True if val is provided (not None/empty after stripping)."""
    if val is None:
        return False
    s = str(val).strip()
    return s != "" and s.lower() != "none"

def optional_chunk_file(maybe_path, chunk):
    """
    Return a resolved (absolute) file path if provided in config (supports abs/rel);
    otherwise None so we can omit it from inputs and CLI.
    """
    if not is_set(maybe_path):
        return None
    return resolve_rel_or_abs(chunk_dir_from_chunk(chunk), maybe_path)


############################################
# define chunking scheme (discover only real folders)
############################################

# match directories like: chr{chrom}-{start}-{end}
_pat = re.compile(r"^chr([^-/]+)-(\d+)-(\d+)$")

# set the list of chunks we'll have
CHUNKS = []
# iterate the items in the directory
with os.scandir(CHUNK_BASE) as it:
    # check each entry
    for entry in it:
        # check if it is a directory and fits the pattern of a chunk
        if entry.is_dir() and _pat.match(entry.name):
            CHUNKS.append(entry.name)

# sort the chunks so it is easier to reproduce
CHUNKS.sort()

# warn if no chunks could be found
if not CHUNKS:
    print(f"[WARN] No chunk directories found under: {CHUNK_BASE}")


############################################
# default target
############################################

# expand through all of the chunks
rule all:
    input:
        # each chunk result
        expand(f"{RESULTS_BASE}/{{chunk}}/result.tsv.gz", chunk=CHUNKS),
        # merge results
        f"{RESULTS_BASE}/merged/results_fdr.tsv.gz"


############################################
# run each chunk
############################################

# the rule to run an interaction analysis
rule run_interaction:
    input:
        # grab all the expression and accessibility files
        expression = lambda wc: resolve_rel_or_abs(chunk_dir_from_chunk(wc.chunk), config["expression_filename"]),
        accessibility = lambda wc: resolve_rel_or_abs(chunk_dir_from_chunk(wc.chunk), config["accessibility_filename"]),
        # fetch confinement
        confinement = config["confinement"], 
        # fetch sample mapping
        smf = config["smf"], 
        # and the covariates
        covariates = config["covariates"]
    output:
        tsv = f"{RESULTS_BASE}/{{chunk}}/result.tsv.gz"
    params:
        in_dir   = lambda wc: chunk_dir_from_chunk(wc.chunk),
        outdir   = lambda wc: chunk_outdir_from_chunk(wc.chunk),
        fixed_effects     = config["fixed_effects"],
        random_effects    = config["random_effects"],
        barcode_column    = config["barcode_column"],
        # genotype prefix (precomputed in Python)
        genotype_prefix   = lambda wc: GENO_TMPL.format(chrom=chrom_from_chunk(wc.chunk)),
        # get the flags for expression or accessibility
        EXPR_FLAG = "--expression_gausnorm"    if bool(config.get("expression_gausnorm", False)) else "",
        ACC_FLAG  = "--accessibility_gausnorm" if bool(config.get("accessibility_gausnorm", False)) else "",
        EXPR_FLAG_BC = "--expression_boxcox"    if bool(config.get("expression_boxcox", False)) else "",
        ACC_FLAG_BC  = "--accessibility_boxcox" if bool(config.get("accessibility_boxcox", False)) else "",
        # optional CLI args for expression/accessibility (only if provided)
        expr_arg = lambda wc: (
            f'--expression_file "{optional_chunk_file(config.get("expression_filename"), wc.chunk)}"'
            if optional_chunk_file(config.get("expression_filename"), wc.chunk) else ""
        ),
        acc_arg = lambda wc: (
            f'--accessibility_file "{optional_chunk_file(config.get("accessibility_filename"), wc.chunk)}"'
            if optional_chunk_file(config.get("accessibility_filename"), wc.chunk) else ""
        ),
        # R invocation
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
            --covariates_file "{input.covariates}" \
            --fixed_effects "{params.fixed_effects}" \
            --random_effects "{params.random_effects}" \
            --barcode_column "{params.barcode_column}" \
            {params.EXPR_FLAG} \
            {params.ACC_FLAG} \
            {params.EXPR_FLAG_BC} \
            {params.ACC_FLAG_BC} \
            {params.expr_arg} \
            {params.acc_arg}


        # check if output exists (your R writes an empty gz when nothing to test)
        test -s "{output.tsv}"
        """


############################################
# write merged results using R script
############################################

rule merge_results:
    input:
        # all chunk result files (ensures rule waits for them)
        expand(f"{RESULTS_BASE}/{{chunk}}/result.tsv.gz", chunk=CHUNKS)
    output:
        merged = f"{RESULTS_BASE}/merged/results_fdr.tsv.gz"
    params:
        # directory containing all chunk subdirectories
        in_dir = RESULTS_BASE,
        rscript = MERGE_SCRIPT_LOC,
        rcmd = R_COMMAND
    shell:
        r"""
        mkdir -p "$(dirname {output.merged})"

        {params.rcmd} "{params.rscript}" \
            --in "{params.in_dir}" \
            --out "{output.merged}"
        """
