######################################################
# Snakemake pipeline for hybrid interaction analysis #
######################################################

import os
import re


############################################
# load default config
############################################

configfile: "config.yaml"


############################################
# configuration and path setup
############################################

CHUNK_BASE   = config["chunk_base"]
RESULTS_BASE = config["results_base"]
GENO_TMPL    = config["genotype_prefix_template"]

RSCRIPT_LOC  = config["script_loc"]
R_COMMAND    = config["r_command"]


############################################
# to store chunk format
############################################

class WCObject:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
    
    def __getitem__(self, key):
        return getattr(self, key)


############################################
# define helpers
############################################

# define chunk_id (you were calling it, but it wasn't defined)
def chunk_id(wc):
    return f"chr{wc.chrom}-{wc.start}-{wc.end}"

# define chunk_dir (also used but missing)
def chunk_dir(wc):
    return f"{CHUNK_BASE}/{chunk_id(wc)}"

# autodetect absolute vs relative paths
def resolve_rel_or_abs(base_dir, maybe_rel_path):
    """Return absolute path unchanged; join relative path to base_dir."""
    if os.path.isabs(maybe_rel_path):
        return maybe_rel_path
    return os.path.join(base_dir, maybe_rel_path)

def chunk_outdir(wc):
    return f"{RESULTS_BASE}/{chunk_id(wc)}"

def chunk_expression(wc):
    return resolve_rel_or_abs(chunk_dir(wc), config["expression_filename"])

def chunk_accessibility(wc):
    return resolve_rel_or_abs(chunk_dir(wc), config["accessibility_filename"])


############################################
# define chunking scheme
############################################

# grab directories that match chr{chrom}-{start}-{end}
pat = re.compile(r"^chr([^-/]+)-(\d+)-(\d+)$")

# create list of cnunk diretoreis
CHUNK_DIRS = []
# check each entry
with os.scandir(CHUNK_BASE) as it:
    for entry in it:
        # check if even a directory
        if entry.is_dir():
            # then check if matches the pattern
            m = pat.match(entry.name)
            if m:
                # then add to the list
                CHUNK_DIRS.append(entry.name)

# sort so that it is easier to check
CHUNK_DIRS.sort()

# extract all the values
CHROMS = [d.split("-")[0].replace("chr", "") for d in CHUNK_DIRS]
STARTS = [d.split("-")[1] for d in CHUNK_DIRS]
ENDS   = [d.split("-")[2] for d in CHUNK_DIRS]

# make object
WC = WCObject(chrom=CHROMS, start=STARTS, end=ENDS)

rule all:
    input:
        # each chunk
        expand(f"{RESULTS_BASE}/chr{{chrom}}-{{start}}-{{end}}/result.tsv.gz",
               chrom=WC.chrom, start=WC.start, end=WC.end, zip = True),
        # merged chunk results
        f"{RESULTS_BASE}/all_results.tsv.gz"


############################################
# run each chunk
############################################

rule run_interaction:
    input:
        expression = chunk_expression,
        accessibility = chunk_accessibility,
        confinement = config["confinement"],
        smf = config["smf"],
        covariates = config["covariates"],
        # OPTIONAL stricter tracking of PLINK files:
        genotype_bed = lambda wc: GENO_TMPL.format(chrom=wc.chrom) + ".bed",
        genotype_bim = lambda wc: GENO_TMPL.format(chrom=wc.chrom) + ".bim",
        genotype_fam = lambda wc: GENO_TMPL.format(chrom=wc.chrom) + ".fam"
    output:
        tsv = f"{RESULTS_BASE}/chr{{chrom}}-{{start}}-{{end}}/result.tsv.gz"
    params:
        outdir = chunk_outdir,
        fixed_effects = config["fixed_effects"],
        random_effects = config["random_effects"],
        interaction_terms = config["interaction_terms"],
        barcode_column = config["barcode_column"],
        genotype_prefix = lambda wc: GENO_TMPL.format(chrom=wc.chrom),
        EXPR_FLAG = "--expression_gausnorm" if bool(config["expression_gausnorm"]) else "",
        ACC_FLAG  = "--accessibility_gausnorm" if bool(config["accessibility_gausnorm"]) else "",
        rscript = RSCRIPT_LOC,
        rcmd = R_COMMAND,
        in_dir = lambda wc: chunk_dir(wc)
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
            --interaction_terms "{params.interaction_terms}" \
            --barcode_column "{params.barcode_column}" \
            --genotype_loc "{params.genotype_prefix}" \
            {params.EXPR_FLAG} \
            {params.ACC_FLAG}

        # Verify output exists
        test -s "{output.tsv}"
        """


############################################
# write merged results
############################################

rule merge_results:
    input:
        # discover all result files produced by per-chunk rules
        expand(
            f"{RESULTS_BASE}/chr{{chrom}}-{{start}}-{{end}}/result.tsv.gz",
            chrom=WC.chrom,
            start=WC.start,
            end=WC.end, 
			zip = True
        )
    output:
        merged = f"{RESULTS_BASE}/all_results.tsv.gz"
    run:
        import os, gzip, shutil, tempfile

        infiles = list(input)

        # remove fiFles that are missing or zero-size
        nonempty = [f for f in infiles if os.path.exists(f) and os.path.getsize(f) > 0]

        # check if output dir exists
        os.makedirs(os.path.dirname(output.merged), exist_ok=True)

        # uncompress before the merge
        tmp_path = tempfile.NamedTemporaryFile(delete=False).name
        # keep track of whether we wrote header
        wrote_header = False
        # open file to write
        with open(tmp_path, "wt") as tmp:
            for f in nonempty:
                with gzip.open(f, "rt") as fin:
                    for i, line in enumerate(fin):
                        if i == 0:
                            # Only write the first header
                            if not wrote_header:
                                tmp.write(line)
                                wrote_header = True
                            # Skip header of later files
                        else:
                            tmp.write(line)

        # create empty file if nothing was done
        if not wrote_header:
            with gzip.open(output.merged, "wt") as fout:
                pass
        else:
            # compress merged file
            with open(tmp_path, "rb") as fin, gzip.open(output.merged, "wb") as fout:
                shutil.copyfileobj(fin, fout)

        os.unlink(tmp_path)
