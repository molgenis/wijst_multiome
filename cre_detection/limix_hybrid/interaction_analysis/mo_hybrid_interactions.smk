import os

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


############################################
# define chunking scheme
############################################

WC = glob_wildcards(f"{CHUNK_BASE}/chr{{chrom}}-{{start}}-{{end}}")

def chunk_id(wc):
    return f"chr{wc.chrom}-{wc.start}-{wc.end}"

def chunk_dir(wc):
    return f"{CHUNK_BASE}/{chunk_id(wc)}"

def chunk_outdir(wc):
    return f"{RESULTS_BASE}/{chunk_id(wc)}"


############################################
# extract chunk info and define rules for each chunk
############################################

rule all:
    input:
        expand(f"{RESULTS_BASE}/chr{{chrom}}-{{start}}-{{end}}/result.tsv.gz",
               chrom=WC.chrom, start=WC.start, end=WC.end)

def chunk_expression(wc):
    return os.path.join(chunk_dir(wc), config["expression_filename"])

def chunk_accessibility(wc):
    return os.path.join(chunk_dir(wc), config["accessibility_filename"])


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
        # Create per-chrom prefix (no extension) for R script
        genotype_prefix = GENO_TMPL,
        rscript = "mo_hybrid_cre_interaction.R"
    threads: 2
    shell:
        r"""
        mkdir -p "{params.outdir}"

        # boolean flags
        EX = "--expression_gausnorm" if "{config[expression_gausnorm]}" == "True" else ""
        AC = "--accessibility_gausnorm" if "{config[accessibility_gausnorm]}" == "True" else ""

        Rscript "{params.rscript}" \
            --in "{chunk_dir(wildcards)}" \
            --out "{params.outdir}" \
            --confinement "{input.confinement}" \
            --smf_loc "{input.smf}" \
            --covariates_file "{input.covariates}" \
            --fixed_effects "{params.fixed_effects}" \
            --random_effects "{params.random_effects}" \
            --interaction_terms "{params.interaction_terms}" \
            --barcode_column "{params.barcode_column}" \
            --genotype_loc "{params.genotype_prefix.format(chrom=wildcards.chrom)}" \
            $EX \
            $AC

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
            end=WC.end
        )
    output:
        merged = "results/merged/all_results.tsv.gz"
    run:
        import os, gzip, shutil, tempfile

        infiles = list(input)

        # Filter out files that are missing or zero-size
        nonempty = [f for f in infiles if os.path.exists(f) and os.path.getsize(f) > 0]

        # Temporary uncompressed file for merging
        tmp_path = tempfile.NamedTemporaryFile(delete=False).name

        wrote_header = False

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

        # If nothing was written at all (all chunks empty), create an empty gz
        if not wrote_header:
            with gzip.open(output.merged, "wt") as fout:
                pass
        else:
            # Compress merged file
            with open(tmp_path, "rb") as fin, gzip.open(output.merged, "wb") as fout:
                shutil.copyfileobj(fin, fout)

        os.unlink(tmp_path)
