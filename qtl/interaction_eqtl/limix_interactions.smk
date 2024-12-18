from pathlib import Path
import re
import os
CHROM = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']

# read configuration
configfile: "./mo_interaction_template.yaml"
includeDir = config["top_dir"]
celltypes=config["celltypes"][0]
print(f"\n\nRunning pipeline for {celltypes} cells\n")
scripts_folder = config['script_folder']
wp3_image_loc = config['wp3_image_loc']
limix_image_loc = config['limix_image_loc']
limix_loc = config['limix_loc']
genotype_loc = config['genotype_loc']

# make the limix path
limix_path="singularity exec --bind "+includeDir+" "+limix_image_loc+" python "+limix_loc

##QTL mapping variables.
# phenotype data
phenotypeFile  = config['phenotype_loc']
if config['phenotype_prepend'] is None:
    formatted_path = f"{phenotypeFile}%s{config['phenotype_append']}" % celltypes
else:
    formatted_path = f"{phenotypeFile}{config['phenotype_prepend']}%s{config['phenotype_append']}" % celltypes
phenotypeFile = Path(formatted_path)

# genotype data
# genotype_prepend = config['genotype_prepend']
# if genotype_prepend is None:
#     genotype_prepend = ''
genotypeFile = config["genotype_loc"]
if config['genotype_prepend'] is None:
    formatted_path = ''.join([genotypeFile, '{chrom}'])  # using {chrom} if genotype is split by chromosome
else:
    formatted_path = ''.join([genotypeFile, config['genotype_prepend'], '{chrom}'])  # using {chrom} if genotype is split by chromosome
genotypeFile = formatted_path


# covariates
covariateFile = config['covariates_loc']
if config['covariates_prepend'] is None:
    formatted_path = f"{covariateFile}%s{config['covariates_append']}" % celltypes
else:
    formatted_path = f"{covariateFile}{config['covariates_prepend']}%s{config['covariates_append']}" % celltypes
covariateFile = Path(formatted_path)

outputFolder=(config["out_folder"] + '%s')  % celltypes
kinshipFile= config["kinship_loc"]
chunkFile = config['chunking_file_loc']
# NOTE: this is currently not celltype specific
annoFile = config['limix_annotation_loc']
if config['limix_annotation_loc'] is None:
    annoFile = Path(annoFile, config['limix_annotation_append'])
else:
    annoFile = Path(annoFile, config['limix_annotation_prepend'] + config['limix_annotation_append'])

# sample mapping file
sampleMappingFile = config['smf_loc']
# get variant-feature file for each cell type
featureVariantFilterFile = config['variant_feature_confinement_loc']
if config['variant_feature_confinement_prepend'] is None:
    formatted_path = f"{featureVariantFilterFile}%s{config['variant_feature_confinement_append']}" % celltypes
else:
    formatted_path = f"{featureVariantFilterFile}{config['variant_feature_confinement_prepend']}%s{config['variant_feature_confinement_append']}" % celltypes
featureVariantFilterFile = Path(formatted_path)

# perform chunked analysis based on the chunks in the chunking file
chunk_chrom, chunk_start, chunk_end=[], [], []
with open(chunkFile) as fp:
    for line in fp:
        re_match=re.match(r"([A-Za-z0-9]+):(\d+)-(\d+)", line.strip())
        chunk_chrom.append(re_match[1])
        chunk_start.append(re_match[2])
        chunk_end.append(re_match[3])

# expand to get all chunks, which will have a .finished file when done
qtlChunks=expand(Path(outputFolder, '{iet}', "qtl", "{chrom}_{start}_{end}.finished"), zip, chrom=chunk_chrom, start=chunk_start, end=chunk_end, allow_missing=True)

rule all:
    input:
        expand(qtlChunks,iet=config["interaction_terms"]),
        expand(Path(outputFolder,"{iet}","iqtl_results_all.txt.gz"), ct=celltypes, iet=config["interaction_terms"])

    output:
        touch(expand(Path(outputFolder, "{iet}", "done.txt"), iet=config["interaction_terms"]))

# run QTL for a chunk
rule run_qtl_mapping:
    input:
        af = annoFile,
        pf = phenotypeFile,
        smf = sampleMappingFile,
        cf = covariateFile,
        kf = kinshipFile,
        fvf = featureVariantFilterFile
        #rf = config["randomeff_files"] if config["randomeff_files"]!='' else []
    output:
        #touch(outputFolder/"{ct}"/"qtl"/"{chrom}_{start}_{end}.finished")
        touch(Path(outputFolder, '{iet}', 'qtl', '{chrom}_{start}_{end}.finished'))
    priority:10
    params:
        od = str(Path(outputFolder, "{iet}", "qtl"))+"/",
        #gen = lambda wildcards: Path(f"{config['genotype_loc']}{genotype_prepend}{wildcards.chrom}{config['genotype_append']}"),
        gen = genotypeFile,
        np = config["numberOfPermutations"],
        maf = config["minorAlleleFrequency"],
        hwe = config["hardyWeinbergCutoff"],
        w = config["windowSize"]
    shell:
        (limix_path + "run_interaction_QTL_analysis.py "
            " --bgen {params.gen} "
            " -af {input.af} "
            " -cf {input.cf} "
            " -pf {input.pf} "
            " -rf {input.kf} "
            " -fvf {input.fvf} "
            " -smf {input.smf} "
            " -od {params.od}/ "
            " -gr {wildcards.chrom}:{wildcards.start}-{wildcards.end} "
            " -np {params.np} "
            " -maf {params.maf} "
            " -t -gm gaussnorm "
            " -it {wildcards.iet} "
            " -w {params.w} "
            " -hwe {params.hwe} "
            " -rs .8 ")

# post-processing on merged chunks
rule all_qtl:
    input:
        qtlChunks
    output:
        Path(outputFolder, "{iet}", "iqtl_results_all.txt.gz")
    params:
        idir = str(Path(outputFolder, "{iet}", "qtl"))+"/",
        odir = str(Path(outputFolder, "{iet}"))+"/"
    shell:
        (limix_path + "post_processing/minimal_interaction_postprocess.py "
            " -id {params.idir} "
            " -od {params.odir} "
            "-sfo -wc ")
