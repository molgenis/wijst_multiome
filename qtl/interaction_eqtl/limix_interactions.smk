from pathlib import Path
import re
import os
CHROM = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']

# read configuration
configfile: "./mo_interaction_template.yaml"
includeDir = config["top_dir"]
celltypes=config["celltypes"]
print(f"\n\nRunning pipeline for {celltypes[0]} cells\n")
scripts_folder = config['script_folder']
wp3_image_loc = config['wp3_image_loc']
limix_image_loc = config['limix_image_loc']
limix_loc = config['limix_loc']
genotype_loc = config['genotype_loc']

# make the limix path
limix_path="singularity exec --bind "+includeDir+" "+limix_image_loc+" python "+limix_loc"

##QTL mapping variables.
phenotypeFile = (config['phenotype_loc'] + config['phenotype_prepend'] + '%s' + config['phenotype_append']) % celltypes[0].split("/")[1]    # use {ct} to indicate celltype
# genotypes as split by chromosome and in bgen format
genotypeFile= config["genotype_loc"] + config['genotype_prepend']  + '{chrom}' + config['genotype_append'] # using {chrom} if genotype is splitted by chromosome
# covariate files are per cell type
covariateFile= (config['covariates_loc'] + config['covariates_prepend'] + '%s' + config['covariates_append']) % celltypes[0].split("/")[1]
# each cell type has its own output folder
outputFolder=(config["out_folder"]+ '%s')  % celltypes[0].split("/")[1]
kinshipFile= config["kinship_loc"]
chunkFile = config['chunking_loc']
# NOTE: this is currently not celltype specific
annoFile = config['limix_annotation_loc']+ config['limix_annotation_prepend'] + config['limix_annotation_append']
sampleMappingFile = config['smf_loc']
# get variant-feature file for each cell type
featureVariantFilterFile = (config['variant_feature_confinement_loc'] + config['variant_feature_confinement_prepend'] + '%s' + config['variant_feature_confinement_append']) % celltypes[0].split("/")[1]

# perform chunked analysis based on the chunks in the chunking file
chunk_chrom, chunk_start, chunk_end=[], [], []
with open(chunkFile) as fp:
    for line in fp:
        re_match=re.match(r"([A-Za-z0-9]+):(\d+)-(\d+)", line.strip())
        chunk_chrom.append(re_match[1])
        chunk_start.append(re_match[2])
        chunk_end.append(re_match[3])

# expand to get all chunks, which will have a .finished file when done
qtlChunks=expand(outputFolder/"{iet}"/"qtl"/"{chrom}_{start}_{end}.finished", zip, chrom=chunk_chrom, start=chunk_start, end=chunk_end, allow_missing=True)


rule all:
    input:
        expand(qtlChunks,iet=config["interaction_terms"]),
        expand(outputFolder/"{ct}/{iet}/iqtl_results_all.txt.gz", ct=celltypes, L2=config["interaction_terms"])
        
    output:
        touch(expand(outputFolder/"{iet}"/"done.txt", iet=config["interaction_terms"]))

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
        touch(outputFolder/'{iet}/qtl/{chrom}_{start}_{end}.finished')
    priority:10
    params:
        od = str(outputFolder/"{iet}"/"qtl")+"/",
        gen  = genotypeFile,
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
        outputFolder/"{ct}/{iet}/iqtl_results_all.txt.gz"
    params:
        idir = str(outputFolder/"{ct}/{iet}/qtl")+"/",
        odir = str(outputFolder/"{ct}/{iet}")+"/"
    shell:
        (limix_path + "post_processing/minimal_interaction_postprocess.py "
            " -id {params.idir} "
            " -od {params.odir} "
            "-sfo -wc ")