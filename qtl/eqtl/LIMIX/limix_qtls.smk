from pathlib import Path
import re
import os

configfile: "./qtl_template.yaml"


includeDir = config["top_dir"]
celltypes=config["celltypes"]
outputFolder=Path(config["out_folder"])
scripts_folder = config['script_folder']
# path to images
limix_image_loc = config['limix_image_loc']
wp3_image_loc = config['wp3_image_loc']
# paste together command to image
limix_path="singularity exec --bind "+includeDir+" "+limix_image_loc+" python /limix_qtl/Limix_QTL/"


##QTL mapping variables.
# phenotype data
phenotypeFile  = config['phenotype_loc']
if config['phenotype_prepend'] is None:
    formatted_path = f"{phenotypeFile}{ct}{config['phenotype_append']}"
else:
    formatted_path = f"{phenotypeFile}{config['phenotype_prepend']}{ct}{config['phenotype_append']}"
phenotypeFile = Path(formatted_path)

# genotype data
genotypeFolder = ["genotype_loc"]
genotypeFile = config["genotype_loc"]
if config['genotype_prepend'] is None:
    formatted_path = ''.join([genotypeFile, '{chrom}'])  # using {chrom} if genotype is split by chromosome
else:
    formatted_path = ''.join([genotypeFile, config['genotype_prepend'], '{chrom}'])  # using {chrom} if genotype is split by chromosome
genotypeFile = formatted_path


# covariates
covariateFile = config['covariates_loc']
if config['covariates_prepend'] is None:
    formatted_path = f"{covariateFile}{ct}{config['covariates_append']}"
else:
    formatted_path = f"{covariateFile}{config['covariates_prepend']}{ct}{config['covariates_append']}"
covariateFile = Path(formatted_path)

outputFolder=Path(config["out_folder"])
kinshipFile= config["kinship_loc"]
chunkFile = config['chunking_file_loc']
# NOTE: this is currently not celltype specific
annoFile = config['limix_annotation_loc']
sampleMappingFile = config['smf_loc']

chunk_chrom, chunk_start, chunk_end=[], [], []
with open(chunkFile) as fp:
    for line in fp:
        re_match=re.match(r"([A-Za-z0-9]+):(\d+)-(\d+)", line.strip())
        chunk_chrom.append(re_match[1])
        chunk_start.append(re_match[2])
        chunk_end.append(re_match[3])

qtlChunks=expand(outputFolder/"{ct}"/"qtl"/"{chrom}_{start}_{end}.finished", zip, chrom=chunk_chrom, start=chunk_start, end=chunk_end, allow_missing=True)

wildcard_constraints:
    end = "([0-9]+)",
    start = "(-?[0-9]+)",
    num = "([0-9]+)",
    chrom = "[0-9]{1,2}|X|Y|MT"


rule all:
    input:
        expand(outputFolder/"{ct}"/"top_qtl_results_all.txt.gz", ct=celltypes),
        expand(outputFolder/"{ct}"/"qtl_results_all.txt.gz", ct=celltypes),
        expand(outputFolder/"{ct}"/"qtl.h5.tgz", ct=celltypes),
        expand(outputFolder/"{ct}"/"qtl.annotation.tgz", ct=celltypes),
        expand(outputFolder/"{ct}"/"qtl.permutations.tgz", ct=celltypes),
#        outputFolder/"LDMatrices.tgz"
    output:
        touch(expand(outputFolder/"{ct}"/"done.txt", ct=celltypes))


rule run_qtl_mapping:
    input:
        af = annoFile,
        pf = phenotypeFile,
        smf = sampleMappingFile,
        cf = covariateFile,
        #rf = config["randomeff_files"] if config["randomeff_files"]!='' else []
    output:
        touch(outputFolder/"{ct}"/"qtl"/"{chrom}_{start}_{end}.finished")
    resources:
        memory = "12000",
        time = "9:59:00"
    params:
        cmd = limix_path,
        od = str(outputFolder/"{ct}"/"qtl")+"/",
        gen  = genotypeFile,
        np = config["numberOfPermutations"],
        maf = config["minorAlleleFrequency"],
        hwe = config["hardyWeinbergCutoff"],
        w = config["windowSize"]
    shell: """{params.cmd}run_QTL_analysis_metaAnalysis.py \
             --bgen {params.gen} \
             -af {input.af} \
             -cf {input.cf} \
             -pf {input.pf} \
             -smf {input.smf} \
             -od {params.od} \
             -gr {wildcards.chrom}:{wildcards.start}-{wildcards.end} \
             -np {params.np} \
             -maf {params.maf} \
             -c -gm gaussnorm \
             -w {params.w} \
             -hwe {params.hwe} \
             -rs .8 """

rule top_feature:
    resources:
        memory = "24000",
        time = "12:59:00"
    input:
        qtlChunks
    output:
        temp(outputFolder/"{ct}"/"top_qtl_results_all.txt")
    priority:10
    params:
        idir = str(outputFolder/"{ct}"/"qtl")+"/",
        odir = str(outputFolder/"{ct}")+"/"
    shell:
        (limix_path + "post_processing/minimal_postprocess.py "
            " -id {params.idir} "
            " -od {params.odir} "
            "-tfb "
            "-sfo ")

rule all_qtl:
    resources:
        memory = "64000",
        time = "5:59:00"
    input:
        qtlChunks
    output:
        outputFolder/"{ct}"/"qtl_results_all.txt.gz"
    priority:10
    params:
        idir = str(outputFolder/"{ct}"/"qtl")+"/",
        odir = str(outputFolder/"{ct}")+"/"
    shell:
        (limix_path + "post_processing/minimal_postprocess.py "
            " -id {params.idir} "
            " -od {params.odir} "
            "-wc "
            "-sfo ")

rule compress_qtl:
    resources:
        memory = "64000",
        time = "5:59:00"
    input:
        topQtl = outputFolder/"{ct}"/"top_qtl_results_all.txt",
        qtlChunks = qtlChunks
    output:
        h5=outputFolder/"{ct}"/"qtl.h5.tgz",
        anno=outputFolder/"{ct}"/"qtl.annotation.tgz",
        perm=outputFolder/"{ct}"/"qtl.permutations.tgz",
        gzTopQtl=outputFolder/"{ct}"/"top_qtl_results_all.txt.gz"
    params:
        idir = str(outputFolder/"{ct}"/"qtl")
    shell: "tar -czf {output.h5} {params.idir}/*.h5 && find {params.idir} | grep '\.pickle.gz$' > {params.idir}/files.txt &&  tar cf {output.perm} -T {params.idir}/files.txt && tar -cf {output.anno} {params.idir}/*.txt.gz && gzip {input.topQtl} "
    ##&& rm -rf {params.idir}

rule make_temporary_files:
    resources:
        memory = "64000",
        time = "5:59:00"
    input:
        bgen_file=genotypeFile+".bgen"
    output:
        temp(genotypeFile+"{chrom}.bgen.z"),
        temp(genotypeFile+"{chrom}.bgen.sample"),
        temp(genotypeFile+"{chrom}.bgen_master.txt")
    shell:
        """
        singularity exec --bind {includeDir} {wp3_image_loc} python {scripts_folder}/make_z.py {input.bgen_file} {genotypeFolder}
        """

rule create_bdose_file_by_chr:
    resources:
        memory = "64000",
        time = "5:59:00"
    input:
        genotypeFile+"{chrom}.bgen",
        genotypeFile+"{chrom}.bgen.bgi",
        genotypeFile+"{chrom}.bgen.z",
        genotypeFile+"{chrom}.bgen.sample",
        genotypeFile+"{chrom}.bgen_master.txt"
    output:
        temp(genotypeFile+"{chrom}.bgen.bdose"),
        temp(genotypeFile+"{chrom}.bgen.bdose.bdose.tmp0"),
        temp(genotypeFile+"{chrom}.bgen.bdose.meta.tmp0")
    shell:
        """
        singularity exec --bind {includeDir} {wp3_image_loc} /tools/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 --in-files {genotypeFile}{wildcards.chrom}.bgen_master.txt --write-bdose --bdose-version 1.1
        """
