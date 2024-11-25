################################################################################
# Draft nuclear sediment pipeline
#
# Aurore GALTIER, 25/11/2024
################################################################################

from snakemake.utils import R
import pandas as pd
import glob

#where to run the script
workdir: "pipeline_path"

#to have errors output
if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")

#define project folder
project = config["project"]

#define summary folder
summary = config["summary"]

indexlibid_df = pd.read_table(config["samples_info"], comment="#").set_index(["indexlibid"], drop=False).sort_index()
INDEXLIBID = indexlibid_df.index.unique()

def get_seqrun(wildcards):
    return (indexlibid_df.loc[(wildcards.indexlibid), "seqrun"])

def get_bam(wildcards):
    return (indexlibid_df.loc[(wildcards.indexlibid), "bam"])


probeset_df = pd.read_table(config["probeset_info"], comment="#").set_index(["probeset"], drop=False).sort_index()
PROBESET = probeset_df.index.unique()

def get_ref(wildcards):
    return (probeset_df.loc[(wildcards.probeset), "path_to_ref"])


def get_bed(wildcards):
    return (probeset_df.loc[(wildcards.probeset), "path_to_bed"])



##############################################
#map_all
##############################################

mapped_bams = expand('{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/{indexlibid}.bam',
                             indexlibid=INDEXLIBID,
                             probeset=PROBESET,
                             seqrun=[get_seqrun(wildcards) for _ in INDEXLIBID])
print(mapped_bams)

rule map_all:
    input:
        bam = get_bam,
        ref = get_ref
    output: mapped_bams
    threads: 8
    conda: "envs/align.yaml"
    shell:
        """
        print('hey! mapping bam')
        bwa bam2bam -t {threads} -g {input.ref} -n 0.01 -o 2 -l 16500 --only-aligned -f {input.bam} | samtools sort -@30 -o {output}
        """

