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

##############################################
#processing
##############################################

processed_bam = expand('{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/deam/{indexlibid}.uniq.L35MQ25.deam53x3.bam',
                             indexlibid=INDEXLIBID,
                             probeset=PROBESET,
                             seqrun=[get_seqrun(wildcards) for _ in INDEXLIBID])
print(processed_bam)

rule processing:
    input: processed_bam
    run:
        print('hey! Running target filtering, deamination, quality filtering and duplicate removal :)')
        pass


rule filter_bam_by_control_sites:
    input:
        control_sites=get_bed,
        bam=mapped_bams,
    output:
        sites_bam="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/{indexlibid}.bam",
    threads: 1
    conda: "envs/processing.yaml"
    shell: """
        bedtools intersect -a {input.bam} -b {input.control_sites} > {output.sites_bam} 
    """

rule uniq_q25_l35:
    input:
        bam = "{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/{indexlibid}.bam",
    output:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/{indexlibid}.uniq.L35MQ25.bam",
    shell:
            """
            bam-rmdup -o {output.bam} -q 25 -l 35 {input.bam}
            """

rule deam_filter:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/{indexlibid}.uniq.L35MQ25.bam"
    output: bam=processed_bam,
        stats="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/deam/{indexlibid}.uniq.L35MQ25.deam53x3.stats"
    conda: "envs/processing.yaml"
    shell: """
    /home/mmeyer/perlscripts/solexa/analysis/filterBAM.pl  -p5 0,1,2 -p3 0,-1,-2 -suffix deam53x3 {input.bam} &> {output.stats}
    """

samtools view -h {input.bam} | python pmdtools.0.60.py --customterminus 0,-1 --header | samtools view -Sb - > {
    output.bam}
--terminal --stats

##############################################
#kraken
##############################################
rule kraken:
    input:
        "{base}.byread",
        "{base}.kraken_spc"
    run:
        print('hey! Running Kraken')
        pass

rule bam_to_fasta:
    input:
        split=get_bam,
        bam="{base}/{base_name}.bam"
    output:
        bam="{base}/{base_name}.bam.fa.gz",
        split="{base}/kraken_split/{base_name}.bam.fa.gz"
    threads: 1
    shell: """
     bam2fastx -a -A -Q {input.bam} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.bam},
     bam2fastx -a -A -Q {input.split} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.split}
    """

rule generic_kraken:
    input: fa="{base}.fa.gz"
    output: kraken="{base}.kraken",
    threads: 1
    shell: """

     if [ $(hostname) != "bionc13" ] ; then 
         echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
         sleep 30
         stop
     fi

    echo THIS FAILS IF THE KRAKEN DB IS NOT COPIED - fix it

     nproc=$(time kraken --threads {threads} --db /mnt/ramdisk/refseqReleaseKraken \
     --output {output.kraken} {input.fa} 2>&1 | grep 'processed' | cut -f1 -d' ')
     echo "count seqs {input.fa}"
     nseq=$(gunzip -c {input.fa} | grep -c -e '>' -e '@') || echo "no seqs found $nseq"

     if [ ! $nseq -eq $nproc ] ; then echo KRAKEN RUN FAILED; stop; fi
    """

rule generic_kraken_summary:
    input: kraken="{base}.kraken"
    output: phylo="{base}.kraken_phylo"
    output:

threads: 1
shell: """

    if [ $(hostname) != "bionc13" ] ; then 
         echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
         sleep 30
         stop
    fi

    echo 'making report..'
    ## translate file has info on each read (but not in an easy-to-digest way)
    ## phylo file has a summary for each level in the taxonomy
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} > {output.phylo}

    """

rule generic_kraken_translate:
    input: kraken="{base}.kraken"
    output: translate="{base}.translate"
    threads: 1
    shell: """

    if [ $(hostname) != "bionc13" ] ; then 
         echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
         sleep 30
         stop
    fi

    echo 'making report..'
    ## translate file has info on each read (but not in an easy-to-digest way)
    ## phylo file has a summary for each level in the taxonomy
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} --translate {output.translate} > /dev/null

    """

rule generic_kraken_byread:
    input: translate="{base}.translate"
    output: byread="{base}.byread"
    threads: 1
    shell: """

     echo 'parsing translate file..'
     python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/parse_kraken_translate_file.py --translate-file {input.translate} > {output.byread}

    """

### summarize the major species groups in the kraken phylogeny file

rule generic_kraken_spc_summary:
    input: spc_krak="{base}.kraken_phylo",
        spc_groups="/mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/data/major_species_groups.txt"
    output: spc_summary="{base}.kraken_spc"
    shell: """

    grep -f {input.spc_groups} {input.spc_krak} | awk '{{print $6,$2}}' > {output.spc_summary}

    n=$(awk '{{s += $2}} END {{print s}}' {output.spc_summary})
    echo "$n : " {input.spc_krak}

    m=$(grep '  Mammalia$' {input.spc_krak} | awk '{{print $2-'$n'}}')
    echo "$n : $m : " {input.spc_krak}

    b=$((grep '  Bacteria$' {input.spc_krak} || echo 0 0) | awk '{{print $2}}')
    echo "$n : $m : $b" {input.spc_krak}

    s=$((grep '  Sauropsida$' {input.spc_krak} || echo 0 0) | awk '{{print $2}}')
    echo "$n : $m : $b : $s :" {input.spc_krak}

    r=$(grep '	root$' {input.spc_krak} | awk '{{print $2-'$n'-'$m'-'$b'-'$s'}}')
    echo "$n : $m : $b : $s : $r :" {input.spc_krak}

    u=$(grep '	unclassified$' {input.spc_krak} | awk '{{print $2}}')
    echo "$n : $m : $b : $s : $r : $u : " {input.spc_krak}

    echo Mammalia $m >> {output.spc_summary}
    echo Bacteria $b >> {output.spc_summary}
    echo Sauropsida $s >> {output.spc_summary}
    echo root $r >> {output.spc_summary}
    echo unclassified $u >> {output.spc_summary}

"""

