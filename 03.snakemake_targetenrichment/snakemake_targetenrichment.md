---
title: PGx analysis using PacBio HiFi target enrichment workflow
authors: Xingliang LIU, Wilson CHENG, Tingting ZHU
date: Monday, November 27, 2023
duration: 60 minutes
---

### Learning Objectives:

* Understand key steps and main output of PacBio HiFi target enrichment workflow
* Understand how to run PacBio HiFi target enrichment workflow using Slurm job management system (JMS) on a high-performing computer (HPC)

### Downsampling

Download downsampled HiFi enrichment dataset from PacBio [DATASET](https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/).

```bash
mkdir -p ~/CUMED_BFX_workshop/03.snakemake_targetenrichment/data/cyp2d6/downsampled
cd ~/CUMED_BFX_workshop/03.snakemake_targetenrichment/data/cyp2d6/downsampled

curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA02016.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA07439.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA09301.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA12244.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA16654.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA16688.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17020.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17039.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17073.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17114.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17209.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17211.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17214.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17215.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17217.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17226.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17227.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17232.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17244.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17276.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17282.bam
curl -O https://downloads.pacbcloud.com/public/dataset/HiFi_amplicon_CYP2D6/downsampled/NA17300.bam
```

### Setting up run environment

```bash
mamba create -n snakemake_targetenrichment \
    --channel conda-forge \
    --channel bioconda \
    python=3.9 snakemake=7.32.4 pysam=0.22.0 lockfile=0.12.2 git=2.42.0 pbmarkdup=1.0.3 mosdepth=0.3.5 

# for lockfile
sudo apt update -y
sudo apt install procmail -y
```

### Preparing de-duplicated HiFi alignment

```bash
# 1. pbmarkdup mark PCR duplicates
# 2. pbmm2 HiFi alignment (sorted)

cd ~/CUMED_BFX_workshop/03.snakemake_targetenrichment/data/cyp2d6/downsampled

SMRT_ROOT=/home/xliu/softwares/pacbio/smrtlinkv13_ea_smrttools
for bam in ~/CUMED_BFX_workshop/03.snakemake_targetenrichment/data/cyp2d6/downsampled/*.bam; do
    echo $bam ...
    pbmarkdup --log-level INFO \
    -j 16 \
    --cross-library \
    $bam ${bam%.bam}.pbmarkdup.bam

    samtools index ${bam%.bam}.pbmarkdup.bam
    $SMRT_ROOT/smrtcmds/bin/pbindex ${bam%.bam}.pbmarkdup.bam

    ~/softwares/pacbio/smrtlinkv13_ea_smrttools/smrtcmds/bin/pbmm2 align \
    	--num-threads 16 \
    	--sort-memory 4G \
    	--preset HIFI \
    	--log-level INFO \
    	--sort \
        --bam-index BAI \
    	~/CUMED_BFX_workshop/01.wdl_humanwgs/dataset/GRCh38/human_GRCh38_no_alt_analysis_set.fasta.mmi \
    	${bam%.bam}.pbmarkdup.bam \
        ${bam%.bam}.pbmarkdup.GRCh38.sorted.bam
done
```

Calculate coverage of CYP2D6 target region.


```bash
mamba activate cumed_workshop

# https://www.twistbioscience.com/resources/data-files/twist-alliance-long-read-pgx-panel-bed-file
echo -e "chr22\t42120672\t42139664\tCYP2D6" > ~/CUMED_BFX_workshop/03.snakemake_targetenrichment/cyp2d6.bed

cd /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/data/cyp2d6/downsampled

for bam in *.pbmarkdup.GRCh38.sorted.bam; do
    echo $bam ...

    mosdepth \
        -t 4 \
        -b ~/CUMED_BFX_workshop/03.snakemake_targetenrichment/cyp2d6.bed \
        -n \
        -x \
        ${bam%.bam}.mosdepth \
        $bam
done

{
    echo -e "sample\tchr\tstart\tend\tcoverage"
    for f in *.pbmarkdup.GRCh38.sorted.mosdepth.regions.bed.gz; do
        printf "${f%.pbmarkdup.GRCh38.sorted.mosdepth.regions.bed.gz}\t";
        zcat $f
    done
} | column -t
```

|sample|chr|start|end|coverage|
|---|---|---|---|---|
|NA02016|chr22|42120672|42139664|25.49|
|NA07439|chr22|42120672|42139664|24.91|
|NA09301|chr22|42120672|42139664|27.74|
|NA12244|chr22|42120672|42139664|20.29|
|NA16654|chr22|42120672|42139664|23.24|
|NA16688|chr22|42120672|42139664|30.13|
|NA17020|chr22|42120672|42139664|18.90|
|NA17039|chr22|42120672|42139664|25.69|
|NA17073|chr22|42120672|42139664|25.97|
|NA17114|chr22|42120672|42139664|19.79|
|NA17209|chr22|42120672|42139664|21.83|
|NA17211|chr22|42120672|42139664|21.74|
|NA17214|chr22|42120672|42139664|18.65|
|NA17215|chr22|42120672|42139664|24.00|
|NA17217|chr22|42120672|42139664|26.07|
|NA17226|chr22|42120672|42139664|22.67|
|NA17227|chr22|42120672|42139664|23.90|
|NA17232|chr22|42120672|42139664|19.40|
|NA17244|chr22|42120672|42139664|33.29|
|NA17276|chr22|42120672|42139664|21.98|
|NA17282|chr22|42120672|42139664|19.55|
|NA17300|chr22|42120672|42139664|18.63|


### Running pipeline using Snakemake + Slurm

PacBio HiFi target enrichment snakemake workflow is usually run with Slurm and due to slim computing resouces available on workshop servers, it will take time to finish. Therefore, we will not run this session but providing scripts for attendees' reference.

```bash
conda activate snakemake_targetenrichment
git clone https://github.com/PacificBiosciences/HiFiTargetEnrichment.git workflow

cd workflow
git show HEAD
# commit 74fe95f48aee91a50a1e5f650dfd964642d3050f
# Author: jrharting <jharting@pacificbiosciences.com>
# Date:   Thu Mar 23 13:31:06 2023 -0700
cd .. 

mkdir -p reference annotation cluster_logs tmp

# run pipeline
# skip demultiplexing
# update workflow/config.yaml
# cp workflow/config.yaml workflow/config.yaml.CYP2D6_PGx_SQIIe
# update workflow/cluster.yaml
# run just variant calling and phasing for a set of bams following demux/markdup/mapping on SL
# <hifi_reads> can be a directory of bams or a textfile with one bam path per line (fofn)
# update workflow/run_snakemake_SLmapped.sh
#   remove: --singularity-args '--nv ', as not running with GPU
batch_name="CYP2D6_PGx_SQIIe"
sbatch workflow/run_snakemake_SLmapped.sh $batch_name /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/data /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/cyp2d6.bed /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/cyp2d6.bed
```

<!--
```bash
# If pangu (CYP2D6 genotyper) returned too many alleles, pharmcat will fail

# batches/CYP2D6_PGx_SQIIe/logs/pharmcat/run_pharmcat/NA02016.log
# ========
# /bin/bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
# org.pharmgkb.pharmcat.reporter.BadOutsideCallException: Line 1: Too many alleles specified in *2x2/*17/*79
# --
# Saving reporter JSON results to batches/CYP2D6_PGx_SQIIe/NA16654/pharmcat/NA16654.preprocessed.filtered.report.json
# Done.
# batches/CYP2D6_PGx_SQIIe/logs/pharmcat/run_pharmcat/NA16688.log
# ========
# /bin/bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
# org.pharmgkb.pharmcat.reporter.BadOutsideCallException: Line 1: Too many alleles specified in *2/[*36 + *10]/*39x2
# --
# 	at org.pharmgkb.pharmcat.Pipeline.call(Pipeline.java:275)
# 	at org.pharmgkb.pharmcat.PharmCAT.main(PharmCAT.java:166)
# batches/CYP2D6_PGx_SQIIe/logs/pharmcat/run_pharmcat/NA17020.log
# ========
# /bin/bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
# org.pharmgkb.pharmcat.reporter.BadOutsideCallException: Line 1: Too many alleles specified in *1/*5/*10
# --
# batches/CYP2D6_PGx_SQIIe/logs/pharmcat/run_pharmcat/NA17244.log
# ========
# /bin/bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
# chr22:42128211
# 	Discarded genotype at this position because REF in VCF (CG) does not match expected reference (C)
# org.pharmgkb.pharmcat.reporter.BadOutsideCallException: Line 1: Too many alleles specified in *2x3/*4x3/unknown
# --
# batches/CYP2D6_PGx_SQIIe/logs/pharmcat/run_pharmcat/NA17300.log
# ========
# /bin/bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
# chr22:42128211
# 	Discarded genotype at this position because REF in VCF (CG) does not match expected reference (C)
# org.pharmgkb.pharmcat.reporter.BadOutsideCallException: Line 1: Too many alleles specified in *1/*5/*6/unknown

conda activate snakemake_targetenrichment
batch_name="CYP2D6_PGx_SQIIe_withQC"
# update workflow/config.yaml
# cp workflow/config.yaml workflow/config.yaml.CYP2D6_PGx_SQIIe_withQC
# sbatch workflow/run_snakemake_SLmapped.sh $batch_name /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/data /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/cyp2d6.bed /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/cyp2d6.bed
sbatch workflow/run_snakemake_SLmapped.sh $batch_name /home/xliu/CUMED_BFX_workshop/03.snakemake_targetenrichment/data /home/xliu/CUMED_BFX_workshop/03.snakemake_targetenrichment/cyp2d6.bed /home/xliu/CUMED_BFX_workshop/03.snakemake_targetenrichment/cyp2d6.bed
```
-->

### Step-by-step execution of pipeline key steps (hacked run)

Use NA02016 sample as example.

#### pbsv

```bash
mamba env create -n TE.pbsv --file ../snakemake_run/workflow/rules/envs/pbsv.yaml
mamba activate TE.pbsv

pbsv --version
# pbsv 2.8.0 (commit SL-release-10.2.0-8-g20b8f57)

mkdir -p pbsv/svsig
time pbsv discover --hifi --log-level INFO --tandem-repeats /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/workflow/aux/tandem_repeats.bed /home/ubuntu/cumed_bfx_workshop/03.snakemake_targetenrichment/data/pbmarkdup.grch38.sorted.bam/na02016.pbmarkdup.grch38.sorted.bam pbsv/svsig/NA02016.GRCh38_noalt.svsig.gz

# real	0m0.192s
# user	0m0.083s
# sys	0m0.009s

mkdir -p pbsv
time pbsv call --hifi -m 20 --log-level INFO --num-threads 4 /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/reference/human_GRCh38_no_alt_analysis_set.fasta pbsv/svsig/NA02016.GRCh38_noalt.svsig.gz pbsv/NA02016.GRCh38_noalt.pbsv.vcf

# real	2m7.686s
# user	0m12.493s
# sys	0m4.742s
```

#### deepvariant

```bash
# gcr.io/deepvariant-docker/deepvariant:1.4.0
singularity pull deepvariant.sif docker://gcr.io/deepvariant-docker/deepvariant:1.4.0

singularity shell deepvariant.sif

mkdir -p deepvariant/examples
# similar to WGS Variant Pipeline, use parallel to speed up DV make_examples

# pipeline specific paramter:
# --vsc_min_fraction_indels: Indel alleles occurring at least this fraction of
#     all counts in our AlleleCount will be advanced as candidates.
#     (default: '0.06')
# pipeline value: 0.12

# add --regions /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/cyp2d6.bed to confine the number of examples to target region: CYP2D6, this will speed up the process
time seq 0 3 | parallel --jobs 4 /opt/deepvariant/bin/make_examples \
    --add_hp_channel \
    --alt_aligned_pileup=diff_channels \
    --min_mapping_quality=1 \
    --parse_sam_aux_fields \
    --partition_size=25000 \
    --max_reads_per_partition=600 \
    --phase_reads \
    --pileup_image_width 199 \
    --norealign_reads \
    --sort_by_haplotypes \
    --track_ref_reads \
    --vsc_min_fraction_indels 0.12 \
    --mode calling \
    --ref /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/reference/human_GRCh38_no_alt_analysis_set.fasta \
    --reads /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/data/pbmarkdup.GRCh38.sorted.bam/NA02016.pbmarkdup.GRCh38.sorted.bam \
    --examples deepvariant/examples/NA02016.GRCh38_noalt.examples.tfrecord@4.gz \
    --gvcf deepvariant/examples/NA02016.GRCh38_noalt.gvcf.tfrecord@4.gz \
    --regions /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/cyp2d6.bed \
    --task {}

# real	0m5.039s
# user	0m13.199s
# sys	0m1.611s

time /opt/deepvariant/bin/call_variants \
	--outfile deepvariant/NA02016.GRCh38_noalt.call_variants_output.tfrecord.gz \
	--examples deepvariant/examples/NA02016.GRCh38_noalt.examples.tfrecord@4.gz \
	--checkpoint "/opt/models/pacbio/model.ckpt"

# real	0m9.520s
# user	0m8.811s
# sys	0m2.410s

time /opt/deepvariant/bin/postprocess_variants \
	--ref /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/reference/human_GRCh38_no_alt_analysis_set.fasta \
	--infile deepvariant/NA02016.GRCh38_noalt.call_variants_output.tfrecord.gz \
	--outfile deepvariant/NA02016.GRCh38_noalt.deepvariant.vcf.gz \
	--nonvariant_site_tfrecord_path deepvariant/examples/NA02016.GRCh38_noalt.gvcf.tfrecord@4.gz \
	--gvcf_outfile deepvariant/NA02016.GRCh38_noalt.deepvariant.g.vcf.gz

# real	0m2.357s
# user	0m2.072s
# sys	0m0.495s

mamba env create -n TE.bcftools --file ../snakemake_run/workflow/rules/envs/bcftools.yaml
mamba activate TE.bcftools
time bcftools stats --threads 4 --fasta-ref /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/reference/human_GRCh38_no_alt_analysis_set.fasta --apply-filters PASS -s NA02016 deepvariant/NA02016.GRCh38_noalt.deepvariant.vcf.gz > deepvariant/NA02016.GRCh38_noalt.deepvariant.vcf.stats.txt

# real	0m0.050s
# user	0m0.021s
# sys	0m0.016s
```

#### whatshap

```bash
mamba env create -n TE.whatshap --file ../snakemake_run/workflow/rules/envs/whatshap.yaml
mamba activate TE.whatshap

mkdir -p whatshap
time whatshap phase --indels \
    --output whatshap/NA02016.GRCh38_noalt.deepvariant.phased.vcf.gz \
    --reference /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/reference/human_GRCh38_no_alt_analysis_set.fasta \
    deepvariant/NA02016.GRCh38_noalt.deepvariant.vcf.gz \
    /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/data/pbmarkdup.GRCh38.sorted.bam/NA02016.pbmarkdup.GRCh38.sorted.bam

# real	0m0.903s
# user	0m0.733s
# sys	0m0.282s

tabix --version
# tabix (htslib) 1.14
# Copyright (C) 2021 Genome Research Ltd.
tabix -p vcf whatshap/NA02016.GRCh38_noalt.deepvariant.phased.vcf.gz

time whatshap stats \
    --gtf whatshap/NA02016.GRCh38_noalt.deepvariant.phased.gtf \
    --tsv whatshap/NA02016.GRCh38_noalt.deepvariant.phased.tsv \
    --block-list whatshap/NA02016.GRCh38_noalt.deepvariant.phased.blocklist \
    --chr-lengths /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/reference/human_GRCh38_no_alt_analysis_set.chr_lengths.txt \
    whatshap/NA02016.GRCh38_noalt.deepvariant.phased.vcf.gz

# Found 1 sample(s) in input VCF
# Reporting results for sample NA02016
# Phasing statistics for sample NA02016 from file whatshap/NA02016.GRCh38_noalt.deepvariant.phased.vcf.gz
# ---------------- Chromosome chr22 ----------------
#      Variants in VCF:      107
#         Heterozygous:       42 (      38 SNVs)
#               Phased:       34 (      32 SNVs)
#             Unphased:        8 (not considered below)
#           Singletons:        0 (not considered below)
#               Blocks:        2
# 
# Block sizes (no. of variants)
#    Median block size:       17.00 variants
#   Average block size:       17.00 variants
#        Largest block:       27    variants
#       Smallest block:        7    variants
# 
# Block lengths (basepairs)
#       Sum of lengths:    12726    bp
#  Median block length:     6363.00 bp
# Average block length:     6363.00 bp
#        Longest block:     9351    bp
#       Shortest block:     3375    bp
#            Block N50:        0    bp
# 
# real	0m0.478s
# user	0m0.463s
# sys	0m0.231s

time whatshap haplotag --tag-supplementary --ignore-read-groups \
    --output whatshap/NA02016.GRCh38_noalt.deepvariant.haplotagged.bam \
    --reference /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/reference/human_GRCh38_no_alt_analysis_set.fasta \
    whatshap/NA02016.GRCh38_noalt.deepvariant.phased.vcf.gz /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/data/pbmarkdup.GRCh38.sorted.bam/NA02016.pbmarkdup.GRCh38.sorted.bam

# Found 1 sample(s) in input VCF
# Keeping 1 sample(s) for haplo-tagging
# Found 1 sample(s) in BAM file
# Reading alignments and detecting alleles ...
# Found 63 reads covering 34 variants
# 
# == SUMMARY ==
# Total alignments processed:                      1004
# Alignments that could be tagged:                   67
# Alignments spanning multiple phase sets:            0
# haplotag - total processing time: 0.40193605422973633
# 
# real	0m0.884s
# user	0m0.803s
# sys	0m0.287s

mamba env create -n TE.samtools --file ../snakemake_run/workflow/rules/envs/samtools.yaml
mamba activate TE.samtools

samtools --version
# samtools 1.13
# Using htslib 1.18
# Copyright (C) 2021 Genome Research Ltd.
# 
# Samtools compilation details:
#     Features:       build=configure curses=yes 
#     CC:             /opt/conda/conda-bld/samtools_1625833129860/_build_env/bin/x86_64-conda-linux-gnu-cc
#     CPPFLAGS:       -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /home/ubuntu/miniforge3/envs/TE.samtools/include
#     CFLAGS:         -Wall -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/ubuntu/miniforge3/envs/TE.samtools/include -fdebug-prefix-map=/opt/conda/conda-bld/samtools_1625833129860/work=/usr/local/src/conda/samtools-1.13 -fdebug-prefix-map=/home/ubuntu/miniforge3/envs/TE.samtools=/usr/local/src/conda-prefix
#     LDFLAGS:        -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,/home/ubuntu/miniforge3/envs/TE.samtools/lib -Wl,-rpath-link,/home/ubuntu/miniforge3/envs/TE.samtools/lib -L/home/ubuntu/miniforge3/envs/TE.samtools/lib
#     HTSDIR:         
#     LIBS:           
#     CURSES_LIB:     -ltinfow -lncursesw
# 
# HTSlib compilation details:
#     Features:       build=configure libcurl=yes S3=yes GCS=yes libdeflate=yes lzma=yes bzip2=yes plugins=yes plugin-path=/home/ubuntu/miniforge3/envs/TE.samtools/libexec/htslib htscodecs=1.5.1
#     CC:             /opt/conda/conda-bld/htslib_1690291177561/_build_env/bin/x86_64-conda-linux-gnu-cc
#     CPPFLAGS:       -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /home/ubuntu/miniforge3/envs/TE.samtools/include
#     CFLAGS:         -Wall -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/ubuntu/miniforge3/envs/TE.samtools/include -fdebug-prefix-map=/opt/conda/conda-bld/htslib_1690291177561/work=/usr/local/src/conda/htslib-1.18 -fdebug-prefix-map=/home/ubuntu/miniforge3/envs/TE.samtools=/usr/local/src/conda-prefix -fvisibility=hidden
#     LDFLAGS:        -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/home/ubuntu/miniforge3/envs/TE.samtools/lib -Wl,-rpath-link,/home/ubuntu/miniforge3/envs/TE.samtools/lib -L/home/ubuntu/miniforge3/envs/TE.samtools/lib -fvisibility=hidden -rdynamic
# 
# HTSlib URL scheme handlers present:
#     built-in:	 preload, data, file
#     Google Cloud Storage:	 gs+http, gs+https, gs
#     Amazon S3:	 s3+https, s3+http, s3
#     S3 Multipart Upload:	 s3w, s3w+https, s3w+http
#     libcurl:	 imaps, pop3, gophers, http, smb, gopher, sftp, ftps, imap, smtp, smtps, rtsp, scp, ftp, telnet, mqtt, https, smbs, tftp, pop3s, dict
#     crypt4gh-needed:	 crypt4gh
#     mem:	 mem

samtools index -@ 4 whatshap/NA02016.GRCh38_noalt.deepvariant.haplotagged.bam
```

<!--
#### pangu

```bash
# force to use python 3.9.18 instead of >=3.9 as it seems python 3.10+ has some problems with pangu 0.2.2
mamba env create -n TE.pangu --file ../snakemake_run/workflow/rules/envs/pangu.yaml
mamba activate TE.pangu

mkdir -p pangu/NA02016/
time pangu -m capture -p pangu/NA02016/ ../snakemake_run/batches/CYP2D6_PGx_SQIIe/NA02016/whatshap/NA02016.GRCh38_noalt.deepvariant.haplotagged.bam

# awk 'BEGIN {{OFS="\t"}} !($2 ~ /\//) {{$2=$2"/[]"}} 1' pangu/NA02016_pharmcat.tsv > pangu/NA02016_pharmcat_fix.tsv
```
-->

#### pharmcat

```bash
singularity pull pharmcat.sif docker://pgkb/pharmcat:2.3.0

singularity shell pharmcat.sif

mkdir -p pharmcat
/pharmcat/pharmcat_vcf_preprocessor.py \
    --missing-to-ref \
    -vcf whatshap/NA02016.GRCh38_noalt.deepvariant.phased.vcf.gz \
    -refFna /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/reference/human_GRCh38_no_alt_analysis_set.fasta \
    -refVcf /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/workflow/aux/pharmcat_positions.vcf.bgz \
    -bf NA02016 \
    -o pharmcat

# mamba env create -n TE.samtools --file ../snakemake_run/workflow/rules/envs/samtools.yaml
mamba activate TE.samtools

bedtools coverage \
    -sorted \
    -g /home/ubuntu/CUMED_BFX_workshop/03.snakemake_targetenrichment/snakemake_run/reference/human_GRCh38_no_alt_analysis_set.chr_lengths.txt \
    -f 1 \
    -header \
    -mean \
    -a pharmcat/NA02016.preprocessed.vcf.bgz \
    -b ~/CUMED_BFX_workshop/03.snakemake_targetenrichment/data/pbmarkdup.GRCh38.sorted.bam/NA02016.pbmarkdup.GRCh38.sorted.bam \
    | ( sed  -u '/^#CHROM/q' ; awk '$11 >= 10' | cut -f1-10 ) > pharmcat/NA02016.preprocessed.filtered.vcf

singularity shell pharmcat.sif

time /pharmcat/pharmcat \
    -vcf pharmcat/NA02016.preprocessed.filtered.vcf \
    -reporterJson \
    -o pharmcat/
```