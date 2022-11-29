#!/usr/bin/bash

set -euo pipefail
shopt -s nullglob

# DESCRIPTION: perform the following steps:
# - Follow GTEx pipeline:
#   - mark duplicates
# - convert bam to filtered, stranded readcount file at positions of interest
#
# Based on bam2realigned_filtered_readcount.bash, except that pipeline *starts* after alignment.  This is because the HTCF cluster doesn't have many >30G nodes.  Therefore, perform the alignment stages with a high mem node, then exit job so can make room for other samples to use that node.  Finish the rest of the pipeline on a (lower mem) node.


input_bam=`realpath $1`

if [ $# -eq 2 ]; then
    output_dir=$2
    mkdir -p $output_dir
    output_dir=`realpath $output_dir`
    cd $output_dir
else
    output_dir=`pwd`
fi

fastq_dir=`realpath "fastqs"`
mkdir -p $fastq_dir

stranded_readcounts_dir=`realpath "stranded_readcounts"`
mkdir -p $stranded_readcounts_dir

raw_bam_dir=`realpath "star_bams"`
mkdir -p $raw_bam_dir

mrk_dup_dir=`realpath "star_bams.mrk_dup"`
mkdir -p $mrk_dup_dir

filtered_bam_dir=`realpath "star_bams.mrk_dup.filtered"`
mkdir -p $filtered_bam_dir

ref_fasta="/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/data/reference/Homo_sapiens_assembly38.fasta"
positions_of_interest_bed="/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/data/reference/gencode.v26.GRCh38.genes.exons_only.sort.merge.20180420.bed"


input_bam_bn=${input_bam##*/} # Ex: EB_TWPID4924_B.bam
sample_name=${input_bam_bn%.bam} # Ex: EB_TWPID4924_B

# Output files
r1_fastq="${fastq_dir}/${sample_name}_1.fastq.gz"
r2_fastq="${fastq_dir}/${sample_name}_2.fastq.gz"
realigned_bam="${raw_bam_dir}/${sample_name}.Aligned.sortedByCoord.out.bam"
markdup_bam="${mrk_dup_dir}/${sample_name}.Aligned.sortedByCoord.out.md.bam"
markdup_bai="${markdup_bam}.bai"
filtered_bam="${filtered_bam_dir}/${sample_name}.Aligned.sortedByCoord.out.md.qc_filtered.bam"
filtered_bai="${filtered_bam}.bai"
output_readcount="${stranded_readcounts_dir}/${sample_name}_stranded_readcounts.filtered_alignments.poi.bed.gz"
finished_fn="${stranded_readcounts_dir}/.finished_readcount.${sample_name}"


# Programs
RUN_SAM2FASTQ="/home/nrockweiler/google_cloud_platform/run_bcall/src/gtex-pipeline/rnaseq/src/run_SamToFastq.py"
RUN_MARKDUP="/home/nrockweiler/google_cloud_platform/run_bcall/src/gtex-pipeline/rnaseq/src/run_MarkDuplicates.py"
SAMTOOLS="/opt/htcf/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/samtools-1.9-jjq5nuaykk2hcqw2v7qimvyvdp7fpvvr/bin/samtools"
SAMCLIP="/home/nrockweiler/google_cloud_platform/run_bcall/src/bam2filtered_readcounts/samclip.pl"
RUN_STAR="/home/nrockweiler/google_cloud_platform/run_bcall/src/gtex-pipeline/rnaseq/src/run_STAR.py"
RUN_MARKDUP="/home/nrockweiler/google_cloud_platform/run_bcall/src/gtex-pipeline/rnaseq/src/run_MarkDuplicates.py"
MPILEUP2STRANDED_READCOUNTS="/home/nrockweiler/google_cloud_platform/run_bcall/src/bam2filtered_readcounts/mpileup2stranded_readcounts.py"
BGZIP="/opt/htcf/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/htslib-1.9-5ezijoygivagniwt7c2lya6boizwsmmc/bin/bgzip"

view_opt="-h -q 20 -F 1796"
mpileup_opt="-O -q 20 -Q 20 -d 1000000000" # NOTE: by default, mpileup will remove marked duplicates.  Additionally, the -Q 20 filter will effectively filter for uniquely mapped reads


# Estimate the read length from the bam.  (The samples have differetn read lengths.)
num_lines=50000
# I get a pipe error 141.  Not sure why.  Temporarily turn off pipefail.
set +o pipefail
read_length=$($SAMTOOLS view $input_bam | head -n $num_lines | awk '{print length($10)}' | sort | uniq -c | sort -k2,2nr | head -n 1 | awk '{print $2}')
set -o pipefail


#
# Mark duplicates
#

if [ -s $filtered_bai ]; then
    >&2 echo "filtered bam already existed.  not going to filter  again"

    if [ -f $finished_fn ]; then
        >&2 echo "stranded readcounts finished file already existed.  not going to rerun filter"
    else
        # Filter bam
        # At positions of interest, filter alignments for good quality alignments.  Output the information in bam and stranded readcount formats
        $SAMTOOLS mpileup $mpileup_opt -f $ref_fasta -l $positions_of_interest_bed $filtered_bam | \
            python3.5 $MPILEUP2STRANDED_READCOUNTS -s $sample_name --output-format bed --read-length $read_length | \
            $BGZIP > $output_readcount
    fi
elif [ -s $markdup_bai ]; then

    # Filter bam
    # At positions of interest, filter alignments for good quality alignments.  Output the information in bam and stranded readcount formats
    
    $SAMTOOLS view $view_opt -L $positions_of_interest_bed $markdup_bam | \
        perl $SAMCLIP --ref $ref_fasta - 2> "${sample_name}.samclip_stats.poi.txt" | \
        $SAMTOOLS view -b | \
        tee $filtered_bam | \
        $SAMTOOLS mpileup $mpileup_opt -f $ref_fasta -l $positions_of_interest_bed - | \
        python3.5 $MPILEUP2STRANDED_READCOUNTS -s $sample_name --output-format bed --read-length $read_length | \
        $BGZIP > $output_readcount
    
        $SAMTOOLS index $filtered_bam

else

    python3 $RUN_MARKDUP $realigned_bam $sample_name -o $mrk_dup_dir
    # From XXX.bam, creates  XXX.md.bam
    
    # Index markdup bam
    $SAMTOOLS index $markdup_bam

    # Filter bam
    # At positions of interest, filter alignments for good quality alignments.  Output the information in bam and stranded readcount formats
    
    $SAMTOOLS view $view_opt -L $positions_of_interest_bed $markdup_bam | \
        perl $SAMCLIP --ref $ref_fasta - 2> "${sample_name}.samclip_stats.poi.txt" | \
        $SAMTOOLS view -b | \
        tee $filtered_bam | \
        $SAMTOOLS mpileup $mpileup_opt -f $ref_fasta -l $positions_of_interest_bed - | \
        python3.5 $MPILEUP2STRANDED_READCOUNTS -s $sample_name --output-format bed --read-length $read_length | \
        $BGZIP > $output_readcount
    
        $SAMTOOLS index $filtered_bam
fi

touch $finished_fn

#
# Clean up
#

# Toss intermediate files
if [ -f "$filtered_bam" ]; then
    rm -f $realigned_bam "${realigned_bam}.bai"
    rm -f $markdup_bam "${markdup_bam}.bai"
fi

