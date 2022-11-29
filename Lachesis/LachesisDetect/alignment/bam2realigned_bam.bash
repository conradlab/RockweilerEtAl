#!/usr/bin/bash

set -euo pipefail
shopt -s nullglob

# DESCRIPTION: perform the following steps:
# - Follow GTEx pipeline:
#   - Convert to fastq (star does not play well with our bams)
#   - Align
# - Submit downstream script (realigned_bam2filtered_readcount.bash)
# Based on bam2realigned_filtered_readcount.bash, except that pipeline *ends( after alignment.  This is because the HTCF cluster doesn't have many >30G nodes.  Therefore, perform the alignment stages with a high mem node, then exit job so can make room for other samples to use that node.  Finish the rest of the pipeline on a (lower mem) node.

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
r_unpaired_fastq="${fastq_dir}/${sample_name}_unpaired.fastq.gz"
realigned_unsorted_bam="${raw_bam_dir}/${sample_name}.Aligned.toTranscriptome.out.bam"
realigned_bam="${raw_bam_dir}/${sample_name}.Aligned.sortedByCoord.out.bam"
realigned_bai="${realigned_bam}.bai"
markdup_bam="${mrk_dup_dir}/${sample_name}.Aligned.sortedByCoord.out.md.bam"
filtered_bam="${filtered_bam_dir}/${sample_name}.Aligned.sortedByCoord.out.md.qc_filtered.bam"
output_readcount="${stranded_readcounts_dir}/${sample_name}_stranded_readcounts.filtered_alignments.poi.bed.gz"
finished_alignment_fn="${raw_bam_dir}/.finished_alignment.${sample_name}"
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

# If the final file already exists, don't bother re-running alignment
if [ -s $finished_alignment_fn ]; then
    >&2 echo "realigned bam already existed.  not running run_STAR.py"
else
    #
    # Convert bam to fastq
    #
    
    # Estimate the read length from the bam.  (The samples have differetn read lengths.)
    num_lines=50000
    # I get a pipe error 141.  Not sure why.  Temporarily turn off pipefail.
    set +o pipefail
    read_length=$($SAMTOOLS view $input_bam | head -n $num_lines | awk '{print length($10)}' | sort | uniq -c | sort -k2,2nr | head -n 1 | awk '{print $2}')
    set -o pipefail
    
    # Given a read length of N, use an index of N - 1
    overhang="$(($read_length - 1))"
    star_index="/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/data/reference/hg38_star_index/overhang_${overhang}"
    
    
    # Convert bam to fastq
    python3 $RUN_SAM2FASTQ $input_bam -p $sample_name -o $fastq_dir
    # From sample XXX, creates XXX_1.fastq.gz, XXX_2.fastq.gz
    
    
    #
    # Align
    #
    
    python3 $RUN_STAR $star_index $r1_fastq $r2_fastq --prefix $sample_name -o $raw_bam_dir
    # From sample XXX, creates XXX.Aligned.sortedByCoord.out.bam (and index)

    touch $finished_alignment_fn
fi

if [ -s $finished_fn ]; then    
    >&2 echo "stranded readcount already existed.  not running realigned_bam2filtered_readcount"
else
    # Submit downstream process
    SUBMIT_2READCOUNT="/home/nrockweiler/twins_uk/src/alignment/realigned_bam2filtered_readcount.slurm.bash"
    echo "$SUBMIT_2READCOUNT $input_bam $output_dir"
    sbatch $SUBMIT_2READCOUNT $input_bam $output_dir

fi

# Toss intermediate files
if [ -f "$realigned_bai" ]; then
    rm -f $r1_fastq $r2_fastq $r_unpaired_fastq
    rm -f $realigned_unsorted_bam "${realigned_unsorted_bam}.bai"
fi

