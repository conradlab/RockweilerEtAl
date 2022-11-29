#!/usr/bin/bash

set -euo pipefail

# Usage
# bash run_rnaseq_qc.bash <bam> <output_dir>
#
# DESCRIPTION
# Run rna-seqc (v2.3.4) on a given bam.  Sample name is inferred from the bam prefix.

bam=`realpath $1`
output_dir=$2

bam_bn=${bam##*/} # Ex: EB_TWPID4924_B.bam
sample_name=${bam_bn%.bam} # Ex: EB_TWPID4924_B

RNASEQC="/home/nrockweiler/software/rnaseqc.v2.3.4/rnaseqc.v2.3.4.linux"
COLLAPSED_GTF="/home/nrockweiler/mut/mut_corona/data/20161103_reference/gencodev26/gencode.v26.GRCh38.genes.gtf"



finished_rnaseqc_fn="${output_dir}/.finished_rnaseqc.${sample_name}"


if [ -s $finished_rnaseqc_fn ]; then

    >&2 echo "rnaseqc already run.  not rerunning rnaseqc."

else

    $RNASEQC $COLLAPSED_GTF $bam $output_dir -s $sample_name
    
    touch $finished_rnaseqc_fn

fi
