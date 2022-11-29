#!/bin/bash

#SBATCH --array=1-17

#SBATCH -o bcall_mk_error_model.o%A.%a
#SBATCH -e bcall_mk_error_model.e%A.%a
#SBATCH --mem=16G

# DESCRIPTION
# - Submit a slurm array of jobs to generate the bcall error models.  Output files are created in the CWD.
#
# USAGE
# $ sbatch bcall_mk_error_model.slurm.bash <subset_fn_prefix>
#
# <subset_fn_prefix> = Prefix of subset filenames.  IMPORTANT: include trailing "." (if applicable).  Script will look for files matching "<subset_prefix>.<array_index>.txt."
# NOTE
# - If resubmit in the same CWD, the original output files will be overwritten.

ID=${SLURM_ARRAY_TASK_ID}

BCALL="/home/nrockweiler/software/bcall/build/bcall_stranded"

# Input files
subset_prefix=$1 # E.g., readcount_subset.
subset_fn="${subset_prefix}${ID}.txt"
subset_bn=${subset_fn##*/}
subset_bn_no_ext=${subset_bn%.txt}

# Output files
dump_fn="${subset_bn_no_ext}.dump"
total_counts_fn="${subset_bn_no_ext}.total_counts.txt"

$BCALL prior-dump $subset_fn $dump_fn > $total_counts_fn
