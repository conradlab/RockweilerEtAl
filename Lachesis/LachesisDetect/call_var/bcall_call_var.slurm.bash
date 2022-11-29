#!/bin/bash

#SBATCH --array=1-34

#SBATCH -o bcall_call_var.o%A.%a
#SBATCH -e bcall_call_var.e%A.%a
#SBATCH --mem=16G

# DESCRIPTION
# - Submit a slurm array of jobs to call variants using bcall.  Each array job runs bcall on a subset of the samples.  Output files are created in the CWD.
#
# USAGE
# $ sbatch bcall_call_var.slurm.bash <merged_dump_fn> <subset_prefix>
#
# <merged_dump_fn> = Path to merged dump file.  File is generated from bcall_merge_error_models.slurm.bash.
# <subset_fn_prefix> = Prefix of subset filenames.  IMPORTANT: include trailing "." (if applicable).  Script will look for files matching "<subset_prefix>.<array_index>.txt."
#
# NOTE
# - If resubmit in the same CWD, the original output files will be overwritten.

if [ "$#" -ne 2 ]; then
    >&2 echo "ERROR: incorrect number of arguments."
    exit 1
fi

merged_dump_fn=$1
subset_prefix=$2 # E.g., readcount_subset.

ID=${SLURM_ARRAY_TASK_ID}

BCALL="/home/nrockweiler/software/bcall/build/bcall_stranded"

# Input files
subset_fn="${subset_prefix}${ID}.txt"
subset_bn=${subset_fn##*/}
subset_bn_no_ext=${subset_bn%.txt}

# Output files
var_calls_fn="${subset_bn_no_ext}.var_calls.txt"

$BCALL call-using-merged $subset_fn $merged_dump_fn > $var_calls_fn
