#!/bin/bash

#SBATCH -o bcall_merge_error_models.o%j
#SBATCH -e bcall_merge_error_models.e%j
#SBATCH --mem=32G

# DESCRIPTION
# - Submit a slurm job to merge the error models generated from bcall_mk_error_model.slurm.bash.
#
# USAGE
# $ sbatch bcall_merge_error_models.slurm.bash <list_of_prior_dumps_fn> <output_prefix>
#
# <list_of_prior_dumps_fn> = Path to file containing list of dump filenames.  Format: Col1: name of the subset, col2: path to the dump file (generated from bcall_mk_error_model.slurm.bash)
# <output_prefix> = Prefix to add to the output file (a .dump file)

if [ "$#" -ne 2 ]; then
    >&2 echo "ERROR: incorrect number of arguments."
    exit 1
fi

BCALL="/home/nrockweiler/software/bcall/build/bcall_stranded"

list_of_prior_dumps_fn=$1
output_prefix=$2

# Output files
merged_dump_fn="${output_prefix}.dump"
total_counts_fn="${output_prefix}.total_counts.txt"

# Run bcall
$BCALL prior-merge $list_of_prior_dumps_fn $merged_dump_fn

# Print plain-text
$BCALL prior-print $merged_dump_fn > $total_counts_fn 
