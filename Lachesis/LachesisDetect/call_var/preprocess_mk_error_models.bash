#!/usr/bin/bash

set -ue
set -o pipefail
shopt -s extglob

# USAGE: preprocess_mk_error_models.bash <set id> <project> [<num_samples_per_subset>]
#
# DESCRIPTION: preprocess and kick off make error model jobs.  The following steps are performed:
# - given a set ID, divide the list of samples in that set into smaller groups 
# - For each sample subset, generate a file of readcount files (needed for generating the error model)
# - Check if all the readcount files exist
# - If all files exist, submit the jobs to generate the error model as a slurm array

num_samples_per_subset=100

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]
then
    >&2 "ERROR: incorrect number of options"
    exit
elif [ "$#" -eq 3 ]
then
    num_samples_per_subset=$3
fi

set_id=$1
project=$2


if [ $project == 'gtex' ]; then
    READCOUNT_DIR="/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/results/stranded_readcount_files.transcriptome"
    SAMPLE_SET_DIR="/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/data/reference/sets/by_name_w_stragglers" # Use the stragglers since the samples were reorganized using this convention.
    ANALYSIS_DIR="/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/results/20190221_mk_error_mdl_v2"

    sample_set_fn="${SAMPLE_SET_DIR}/list_of_gtexv8_sample_ids.by_name_w_stragglers.${set_id}.20180919.txt"
    readcount_subdir="${READCOUNT_DIR}/${set_id}"
elif [ $project == 'twins_uk' ]; then
    readcount_subdir="/scratch/dclab/twins/alignments/${set_id}/stranded_readcounts"
    SAMPLE_SET_DIR="/scratch/dclab/nrockweiler/twins_uk/data/reference"
    ANALYSIS_DIR="/scratch/dclab/nrockweiler/twins_uk/results/20191120_mk_error_mdl"

    sample_set_fn="${SAMPLE_SET_DIR}/list_of_samples.${set_id}.ega.*.txt"
else
    >&2 echo "ERROR: unknown project $project"
    exit
fi

READCOUNT_SUFFIX="_stranded_readcounts.filtered_alignments.poi.bed.gz"


analysis_subdir="${ANALYSIS_DIR}/${set_id}"


mkdir -p $analysis_subdir
cd $analysis_subdir

# Create the list of readcounts file
readcount_subset_prefix="readcount_subset.${set_id}."
awk -v dir="$readcount_subdir" -v suffix="$READCOUNT_SUFFIX" 'BEGIN{OFS="\t"; FS="\t"} {print $1, dir "/" $1 suffix}' $sample_set_fn | split -l $num_samples_per_subset -a 2 --numeric-suffixes=1 --additional-suffix=".txt" - $readcount_subset_prefix

# Remove the 0 padding
rename -f s/\.$set_id.0/.$set_id./ ${readcount_subset_prefix}*.txt

# Calculate the number of subsets
num_readcount_subsets=`ls -1 ${readcount_subset_prefix}+([0-9]).txt | wc -l`

# Sanitize for the near 100% of samples that are available
python3 /home/nrockweiler/google_cloud_platform/run_bcall/src/firecloud/check_have_all_files.py -f $sample_set_fn -sd $readcount_subdir

exit_code=$?
if [ $exit_code -eq 0 ]; then
    # Submit mk error model
    sbatch --array=1-$num_readcount_subsets /home/nrockweiler/mut/mut_corona/src/call_var/bcall_mk_error_model.slurm.bash "readcount_subset.${set_id}." # command line sbatch arguments trump values defined in job script (https://help.rc.ufl.edu/doc/Using_Variables_in_SLURM_Jobs)
else
    # Note: since set -e, won't ever get here.  However, check_have_all_files.py will error out.
    >&2 echo "ERROR: missing samples.  Mk error model job NOT submitted"
fi

