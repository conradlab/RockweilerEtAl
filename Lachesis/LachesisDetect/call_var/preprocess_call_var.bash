#!/usr/bin/bash

set -ue
set -o pipefail
shopt -s extglob

# USAGE: preprocess_call_var.bash <set id> <project>
#
# DESCRIPTION: preprocess and kick off call var jobs.  The following steps are performed:
# - Given a set ID, make the analysis subdir
# - Divide the list of sample in that set into smaller groups (N = 50 samples per subset)
# - For each sample subset, gerneate a file of readcount files
# - Submit the call var job as a slurm array (each readcount subset is sumitted as a job in the array) 

set_id=$1
project=$2

if [ $project == 'gtex' ]; then
    error_dump_fn="/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/results/20190221_mk_error_mdl_v2/merged_dumps.all_sets.20190306.dump"
    ANALYSIS_DIR="/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/results/20190301_call_var_v2"
    SAMPLE_SET_DIR="/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/data/reference/sets/by_name_w_stragglers" # Use the stragglers since the samples were reorganized using this convention.
    READCOUNT_DIR="/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/results/stranded_readcount_files.transcriptome"
    readcount_subdir="${READCOUNT_DIR}/${set_id}"
    sample_set_fn="${SAMPLE_SET_DIR}/list_of_gtexv8_sample_ids.by_name_w_stragglers.${set_id}.20180919.txt"
elif [ $project == 'twins_uk' ]; then

    error_mdl_type=$3
    if [ $error_mdl_type == "tissue_specific" ]; then
        error_dump_fn="/scratch/dclab/nrockweiler/twins_uk/results/20191120_mk_error_mdl/${set_id}/merged_dumps.${set_id}.*.dump"
    elif [ $error_mdl_type == "combined_tissues" ]; then
        error_dump_fn="/scratch/dclab/nrockweiler/twins_uk/results/20191120_mk_error_mdl/combined_tissues/merged_dumps.combined_tissues.20191127.dump"
    else
        >&2 echo "ERROR: unknown error model type '$error_mdl_type'"
        exit
    fi

    ANALYSIS_DIR="/scratch/dclab/nrockweiler/twins_uk/results/20191125_call_var/${error_mdl_type}"
    SAMPLE_SET_DIR="/scratch/dclab/nrockweiler/twins_uk/data/reference"
    sample_set_fn="${SAMPLE_SET_DIR}/list_of_samples.${set_id}.ega.*.txt"
    readcount_subdir="/scratch/dclab/twins/alignments/${set_id}/stranded_readcounts"
fi



READCOUNT_SUFFIX="_stranded_readcounts.filtered_alignments.poi.bed.gz"
N=50 # # of samples to process in each subset


analysis_subdir="${ANALYSIS_DIR}/${set_id}"
readcount_subset_prefix="readcount_subset.${set_id}."

mkdir -p $analysis_subdir
cd $analysis_subdir

# Create the list of readcounts file
awk -v dir="$readcount_subdir" -v suffix="$READCOUNT_SUFFIX" 'BEGIN{OFS="\t"; FS="\t"} {print $1, dir "/" $1 suffix}' $sample_set_fn | split -l $N -a 2 --numeric-suffixes=1 --additional-suffix=".txt" - $readcount_subset_prefix

# Remove the 0 padding
rename -f s/\.$set_id.0/.$set_id./ $readcount_subset_prefix*

num_subsets=`ls -1 ${readcount_subset_prefix}+([0-9]).txt | wc -l`

# Submit call var job
sbatch --array=1-$num_subsets /home/nrockweiler/mut/mut_corona/src/call_var/bcall_call_var.slurm.bash $error_dump_fn $readcount_subset_prefix

