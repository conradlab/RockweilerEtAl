#!/usr/bin/bash

set -ue
set -o pipefail

# USAGE: preprocess_merge_error_models.bash <set id>
#
# DESCRIPTION: preprocess and kick off merge error model jobs.  The following steps are performed:
# - create list of prior dumps file
# - submit job to merge the error models

set_id=$1
main_analysis_dir=$(readlink -m $2)

# For GTEx, main_analysis_dir = "/scratch/dclab/nrockweiler/google_cloud_platform/run_bcall/results/20190221_mk_error_mdl_v2"

DATE=`date +%Y%m%d`


analysis_subdir="${main_analysis_dir}/${set_id}"
mkdir -p $analysis_subdir
cd $analysis_subdir

# Create list of prior dumps
list_of_prior_dumps_fn="list_of_prior_dumps.${set_id}.${DATE}.txt"
find . -name  "readcount_subset*dump" | sort -k1,1V | awk 'BEGIN{FS="/"; OFS="\t"} {print $2, $0}' > $list_of_prior_dumps_fn

# Submit merge error model
dump_prefix="merged_dumps.${set_id}.${DATE}"
sbatch /home/nrockweiler/mut/mut_corona/src/call_var/bcall_merge_error_models.slurm.bash $list_of_prior_dumps_fn $dump_prefix
