#! /home/nrockweiler/software/R-3.6.0/bin/Rscript

#
# IMPORTANT
#
# This code was used to estimate the percent of donor variance explained by the GLM on sample mutation burden using rpt() from the rptR package.  
# This function call takes several hours and be parallelized.  Therefore, this part of the analysis was executed on a high-performance computing cluster.  This script is provided to show how the data was analyzed but is not intended to be run in a different computing environment (i.e., there are hardcoded paths).
options(warn=1)

library("rptR")

source("/home/nrockweiler/mut/mut_corona/src/utils_v8.R")
source("/home/nrockweiler/mut/mut_corona/src/basic_utils.R")
source("/home/nrockweiler/mut/mut_corona/src/plot_utils.R")

local = 4

tissue_palette = create_tissue_palette(local)
today = get_todays_date()


#nboots = c(10) # DEBUG
nboots = c(500, 1000)
ncores = 4

remove_sex_tissues = FALSE
metadata = load_master_sample_metadata(local)
metadata = metadata[metadata$is_hypermutated == 0,]
sex_tissues = get_sex_specific_tissues()

if (remove_sex_tissues) {
  metadata = metadata[! metadata$SMTSD %in% sex_tissues,]
}


metadata$RACE = factor(metadata$RACE, levels=c("White", "Black or African American", "Asian", "American Indian or Alaska Native", "Unknown"))

#coi = c("num_em_var", "num_mbp_w_ge20X", "SMRIN", "SMTSD", "SUBJID", "SMNABTCH", "AGE", "SEX", "RACE", "genotype_data", "avg_simulated_var_power")
#metadata = metadata[,coi] # Try to avoid the extra overhead of QR decomposition of input data for irrelevant columns (https://stackoverflow.com/questions/37090722/lme4lmer-reports-fixed-effect-model-matrix-is-rank-deficient-do-i-need-a-fi)


# Analyze all ancestries combined
rpt_combined_races = vector("list", length=0)
print("combined ancestries")
for (nboot in nboots) {
    print(nboot)

    x = tryCatch({
        y = rpt(log10(num_em_var + 1) ~ num_mbp_w_ge20X + SMRIN + SMTSD + (1|SUBJID) + (1|SMNABTCH) + AGE + SEX + RACE + genotype_data + avg_simulated_var_power + SEX:SMTSD + AGE:SMTSD + RACE:SMTSD,
            data=metadata,
            nboot=nboot,
            grname="SUBJID",
            parallel=TRUE,
            ncores=ncores)

        y
    },
    error=function(cond) {
        message(sprintf("problem with rpt: %s", cond))

        NULL
    })

    rpt_combined_races[[as.character(nboot)]] = x
}

# Analyze each ancestry separately
races = levels(metadata$RACE)
rpt_individual_races = vector("list", length=0)
for (r in races) {
    print(r)
    metadata_subset = metadata[metadata$RACE == r,]
    #metadata_subset$RACE = NULL

    rpt_individual_races[[r]] = vector("list", length=0)
    for (nboot in nboots) {
        print(nboot)

        # Since this is an ancestry-specifc analysis, drop RACE effect and interactions
        if (r == "American Indian or Alaska Native") {
            # There are no females for this ancestry.  Drop SEX effect and interactions
            x = tryCatch({
                    y = rpt(log10(num_em_var + 1) ~ num_mbp_w_ge20X + SMRIN + SMTSD + (1|SUBJID) + (1|SMNABTCH) + AGE + genotype_data + avg_simulated_var_power + AGE:SMTSD,
                        datatype="Gaussian",
                        data=metadata_subset,
                        nboot=nboot,
                        grname="SUBJID", 
                        parallel=TRUE,
                        ncores=ncores)

                    y
                },
            error=function(cond) {
                message(sprintf("problem with rpt: %s", cond))

                NULL
            })
        } else {
            x = tryCatch({
                y = rpt(log10(num_em_var + 1) ~ num_mbp_w_ge20X + SMRIN + SMTSD + (1|SUBJID) + (1|SMNABTCH) + AGE + SEX + genotype_data + avg_simulated_var_power + SEX:SMTSD + AGE:SMTSD,
                    datatype="Gaussian",
                    data=metadata_subset,
                    nboot=nboot,
                    grname="SUBJID", 
                    parallel=TRUE,
                    ncores=ncores)

                y
            },
            error=function(cond) {
                message(sprintf("problem with rpt: %s", cond))

                NULL
            })

            rpt_individual_races[[r]][[as.character(nboot)]] = x
        }
    }
}



rdata_fn = sprintf("repeat_estimates.%s.Rdata", today)
save(list=ls(all.names=TRUE), file=rdata_fn)
