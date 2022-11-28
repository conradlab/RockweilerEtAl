suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))

create_tissue_palette <- function(local) {
    # Creates a ggplot color and fill scale for consistent coloring of tissues
    
    tissue_colors = load_tissue_colors(local)
    mycolors = tissue_colors$color_code
    names(mycolors) = tissue_colors$SMTSD
    fill_scale <- scale_fill_manual(name = "", values=mycolors)
    col_scale <- scale_color_manual(name = "", values=mycolors)
    
    return(list(fill_scale=fill_scale, col_scale=col_scale))
}

load_tissue_colors <- function(local) {
    metadata_dir = get_metadata_dir(local)
    tissue_colors_fn = file.path(metadata_dir, "tissue2color.txt")
    
    tissue_colors_df = load_df(tissue_colors_fn)
    return(tissue_colors_df)
}

get_metadata_dir <- function(local) {
    
    metadata_dir = "../data"
    
    return(metadata_dir)
}

load_df <- function(fn, header=TRUE) {
    colors_df = data.frame(read.table(fn, header=header, stringsAsFactors=FALSE, comment.char="", sep="\t", quote=""))
    return(colors_df)
}

get_todays_date <- function() {
    today = sprintf("%s", format(Sys.time(), "%Y%m%d"))
    return(today)
}

load_master_sample_metadata <- function(local) {
    metadata_dir = get_metadata_dir(local)
    
    master_sample_metadata_fn = file.path(metadata_dir, "merged_sample_metadata.filtered.20220422.tsv")
    
    master_sample_metadata_df = load_df(master_sample_metadata_fn)
    return(master_sample_metadata_df)
}

sanitize_var_for_fn <- function(var) {
    var = tolower(gsub(">|-| |\\(|\\)|\\|", "_", var)) # Replace ">", "-", " ", "(", ")" and "|"
    var = gsub("_+","_", var)
    var = gsub("_$","", var)
    return(var)
}

write_new_file_message <- function(fn) {
    writeLines(sprintf("INFO: created %s", fn), con=stderr())
}

my_theme <- function () {
    # Inspired by http://joeystanley.com/blog/custom-themes-in-ggplot2
    theme_bw(base_size=12) %+replace% 
        theme(
            plot.title = element_text(hjust = 0.5) # Center title
        )
}

timestamp_figure <- function(line=-1,side=1,outer=TRUE,adj=0, cex=0.5) {
    mtext(format(Sys.time(), "%Y%m%d %H:%M"), cex=cex, line=line, side=side, adj=adj, outer=outer)
}

timestamp_ggfigure <- function(p) {
    footnote = format(Sys.time(), "%Y%m%d %H:%M") # From http://statmodeling.com/best-way-to-add-a-footnote-to-a-plot-created-with-ggplot2.html
    #grid.newpage()
    g <- arrangeGrob(p, bottom = textGrob(footnote, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontface = "italic", fontsize = 8)))
    
    return(g)
}

create_vartype_palette <- function(local) {
    # Creates a ggplot color and fill scale for consistent coloring of variant types
    vartype_colors_df = load_vartype_colors(local)
    
    mycolors = vartype_colors_df$color_code
    names(mycolors) = vartype_colors_df$X.var_type
    fill_scale <- scale_fill_manual(name = "", values=mycolors)
    col_scale <- scale_color_manual(name = "", values=mycolors)
    
    return(list(fill_scale=fill_scale, col_scale=col_scale))
}

load_vartype_colors <- function(local) {
    metadata_dir = get_metadata_dir(local)
    vartype_colors_fn = file.path(metadata_dir, "vartype2color.txt")
    
    vartype_colors_df = load_df(vartype_colors_fn)
    return(vartype_colors_df)
}

generate_compressed_var <- function () {
    ### Create a lookup table for variant to compressed variant and a lookup table for compressed variant to variant ###
    
    NUCS = c('C', 'T', 'G', 'A') # Go in this out-of-alphabetical order so the reference base is is always a pyrimidine (and matches literature)
    compressed2uncompressed = list()
    uncompressed2compressed = list()
    for (i in NUCS) {
        for (j in NUCS) {
            # Since these are var, A>A, C>C, etc. are NOT valid variants
            if (i != j) {
                variant = paste(i, j, sep=">")
                
                variant_compliment = seq2seq_compliment(variant)
                if (variant_compliment %in% compressed2uncompressed) {
                    compressed2uncompressed[[variant_compliment]] = c(compressed2uncompressed[[variant_compliment]], variant)
                    uncompressed2compressed[[variant]] = variant_compliment
                } else {
                    compressed2uncompressed[[variant]] = variant
                    uncompressed2compressed[[variant]] = variant
                }
            }
        }
    }
    
    # Sort in alphabetical order
    compressed2uncompressed = compressed2uncompressed[order(names(compressed2uncompressed))]
    uncompressed2compressed = uncompressed2compressed[order(names(uncompressed2compressed))]
    
    return(list(compressed2uncompressed=compressed2uncompressed, uncompressed2compressed=uncompressed2compressed))
}

seq2seq_compliment <- function(seq) {
    
    seq_compliment = chartr("ACGT", "TGCA", seq)
    
    return(seq_compliment)
}