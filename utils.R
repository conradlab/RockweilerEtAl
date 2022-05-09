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

timestamp_ggfigure <- function(p) {
    footnote = format(Sys.time(), "%Y%m%d %H:%M") # From http://statmodeling.com/best-way-to-add-a-footnote-to-a-plot-created-with-ggplot2.html
    #grid.newpage()
    g <- arrangeGrob(p, bottom = textGrob(footnote, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontface = "italic", fontsize = 8)))
    
    return(g)
}