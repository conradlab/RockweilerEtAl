suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(RColorBrewer))

if (packageVersion("igraph") < "1.1.1") {
    stop("Must upgrade igraph to >= 1.1.1")
}

#
# Functions
# 
vertex2index <- function (g, v) {
    index = which(V(g)$name == v)
    
    return(index)
}

index2vertex <- function(g, index) {
    v = V(g)[index]
    
    return(v)
}

get_parent_vertex <- function(g, child_vertex) {
    my_edge = get.edges(g, E(g)[to(child_vertex)])
    if (length(my_edge)) {
        parent_index = my_edge[1]
        parent = index2vertex(g, parent_index)
    } else {
        parent = FALSE # There is no parent
    }
    return(parent)
}

update_parent <- function(g, v_index, ref, root_index) {
    # Get child vertex
    v = index2vertex(g, v_index)
    # Get parent vertex
    parent_v = get_parent_vertex(g, v)
    
    # Get the parent's current variant value
    a = get.vertex.attribute(g, "VARIANT", index=parent_v) # If this is the first time the parent's attribute has been queried, it will be NA.
    
    # Get the uncertainty value for parent and child
    child_uncertainty = get.vertex.attribute(g, "VARIANT_CONTAINS_NA", index=v)
    parent_uncertainty = get.vertex.attribute(g, "VARIANT_CONTAINS_NA", index=parent_v)
    
    # Get the child's variant value
    v_a = get.vertex.attribute(g, "VARIANT", index=v)
    
    # Update the parent's variant value with the child's value
    # If the parent's allele has not been assigned yet
    
    # Also update the uncertainty
    # If parent is not ref, move uncertainty upwards
    if (is.na(a)) {
        new_value = v_a
        # Else, 1) is root or 2) one of the other parent's children has set the value.  If 1: keep as is.  If 2: If there is consensus, keep as is, if unknown, keep as is, else set to confict (i.e., ref)
        if (child_uncertainty) {
            parent_uncertainty = 1
        }
    } else if (parent_v == root_index) {
        new_value = a
        parent_uncertainty = 0
    } else {
        if (is.na(v_a)) { # Need to check if is.na first, otherwise operators like == on na will return na.
            new_value = a
            if (child_uncertainty) {
                parent_uncertainty = 1
            }
        } else if (a == v_a) {
            new_value = a
            if (child_uncertainty) {
                parent_uncertainty = 1
            }
        } else {
            new_value = ref
            parent_uncertainty = 0
        }
    }
    
    g = set.vertex.attribute(g, "VARIANT", index=parent_v, value=new_value)
    g = set.vertex.attribute(g, "VARIANT_CONTAINS_NA", index=parent_v, value=parent_uncertainty)
    
    return(g)
}


update_parent_unweighted <- function(g, v_index, ref, root_index) {
    # Get child vertex
    v = index2vertex(g, v_index)
    # Get parent vertex
    parent_v = get_parent_vertex(g, v)
    
    # Get the parent's current variant value
    a = get.vertex.attribute(g, "VARIANT", index=parent_v) # If this is the first time the parent's attribute has been queried, it will be NA.
    
    # Get the child's variant value
    v_a = get.vertex.attribute(g, "VARIANT", index=v)
    
    # Update the parent's variant value with the child's value
    # If the parent's allele has not been assigned yet
    if (is.na(a)) {
        new_value = v_a
        
        # Else, 1) is root or 2) one of the other parent's children has set the value.  If 1: keep as is.  If 2: If there is consensus, keep as is, if unknown, keep as is, else set to confict (i.e., ref)
    } else if (parent_v == root_index) {
        new_value = a
    } else {
        if (is.na(v_a)) { # Need to check if is.na first, otherwise operators like == on na will return na.
            new_value = a
        } else if (a == v_a) {
            new_value = a
        } else {
            new_value = ref
        }
    }
    
    g = set.vertex.attribute(g, "VARIANT", index=parent_v, value=new_value)
    
    return(g)
}

add_variant_to_tree_OLD <- function(g, root_index, ref, vartype, vartype_colors_palette, alt="T") {
    b = bfs(g, root_index, dist=TRUE) # TODO: could save this since the BFS is the same for all variants
    vertex_order = sort(b$dist, decreasing=TRUE)
    
    vertex_order_names = names(vertex_order)
    
    for (i in 2:length(vertex_order)) {
        level_previous = vertex_order[i-1]
        v_previous_index = vertex_order_names[i-1]
        
        level_current = vertex_order[i]
        v_current_index = vertex_order_names[i]
        
        g = update_parent(g, v_previous_index, ref, root_index)
    }
    
    # Now that the tree has all the internal variants added, construct the edge types
    g = add_edge_annotation(g, vartype, ref, vartype_colors_palette, b, alt=alt)
    return(g)
}

add_variant_to_tree <- function(g, ancestor_paths, ref, vartype, vartype_colors_palette, alt="T") {
    
    var_vertexes = V(g)[which(V(g)$VARIANT == alt)]
    num_var_vertexes = length(var_vertexes)
    
    if (num_var_vertexes < 2) {
        stop("Variant is not recurrent.  All variants must be recurrent")
    }
    
    # Initialize mrca to the first vertex
    mrca = var_vertexes[1]
    mrca_ancestors = ancestor_paths[[mrca]]
    for (i in 2:num_var_vertexes) {
        w = var_vertexes[i]
        w_ancestors = ancestor_paths[[w]]
        
        mrca = get_mrca(mrca_ancestors, w_ancestors)
        
        # Update
        mrca_ancestors = ancestor_paths[[mrca]]
    }
    
    # Set all nodes from mrca up to, but not including the root (i.e., [MRCA, root)) to alt
    num_mrca_ancestors = length(mrca_ancestors)
    if (num_mrca_ancestors > 1) {
        V(g)[mrca_ancestors[1:(num_mrca_ancestors - 1)]]$VARIANT = alt
    }
    
    g = add_edge_annotation(g, mrca_ancestors, vartype, ref, vartype_colors_palette, alt=alt)
    return(list(g=g, mrca=mrca))
}

get_root <- function(g, root_name) {
    root_index = which(names(V(g)) == root_name)
    
    return(root_index)
}


get_mrca <- function(path1, path2) {
    # Paths are sorted from leaf to root
    mrca = intersection(path1, path2)[1]
    
    return(mrca)
}

get_leaves <- function(g) {
    leaf_indexes = which(degree(g, V(g), mode="out") == 0)
    leaf_names = names(leaf_indexes)
    
    return(list(leaf_names=leaf_names, leaf_indexes=leaf_indexes))
}

get_children <- function(g, v) {
    children = neighbors(g, v, mode=c("out"))
    return(children)
}

get_profiled_tissue_indexes <- function(g, profiled_tissues) {
    # Find the leaves
    tmp = get_leaves(g)
    leaf_indexes = tmp$leaf_indexes
    leaf_names = tmp$leaf_names
    
    profiled_tissue_indexes = unname(leaf_indexes[leaf_names %in% profiled_tissues])
    
    return(profiled_tissue_indexes)
}


df2uniq_vartype_df <- function(var_df, var_bin_df) { # FIXME delete me?
    var_cat = apply(var_df, 1, function(x) {paste(x, collapse=" ")})
    counts = table(var_cat)
    
    # Expand concatentated column back to dataframe
    var_strs = rownames(counts)
    vars = strsplit(var_strs, " ", fixed=TRUE)
    
    var_compressed_df = as.data.frame(t(as.data.frame(strsplit(var_strs, " ", fixed=TRUE), stringsAsFactors=FALSE)), stringsAsFactors=FALSE)
    rownames(var_compressed_df) = NULL
    colnames(var_compressed_df) = colnames(var_df)
    
    var_compressed_df$count = as.numeric(unname(counts)) # w/o the as.numeric, the values are of type "table"
    
    # Convert "NA" back to NA
    var_compressed_df[var_compressed_df=="NA"] = NA
    
    return(var_compressed_df)
}

add_all_variants_to_tree <- function(donor_sub_g, ancestor_paths, donor, var_compressed_counts, var_compressed_mut_vectors, vartype_colors_palette, base_colors_df, core_vartypes, suffix_sans_extension="", plot_individual_var_tree=FALSE, plot_random_var_tree=0.01, save_file=TRUE, prefix="", vertex_size=5, print_var_info=FALSE, alt="T", var_type_index=NULL) {
    num_edges = length(E(donor_sub_g))
    g_consensus = donor_sub_g
    
    # For each edge, set each of the vartype counters, e.g., T>C to 0
    for (t in core_vartypes) {
        g_consensus = set.edge.attribute(g_consensus, t, value=0)
    }
    
    # Find the root
    root_index = which(names(V(donor_sub_g)) == "Zygote") # An alternative could be in == 0
    
    # Find the leaves
    leaf_indexes = which(degree(donor_sub_g, V(donor_sub_g), mode="out") == 0)
    leaf_names = names(leaf_indexes)
    
    num_used_vartypes = ncol(var_compressed_counts) - 1 # -1 since there is 1 extra column for $count (total count)
    used_vartypes = colnames(var_compressed_counts)[1:num_used_vartypes]
    
    # Create a subset of just the tissues (no reference)
    num_tissues_profiled = ncol(var_compressed_mut_vectors) - 1 # -1 since there is 1 extra column: $ref
    only_var_df = var_compressed_mut_vectors[,1:num_tissues_profiled, drop=FALSE]
    profiled_tissues_only_var_df = only_var_df
    profiled_tissues = colnames(profiled_tissues_only_var_df)
    profiled_tissue_indexes = unname(leaf_indexes[leaf_names %in% profiled_tissues])
    
    # Set the unprofiled tissues to NA
    unprofiled_tissues = leaf_names[! leaf_names %in% colnames(only_var_df)]
    only_var_df[,unprofiled_tissues] = NA

    # Reorder only_var_df to have the same order as the graph
    only_var_df = only_var_df[,match(leaf_names, colnames(only_var_df))]

    num_variants = sum(var_compressed_counts$count)
    
    var_compressed_mrca = data.frame(var_compressed=rownames(var_compressed_counts))
    var_compressed_mrca$mrca = NA
 
    var_type_index_to_use = var_type_index
    for (j in 1:nrow(var_compressed_counts)) {
        
        my_g = donor_sub_g
        
        ref = "C" # Dummy reference
        alt = "T" # Dummy alt
        vartype = sprintf("%s>%s", ref, alt) # Dummy vartype
        
        # Set the vertex variant
        #V(my_g)$VARIANT = ref
        # Set root to refernece
        V(my_g)[root_index]$VARIANT = ref
        
        # Set all the leaves
        V(my_g)[leaf_indexes]$VARIANT = unlist(only_var_df[j,])
        
        # Add variant to tree
        tmp = add_variant_to_tree(my_g, ancestor_paths, ref, vartype, vartype_colors_palette, alt=alt)
        my_g = tmp$g
        mrca = tmp$mrca # Will store the vertex index
        
        var_compressed_mrca$mrca[j] = mrca
        
        if (is.null(var_type_index)) {
            var_type_index_to_use = j
        }
        
        # Plot the variant tree if option specified, or for a random set of the variants (default 1%)
        if ((plot_individual_var_tree) | (runif(1) <= plot_random_var_tree)) {
            title = sprintf("%s: var type %.0f", donor, var_type_index_to_use)
            
            # Since the mutation vectors are getting pretty long, don't put them on the plot
            #subt = paste(c(sprintf("Ref: %s", ref), paste(colnames(profiled_tissues_only_var_df), profiled_tissues_only_var_df[j,], sep=": ")), collapse="\n")
            var_prefix = paste(prefix, var_type_index_to_use, sep=".")
            # Really, this is a var_type graph
            plot_g(my_g, "edge_weight", "categorical_and_numeric", var_prefix, title, suffix_sans_extension=suffix_sans_extension, profiled_tissue_indexes=profiled_tissue_indexes, mrca=mrca, save_file=save_file, vertex_size=vertex_size, set_vertex_color=TRUE, base_colors_df=base_colors_df)
        }
        
        # Update the counters
        for (k in 1:num_used_vartypes) {
            vartype = used_vartypes[k]
            var_freq = var_compressed_counts[j,k]
            counter = get.edge.attribute(g_consensus, vartype)
            counter = counter + (E(my_g)$edge_weight * var_freq)
            g_consensus = set.edge.attribute(g_consensus, vartype, value=counter)
        }
    }
    
    # Calculate the total # of variants for each edge
    # Also calculate the max (used for setting global width/color scales)
    all_variant_types_counter = rep(0, times=num_edges)
    max_num_var = 0
    for (vartype_i in core_vartypes) {
        variant_type_counter = get.edge.attribute(g_consensus, vartype_i)
        
        all_variant_types_counter = all_variant_types_counter + variant_type_counter
        
        # Update max
        max_num_var_i = max(variant_type_counter)
        
        if (max_num_var_i > max_num_var) {
            max_num_var = max_num_var_i
        }
    }
    E(g_consensus)$num_var = all_variant_types_counter

    # Print how many variants
    if (print_var_info) {
        writeLines(sprintf("There were %.0f variants total", num_variants), con=stderr())
    }

    return(list(g_consensus=g_consensus, max_num_var=max_num_var, profiled_tissue_indexes=profiled_tissue_indexes, var_compressed_mrca=var_compressed_mrca))
}

add_edge_annotation <- function(g, mrca_ancestors, vartype, ref, vartype_colors_palette, alt="T") {
    
    # Add the variant change and color to each edge
    # For now, only want to change the unknown edges to purple and the mrca's incident edge to the variant type
    E(g)$edge_color = "black"
    for (e in E(g)) {
        
        v_indexes = get.edges(g, e)
        parent_index = v_indexes[1]
        child_index = v_indexes[2]
        
        parent_allele = V(g)[parent_index]$VARIANT
        child_allele = V(g)[child_index]$VARIANT
        
        if (is.na(parent_allele) | is.na(child_allele)) {
            variant_change = "unknown"
            E(g)[e]$edge_color = vartype_colors_palette[variant_change]
        } else {
            variant_change = sprintf("%s>%s", parent_allele, child_allele) # For edges there there is no change, will get values like A>A.
        }
    }
    
    # Add the edge weights
    E(g)$edge_weight = 0
    
    # New: mutation is added from MRCA to root with weight = 1/N, where N = number of edges
    num_mrca_ancestors = length(mrca_ancestors)
    num_mutation_edges = num_mrca_ancestors - 1
    edge_weight = 1 / num_mutation_edges
    
    for (mrca_ancestor in mrca_ancestors[1:(num_mrca_ancestors - 1)]) {
        # Look up the incidence edge
        edge_w_var = incident(g, mrca_ancestor, mode="in")
        E(g)[edge_w_var]$edge_weight = edge_weight
        E(g)[edge_w_var]$edge_color = vartype_colors_palette[vartype]
    }
    
    return(g)
}


set_edge_attributes_for_plot <- function(g, type, category, value_global_max=NULL, scientific_notation=TRUE, set_labels=TRUE, set_edge_color=TRUE, set_na_colors=FALSE, set_widths=TRUE, diverging_palette=FALSE, na_color="grey", alt_edge_label_attr=NULL, use_log_trans=FALSE, pseudo_count=1) {
    # Set edge color based on how many variants
    # Inspired by http://stackoverflow.com/questions/22255465/assign-colors-to-a-range-of-values.
    
    # Another possibility might be to use colorbrewer palette and then assign color based on percentile of value in data
    if (category == "numeric") {
        g = scale_edges(g, type, value_global_max=value_global_max, scientific_notation = scientific_notation, set_labels=set_labels, set_colors=set_edge_color, set_na_colors = set_na_colors, set_widths = set_widths, diverging_palette=diverging_palette, na_color=na_color, alt_edge_label_attr=alt_edge_label_attr, use_log_trans=use_log_trans, pseudo_count=pseudo_count)
    
    } else if (category == "categorical_and_numeric") {
        
        g = scale_edges(g, type, value_global_max=value_global_max, set_widths=TRUE, set_labels=TRUE, scientific_notation = scientific_notation, set_labels=set_labels, set_colors=set_edge_color, set_na_colors = set_na_colors, set_widths = set_widths, diverging_palette=diverging_palette, na_color=na_color, alt_edge_label_attr=alt_edge_label_attr, use_log_trans=use_log_trans, pseudo_count=pseudo_count)
        
    } else if (category == "em_consensus_var") {
        g = scale_edges(g, type, value_global_max=value_global_max, scientific_notation = scientific_notation, set_labels=set_labels, set_colors=set_edge_color, set_na_colors = set_na_colors, set_widths = set_widths, diverging_palette=diverging_palette, na_color=na_color, alt_edge_label_attr=alt_edge_label_attr, use_log_trans=use_log_trans, pseudo_count=pseudo_count) #sets the width and color
        # Update the edge with the binom p param
        cluster_indexes = c(2,7,1)
        orig_edge_labels = E(g)$edge_label[cluster_indexes] # The numeric weights have already been transformed to pretty strings
        binom_p_labels = as.numeric(E(g)$binom_p[cluster_indexes])
        
        updated_edge_labels = sprintf("%s\np = %.2f", orig_edge_labels, binom_p_labels)
        E(g)$edge_label[cluster_indexes] = updated_edge_labels
    } else if (category == "leave_edges_alone") {
        # Do nothing.  Don't update widths
    } else {
        # no need to set edge colors since they're already defined
        E(g)$edge_width = 8
        E(g)$edge_label = ""
    }

    if (type == "tissue_tree") {
        # Label edges with the edge number
        if (set_labels) {
            E(g)$edge_label = 1:length(E(g))
        }
    }
    return(g)
}


scale_edges <- function(g, type, set_colors=TRUE, set_widths=TRUE, set_labels=TRUE, value_global_max=NULL, use_log_trans=FALSE, print_nonzero_weights_only=TRUE, scientific_notation=TRUE, set_na_colors=FALSE, diverging_palette=FALSE, na_color="grey", alt_edge_label_attr=NULL, pseudo_count=1) {
    values = get.edge.attribute(g, type)
    #print(edge_attr(g))
    raw_values = values
    if (use_log_trans) {
        values = log10(values + pseudo_count)
        warning("WARNING: log transformation created -Inf values.  Setting to NA")
        values[values == -Inf] = NA
    }
    
    # Start min at 0 or if attribute is neg, min value
    # If diverging, must make range symmetric so 0 is in middle
    if (diverging_palette) {
        # Use local scale, else, will use global scale
        if (is.null(value_global_max)) {
            value_global_max = max(abs(values), na.rm=TRUE)
        }
        value_global_min = -1 * value_global_max
    } else {
        # Use local scale, else, will use global scale
        if (is.null(value_global_max)) {
            value_global_max = max(values, na.rm=TRUE)
        }
        # If values dip to negatives, use neg value; else ues 0
        if (min(values, na.rm=TRUE) < 0) {
            value_global_min = min(values, na.rm=TRUE)
        } else {
            value_global_min = 0
        }
    }
    
    range_of_values = value_global_max - value_global_min
    
    if (range_of_values) { # All values are not the same
        values_scaled = (values - value_global_min) / range_of_values
    } else {
        values_scaled = values - value_global_min 
    }
    
    # Edge color
    if (set_colors) {
        #print("raw")
        #print(raw_values)
        #print(raw_values[7])
        #print("scaled")
        #print(values_scaled)
        if (diverging_palette) {
            f_discrete = brewer.pal(11, 'RdYlBu') 
            # Convert to continuous
            # Inspired by https://sites.google.com/site/seascapemodelling/home/r-code/customising-plots/continous-colour-brewer-palettes
            f = colorRamp(f_discrete)
            
        } else {
            f <- colorRamp(viridis(100))
        }
        E(g)$edge_color = NA

        # Handle non NA values
        non_na_edges = which(! is.na(values_scaled))
        colors <- rgb(f(values_scaled[non_na_edges])/255)
        E(g)$edge_color[non_na_edges] = colors
        E(g)$values_scaled = values_scaled
        E(g)$values_used = values
        
        if (set_na_colors) {
            # In older code, I guess weight_unknown was set for NA weights.  If that attribute is set use it, else use the NA values
            if (is.null(E(g)$weight_unknown)) {
                na_edges = which(is.na(values_scaled))
                E(g)[na_edges]$edge_color = na_color
            } else {
                edges_of_interest = which(E(g)$weight_unknown == TRUE)
                E(g)[edges_of_interest]$edge_color = na_color
            }
        }
    }
    # Edge width
    if (set_widths) {
        min_width = 1
        max_width = 40
        E(g)$edge_width = unname(quantile(seq(min_width, max_width), values_scaled)) # Map the range [0,1] to [1,12]
    }
    
    # Edge label
    if (set_labels) {
        if (is.null(alt_edge_label_attr)) {
            if (print_nonzero_weights_only) {
                if (scientific_notation) {
                    edge_labels_pretty = sapply(raw_values, function(x) {if (x == 0 | is.na(x)) {""} else {sprintf("%.2E", x)}})
                } else {
                    edge_labels_pretty = sapply(raw_values, function(x) {if (x == 0 | is.na(x)) {""} else {sprintf("%.2f", x)}})
                }
            } else {
                if (scientific_notation) {
                    edge_labels_pretty = sprintf("%.2E", raw_values)
                } else {
                    edge_labels_pretty = sprintf("%.2f", raw_values)
                }
            }
        } else {
            # Use an alternative edge label
            edge_labels_pretty = get.edge.attribute(g, alt_edge_label_attr)
        }
        E(g)$edge_label = edge_labels_pretty
    }
    #print(edge_labels_pretty)
    #stop()
    return(g)
}

load_tissue_tree <- function(local, name="small", vertex_fn=NULL, edge_fn=NULL) {
    metadata_repo_dir = get_metadata_dir(local)

    # Load the tissue abbreviations
    tissue_abbrev_fn = file.path(metadata_repo_dir, "tissueSiteDetail.20170324.tsv")
    tissue_abbrev_df = read.table(tissue_abbrev_fn, header=TRUE, stringsAsFactors=FALSE, sep="\t")
    
    # Create tree
    if (is.null(vertex_fn) & is.null(edge_fn)) {
        if (name == "full.no_gtex_blacklist.no_sex.no_cell_lines") { 
            # No GTEx blacklisted tissues
            # No sex-specific tissues (e.g., prostate)
            # No cell lines
            vertex_fn = file.path(metadata_repo_dir, "trees", "tissue_lineage_vertexes.no_cycles.no_gtex_blacklist.no_sex.no_cell_lines.tsv")
            edge_fn = file.path(metadata_repo_dir, "trees", "tissue_lineage_edges.no_cycles.no_gtex_blacklist.no_sex.no_cell_lines.tsv")
        } else if (name == "germ_layers") {
            # No GTEx blacklisted tissues
            # No sex-specific tissues (e.g., prostate)
            # No cell lines
            # Tissues are divided by germ layer
            vertex_fn = file.path(metadata_repo_dir, "trees", "tissue_lineage_vertexes.germ_layers.tsv")
            edge_fn = file.path(metadata_repo_dir, "trees", "tissue_lineage_edges.germ_layers.tsv")
        } else {
            stop(sprintf("unknown tree name '%s'", name))
        }
    } else if (is.null(vertex_fn) | is.null(edge_fn)) {
        stop("load_tissue_tree() requires either both or neigther vertex and edge files to be supplied")
    }
    
    vertex_df = read.table(vertex_fn, header=FALSE, stringsAsFactors=FALSE, sep="\t")
    edge_df = read.table(edge_fn, header=FALSE, stringsAsFactors=FALSE, sep="\t")
    
    tissue_g = graph_from_data_frame(edge_df, directed=TRUE, vertices=vertex_df)
    
    # Add the tissue abbreviation
    # First look up the tissues with no abbreviation and add colors for them
    v_tissue_names = V(tissue_g)$name
    
    tissues_w_no_abbrev = v_tissue_names[! v_tissue_names %in% tissue_abbrev_df$tissue_site_detail]
    tissues_w_no_abbrev_df = data.frame(tissue_site_detail=tissues_w_no_abbrev,
                                        tissue_site_detail_abbr=tissues_w_no_abbrev, 
                                        tissue_site_detail_id=tissues_w_no_abbrev, 
                                        tissue_site=tissues_w_no_abbrev, 
                                        tissue_color_hex="D3D3D3", 
                                        tissue_color_rgb="211,211,211")
    
    # Combine
    tissue_abbrev_df = rbind(tissue_abbrev_df, tissues_w_no_abbrev_df)
    
    # Update the germlevel colors
    tissue_abbrev_df[which(tissue_abbrev_df$tissue_site_detail == "Ectoderm"), c("tissue_color_hex", "tissue_color_rgb")] = c("3498db", "52,152,219")
    tissue_abbrev_df[which(tissue_abbrev_df$tissue_site_detail == "Endoderm"), c("tissue_color_hex", "tissue_color_rgb")] = c("f1c40f", "241,196,15")
    tissue_abbrev_df[which(tissue_abbrev_df$tissue_site_detail == "Mesoderm"), c("tissue_color_hex", "tissue_color_rgb")] = c("e74c3c", "231,76,60")
    
    index_order = match(v_tissue_names, tissue_abbrev_df$tissue_site_detail)
    V(tissue_g)$abbrev = tissue_abbrev_df$tissue_site_detail_abbr[index_order]
    
    # Add color
    # Convert hex color to real hex color
    tissue_abbrev_df$tissue_color_hex = paste0("#", tissue_abbrev_df$tissue_color_hex)
    V(tissue_g)$tissue_color = tissue_abbrev_df$tissue_color_hex[index_order]
    
    return(list(tissue_g=tissue_g, tissue_abbrev_df=tissue_abbrev_df))
}


load_donor_var <- function(var_fn) {
    PURINE = get_purines()
    donor = names(var_fn)
    # Step 1: Load a donor's variants
    tmp = read.table(var_fn, stringsAsFactors=FALSE, sep="\t", comment.char="", header=TRUE, row.names=1, check.names=FALSE, colClasses=c("character")) # IMPORTANT: include colClasses argument to prevent R from reading T as TRUE
    colnames(tmp) = gsub("^#", "", colnames(tmp)) # The text file use "^#"" on the first column name.  Now that we're in R space with defined headers, remove this extra annotation
    var_df = as.data.frame(t(tmp), stringsAsFactors=FALSE)
    
    # Add the donor name to the variant infomration since variants from multiple variants will be combined later
    rownames(var_df) = paste(donor, rownames(var_df), sep = ":")
    
    vaf_df = var_df[,1:2, drop=FALSE] # Rows 1 and 2 contain the em and bbcall VAFs, respectivively
    var_df = var_df[,3:ncol(var_df), drop=FALSE] # Rows 3+ contains the variant infomration
    
    donor_tissues = colnames(var_df)
    # Add reference and variant fields
    var_df$ref = sapply(rownames(var_df), function(y) {strsplit(y, ":", fixed=TRUE)[[1]][4]})
    
    # Switch purine-based (no pun intended!) variant calls to pyrimidine-based
    roi = var_df$ref %in% PURINE
    if (any(roi)) { # R doesn't like it if try subset 0 rows (i.e., all FALSES)
        var_df[roi,] = apply(var_df[roi,], 2, seq2seq_compliment)
    }
    
    return(list(var_df=var_df, vaf_df=vaf_df, donor_tissues=donor_tissues))
}


tissue_tree2donor_tree <- function(donor_tissues, tissue_g, tissue_abbrev_df, reserved_tissues) {
    # Step 1: Find the subset of the tissue tree that matches this donor's tissues
    donor_g = tissue_g
    
    # Given a list of donor tissues, create a unique list of the parent vertexes.  Then create a subgraph of these vertexes
    donor_vertexes = c()
    for (my_t in donor_tissues) {
        #print(my_t)
        # Get parent vertex
        parent_v = get_parent_vertex(donor_g, my_t)
        while(parent_v) {
            parent_v = get_parent_vertex(donor_g, my_t)
            
            if (! parent_v) { # Don't add the FALSE elements to donor_vertexes
                break
            }
            donor_vertexes = c(donor_vertexes, parent_v)
            my_t = parent_v
        }
    }
    
    # Create a unique list
    donor_vertexes = unique(c(donor_tissues, names(donor_vertexes)))
    
    # Make a subtree of just these
    # Find the vertex indexes that correspond to these vertex names
    donor_vertex_ids = sapply(donor_vertexes, function(x) {which(V(donor_g)$name == x)})
    donor_sub_g = induced_subgraph(donor_g, donor_vertex_ids, impl="auto")
    
    #plot_tissue_g(donor_sub_g, "")
    
    # Step 2: compress internal nodes to improve clarity
    # Find the root
    root_index = which(names(V(donor_sub_g)) == "Zygote")
    # Find the distance to root for each vertex
    b = bfs(donor_sub_g, root_index, dist=TRUE)
    
    current_vertex_depths = sort(b$dist, decreasing=TRUE)
    max_depth = as.numeric(current_vertex_depths[1])
    
    for (i in max_depth:0) { # Root is at level 0
        #print(sprintf("\n\nOn level %d", i))

        # Find all the vertexes at this level
        v_names = names(which(current_vertex_depths == i))
        
        # Try to compress the each vertex
        for (v_name in v_names) {
            #print(sprintf("Looking at vertex %s", v_name))
            v_index = which(V(donor_sub_g)$name == v_name)
            
            # Find it's parent
            parent_v = get_parent_vertex(donor_sub_g, v_index)
            parent_index = as.numeric(parent_v)
            
            if (parent_v != FALSE) { # the != was needed otherwise R compained that the string could not be interpreted as a logical
                parent_v = names(parent_v) # Get just the anme
                in_degree = degree(donor_sub_g, parent_v, mode="in")
                out_degree = degree(donor_sub_g, parent_v, mode="out")
                
                # Determine if parent is a reserved tissue (i.e., GTEX or germlayer.  Note: some GTEx tissues are internal nodes which may have in degree = out degree = 1 and we don't want to delete these)
                if (parent_v %in% reserved_tissues) {
                    p_is_not_gtex_tissue = FALSE
                } else {
                    p_is_not_gtex_tissue = TRUE
                }
                
                if (((in_degree == 1 ) & (out_degree == 1)) & (p_is_not_gtex_tissue)) {
                    # If the above conditions are true, the parent index can be safely deleted.  (igraph will automatically remove the edges to/from the parent)
                    
                    # Find the grandparent
                    grandparent_v = names(get_parent_vertex(donor_sub_g, parent_v))  # gp will be defined since in_degree of parent is 1
                    
                    # Delete the parent
                    donor_sub_g = delete_vertices(donor_sub_g, parent_index)
                    #print(sprintf("Just deleted parent (%s)", parent_v))
                    
                    # Add an edge between the child and the grandparent (gp)
                    # Since the graph has changed look up the indexes
                    grandparent_index = as.numeric(vertex2index(donor_sub_g, grandparent_v))
                    v_index = as.numeric(vertex2index(donor_sub_g, v_name))
                    
                    donor_sub_g = add_edges(donor_sub_g, c(grandparent_index, v_index))
                    
                    # Update the current_vertex_depths
                    # Delete the child index name
                    index_to_delete = which(names(current_vertex_depths) == v_name)
                    current_vertex_depths = current_vertex_depths[-index_to_delete]
                    
                    # Rename the parent index name to the child index name
                    index_to_update = which(names(current_vertex_depths) == parent_v)
                    names(current_vertex_depths)[index_to_update] = v_name

                }
            } # Else at the root of the tree can merge anyhigher up so end
            #plot_tissue_g(donor_sub_g, "")
        }
    }
   
    # Determine the ~width and height of the tree (to be used to scale size of figures)
    height = max(current_vertex_depths) + 1 # Assign height to the # of levels in the tree
    width = max(table(current_vertex_depths)) # Assign width to the # of nodes in the level that has the most nodes 
    return(list(donor_sub_g=donor_sub_g, height=height, width=width))
}


tissue_tree2donor_tree_recursive <- function(donor_tissues, tissue_g, tissue_abbrev_df, reserved_tissues) {
    # Step 1: Find the subset of the tissue tree that matches this donor's tissues
    donor_g = tissue_g
    
    # Given a list of donor tissues, create a unique list of the parent vertexes.  Then create a subgraph of these vertexes
    donor_vertexes = c()
    for (my_t in donor_tissues) {
        #print(my_t)
        # Get parent vertex
        parent_v = get_parent_vertex(donor_g, my_t)
        while(parent_v) {
            parent_v = get_parent_vertex(donor_g, my_t)
            
            if (! parent_v) { # Don't add the FALSE elements to donor_vertexes
                break
            }
            donor_vertexes = c(donor_vertexes, parent_v)
            my_t = parent_v
        }
    }
    
    # Create a unique list
    donor_vertexes = unique(c(donor_tissues, names(donor_vertexes)))
    
    # Make a subtree of just these
    # Find the vertex indexes that correspond to these vertex names
    donor_vertex_ids = sapply(donor_vertexes, function(x) {which(V(donor_g)$name == x)})
    donor_sub_g = induced_subgraph(donor_g, donor_vertex_ids, impl="auto")
    
    plot_tissue_g(donor_sub_g, "")
    
    # Step 2: compress internal nodes to improve clarity
    # Find the root
    root_index = which(names(V(donor_sub_g)) == "Zygote")
    
    # Find all the leaves
    leaves = V(donor_sub_g)[degree(donor_sub_g, mode="out") == 0]
    leaf_names = names(leaves) # Hold on your hats here, the graph will be changing so get the vertex names since the indexes will be modified below
    
    # THIS DOESN'T WORK BECAUSE ANCESTORS OF CHILDREN CAN BE COMPRESSED WHILE CHILDREN CANNOT BE FURTHER COMPRESSED
    
    v_name = leaf_names[4]
    for (v_name in leaf_names) {
        # Look up the vertex's index
        v_index = which(V(donor_sub_g)$name == v_name)
        
        donor_sub_g = try_to_compress(donor_sub_g, v_name)
                    
        plot_tissue_g(donor_sub_g, "")
    }
    
    return(donor_sub_g)
}

try_to_compress <- function(g, v_name) {
    print(sprintf("Working on %s", v_name))
    # Look up the vertex's index
    v_index = which(V(g)$name == v_name)
    
    # Find it's parent
    parent_v = get_parent_vertex(g, v_index)
    parent_index = as.numeric(parent_v)
     
    
    if (parent_v != FALSE) { # the != was needed otherwise R compained that the string could not be interpreted as a logical
        parent_v = names(parent_v) # Get just the anme
        in_degree = degree(g, parent_v, mode="in")
        out_degree = degree(g, parent_v, mode="out")
        
        # Determine if parent is a gtex tissue.  (Some gtex tissues are internal nodes which may have in degree = out degree = 1 and we don't want to delete these)
        if (parent_v %in% reserved_tissues) {
            p_is_not_gtex_tissue = FALSE
        } else {
            p_is_not_gtex_tissue = TRUE
        }
        
        if (((in_degree == 1 ) & (out_degree == 1)) & (p_is_not_gtex_tissue)) {
            # If the above conditions are true, the parent index can be safely deleted.  (igraph will automatically remove the edges to/from the parent)
            
            # Find the grandparent
            grandparent_v = names(get_parent_vertex(g, parent_v))  # gp will be defined since in_degree of parent is 1
            
            # Delete the parent
            g = delete_vertices(g, parent_index)
            print(sprintf("Just deleted parent (%s)", parent_v))
            
            # Add an edge between the child and the grandparent (gp)
            # Since the graph has changed look up the indexes
            grandparent_index = as.numeric(vertex2index(g, grandparent_v))
            v_index = as.numeric(vertex2index(g, v_name))
            
            g = add_edges(g, c(grandparent_index, v_index))
            plot_tissue_g(g, "")
            print(sprintf("Want to compress current node again (%s)", v_name))
            g = try_to_compress(g, v_name)
        } else { # The child node can't be merged anymore, but perhaps its parent can be
            print(sprintf("The child (%s) node can't be merged anymore, but perhaps its parent (%s) can be", v_name, parent_v))
            g = try_to_compress(g, parent_v)
        }
    } # Else at the root of the tree can merge anyhigher up so end
    
    return(g)
}

# Build index of leaf to root paths for all leaves
get_ancestor_paths <- function(g) {
    num_vertexes = length(V(g))
    ancestor_paths = vector("list", num_vertexes)
    for (v in V(g)) {
        v_ancestors = subcomponent(g, v, mode="in")
        ancestor_paths[[v]] = v_ancestors
    }
    return(ancestor_paths)
}

get_node2germlayer <- function(g, ancestor_paths, germlayers=c("Mesoderm", "Endoderm", "Ectoderm")) {
    num_vertexes = length(V(g))
    node2germlayer = vector("character", num_vertexes)
    for (path in ancestor_paths) {
        v = path[1]
        hits = path[(names(path) %in% germlayers)]
        if (length(hits) == 0) {
            # If no hit, report what node is (either gastrual or zygote)
            node2germlayer[v] = names(v)
        } else if (length(hits) == 1) {
            node2germlayer[v] = names(hits)[1]
        } else {
            writeLines("ERROR: got multiple germlayers for a node.  This should not be possible", con=stderr())
            stop()
        }
    }
    return(node2germlayer)
}

get_descendants <- function(g, v) {
    v_descendants = subcomponent(g, v, mode="out") # Includes itself
    
    return(v_descendants)
}

get_descendant_leaves <- function(g, v, leaves) {
    v_descendants = get_descendants(g, v)
    #descendant_leaf_indexes = intersect(names(v_descendants), leaves)
    descendant_leaf_indexes = intersect(v_descendants, leaves)
    
    return(descendant_leaf_indexes)
}

mk_descendent_trees <- function(tissue_g, return_lookup=FALSE) {
    # descendent_tree: key = mutation vector string; value = tree with edge marked
    # mutation_vector_str2edge = dataframe 
    
    
    mutation_prob_attr = "mutation_prob"
    # By definition, in a descendent tree, the variant is mapped to a single edge, and thus has an edge weight of 1
    edge_weight = 1
    
    leaves = get_leaves(tissue_g)
    leaf_names_alpha_sorted = sort(leaves$leaf_names)
    leaf_edge_indexes = as.numeric(incident_edges(tissue_g, leaves$leaf_indexes, mode="in"))
    non_leaf_edge_indexes = setdiff(E(tissue_g), leaf_edge_indexes)
    
    num_descendent_trees = length(non_leaf_edge_indexes)
    descendent_trees = list()
    mutation_vector_str2edge = data.frame(
        mutation_vector=character(length=num_descendent_trees),
        head_node=character(length=num_descendent_trees),
        tail_node=character(length=num_descendent_trees),
        edge_index=character(length=num_descendent_trees),
        pretty_edge=character(length=num_descendent_trees),
        stringsAsFactors=FALSE)
    
    for (i in 1:num_descendent_trees) {
        edge = non_leaf_edge_indexes[i]
        g = tissue_g
        g = set.edge.attribute(g, name=mutation_prob_attr, value = 0)
        g = set.edge.attribute(g, name=mutation_prob_attr, index=edge, value=edge_weight)
        
        # Calculate the mutation vector (i.e., the alphabetically sorted list of tissues with the variant status)
        head_node = head_of(tissue_g, edge)
        tail_node = tail_of(tissue_g, edge)
        toi_indexes = get_descendant_leaves(tissue_g, head_node, leaves$leaf_indexes)
        toi = sort(names(leaves$leaf_indexes)[match(toi_indexes, leaves$leaf_indexes)])
        mutation_vector = rep(0, length(leaf_names_alpha_sorted))
        names(mutation_vector) = leaf_names_alpha_sorted
        mutation_vector[toi] = 1
        mutation_vector_str = paste(mutation_vector, collapse=",")
            
        descendent_trees[[mutation_vector_str]] = g
        
        mutation_vector_str2edge$mutation_vector[i] = mutation_vector_str
        mutation_vector_str2edge$head_node[i] = V(tissue_g)$name[head_node]
        mutation_vector_str2edge$tail_node[i] = V(tissue_g)$name[tail_node]
        mutation_vector_str2edge$edge[i] = edge
        mutation_vector_str2edge$pretty_edge[i] = sprintf("%s->%s", V(tissue_g)$name[tail_node], V(tissue_g)$name[head_node])
    }
    
    # Double-check that the mutation vectors are all unique
    # If they're not unique, then the tree hasn't been pruned enough, e.g., there might be parent with only 1 child
    mutation_vector_freq = as.data.frame(table(names(descendent_trees)))
    if (max(mutation_vector_freq$Freq) > 1) {
        stop("non unique mutation vectors; tree most likely needs to be pruned, e.g., is there a parent with only 1 child?")
    }
    
    if (return_lookup) {
        # This is mainly a backwards compatible effort.  Added mutation_vector_str2edge at a later date.
        return(list(descendent_trees=descendent_trees, mutation_vector_str2edge=mutation_vector_str2edge))
    } else {
        return(descendent_trees) # List of trees; each tree has mutation_prob_attr set to 1 where the mutation occurred
    }
}

mk_donor_g_from_var_df_wrapper <- function(var_compressed_counts, var_compressed_mut_vectors, num_var_list, donor_id, tissue_g, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension=NULL, plot_total_num_var_tree=TRUE, plot_num_var_type_tree=TRUE, plot_individual_var_tree=FALSE, save_file=FALSE, vertex_size=5, toi=NULL, print_var_info=FALSE, var_type_index=NULL) {
    
    # Build index of leaf to root paths for all leaves
    ancestor_paths = get_ancestor_paths(tissue_g)
    
    tmp = add_all_variants_to_tree(tissue_g, ancestor_paths, donor_id, var_compressed_counts, var_compressed_mut_vectors, vartype_colors_palette, base_colors_df, core_vartypes=core_vartypes, plot_individual_var_tree=plot_individual_var_tree, suffix_sans_extension=suffix_sans_extension, prefix=donor_id, save_file=save_file, vertex_size=vertex_size, print_var_info=print_var_info, var_type_index=var_type_index)
    g_consensus = tmp$g_consensus
    max_num_var = tmp$max_num_var
    profiled_tissue_indexes = tmp$profiled_tissue_indexes
    var_compressed_mrca = tmp$var_compressed_mrca
    
    # Total number of variants on tree
    if (plot_total_num_var_tree) {
        plot_g(g_consensus, type="num_var", category="numeric", title=sprintf("%s: total # var consensus tree", donor_id), prefix=donor_id, suffix_sans_extension=suffix_sans_extension, profiled_tissue_indexes=profiled_tissue_indexes, save_file=save_file, vertex_size=vertex_size)
    }
    
    # Number of each type of variant, e.g., C>T on tree
    if (plot_num_var_type_tree) {
        plot_g_num_var_type(g_consensus, core_vartypes, max_num_var, donor_id, profiled_tissue_indexes, num_var_list, prefix=donor_id, suffix_sans_extension=suffix_sans_extension, save_file=save_file, vertex_size=vertex_size)
    }
    
    return(list(g_consensus=g_consensus, var_compressed_mrca=var_compressed_mrca))
}

# Copied from https://gist.githubusercontent.com/ChrKoenig/2289006de68f26f04e3e7fa1eb7a683c/raw/a80e2fbf6b0fed715bff2ce0d50c8be6d661b286/graph_to_newick.R.  This function was mentioned in this SO post (https://stackoverflow.com/questions/40047248/layout-igraph-object-as-fan-r) for converting igraph to newick
# Initially, I didn’t think it was producing the correct output on my data. However, I later realized that some of my node names contain "()" which is most likely why the visualization wasn’t working and led me to the probably incorrect conclusion that it might not be working. Since the one I wrote seems to be working, I’ll continue to use it. Plus, I understand how it works.
graph_to_newick = function(graph, root) {
    ##### Define callback functions for DFS
    # 1. function to be called whenever a vertex is visited
    f.in <- function(graph, data, extra) {
        curr_node = extra$names[data['vid']+1] # Get vertex name (Add 1 to index because igraph uses 0-based indexing)
        prev_node = extra$order.in[which(extra$order.in == curr_node)-1]
        next_node = extra$order.out[which(extra$order.out == curr_node)+1]
        if(length(extra$distances[prev_node]) == 0){ # first node / root
            cat("")
        } else{
            if(extra$distances[prev_node] < extra$distances[curr_node]){cat("(")}
        }
        FALSE # Do not terminate
    }
    
    # 2. function to be called whenever the subtree of a vertex is completed
    f.out <- function(graph, data, extra) {
        curr_node = extra$names[data['vid']+1] # Get vertex name (Add 1 to index because igraph uses 0-based indexing)
        prev_node = extra$order.in[which(extra$order.in == curr_node)-1]
        next_node = extra$order.out[which(extra$order.out == curr_node)+1]
        if(length(extra$distances[prev_node]) == 0){ # first node / root
            cat("")
        } else if(is.na(next_node)){ # last / root
            cat(curr_node)
            cat(";")
        } else{
            cat(curr_node)
            if(extra$distances[next_node] < extra$distances[curr_node]){cat(")")}
            if(extra$distances[next_node] == extra$distances[curr_node]){cat(",")}
            if(extra$distances[next_node] > extra$distances[curr_node]){cat(",")}
        }
        FALSE # Do not terminate
    }
    
    ##### Instructions
    # Obtain extra arguments (distance from root, order of node visits etc.) by running a test DFS without callback functions
    tmp_dfs = graph.dfs(graph, root = which(names(V(graph)) == root), dist = T, order = T, order.out = T) 
    # Organize results as extra arguments for the actual DFS
    extra = list(names  = names(V(graph)), order.in = names(tmp_dfs$order), order.out = names(tmp_dfs$order.out), 
                 distances = tmp_dfs$dist, maxdist = max(tmp_dfs$dist))
    # Run DFS with callback functions
    newick_string = capture.output(graph.dfs(graph, root = which(names(V(graph)) == root), 
                                             dist = T, father = T, in.callback = f.in, out.callback = f.out, extra = extra))[1]
    # Tidy up
    newick_string = gsub(pattern = "\\$root", paste(root, ";", sep  = ""), newick_string)
    return(newick_string)
}

igraph2newick <- function(g, root_name) {
    root_index = which(names(V(g)) == root_name)
    # IMPORTANT: scrub tissue names.  Some gtex samples use () in them which is used for the newick format!!
    #names(V(g)) = unlist(lapply(names(V(g)), sanitize_var_for_fn))
    g = set.vertex.attribute(g, "name", value=unlist(lapply(names(V(g)), sanitize_var_for_fn)))
    
    b = bfs(g, root_index, dist=TRUE)
    current_vertex_depths = sort(b$dist, decreasing=TRUE)
    max_depth = as.numeric(current_vertex_depths[1])
    
    for (i in max_depth:1) { # Root is at level 0
        
        # Find all the vertexes at this level
        v_names = names(which(current_vertex_depths == i))
        
        parent_names = c()
        for (v_name in v_names) {
            # Find it's parent
            v_index = which(V(g)$name == v_name)
            parent_v = get_parent_vertex(g, v_index)
            parent_index = as.numeric(parent_v)
            
            if (is_na_or_false(V(g)[parent_index]$children)) {
                # If the child has children, use the string version (containing it and its children)
                if (is_na_or_false(V(g)[v_index]$formatted)) {
                    # Child does not have any children, use the v_name
                    vertex_formatted = v_name
                } else {
                    # Child has children, use the string version
                    vertex_formatted = V(g)[v_index]$formatted
                }
                V(g)[parent_index]$children = vertex_formatted
                
            } else {
                # If the child has children, use the string version (containing it and its children)
                if (is_na_or_false(V(g)[v_index]$formatted)) {
                    # Child does not have any children, use the v_name
                    vertex_formatted = v_name
                } else {
                    # Child has children, use the string version
                    vertex_formatted = V(g)[v_index]$formatted
                }
                current_children = sprintf("%s,%s", V(g)[parent_index]$children, vertex_formatted)
                V(g)[parent_index]$children = current_children
            }
            parent_names = c(parent_names, names(parent_v))
        }
        
        # Format
        for (parent_name in unique(parent_names)) {
            # Find the grandchild
            if (i==1) {
                V(g)[parent_name]$formatted = sprintf("(%s)%s;", V(g)[parent_name]$children, parent_name)
            } else {
                V(g)[parent_name]$formatted = sprintf("(%s)%s", V(g)[parent_name]$children, parent_name)
            }
        }
    }
    
    newick_str = V(g)[root_index]$formatted
    
    return(newick_str)
}


#
# Plotting functions
#

# Plot the # of variants along each edge on a tree.  Edge width and color are proportional to # of variants along edge.  
# If global max is given, the width and color scales will be defined by the global max.  Otherwise, the scales will be determined using the max of the tree's data.
plot_g_num_var <- function(g, title=NULL, attribute_name=NULL, value_global_max=NULL, prefix="", suffix_sans_extension="", save_file=FALSE, height=NULL, width=NULL, vertex_size=5) {
    
    if (! is.null(attribute_name)) {
        g = set_edge_attributes_for_plot(g, attribute_name, value_global_max, scientific_notation = TRUE)
    }
    if (save_file) {
        attribute_name_pretty = attribute_name
        if (prefix != "") {
            prefix = paste0(prefix, ".")
        }
        if (suffix_sans_extension != "") {
            suffix_sans_extension = paste0(".", suffix_sans_extension)
        }
        plot_fn = sprintf("%s%snum_var_tree%s.png", prefix, attribute_name_pretty, suffix_sans_extension)
        
        plot_fn = sanitize_var_for_fn(plot_fn)
        pdf(file=plot_fn, height=40, width=40)
    }

    plot.igraph(g, 
        layout=layout_as_tree, 
        vertex.color=V(g)$tissue_color, 
        vertex.size=vertex_size,
        vertex.label=V(g)$abbrev,
        vertex.label.cex=0.8, 
        vertex.label.color="black",
        vertex.frame.color=V(g)$tissue_color,
        vertex.label.family="sans",
        vertex.label.font=1,
        edge.color=E(g)$edge_color,
        edge.width=E(g)$edge_width,
        edge.arrow.mode=0, # Turn off arrows since edge widths don't always play well with arrow size
        edge.label=E(g)$edge_label,
        edge.label.color="black",
        edge.label.family="sans",
        edge.label.font=1,
        main=title)

    if (save_file) {
        dev.off()
        write_new_file_message(plot_fn)
    }
}


# Add the triangle shape
# Code copied from igraph's shape() function
mytriangle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
    }
    vertex.size <- 1/200 * params("vertex", "size")
    
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    
    symbols(x=coords[,1], 
            y=coords[,2], 
            bg=vertex.color,
            stars=cbind(vertex.size, 
            vertex.size, 
            vertex.size),
            add=TRUE, 
            inches=FALSE)
}

add_additional_graph_shapes <- function() {

    add_shape("triangle", clip=shapes("circle")$clip,
              plot=mytriangle)
}

plot_g <- function(g, type, category, prefix=NULL, title=NULL, value_global_max=NULL, subt=NULL, suffix_sans_extension="", profiled_tissue_indexes=NULL, mrca=NULL, save_file=TRUE, vertex_size=5, set_vertex_color=FALSE, base_colors_df=NULL, fig_width=40, fig_height=40, vertex_label_cex=0.8, edge_label_cex=0.8, output_dir=NULL, scientific_notation=TRUE, set_edge_labels=TRUE, set_vertex_labels=TRUE, set_edge_color=TRUE, set_na_colors=FALSE, file_format="pdf", mk_plot=TRUE, transparent_bg=FALSE, set_widths=TRUE, diverging_palette=FALSE, na_color="grey", alt_edge_label_attr=NULL, edge_label_color="red", use_log_trans=FALSE, pseudo_count=1) {
    # subt = subtitle
    # save_file = true, then save plot to pdf; else display plot
    if (is.null(output_dir)) {
        output_dir = getwd()
    }
    g = set_edge_attributes_for_plot(g, type, category, value_global_max, scientific_notation = scientific_notation, set_labels=set_edge_labels, set_edge_color=set_edge_color, set_na_colors=set_na_colors, set_widths=set_widths, diverging_palette=diverging_palette, na_color=na_color, alt_edge_label_attr=alt_edge_label_attr, use_log_trans=use_log_trans, pseudo_count=pseudo_count)

    if (suffix_sans_extension != "") {
        suffix_sans_extension = paste0(".", suffix_sans_extension)
    }
    # Use squares with black borders to represent profiled tissues
    V(g)$vertex_border_color = V(g)$tissue_color
    V(g)$vertex_fill_color = V(g)$tissue_color
    V(g)$vertex_shape = "circle"
    V(g)$label_color = "black"
    if (! is.null(profiled_tissue_indexes)) {
        V(g)[profiled_tissue_indexes]$vertex_border_color = "black"
        V(g)[profiled_tissue_indexes]$vertex_shape = "square"
    }
    
    # Annotate the mrca
    if (! is.null(mrca)) {
        V(g)[mrca]$vertex_shape = "sphere"
    }
    
    if (set_vertex_color) {
        if (is.null(base_colors_df)) {
            base_colors_df = load_base_colors()
        }
        bases = V(g)$VARIANT
        bases[is.na(bases)] = "unknown"
        
        colors = base_colors_df$color_code[match(bases, base_colors_df$base)] # https://stackoverflow.com/questions/35636315/replace-values-in-a-dataframe-based-on-lookup-table

        V(g)$vertex_border_color = colors
        V(g)$vertex_fill_color = colors
        #V(g)$label_color = border_colors

    }

    if (set_vertex_labels == FALSE) {
        V(g)$abbrev = ""
    }
    if (save_file) {
        if (is.null(prefix)) {
            plot_fn = sprintf("%s%s", type, suffix_sans_extension)    
        } else {
            plot_fn = sprintf("%s.%s%s", prefix, type, suffix_sans_extension)
        }
        plot_fn = sanitize_var_for_fn(plot_fn)
        plot_fn = file.path(output_dir, plot_fn)
        if (file_format == "pdf") {
            plot_fn = sprintf("%s.pdf", plot_fn)
            pdf(file=plot_fn, width=fig_width, height=fig_height)
        } else {
            plot_fn = sprintf("%s.tiff", plot_fn)
            tiff(plot_fn, units="in", width=fig_width, height=fig_height, res=300)
        }
    }
    if (mk_plot) {
        if (transparent_bg) {
            par(bg=NA)
        }
        plot.igraph(g, 
                    layout=layout_as_tree, 
                    vertex.color=V(g)$vertex_fill_color, 
                    vertex.size=vertex_size,
                    vertex.label=V(g)$abbrev,
                    vertex.label.cex=vertex_label_cex, 
                    vertex.frame.color=V(g)$vertex_border_color,
                    vertex.shape=V(g)$vertex_shape,
                    vertex.label.family="sans",
                    vertex.label.font=1,
                    vertex.label.color = V(g)$label_color,
                    edge.color=E(g)$edge_color,
                    edge.width=E(g)$edge_width,
                    edge.arrow.mode=0, # Turn off arrows since edge widths don't always play well with arrow size
                    edge.label=E(g)$edge_label,
                    edge.label.color=edge_label_color,
                    edge.label.family="sans",
                    edge.label.font=1,
                    edge.label.cex=edge_label_cex,
                    main=title,
                    sub=subt)
        timestamp_figure(cex=0.7)
        if (save_file) {
            dev.off()
            write_new_file_message(plot_fn)
        }
    }
    
    return(g)
}

plot_g_var_type <- function(g, title=NULL, attribute_name=NULL, value_global_max=NULL, subt=NULL, suffix_sans_extension="", prefix="") {
    
    if (! is.null(attribute_name)) {
        g = set_edge_attributes_for_plot(g, attribute_name, value_global_max, scientific_notation = TRUE)
    } else {
        attribute_name = ""
    }

    attribute_name_pretty = attribute_name
    if (prefix != "") {
        prefix = paste0(prefix, ".")
    }
    if (suffix_sans_extension != "") {
        suffix_sans_extension = paste0(".", suffix_sans_extension)
    }
    plot_fn = sprintf("%s%s.var_type_tree%s.png", prefix, attribute_name_pretty, suffix_sans_extension)
    plot_fn = sanitize_var_for_fn(plot_fn)
    png(plot_fn, width = 700, height = 700) # I couldn't get pdf to work with the sans serif font.

    plot <- plot.igraph(g, 
        layout=layout_as_tree, 
        vertex.color=V(g)$tissue_color, 
        vertex.size=25,
        vertex.label=V(g)$abbrev,
        vertex.label.cex=0.8, 
        vertex.label.color="black",
        vertex.frame.color=V(g)$tissue_color,
        vertex.label.family="sans",
        vertex.label.font=1,
        edge.color=E(g)$edge_color,
        edge.width=2,
        main=title,
        sub=subt)

    dev.off()
    write_new_file_message(plot_fn)
}

# Plot the GTEx tissue tree
plot_tissue_g <- function(g, title=NULL) {
    plot.igraph(g, 
        layout=layout_as_tree, 
        vertex.color=V(g)$tissue_color, 
        vertex.label=V(g)$abbrev,
        vertex.size=5, 
        vertex.label.cex=0.8, 
        vertex.label.color="black",
        vertex.frame.color=V(g)$tissue_color,
        vertex.label.family="sans",
        vertex.label.font=1,
        main=title)
}


plot_g_num_var_type <- function(g_consensus, core_vartypes, max_num_var, profiled_tissue_indexes, num_var_list, donor=NULL, prefix=NULL, suffix_sans_extension=NULL, save_file=FALSE, vertex_size=5, fig_height=40, fig_width=40, vertex_label_cex=0.8) {
    # Use the global max for setting the color and edge width scales since we will want to compare different variant types
    for (variant_type in core_vartypes) {
        if (is.null(donor)) {
            my_title = sprintf("%s consensus tree\nn = %.0f variants", variant_type, num_var_list[variant_type])
        } else {
            my_title = sprintf("%s: %s consensus tree\nn = %.0f variants", donor, variant_type, num_var_list[variant_type])
        }
        plot_g(g_consensus, variant_type, "numeric", prefix=prefix, title=my_title, value_global_max=max_num_var, suffix_sans_extension=suffix_sans_extension, profiled_tissue_indexes=profiled_tissue_indexes, save_file=save_file, vertex_size=vertex_size, fig_height=fig_height, fig_width=fig_width, vertex_label_cex=vertex_label_cex)
    }
}


mk_complete_mutation_vector <- function(i, toi, possible_profile_vectors, recurrent_only=TRUE) {
    
    profiled_tissues = toi[possible_profile_vectors[i,] == 1]
    num_profiled_tissues = length(profiled_tissues)
    
    var_df = expand.grid(rep(list(0:1), num_profiled_tissues))
    # Remove the all ref row
    var_df = as.data.frame(var_df[apply(var_df, 1, sum) > 0,]) # as.data.frame needed when there is only 1 column in the df and R converts it to a vector
    
    # Filter for recurrent rows (if applicable)
    if (recurrent_only) {
        var_df = as.data.frame(var_df[rowSums(var_df) > 1,])
        if (nrow(var_df) == 0) {
            var_df = NULL
        }
    }
    
    if (! is.null(var_df)) {
        colnames(var_df) = profiled_tissues
        
        # Add the ref on the end
        var_df$ref = 0
        
        # Add dummy variant alleles.  Use C>T mutation
        var_df[var_df == 0] = "C"
        var_df[var_df == 1] = "T"
    }
    
    return(var_df)
}

visualize_o_m_p_trees <- function(info_df_list, prefix=NULL, suffix_sans_extension="", save_file=FALSE) { # O M P = observed, missing, possible
    info_df_melt = info_df_list$info_df_melt
    info_df_t = info_df_list$info_df_t
    
    # Edge weight distribution
    plot_edge_weight_distribution(info_df_melt, save_file=save_file, prefix=prefix, suffix_sans_extension=suffix_sans_extension, fig_width=6, fig_height=4)
    
    # Missing vs. observed
    plot_edge_type1_vs_edge_type2(info_df_t, "observed", "missing", "observed weight", "missing weight", save_file=save_file, prefix=prefix, suffix_sans_extension=suffix_sans_extension, fig_width=4, fig_height=4)
    plot_edge_weight_ratio(info_df_t, "observed", "missing", "observed weight", "missing weight", save_file=save_file, prefix=prefix, suffix_sans_extension=suffix_sans_extension, fig_width=6, fig_height=4)
    
    # Possible vs. observed
    plot_edge_type1_vs_edge_type2(info_df_t, "observed", "possible", "observed weight", "possible weight", save_file=save_file, prefix=prefix, suffix_sans_extension=suffix_sans_extension, fig_width=4, fig_height=4)
    plot_edge_weight_ratio(info_df_t, "observed", "possible", "observed weight", "possible weight", save_file=save_file, prefix=prefix, suffix_sans_extension=suffix_sans_extension, fig_width=6, fig_height=4)
}

plot_edge_weight_distribution <- function(info_df_melt, save_file=FALSE, prefix=NULL, suffix_sans_extension="", fig_width=6, fig_height=4) {
    
    # OMP weights on each edge
    p = ggplot(info_df_melt, aes(edge, weight)) +
        geom_col(aes(fill=type), position="dodge") +
        theme_bw() +
        xlab("Edge") +
        ylab("Weight") +
        ggtitle("Edge weights") +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=12),
              strip.text=element_text(size=12),
              title=element_text(size=12),
              plot.title = element_text(hjust = 0.5)
        )
    
    
    if (save_file) {
        p2 = timestamp_ggfigure(p)
        if (is.null(prefix)) {
            plot_fn = sprintf("edge_weight_distribution.%s.pdf", suffix_sans_extension)    
        } else {
            plot_fn = sprintf("%s.edge_weight_distribution.%s.pdf", prefix,  suffix_sans_extension)
        }
        
        save_fig(p2, plot_fn, width=fig_width, height=fig_height)
    } else {
        grid.newpage()
        p2 = timestamp_ggfigure(p)
        grid.draw(p2)
    }
}

plot_edge_type1_vs_edge_type2 <- function(df, x_var, y_var, my_x_lab, my_y_lab, save_file=FALSE, prefix=NULL, suffix_sans_extension="", fig_width=6, fig_height=4) {
    point_labels = seq(1, nrow(df))
    
    p = ggplot(df, aes_string(x=x_var, y=y_var, label=point_labels)) +
        geom_point(size=3, col="blue", alpha=0.7) +
        geom_abline(slope=1, intercept=0, linetype=3) +
        geom_text(hjust = 0, nudge_x = 0.01) +
        theme_bw() +
        xlab(my_x_lab) +
        ylab(my_y_lab) +
        ggtitle(sprintf("%s vs. %s", my_y_lab, my_x_lab)) +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=12),
              strip.text=element_text(size=12),
              title=element_text(size=12),
              plot.title = element_text(hjust = 0.5)
        )
    
    if (save_file) {
        p2 = timestamp_ggfigure(p)
        y_var_for_fn = sanitize_var_for_fn(y_var)
        x_var_for_fn = sanitize_var_for_fn(x_var)
        if (is.null(prefix)) {
            plot_fn = sprintf("%s_vs_%s.%s.pdf", y_var_for_fn, x_var_for_fn, suffix_sans_extension)    
        } else {
            plot_fn = sprintf("%s.%s_vs_%s.%s.pdf", prefix, y_var_for_fn, x_var_for_fn, suffix_sans_extension)
        }
        
        save_fig(p2, plot_fn, width=fig_width, height=fig_height)
    } else {
        grid.newpage()
        p2 = timestamp_ggfigure(p)
        grid.draw(p2)
    }
}


plot_edge_weight_ratio <- function(df, x_var, y_var, x_var_pretty, y_var_pretty, save_file=FALSE, prefix=NULL, suffix_sans_extension="", fig_width=6, fig_height=4) {
    
    ratio_df = data.frame(edge=seq(1, nrow(df)), ratio=df[,x_var] / df[,y_var])
    
    ratio_name = sprintf("%s / %s", x_var_pretty, y_var_pretty)
    
    p = ggplot(ratio_df, aes(edge, ratio)) +
        geom_col() +
        geom_hline(yintercept=1, linetype=3) +
        theme_bw() +
        xlab("edge") +
        ylab(ratio_name) +
        ggtitle(sprintf("Edge %s ratios", ratio_name)) +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=12),
              strip.text=element_text(size=12),
              title=element_text(size=12),
              plot.title = element_text(hjust = 0.5)
        )
    
    if (save_file) {
        p2 = timestamp_ggfigure(p)
        y_var_for_fn = sanitize_var_for_fn(y_var)
        x_var_for_fn = sanitize_var_for_fn(x_var)
        if (is.null(prefix)) {
            plot_fn = sprintf("%s_vs_%s_ratio.%s.pdf", x_var_for_fn, y_var_for_fn, suffix_sans_extension)   
        } else {
            plot_fn = sprintf("%s.%s_vs_%s_ratio.%s.pdf", prefix, x_var_for_fn, y_var_for_fn, suffix_sans_extension) 
        }
        
        save_fig(p2, plot_fn, width=fig_width, height=fig_height)
    } else {
        grid.newpage()
        p2 = timestamp_ggfigure(p)
        grid.draw(p2)
    }
}

plot_recurrence_vs_time_plot_cmd <- function(merge_df, save_file=FALSE, prefix=NULL, suffix_sans_extension="", fig_width=6, fig_height=4, log_scale=TRUE, my_title="EM VAF vs. distance from tree root") {
    # Use em VAF
    
    max_x = max(merge_df$dist_from_root, na.rm = TRUE) # There can be NAs if all TP variants were blacklisted
    n_x = seq(1,max_x)
    merge_df$my_x_bins = factor(merge_df$dist_from_root, levels=n_x)
    
    n_per_bin_table = table(merge_df$my_x_bins)
    n_per_bin = unname(n_per_bin_table)
    
    max_y = max(merge_df$median_vaf_em) * 1.05 # Add fudge factor so the # of samples label is not cutoff
    
    p = ggplot(merge_df, aes(my_x_bins, median_vaf_em)) +
        geom_boxplot(fill="grey", outlier.shape=NA) +
        geom_jitter(alpha=0.1) +
        #geom_violin(fill="grey") + # Redo the cut_width so the factor labels are nce
        theme_bw() +
        xlab("Distance from root") +
        ylab("EM vaf") +
        ggtitle(my_title) +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=12),
              strip.text=element_text(size=12),
              title=element_text(size=12),
              plot.title = element_text(hjust = 0.5)
        ) +
        annotate("text",
                 x=n_x,
                 y=max_y,
                 label=sprintf("%.0f", n_per_bin),
                 size=4) +
        expand_limits(y=max_y)
    if (log_scale) {
        p = p +
            scale_y_log10() +
            annotation_logticks(sides="l")
    }
    
    if (save_file) {
        p2 = timestamp_ggfigure(p)
        if (is.null(prefix)) {
            plot_fn = sprintf("em_vaf_vs_dist_from_root.%s.pdf", suffix_sans_extension)    
        } else {
            plot_fn = sprintf("%s.em_vaf_vs_dist_from_root.%s.pdf", prefix,  suffix_sans_extension)
        }
        
        save_fig(p2, plot_fn, width=fig_width, height=fig_height)
    } else {
        grid.newpage()
        p2 = timestamp_ggfigure(p)
        grid.draw(p2)
    }
}

plot_recurrence_vs_embyronic_stage_plot_cmd <- function(merge_df, save_file=FALSE, prefix=NULL, suffix_sans_extension="", fig_width=6, fig_height=4, log_scale=TRUE) {
    # Use em VAF
    
    # max_x = max(merge_df$ts_start, na.rm = TRUE) # There can be NAs if all TP variants were blacklisted
    # n_x = seq(1,max_x)
    # merge_df$my_x_bins = factor(merge_df$ts_start, levels=n_x)
    
    n_per_bin_table = table(merge_df$ts_start)
    n_x = names(n_per_bin_table)
    n_per_bin = unname(n_per_bin_table)
    
    max_y = max(merge_df$median_vaf_em) * 1.05
    p = ggplot(merge_df, aes(factor(ts_start), median_vaf_em)) +
        # geom_boxplot(fill="grey") +
        geom_boxplot(fill="grey", outlier.shape=NA) +
        geom_jitter(alpha=0.05) +
        # geom_violin(fill="grey") + # Redo the cut_width so the factor labels are nce
        theme_bw() +
        xlab("Embryonic stage") +
        ylab("EM vaf") +
        ggtitle("EM VAF vs. embyronic stage") +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=12),
              strip.text=element_text(size=12),
              title=element_text(size=12),
              plot.title = element_text(hjust = 0.5)
        ) +
        annotate("text",
                 x=n_x,
                 y=max_y,
                 label=sprintf("%.0f", n_per_bin),
                 size=4) +
    expand_limits(y=max_y)
    
    if (log_scale) {
        p = p +
            scale_y_log10() +
            annotation_logticks(sides="l")
    }
    
    if (save_file) {
        p2 = timestamp_ggfigure(p)
        if (is.null(prefix)) {
            plot_fn = sprintf("em_vaf_vs_embyronic_stage.%s.pdf", suffix_sans_extension)    
        } else {
            plot_fn = sprintf("%s.em_vaf_vs_embyronic_stage.%s.pdf", prefix,  suffix_sans_extension)
        }
        
        save_fig(p2, plot_fn, width=fig_width, height=fig_height)
    } else {
        grid.newpage()
        p2 = timestamp_ggfigure(p)
        grid.draw(p2)
    }
}

plot_recurrence_vs_time <- function(master_vaf_df, var_compressed_mrca, tissue_g, save_file=FALSE, prefix=NULL, suffix_sans_extension="", fig_width=6, fig_height=4) {
    # Add the vertex name to the mrca
    var_compressed_mrca$mrca_name = names(V(tissue_g)[var_compressed_mrca$mrca])
    # Distance from zygote
    root_index = get_root(tissue_g, "Zygote")
    V(tissue_g)$dist_from_root = bfs(tissue_g, root_index, dist=TRUE)$dist
    var_compressed_mrca$dist_from_root = V(tissue_g)[var_compressed_mrca$mrca]$dist_from_root
    
    master_vaf_df$row_names = rownames(master_vaf_df)
    
    merge_df = merge(master_vaf_df, var_compressed_mrca, by.x="bin_str", by.y="var_compressed", all.x=TRUE)
    rownames(merge_df) = merge_df$row_names # Note: CANNOT use the same rownames as master_vaf_df because merge() does NOT preserve order
    merge_df$row_names = NULL
    
    max_x = max(merge_df$dist_from_root, na.rm = TRUE) # There can be NAs if all TP variants were blacklisted
    n_x = seq(1,max_x)
    merge_df$my_x_bins = factor(merge_df$dist_from_root, levels=n_x)
 
    n_per_bin_table = table(merge_df$my_x_bins)
    n_per_bin = unname(n_per_bin_table)
    # Orginally, DC had requested to go from newest to oldest, but I'm not sure how to reverse the factor levels that well along with the annotation labels.
    
    # Use bbcall VAF
    p = ggplot(merge_df, aes(my_x_bins, median_vaf_bbcall)) +
        geom_boxplot(fill="grey") + # Redo the cut_width so the factor labels are nce
        theme_bw() +
        xlab("Distance from root") +
        ylab("Median bbcall vaf") +
        ggtitle("Median bbcall VAF vs. distance from tree root") +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=12),
              strip.text=element_text(size=12),
              title=element_text(size=12),
              plot.title = element_text(hjust = 0.5)
        ) +
        scale_y_log10() +
        annotation_logticks(sides="l") +
        expand_limits(y=1) +
        annotate("text",
                 x=n_x,
                 y=1,
                 label=sprintf("%.0f", n_per_bin),
                 size=4)
    
    if (save_file) {
        p2 = timestamp_ggfigure(p)
        if (is.null(prefix)) {
            plot_fn = sprintf("bbcall_vaf_vs_dist_from_root.%s.pdf", suffix_sans_extension)    
        } else {
            plot_fn = sprintf("%s.bbcall_vaf_vs_dist_from_root.%s.pdf", prefix,  suffix_sans_extension)
        }
        
        save_fig(p2, plot_fn, width=fig_width, height=fig_height)
    } else {
        grid.newpage()
        p2 = timestamp_ggfigure(p)
        grid.draw(p2)
    }
    
    # Use em VAF
    plot_recurrence_vs_time_plot_cmd(merge_df, save_file, prefix, suffix_sans_extension, fig_width, fig_height)
    return(merge_df)
    
}

plot_em_vaf_vs_bbcall_vaf <- function(master_vaf_df, save_file=FALSE, prefix=NULL, suffix_sans_extension="", fig_width=6, fig_height=4) {
   
    p = ggplot(master_vaf_df, aes(x=median_vaf_bbcall, y=median_vaf_em)) +
        geom_point(alpha=0.7) +
        xlab("Median bbcall vaf") +
        ylab("EM vaf") +
        ggtitle("EM VAF vs. median bbcall VAF") +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=12),
              strip.text=element_text(size=12),
              title=element_text(size=12),
              plot.title = element_text(hjust = 0.5)
        ) +
        scale_y_log10() +
        scale_x_log10() +
        annotation_logticks(sides="lb") +
        geom_abline(slope=1, intercept = 0, color="red", linetype=3)
    
    if (save_file) {
        p2 = timestamp_ggfigure(p)
        if (is.null(prefix)) {
            plot_fn = sprintf("em_vaf_vs_bbcall_vaf.%s.pdf", suffix_sans_extension)    
        } else {
            plot_fn = sprintf("%s.em_vaf_vs_bbcall_vaf.%s.pdf", prefix,  suffix_sans_extension)
        }
        
        save_fig(p2, plot_fn, width=fig_width, height=fig_height)
    } else {
        grid.newpage()
        p2 = timestamp_ggfigure(p)
        grid.draw(p2)
    }
}

combine_o_m_p_trees <- function(g_observed, missing_tree_weights, possbile_tree_weights, num_edges) {
    info_df = as.data.frame(rbind(E(g_observed)$num_var, missing_tree_weights, possbile_tree_weights)) # O M P = observed, missing, possible
    colnames(info_df) = 1:num_edges
    info_df$type = c("observed", "missing", "possible")
    info_df$type = factor(info_df$type, levels=info_df$type)
    
    # Melt
    info_df_melt = melt(info_df, id.vars=c("type"))
    colnames(info_df_melt) = c("type", "edge", "weight")
    info_df_melt$edge = factor(info_df_melt$edge)
    
    # Transpose
    info_df_t = info_df[,1:num_edges]
    rownames(info_df_t) = info_df$type
    info_df_t = as.data.frame(t(info_df_t))
    
    return(list(info_df=info_df, info_df_melt=info_df_melt, info_df_t=info_df_t))
}

create_observed_tree <- function(var_compressed_counts, var_compressed_mut_vectors, num_var_list, toi, tissue_g, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension, vertex_size, plot_individual_var_tree=FALSE, save_file=FALSE) {
    donor_id = "dummy"
    
    # Note: this wrapper will compress variants
    tmp = mk_donor_g_from_var_df_wrapper(var_compressed_counts, var_compressed_mut_vectors, num_var_list, donor_id, tissue_g, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension=suffix_sans_extension, plot_total_num_var_tree=FALSE, plot_num_var_type_tree=FALSE, plot_individual_var_tree=plot_individual_var_tree, save_file=save_file, vertex_size=vertex_size, toi=toi, print_var_info=TRUE)
    g_observed = tmp$g_consensus
    var_compressed_mrca = tmp$var_compressed_mrca
    
    # Normalize so that the total weight is 1
    for (e_attr in names(edge_attr(g_observed))) {
        value = get.edge.attribute(g_observed, e_attr)
        sum_value = sum(value)
        if (sum_value != 0) {
            norm_value = value / sum_value
        } else { # Otherwise, just keep the weights as is (all 0s)
            norm_value = value
        }
        g_observed = set.edge.attribute(g_observed, e_attr, value=norm_value)
    }
    
    return(list(g_observed=g_observed, var_compressed_mrca=var_compressed_mrca))
}

enumerate_all_missing_and_possible_weights <- function(var_df_list, num_edges, toi, tissue_g, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension) {
    
    num_possible_profile_vectors = length(var_df_list) # FIXME IS THIS RIGHT?
    donor_id = "dummy"
    
    missing_mutation_case_df = as.data.frame(matrix(NA, nrow=num_possible_profile_vectors, ncol=num_edges))
    accessible_mutation_case_df = as.data.frame(matrix(NA, nrow=num_possible_profile_vectors, ncol=num_edges))
    
    for (i in 1:num_possible_profile_vectors) {
        
        profiling_case = var_df_list[[i]]
        
        missing_mutation_case_weights = rep(0, num_edges)
        accessible_mutation_case_weights = rep(0, num_edges)
        
        for (j in 1:nrow(profiling_case)) {
            mutation_case = profiling_case[j,]
            num_var_list = colSums(var_df2compressed_var_df(mutation_case)$var_compressed_counts)
            mutation_g = mk_donor_g_from_var_df_wrapper(mutation_case, num_var_list, donor_id, tissue_g, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension=suffix_sans_extension, plot_total_num_var_tree=FALSE, plot_num_var_type_tree=FALSE, plot_individual_var_tree=FALSE, save_file=FALSE, vertex_size=25, toi=toi)
            
            # Missing tree
            inaccessible_mut_indexes = which(E(mutation_g)$num_var == 0)
            
            inaccessible_mut_weights = rep(0, num_edges)
            inaccessible_mut_weights[inaccessible_mut_indexes] = 1 / length(inaccessible_mut_indexes)
            
            missing_mutation_case_weights = missing_mutation_case_weights + inaccessible_mut_weights
            
            # Accessible tree
            accessible_mutation_case_weights = accessible_mutation_case_weights + E(mutation_g)$num_var
            
        }
        
        # Normalize case weights
        missing_mutation_case_df[i,] = missing_mutation_case_weights / sum(missing_mutation_case_weights) # Use graphs or just an array?  I think an array should be sufficient
        
        accessible_mutation_case_df[i,] = accessible_mutation_case_weights / sum(accessible_mutation_case_weights)
    }
    
    return(list(missing_mutation_case_df=missing_mutation_case_df, accessible_mutation_case_df=accessible_mutation_case_df))
}


create_missing_and_possible_trees <- function(donor_var_df, var_df_list, num_edges, toi, num_tissues, tissue_g, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension) {
    
    # First, calculate the missing and possible weights for all possibe mutations
    tmp = enumerate_all_missing_and_possible_weights(var_df_list, num_edges, toi, tissue_g, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension)
    missing_mutation_case_df = tmp$missing_mutation_case_df
    accessible_mutation_case_df = tmp$accessible_mutation_case_df
    
    # Next, calculate the missing and possible weights for just the donor mutations found
    donor_num_var = nrow(donor_var_df)
    missing_tree_weights = rep(0, num_edges)
    possbile_tree_weights = rep(0, num_edges)
    
    # FIXME: to increase efficiency, compress the variants and just add K*weights

    
    for (k in 1:donor_num_var) {
        var = unlist(donor_var_df[k,1:num_tissues]) # Unlist to avoid removing nownames which would mess with the comparisons below
        # Look up which mutation case this is
        index = which(apply(var_all_possible_df, 1, function(x) {identical(unlist(x[1:num_tissues]), var)}))
        
        if (length(index) != 1) {
            print(var)
            stop("couldn't find current mutation in list of all possible mutations")
        }
        case = var_all_possible_df$case[index]
        
        missing_tree_weights = missing_tree_weights +  missing_mutation_case_df[case,]
        possbile_tree_weights = possbile_tree_weights +  accessible_mutation_case_df[case,]
    }
    
    # Normalize
    missing_tree_weights = missing_tree_weights / sum(missing_tree_weights)
    possbile_tree_weights = possbile_tree_weights / sum(possbile_tree_weights)
    
    return(list(missing_tree_weights=missing_tree_weights, possbile_tree_weights=possbile_tree_weights))   
}

create_missing_and_possible_trees <- function(donor_var_df, var_df_list, num_edges, toi, num_tissues, tissue_g, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension) {
    
    # First, calculate the missing and possible weights for all possibe mutations
    tmp = enumerate_all_missing_and_possible_weights(var_df_list, num_edges, toi, tissue_g, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension)
    missing_mutation_case_df = tmp$missing_mutation_case_df
    accessible_mutation_case_df = tmp$accessible_mutation_case_df
    
    # Next, calculate the missing and possible weights for just the donor mutations found
    donor_num_var = nrow(donor_var_df)
    missing_tree_weights = rep(0, num_edges)
    possbile_tree_weights = rep(0, num_edges)
    
    for (k in 1:donor_num_var) {
        var = unlist(donor_var_df[k,1:num_tissues]) # Unlist to avoid removing nownames which would mess with the comparisons below
        var_freq = donor_var_df$count[k]
        # Look up which mutation case this is
        index = which(apply(var_all_possible_df, 1, function(x) {identical(unlist(x[1:num_tissues]), var)}))
        
        if (length(index) != 1) {
            print(var)
            stop("couldn't find current mutation in list of all possible mutations")
        }
        case = var_all_possible_df$case[index]
        
        missing_tree_weights = missing_tree_weights +  missing_mutation_case_df[case,]*var_freq
        possbile_tree_weights = possbile_tree_weights +  accessible_mutation_case_df[case,]*var_freq
    }
    
    # Normalize
    missing_tree_weights = missing_tree_weights / sum(missing_tree_weights)
    possbile_tree_weights = possbile_tree_weights / sum(possbile_tree_weights)
    
    return(list(missing_tree_weights=missing_tree_weights, possbile_tree_weights=possbile_tree_weights))   
}


enumerate_all_profiling_vectors <- function(num_tissues, toi) {
    
    possible_profile_vectors = expand.grid(rep(list(0:1), num_tissues)) # Copied from https://stackoverflow.com/questions/18705153/generate-list-of-all-possible-combinations-of-elements-of-vector
    colnames(possible_profile_vectors) = toi
    # Remove the all = 0 
    possible_profile_vectors = possible_profile_vectors[rowSums(possible_profile_vectors) != 0,]
    rownames(possible_profile_vectors) = 1:nrow(possible_profile_vectors)
    
    return(possible_profile_vectors)
}

visualize_profiling_vectors <- function(possible_profile_vectors, tissue_g, toi) {
    
    for (i in 1:nrow(possible_profile_vectors)) {
        profiled_tissues = toi[possible_profile_vectors[i,] == 1]
        num_profiled_tissues = length(profiled_tissues)
        profiled_tissue_indexes = get_profiled_tissue_indexes(tissue_g, profiled_tissues)
        my_title = sprintf("Profile vector: case %.0f\n(%.0f tissues profiled)", i, num_profiled_tissues)
        plot_g(tissue_g, type="tissue_tree", category="categorical", title=my_title, save_file=FALSE, vertex_size=25, profiled_tissue_indexes=profiled_tissue_indexes)
    }
}


enumerate_all_mutation_vectors <- function(possible_profile_vectors, tissue_g, toi, num_tissues, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension, plot_individual_var_tree=FALSE) {
    
    # Make all possible mutation vectors & count the # of variants
    # TODO: Could probably count the # of mutations with an equation
    num_possible_profile_vectors = nrow(possible_profile_vectors)
    var_df_list = vector("list", num_possible_profile_vectors)
    num_var = 0
    for (case_i in 1:num_possible_profile_vectors) {
        var_df = mk_complete_mutation_vector(case_i, toi, possible_profile_vectors, recurrent_only=TRUE)
        var_df_list[[case_i]] = var_df

        if (! is.null(var_df)) {
            num_var = num_var + nrow(var_df)
        }
    }
    
    # For historical purposes, I think the mutation vector generation is based on which tissues were sampled, thus the # of columns names can change.  Convert the list to a df.  Use NA if not profiled.
    # TODO: Might want to just make mk_complete_mutation_vector output a df
    var_all_possible_df = as.data.frame(matrix(NA, nrow=num_var, ncol=(num_tissues+1)))
    colnames(var_all_possible_df) = c(toi, "case")
    
    current_row = 1
    for (case_i in 1:num_possible_profile_vectors) {
        
        var_df = var_df_list[[case_i]]
        
        if (! (is.null(var_df))) {
            num_var_in_case = nrow(var_df)
            final_row = current_row + num_var_in_case - 1
            
            # Update the df of all possible variants
            # Find the column indexes that weren't profiled
            profiled_tissue_indexes = which(possible_profile_vectors[case_i,] != 0)
            
            var_all_possible_df[current_row:final_row, profiled_tissue_indexes] = var_df[,1:(ncol(var_df) - 1)] # -1 since var_df contains ref column
            # Update the case column
            var_all_possible_df[current_row:final_row, "case"] = case_i
            
            # Update
            current_row = final_row + 1
        }
    }
    
    # Create the graph for each mutation vector
    # Note: since I also want to recover the weights for each possible mutation vector, and mk_donor_g_from_var_df_wrapper returns the consensus, run function on each variant individually.  Knowing the weights for each possible mutation will be helpful for unit testing
    edge_names = as.character(seq(1:length(E(tissue_g))))
    var_all_possible_df[,edge_names] = NA
    
    #g_possible_list = vector("list", num_possible_profile_vectors)
    row_counter = 1
    for (case_i in 1:num_possible_profile_vectors) {
        var_df = var_df_list[[case_i]]
        
        if (! is.null(var_df)) {
            donor_id = sprintf("Case %.0f", case_i)
            
            for (var_j in 1:nrow(var_df)) {
                v = var_df[var_j,]
                g = mk_donor_g_from_var_df_wrapper(v, var_df_list, donor_id, tissue_g, tissue_abbrev_df, vartype_colors_palette, base_colors_df, core_vartypes, reserved_tissues, suffix_sans_extension=suffix_sans_extension, plot_total_num_var_tree=FALSE, plot_num_var_type_tree=FALSE, plot_individual_var_tree=plot_individual_var_tree, save_file=FALSE, vertex_size=25, toi=toi, var_type_index=var_j)
                #g_possible_list[[case_i]] = g
            
                edge_weights = edge_attr(g)$num_var
                var_all_possible_df[row_counter,edge_names] = edge_weights
                row_counter = row_counter + 1
            }
        }
    }
    
    return(list(var_df_list=var_df_list, var_all_possible_df=var_all_possible_df))
}


var_df2compressed_var_df <- function(df, vaf_df) {
    sample_names = colnames(df) # includes ref
    
    bin_df = df
    bin_df[] = "C"
    
    # Set the NA positions
    bin_df[is.na(df)] = NA

    bin_df[df != df[,"ref"]] = "T" # Inspired by https://stackoverflow.com/questions/26724139/by-row-replace-values-equal-to-value-in-specified-column
    
    colnames(bin_df) = paste("bin", colnames(bin_df), sep="_")
    bin_df$bin_str = do.call(paste, bin_df[])
    
    # Create alt column
    if (ncol(df) == 2) { # 1 tissue and 1 ref
        df$alt = tissue_alleles = df[,1]
    } else {
        # Remove an NA alleles and remove ref
        df$alt = apply(df, 1, function(x) {
	    y = unique(x[! is.na(x)])
            y[y != x["ref"]]
        })
    }
    
    df$vartype_str = do.call(paste, c(df[,c("ref", "alt")], sep=">"))
    
    # Combine (row order is the same)
    df_combined = cbind(df, bin_df)
    vaf_df_combined = cbind(vaf_df, bin_str=bin_df$bin_str, my_vartype_str=df$vartype_str)
    
    # Don't want to break earlier stuff; not sure what vaf_df_combined, but now need it as a df
    vaf_df_combined_real_df = data.frame(vaf_df, bin_str=bin_df$bin_str, my_vartype_str=df$vartype_str)
    
    var_compressed = as.data.frame.matrix(table(df_combined$bin_str, df_combined$vartype_str)) # Important: need as.data.frame.matrix; as.data.frame reports a melted df
    
    # Compress on the variant types 
    uncomp2comp = generate_compressed_var()$uncompressed2compressed
    uncomp_var_used = colnames(var_compressed)
    comp_var_used = unname(unlist(uncomp2comp[uncomp_var_used]))
    var_compressed_max = var_compressed
    colnames(var_compressed_max) = comp_var_used
    # Combine by columns
    if (ncol(var_compressed_max) > 1) {
        var_compressed_max = data.frame(do.call(cbind, by(t(var_compressed_max),INDICES=names(var_compressed_max),FUN=colSums)), stringsAsFactors = FALSE, check.names=FALSE) # Inspired by https://stackoverflow.com/questions/8961063/combining-multiple-identically-named-columns-in-r
    }
    
    # Add on the total # of variants
    var_compressed$count = rowSums(var_compressed)
    var_compressed_max$count = rowSums(var_compressed_max)
    
    # Reinflate the binary vector
    bin_df_compressed = as.data.frame(do.call(rbind, strsplit(rownames(var_compressed), " ")), stringsAsFactors = FALSE, check.names=FALSE)
    # Convert "NA" back to NA
    bin_df_compressed[bin_df_compressed == "NA"] = NA
    colnames(bin_df_compressed) = sample_names
    
    #var_compressed = cbind(var_compressed, bin_df_compressed)
    #var_compressed_max = cbind(bin_df_compressed, var_compressed_max)
    
    return(list(var_compressed_counts=var_compressed_max, var_compressed_mut_vectors=bin_df_compressed, vaf_df_combined=vaf_df_combined, vaf_df_combined_real_df=vaf_df_combined_real_df))
}

master_var2mega_df <- function(master_var_df, remove_ref=FALSE, set_na_to_ref=FALSE, remove_pan_tissue_var=TRUE) {
    mega_df = master_var_df
    mega_missing_val = "-"
    
    # MOve ref to the first species
    # This is helpful 1) to access the ref col/row and 2) in MEGA, the first species is used as the reference in the dot notation
    cols = 1:ncol(mega_df)
    ref_index = which(colnames(mega_df) == "ref")
    nonref_indexes = rows[which(cols != ref_index)]
    new_col_order = c(ref_index, nonref_indexes)
    
    mega_df = mega_df[,new_col_order]
    
    # Remove sites that are no longer recurrent (this can happen if master_var_df was subset to a list of toi)
    roi_recurrent = apply(mega_df, 1, function(x) {
        if (sum(x != x[["ref"]], na.rm=TRUE) >= 2) { # Since multialleic sites are not possible, finding at leasdt 2 differences means that the same variant, e.g., T>C is found in at leasdt 2 tissues
            return(TRUE)
        } else{
            return(FALSE)
        }
    })
    
    mega_df = mega_df[roi_recurrent,]
    
    # Remove pan-tissue variants (i think this might be only helpful if including the reference)
    if (remove_pan_tissue_var) {
        num_tissues = ncol
        roi_not_pan_tissue = apply(mega_df, 1, function(x) {
            if (sum(x == x[["ref"]], na.rm=TRUE) >= 2) { # Find the ref at least 2X (once in the ref, and once in another tissue; this means it is not a pan tissue site)
                return(TRUE)
            } else{
                return(FALSE)
            }
        })
        
        mega_df = mega_df[roi_not_pan_tissue,]
    }
    
    # This transposed the data?? FIXME tissue x site
    # Assumption: any NA site is ref
    # Rename NA sites to the ref
    if (set_na_to_ref) {
        mega_df = as.data.frame(apply(mega_df, 1, function(x) {x[which(is.na(x))] = x[["ref"]]; return(x)}))
    } else {
        mega_df = as.data.frame(apply(mega_df, 1, function(x) {x[which(is.na(x))] = mega_missing_val; return(x)}))
    }
    
    # Remove ref row
    if (remove_ref) {
        mega_df = mega_df[which(rownames(mega_df) != "ref"),]
    }
    
    return(mega_df)
}




load_time_dim_data <- function(species) {
    if (species == "human") {
        embryonic_staging_fn = "/home/nicole/projects/sommut/sommut/data/v8/embryonic_staging_of_tissue_development.human.20200316.txt"
        embryonic_time_df = read.table(embryonic_staging_fn, header=TRUE, stringsAsFactors = FALSE, sep="\t", quote="")
        
        # Rename Carnegie stage to generic stage
        colnames(embryonic_time_df)[colnames(embryonic_time_df) == "carnegie_stage_start"] = "stage_start"
    } else if (species == "mouse") {
        embryonic_staging_fn = "/home/nicole/projects/sommut/sommut/data/embryonic_staging_of_tissue_development.tsv"
        embryonic_staging_df = read.table(embryonic_staging_fn, header=TRUE, stringsAsFactors = FALSE, sep="\t", quote="")
        
        embryonic_staging2dpc_fn = "/home/nicole/projects/sommut/sommut/data/embryonic_staging2dpc.txt"
        embryonic_staging2dpc_df = read.table(embryonic_staging2dpc_fn, header=TRUE, stringsAsFactors = FALSE, sep="\t", quote="")
        
        embryonic_time_df = merge(embryonic_staging_df, embryonic_staging2dpc_df, by.x="ts_start", by.y="ts_stage")
        
        # Rename Theiller stage to generic stage
        colnames(embryonic_staging_df)[colnames(embryonic_staging_df) == "ts_stage"] = "stage_start"
    } else {
        writeLines(sprintf("ERROR: unknown species %s", species), con=stderr())
        stop()
    }
    
    return(embryonic_time_df)
}

add_time_dim_to_vertexes <- function(donor_g, embryonic_time_df) {
    # Add the timing to the vertex information
    # Subset the vertex info to just the voi; reorder the timing data to match the order of the vertexes
    voi = names(V(donor_g))
    tmp = embryonic_time_df[embryonic_time_df$tissue %in% voi,]
    tmp = tmp[match(voi, tmp$tissue),]
    
    donor_g = set.vertex.attribute(donor_g, "stage_start", value=tmp$stage_start)
    donor_g = set.vertex.attribute(donor_g, "dpc", value=tmp$dpc)
    
    # FIXME
    # Currently, don't have timing for tissues.  Set to +1 max
    max_time = max(V(donor_g)$dpc, na.rm=TRUE)
    V(donor_g)$dpc[is.na(V(donor_g)$dpc)] = max_time + 1
    
    return(donor_g)
}

plot_g_w_time_dim <- function(donor_g, vertex_label_cex=0.8, edge_label_cex=0.8, vertex_size = 5) {
    # Original layout (bottom to top)
    vertex_coordinates = layout_as_tree(donor_g, flip.y = FALSE)
    # Adjust the y-values to map onto the new timing axis
    vertex_coordinates[,2] = V(donor_g)$dpc # Using DPC
    
    # Create new labels for y-axis
    range(vertex_coordinates[,2], na.rm = TRUE)
    pretty_breaks = pretty(vertex_coordinates[,2], n=4)
    pretty_range = range(pretty_breaks)
    pretty_min = min(pretty_breaks)
    pretty_max = max(pretty_breaks)
    
    # Normaize the y-axis myself since we need to use the pretty_range (which might be a little larger than the range of just the data)
    plot_min = -1
    plot_max = 1
    slope = (plot_max - plot_min)/(pretty_max - pretty_min)
    vertex_coordinates[,2] = slope * vertex_coordinates[,2] + plot_min
    
    pretty_breaks_transformed = slope * pretty_breaks + plot_min
    
    # Add the edge and vertex attributes
    g = plot_g(donor_g, type="mutation_prob", category="numeric", save_file=FALSE, vertex_size=vertex_size, set_vertex_labels=TRUE, set_edge_labels=TRUE, value_global_max=1, set_na_colors=TRUE, mk_plot=FALSE)

    plot.igraph(g, 
                layout=vertex_coordinates, 
                
                vertex.color=V(g)$vertex_fill_color, 
                vertex.size=vertex_size,
                vertex.label=V(g)$abbrev,
                vertex.label.cex=vertex_label_cex, 
                vertex.frame.color=V(g)$vertex_border_color,
                vertex.shape=V(g)$vertex_shape,
                vertex.label.family="sans",
                vertex.label.font=1,
                vertex.label.color = V(g)$label_color,
                edge.color=E(g)$edge_color,
                edge.width=E(g)$edge_width,
                edge.arrow.mode=0, # Turn off arrows since edge widths don't always play well with arrow size
                edge.label=E(g)$edge_label,
                edge.label.color="red",
                edge.label.family="sans",
                edge.label.font=1,
                edge.label.cex=edge_label_cex,
                
                axes=FALSE)
    axis(2, at=pretty_breaks_transformed, labels=pretty_breaks)
    timestamp_figure(cex=0.7)
}

shuffle_mutation_type <- function(var_types) {
    if (length(var_types) > 1) {
        var_types_shuffled = sample(var_types)
    } else {
        var_types_shuffled = var_types
    }
    return(var_types_shuffled)
}

shuffle_mutation_vector <- function(x) {
    # Only shuffle the ref and alt; keep NAs as NAs; keep the ref
    #x = as.character(x)
    fixed_indexes = is.na(x)
    x_unfixed = x[! fixed_indexes]
    if (length(x_unfixed) > 1) {
        x_unfixed_shuffled = sample(x_unfixed)
    } else {
        x_unfixed_shuffled = x_unfixed
    }
    
    # Put humpty dumpty back together again
    x_shuffled = x
    x_shuffled[! fixed_indexes] = x_unfixed_shuffled
    return(x_shuffled)
}

g2root_distances <- function(g, root_name="Zygote") {
    distance_from_root = data.frame(bfs(g, root_name, dist=TRUE)$dist)
    colnames(distance_from_root) = c("distance")
    distance_from_root$edge = rownames(distance_from_root)
    rownames(distance_from_root)= NULL
    # Sort on distance, then on edge (alphabetical)
    distance_from_root = distance_from_root[with(distance_from_root, order(distance, edge)),]
    
    return(distance_from_root)
}
