
load_mutation_type_data <- function(mutation_type_dist_fn, mutation_vector_str2edge, edge_indexes_to_exclude, min_raw_mutations_per_edge=5) {
    
    mutation_type_dist = readRDS(mutation_type_dist_fn)
    
    # Remove the edges that have low power
    # Well, for now this is a moot point since no variants were observed on these edges anyways.
    mutation_type_dist = mutation_type_dist[! mutation_type_dist$edge %in% edge_indexes_to_exclude,]
    
    mutation_type_dist = merge(mutation_type_dist, mutation_vector_str2edge[,c("head_node", "tail_node", "pretty_edge", "edge_index")], by.x="edge", by.y="head_node")
    
    mutation_type_dist_filtered = mutation_type_dist[mutation_type_dist$total_count >= min_raw_mutations_per_edge,]
    nrow_before = nrow(mutation_type_dist)
    nrow_after = nrow(mutation_type_dist_filtered)
    nrow_removed = nrow_before - nrow_after
    
    writeLines(sprintf("Requiring at least %d mutations on an edge removed %d/%d (%.1f%%) edges", min_raw_mutations_per_edge, nrow_removed, nrow_before, nrow_removed/nrow_before *100))
    
    mutation_type_dist_melt = reshape2::melt(mutation_type_dist_filtered, measure.vars=compressed_var)
    
    return(list(mutation_type_dist=mutation_type_dist, mutation_type_dist_filtered=mutation_type_dist_filtered, mutation_type_dist_melt=mutation_type_dist_melt))
}


load_data <- function(tree_name, consensus_g_fn, local, edge_indexes_to_exclude) {

    tmp = load_tissue_tree(local, name=tree_name)
    tissue_g = tmp$tissue_g
    num_edges = length(E(tissue_g))
    
    leaves = get_leaves(tissue_g)
    leaf_edge_indexes = as.numeric(incident_edges(tissue_g, leaves$leaf_indexes, mode="in"))
    non_leaf_edge_indexes = setdiff(E(tissue_g), leaf_edge_indexes)
    
    mutation_vector_str2edge = mk_descendent_trees(tissue_g, return_lookup=TRUE)$mutation_vector_str2edge
    mutation_vector_str2edge$edge_index = as.numeric(mutation_vector_str2edge$edge)
    mutation_vector_str2edge$edge = NULL
    
    # Remove edges with low power (see compare_obs_vs_expected_trees.Rmd for how indexes were determined)
    
    edges_to_exclude = mutation_vector_str2edge$head_node[mutation_vector_str2edge$edge %in% edge_indexes_to_exclude] # Yes, use head_of to get the *incident* node, e.g., head of A->B will return B
    mutation_vector_str2edge = mutation_vector_str2edge[! mutation_vector_str2edge$head_node %in% edges_to_exclude,]
    
    consensus_g = readRDS(consensus_g_fn)
    
    E(consensus_g)$edge_width = 3
    E(consensus_g)$edge_width[non_leaf_edge_indexes] = 25
    
    # Set weight to NA for edges with lower power
    E(consensus_g)$weight_unknown[edge_indexes_to_exclude] = TRUE
    
    # Set germlayer nodes to gray (like the rest of the internal nodes)
    germlayer_nodes = which(names(V(consensus_g)) %in% c("Ectoderm", "Mesoderm", "Endoderm"))
    V(consensus_g)$tissue_color[germlayer_nodes] = "#D3D3D3"
    
    # Calculate the distance to root for each vertex
    distance_from_root = g2root_distances(tissue_g)
    
    return(list(tissue_g=tissue_g,
                consensus_g=consensus_g,
                leaves=leaves,
                non_leaf_edge_indexes=non_leaf_edge_indexes,
                distance_from_root=distance_from_root,
                mutation_vector_str2edge=mutation_vector_str2edge))
}

mk_color_bar <- function(g, fn, log_trans=FALSE, percent=FALSE, save_file=TRUE, use_custom_range=FALSE, my_min=NULL, my_max=NULL, ntick_labels = 2) {
    if (use_custom_range) {
        if (is.null(my_min) | is.null(my_max)) {
            stop("ERROR: my_min and my_max must be specified")
        }
    } else {
        my_max = max(E(g)$values_used, na.rm=TRUE)
        my_min = min(E(g)$values_used, na.rm=TRUE)
    }
    
    
    f <- colorRamp(viridis(100))
    
    colors <- rgb(f(seq(0,1,by=0.005))/255)
    lut = colors
    scale = (length(lut)-1)/(my_max-my_min)
    
    if (save_file) {
        pdf(fn, width=2, height=7)
    }
    #par(bg=NA)
    par(las=1) # Execute for vertical color ramp
    plot(c(0,10), c(my_min,my_max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
    for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + my_min
        rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
    
    yticks_val <- pretty_breaks(n=ntick_labels)(seq(my_min, my_max, length.out=ntick_labels))
    
    if (log_trans) {
        if (percent) {
            yticks_lab = sprintf("%.1f%%", 10^(yticks_val)*100)
        } else {
            yticks_lab = sprintf("%.1f", 10^(yticks_val))
        }
    } else {
        if (percent) {
            yticks_lab = sprintf("%.1f%%", yticks_val*100)
        } else {
            yticks_lab = sprintf("%.1f", yticks_val)
        }
    }
    axis(2, at=yticks_val, lab=yticks_lab, cex.axis=1)
    
    timestamp_figure()
    
    if (save_file) {
        dev.off()
        write_new_file_message(fn)
    }
}


compare_mutation_type_distribution <- function(name_1, name_2, data, compressed_var, test_name) {
    
    if (name_1 %in% data$edge & name_2 %in% data$edge) {
        
        data_1 = as.numeric(data[data$edge == name_1, compressed_var])
        data_2 = as.numeric(data[data$edge == name_2, compressed_var])
        
        # If both data 1 and data 2 have 0 counts for factor, remove it.  It will cause div by 0 issus and could argue that if category not observed nor expected, then it's not informative so remove it
        vartype_to_remove = data_1 == 0 & data_2 == 0
        data_1 = data_1[! vartype_to_remove]
        data_2 = data_2[! vartype_to_remove]
        
        # Expected counts of 0 are not fun
        # Maybe 1) drop that factor levels from data1, calculate p-value, 2) repeat with data2, and take the min p-value or 3) drop from both?
        # Maybe combined C* and T*.  Is there any better way to group the data?
        if (any(data_1 == 0) & any(data_2 == 0)) {
            print("Zeros in both groups")
            # Drop the 0's from data_1
            data_1_sanitized = data_1[data_1 != 0]
            data_2_sanitized = data_2[data_1 != 0] # Drop the same data_1 levels
            tmp1 = perform_test(data_2_sanitized, data_1_sanitized, test_name, pick_expected_group=FALSE)
            
            # Drop the 0's from data_2
            data_1_sanitized = data_1[data_2 != 0] # Drop the same data_2 levels
            data_2_sanitized = data_2[data_2 != 0] 
            tmp2 = perform_test(data_1_sanitized, data_2_sanitized, test_name, pick_expected_group=FALSE)
            
            # Report the minimum p-value
            # min p-value: would throw out data where Obs=0, exp=22; min p-value would keep it?  FIXME
            if (tmp1[1] < tmp2[1]) {
                p_val = tmp1[1]
                p_val.chisq = tmp1[2]
            } else {
                p_val = tmp2[1]
                p_val.chisq = tmp2[2]
            }
            
        } else {
            
            tmp = perform_test(data_1, data_2, test_name)
            p_val = tmp[1]
            p_val.chisq = tmp[2]
        }
    } else {
        p_val = NA
        p_val.chisq = NA
    }
    return(c(p_val, p_val.chisq))
}


perform_test <- function(data_1, data_2, test_name, pick_expected_group=TRUE) {
    # Use the class with more observations  as the expected (may be estimated better, fewer changes of having 0, fewer computations for xmulti)
    if (pick_expected_group) {
        if (sum(data_1) > sum(data_2)) {
            var_1 = data_2
            var_2 = data_1
        } else {
            # Data 2 is used for expected
            var_1 = data_1
            var_2 = data_2
        }
    } else {
        var_1 = data_1
        var_2 = data_2
    }
    
    if (length(var_1) >= 2) { # Need 2+ groups
        if (test_name == 'exact_multinomial') {
            p_val = xmulti(var_1, var_2, detail=0)$pProb
        } else if (test_name == 'mc_multinomial') {  
            # Got an error:
            # This operation could take more than several minutes.
            # The monte carlo version, "xmonte" is recommended for this case. 
            
            # The error message suggests to use xmonte:
            # Use xmonte to compute a P value to test whether a set of counts fits a specific multinomial distribution. It does this by examining a large number of random outcomes and finding the probability of those cases which deviate from the expectation by at least as much as the observed.
            # FIXME why is there no pProb?
            p_val = xmonte(var_1, var_2, detail=0)$pLLR
        }
        
        # Also perform chi-sequare
        #print(var_1)
        #print(var_2)
        p_val.chisq = chisq.test(var_1, p=var_2, rescale.p=TRUE)$p.value
    } else {
        p_val = NA
        p_val.chisq = NA
    }
    
    return(c(p_val, p_val.chisq))
}

compare_c_t_collapsed_mutation_type_distributions_in_tree <- function(compressed_var, mutation_type_dist) {
    c_var = compressed_var[grep("^C>", compressed_var)]
    t_var = compressed_var[grep("^T>", compressed_var)]
    
    mutation_type_dist_compressed = mutation_type_dist
    mutation_type_dist_compressed$c_var = apply(mutation_type_dist_compressed[,c_var], 1, sum)
    mutation_type_dist_compressed$t_var = apply(mutation_type_dist_compressed[,t_var], 1, sum)
    
    x = compare_mutation_type_distributions_in_tree(mutation_type_dist_compressed, mutation_vector_str2edge, "mc_multinomial", leaves$leaf_names, c("c_var", "t_var"))
}

compare_mutation_type_distributions_in_tree <- function(mutation_type_dist, mutation_vector_str2edge, test_name, leaf_names, compressed_var, sig_threshold=0.05) {
    # Comparisons of interest:
    # Parent vs child
    # Child vs child
    
    # Get internal edges
    # For each internal node (parent)
    #    get children
    #    for each child
    #        compare child, parent
    #    for each child combination
    #        compare child i, child j
    
    # Counts rounded to nearest integer
    mutation_type_dist_counts = mutation_type_dist
    mutation_type_dist_counts[,compressed_var] = round(mutation_type_dist_counts[,compressed_var]*mutation_type_dist_counts$total_count)
    
    compare_mut_type_results_df = data.frame(edge1=character(0),
                                             edge2=character(0), 
                                             p_value.multinomial=numeric(0),
                                             p_value.chisq=numeric(0), 
                                             comparison_type=character(0),
                                             comparison_group=character(0),
                                             stringsAsFactors = FALSE
    )
    
    for (parent in mutation_vector_str2edge$head_node) {
        children = names(get_children(tissue_g, parent))
        
        # Remove comparisons involving leaves
        children = setdiff(children, leaf_names)
        # print("THESE ARE THE CHILDREN")
        # print(children)
        if (length(children) >= 2) {
            combinations = t(combn(children, 2))
            
            num_comparisons = length(children) + nrow(combinations)
            tmp = data.frame(edge1=character(num_comparisons),
                             edge2=character(num_comparisons),
                             p_value.multinomial=numeric(num_comparisons),
                             p_value.chisq=numeric(num_comparisons),
                             comparison_type=character(num_comparisons),
                             comparison_group=character(num_comparisons),
                             stringsAsFactors = FALSE)
            i = 1
            
            # Compare child, parent
            # print("Parent vs. child")
            for (child in children) {
                # print(parent)
                # print(child)
                p_vals = compare_mutation_type_distribution(child, parent, mutation_type_dist_counts, compressed_var, test_name)
                
                tmp[i,] = c(parent, child, p_vals, "parent_vs_child", parent)
                i = i+1
            }
            
            # Compare children
            # print("Child vs. child")
            for (k in 1:nrow(combinations)) {
                child_i = combinations[k,1]
                child_j = combinations[k,2]
                # print(child_i)
                # print(child_j)
                
                p_vals = compare_mutation_type_distribution(child_i, child_j, mutation_type_dist_counts, compressed_var, test_name)
                tmp[i,] = c(child_i, child_j, p_vals, "child_vs_child", parent)
                i = i+1
            }
            
            compare_mut_type_results_df = rbind(compare_mut_type_results_df, tmp)
        }
    }
    
    compare_mut_type_results_df$q_value.multinomial = p.adjust(compare_mut_type_results_df$p_value.multinomial, method="fdr")
    
    # Print table of how many comparisons were/n't significant
    writeLines(sprintf("\n%d/%d of edge vs. edge comparisons were significant", sum(compare_mut_type_results_df$q_value.multinomial <= sig_threshold, na.rm=TRUE), nrow(compare_mut_type_results_df)))
    
    # Print list of comparisons that were significant
    writeLines("\nThe following edge vs. edge comparisons were significant")
    write.table(compare_mut_type_results_df[compare_mut_type_results_df$q_value.multinomial < sig_threshold & ! is.na(compare_mut_type_results_df$q_value.multinomial), c("edge1", "edge2", "comparison_type")], sep="|", quote=FALSE, row.names = FALSE)
    
    return(compare_mut_type_results_df)
}

test_if_equal_counts <- function(mutation_type_dist, compressed_var, sig_threshold=0.05, test_name='mc_multinomial') {
    # Counts rounded to nearest integer
    mutation_type_dist_counts = mutation_type_dist
    mutation_type_dist_counts[,compressed_var] = round(mutation_type_dist_counts[,compressed_var]*mutation_type_dist_counts$total_count)
    mutation_type_dist_counts$total_count_round = apply(mutation_type_dist_counts[,compressed_var], 1, sum)
    # Step 1: test if distribution is different than uniform

    mutation_type_dist_counts$p_value.chisq = apply(mutation_type_dist_counts[,compressed_var], 1, function(x){
        chisq.test(x)$p.value # Test if counts have a uniform distribution
    })
    
    expected_p = rep(1/length(compressed_var), length(compressed_var))
    if (test_name == 'exact_multinomial') {
        mutation_type_dist_counts$p_value.multinomial = apply(mutation_type_dist_counts[,compressed_var], 1, function(x){
            xmulti(obs=x, expr=expected_p, detail=0)$pProb
        })
    } else if (test_name == 'mc_multinomial') {
        mutation_type_dist_counts$p_value.multinomial = apply(mutation_type_dist_counts[,compressed_var], 1, function(x){
            xmonte(obs=x, expr=expected_p, detail=0)$pLLR
        })
    }
    
    # Step 2: post-hoc tests: test which mutation types are significant
    # Will test for significance with p-values rather than q-values because want to have low bar to get to next stage.  The next stage takes into account multiple tests so any borderline cases should be taken care of.
    mutation_type_dist_counts$is_sig.multinomial = mutation_type_dist_counts$p_value.multinomial <= sig_threshold
    step2_df = mutation_type_dist_counts[mutation_type_dist_counts$is_sig.multinomial,]
    failed_df = mutation_type_dist_counts[! mutation_type_dist_counts$is_sig.multinomial,]
    
    
    step2_melt_df = reshape2::melt(step2_df, measure.vars = compressed_var)
    
    
    # mapply for the win!
    step2_melt_df$p_value.binomial = mapply(my_binom.test,
                                            x=step2_melt_df$value,
                                            n=step2_melt_df$total_count_round,
                                            p=1/length(compressed_var)
    )
    
    step2_melt_df$q_value.binomial = p.adjust(step2_melt_df$p_value.binomial, method="fdr")
    
    step2_melt_df$is_sig.binomial = step2_melt_df$q_value.binomial <= sig_threshold
    
    # Merge with the edges that had a uniform distribution
    if (nrow(failed_df) > 0) {
        failed_melt_df = reshape2::melt(failed_df, measure.vars = compressed_var)
        failed_melt_df$p_value.binomial = NA
        failed_melt_df$q_value.binomial = NA
        failed_melt_df$is_sig.binomial = NA
        
        mutation_type_dist_counts_melt = rbind(step2_melt_df, failed_melt_df)
    } else {
        mutation_type_dist_counts_melt = step2_melt_df
    }
    
    # Print table of how many comparisons were/n't significant
    writeLines(sprintf("\n%d/%d edges had vartype distribution != uniform", nrow(step2_df), nrow(mutation_type_dist_counts)))
    writeLines(sprintf("\nOf the edges that did have a non-uniform vartype distribution, %d/%d vartypes were different than expected", sum(step2_melt_df$is_sig.binomial), nrow(step2_melt_df)))
    
    # Calculate the breakdown of var types that were signifificantly different from random
    writeLines("\nBreakdown of vartypes that were significant:")
    num_sig_vartypes_df = data.frame(prop.table(table(mutation_type_dist_counts_melt$variable[mutation_type_dist_counts_melt$is_sig.binomial])))
    write.table(num_sig_vartypes_df, sep="|", quote=FALSE, row.names = FALSE)
    
    # Print # of var types that were significantly different from uniform
    writeLines("\nFor edges with at least 1 sig var type, # vartypes that were significant per edge:")
    num_sig_vartypes_df = aggregate(is_sig.binomial ~ edge, mutation_type_dist_counts_melt, sum)
    num_sig_vartypes_df = num_sig_vartypes_df[num_sig_vartypes_df$is_sig.binomial > 0,]
    colnames(num_sig_vartypes_df) = c("edge", "num_sig_vartypes")
    
    write.table(num_sig_vartypes_df, sep="|", quote=FALSE, row.names = FALSE)
    
    return(mutation_type_dist_counts_melt)
}

test_if_predominant_type_exists <- function(mutation_type_dist, compressed_var, sig_threshold=0.05) {
    # Counts rounded to nearest integer
    mutation_type_dist_counts = mutation_type_dist
    mutation_type_dist_counts[,compressed_var] = round(mutation_type_dist_counts[,compressed_var]*mutation_type_dist_counts$total_count)
    mutation_type_dist_counts$total_count_round = apply(mutation_type_dist_counts[,compressed_var], 1, sum)
    
    # Test if the most common mutation type is more common than the next most common mutation type
    # FIXME Do I need to do this with a beta binom test inscase the p parameter is not estimated well?
    mutation_type_dist_counts$p_value = NA
    mutation_type_dist_counts$predominant_type = NA
    mutation_type_dist_counts$at_least_2_var_types = NA
    for (i in 1:nrow(mutation_type_dist_counts)) {
        # Find the max mutation type
        sorted_counts = sort(mutation_type_dist_counts[i, compressed_var], decreasing=TRUE)
        total_counts = mutation_type_dist_counts$total_count_round[i]
        p = as.numeric(sorted_counts[2]/total_counts)
        mutation_type_dist_counts$p_value[i] = binom.test(x=as.numeric(sorted_counts[1]), 
                                                          n=total_counts, 
                                                          p=p,
                                                          alternative=c("greater"))$p.value
        mutation_type_dist_counts$predominant_type[i] = names(sorted_counts[1])
        # New: 20210430 test is only really applicable if have at least 2 var types
        if (p == 0) {
            mutation_type_dist_counts$at_least_2_var_types[i] = FALSE
        } else {
            mutation_type_dist_counts$at_least_2_var_types[i] = TRUE
        }
    }
    
    mutation_type_dist_counts$p_value[mutation_type_dist_counts$at_least_2_var_types == FALSE] = NA
    
    mutation_type_dist_counts$q_value = p.adjust(mutation_type_dist_counts$p_value, method="fdr") # the NA wont contribute to the # of hypothesis tests
    
    mutation_type_dist_counts$is_sig = mutation_type_dist_counts$q_value <= sig_threshold
    
    writeLines(sprintf("%d/%d edges had a significant predominant edge", sum(mutation_type_dist_counts$is_sig, na.rm=TRUE), nrow(mutation_type_dist_counts)))
    
    return(mutation_type_dist_counts)
}

my_binom.test <- function(x,n,p) {
    # Calls binom.test and returns just the p-value (to play nice with mapply)
    results = binom.test(x=x, n=n, p=p)
    p_value = results$p.value
    return(p_value)
}

plot_predominant_mutation_type <- function(predominant_edge_df, g, vartype_palette_df, na_color="gray", non_sig_color="#051E9D") {
    
    # Add edge colors to df
    predominant_edge_df = merge(predominant_edge_df, vartype_palette_df, by.x=c("predominant_type"), by.y=c("X.var_type"))
    
    # Set non-significant edges to have thin widths.  Using a different color made it difficult to draw attention to it (especitally since gray and black are already taken)
    sig_edge_indexes = predominant_edge_df$edge_index[predominant_edge_df$is_sig & ! is.na(predominant_edge_df$is_sig)]
    na_edge_indexes = predominant_edge_df$edge_index[is.na(predominant_edge_df$is_sig)]
    
    E(g)$edge_width = 3
    E(g)$edge_width[sig_edge_indexes] = 25
    
    # Set NA edges (not enough mutations to do a test) to have dashed
    E(g)$lty = rep(1, length(E(g)))
    E(g)$lty[na_edge_indexes] = 5 # Dashed
    
    # Add edge colors to graph
    E(g)$edge_color = na_color
    E(g)$edge_color[sig_edge_indexes] = predominant_edge_df$color_code[predominant_edge_df$is_sig & ! is.na(predominant_edge_df$is_sig)]
    
    
    consensus_g_raw_weights = plot_g(g, type="predominant_type", category="leave_edges_alone", save_file=FALSE, set_na_colors = FALSE, set_widths=FALSE, set_edge_labels=FALSE, set_vertex_labels = FALSE)
}

plot_blank_tissue_tree <- function(tissue_g, mutation_type_dist_filtered, non_leaf_edge_indexes, save_file=TRUE, edge_width=3) {
    E(tissue_g)$edge_width = edge_width
    
    # Update the edges that were NA 
    unknown_edge_indexes = setdiff(non_leaf_edge_indexes, mutation_type_dist_filtered$edge_index)
    E(tissue_g)$edge_width[unknown_edge_indexes] = 25
    
    # Set germlayer nodes to gray (like the rest of the internal nodes)
    germlayer_nodes = which(names(V(tissue_g)) %in% c("Ectoderm", "Mesoderm", "Endoderm"))
    V(tissue_g)$tissue_color[germlayer_nodes] = "#D3D3D3"
    
    fn = file.path(fig_dir, sprintf("mutation_type_distribution.tree.%s.%s.pdf", tree_name, today))
    if (save_file) {
        pdf(fn, height=20, width = 20)
    }
    x = plot_g(tissue_g, "tissue_tree", "leave_edges_alone", set_edge_labels=FALSE, set_vertex_labels=FALSE, set_edge_color=TRUE, set_na_colors=FALSE, save_file = FALSE)
    
    if (save_file) {
        dev.off()
        write_new_file_message(fn)
    }
}


plot_mutation_type_within_edge <- function(compare_mut_type_within_edge_df) {
    facet_pretty_names = c("TRUE"="Vartype distribution is\nNOT uniform", 
                           "FALSE"="Vartype distribution is\nuniform")
    p = ggplot(compare_mut_type_within_edge_df, aes(x=pretty_edge, y=value)) +
        geom_bar(aes(fill=variable, col=is_sig.binomial), position = "fill", stat="identity", size=1) +
        facet_grid(~is_sig.multinomial, labeller = as_labeller(facet_pretty_names)) +
        coord_flip() +
        vartype_palette$fill_scale +
        scale_color_hue(name="Is vartype different\nfrom uniform?") +
        scale_y_continuous(labels=scales::percent) +
        xlab("Edge") +
        ylab("Percent of mutations") +
        my_theme() +
        theme_classic()
    grid.newpage()
    p2 = timestamp_ggfigure(p)
    grid.draw(p2)
}

plot_mutation_type_distribution_comparisons <- function(mutation_type_dist_melt, distance_from_root, tree_name, today, fig_dir, save_file=TRUE, width=10, height=4) {
    
    tmp = mutation_type_dist_melt
    tmp$parent=tmp$tail_node
    tmp$is_parent = "children"
    tmp$is_parent[tmp$parent == tmp$edge] = "parent"
    tmp$is_parent = factor(tmp$is_parent, levels=c("parent", "children"))
    
    # Order edges based on distance from top  Therefore, parents will always be before children
    tmp$edge = factor(tmp$edge, levels=distance_from_root$edge)
    
    # Time to do some surgery!  Need to repeat data that is both a parent and a child
    parents = unique(tmp$parent)
    children = unique(tmp$edge)
    
    teen_parents = intersect(parents, children)
    
    for (teen_parent in teen_parents) {
        teen_parent_df = tmp[tmp$edge == teen_parent,]
        teen_parent_df$parent = teen_parent
        teen_parent_df$is_parent = "parent"
        
        tmp = rbind(tmp, teen_parent_df)
    }
    # Exclude the zygote parent (non informative because don't have zygote distribution since it's the first node)
    tmp = tmp[tmp$parent != "Zygote",]
    
    tmp$parent = factor(tmp$parent, levels=distance_from_root$edge)
    
    fn = file.path(fig_dir, sprintf("mutation_type_distribution.compare_differences.%s.%s.pdf", tree_name, today))
    if (save_file) {
        pdf(fn, width=width, height=height)
    }
    p = ggplot(tmp, aes(x=edge, y=value)) +
        geom_bar(aes(fill=variable), position = "fill", stat="identity") +
        facet_grid(~parent, scales="free_x", space="free") +
        vartype_palette$fill_scale +
        scale_y_continuous(labels=scales::percent) +
        xlab("Edge") +
        ylab("Percent of mutations") +
        my_theme() +
        theme_classic() +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
    grid.newpage()
    p2 = timestamp_ggfigure(p)
    grid.draw(p2)
    
    if (save_file) {
        dev.off()
        write_new_file_message(fn)
    }
}















