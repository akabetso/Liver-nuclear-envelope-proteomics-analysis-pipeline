####### Here you find functions to use in the main scripts
set.seed(42)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(dplyr)
library(ggplot2)

### aggregate timepoints replicates
levels <- c("CT0", "CT4", "CT8", "CT12", "CT16", "CT20")

### aggregate function
aggregate_timepoints <- function(mat, levels) {
    sample_info <- data.frame(sample = colnames(mat)) %>%
        mutate(
        timepoint = gsub("(_\\d+)$", "", sample), 
        timepoint = factor(timepoint, levels = levels)
        )

    # Average replicates per timepoint
    mat_avg <- mat %>%
        as.data.frame() %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        left_join(sample_info, by = "sample") %>%
        group_by(timepoint) %>%
        summarise(across(where(is.numeric), mean)) %>%
        arrange(timepoint) %>%             
        column_to_rownames("timepoint") %>%
        t()
}

aggregated_timepoints <- function(mat) {
    aggregated <- aggregate_timepoints(mat, levels)
    cat("\nTimepoints replicates aggregatted \n-------------------------------------------")

    return(aggregated)
}



###########################################################
# Plot cycling genes
###########################################################

calculate_phase <- function(y, time, period = 23.6) {
    time_rad <- 2 * pi * time / period
    cos_term <- cos(time_rad)
    sin_term <- sin(time_rad)

    model <- lm(y ~ cos_term + sin_term)
    coefficients <- coef(model)

    # get coeff for cosine and sine terms
    # remember : y = c + A*cos(wt) + B*sin(wt)
    A <- coefficients["cos_term"]
    B <- coefficients["sin_term"]

    # get amplitude and phase
    amplitude <- sqrt(A^2 + B^2)
    phase_rad <- atan2(B, A) # phase in rads

    # Convert phase to hours
    phase_hours <- (phase_rad * period) / (2 * pi)
    if (phase_hours < 0) {
        phase_hours <- phase_hours + period
    }
    
    return(list(phase = unname(phase_hours), amplitude = unname(amplitude)))
}

### calculate the phase for all cycling proteins fdr < cyclers_0.05

get_phase_param <- function(mat) {
    timepoints = seq(0, 20, by = 4)
    result_list <- lapply(rownames(mat), function(gene_name) {
        prot_expr <- as.numeric(mat[gene_name, ])
        phase_result <- calculate_phase(prot_expr, timepoints)
        return(data.frame(
            gene_name = gene_name,
            phase = phase_result$phase,
            amplitude = phase_result$amplitude,
            stringsAsFactors = FALSE
        ))
    })
    phase_cycle <- do.call(rbind, result_list)
    rownames(phase_cycle) <- NULL
    return(phase_cycle)
}


plot_phase_cycles_heatmap <- function(phase_cosine, c_mat, path, title, showrownames) {
    #phase_cosine <- get_phase_param(c_mat)
    #print(head(phase_cosine))
    order_gene <- order(phase_cosine$phase)
    mat_cm <- c_mat[order_gene, ]
    phase_order <- phase_cosine$phase[order_gene]

    mat_scaled <- t(scale(t(mat_cm)))
    col_fun <- colorRamp2(
        c(-2, 2),
        c("blue", "yellow")
    )

    # Phase annotation (left side)
    phase_col_fun <- colorRamp2(
        c(0, 12, 24),
        c("cyan", "red", "green")
    )

    ha_left <- rowAnnotation(
        Phase = phase_order,
        col = list(Phase = phase_col_fun),
        show_annotation_name = FALSE,
        annotation_legend_param = list(
            Phase = list(
                at = c(0, 6, 12, 18, 24),
                labels = c("0", "6", "12", "16", "20"),
                title_gp = gpar(fontsize = 24),
                labels_gp = gpar(fontsize = 18),
                title_gap = unit(1, "cm")
            )
        )
    )
    timepoints <- seq(0, 20, by = 4)
    day_night <- ifelse(timepoints %% 24 < 12, "Day", "Night")

    ha_top <- HeatmapAnnotation(
        Period = day_night,
        col = list(Period = c("Day" = "white", "Night" = "black")),
        show_annotation_name = FALSE,
        show_legend = FALSE,
        annotation_legend_param = list(
          Period = list(
            title_gp = gpar(fontsize = 24),
            labels_gp = gpar(fontsize = 18),
            title_gap = unit(1, "cm")
          )
        )
    )

    ht <- Heatmap(
        mat_scaled,
        name = "Z-score",
        col = col_fun,
        
        # Remove clustering - genes ordered by phase
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        
        # Remove row/column names from heatmap body
        show_row_names = showrownames,
        show_column_names = TRUE,
        column_names_side = "bottom",     
        column_names_rot = 45,
        
        # Add annotations
        left_annotation = ha_left,
        top_annotation = ha_top,
        
        # Appearance
        row_title = paste(nrow(c_mat), "Rhythmic Genes"),
        row_title_side = "left",
        column_title = title,
        column_title_side = "top",    
        column_names_gp   = gpar(fontsize = 24),   
        column_title_gp = gpar(fontsize = 24),
        row_title_gp = gpar(fontsize = 24),
        
        # Border
        border = TRUE,
        
        # Heatmap legend
        heatmap_legend_param = list(
            title = "Expression (Z-score)",
            at = c(-2, -1, 0, 1, 2),
            labels = c("-2", "-1", "0", "1", "2"),
            legend_height = unit(5, "cm"),
            title_gp = gpar(fontsize = 24),
            labels_gp = gpar(fontsize = 18),
            title_gp = gpar(lineheight = 1.5)
        )
    )

    png(path, width = 12, height = 18, res = 300, units = "in")  
    draw(ht,
        merge_legend = FALSE,             
        heatmap_legend_side = "right",    
        annotation_legend_side = "right",  
        legend_gap = unit(2, "cm")         
    ) 
    dev.off() 
}


#################### Make a circular plot of phase and amplitude
phase_amplitude_circular_plot <- function(phase_res, title, path) {
    period = 23.6
    phase_res$phase_rad <- 2 * pi * phase_res$phase / 23.6
    phase_res_binning <- phase_res %>%
        mutate(bin = cut(
            phase_rad,
            breaks = seq(0, 2*pi + 0.107, by = 2*pi/period),
            include.lowest = TRUE,
            right = TRUE
        )) %>%
        group_by(bin) %>%
        summarise(
            amp_mean = sum(amplitude),
            #phase_mid = mean(as.numeric(sub("\\((.+),(.+)\\]", "\\1", bin))) + pi/period#, # bin center
        ) %>%
        mutate(
            # Extract bin boundaries safely
            lower_bound = as.numeric(gsub(".*?(\\d+\\.?\\d*).*", "\\1", bin)),
            upper_bound = as.numeric(gsub(".*,(\\d+\\.?\\d*).*", "\\1", bin)),
            # Calculate midpoint
            phase_mid = (lower_bound + upper_bound) / 2
        ) %>%
        dplyr::select(-lower_bound, -upper_bound) 
    print(phase_res_binning)

    plot_interest <- ggplot(phase_res_binning, aes(x = phase_mid, y = amp_mean)) +
            geom_col(fill = "steelblue", color = "red", alpha = 0.7, width = 2*pi/period * 0.99) +
            #geom_text(aes(label = label, y = amp_mean + 0.5), size = 3, color = "red") +
            coord_polar(start = -pi/2, direction = -1) +
            annotate("rect", xmin = 0, xmax = 2*pi*(6/period), ymin = 0, 
                    ymax = Inf, alpha = 0.1, fill = "black") + 
            annotate("rect", xmin = 2*pi*(18/period), xmax = 2*pi, ymin = 0, 
                    ymax = Inf, alpha = 0.6, fill = "black") + 
            annotate("rect", xmin = 2*pi*(12/period), xmax = 2*pi*(18/period), ymin = 0, 
                ymax = Inf, alpha = 0.6, fill = "black") + 
            annotate("rect", xmin = 2*pi*(6/period), xmax = 2*pi*(12/period), ymin = 0, 
                ymax = Inf, alpha = 0.1, fill = "black") + 
            scale_x_continuous(
                    limits = c(0, 2*pi),
                    breaks = seq(0, 2*pi, by = 2 * pi / period),
                    labels = c(0, "", 2, "", 4, "", 6, "", 8, "", 10, "", 12, "", 14, "", 
                                16, "", 18, "", 20, "", 22, "")) +
                    labs(title = title,
                    x = paste("Time (hours, 23.6)", "\n", "dark (night) | white (day)"),
                    y = "Amplitude\n(sum per bin)") + theme_minimal() +
                    theme(
                    axis.title.y = element_text(vjust = 0.95, angle = 0, size = 18),
                    axis.text.x = element_text(size = 18, color = "black"),
                    plot.title = element_text(size = 18),
                    legend.title = element_text(size = 18),
                    legend.text = element_text(size = 18)
                )
    ggsave(path, plot_interest, bg = "white")
}

################## Cosine line plots 
cosine_line_plot <- function(gene_name, jtk_res, mat, timepoints, path, period = 23.6) {
    data_to_plot <- data.frame(tpts = timepoints)
    protein_expr <- as.numeric(mat[gene_name, ])
    time_rad <- 2 * pi * timepoints / 23.6
    cos_term <- cos(time_rad)
    sin_term <- sin(time_rad)

    model <- lm(protein_expr ~ cos_term + sin_term)
    fitted_vals <- fitted(model)

    data_to_plot$observed = protein_expr
    data_to_plot$fitted = fitted_vals

    q_val <- jtk_res[rownames(jtk_res) == gene_name, ]

    plot <- ggplot(data_to_plot) +
        #geom_point(aes(x = time_points, y = observed, color = "Observed"), size = 3) +
        geom_line(aes(x = tpts, y = fitted, color = "Fitted"), linewidth = 1) +
        scale_color_manual(
            name = "Data Type",
            values = c("Fitted" = "yellow"),
            labels = c("Fitted" = paste0("Fitted\n(p = ", round(as.numeric(q_val$ADJ.P), 2), ")"))) +
        labs(title = paste(gene_name),
                x = "Circadian time (hours)", 
                y = "Expression") +
        theme_minimal()
    ggsave(path, plot, bg = "white", width = 8, height = 4)
}

################## Cosine line plots for per and arntl proteins accross conditions
cosine_line_plot_2 <- function(genes, jtk_res, mat, timepoints, path, title, period = 23.6) {
    data_to_plot <- data.frame(tpts = timepoints)
    for (gene in genes) {
        protein_expr <- as.numeric(mat[gene, , drop = FALSE])
        time_rad <- 2 * pi * timepoints / 23.6
        cos_term <- cos(time_rad)
        sin_term <- sin(time_rad)

        model <- lm(protein_expr ~ cos_term + sin_term)
        fitted_vals <- fitted(model)
        fitted <- paste0("fitted_", gene)

        #data_to_plot[["condition"]] = name
        data_to_plot[[fitted]] = fitted_vals
    }

    p_val_1 <- jtk_res[rownames(jtk_res) == genes[[1]], ]
    p_val_2 <- jtk_res[rownames(jtk_res) == genes[[2]], ]
    p_val_3 <- jtk_res[rownames(jtk_res) == genes[[3]], ]
    p_val_4 <- jtk_res[rownames(jtk_res) == genes[[4]], ]
    p_val_5 <- jtk_res[rownames(jtk_res) == genes[[5]], ]

    plot <- ggplot(data_to_plot) +
        # Fixed: linetype (not line_type), and it should map to a variable or be a constant string
        geom_line(aes(x = tpts, y = fitted_Clock, color = "Clock"), linewidth = 1) +
        geom_line(aes(x = tpts, y = fitted_Bmal1, color = "Bmal1"), linewidth = 1) +
        geom_line(aes(x = tpts, y = fitted_Per1, color = "Per1"), linewidth = 1) +
        geom_line(aes(x = tpts, y = fitted_Per2, color = "Per2"), linewidth = 1) +
        geom_line(aes(x = tpts, y = fitted_Per3, color = "Per3"), linewidth = 1) +
        
        # Color scale for conditions
        scale_color_manual(
            name = "Clock - Bmal - Per complex",
            values = c("Per1" = "blue", "Clock" = "green", "Bmal1" = "red", "Per2" = "yellow", "Per3" = "orange"),
            labels = c(
            "Clock" = paste0("Clock (p-value =", round(p_val_1$ADJ.P, 2), ")"),
            "Bmal1" = paste0("Bmal1 (p-value =", round(p_val_2$ADJ.P, 2), ")"),
            "Per1" = paste0("Per1 (p-value = ", round(p_val_3$ADJ.P, 2), ")"),
            "Per2" = paste0("Per2 (p-value = ", round(p_val_4$ADJ.P, 2), ")"),
            "Per3" = paste0("Per3 (p-value = ", round(p_val_5$ADJ.P, 2), ")")
            )
        ) +
        
        labs(title = title,
            x = "Circadian time (hours)", 
            y = "Expression") +
        theme_minimal()
    ggsave(path, plot, bg = "white", width = 8, height = 4)
}


###### This function is repeated just to change one parameter, add k clusters for 
###### studying elements in cycling heatmap plot. I'm too lazy to modify the above function
###### to include this functionalidty so I copy and modigy it here.
plot_phase_cycles_heatmap_k_clusters <- function(phase_cosine, c_mat, path, title, showrownames, k_clusters = 4) {

    # Order genes by phase
    order_gene   <- order(phase_cosine$phase)
    mat_cm       <- c_mat[order_gene, ]
    phase_order  <- phase_cosine$phase[order_gene]
    mat_scaled   <- t(scale(t(mat_cm)))

    # ── Define clusters by splitting phase range into k equal bins ──
    phase_breaks <- seq(min(phase_order), max(phase_order),
                        length.out = k_clusters + 1)

    cluster_labels <- cut(
        phase_order,
        breaks        = phase_breaks,
        labels        = paste0("C", 1:k_clusters),
        include.lowest = TRUE
    )

    row_split_factor <- factor(cluster_labels,
                               levels = paste0("C", 1:k_clusters))

    # ── Color functions ──
    col_fun <- colorRamp2(c(-2, 2), c("blue", "yellow"))

    phase_col_fun <- colorRamp2(
        c(0, 12, 24),
        c("cyan", "red", "green")
    )

    # ── Annotations ──
    ha_left <- rowAnnotation(
        Phase = phase_order,
        col   = list(Phase = phase_col_fun),
        show_annotation_name = FALSE,
        annotation_legend_param = list(
            Phase = list(
                at        = c(0, 6, 12, 18, 24),
                labels    = c("0", "6", "12", "18", "24"),
                title_gp  = gpar(fontsize = 24),
                labels_gp = gpar(fontsize = 18),
                title_gap = unit(1, "cm")
            )
        )
    )

    timepoints <- seq(0, 20, by = 4)
    day_night  <- ifelse(timepoints %% 24 < 12, "Day", "Night")

    ha_top <- HeatmapAnnotation(
        Period = day_night,
        col    = list(Period = c("Day" = "white", "Night" = "black")),
        show_annotation_name = FALSE,
        show_legend          = FALSE,
        annotation_legend_param = list(
            Period = list(
                title_gp  = gpar(fontsize = 24),
                labels_gp = gpar(fontsize = 18),
                title_gap = unit(1, "cm")
            )
        )
    )

    # ── Heatmap ──
    ht <- Heatmap(
        mat_scaled,
        name             = "Z-score",
        col              = col_fun,

        cluster_rows     = FALSE,   # preserve phase order
        cluster_columns  = FALSE,

        row_split        = row_split_factor,  # visual cluster separation
        row_gap          = unit(3, "mm"),

        show_row_names   = showrownames,
        show_column_names = TRUE,
        column_names_side = "bottom",
        column_names_rot  = 45,

        left_annotation  = ha_left,
        top_annotation   = ha_top,

        row_title        = paste0("C", 1:k_clusters),
        row_title_side   = "left",
        row_title_gp     = gpar(fontsize = 14, fontface = "bold"),

        column_title     = title,
        column_title_side = "top",
        column_names_gp  = gpar(fontsize = 24),
        column_title_gp  = gpar(fontsize = 24),

        border = TRUE,

        heatmap_legend_param = list(
            title         = "Expression (Z-score)",
            at            = c(-2, -1, 0, 1, 2),
            labels        = c("-2", "-1", "0", "1", "2"),
            legend_height = unit(5, "cm"),
            title_gp      = gpar(fontsize = 24),
            labels_gp     = gpar(fontsize = 18)
        )
    )

    # ── Save plot ──
    png(path, width = 12, height = 18, res = 300, units = "in")
    draw(ht,
        merge_legend          = FALSE,
        heatmap_legend_side   = "right",
        annotation_legend_side = "right",
        legend_gap            = unit(2, "cm")
    )
    dev.off()

    # ── Return cluster assignments for downstream use ──
    cluster_df <- data.frame(
        Gene    = rownames(mat_cm),
        Phase   = phase_order,
        Cluster = cluster_labels
    )

    return(invisible(cluster_df))
}