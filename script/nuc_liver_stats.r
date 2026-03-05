
##############################################################################################################
#################### LIVER NUCLEAR ENVELOPE PROTEOME ANALYSIS
set.seed(42)
library("dplyr")
library(tidyverse)
library(ggplot2)
library(VIM)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(limma)
library(tidyverse)
library(ggrepel)
library(gtools)
library(plotrix)
library(imputeLCMD)
library(ReactomePA)
library(enrichplot)
library(httr)
library(jsonlite)
library(igraph)
library(ggraph)
library(ggVennDiagram)
jtk_module <- "script/rythmicity_jtk_cycle.r"


getwd()
list.files("data")
# load
data <- read.delim("data/proteins.txt", stringsAsFactors = FALSE)
colnames(data)

### Drop values that are contaminants
any(data$Contaminant == "VRAI")
table(data$Contaminant == "VRAI") # 51 contaminants found

### Drop the contaminants rows
data <- data[data$Contaminant != "VRAI", ]


# extract column intensities name
head(data$Gene.Symbol)
extract_sample_cols <- function(data) {
    colnames(data)[grepl("Abundances..Normalized", colnames(data))]
}
extract_sample_cols <- extract_sample_cols(data)
new_data <- data[, extract_sample_cols]
colnames(new_data)
head(data$Abundances..Normalized...F2..Sample..CT0)

# clean extracted columns names
### 28 columns ordered according to replicates, 4 each per condition.
state <- 1
for (col_names in colnames(new_data)){
    state <- if (state > 4) 1 else state
    new_name <- sub(".*\\.\\.", "", col_names)
    new_name <- paste0(new_name, "_", state)
    state <- state + 1
    colnames(new_data)[colnames(new_data) == col_names] <- new_name
}
colnames(new_data)
table(is.na(new_data)) # I do not have na's yet but I have missing values or empty spaces
tail(new_data)


### insert gene symbol in the begining
sub_list <- data[, "Gene.Symbol"]
proteins <- cbind(sub_list, new_data) # Complete data cleaning with column needed for downstream analysis
colnames(proteins)[colnames(proteins) == "sub_list"] <- "Gene_names"
colnames(proteins)
tail(proteins)


### Check for duplicated genes if any make unique!
proteins$Gene_names %>% duplicated() %>% any()


### Make a table of duplicated genes
proteins %>% group_by(Gene_names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
# # A tibble: 6 × 2
#   Gene_names frequency
#   <chr>          <int>
# 1 ""                31
# 2 "Cux1"             2
# 3 "H2-D1"            2
# 4 "H2-Eb1"           2
# 5 "Macf1"            2
# 6 "Tmpo"             2

### find if missing names have enseml id.
data$Ensembl.Gene.ID %>% duplicated() %>% any()
table(data$Gene.Symbol == "") 
# 31 genes have no name and even more ensemble id has no name



# Make unique names using the annotation in the "Gene.names" 
# column as primary names and the annotation in "Protein.IDs" 
# as name for those that do not have an gene name.

new_names <- as.character(proteins$Gene_names)
blank_mask <- is.na(new_names) | new_names == "" 
length(which(blank_mask == TRUE)) 
new_names[blank_mask] <- paste0("Unknown_", seq_len(sum(blank_mask)))
new_names <- make.unique(new_names, sep = "_")
proteins$Gene_names <- new_names 

# sanity check: should be FALSE
any(duplicated(proteins$Gene_names)) 
tail(proteins)





                ###############################################
                ######## Generate a SummarizedExperiment object
                ###############################################




################################################ drop perko for wildtype timepoints
# ### PerKO has only one timepoint and is not effective for biological rythms.
# ### Keep three replicates to match perko for equal use in the jtk algo
# drops <- c("Perko_1", "Perko_2", "Perko_3", "CT0_2", "CT8_4", "CT4_2", "CT12_3", "CT16_3", "CT20_1")
# proteins <- proteins[ , !(colnames(proteins) %in% drops)]
# head(proteins)

### Transform wide df to long df format. 
proteins_long <- proteins %>%
    pivot_longer(
        cols = -c(Gene_names),
        names_to = "label",
        values_to = "value"
        ) %>%
    separate(label, into = c("condition", "replicate"), sep = "_") %>%
    mutate(replicate = as.integer(replicate))
head(proteins_long)

### make a summary plot of current state, shows the missing diversity per samples.
# make numerical the value column
proteins_long <- proteins_long %>%
    mutate(value = as.numeric(str_replace(value, ",", ".")))

# create sample id and coerce numeric safely
### remeber, moving from wide to long inserts na in the empty spaces!
proteins_long <- proteins_long %>%
    mutate(sample = paste0(condition, "_", replicate))
head(proteins_long)
table(is.na(proteins_long$value)) 
 

# define detection (presence) rule, 1 for value (presence), 0 for absence(na)
proteins_long <- proteins_long %>%
    mutate(detected = if_else(!is.na(value), 1L, 0L))
tail(proteins_long)

# per-protein detection counts across ALL samples (how many samples/replicates each protein is detected in)
protein_detect <- proteins_long %>%
    group_by(Gene_names) %>%
    summarise(
        detected_in_samples = sum(detected, na.rm = TRUE),
        total_sample    = n_distinct(sample),
        detected_frac = detected_in_samples / total_sample,
        .groups = "drop") 
head(protein_detect)


# frequency table: how many proteins are detected in exactly k samples
### Normally, we would expect the same genes to be expressed in all samples, so the idea
### is to find how many sample have unique expressions. 
### We are projecting the detection abundance per sample counts. For example in 27 samples 3000
### unique proteins were detected or in 1 sample 122 unique proteins were detected. 122 proteins
### detected in 1 sample only is plausibly technical noise. 
which(protein_detect$detected_in_samples == 1)
freq_table <- protein_detect %>%
    group_by(detected_in_samples) %>%
    summarise(sample_dec_count = length(detected_in_samples))
head(freq_table, n = 20)
sum(freq_table$sample_dec_count)
nrow(proteins)

# quick barplot: number of proteins detected in k replicates
plot <- ggplot(freq_table, aes(x = factor(detected_in_samples), y = sample_dec_count)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = sample_dec_count), vjust = -0.5, color = "black", size = 4) +
    labs(x = "K Samples",
        y = "Frequency (proteins)",
        title = "Sample occurence to proteins frequency distribution from raw dataset") +
    theme_minimal()
ggsave("sample_occurence_freq.png", plot, bg = "white", width = 12, height = 12)

# per-replicate summary: how many proteins detected in each sample (useful QC)
per_sample_counts <- proteins_long %>%
    group_by(sample, condition) %>%
    summarise(n_proteins = sum(detected, na.rm = TRUE), .groups = "drop")

pp <- ggplot(per_sample_counts, aes(x = sample, y = n_proteins, fill = condition)) +
    geom_col() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(x = "Sample (condition_replicate)", y = "Frequency (proteins detected)",
        title = "Proteins replicate frequency distribution from raw dataset")
ggsave("prot_replicate_occurence.png", pp, bg = "white", width = 12, height = 12)




#######################################################################
# Firts apply median center normalization before impuation
#######################################################################

### This method helps to correct for sample loading differences by
### centering each samples distribution to the median
proteins_wide <- proteins_long %>%
            dplyr::select(Gene_names, value, sample) %>%
            pivot_wider(names_from = sample, values_from = value) %>%
            column_to_rownames("Gene_names")
head(proteins)


### Convert to matrix for calculations
proteins_mat <- as.matrix(proteins_wide)
head(proteins_mat)
str(proteins_mat)

### Apply log 2 transformation
### There're many transformation that can be done normalise the data(VSN, logtransform...). 
### Data normalization involves transforming variables to bring them to a 
### common scale or comparable state, making them suitable for statistical analysis.
### it Remove technical variation, meet statistical assumptions, enable fair comparisons
### and improve interpretability.
log_intensities <- log2(proteins_mat)
tail(log_intensities)
table(is.na(log_intensities)) 

### Calculate the median for each sample to correct for intensity biassness or
### proteins loading defferences.
### Each sample's protein intensities are shifted so they all have the same median 
### value and removes systematic offset between samples caused by technical factors, 
### aligns distributions so differences reflect biological variation

### First, get the median intensity for each sample.
median_intensities <- apply(log_intensities, 2, median, na.rm = TRUE)

### Make a subtraction of each sample's median to center at zero
normalized_intensities <- sweep(log_intensities, 2, median_intensities, "-")
head(normalized_intensities)
tail(normalized_intensities)

### Check median centering.
check_median <- apply(normalized_intensities, 2, median, na.rm = TRUE)

png("sample_loading_diff.png", width = 8, height = 12, res = 300, units = "in") 
boxplot(log_intensities, main = "Sample loading differences", 
      ylab = "Log2 Intensity", las = 2, cex.axis = 0.7,
      xlab = "Samples")
dev.off()    

# After normalization
png("sample_loading_correction.png", width = 8, height = 12, res = 300, units = "in") 
boxplot(normalized_intensities, main = "Sample loading correction", 
        ylab = "Log2 Intensity", las = 2, cex.axis = 0.7,
        xlab = "Samples")
dev.off()



##############################################################################
# Imputation. First clean values that are less than 50% in each sample
##############################################################################


### The output of mass spectrometry has a lot of missing values described as technical
### limititation due to LOD or simply technical noise for random variables.
### Many methods exist for imputation, but based on the distribution of prots intensities,
### the missing values are mostly driven by LOD so we decide to use Left centered inputation
### methods for this but we tested other methods and those that worked well based on reference
### data is the one we are using here.

### Convert median normalization to long format.
median_normalized_long <- as.data.frame(normalized_intensities) %>%
    rownames_to_column("Gene_names") %>%
    pivot_longer(-Gene_names, names_to = "sample", values_to = "med_norm_vals") %>%
    mutate(sample_temp = sample) %>%
    separate(sample_temp, into = c("condition", "replicate"), sep = "_") %>%
    mutate(replicate = as.integer(replicate))
head(median_normalized_long)

### Now get the detection per condition
median_normalized_long <- median_normalized_long %>%
  mutate(detected = if_else(!is.na(med_norm_vals), 1L, 0L))
print(median_normalized_long, n = 1000, width = Inf)
print(median_normalized_long, n = 10, width = Inf)
table(is.na(median_normalized_long))
table(is.na(proteins_long))
n_distinct(proteins_long)
n_distinct(median_normalized_long)

### Keep only proteins with <= 50% missing values
### Group all the proteins that are in each condition, count those proteins and tell how many
### samples contain each protein in one condition, repeat for all. 
### Since I have 4 replicates per condition, each protein should be found in atleast two samples 
### to be 50% detected in that condition. So, only keep detection >1 in each condition.
proteins_filtered <- median_normalized_long %>%
        group_by(Gene_names, condition) %>%
        mutate(
            n_detections = sum(detected, na.rm = TRUE)
            ) %>%
            filter(n_detections > 1) %>%
        ungroup() 

tail(median_normalized_long)
print(tail(proteins_filtered, 27), n = Inf)
table(proteins_filtered$n_detections)

table(is.na(proteins_filtered))
n_distinct(proteins_filtered)
n_distinct(median_normalized_long)

#### Observe prots per sample for all 24 sample
protein_detect <- proteins_filtered %>%
    group_by(Gene_names) %>%
    summarise(
        detected_in_samples = sum(detected, na.rm = TRUE),
        total_sample = 24,
        detected_frac = detected_in_samples / total_sample,
        .groups = "drop")

table(is.na(proteins_filtered))
print(tail(protein_detect, 27), n = Inf)
tail(protein_detect)
protein_detect[protein_detect$Gene_names == "Myh1", ]

# frequency table: how many proteins are detected in exactly k replicates
freq_table <- protein_detect %>%
  count(detected_in_samples) 

plot <- ggplot(freq_table, aes(x = factor(detected_in_samples), y = n)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = n), vjust = -0.5, color = "black", size = 4) +
    labs(x = "Detections per replicates",
        y = "Number of proteins",
        title = "Sample occurence to proteins frequency distribution corrected") +
    theme_minimal()
ggsave("sample_occurence_freq_correct.png", plot, bg = "white")

### Now convert back to wide and impute
proteins_filtered_mat <- proteins_filtered %>%
    dplyr::select(Gene_names, sample, med_norm_vals) %>%
    pivot_wider(names_from = sample, values_from = med_norm_vals) %>%
    column_to_rownames("Gene_names") %>% 
    as.matrix()
head(proteins_filtered_mat)
tail(proteins_filtered_mat)
proteins_filtered_mat[rownames(proteins_filtered_mat) == "Arel1", ]

# Create the histogram
png("filtered_prots_intensity_distr.png", width = 12, height = 12, res = 300, units = "in")
hist(proteins_filtered_mat, 
    main = "Filtered median normalized log transformed intensity distribution",
    breaks = 50, 
    col = "lightblue")
dev.off()

############################### Imputations #############################################

# probabilistic imputation
### q is the quantile used for estimating missing values. tune.sigma is the standard deviatin 
### multiplier when adding random noise. 
### if low missing values %, keep q low to inpute from low abundance. 
### higher tune.sigma creates more noise and uniformity, avoids spikes, lower creates more 
### outliers, more spikes and less uniformity.
### The parameters below is chosen after varying many others and comparing the normal distributions
### The line below shows the na percentage. helpful for param selection, lower q for lower NAs.
set.seed(42)
(sum(is.na(proteins_filtered_mat)) / (nrow(proteins_filtered_mat) * ncol(proteins_filtered_mat))) * 100
mini_prob <- impute.MinProb(proteins_filtered_mat, q = 0.03, tune.sigma = 2)

# res_qrilc <- impute.QRILC(filtered_proteins_mat, tune.sigma = 1)
# qrilc <- if(!is.null(res_qrilc$dataSet.imputed)) res_qrilc$dataSet.imputed else res_qrilc[[1]]
    
png("imputaion_to_na_dist.png", width = 12, height = 12, res = 300, units = "in")
tail(mini_prob)
par(mfrow = c(1, 2))
hist(proteins_filtered_mat, main="Before imputation (with NAs)")
hist(mini_prob, main="After MinProb imputation")
dev.off()

### plot missing values proportions
png("missing_vals_prop.png", width = 12, height = 12, res = 300, units = "in")
aggr(proteins_filtered_mat, 
    numbers = TRUE,           
    prop = TRUE,               
    cex.axis = 0.7,            
    col = c("skyblue", "red"))
dev.off()


### Plot imputation diagnostics
### Observe the normal distribution of the raw matrix, the imputed matrix and the area imputed 
### or desity of impuation.
visualize_imputation <- function(data, imput_data, save_path) {
    if (!all(dim(data) == dim(imput_data))) stop("Matrices must match.")
    df_orig <- as.data.frame(data) %>% rownames_to_column("Gene")
    head(df_orig)
    df_imp  <- as.data.frame(imput_data) %>% rownames_to_column("Gene")
    df_long <- df_orig %>%
        pivot_longer(-Gene, names_to = "sample", values_to = "orig") %>%
        left_join(df_imp %>% pivot_longer(-Gene, names_to = "sample", values_to = "imputed"),
            by = c("Gene", "sample")) %>%
        mutate(is_imputed = is.na(orig), inputaion = ifelse(is_imputed, "Imputed_matrix", "Imputed_matrix"),
            observation = ifelse(!is_imputed, "Observed_matrix", "NA"),
            other = ifelse(is_imputed, "Missing_value_concentration", "NA"),
            imputed_NAs = ifelse(is_imputed, imputed, NA_real_))

    pp <- ggplot(df_long) +
        geom_density(aes(x = orig, fill = "Observed_matrix"), na.rm = TRUE, alpha = 0.6) +
        geom_density(aes(x = imputed, fill = "Imputed_matrix"), alpha = 0.6) +
        geom_density(aes(x = imputed_NAs, fill = "Missing_value_concentration"), na.rm = TRUE, alpha = 0.6) +
        scale_fill_manual(
            name = "Data Source",
            values = c(Observed_matrix = "steelblue", Imputed_matrix = "tomato", Missing_value_concentration = "green")) +
        labs(x = "Intensities", y = "Density", title = "Imputation density plot (diagnostic)")
    ggsave(save_path, pp, bg = "white", width = 12, height = 12)
    return(df_long)
}
df_long <- visualize_imputation(proteins_filtered_mat, mini_prob, "imputation_density.png")
head(df_long)
tail(df_long)



### Using PCA and Pearson correlation : Capture the similarities and diversities between samples, equally capture outliers.
### Plot PCA
plot_pca <- function(mat_in, sample_info, title){
    # PCA on samples
    pca <- prcomp(t(mat_in), scale. = TRUE, center = TRUE)
    pc_var <- (pca$sdev^2) / sum(pca$sdev^2) * 100 

    df <- as.data.frame(pca$x[,1:2]) %>%
        rownames_to_column("sample") %>%
        left_join(sample_info, by = "sample") %>%
        mutate(replicate = as.factor(replicate)) # Factorized for ggplot as shape takes only categorized data

  ggplot(df, aes(x = PC1, y = PC2, color = condition, shape = replicate)) +
    geom_point(size=4) +
    labs(
        title = paste0("PCA: ", title),
        x = paste0("PC1 (", round(pc_var[1],1), "%)"),
        y = paste0("PC2 (", round(pc_var[2],1), "%)")
    ) +
    theme_minimal()
}

sample_info <- as.data.frame(mini_prob) %>%
    rownames_to_column("Gene_names") %>%
    pivot_longer(-Gene_names, names_to = "sample", values_to = "med_norm_vals") %>%
    mutate(sample_temp = sample) %>%
    separate(sample_temp, into = c("condition", "replicate"), sep = "_")

png("pca_mini_prob.png", width = 800, height = 600)
plot_pca(mini_prob, sample_info, "MinProb")
dev.off()

### Plot pearson correlation matrix
plot_corr_matrix <- function(mat_in, title){
    cor_mat <- cor(mat_in, use="pairwise.complete.obs", method="pearson")
    
    pheatmap(
        cor_mat,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        main = paste0("Correlation matrix: ", title),
        display_numbers = TRUE,
        fontsize_number = 8
    )
}

png("pcorr_mini_prob.png", width = 12, height = 12, res = 300, units = "in")
plot_corr_matrix(mini_prob, "MinProb")
dev.off()


#################### Filter Outliers #####################################

### In PCA, we expect samples to separate by biological condition. 
### Replicates within same group should cluster together.
### Good separation indicates: Strong biological signal, Low technical noise, Successful 
### experimental design.
###
### Poor separation (mixed samples) suggests: Weak biological effect, High technical 
### variability and batch effects dominating biology.

### It is possible that extreme biological outliers may represent genuine rare events 
### or pathological states, but having outliers within replicates is plausibly noise.

### Outliers can be identified via PCA distance, robust Mahalanobis distance, median 
### absolute deviation, etc... and filtered to reduce batch effect within replicates 
### sample seperation are truely biolocally driven.


### Drop sample outliers
drop_sample_outliers <- function(cols, mat){
    mat <- mat[ , !(colnames(mat) %in% cols)]
    return(mat)
}

col_mini_prob <- c("CT0_2", "CT8_4", "Perko_3")
prots_with_no_outlier <- drop_sample_outliers(col_mini_prob, mini_prob)
head(prots_with_no_outlier)
png("pca_corrected_outlier.png", width = 800, height = 600)
plot_pca(prots_with_no_outlier, sample_info, "MinProb")
dev.off()

png("pcorr_outlier_dropped.png", width = 12, height = 12, res = 300, units = "in")
plot_corr_matrix(prots_with_no_outlier, "MinProb")
dev.off()


### Pheatmap visualization
pheatmap_plot <- function(mat, title, show_rownames, file_name) {
        # Parse sample info from column names 
        sample_info <- data.frame(sample = colnames(mat)) %>%
            mutate(
            timepoint = gsub("(_\\d+)$", "", sample), # e.g. CT0_1 → CT0
            timepoint = factor(timepoint,
                levels = c("CT0","CT4","CT8","CT12","CT16","CT20", "Perko")
                )
            )

        # Average replicates per timepoint (factor order)
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

        # Heatmap 
        png_f <- pheatmap(
            mat_avg,
            scale = "row",
            clustering_method = "ward.D2",
            clustering_distance_rows = "manhattan",
            cutree_rows = 4,
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            show_rownames = show_rownames,
            main = title,
            color = colorRampPalette(c("green", "gray", "firebrick3"))(10),
            fontsize_number = 12
        )

        png(file_name, width = 8, height = 12, res = 300, units = "in")  
        print(png_f)
        prots_clusters <- cutree(png_f$tree_row, k = 4)
        dev.off()       
        
        return (prots_clusters)
}

### make plots
title <- "Hierachichal clustering of proteins expression (Averaged replicates)"
prots_clusters1 <- pheatmap_plot(prots_with_no_outlier, title, FALSE, "pheatmap_mini_probb.png")
str(prots_clusters1)

########## Run GO enrichment on clusters
go_enrichment <- function(genes_set, clust) {
    prots_hc <- names(genes_set[genes_set == clust])
    enrich <- enrichGO(gene = prots_hc,
        OrgDb = org.Mm.eg.db,  
        keyType = "SYMBOL",
        ont = "BP")
    cat(head(enrich$Description, 20), sep = "\n")
}
go_enrichment(prots_clusters1, 1)
go_enrichment(prots_clusters1, 2)
go_enrichment(prots_clusters1, 3)
go_enrichment(prots_clusters1, 4)


sample_hc <- function(mat, path) {
    pp <- pheatmap(mat,
        cutree_cols = 7,
        clustering_method = "ward.D2",
        clustering_distance_rows = "manhattan",
        cutree_rows = 4,
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(10),
        scale = "row",
        main = "Hierachichal clustering of proteins expression",
        fontsize_row = 4,  
        show_rownames = FALSE,  
        gaps_col = c(3, 7, 10, 14, 18, 22)
        )
    ggsave(path, pp, bg = "white", width = 12, height = 12)
}

levels <- c("CT0_1", "CT0_3", "CT0_4", 
    "CT4_1", "CT4_2", "CT4_3", "CT4_4", 
    "CT8_1", "CT8_2", "CT8_3", 
    "CT12_1", "CT12_2", "CT12_3", "CT12_4", 
    "CT16_1", "CT16_2", "CT16_3", "CT16_4", 
    "CT20_1", "CT20_2", "CT20_3", "CT20_4", 
    "Perko_1", "Perko_2")
proteins <- prots_with_no_outlier
proteins <- proteins[, levels]
sample_hc(proteins, "replicate_pheatmap.png")

####################### Subset replicates to enforce biologically driven ######
### The dataset were acquiered in one batch, no meta data is avalaible, checking for 
### batch effects becomes difficult. But to ensure data is biologically driven, we
### check for noise and technical variations and limit as much as possible. From unsupervised
### PCA, the're is one outlier which can skew the data, extracting this outlier does'nt provide
### enough information on biological driving factors because replicates are fuzilly clustered. 
### While this can come from variability accross samples (mouses have different physiology), 
### two possible down stream is to average the replicates (done above) and/or to subset 
### a balanced set of highly correlated samples to observed biological principles underlying the 
### observation. 

################# THe averaging perfomed above provides even more better pattern in proteins expresion
###### profiles compared to replicates subset. So we staty with this approach. the gaol is 
###### identifying patterns. ##############################################################



#### Get only the top differentially expressded genes to observe patterns
protein_vars <- apply(prots_with_no_outlier, 1, var)
top_proteins <- names(sort(protein_vars, decreasing = TRUE)[1:500])
mat_subset <- prots_with_no_outlier[top_proteins, ]
tt <- pheatmap(mat_subset,
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_method = "ward.D2",
    show_rownames = FALSE,
    fontsize = 10,
    cutree_rows = 3,
    main = "Top 500 Variable Proteins",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(10))
ggsave("top_variable_prots.png", tt, bg = "white", width = 12, height = 12)
prots_clusters2 <- cutree(tt$tree_row, k = 3)
go_enrichment(prots_clusters2, 1)
go_enrichment(prots_clusters2, 2)
go_enrichment(prots_clusters2, 3)



####################################################################
#       grep protein of interest and viz
####################################################################

plot_protein_of_interest <- function(mat) { 
    prots_interest <- character(0)
    grep_list <- c("Nup", "Per", "H3")  
    for (idx in seq_len(length(grep_list))) {
        val <- grep(grep_list[idx], rownames(mat), value = TRUE)
        prots_interest <- c(prots_interest, val)
    }

    mat_long <- as.data.frame(mat) %>%
        tibble::rownames_to_column("Gene_names") %>%
        pivot_longer(-Gene_names, names_to = "sample", values_to = "value") 

    mat_interest_long <- mat_long %>% 
        filter(Gene_names %in% prots_interest)

    mat_interest_wide <- mat_interest_long %>%
        dplyr::select(Gene_names, sample, value) %>%
        pivot_wider(names_from = sample, values_from = value) %>%
        column_to_rownames("Gene_names") %>%
        as.matrix()
    
    return(mat_interest_wide)
}

############# subset and plot proteins of interest 
prots_interest_mat <- plot_protein_of_interest(prots_with_no_outlier)
title <- "Hierachichal clustering of proteins of interest"
prots_clusters3 <- pheatmap_plot(prots_interest_mat, title, TRUE, "hcl_prots_interest.png")
go_enrichment(prots_clusters3, 1)
go_enrichment(prots_clusters3, 2)
go_enrichment(prots_clusters3, 3)
go_enrichment(prots_clusters3, 4)


################### t-test between CT20 and Perko ##################################

##### averaged the matrix
sample_info <- data.frame(sample = colnames(proteins)) %>%
    mutate(
    timepoint = gsub("(_\\d+)$", "", sample),
    timepoint = factor(timepoint, levels = c("CT0","CT4","CT8","CT12","CT16","CT20", "Perko")))

proteins_avg <- proteins %>%
    t() %>% 
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(sample_info, by = "sample") %>%
    group_by(timepoint) %>%
    summarise(across(where(is.numeric), mean)) %>%
    arrange(timepoint) %>%             
    column_to_rownames("timepoint") %>%
    t() 


head(proteins)
head(proteins_avg)

t_test <- function(dt,grp1,grp2){
    x <- dt[ , grp1] 
    y <- dt[ , grp2]
    ttest <- t.test(x, y)
    results <- tibble(p_val = ttest$p.value, 
        log2FC = mean(x) - mean(y))
    return(results)
}


################# Calculate foldchange and pvalue with limma ####################
set.seed(42)
ct20s  <- c("CT20_1", "CT20_2", "CT20_3", "CT20_4")
perkos <- c("Perko_1", "Perko_2")

mat <- as.data.frame(proteins) %>%
    dplyr::select(all_of(c(ct20s, perkos))) %>% as.matrix()

group <- factor(c(rep("CT20", length(ct20s)), rep("Perko", length(perkos))))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(mat, design)
contrast <- makeContrasts(CT20 - Perko, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
head(fit2$coefficients)
fit_results <- topTable(fit2, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("Gene_names") 
head(fit_results)

####################### Make a volcano plot ####################################
plot_volcano <- function(fit_results, method, path, threshold) {
    fit_results$diffexpressed <- "NO"
    fit_results$diffexpressed[fit_results$logFC > 1.5 & fit_results[[method]] < threshold] <- "UP"
    fit_results$diffexpressed[fit_results$logFC < -1.5 & fit_results[[method]] < threshold] <- "DOWN"
    fit_results$label <- ifelse(fit_results$Gene_names %in% head(fit_results[order(
        fit_results[[method]]), "Gene_names"], 80), fit_results$Gene_names, NA)

    pp <- ggplot(fit_results, aes(logFC, -log10(fit_results[[method]]), col = diffexpressed, label = label)) + 
        geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = "dashed") +
        geom_hline(yintercept = -log10(threshold), col = "gray", linetype = "dashed") +
        geom_point(size = 2) +
        scale_color_manual(values = c("skyblue", "grey", "red"),
            labels = c("Downregulated", "Not significant", "Upregulated")) +
        labs(color = "severe", x = expression("log"[2]*"FC"), y = expression("-log"[10]*p-value)) +
        ggtitle(paste("CT20 vs  PerKO (", method, " < ", threshold, " logfc threshold 1.5)")) + 
        geom_text_repel(max.overlaps = Inf)
    ggsave(path, bg = "white", width = 12, height = 12)
}
plot_volcano(fit_results, "P.Value", "pval_volcano.png", 0.05)
plot_volcano(fit_results, "adj.P.Val", "bh_volcano.png", 0.3)
which(!is.na(fit_results$label))
################# Find GO terms with the top Genes between CT20 Perko ############
get_meta_pathways <- function(fit_results, regulation, threshold_top_variable) {
    fit_results$diffexpresssed <- "NO"
    fit_results$diffexpresssed[fit_results$logFC > 1.5 & fit_results$P.Value < 0.05] <- "UP"
    fit_results$diffexpresssed[fit_results$logFC < -1.5 & fit_results$P.Value < 0.05] <- "DOWN"
    fit_results$label <- ifelse(fit_results$Gene_names %in% head(fit_results[order(
        fit_results$P.Value), "Gene_names"], threshold_top_variable), fit_results$Gene_names, NA)

    up_genes <- fit_results[fit_results$diffexpresssed == regulation, ]
    labeled_gene_pos <- which(!is.na(up_genes$label))
    gene_names <- up_genes$Gene_names[labeled_gene_pos]
    go_enrichment <- function(genes_set) {
        enrich <- enrichGO(gene = genes_set,
            OrgDb = org.Mm.eg.db,  
            keyType = "SYMBOL",
            ont = "BP")
        cat(head(enrich$Description, 10), sep = "\n")
    }
    go_enrichment(gene_names)


    gene_df <- bitr(gene_names, 
                    fromType = "SYMBOL",      
                    toType = "ENTREZID",      
                    OrgDb = org.Mm.eg.db)
    kegg_results <- enrichKEGG(
        gene = gene_df$ENTREZID,
        organism = "mmu",          
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    
    head(kegg_results)
    kegg_results[ , c("category", "subcategory", "Description")]
}
get_meta_pathways(fit_results, "DOWN", 50)
get_meta_pathways(fit_results, "UP", 80)
colnames(fit_results)



#########################################################################################
# Get biological rythmns significance with JTK algorithm developed by Houghes Lab
#########################################################################################

### In order to run rythmicity analysis, we have to drop the PerKO samples because 
### rythmicity follows timepoint and timepoints run from ct0 to ct20. The algorithm 
### requires balanced experimental conditions only so I have to balance the samples    
#### replicates by taking out less correlated replicates.

col_names <- c("Perko_1", "Perko_2", "CT20_3", "CT16_3", "CT12_3", "CT4_2")
prots_timepoints <- proteins[ , !colnames(proteins) %in% col_names]
head(prots_timepoints)

### Save wildtype, then rerun for perko and proceed downstream parallely 
write.table(prots_timepoints, "jtk/proteins.txt", sep = "\t")

### run JTK cycle, will read from the file saved above.
run_jtk_module <- function(module_file) {
    source(module_file, local = FALSE)
}
run_jtk_module(jtk_module)

#### Read from jtk saved text output
jtk_res <- read.csv("jtk/jtk_nucProts.txt", header = TRUE, sep = "\t")
str(jtk_res)

### aggregate replicates in timepoints
path <- "script/utilities.r"
source(path)
avg_timepoints <- aggregated_timepoints(prots_timepoints)
head(avg_timepoints)
str(avg_timepoints, 2)


### plot cycling genes heatmap
phase_heatmap <- function(mat, jtk_res, method, threshold, showrownames, nup_genes = NULL, k_clusters = NULL) {
    cycle_jtk <- subset(jtk_res, jtk_res[[method]] <= threshold)
    name_def <- "_"
    c_gene <- rownames(cycle_jtk)
    if (!is.null(nup_genes)) {
        cycle_jtk <- subset(jtk_res, jtk_res[[method]] <= threshold & rownames(jtk_res) %in% nup_genes)
        cyclers_mat <- mat[rownames(mat) %in% rownames(cycle_jtk), , drop = FALSE]
        name_def <- "_nups"}
    else {cyclers_mat <- mat[rownames(mat) %in% c_gene, , drop = FALSE]}
    phase_cosine <- get_phase_param(cyclers_mat)
    title <- paste("JTK cyclers_", name_def, "-", method, " < 0.05 -", nrow(cycle_jtk), "proteins")
    filename <- paste("jtk/cyclers_", name_def, "_", method, ".png")
    plot_phase_cycles_heatmap(phase_cosine, cyclers_mat, filename, title, showrownames)

    ### Make a pie plot to show period proportions
    label <- round(as.numeric(names(table(cycle_jtk$PER))), 1)
    png(filename = paste("jtk/pc_", name_def, "_", method, ".png"), width = 800, height = 600, bg = "white") 
    pie3D(table(cycle_jtk$PER), labels = paste(label, "hours \n", unique(table(cycle_jtk$PER)), "prots"),
        main = paste0("Rhythmic nuclear ", name_def, " proteins period counts ", method, " < ", threshold), 
        explode=0.25, radius=.7, labelcex = 1.5,  start=2, col = c("blue", "green", "red", "orange"))
    dev.off()

    if(!is.null(k_clusters)) {
        filename <- paste("jtk/cyclers_", name_def, "_kclusters", method, ".png")
        kclus <- plot_phase_cycles_heatmap_k_clusters(phase_cosine, cyclers_mat, filename, title, showrownames, k_clusters)
        return (invisible(kclus))
    }  
}
k_clusters <- phase_heatmap(avg_timepoints, jtk_res, "ADJ.P", 0.05, FALSE, nup_genes = NULL, 3)
phase_heatmap(avg_timepoints, jtk_res, "BH.Q", 0.05, FALSE)

############# Get the meta processes of all cycling proteins
go_enrichment_all_cycle_prots <- function(genes_set, clust) {
    gene_clust <- subset(genes_set, Cluster == clust)$Gene
    enrich <- enrichGO(gene = gene_clust,
        OrgDb = org.Mm.eg.db,  
        keyType = "SYMBOL",
        ont = "BP")
    pp <- enrich@result %>%
        filter(p.adjust < 0.05) %>%
        slice_min(p.adjust, n = 15) %>%
        mutate(Description = fct_reorder(Description, Count)) %>%
        ggplot(aes(x = Count, y = Description, fill = p.adjust)) +
        geom_col() +
        scale_fill_gradient(low = "firebrick3", high = "steelblue") +
        labs(
            title = paste("Cluster", clust),
            x     = "Gene Count",
            y     = NULL,
            fill  = "adj. p-value"
        ) +
        theme_bw() +
        theme(
            plot.title  = element_text(size = 24, face = "bold"),
            axis.text.y = element_text(size = 16)
        )
    ggsave(paste0("jtk/bp_all_cyclers_", clust, ".png"), pp, bg = "white", width = 12, height = 12)

    gene_df <- bitr(gene_clust, 
                    fromType = "SYMBOL",      
                    toType = "ENTREZID",      
                    OrgDb = org.Mm.eg.db)
    kegg_results <- enrichKEGG(
        gene = gene_df$ENTREZID,
        organism = "mmu",          
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    
    head(kegg_results)
    kegg_results[ , c("category", "subcategory", "Description")]
}
go_enrichment_all_cycle_prots(k_clusters, "C1")
go_enrichment_all_cycle_prots(k_clusters, "C2")
go_enrichment_all_cycle_prots(k_clusters, "C3")

###### Make a reactome analysis of the all cycling proteins

cycling_genes <- subset(jtk_res, ADJ.P <= 0.05)
gene_df <- bitr(rownames(cycling_genes), 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Mm.eg.db)

gene_amp_df <- gene_df %>% 
    left_join(jtk_res %>% 
        rownames_to_column("SYMBOL") %>% 
        dplyr::select(SYMBOL, AMP), by = "SYMBOL")
# 22 Genes were not mapped, so keep track of this!

reactome_list <- gene_amp_df %>%
    arrange(desc(AMP)) %>%
    { setNames(.$AMP, .$ENTREZID) } 

y <- gsePathway(
    reactome_list, 
    organism = "mouse",
    minGSSize = 5,     
    maxGSSize = 500, pvalueCutoff=1)
res <- as.data.frame(y)
head(res)
head(jtk_res)


y2 <- pairwise_termsim(y)  
pp <- emapplot(y2, showCategory = 30)
ggsave("jtk/reactome_gse_pw.png", pp, bg = "white", width = 12, height = 12)

# readmoe here on reactome : 
# https://bioconductor.statistik.uni-dortmund.de/packages/3.1/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
de <- names(reactome_list)
x <- enrichPathway(gene=de,pvalueCutoff = 0.05, organism = "mouse", readable=T)
x2 <- pairwise_termsim(x)
p1 <- emapplot(x2, showCategory = 30)
ggsave("jtk/reactome_enrich_pw.png", p1, bg = "white", width = 12, height = 12)

############## Search for transcription factors driving the genes of interest

#### Get TF for cycle proteins clusters
chea_TF_cycle_prots_clusters <- function(gene_set, path, clust) {
    url = "https://maayanlab.cloud/chea3/api/enrich/"
    encode = "json"
    payload = list(query_name = "myQuery", gene_set = gene_set)
    response = POST(url = url, body = payload, encode = encode)
    json = content(response, "text")
    results = fromJSON(json)

    top_tfs <- results$`Integrated--meanRank` %>%
        as.data.frame() %>%
        mutate(Rank = as.numeric(Rank)) %>%
        slice_min(Rank, n = 20) %>%
        mutate(TF = fct_reorder(TF, -Rank))

    tt <- ggplot(top_tfs, aes(x = Rank, y = TF, fill = Rank)) +
        geom_col() +
        scale_fill_gradient(low = "firebrick3", high = "steelblue") +
        labs(
            title = paste("Top 20 Transcription Factors (ChEA3) - Cluster ", clust),
            x     = "Mean Rank (lower = better)",
            y     = NULL
        ) +
        theme_bw() +
        theme(
            plot.title  = element_text(face = "bold", size = 24),
            axis.text.y = element_text(size = 16),
            axis.text.x = element_text(size = 16),
            legend.text = element_text(size = 16),  
            axis.title.x = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 16, face = "bold")
        )
    ggsave(path, tt, bg = "white", width = 12, height = 12)
    
    return (results)
}
cluster1 <- subset(k_clusters, Cluster == "C1")$Gene
chea_clust1_res <- chea_TF_cycle_prots_clusters(cluster1, "jtk/cyclers_clus1_cheaTF.png", 1)
cluster2 <- subset(k_clusters, Cluster == "C2")$Gene
chea_clust2_res <- chea_TF_cycle_prots_clusters(cluster2, "jtk/cyclers_clus2_cheaTF.png", 2)
cluster3 <- subset(k_clusters, Cluster == "C3")$Gene
chea_clust3_res <- chea_TF_cycle_prots_clusters(cluster3, "jtk/cyclers_clus3_cheaTF.png", 3)


### Plot the TF gene network
chea_TF_gene_network_clust <- function(chea_clust_res, path, clust) {
    tf_genes <- chea_clust_res$`Integrated--meanRank` %>%
        as.data.frame() %>%
        mutate(Rank = as.numeric(Rank)) %>%
        slice_min(Rank, n = 5) %>%      
        rowwise() %>%
        mutate(genes = list(strsplit(Overlapping_Genes, ",")[[1]])) %>%
        unnest(genes) %>%
        dplyr::select(TF, genes)

    g <- graph_from_data_frame(tf_genes, directed = TRUE)
    V(g)$type <- ifelse(V(g)$name %in% tf_genes$TF, "TF", "Gene")

    tt <- ggraph(g, layout = "fr") +
        geom_edge_link(alpha = 0.3, color = "gray60") +
        geom_node_point(aes(color = type, size = type)) +
        geom_node_text(aes(label = name, color = type),
                    repel = TRUE, size = 4) +
        scale_color_manual(values = c("TF" = "firebrick3", "Gene" = "steelblue")) +
        scale_size_manual(values  = c("TF" = 6, "Gene" = 3)) +
        labs(title = paste("Top 5 TF-Gene Regulatory Network - Cluster ", clust)) +
        theme_graph()
    ggsave(path, tt, bg = "white", width = 12, height = 12)
}
chea_TF_gene_network_clust(chea_clust1_res, "jtk/cyclers_clus1_cheaTF_network.png", 1)
chea_TF_gene_network_clust(chea_clust2_res, "jtk/cyclers_clus2_cheaTF_network.png", 2)
chea_TF_gene_network_clust(chea_clust3_res, "jtk/cyclers_clus3_cheaTF_network.png", 3)


#### Get list of proteins of prots of interest and make phase heatmap 
prots_interest <- character(0)
grep_list <- c("Nup", "Per", "H3")  
for (idx in seq_len(length(grep_list))) {
    val <- grep(grep_list[idx], rownames(avg_timepoints), value = TRUE)
    prots_interest <- c(prots_interest, val)
}
phase_heatmap(avg_timepoints, jtk_res, "ADJ.P", 0.05, TRUE, prots_interest)

############ Make a phase amplitude circular plot for cycling proteins
cyclers_polar_plot <- function(mat, jtk_res, method, threshold, nups = NULL){
    def_name <- ""
    cycle_jtk <- subset(jtk_res, jtk_res[[method]] <= threshold)
    c_gene <- rownames(cycle_jtk)
    cyclers_mat <- mat[rownames(mat) %in% c_gene, , drop = FALSE]
    if (!is.null(nups)) {
        cycle_jtk <- subset(jtk_res, jtk_res[[method]] <= threshold & rownames(jtk_res) %in% nups)
        cyclers_mat <- mat[rownames(mat) %in% rownames(cycle_jtk), , drop = FALSE]
        def_name <- "nups"
    }
    phase_cosine <- get_phase_param(cyclers_mat)
    title <- paste("Rhthmic", def_name, "proteins phase/Amplitude \nditribution - ", method, " < ", threshold, 
        nrow(cycle_jtk), "proteins")
    filename <- paste0("jtk/amp_dist_", def_name, method, ".png")
    phase_amplitude_circular_plot(phase_cosine, title, filename)
}
cyclers_polar_plot(avg_timepoints, jtk_res, "ADJ.P", 0.05)
cyclers_polar_plot(avg_timepoints, jtk_res, "BH.Q", 0.05)
cyclers_polar_plot(avg_timepoints, jtk_res, "ADJ.P", 0.05, prots_interest)

### check the meta biological processes of cyling proteins of interest.
get_pathways_prots_interest <- function(genes) {
    enrich <- enrichGO(gene = genes,
        OrgDb = org.Mm.eg.db,  
        keyType = "SYMBOL",
        ont = "BP")
    cat(head(enrich$Description, 10), sep = "\n")


    gene_df <- bitr(genes, 
                    fromType = "SYMBOL",      
                    toType = "ENTREZID",      
                    OrgDb = org.Mm.eg.db)
    kegg_results <- enrichKEGG(
        gene = gene_df$ENTREZID,
        organism = "mmu",          
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    
    head(kegg_results)
    kegg_results[ , c("category", "subcategory", "Description")]
}
get_pathways_prots_interest(c("Nup37", "Nup35", "Nup62"))


############# make a line plot for reference proteins
clock_transcripts <- c("Per1", "Cry2", "Clock", "Bmal1", "Per2", "Per3", "Nr1d1", "Nr1d2", "Dbp", "Hlf", "Tef", "Nfil3")
for (clock_gene in clock_transcripts) {
    cosine_line_plot(clock_gene, jtk_res, avg_timepoints, seq(0, 20, by = 4), 
        paste0("jtk/", clock_gene, ".png"), "Dynamics")
}

#### Repeat the same line plot from reference proteins to proteins of interest.
prots_interest <- character(0)
grep_list <- c("Nup", "Per", "H3")  
for (idx in seq_len(length(grep_list))) {
    val <- grep(grep_list[idx], rownames(avg_timepoints), value = TRUE)
    prots_interest <- c(prots_interest, val)
}
for (clock_gene in prots_interest) {
    cosine_line_plot(clock_gene, jtk_res, avg_timepoints, seq(0, 20, by = 4), 
        paste0("jtk/nups/", clock_gene, ".png"), "Dynamics")
}

############# make another line plot for reference proteins for observing phase diff between per and clock
#### Make a cosine line plot for nups proteins

path <- "script/utilities.r"
source(path)
clock_transcripts <- c("Clock", "Bmal1", "Per1", "Per2", "Per3")
cosine_line_plot_2(clock_transcripts, jtk_res, avg_timepoints, seq(0, 20, by = 4), 
    paste0("jtk/line_plot_feedback_complex", ".png"), "Clock-Bmal-Per Feedback Complex")

########### To end the analysis, download all liver cycling genes from literature and 
########### cross check detection with ours using a venn diagram. rythmic genes downloaded from :
########### https://pmc.ncbi.nlm.nih.gov/articles/PMC3694775/#S7
literature_cycling_genes <- read.delim("data/rythmic_genes_literature.txt", stringsAsFactors = FALSE)
colnames(literature_cycling_genes)
literature_cycling_genes <- literature_cycling_genes["Gene.names"]
head(literature_cycling_genes)
nrow(literature_cycling_genes)
nrow(cycling_genes)

cross_check_cyclers <- list(Takahashi = literature_cycling_genes$Gene.names, 
    Kiran = rownames(cycling_genes))
head(cross_check_cyclers)

venn_plot <- ggVennDiagram(cross_check_cyclers, set_size = 3.5, label_size = 4,) +
    ggplot2::ggtitle("Cross check cycling proteins\n found in litereature and this dataset") +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::theme(
        plot.title  = ggplot2::element_text(hjust = 0.5, size = 20),
        plot.margin = margin(5, 5, 5, 5)
    )

ggsave("jtk/cross_check_cycling_proteins.png", venn_plot, width = 8, height = 6, bg = "white")



############################# END OF ANALYSIS #########################################
## important sites :
# https://bioconductor.statistik.uni-dortmund.de/packages/3.1/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
# https://maayanlab.cloud/chea3/
# https://ab604.github.io/docs/bspr_workshop_2018/transform.html
# https://www.cell.com/cell/fulltext/S0092-8674(19)31003-7#mmc3)?and
# https://github.com/41ison/SAM-for-proteomics
