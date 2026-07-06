library(readxl)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(ggvenn)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Dm.eg.db)
library(stringr)

# uploaded the table
rowData <- as.data.frame(read_excel("data.xlsx"))
# deleted the first 2 lines with technical information
data <- rowData[-(1:2), ]
# row name = the text in the first column
rownames(data) <- make.unique(as.character(data[, 1]))
# removed the first column
data <- data[ ,-1]

# preparation of metadata
coldata <- data.frame(
  row.names = c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"),
  Genotype = factor(c("Het","Het","Het","Het", "Hom","Hom","Hom","Hom", "Het","Het","Het","Het", "Hom","Hom","Hom","Hom")),
  Diet     = factor(c("Complete","Complete","Complete","Complete", "Complete","Complete","Complete","Complete", "HalfTyr","HalfTyr","HalfTyr","HalfTyr", "HalfTyr","HalfTyr","HalfTyr","HalfTyr"))
)

# checking whether columns in data and coldata match
all(colnames(data) == rownames(coldata))

# converted the table from a text format to a numeric
data_numeric <- matrix(as.numeric(as.matrix(data)), 
                       nrow = nrow(data), 
                       dimnames = dimnames(data))

# validations and synchronization of countData with colData, initialization of result matrices, Writing of the design matrix (~ Genotype + Diet + Genotype:Diet)
dds <- DESeqDataSetFromMatrix(
  countData = data_numeric,
  colData   = coldata,
  design    = ~ Genotype + Diet + Genotype:Diet
)

# Size factor estimation, variance estimation, GLM construction and Wald test
dds <- DESeq(dds)

# we save intermediate results Hom vs Het (100% Tyr) 
res_disease_100 <- results(dds, contrast = c("Genotype", "Hom", "Het"), alpha = 0.05)
# we save intermediate results Hom vs Het (50% Tyr) 
res_disease_50 <- results(dds, contrast = list(c("Genotype_Hom_vs_Het", "GenotypeHom.DietHalfTyr")), alpha = 0.05)
# we save intermediate results Hom 50% Tyr vs Hom 100% Tyr
res_diet_hom <- results(dds, contrast = list(c("Diet_HalfTyr_vs_Complete", "GenotypeHom.DietHalfTyr")), alpha = 0.05)
# we save intermediate results Het 50% Tyr vs Het 100% Tyr
res_diet_het <- results(dds, contrast = c("Diet", "HalfTyr", "Complete"), alpha = 0.05)



# ____________________________Volcano plot______________________________________
make_volcano <- function(res_data, plot_title) {

# we write down the 10 most statistically significant genes
top_genes <- head(rownames(res_data[order(res_data$padj), ]), 10)

# Volcano plot preparation
p <- EnhancedVolcano(
  res_data,
  lab = rownames(res_data),              # take the names of the genes from the names of the lines
  x = 'log2FoldChange',                  # The X-axis is the magnitude of the effect
  y = 'padj',                            # Y-axis - adjusted p-value (FDR)
  selectLab = top_genes,
  pCutoff = 0.05,                        # Horizontal line of significance
  FCcutoff = 1.0,                        # Vertical lines of effect power (doubled)
  pointSize = 0.5,                       # The size of gene points
  labSize = 2.0,                         # The size of the text with the names of the top genes
  title = plot_title,
  subtitle = 'Volcano plot',
  legendPosition = 'right'
)

p_final <- p + 
  coord_cartesian(
    xlim = c(-4, 4),   # zoom in on the X axis
    ylim = c(0, 30)    # zoom in on the Y axis
  ) +
  # we make a smaller grid
  scale_x_continuous(breaks = seq(-4, 4, by = 1)) +
  scale_y_continuous(breaks = seq(0, 30, by = 5)) 

return(p_final)
}

# save the Volcano plot option
# Disease effect 100% Tyr
p1 <- make_volcano(res_disease_100, 'Disease Effect: Hom vs Het (100% Tyr)')
ggsave("Volcano_Disease_100.png", plot = p1, width = 12, height = 9, dpi = 300)
# Disease effect 50% Tyr
p2 <- make_volcano(res_disease_50, 'Disease Effect: Hom vs Het (50% Tyr)')
ggsave("Volcano_Disease_50.png", plot = p2, width = 12, height = 9, dpi = 300)
# Diet effect on Het
p3 <- make_volcano(res_diet_het, 'Diet Effect: 50% vs 100% Tyr (Het)')
ggsave("Volcano_Diet_Het.png", plot = p3, width = 12, height = 9, dpi = 300)
# Diet effect on Hom
p4 <- make_volcano(res_diet_hom, 'Diet Effect: 50% vs 100% Tyr (Hom)')
ggsave("Volcano_Diet_Hom.png", plot = p4, width = 12, height = 9, dpi = 300)



# ____________________________MA plot___________________________________________
# Variance Stabilizing Transformation
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd) # Pure matrix of normalized values

make_ma_plot <- function(res_data, plot_title) {
  
  # convert the results into a dataframe
  df <- as.data.frame(res_data)
  
  # We create a column for color (Up, Down, NS)
  df$Significance <- "NS"
  df$Significance[df$padj < 0.05 & df$log2FoldChange > 1] <- "Up"
  df$Significance[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down"
  
  # Defense against NA in padj
  df$Significance[is.na(df$Significance)] <- "NS"
  
  # count the exact number for the legend
  n_up <- sum(df$Significance == "Up")
  n_down <- sum(df$Significance == "Down")
  n_ns <- sum(df$Significance == "NS")
  
  # Convert into a factor with clear levels
  df$Significance <- factor(df$Significance, levels = c("Up", "Down", "NS"))
  df$GeneName <- rownames(df)
  
  # find the TOP-10 genes for signatures
  top_genes_df <- df[!is.na(df$padj), ]
  top_genes <- head(rownames(top_genes_df[order(top_genes_df$padj), ]), 10)
  df$GeneName <- rownames(df)
  
  # build a plot
  # logarithmize the X axis
  p <- ggplot(df, aes(x = baseMean, y = log2FoldChange, color = Significance)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_x_log10(
      breaks = c(1, 10, 100, 1000, 10000, 100000),
      labels = c("1", "10", "100", "1K", "10K", "100K")
    ) +
    coord_cartesian(ylim = c(-4, 4)) +
    scale_color_manual(
      values = c("Up" = "#B31B21", "Down" = "#1465B2", "NS" = "darkgray"),
      labels = c("Up" = paste0("Up (", n_up, ")"), 
                 "Down" = paste0("Down (", n_down, ")"), 
                 "NS" = paste0("NS (", n_ns, ")"))
    ) +
    
    # Dashed threshold lines
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = 0, color = "black") +
    
    # Axis labels and headings
    labs(
      title = plot_title,
      subtitle = paste0("Up: ", sum(df$Significance == "Up"), 
                        " | Down: ", sum(df$Significance == "Down")),
      x = "Mean of Normalized Counts (Log10 scale)",
      y = "Log2 Fold Change"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    )

  # add the signatures of the TOP-10 genes
  library(ggrepel)
  df_top <- df[df$GeneName %in% top_genes & df$Significance != "NS", ]
  if(nrow(df_top) > 0) {
    p <- p + geom_label_repel(
      data = df_top,
      aes(label = GeneName),
      size = 2.5,
      label.padding = unit(0.15, "lines"),
      box.padding = unit(0.3, "lines"),
      point.padding = unit(0.2, "lines"),
      color = "black", fill = "white", fontface = "bold",
      show.legend = FALSE
    )
  }
  
  return(p)
}

# MA plot for Hom vs Het (100% Tyr)
ma1 <- make_ma_plot(res_disease_100, 'MA Plot: Disease Effect - Hom vs Het (100% Tyr)')
ggsave("MA_Disease_100.png", plot = ma1, width = 11, height = 8, dpi = 300, bg = "white")

# MA plot for Hom vs Het (50% Tyr)
ma2 <- make_ma_plot(res_disease_50, 'MA Plot: Disease Effect - Hom vs Het (50% Tyr)')
ggsave("MA_Disease_50.png", plot = ma2, width = 11, height = 8, dpi = 300, bg = "white")

# MA plot for Het 50% Tyr vs Het 100% Tyr
ma3 <- make_ma_plot(res_diet_het, 'MA Plot: Diet Effect - Het 50% vs 100% Tyr')
ggsave("MA_Diet_Het.png", plot = ma3, width = 11, height = 8, dpi = 300, bg = "white")

# MA plot for Hom 50% Tyr vs Hom 100% Tyr
ma4 <- make_ma_plot(res_diet_hom, 'MA Plot: Diet Effect - Hom 50% vs 100% Tyr')
ggsave("MA_Diet_Hom.png", plot = ma4, width = 11, height = 8, dpi = 300, bg = "white")



# ____________________Heatmap of significant genes______________________________
vst_mat <- assay(vst(dds, blind = FALSE))

make_heatmap <- function(vst_matrix, gene_list, plot_title, file_name) {
  
  # if there is no list, it is empty or there are less than 2 genes - exit
  if (is.null(gene_list) || length(gene_list) < 2) {
    message("Too few genes (< 2) are listed for: ", plot_title)
    return(invisible(FALSE))
  }
  
  # take from the matrix only those genes that passed the filter
  heatmap_data <- vst_matrix[rownames(vst_matrix) %in% gene_list, ]
  
  # if there are almost no significant genes, the Heatmap is not built
  if (is.null(heatmap_data) || nrow(heatmap_data) < 2) {
    message("Too few genes to build a Heatmap:", plot_title)
    return(invisible(FALSE))
  }
  
  # sorting columns by names
  heatmap_data <- heatmap_data[, order(colnames(heatmap_data))]
  
  # take the metadata from dds for the colored bars on top
  annotation_col <- as.data.frame(colData(dds)[, c("Genotype", "Diet")])
  
  # set colors for groups
  ann_colors <- list(
    Genotype = c(Het = "#1465B2", Hom = "#B31B21"),
    Diet = c(Complete = "#2CA02C", HalfTyr = "#FF7F0E")
  )
  
  # color palette: from blue (low expression) to red (higher)
  display_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  
  # draw and save
  pheatmap(
    heatmap_data,
    cluster_rows = TRUE,           # group similar genes together
    cluster_cols = FALSE,          # shows whether repeats are grouped together
    scale = "row",                 # converts to Z-score (relative change)
    show_rownames = (nrow(heatmap_data) <= 50), # hide the names if there are too many genes
    show_colnames = TRUE,          # show the ciphers
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    color = display_colors,
    filename = file_name,          # saving
    width = 10,
    height = 8
  )
}

genes_dis100 <- head(rownames(res_disease_100[order(res_disease_100$padj), ]), 50)
make_heatmap(vst_mat, genes_dis100, 
             "Heatmap: Top 50 DEGs - Hom vs Het (100% Tyr)", 
             "Heatmap_Disease_100.png")

genes_dis50 <- head(rownames(res_disease_50[order(res_disease_50$padj), ]), 50)
make_heatmap(vst_mat, genes_dis50, 
             "Heatmap: Top 50 DEGs - Hom vs Het (50% Tyr)", 
             "Heatmap_Disease_50.png")

genes_diet_hom <- rownames(subset(res_diet_hom, padj < 0.05 & abs(log2FoldChange) > 1))
make_heatmap(vst_mat, genes_diet_hom, 
             "Heatmap: All DEGs - Diet effect Hom 50% Tyr vs Hom 100% Tyr", 
             "Heatmap_Diet_Hom.png")

genes_diet_het <- rownames(subset(res_diet_het, padj < 0.05 & abs(log2FoldChange) > 1))
make_heatmap(vst_mat, genes_diet_het, 
             "Heatmap: All DEGs - Diet Effect Het 50% Tyr vs Het 100% Tyr", 
             "Heatmap_Diet_Het.png")



# =========================Gene Set Enrichment Analysis=========================
make_gsea_analysis <- function(res_data, plot_title, file_prefix) {
  
  df <- as.data.frame(res_data)
  # remove genes from NA
  df <- df[!is.na(df$log2FoldChange), ]
  
  # calculate the ranking metric (Rank = shift sign * (-log10 of p-value))
  df$ranking_metric <- sign(df$log2FoldChange) * (-log10(df$pvalue))
  # Protection against infinite values if p-value was 0
  df$ranking_metric[is.na(df$ranking_metric) | is.infinite(df$ranking_metric)] <- 
    sign(df$log2FoldChange[is.na(df$ranking_metric) | is.infinite(df$ranking_metric)]) * 300
  
  # sort the genes in descending order of rank
  df_sorted <- df[order(-df$ranking_metric), ]
  
  # form a named vector for clusterProfiler
  gene_list <- df_sorted$ranking_metric
  names(gene_list) <- rownames(df_sorted)
  
  message("Running gseGO for Drosophila: ", plot_title, "... this may take 1-2 minutes.")
  
  # org.Dm.eg.db automatically recognizes names
  gsea_res <- clusterProfiler::gseGO(
    geneList     = gene_list,
    OrgDb        = org.Dm.eg.db, # database
    Ontology     = "BP",         # Biological Process     
    keyType      = "SYMBOL",     
    minGSSize    = 10,          
    maxGSSize    = 500,         
    pvalueCutoff = 0.05,        
    verbose      = FALSE,
    pAdjustMethod = "BH"         # Benjamini-Hochberg correction
  )
  
  # if no path is enriched - leave
  if (is.null(gsea_res) || nrow(gsea_res) == 0) {
    message("No significant biological pathways were found for: ", plot_title)
    return(NULL)
  }
  
  # Graph 1: Dotplot
  # Shows activated (red) and suppressed (blue) processes
  p_dot <- dotplot(gsea_res, showCategory = 20, split = ".sign") + 
    labs(title = paste0("GSEA (Drosophila): ", plot_title)) +
    theme(plot.title = element_text(face = "bold", size = 12))
  
  ggsave(paste0(file_prefix, "_GSEA_Dotplot.png"), plot = p_dot, 
         width = 11, height = 20, dpi = 300, bg = "white")
  
  # Graph 2: Top Pathway. GSEA curve-graph (Enrichment Score) for the top-1 strongest path
  p_gsea <- enrichplot::gseaplot2(gsea_res, geneSetID = 1, title = gsea_res$Description[1])
  
  ggsave(paste0(file_prefix, "_GSEA_TopPath.png"), plot = p_gsea, 
         width = 9, height = 7, dpi = 300, bg = "white")
  
  return(gsea_res)
}

# Hom vs Het (100% Tyr)
gsea_disease_100 <- make_gsea_analysis(
  res_data = res_disease_100, 
  plot_title = "Disease Effect - Hom vs Het (100% Tyr)", 
  file_prefix = "Disease_100"
)

# Hom vs Het (50% Tyr)
gsea_disease_50 <- make_gsea_analysis(
  res_data = res_disease_50, 
  plot_title = "Disease Effect - Hom vs Het (50% Tyr)", 
  file_prefix = "Disease_50"
)

# Hom 50% Tyr vs Hom 100% Tyr
gsea_diet_hom <- make_gsea_analysis(
  res_data = res_diet_hom, 
  plot_title = "Diet Effect - Hom 50% vs 100% Tyr", 
  file_prefix = "Diet_Hom"
)

# Het 50% Tyr vs Het 100% Tyr
gsea_diet_het <- make_gsea_analysis(
  res_data = res_diet_het, 
  plot_title = "Diet Effect - Het 50% vs 100% Tyr", 
  file_prefix = "Diet_Het"
)



# ======================Over-representation analysis (ORA)======================
make_ora_analysis <- function(res_data, plot_title, file_prefix) {
  
  df <- as.data.frame(res_data)
  
  # select only significant genes. We take the standard threshold padj < 0.05
  sig_genes <- rownames(df[which(df$padj < 0.05), ])
  
  # if there are no significant genes or there are less than 2-3 of them, the analysis will not start
  if (length(sig_genes) < 2) {
    message("Too few significant genes for ORA in contrast: ", plot_title)
    return(NULL)
  }
  
  message("Run ORA (enrichGO) for: ", plot_title, ". Number of genes: ", length(sig_genes))
  
  # Running a saturation analysis
  ora_res <- enrichGO(
    gene          = sig_genes,
    OrgDb         = org.Dm.eg.db, # database
    keyType       = "SYMBOL",     # gene names
    pAdjustMethod = "BH",         # Benjamini-Hochberg correction
    pvalueCutoff  = 0.05,         # leave only significant paths
    qvalueCutoff  = 0.2           # additional reliability filter
  )
  
  if (is.null(ora_res) || nrow(ora_res) == 0) {
    message("No significant processes were found in ORA for: ", plot_title)
    return(NULL)
  }
  
  # Visualization
  ora_res_plot <- ora_res
  ora_res_plot@result$Description <- stringr::str_wrap(ora_res_plot@result$Description, width = 40)
  
  p_bar <- barplot(ora_res_plot, showCategory = 12) + 
    labs(title = paste0("ORA Barplot: ", plot_title)) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8, lineheight = 0.8))
  
  ggsave(paste0(file_prefix, "_ORA_Barplot.png"), plot = p_bar, 
         width = 9, height = 8, dpi = 300, bg = "white")
  
  return(ora_res)
}

# Hom vs Het (100% Tyr)
ora_disease_100 <- make_ora_analysis(res_disease_100, "Disease Effect - 100% Tyr", "Disease_100")

# Hom vs Het (50% Tyr)
ora_disease_50 <- make_ora_analysis(res_disease_50, "Disease Effect - 50% Tyr", "Disease_50")

# Hom 50% Tyr vs Hom 100% Tyr
ora_diet_hom <- make_ora_analysis(res_diet_hom, "Diet Effect - Hom 50% vs 100%", "Diet_Hom")

# Het 50% Tyr vs Het 100% Tyr
ora_diet_het <- make_ora_analysis(res_diet_het, "Diet Effect - Het 50% vs 100%", "Diet_Het")



# ============================Specific gene Heatmap=============================
library(KEGGREST)

# Let's start the variance stabilization transform
vst_counts <- DESeq2::vst(dds, blind = FALSE)

# ----------------------group of gene by coded----------------------------------
kegg_pathway <- keggGet("dme00280")
kegg_genes <- kegg_pathway[[1]]$GENE

# extract purely SYMBOL gene names (they go in pairs with ID in KEGG)
valine_genes_raw <- kegg_genes[seq(2, length(kegg_genes), by = 2)]
# remove system signatures, leaving clean gene names
valine_gene_symbols <- sapply(strsplit(valine_genes_raw, ";"), "[", 1)

# filter the normalized VST-data purely for this list of genes
available_genes <- intersect(valine_gene_symbols, rownames(vst_counts))
valine_matrix <- assay(vst_counts)[available_genes, ]

# draw the target heatmap
pheatmap(
  valine_matrix, 
  scale = "row",          # normalize by rows to see shifts in expression
  show_rownames = TRUE,   # show the names of the genes
  cluster_cols = FALSE,   # do not mix the columns
  annotation_col = as.data.frame(colData(dds)[, c("Genotype", "Diet")]), # add group signatures
  main = "Expression of Valine Metabolism Genes",
  filename = "Code_Target_Genes_Heatmap.png",
  width = 7,
  height = 4 
)


# -----------------------------target gene--------------------------------------
# enter the names of the genes manually
my_target_genes <- c("Hibch", "Thor")

# check whether these genes are present in the normalized data
available_genes <- intersect(my_target_genes, rownames(vst_counts))

# Check if the genes are spelled correctly
if(length(available_genes) < length(my_target_genes)) {
  message("WARNING! Some genes were not found in the sample matrix. Check your spelling.")
}

# extract the matrix purely for these genes
valine_matrix <- assay(vst_counts)[available_genes, , drop = FALSE] 
# Note: drop = FALSE is necessary so that R does not break the structure of the matrix when there are so few genes

# Draw a Heatmap
pheatmap(
  valine_matrix, 
  scale = "row",          
  show_rownames = TRUE,   
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  annotation_col = as.data.frame(colData(dds)[, c("Genotype", "Diet")]), 
  main = "Expression of Selected Target Genes",
  filename = "Target_Genes_Heatmap.png",
  width = 7,
  height = 4 
)



#===========================Interaction analysis================================
# extract the results for interaction
res_interaction <- results(dds, name = "GenotypeHom.DietHalfTyr")

# translate into a table and filter significant interaction genes
df_int <- as.data.frame(res_interaction)
sig_interaction_genes <- df_int[which(df_int$padj < 0.05 & abs(df_int$log2FoldChange) > 1), ]

message("The number of genes with an interaction effect: ", nrow(sig_interaction_genes))



# ___________________________Venn diagrams______________________________________
make_venn_plot <- function(venn_data, colors, file_name) {
  
  # build a basic Venn diagram
  p <- ggvenn(
    venn_data, 
    fill_color = colors,       # Colors are passed as an argument
    stroke_size = 0.5,         # Circle line thickness
    set_name_size = 4,         # Font size for group names
    text_size = 4.5            # The size of the numbers inside the circles
  )
  
  # save the file
  ggsave(
    filename = file_name, 
    plot = p, 
    width = 8, 
    height = 8, 
    dpi = 300, 
    bg = "white"
  )
  
  return(p)
}

# extract significant genes
g_dis_100 <- rownames(subset(res_disease_100, padj < 0.05 & abs(log2FoldChange) > 1))
g_dis_50  <- rownames(subset(res_disease_50,  padj < 0.05 & abs(log2FoldChange) > 1))
g_diet_het <- rownames(subset(res_diet_het,   padj < 0.05 & abs(log2FoldChange) > 1))
g_diet_hom <- rownames(subset(res_diet_hom,   padj < 0.05 & abs(log2FoldChange) > 1))

# create a list for comparing genotypes
all_effects <- list(
  "Hom_vs_Het_100" = g_dis_100,
  "Hom_vs_Het_50"  = g_dis_50,
  "Diet_Effect_Het" = g_diet_het,
  "Diet_Effect_Hom" = g_diet_hom
)

# generate a matrix of all possible unique pairs
combinations <- t(combn(names(all_effects), 2))
pairs_matrix <- t(combn(names(all_effects), 2))

# Chart palette
venn_colors <- c("#1465B2", "#B31B21", "#2CA02C", "#FF7F0E", "#9467BD", "#8C564B")

for (i in 1:nrow(pairs_matrix)) {
  name1 <- pairs_matrix[i, 1]
  name2 <- pairs_matrix[i, 2]
  
  # Form a pair for the current chart
  current_pair <- list()
  current_pair[[name1]] <- all_effects[[name1]]
  current_pair[[name2]] <- all_effects[[name2]]
  
  # compose the name of the file and the title of the graph
  file_title <- paste0("Venn_", name1, "_vs_", name2, ".png")
  plot_title <- paste0("Venn: ", name1, " vs ", name2)
  
  # choose a couple of colors from the palette
  col_pair <- c(venn_colors[i], venn_colors[(i + 2) %% 6 + 1])
  
  # call the method
  make_venn_plot(current_pair, col_pair, file_title)
}



#======================Candidate genes==========================================
# pull out the leaders for the main effect of the disease on the diet
top_candidates <- as.data.frame(res_disease_50)

# filtering
top_candidates <- top_candidates[which(top_candidates$padj < 0.01 & abs(top_candidates$log2FoldChange) > 1.5), ]

# sorting
top_candidates <- top_candidates[order(top_candidates$padj), ]

# show the first ten
head(top_candidates, 10)
