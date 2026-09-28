library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)

monkeycount_data <- read_csv(".../monkeyTSC.com.data.csv")
monkey_TSC_tpmdata <- read_csv(".../monkeyTSC.com.tpm.data.csv")

monkeycount_data <- as.data.frame(monkeycount_data)
rownames(monkeycount_data) <- monkeycount_data[, 1]
monkeycount_data <- monkeycount_data[, -1]


group_info <- data.frame(
  sample = c("RES_2_TS_1", "RES_2_TS_2", "RES_2_TS_3", "RES_2_TS_4", "RES_2_TS_5",
             "RES_N1_TS_1", "RES_N1_TS_2", "RES_N1_TS_3", "RES_N1_TS_4", "RES_N1_TS_5",
             "RTS_1_1", "RTS_1_2", "RTS_2_1", "RTS_2_2", "RTS_3_1", "RTS_3_2", "RTS_4", "RTS_5"),
  group = c("M_P_dECM", "M_P_dECM", "m_Col_IV", "M_T_dECM", "M_T_dECM",
            "M_P_dECM", "M_P_dECM", "m_Col_IV", "M_T_dECM", "M_T_dECM",
            "M_T_dECM", "M_T_dECM", "M_T_dECM", "M_T_dECM", "m_Col_IV", "m_Col_IV", "M_P_dECM", "M_P_dECM")
)
print(group_info)

group_info$major_group <- ifelse(grepl("^RES", group_info$sample), "RES",
                                        ifelse(grepl("^RTS", group_info$sample), "RTS", NA))
group_info$major_group <- factor(group_info$major_group, levels = c("RES", "RTS"))
print(group_info)

dge <- DGEList(counts = monkeycount_data)

cpm_data <- cpm(dge)
keep <- rowSums(cpm_data > 1) >= 3  
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)
print(dge$samples$norm.factors) 

v <- voom(dge, plot = TRUE, normalize.method = "quantile")
expr_normalized <- v$E

print(dim(expr_normalized))


#### heatmap ####

group_info$new_group <- ifelse(
  str_detect(group_info$sample, "^RES"), 
  paste0("iTSCs_", group_info$group),
  paste0("TSCs_", group_info$group)
)

print(group_info)

new_groups <- unique(group_info$new_group)

new_expression_matrix <- matrix(nrow = nrow(monkey_TSC_tpmdata), 
                                ncol = length(new_groups))
colnames(new_expression_matrix) <- new_groups
rownames(new_expression_matrix) <- rownames(monkey_TSC_tpmdata)

for(grp in new_groups) {
  samples_in_group <- group_info$sample[group_info$new_group == grp]

  actual_colnames <- colnames(monkey_TSC_tpmdata)
  
  valid_samples <- samples_in_group[samples_in_group %in% actual_colnames]

  if(length(valid_samples) == 0) {
    valid_samples <- actual_colnames[sapply(samples_in_group, function(x) any(grepl(x, actual_colnames)))]
  }

  if(length(valid_samples) > 0) {
    if(length(valid_samples) > 1) {
      new_expression_matrix[, grp] <- rowMeans(monkey_TSC_tpmdata[, valid_samples, drop = FALSE])
    } else {
      new_expression_matrix[, grp] <- monkey_TSC_tpmdata[, valid_samples]
    }
    print(paste("success", grp, "n=:", length(valid_samples)))
  } else {
    print(paste("warn", grp, "no sample"))
  }
}

new_expression_matrix <- new_expression_matrix[, colSums(is.na(new_expression_matrix)) < nrow(new_expression_matrix)]

print(dim(new_expression_matrix))
print(colnames(new_expression_matrix))
print(head(new_expression_matrix))

genes_to_plot <- c("GATA3","GATA2","KRT7","TFAP2C","EPAS1","EGFR","KRT18","TP63","TEAD4","ELF5", "COL4A4", "LAMA2","ECM1",
                   "CDX2","CLDN4","NR2F2","TEAD1","DLX3","KLF5","GRHL1","ZFHX3", "CDH1", "EPCAM", #TSC
                   "SALL4","DPPA2","DNMT3B","POU5F1","TDGF1","DPPA4",  #ESC
                   "MMP2", "SP100", "FN1","ITGA5","IGF2","FSTL1","PAPPA", #EVT
                   "CGA","CYP19A1","SDC1","INHA","LHCGR","SP6",  #STB
                   "BNC1","CDH2", "EVC","LYPD6B","ISL1", "SALL2") #AMN 

rownames_new <- rownames(new_expression_matrix)

missing_genes <- genes_to_plot[!genes_to_plot %in% rownames_new]
print(missing_genes)

selected_genes <- monkey_TSC_tpmdata[genes_to_plot, ]  
selected_genes <- na.omit(selected_genes)


library(tidyverse)
library(ggplot2)
library(ggh4x) 
library(RColorBrewer)

selected_genes <- new_expression_matrix[genes_to_plot, ]  
selected_genes <- na.omit(selected_genes)

gene_groups <- data.frame(
  GeneClass = c(rep("TSC", 23), 
                rep("ESC", 6),
                rep("EVT", 7),      
                rep("STB", 6),  
                rep("AMN", 6))
)
rownames(gene_groups) <- rownames(selected_genes)

sample_annotation_df <- data.frame(
  Group = str_remove(colnames(selected_genes), "^(iTSCs|TSCs)_")
)
rownames(sample_annotation_df) <- colnames(selected_genes)

print(table(sample_annotation_df$Group))
print(sample_annotation_df)

heatmap_data <- selected_genes %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  mutate(Expression = log2(Expression + 1)) %>%
  left_join(gene_groups %>% rownames_to_column("Gene"), by = "Gene") %>%
  left_join(sample_annotation_df %>% rownames_to_column("Sample"), by = "Sample") %>%
  mutate(
    Gene = factor(Gene, levels = rev(genes_to_plot[genes_to_plot %in% unique(Gene)])),
    Sample = factor(Sample, levels = colnames(selected_genes)),
    GeneClass = factor(GeneClass, levels = c("TSC", "ESC", "EVT", "STB", "AMN")),
    Group = factor(Group, levels = c("M_P_dECM", "m_Col_IV", "M_T_dECM"))
  ) %>%
  group_by(Gene) %>%
  mutate(Scaled_Expression = scale(Expression)[,1]) %>%
  ungroup()

strips <- strip_themed(
  background_y = elem_list_rect(
    fill = c("#EC7A05","#9AC2CD99","#82B7F590","#5A77D190", "#FAD77B")  
  ),
  background_x = elem_list_rect(
    fill = c("#5A77D1", "#9AC2CD","#6B9B3C")  
  )
)


p1 <- ggplot(heatmap_data, aes(x = Sample, y = Gene, fill = Scaled_Expression)) +
  geom_tile(color = "grey80", linewidth = 0.2) +
  facet_grid2(GeneClass ~ Group, 
              scales = "free", space = "free",
              switch = "y", strip = strips) +
  scale_fill_gradientn(
    colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
    name = "Z-score",
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_x_discrete(expand = c(0, 0), position = "bottom") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Gene Expression Heatmap (All Genes)") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text.x = element_text(color = "white", size = 15, face = "bold"),
    strip.text.y = element_text(color = "white", size = 16, face = "bold"),
    strip.placement = "outside",
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(face = "bold"),
    panel.spacing = unit(0.1, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 30, 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)  # 标题加大并加粗
  )

print(p1)

heatmap_data_tsc <- selected_genes %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  left_join(gene_groups %>% rownames_to_column("Gene"), by = "Gene") %>%
  filter(GeneClass == "TSC") %>%
  pivot_longer(cols = -c(Gene, GeneClass), names_to = "Sample", values_to = "Expression") %>%
  # 在pivot_longer之后添加log2转换
  mutate(Expression = log2(Expression + 1)) %>%
  left_join(sample_annotation_df %>% rownames_to_column("Sample"), by = "Sample") %>%
  mutate(
    Gene = factor(Gene, levels = rev(genes_to_plot[genes_to_plot %in% unique(Gene)])),
    Sample = factor(Sample, levels = colnames(selected_genes)),
    GeneClass = factor(GeneClass, levels = "TSC"),
    Group = factor(Group, levels = c("M_P_dECM", "m_Col_IV", "M_T_dECM"))
  )

heatmap_data_tsc <- heatmap_data_tsc %>%
  group_by(Gene) %>%
  mutate(Scaled_Expression = scale(Expression)[,1]) %>%
  ungroup()

strips_tsc <- strip_themed(
  background_y = elem_list_rect(fill = "#EC7A05"),
  background_x = elem_list_rect(fill = c("#5A77D1", "#9AC2CD","#6B9B3C"))
)

p2 <- ggplot(heatmap_data_tsc, aes(x = Sample, y = Gene, fill = Scaled_Expression)) +
  geom_tile(color = "black", linewidth = 0.2) +
  facet_grid2(GeneClass ~ Group, 
              scales = "free", space = "free",
              switch = "y", strip = strips_tsc) +
  scale_fill_gradientn(
    colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
    name = "Z-score",
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_x_discrete(expand = c(0, 0), position = "bottom") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "TSC Genes Expression Heatmap") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text.x = element_text(color = "white", size = 14, face = "bold"),
    strip.text.y = element_text(color = "white", size = 15, face = "bold"),
    strip.placement = "outside",
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(face = "bold"),
    panel.spacing = unit(0.1, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 30, 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)  # 标题加大并加粗
  )

print(p2)

heatmap_data_no_tsc <- selected_genes %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  left_join(gene_groups %>% rownames_to_column("Gene"), by = "Gene") %>%
  filter(GeneClass != "TSC") %>%
  pivot_longer(cols = -c(Gene, GeneClass), names_to = "Sample", values_to = "Expression") %>%
  mutate(Expression = log2(Expression + 1)) %>%
  left_join(sample_annotation_df %>% rownames_to_column("Sample"), by = "Sample") %>%
  mutate(
    Gene = factor(Gene, levels = rev(genes_to_plot[genes_to_plot %in% unique(Gene)])),
    Sample = factor(Sample, levels = colnames(selected_genes)),
    GeneClass = factor(GeneClass, levels = c("ESC", "EVT", "STB", "AMN")),
    Group = factor(Group, levels = c("M_P_dECM", "m_Col_IV", "M_T_dECM"))
  )

heatmap_data_no_tsc <- heatmap_data_no_tsc %>%
  group_by(Gene) %>%
  mutate(Scaled_Expression = scale(Expression)[,1]) %>%
  ungroup()

strips_no_tsc <- strip_themed(
  background_y = elem_list_rect(
    fill = c("#9AC2CD99","#82B7F590","#5A77D190", "#FAD77B")
  ),
  background_x = elem_list_rect(
    fill = c("#5A77D1", "#9AC2CD","#6B9B3C")
  )
)

p3 <- ggplot(heatmap_data_no_tsc, aes(x = Sample, y = Gene, fill = Scaled_Expression)) +
  geom_tile(color = "black", linewidth = 0.2) +
  facet_grid2(GeneClass ~ Group, 
              scales = "free", space = "free",
              switch = "y", strip = strips_no_tsc) +
  scale_fill_gradientn(
    colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
    name = "Z-score",
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_x_discrete(expand = c(0, 0), position = "bottom") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Non-TSC Genes Expression Heatmap") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text.x = element_text(color = "white", size = 14, face = "bold"),
    strip.text.y = element_text(color = "white", size = 15, face = "bold"),
    strip.placement = "outside",
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(face = "bold"),
    panel.spacing = unit(0.1, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 30, 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)  # 标题加大并加粗
  )

print(p3)




##### violin #####

expr_matrix <- as.matrix(monkey_TSC_tpmdata)

target_genes <- c("KRT18", "GATA3", "NR2F2", "TFAP2C", "TEAD4", "CDX2", "FN1", "SP100", "ELF5", "COL4A4", "SDC1", "DNMT3B", "MMP2", "FSTL1", "ITGA5", "EPCAM")

target_genes <- c("CDH5", "CGA", "CGB7", "CGB8", "CYP19A1", "DLX3", "ERVW-1", "GATA2", "GATA3", "GCM1", "HSD3B1", "IGF2", "INSL4",
                  "KRT18", "KRT19", "KRT23", "KRT7", "KRT8", "LGALS13", "MMP14", "MUC15", "PSG2", "PSG3", "PSG4", "PSG5", "PSG6", 
                  "S100P", "SLC22A11-1", "SLC22A11-2", "TBX3", "TFAP2C", "VGLL1", "XAGE3")


df <- data.frame(
  sample = c("RTS_1_1", "RTS_1_2", "RTS_2_1", "RTS_2_2", "RTS_3_1", "RTS_3_2",  "RTS_4", "RTS_5"),
  group = c("M_T_dECM", "M_T_dECM", "M_T_dECM", "M_T_dECM", "m_Col_IV", "m_Col_IV", "M_P_dECM", "M_P_dECM")
)

group_info <- df

group_order <- c("M_P_dECM", "m_Col_IV")

group_colors <- c("M_P_dECM" = "#e14c33", "m_Col_IV" = "#4b89bf")

target_expr <- expr_matrix[rownames(expr_matrix) %in% target_genes, ]

library(reshape2)
expr_melt <- melt(target_expr)
colnames(expr_melt) <- c("Gene", "Sample", "Expression")

expr_melt <- expr_melt %>%
  left_join(group_info, by = c("Sample" = "sample"))

expr_melt <- na.omit(expr_melt)
expr_melt <- expr_melt[expr_melt$group %in% group_order, ]

expr_melt$group <- factor(expr_melt$group, levels = group_order)
expr_melt$Gene <- factor(expr_melt$Gene, levels = target_genes)

library(ggpubr)

p <- ggplot(expr_melt, aes(x = group, y = log2(Expression + 1), fill = group)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = 1) +  
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5, color = "black", linewidth = 1) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("M_P_dECM", "m_Col_IV")),
    label = "p.format",  
    label.sep = "=",
    size = 5,
    vjust = -0.5,
    label.y.npc = 0.95,
    color = "black",  
    bracket.size = 1  
  ) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("M_P_dECM", "m_Col_IV")),
    label = "p.signif",  
    size = 5,
    vjust = 1.5,
    label.y.npc = 0.95,
    color = "black",  
    bracket.size = 1  
  ) +
  facet_wrap(~ Gene, ncol = 5) +  
  labs(title = "Expression Distribution of Key Genes",
       x = "Group",
       y = "log2(TPM + 1)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 10, color = "black"),  
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),  
    axis.title.x = element_text(face = "bold", size = 14, color = "black"),  
    axis.title.y = element_text(face = "bold", size = 14, color = "black"),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, color = "black"),  
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12, color = "black"),  
    panel.grid = element_blank(),  
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 1, color = "black"),  
    # 显著性标签加粗
    text = element_text(face = "bold", color = "black") 
  ) +
  scale_fill_manual(values = group_colors)

print(p)

individual_plots <- list()

for (gene in target_genes) {
  gene_data <- expr_melt[expr_melt$Gene == gene, ]
  p_single <- ggplot(gene_data, aes(x = group, y = log2(Expression + 1), fill = group)) +
    geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = 0.5) +  
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5, color = "black", linewidth = 0.5) +  
    stat_compare_means(
      method = "t.test",
      comparisons = list(c("M_P_dECM", "m_Col_IV")),
      label = "p.format",  
      label.sep = "=",
      size = 6,
      vjust = -0.5,
      label.y.npc = 0.95,
      color = "black",  
      bracket.size = 0.5  
    ) +
    stat_compare_means(
      method = "t.test",
      comparisons = list(c("M_P_dECM", "m_Col_IV")),
      label = "p.signif",  
      size = 6,
      vjust = 1.5,
      label.y.npc = 0.95,
      color = "black",  
      bracket.size = 0.5  
    ) +
    labs(title = paste("Expression of", gene),
         x = "",
         y = "log2(TPM + 1)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 17, face = "bold", color = "black"),  
      axis.text.y = element_text(size = 15, face = "bold", color = "black"),  
      axis.title.x = element_text(size = 16, face = "bold", color = "black"),  
      axis.title.y = element_text(size = 16, face = "bold", color = "black"),  
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, color = "black"),  
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.5, color = "black"),  
      text = element_text(face = "bold", color = "black")  
    ) +
    scale_fill_manual(values = group_colors)
  
  individual_plots[[gene]] <- p_single

  ggsave(paste0("RTS_genes_t_group_", gene, ".pdf"), 
         plot = p_single, 
         width = 4, 
         height = 6, 
         dpi = 300)
}



#### vocanol and GO ####
expr_matrix <- as.matrix(expr_normalized)

group_info <- data.frame(
  sample = c("RES_2_TS_1", "RES_2_TS_2", "RES_2_TS_3", "RES_2_TS_4", "RES_2_TS_5",
             "RES_N1_TS_1", "RES_N1_TS_2", "RES_N1_TS_3", "RES_N1_TS_4", "RES_N1_TS_5",
             "RTS_1_1", "RTS_1_2", "RTS_2_1", "RTS_2_2", "RTS_3_1", "RTS_3_2", "RTS_4", "RTS_5"),
  group = c("RES", "RES", "RES", "RES", "RES",
            "RES", "RES", "RES", "RES", "RES",
            "RTS", "RTS", "RTS", "RTS", "RTS", "RTS", "RTS", "RTS")
)

print(group_info)

group <- factor(group_info$group, levels = c("RES", "RTS"))
design <- model.matrix(~0 + group)
colnames(design) <- c("RES", "RTS")

fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(RES - RTS, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(fit2, number = Inf, adjust.method = "BH")
deg_results$Gene <- rownames(deg_results)

deg_results <- deg_results %>%
  mutate(
    log2FC = logFC,
    neg_log10_pval = -log10(P.Value),
    padj = adj.P.Val
  )

fc_threshold <- 0.5  # log2FC
pval_threshold <- 0.05  # P

deg_results <- deg_results %>%
  mutate(
    significance = case_when(
      log2FC >= fc_threshold & P.Value < pval_threshold ~ "Up",
      log2FC <= -fc_threshold & P.Value < pval_threshold ~ "Down",
      TRUE ~ "Not significant"
    )
  )

library(dplyr)

label_data <- deg_results %>%
  filter(significance %in% c("Up", "Down")) %>%
  group_by(significance) %>%
  arrange(P.Value) %>%
  slice_head(n = 5) %>%
  ungroup()


focus_genes <- c("DNMT3B", "TEDA4", "GATA3", "TFAP2C", "KRT18", "SDC1", "FN1", "ELF5")


volcano_plot <- ggplot(deg_results, aes(x = log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = significance, fill = significance), 
             alpha = 0.7, size = 2, shape = 21) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
             linetype = "dashed", color = "black", linewidth = 0.3, alpha = 0.6) +
  geom_hline(yintercept = -log10(pval_threshold), 
             linetype = "dashed", color = "black", linewidth = 0.3, alpha = 0.6) +
  scale_color_manual(
    name = "",
    values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Not significant" = "gray60"),
    labels = c("Down-regulated", "Not significant", "Up-regulated")
  ) +
  scale_fill_manual(
    name = "",
    values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Not significant" = "gray60"),
    labels = c("Down-regulated", "Not significant", "Up-regulated")
  ) +
  labs(
    title = "iTSCs vs TSCs",
    x = expression(Log[2] * " Fold Change"),
    y = expression(-Log[10] * " P-value"),
    caption = paste0("Fold change threshold: ±", fc_threshold, 
                     "; P-value threshold: ", pval_threshold)
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, 
                              margin = margin(b = 15)),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12, face = "bold"),
    legend.key.size = unit(0.8, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 5, l = 0),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    plot.caption = element_text(size = 14, color = "gray40", hjust = 0, face = "bold")
  ) +
  coord_cartesian(
    xlim = c(-max(abs(deg_results$log2FC)) * 1.05, max(abs(deg_results$log2FC)) * 1.05)
  ) +
  geom_text_repel(
    data = label_data,
    aes(x = log2FC, y = neg_log10_pval, label = Gene),
    size = 4.5,
    fontface = "bold",
    family = "Arial", 
    box.padding = 1.5,
    point.padding = 0.2,
    max.overlaps = 8,
    segment.color = "gray40",
    segment.size = 0.3,
    min.segment.length = 0.3,
    force = 3,
    force_pull = 3,
    max.time = 2,
    max.iter = 10000,
    direction = "both",
    seed = 123,
    inherit.aes = FALSE
  ) +
  
  geom_label_repel(
    data = subset(deg_results, Gene %in% focus_genes),
    aes(x = log2FC, y = neg_log10_pval, label = Gene),
    size = 5,
    fontface = "bold.italic",
    family = "Arial",  
    color = "#984EA3",
    fill = "white",
    alpha = 0.8,
    box.padding = 1.2,
    point.padding = 0.3,
    max.overlaps = 10,
    segment.color = "grey40",
    segment.size = 0.4,
    segment.alpha = 1,
    min.segment.length = 0.3,
    force = 3,
    force_pull = 2,
    max.time = 3,
    max.iter = 20000,
    direction = "both",
    nudge_x = 0.3,
    nudge_y = 0.3,
    seed = 456,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = subset(deg_results, Gene %in% focus_genes),
    aes(x = log2FC, y = neg_log10_pval),
    shape = 21, size = 3, fill = "#984EA3", color = "black", stroke = 0.5,
    inherit.aes = FALSE
  )

print(volcano_plot)



library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)

up_genes <- deg_results %>%
  filter(significance == "Up" & adj.P.Val < 0.05 & abs(log2FC) > fc_threshold) %>%
  pull(Gene)

down_genes <- deg_results %>%
  filter(significance == "Down" & adj.P.Val < 0.05 & abs(log2FC) > fc_threshold) %>%
  pull(Gene)

print(paste("up:", length(up_genes)))
print(paste("down:", length(down_genes)))

if(length(up_genes) >= 10) {
  go_up <- enrichGO(
    gene = up_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",  #  "BP", "MF", "CC" 或 "ALL"
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
  )
  
  kegg_up <- enrichKEGG(
    gene = bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )
} else {
  warning("no genes")
}

if(length(down_genes) >= 10) {
  go_down <- enrichGO(
    gene = down_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
  )
  
  kegg_down <- enrichKEGG(
    gene = bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )
} else {
  warning("no genes")
}

if(exists("go_up") && nrow(go_up) > 0) {
  p_go_up <- dotplot(go_up, showCategory = 12, 
                     title = "GO Enrichment - Up-regulated Genes") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold")
    )
  print(p_go_up)
  ggsave("GO_upregulated_dotplot.pdf", plot = p_go_up, width = 9, height = 9)
}

if(exists("kegg_up") && nrow(kegg_up) > 0) {
  p_kegg_up <- dotplot(kegg_up, showCategory = 12, 
                       title = "KEGG Enrichment - Up-regulated Genes") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold")
    )
  
  print(p_kegg_up)
  ggsave("KEGG_upregulated_dotplot.pdf", plot = p_kegg_up, width = 9, height = 9)
}

if(exists("go_down") && nrow(go_down) > 0) {
  p_go_down <- dotplot(go_down, showCategory = 12, 
                       title = "GO Enrichment - Down-regulated Genes") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold")
    )
  print(p_go_down)
  ggsave("GO_downregulated_dotplot.pdf", plot = p_go_down, width = 9, height = 9)
}


if(exists("kegg_down") && nrow(kegg_down) > 0) {
  p_kegg_down <- dotplot(kegg_down, showCategory = 12, 
                         title = "KEGG Enrichment - Down-regulated Genes") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold")
    )
  
  print(p_kegg_down)
  ggsave("KEGG_downregulated_dotplot.pdf", plot = p_kegg_down, width = 9, height = 9)
}

if(!is.null(go_up) && !is.null(go_down) && nrow(go_up) > 0 && nrow(go_down) > 0) {

  common_cols <- intersect(colnames(go_up@result), colnames(go_down@result))

  important_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", 
                      "p.adjust", "qvalue", "geneID", "Count")

  available_cols <- intersect(important_cols, common_cols)

  go_up_df <- go_up@result[, available_cols] %>% 
    mutate(Direction = "Up") %>%
    head(10)
  
  go_down_df <- go_down@result[, available_cols] %>% 
    mutate(Direction = "Down") %>%
    head(10)

  go_combined <- rbind(go_up_df, go_down_df)

  p_combined <- ggplot(go_combined, aes(x = Direction, y = Description, 
                                        size = Count, color = -log10(pvalue))) +
    geom_point(alpha = 0.8) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(P-value)") +
    scale_size_continuous(name = "Gene Count", range = c(3, 8)) +
    labs(
      title = "GO Enrichment Comparison: Up vs Down-regulated Genes",
      x = "Regulation Direction",
      y = "GO Terms"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold")
    )
  
  print(p_combined)
  ggsave("GO_comparison_up_down.pdf", plot = p_combined, width = 12, height = 10)
  
} else {
  message("no result")
}

if(exists("go_up")) write.csv(go_up@result, "GO_upregulated_results.csv", row.names = FALSE)
if(exists("go_down")) write.csv(go_down@result, "GO_downregulated_results.csv", row.names = FALSE)
if(exists("kegg_up")) write.csv(kegg_up@result, "KEGG_upregulated_results.csv", row.names = FALSE)
if(exists("kegg_down")) write.csv(kegg_down@result, "KEGG_downregulated_results.csv", row.names = FALSE)



enrich_BP<-enrichGO(down_genes, org.Hs.eg.db, keyType="SYMBOL",ont ="BP", pvalueCutoff=0.05, pAdjustMethod ="BH", qvalueCutoff =0.05)
enrich_MF<-enrichGO(down_genes, org.Hs.eg.db, keyType="SYMBOL",ont ="MF", pvalueCutoff=0.05, pAdjustMethod ="BH", qvalueCutoff =0.05)
enrich_CC<-enrichGO(down_genes, org.Hs.eg.db, keyType="SYMBOL",ont ="CC", pvalueCutoff=0.05, pAdjustMethod ="BH", qvalueCutoff =0.05)
library(enrichplot)
bp_df<- as.data.frame(enrich_BP)
mf_df<- as.data.frame(enrich_MF)
cc_df<- as.data.frame(enrich_CC)
bp_df$Category <-"Biological Process"  
mf_df$Category <-"Molecular Function"  
cc_df$Category <-"Cellular Component"  
combined_df<- rbind(bp_df, mf_df, cc_df)
combined_df<- combined_df[combined_df$p.adjust <0.05, ]

library(ggplot2)
library(dplyr)
library(stringr)
combined_top <- combined_df %>% 
  group_by(Category) %>%   
  arrange(p.adjust) %>% 
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = str_wrap(Description, width = 60))

library(ggpubr)

combined_top$GeneRatio <- sapply(strsplit(combined_top$GeneRatio, "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})

ggplot(combined_top, aes(x = reorder(Description, -log10(p.adjust)), 
                         y = -log10(p.adjust), 
                         fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  facet_grid(Category ~ ., scales = "free_y", space = "free") +
  theme_pubr() + 
  labs(x = "GO Term", y = "-log10(Adjusted p-value)") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    strip.text = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    text = element_text(face = "bold")
  )





#### venn ####
library(ggVennDiagram)
library(ggplot2)

up_genes <- deg_results %>% filter(significance == "Up") %>% pull(Gene)
down_genes <- deg_results %>% filter(significance == "Down") %>% pull(Gene)
not_sig_genes <- deg_results %>% filter(significance == "Not significant") %>% pull(Gene)

rna_list <- list(
  "Up-regulated" = c(up_genes, not_sig_genes),
  "Down-regulated" = c(down_genes, not_sig_genes)
)


library(ggvenn)

ggvenn(rna_list, 
       columns = c("Up-regulated", "Down-regulated"),
       show_elements = F,
       fill_color = c( "#f06152","#76afda" ),  
       fill_alpha = 0.6,
       stroke_color = "black",
       stroke_size = 0,
       set_name_color = "black",
       set_name_size = 6,
       text_color = "black",
       text_size = 6,
       auto_scale = T
) +
  labs(title = "Gene Expression Patterns Venn Diagram") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))


p <- ggVennDiagram(rna_list, 
                   category.names = c("Up-regulated", "Down-regulated"),
                   set_color = "black",
                   set_size = 5,
                   label = "both",  # "count","percent","both","none"
                   label_alpha = 0,  
                   label_color = "black", 
                   label_size = 5,  
                   label_font = "serif") +
  coord_flip() +
  
  scale_fill_gradient(low = "#76afda", high = "#f06152") +
  
  labs(title = "Gene Expression Patterns Venn Diagram") +
  theme(legend.position = "none",
        text = element_text(family = "serif"),
        plot.title = element_text(hjust = 0.5, size = 15))

print(p)
