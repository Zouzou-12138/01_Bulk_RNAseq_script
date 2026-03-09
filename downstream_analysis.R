suppressMessages({
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)
library(venn)   
library(dplyr)
library(tools)
})

# 创建必要的输出目录
dir.create("./07_AdvancedQC", showWarnings = FALSE)
dir.create("./08_DE", showWarnings = FALSE)
dir.create("./12_PCA", showWarnings = FALSE)
dir.create("./09_GO", showWarnings = FALSE)
dir.create("./10_KEGG", showWarnings = FALSE)  # 新增KEGG目录

# 1. 读取文件
cat("读取表达矩阵...\n")
gene_counts <- read_csv("/media/desk16/iy11913/master_project/01_Bulk_RNAseq/06_GeneExpression/gene_expression_matrix.csv")
gene_TPM <- read_csv("/media/desk16/iy11913/master_project/01_Bulk_RNAseq/06_GeneExpression/gene_expression_TPM.csv")

# 2. 读取分组信息（外部文件）
cat("读取分组信息...\n")
if(!file.exists("sample_group.csv")) {
  # 如果文件不存在，创建一个示例文件
  sample_group_df <- data.frame(
    sample = c("WT-testis-1", "WT-testis-2", "C2-KO-testis-1", "C2-KO-testis-2",
               "C6-KO-testis-1", "C6-KO-testis-2", "C26-DKO-testis-1", "C26-DKO-testis-2"),
    group = c("WT", "WT", "C2_KO", "C2_KO", "C6_KO", "C6_KO", "C26_DKO", "C26_DKO")
  )
  write.csv(sample_group_df, "sample_group.csv", row.names = FALSE)
  cat("已创建示例分组文件: sample_group.csv\n")
}

sample_groups <- read.csv("sample_group.csv")
group_vector <- setNames(sample_groups$group, sample_groups$sample)
group_levels <- unique(sample_groups$group)
cat("分组信息:", paste(group_levels, collapse = ", "), "\n")

# 3. 数据过滤
cat("过滤低表达基因...\n")

# 自动识别样本列（所有数值列）
sample_cols <- which(sapply(gene_counts, is.numeric))
cat("检测到", length(sample_cols), "个数值列（样本）\n")

if(length(sample_cols) == 0) {
  stop("没有找到数值列，请检查数据格式")
}

# 过滤基因
gene_counts <- gene_counts[apply(gene_counts[, sample_cols] > 10, 1, all),]
gene_TPM <- gene_TPM[apply(gene_TPM[, sample_cols] > 1, 1, all), ]

cat("过滤后剩余基因数:", nrow(gene_counts), "\n")

# 4. 计算相关性热图
cat("绘制相关性热图...\n")
exprSet <- as.matrix(as.data.frame(gene_TPM)[,-1])
rownames(exprSet) <- as.data.frame(gene_TPM)[,1]

# 筛选top800变异基因并计算相关性
gene_mad <- apply(exprSet, 1, mad)
top_genes <- names(sort(gene_mad, decreasing = TRUE)[1:min(800, nrow(exprSet))])
exprSet_subset <- exprSet[top_genes, , drop = FALSE]
M <- cor(exprSet_subset, method = "spearman")

# 注释信息（使用外部分组）
group <- group_vector[colnames(exprSet_subset)]
anno <- data.frame(Group = group, row.names = colnames(exprSet_subset))

# 设置分组颜色
group_colors <- c("#4DBBD5", "#E64B35", "#E64B35", "#E64B35")
names(group_colors) <- group_levels
ann_colors <- list(Group = group_colors)

# 绘制并保存相关性热图
p_cor <- pheatmap(M, display_numbers = TRUE, annotation_col = anno, 
                  annotation_colors = ann_colors, cellheight = 40, cellwidth = 40,
                  color = colorRampPalette(c("blue", "white", "red"))(100),
                  main = "Correlation Analysis")
ggsave("./07_AdvancedQC/correlation_heatmap.png", p_cor, width = 10, height = 8, dpi = 300)
print(p_cor)

# 5. 基因表达热图
cat("绘制基因表达热图...\n")
exprSet_raw <- as.matrix(as.data.frame(gene_TPM)[,-1])
rownames(exprSet_raw) <- as.data.frame(gene_TPM)[,1]

# 筛选top2000变异基因并log转换
gene_mad <- apply(exprSet_raw, 1, mad)
top_genes_exp <- names(sort(gene_mad, decreasing = TRUE)[1:min(2000, nrow(exprSet_raw))])
exprSet_log <- log2(exprSet_raw[top_genes_exp, , drop = FALSE] + 1)

# 绘制并保存表达热图
p_exp <- pheatmap(exprSet_log, scale = "row", show_rownames = FALSE,
                  annotation_col = anno, annotation_colors = ann_colors,
                  color = colorRampPalette(c("navy", "white", "firebrick3"))(200),
                  main = "Gene Expression Heatmap (Top 2000 Most Variable Genes)")
ggsave("./07_AdvancedQC/gene_expression_heatmap.png", p_exp, width = 12, height = 10, dpi = 300)
print(p_exp)

# 6. PCA分析
cat("进行PCA分析...\n")
expr_pca <- as.data.frame(gene_TPM)
rownames(expr_pca) <- expr_pca[,1]
expr_pca <- expr_pca[,-1]

expr_log <- log2(expr_pca + 1)
pca_result <- prcomp(t(expr_log), center = TRUE, scale. = TRUE)

pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data$condition <- group_vector[rownames(pca_data)]
pca_data$name <- rownames(pca_data)

percentVar <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) + 
  geom_text_repel(size = 3, max.overlaps = 20) + 
  labs(x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%"),
       title = "PCA Analysis (log2 TPM)") +
  scale_color_manual(values = group_colors) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("12_PCA/PCA_plot.png", p_pca, width = 8, height = 6, dpi = 300)
print(p_pca)

# 7. 差异表达分析
cat("进行差异表达分析...\n")
suppressMessages({
  gene_info <- read_csv("./mouse_gene_info.csv", na = "NA")
})

# 添加gene_name
gene_counts <- gene_counts %>%
  left_join(gene_info[c("GeneID", "Gene_name")], by = c("gene_id" = "GeneID")) %>%
  relocate(Gene_name, .after = 1)

gene_counts_df <- as.data.frame(gene_counts)

# 提取counts数据
counts_data <- gene_counts_df[, 3:10]
rownames(counts_data) <- gene_counts_df$gene_id

# 使用外部分组信息
sample_info <- data.frame(
  group = factor(group_vector[colnames(counts_data)], levels = group_levels)
)
rownames(sample_info) <- colnames(counts_data)

cat("样本分组:\n")
print(sample_info)

# DESeq2分析
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = sample_info, design = ~ group)
dds <- DESeq(dds)

# 自动生成所有非WT vs WT的比较
comparisons <- list()
for(g in group_levels[group_levels != "WT"]) {
  comparisons[[g]] <- c("group", g, "WT")
}

# 对每个比较进行分析
for(comp in names(comparisons)) {
  cat("分析比较:", comp, "vs WT\n")
  
  # 获取结果
  res <- results(dds, contrast = comparisons[[comp]], alpha = 0.05)
  res_df <- as.data.frame(res)
  
  # 添加基因信息
  res_df$Gene_id <- rownames(res_df)
  res_df$Gene_name <- gene_counts_df$Gene_name[match(rownames(res_df), gene_counts_df$gene_id)]
  res_df$regulated <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                          ifelse(res_df$log2FoldChange > 1, "up", "down"), "none")
  
  # 保存DEG列表
  deg <- res_df %>% filter(regulated != "none") %>% arrange(regulated, padj)
  write.csv(deg, file = paste0("./08_DE/DEG_", comp, "_vs_WT.csv"), row.names = FALSE)
  
  # 统计DEG数量
  up_count <- sum(res_df$regulated == "up", na.rm = TRUE)
  down_count <- sum(res_df$regulated == "down", na.rm = TRUE)
  cat("  上调基因:", up_count, " 下调基因:", down_count, "\n")
  
  # 绘制火山图
  plot_data <- res_df %>%
    mutate(diff = case_when(padj < 0.05 & log2FoldChange > 1 ~ "Up",
                             padj < 0.05 & log2FoldChange < -1 ~ "Down",
                             TRUE ~ "Not significant"),
           diff = factor(diff, levels = c("Down", "Not significant", "Up")))
  
  up <- sum(plot_data$diff == "Up", na.rm = TRUE)
  down <- sum(plot_data$diff == "Down", na.rm = TRUE)
  
  top_genes <- plot_data %>% 
    filter(diff != "Not significant") %>% 
    arrange(padj) %>% 
    group_by(diff) %>% 
    slice_head(n = 5)
  
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj), color = diff)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_color_manual(values = c("Down" = "#377EB8", "Not significant" = "grey60", "Up" = "#E41A1C"),
                       labels = c(paste0("Down (", down, ")"), "Not significant", paste0("Up (", up, ")"))) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = top_genes, aes(label = Gene_name), size = 3.5, max.overlaps = Inf) +
    labs(title = paste0(comp, " vs WT"), 
         x = "log2 Fold Change",
         y = "-log10(adjusted p-value)") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(filename = paste0("./08_DE/volcano_", comp, "_vs_WT.png"), plot = p, width = 8, height = 6, dpi = 300)
  print(p)
}

# 8. GO富集分析
cat("\n开始GO富集分析...\n")

# 获取所有差异表达文件
deg_files <- list.files("./08_DE", pattern = "\\.csv$", full.names = TRUE)

# 遍历每个文件
for(file in deg_files) {
  
  # 提取文件名
  file_name <- tools::file_path_sans_ext(basename(file))
  cat("\n处理文件:", file_name, "\n")
  
  # 创建该文件的输出子目录
  output_dir <- paste0("./09_GO/", file_name)
  dir.create(output_dir, showWarnings = FALSE)
  
  # 读取差异分析结果
  df <- read.csv(file)
  
  # 提取上调和下调的基因
  upregulated_genes <- df$Gene_id[df$regulated == "up"]
  downregulated_genes <- df$Gene_id[df$regulated == "down"]
  
  cat("上调基因数量:", length(upregulated_genes), "\n")
  cat("下调基因数量:", length(downregulated_genes), "\n")
  
  # 定义函数进行GO分析
  run_go_analysis <- function(gene_list, direction, file_prefix) {
    
    if(length(gene_list) == 0) {
      cat(direction, ": 没有基因\n")
      return(NULL)
    }
    
    # 进行GO富集分析（三个本体）
    for(ont in c("BP", "CC", "MF")) {
      cat("  分析", ont, "...\n")
      
      ego <- tryCatch({
        enrichGO(gene = gene_list, 
                 OrgDb = org.Mm.eg.db, 
                 keyType = 'ENSEMBL', 
                 ont = ont, 
                 pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05,
                 readable = TRUE)
      }, error = function(e) {
        cat("    GO分析出错:", e$message, "\n")
        return(NULL)
      })
      
      # 保存结果到文件
      if(!is.null(ego) && nrow(ego@result) > 0) {
        # 保存完整结果
        write.csv(as.data.frame(ego@result), 
                  file = paste0(output_dir, "/", file_prefix, "_", direction, "_", ont, "_full.csv"), 
                  row.names = FALSE)
        
        # 保存显著结果
        sig_result <- as.data.frame(ego@result) %>% filter(p.adjust < 0.05)
        if(nrow(sig_result) > 0) {
          write.csv(sig_result, 
                    file = paste0(output_dir, "/", file_prefix, "_", direction, "_", ont, "_significant.csv"), 
                    row.names = FALSE)
          
          # 绘制dotplot
          p <- dotplot(ego, showCategory = 15) + 
            ggtitle(paste0(file_prefix, " - ", direction, " - ", ont)) +
            theme(plot.title = element_text(hjust = 0.5, size = 12))
          
          ggsave(filename = paste0(output_dir, "/", file_prefix, "_", direction, "_", ont, "_dotplot.png"), 
                 plot = p, width = 10, height = 8, dpi = 300)
          
          # 绘制barplot
          p2 <- barplot(ego, showCategory = 15) + 
            ggtitle(paste0(file_prefix, " - ", direction, " - ", ont)) +
            theme(plot.title = element_text(hjust = 0.5, size = 12))
          
          ggsave(filename = paste0(output_dir, "/", file_prefix, "_", direction, "_", ont, "_barplot.png"), 
                 plot = p2, width = 10, height = 8, dpi = 300)
          
          cat("    ", direction, "-", ont, ": 找到", nrow(sig_result), "个显著GO terms\n")
        } else {
          cat("    ", direction, "-", ont, ": 没有显著GO terms\n")
        }
      }
    }
  }
  
  # 对上调和下调基因分别进行GO分析
  run_go_analysis(upregulated_genes, "up", file_name)
  run_go_analysis(downregulated_genes, "down", file_name)
}

cat("\n所有GO分析完成！结果保存在 ./09_GO 目录\n")

# 9. KEGG富集分析（保存到./10_KEGG）
cat("\n开始KEGG富集分析...\n")

# 重新获取差异表达文件
deg_files <- list.files("./08_DE", pattern = "\\.csv$", full.names = TRUE)

# 遍历每个文件
for(file in deg_files) {
  
  # 提取文件名
  file_name <- tools::file_path_sans_ext(basename(file))
  cat("\n处理文件:", file_name, "\n")
  
  # 创建KEGG输出子目录
  output_dir <- paste0("./10_KEGG/", file_name)
  dir.create(output_dir, showWarnings = FALSE)
  
  # 读取差异分析结果
  df <- read.csv(file)
  
  # 提取上调和下调的基因
  upregulated_genes <- df$Gene_id[df$regulated == "up"]
  downregulated_genes <- df$Gene_id[df$regulated == "down"]
  
  cat("上调基因数量:", length(upregulated_genes), "\n")
  cat("下调基因数量:", length(downregulated_genes), "\n")
  
  # 定义函数进行KEGG分析
  run_kegg_analysis <- function(gene_list, direction, file_prefix) {
    
    if(length(gene_list) == 0) {
      cat(direction, ": 没有基因\n")
      return(NULL)
    }
    
    cat("  转换基因ID...\n")
    # 将ENSEMBL ID转换为ENTREZ ID
    gene_entrez <- tryCatch({
      bitr(gene_list, 
           fromType = "ENSEMBL", 
           toType = "ENTREZID", 
           OrgDb = org.Mm.eg.db)
    }, error = function(e) {
      cat("    ID转换出错:", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(gene_entrez) || nrow(gene_entrez) == 0) {
      cat(direction, ": 没有成功转换的基因ID\n")
      return(NULL)
    }
    
    cat("  进行KEGG富集分析...\n")
    # 进行KEGG富集分析
    kk <- tryCatch({
      enrichKEGG(gene = gene_entrez$ENTREZID,
                 organism = 'mmu',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
    }, error = function(e) {
      cat("    KEGG分析出错:", e$message, "\n")
      return(NULL)
    })
    
    # 保存结果到文件
    if(!is.null(kk) && nrow(kk@result) > 0) {
      # 保存完整结果
      write.csv(as.data.frame(kk@result), 
                file = paste0(output_dir, "/", file_prefix, "_", direction, "_KEGG_full.csv"), 
                row.names = FALSE)
      
      # 保存显著结果
      sig_result <- as.data.frame(kk@result) %>% filter(p.adjust < 0.05)
      if(nrow(sig_result) > 0) {
        write.csv(sig_result, 
                  file = paste0(output_dir, "/", file_prefix, "_", direction, "_KEGG_significant.csv"), 
                  row.names = FALSE)
        
        # 绘制dotplot
        p <- dotplot(kk, showCategory = 15) + 
          ggtitle(paste0(file_prefix, " - ", direction, " - KEGG Pathway")) +
          theme(plot.title = element_text(hjust = 0.5, size = 12))
        
        ggsave(filename = paste0(output_dir, "/", file_prefix, "_", direction, "_KEGG_dotplot.png"), 
               plot = p, width = 10, height = 8, dpi = 300)
        
        # 绘制barplot
        p2 <- barplot(kk, showCategory = 15) + 
          ggtitle(paste0(file_prefix, " - ", direction, " - KEGG Pathway")) +
          theme(plot.title = element_text(hjust = 0.5, size = 12))
        
        ggsave(filename = paste0(output_dir, "/", file_prefix, "_", direction, "_KEGG_barplot.png"), 
               plot = p2, width = 10, height = 8, dpi = 300)
        
        # 准备foldChange用于网络图
        fc_vector <- setNames(df$log2FoldChange[match(gene_entrez$ENSEMBL, df$Gene_id)], 
                             gene_entrez$ENTREZID)
        
        # 绘制cnetplot
        tryCatch({
          p3 <- cnetplot(kk, showCategory = 10, foldChange = fc_vector)
          ggsave(filename = paste0(output_dir, "/", file_prefix, "_", direction, "_KEGG_cnetplot.png"), 
                 plot = p3, width = 12, height = 10, dpi = 300)
        }, error = function(e) {
          cat("    cnetplot绘制失败:", e$message, "\n")
        })
        
        # 绘制heatplot
        tryCatch({
          p4 <- heatplot(kk, showCategory = 10, foldChange = fc_vector)
          ggsave(filename = paste0(output_dir, "/", file_prefix, "_", direction, "_KEGG_heatplot.png"), 
                 plot = p4, width = 12, height = 10, dpi = 300)
        }, error = function(e) {
          cat("    heatplot绘制失败:", e$message, "\n")
        })
        
        # 绘制emapplot
        tryCatch({
          kk_pair <- pairwise_termsim(kk)
          p5 <- emapplot(kk_pair, showCategory = 10)
          ggsave(filename = paste0(output_dir, "/", file_prefix, "_", direction, "_KEGG_emapplot.png"), 
                 plot = p5, width = 12, height = 10, dpi = 300)
        }, error = function(e) {
          cat("    emapplot绘制失败:", e$message, "\n")
        })
        
        # 绘制pathway图（如果安装了pathview包）
        if(requireNamespace("pathview", quietly = TRUE)) {
          tryCatch({
            # 对每个显著通路绘制pathway图
            for(i in 1:min(5, nrow(sig_result))) {
              pathway_id <- sig_result$ID[i]
              pathway_name <- gsub("/", "_", sig_result$Description[i])
              pathview::pathview(gene.data = fc_vector,
                                 pathway.id = pathway_id,
                                 species = "mmu",
                                 out.suffix = paste0(file_prefix, "_", direction, "_", pathway_name),
                                 kegg.native = TRUE)
              # 移动生成的文件到输出目录
              files_to_move <- list.files(pattern = paste0(pathway_id, ".*.", "png"))
              for(f in files_to_move) {
                file.rename(f, file.path(output_dir, f))
              }
            }
          }, error = function(e) {
            cat("    pathview绘制失败:", e$message, "\n")
          })
        }
        
        cat("    ", direction, "- KEGG: 找到", nrow(sig_result), "个显著通路\n")
      } else {
        cat("    ", direction, "- KEGG: 没有显著通路\n")
      }
    } else {
      cat("    ", direction, "- KEGG: 没有富集结果\n")
    }
  }
  
  # 对上调和下调基因分别进行KEGG分析
  run_kegg_analysis(upregulated_genes, "up", file_name)
  run_kegg_analysis(downregulated_genes, "down", file_name)
}

cat("\n所有KEGG分析完成！结果保存在 ./10_KEGG 目录\n")
cat("\n全部分析完成！\n")