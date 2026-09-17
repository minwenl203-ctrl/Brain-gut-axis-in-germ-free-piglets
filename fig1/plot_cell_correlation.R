library(Seurat)
library(corrplot)
library(pheatmap)
library(dplyr)

# 1. 加载数据
brain_ileum <- readRDS("06_manual_annotate/Brain_new.rds")
output_dir      <- "./07_plot/"
# 2. 计算每个细胞类型的平均表达（log1p标准化后的数据）
avg_exp <- AverageExpression(brain_ileum,
                             assays = "RNA",
                             slot = "data",      # 已标准化且log1p转换
                             group.by = "manual_celltype",
                             return.seurat = FALSE)$RNA

# 确保是矩阵格式
avg_exp <- as.matrix(avg_exp)

# 查看矩阵：行是基因，列是细胞类型
cat("Average expression matrix dimensions:\n")
print(dim(avg_exp))
cat("Column names (cell types):\n")
print(colnames(avg_exp))

# 3. 过滤低表达基因（至少在一类细胞中平均log1p表达 > 0.1）
keep <- apply(avg_exp, 1, max) > 0.1
avg_exp_filt <- avg_exp[keep, ]

cat("Filtered matrix dimensions:\n")
print(dim(avg_exp_filt))

# 可选：查看过滤后的基因数和细胞类型数
dim(avg_exp_filt)

# 4. 计算相关性矩阵（Pearson 和 Spearman）
cor_pearson <- cor(avg_exp_filt, method = "pearson")
cor_spearman <- cor(avg_exp_filt, method = "spearman")

# 5. 可视化
# 设置颜色条
my_colors <- colorRampPalette(c("#458f8fff", "white", "#6f4b8dff"))(100)

# 5.1 使用 corrplot 绘制混合图（上三角圆形，下三角数字）
pdf(file.path(output_dir, "celltype_corr_pearson_mixed.pdf"), width = 7.75, height = 7)
corrplot.mixed(cor_pearson,
               upper = "circle",
               lower = "number",
               tl.pos = "lt",
               diag = "u",
               tl.col = "black",
               upper.col = my_colors,
               lower.col = "black",
               number.cex = 1.2,
               number.digits = 2,
               mar = c(0, 0, 2, 0))
dev.off()

pdf(file.path(output_dir, "celltype_corr_spearman_mixed.pdf"), width = 7.75, height = 7)
corrplot.mixed(cor_spearman,
               upper = "circle",
               lower = "number",
               tl.pos = "lt",
               diag = "u",
               tl.col = "black",
               upper.col = my_colors,
               lower.col = "black",
               number.cex = 1.2,
               number.digits = 2,
               mar = c(0, 0, 2, 0))
dev.off()

# 5.2 使用 pheatmap 绘制带聚类的热图（显示数值）
pdf(file.path(output_dir, "celltype_corr_pearson_heatmap.pdf"), width = 7, height = 6)
pheatmap(cor_pearson,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = FALSE,
         number_format = "%.2f",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         color = my_colors,
         main = "Pearson correlation between cell types")
dev.off()

pdf(file.path(output_dir, "celltype_corr_spearman_heatmap.pdf"), width = 7, height = 6)
pheatmap(cor_spearman,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = FALSE,
         number_format = "%.2f",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         color = my_colors,
         main = "Spearman correlation between cell types")
dev.off()

cat("Correlation plots generated successfully!\n")
