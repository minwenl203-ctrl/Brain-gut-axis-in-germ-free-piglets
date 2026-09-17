# 程序功能：scRNA-seq 分析通用绘图参数与分组信息
# Functions: shared plotting themes, group orders and color palettes for figures 1-6

# 1. ggplot2 通用主题（类似 example 模板中的 main_theme）
main_theme = theme(panel.background = element_blank(),
                   panel.grid = element_blank(),
                   axis.line.x = element_line(size = .5, colour = "black"),
                   axis.line.y = element_line(size = .5, colour = "black"),
                   axis.ticks = element_line(color = "black"),
                   axis.text = element_text(color = "black", size = 7),
                   legend.position = "right",
                   legend.background = element_blank(),
                   legend.key = element_blank(),
                   legend.text = element_text(size = 7),
                   text = element_text(family = "sans", size = 7))

# 2. 分组顺序（GF_D0 为无菌对照，GF_D30 为无菌 30 天，MMT_D30 为粪菌移植 30 天）
group_order <- c("GF_D0", "GF_D30", "MMT_D30")

# 3. 胶质细胞亚型（fig6 Ro/e 分析用）
subtype_keep <- c("Oligodendrocytes", "Microglia", "Astrocytes")

# 4. 细胞类型颜色（与 07_plot.R 一致，fig3 火山图用）
celltype_cols <- c("Oligodendrocytes" = "#E3EDC0", "Microglia" = "#F2CA8D", "Astrocytes" = "#C1E6F3")

# 5. 相关性热图颜色（fig1 用）
cor_colors <- colorRampPalette(c("#458f8fff", "white", "#6f4b8dff"))(100)
