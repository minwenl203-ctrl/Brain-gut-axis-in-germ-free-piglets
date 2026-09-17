#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 差异基因 TPM 矩阵，分肠段 edgeR 差异分析
# 比较组: (1) D30C vs D0   (2) D30T vs D30C
# 统计各肠段上/下调基因数，绘制双向柱形图：
#   - x 轴上方 = 上调，下方 = 下调
#   - 颜色区分比较组（红 = D30C vs D0，蓝 = D30T vs D30C）
#   - 上方 = 上调，下方 = 下调；坐标轴置于底层，无柱边框
# 输出:
#   1) DEG_summary_by_segment.csv    -- 两比较组 x 肠段 上/下调基因数汇总
#   2) DEG_table_perSegment.csv      -- 全部基因 x 肠段 x 比较组 完整结果
#   3) DEG_count_segments.pdf / png  -- 双向柱形图
# ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
})


FC_CUT  <- 1      # |log2FC| > 1
FDR_CUT <- 0.05   # FDR < 0.05

d <- read.csv("gut_bulkRNA_seq_TPM_20260801.csv", row.names = 1, check.names = FALSE)

# 移除完全重复的列（D30C_IleMuco_B1_1 与 D30C_IleMuco_B1 数值完全相同）
d <- d[, colnames(d) != "D30C_IleMuco_B1_1", drop = FALSE]

# 肠段 -> 各时间点的样本名列前缀（D0 与 D30C/D30T 命名不一致：jejunum/JejMuco, ileum/IleMuco）
seg_map <- data.frame(
  segment = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Colon"),
  D0      = c("duo", "jejunum", "ileum", "cecum", "colon"),
  D30C    = c("duo", "JejMuco", "IleMuco", "cecum", "colon"),
  D30T    = c("duo", "JejMuco", "IleMuco", "cecum", "colon"),
  stringsAsFactors = FALSE
)

# 比较组: ref = 参照组, trt = 处理组; logFC = log2(trt / ref)
comps <- list(
  "D30C vs D0"  = c(ref = "D0",   trt = "D30C"),
  "D30T vs D30C" = c(ref = "D30C", trt = "D30T")
)

all_res <- list()
for (cn in names(comps)) {
  ref <- comps[[cn]][["ref"]]
  trt <- comps[[cn]][["trt"]]
  for (i in seq_len(nrow(seg_map))) {
    seg   <- seg_map$segment[i]
    c_ref <- grep(paste0("^", ref, "_", seg_map[[ref]][i], "_"), colnames(d), value = TRUE)
    c_trt <- grep(paste0("^", trt, "_", seg_map[[trt]][i], "_"), colnames(d), value = TRUE)
    sub   <- d[, c(c_ref, c_trt)]
    grp   <- factor(c(rep(ref, length(c_ref)), rep(trt, length(c_trt))),
                    levels = c(ref, trt))

    # TPM 取整作为计数近似输入 edgeR
    y <- DGEList(counts = round(sub), group = grp)
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y)
    et <- exactTest(y)               # logFC = log2(trt / ref)
    tt <- topTags(et, n = Inf)$table
    tt$gene       <- rownames(tt)
    tt$comparison <- cn
    tt$segment    <- seg
    tt$direction  <- ifelse(tt$FDR < FDR_CUT & tt$logFC >  FC_CUT, "Up",
                     ifelse(tt$FDR < FDR_CUT & tt$logFC < -FC_CUT, "Down", "NS"))
    all_res[[paste(cn, seg)]] <- tt

    cat(sprintf("%-14s %-9s | %s n=%d | %s n=%d | tested=%d | Up=%d | Down=%d\n",
                cn, seg, ref, length(c_ref), trt, length(c_trt), nrow(tt),
                sum(tt$direction == "Up"), sum(tt$direction == "Down")))
  }
}

res <- do.call(rbind, all_res)
res$segment    <- factor(res$segment, levels = seg_map$segment)
res$comparison <- factor(res$comparison, levels = names(comps))

# ---- 汇总统计 ----
sig <- res[res$direction != "NS", ]
tb  <- table(sig$comparison, sig$segment, sig$direction)
summ <- expand.grid(comparison = levels(res$comparison), segment = levels(res$segment))
summ$Up    <- as.integer(tb[, , "Up"])
summ$Down  <- as.integer(tb[, , "Down"])
summ$Total <- summ$Up + summ$Down

write.csv(summ, "DEG_summary_by_segment.csv", row.names = FALSE)
write.csv(res,  "DEG_table_perSegment.csv", row.names = FALSE)

# ---- 双向柱形图：上 = 上调，下 = 下调；颜色 = 比较组 ----
plt <- data.frame(
  comparison = summ$comparison,
  segment    = summ$segment,
  dir        = rep(c("Up", "Down"), each = nrow(summ)),
  n          = c(summ$Up, summ$Down)
)
plt$value <- ifelse(plt$dir == "Up", plt$n, -plt$n)
plt$x     <- as.numeric(plt$segment) +
             ifelse(plt$comparison == names(comps)[1], -0.21, 0.21)

p <- ggplot(plt, aes(x = x, y = value, fill = comparison)) +
  # 坐标轴置于底层：先画轴线与零线，再画柱子
  # 注意 limits 必须包住全部柱子（最左柱左缘 0.59，最右柱右缘 5.41），否则柱子会被裁剪移除
  geom_segment(aes(x = 0.5, xend = 5.5, y = -Inf, yend = -Inf),
               colour = "black", linewidth = 0.4) +
  geom_segment(aes(x = 0.5, xend = 0.5, y = -Inf, yend = Inf),
               colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_col(width = 0.4, colour = NA) +
  scale_fill_manual(values = c("D30C vs D0" = "#6674c5ff", "D30T vs D30C" = "#9fc5a6ff")) +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.15, 0.15))) +
  scale_x_continuous(breaks = seq_along(levels(plt$segment)),
                     labels = levels(plt$segment),
                     limits = c(0.5, 5.5)) +
  labs(x = NULL, y = "Number of DEGs", fill = NULL,
       title = "DEG counts by intestinal segment",
       subtitle = expression("edgeR exactTest, FDR < 0.05 & " * "|" * log[2] * "FC| > 1")) +
  theme_classic(base_size = 14) +
  theme(axis.line        = element_blank(),   # 轴线由底层 geom_segment 绘制
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top",
        axis.text.x   = element_text(size = 12, face = "bold"),
        axis.title.y  = element_text(size = 13))

ggsave("DEG_count_segments.pdf", p, width = 6.5, height = 6)
# ggsave("DEG_count_segments.png", p, width = 8.5, height = 6, dpi = 300)

cat("\n===== 汇总（FDR<0.05 & |log2FC|>1）=====\n")
print(summ)
cat("\nDone.\n")
