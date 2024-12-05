# 安装必要的包
install.packages(c("readr", "ggplot2", "vegan"))  # 安装所需包
library(readr)
library(ggplot2)
library(vegan)

# 1. 读取数据
file_path <- "C:/Users/Lamarck/Desktop/pca.csv"  # 替换为实际文件路径
data <- read_csv(file_path)

# 设置基因名称为行名，并移除第一列
rownames(data) <- data$NAME
data <- data[, -1]  # 移除基因名称列

# 删除全为 0 的行
data <- data[rowSums(data != 0) > 0, ]

# 提取分组信息
group <- ifelse(grepl("^control", colnames(data)), "Control", "Treated")  # 根据列名生成分组信息

# 确保数据为数值型矩阵
data <- as.matrix(data)

# 2. 数据标准化
data_scaled <- t(scale(t(data)))  # 对行（基因）进行标准化

# 移除包含 NaN 的行
data_scaled <- data_scaled[complete.cases(data_scaled), ]

# 3. PCA 分析
pca_result <- prcomp(t(data_scaled), center = TRUE, scale. = TRUE)

# 查看 PCA 解释的方差比例
summary(pca_result)

# 提取 PCA 数据
pca_data <- data.frame(
  Sample = colnames(data),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = group
)

# 4. 组间统计：PERMANOVA
adonis_result <- adonis2(t(data_scaled) ~ group, method = "euclidean")
print(adonis_result)

# 提取 PERMANOVA 结果 P 值
adonis_pval <- adonis_result$aov.tab$`Pr(>F)`[1]

# 5. 绘制 PCA 图：标记样本点 + 95% 置信椭圆 + 组间统计 P 值
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  # 样本点
  geom_point(size = 4, alpha = 0.8) +
  # 样本标记
  geom_text(vjust = -1, size = 3) +
  # 添加 95% 置信椭圆
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group)) +
  # 添加标题和标签
  labs(
    title = "PCA Analysis of Gene Expression",
    subtitle = paste0("PERMANOVA P = ", format(adonis_pval, digits = 3)),  # 添加 P 值
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "% Variance Explained)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "% Variance Explained)")
  ) +
  # 设置主题
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  # 设置颜色和填充
  scale_color_manual(values = c("Control" = "red", "Treated" = "blue")) +
  scale_fill_manual(values = c("Control" = "red", "Treated" = "blue"))
