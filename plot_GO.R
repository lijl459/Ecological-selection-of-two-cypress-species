
library(ggplot2)
library(dplyr)
library(viridis)

# 读取数据文件
enrich_data <- read.csv("combined_go_results.csv")  # 这是包含两物种数据的文件

# 过滤出 p-value 小于 0.05 的行
filtered_data <- enrich_data %>%
  filter(classicFisher < 0.05)

# 计算富集比率
filtered_data$Ratio <- filtered_data$Significant / filtered_data$Annotated

# 对 p-value 进行 -log10 转换，调整颜色渐变的分布
filtered_data$log_pvalue <- log10(filtered_data$classicFisher)

# 设置 p-value 的中间值
mid_pvalue <- median(filtered_data$log_pvalue)


ggplot(filtered_data, aes(x = Significant, y = Term, color = log_pvalue, shape = Species)) +
  geom_point(aes(size = Ratio)) +  # 设置气泡的大小、颜色和形状
  scale_color_viridis(option = "A", direction = -1, begin = 0.8, end = 0, discrete = FALSE) +  # 使用 viridis 渐变色
  scale_size(range = c(2, 8)) +  # 控制气泡的大小范围
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   +
  labs(
    x = "Number of Significant Genes",
    y = NULL,
    size = "Ratio (Significant / Annotated)",
    color = "log10(p-value)",
    shape = "Species"
  ) +
  ggtitle("GO Enrichment for Biologiacal Process")  # 图表标题



ggsave("Go_enrichment_bubble_plot.jpg", width = 12, height = 7, dpi = 300)
ggsave("Go_enrichment_bubble_plot.pdf", width = 12, height = 7)




