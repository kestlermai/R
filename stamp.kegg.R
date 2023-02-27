library(reshape2)
library(doBy)
setwd("C:/Users/11146/Desktop/NAFLD数据2.0/2 Ayesha MK 2020")
#在 KEGG 第二层级统计求和
pathway <- read.delim("pathway_prediction.txt", check.names = FALSE)#列名为样本名，行名为level2通路
pathway2 <- melt(pathway, id = 'level2')#转换为短数据
pathway2 <- summaryBy(value~variable+level2, pathway2, FUN = sum)
#转化为 STAMP 输入文件格式
pathway2 <- dcast(pathway2, level2~variable)
write.table(pathway2, 'stamp.pathway2_level.txt', row.names = FALSE, sep = '\t', quote = FALSE)
#stamp路径不能有中文


group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
names(group)[1] <- 'variable'
tax4fun_pathway <- tax4fun_pathway[c(group$variable, 'Pathway2_level')]
tax4fun_pathway <- melt(tax4fun_pathway, id = 'Pathway2_level')
tax4fun_pathway <- summaryBy(value~variable+Pathway2_level, tax4fun_pathway, FUN = sum)
tax4fun_pathway <- merge(tax4fun_pathway, group, by = 'variable')
#计算均值±标准差
se <- function(x) sd(x) / (length(x))^0.5pathway_stat <- summaryBy(value.sum~group+Pathway2_level, tax4fun_pathway, FUN = c(mean, sd, se))

#可选进行显著性差异分析，这里直接使用 wilcoxon 秩和检验
Pathway2_level <- unique(tax4fun_pathway$Pathway2_level)
for (i in Pathway2_level) {
  tax4fun_pathway_2_i <- subset(tax4fun_pathway, Pathway2_level == i)
  test <- wilcox.test(value.sum~group, tax4fun_pathway_2_i)
  line_t <- which(pathway_stat$Pathway2_level == i & pathway_stat$group == 'treat')
  pathway_stat[line_t,'p_value'] <- test$p.value
  if (test$p.value < 0.05 & test$p.value >= 0.01) {
    pathway_stat[line_t,'sign'] <- '*'
  }
  if (test$p.value < 0.01 & test$p.value >= 0.001) {
    pathway_stat[line_t,'sign'] <- '**'
  }
  if (test$p.value < 0.001) {
    pathway_stat[line_t,'sign'] <- '***'
  }
}

library(ggplot2)

#使用 ggplot2 命令作图
pathway_stat$value.sum.mean <- 100 * pathway_stat$value.sum.mean
pathway_stat$value.sum.sd <- 100 * pathway_stat$value.sum.sd

pathway2_plot <- ggplot(pathway_stat, aes(Pathway2, value.sum.mean, fill = group)) +
  geom_col(position = 'dodge', width = 0.8, colour = 'black', size = 0.05) +  #“dodge 柱状图”样式
  geom_errorbar(aes(ymin = value.sum.mean - value.sum.sd, ymax = value.sum.mean + value.sum.sd), size = 0.05, width = .35, position = position_dodge(width = .8)) +  #添加误差线（均值±标准差）
  scale_fill_manual(values = c('red', 'blue')) +  #颜色填充
  theme(legend.title = element_blank(), legend.position = c(0.9, 0.95)) +  #去除图例标题，调整图例位置
  coord_flip() +  #将横轴和纵轴反转
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent',  color = 'black')) +  #去除默认的背景框
  geom_text(aes(label = sign, y = value.sum.mean + value.sum.sd + 0.5), size = 4, position = position_dodge(0.8)) +  #添加显著性标记“*”
  labs(x = 'KEGG pathway2', y = 'Relative Abundance (%)')  #添加每个坐标轴标题

ggsave('Kegg_pathway2.pdf', pathway2_plot, width = 8, height = 10)
ggsave('Kegg_pathway2.png', pathway2_plot, width = 8, height = 10)