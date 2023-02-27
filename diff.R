setwd("C:/Users/11146/Desktop/NAFLD2.0/12 Cyrielle C 2019")

#设置工作路径session
library(microeco)
library(magrittr)
library(gridExtra)
library(ggplot2)
library(ape)
library(picante)
library(GUniFrac)
library(agricolae)
library(ggpubr)
library(ggdendro)
library(vegan)
library(SpiecEasi)
library(WGCNA)
library(rgexf)
library(randomForest)

group <- read.delim('12group-pre.txt', row.names = 1, check.names = FALSE)
otu_table <- read.delim('otu_table_rare.txt', row.names = 1, check.names = FALSE)
taxonomy <- read.delim('taxonomy.txt', row.names = 1, check.names = FALSE)
tree <- read.tree("rooted_tree.tre")

set.seed(123)#设置随机种子便于重复
theme_set(theme_bw())
taxonomy %<>% tidy_taxonomy

#归一化后的
dataset <- microtable$new(sample_table = group, otu_table = otu_table, tax_table = taxonomy, phylo_tree = tree)
class(dataset)
print(dataset)

#使otu、分类表、树表的行树一致
dataset$tidy_dataset()
print(dataset)

dataset$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
print(dataset)

dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
print(dataset)

dataset$tidy_dataset()
print(dataset)

dataset$sample_sums() %>% range

dataset$rarefy_samples(sample.size = 10000)
dataset$sample_sums() %>% range

dataset$cal_abund()
# return dataset$taxa_abund
class(dataset$taxa_abund)
#导出抽平后的otu表
write.table (dataset[["otu_table"]],file ="otu_table_rare.txt", sep ="\t", row.names = T)

##输出对应的丰度表
write.csv(dataset$taxa_abund$Phylum, '门水平丰度表.csv', quote = FALSE) #差异分析结果输出
write.csv(dataset$taxa_abund$Class, '纲水平丰度表.csv', quote = FALSE)
write.csv(dataset$taxa_abund$Order, '目水平丰度表.csv', quote = FALSE)
write.csv(dataset$taxa_abund$Family, '展科水平丰度表.csv', quote = FALSE)
write.csv(dataset$taxa_abund$Genus, '属水平丰度表.csv', quote = FALSE)
write.csv(dataset$taxa_abund$Species, '种水平丰度表.csv', quote = FALSE)


###绘制细菌构成比堆叠图
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 15,  groupmean = "group")
t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
t3 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
ggsave("堆积图门前15.png",t3,units = "cm",width = 15, height = 12,dpi = 600)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 15,  groupmean = "group")
t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
t3 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
ggsave("堆积图纲前15.png",t3,units = "cm",width = 15, height = 12,dpi = 600)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Order", ntaxa = 15,  groupmean = "group")
t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
t3 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
ggsave("堆积图目前15.png",t3,units = "cm",width = 15, height = 12,dpi = 600)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Family", ntaxa = 15,  groupmean = "group")
t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
t3 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
ggsave("堆积图科前15.png",t3,units = "cm",width = 15, height = 12,dpi = 600)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 15,groupmean = "group")
t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
t3 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
ggsave("堆积图属前15.png",t3,units = "cm",width = 15, height = 12,dpi = 600)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Species", ntaxa = 15,  groupmean = "group")
t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
t3 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
ggsave("堆积图种前15.png",t3,units = "cm",width = 15, height = 12,dpi = 600)

# show 15 taxa at Class level
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 15)
t2 <- t1$plot_box(group = "group")
t2
ggsave("箱式图门水平物种组成前15.png",t2,units = "cm",width = 20, height = 20,dpi = 600)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 15)
t2 <- t1$plot_box(group = "group")
t2
ggsave("箱式图属水平物种组成前15.png",t2,units = "cm",width = 20, height = 20,dpi = 600)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Species", ntaxa = 15)
t2 <- t1$plot_box(group = "group")
t2
ggsave("箱式图种水平物种组成前15.png",t2,units = "cm",width = 20, height = 20,dpi = 600)

#venn图
# merge samples as one community for each group
dataset1 <- dataset$merge_samples(use_group = "group")
# dataset1 is a new microtable object
# create trans_venn object
t1 <- trans_venn$new(dataset1, ratio = NULL)
t2 <- t1$plot_venn()
t2
ggsave("micro-venn.png",t2,units = "cm",width = 15 ,height = 15, dpi = 600)


#1.lefse分析
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "group", alpha = 0.05, lefse_subgroup = NULL)
#未校正P值的lefse
t1 <- trans_diff$new(dataset = dataset, method = "lefse", taxa_level = "Genus", group = "group", alpha = 0.05, lefse_subgroup = NULL, p_adjust_method = NULL)
t1 <- trans_diff$new(dataset = dataset, method = "lefse", taxa_level = "Species", group = "group", alpha = 0.05, lefse_subgroup = NULL, p_adjust_method = NULL)
t1$res_diff #is the LEfSe result
t1$plot_diff_bar(threshold = 2)
write.csv(t1$res_diff, 'lefse属.csv', quote = FALSE)
write.csv(t1$res_diff, 'lefse种.csv', quote = FALSE)
t2 <- t1$plot_diff_bar(threshold = 2)
ggsave("lefse-LDA=2属.png",t2,units = "cm",width = 30, height = 10,dpi = 600) 
ggsave("lefse-LDA=2种.png",t2,units = "cm",width = 30, height = 10,dpi = 600) 


#2.随机森林
t1 <- trans_diff$new(dataset = dataset, method = "rf", group = "group", taxa_level = "Species", p_adjust_method = NULL)
t1 <- trans_diff$new(dataset = dataset, method = "rf", group = "group", taxa_level = "Genus", p_adjust_method = NULL)
t1$res_diff #is the result stored in the object
g1 <- t1$plot_diff_bar(use_number = 1:15)
g2 <- t1$plot_diff_abund(use_number = 1:15, only_abund_plot = FALSE)
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
t3 <- gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 2))
ggsave("随机森林差异菌属前15.png",t3,units = "cm",width = 28, height = 12,dpi = 600)
ggsave("随机森林差异菌种前15.png",t3,units = "cm",width = 28, height = 12,dpi = 600)
write.csv(t1$res_diff, '随机森林分析属.csv', quote = FALSE) #差异分析结果输出
write.csv(t1$res_diff, '随机森林分析种.csv', quote = FALSE) 
