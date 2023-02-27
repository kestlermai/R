setwd("C:/Users/11146/Desktop/NAFLD数据2.0/6 Xiongfeng P 2021")
#install.packages("BiocManager")
#####可视化分析
#设置工作路径session
library(fastcluster)
library(dynamicTreeCut)
library(ggpubr)
library(agricolae)
library(nlme)
library(lattice)
library(permute)
library(vegan)
library(microeco)
library(magrittr)
library(gridExtra)
library(ggplot2)
library(ape)
library(picante)
library(GUniFrac)
library(ggdendro)
library(vegan)
library(SpiecEasi)
library(WGCNA)
library(rgexf)

#归一化，按最小样本序列数
otu_table <- read.delim('otu_table.txt', row.names = 1, check.names = FALSE)
colSums(otu_table)
otu_table_rare = as.data.frame(t(rrarefy(t(otu_table), min(colSums(otu_table)))))
colSums(otu_table_rare)
write.table (otu_table_rare,file ="otu_table_rare.xls", sep ="\t", row.names = T)
#抽平一般在5000以上甚至10000/30000，一般是按最小值进行抽平，但抽平深度太小或太高都会影响结果

####数据导入及处理
group <- read.delim('11group.txt', row.names = 1, check.names = FALSE)
otu_table_rare <- read.delim('otu_table_rare.txt', row.names = 1, check.names = FALSE)
taxonomy <- read.delim('taxonomy.txt', row.names = 1, check.names = FALSE)
tree <- read.tree("rooted_tree.tre")

set.seed(123)#设置随机种子便于重复
theme_set(theme_bw())
taxonomy %<>% tidy_taxonomy

#归一化
dataset <- microtable$new(sample_table = group, otu_table = otu_table_rare, tax_table = taxonomy, phylo_tree = tree)
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

#data normalization
#install.packages("devtools")
#library(devtools)
#devtools::install_github("ChiLiubio/mecodev")
library(mecodev)
dataset_norm <- trans_norm$new(dataset = dataset)
# Centered log-ratio normalization. 中心对数比归一化
dataset <- dataset_norm$norm(method = "CLR")
# Cumulative sum scaling normalization. 累积和缩放标准化 Require metagenomeSeq package to be installed.
#BiocManager::install("metagenomeSeq") 
#用normalization实现分类注释时的biases处理，后续也可以用这个包构建线性模型进行物种差异分析，类似lefse
#安装前需要limma、Biobase、BiocGenerics包
library(metagenomeSeq)
dataset <- dataset_norm$norm(method = "CCS")
# log transformation. 对数变换
dataset <- dataset_norm$norm(method = "log")

write.table (dataset[["otu_table"]],file ="otu_table_norm.txt", sep ="\t", row.names = T)
#####alpha多样性分析
# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = TRUE)
# return dataset$alpha_diversity
class(dataset$alpha_diversity)

# save dataset$alpha_diversity to a directory
dir.create("alpha_diversity2")
dataset$save_alphadiv(dirpath = "alpha_diversity2")

#绘制Alpha多样性分析图
t1 <- trans_alpha$new(dataset = dataset, group = "group")
# return t1$alpha_stat
t1$alpha_stat[1:5, ]
t1$cal_diff(method = "KW")
# return t1$res_alpha_diff
t1$res_alpha_diff[1:5, ]
t1$cal_diff(method = "anova")
# return t1$res_alpha_diff
t1$res_alpha_diff
#t1$plot_alpha(add_letter = TRUE, measure = "PD")

#observed
t1$plot_alpha(pair_compare = TRUE, measure = "Observed", title = "PRJNA728736 2021-Feces")+grids(linetype="dashed")+border("black")
t2 <- t1$plot_alpha(pair_compare = TRUE, measure = "Observed", title = "PRJNA728736 2021-Feces")+grids(linetype="dashed")+border("black")
t3 <- ggpar(t2, font.title = c(10, "black")) #"bold.italic", 
ggsave("alpha-Observed.png",t3,units = "cm",width = 8, height = 11,dpi = 1200)

#ACE
t2 <- t1$plot_alpha(pair_compare = TRUE, measure = "ACE", title = "PRJNA728736 2021-Feces")+grids(linetype="dashed")+border("black")
t3 <- ggpar(t2, font.title = c(10, "black")) #"bold.italic", 
ggsave("alpha-ACE.png",t3,units = "cm",width = 8, height = 11,dpi = 1200)

#Chao1
t2 <- t1$plot_alpha(pair_compare = TRUE, measure = "Chao1", title = "PRJNA728736 2021-Feces")+grids(linetype="dashed")+border("black")
t3 <- ggpar(t2, font.title = c(10, "black")) #"bold.italic", 
ggsave("alpha-Chao1.png",t3,units = "cm",width = 8, height = 11,dpi = 1200)

#Shannon
t2 <- t1$plot_alpha(pair_compare = TRUE, measure = "Shannon", title = "PRJNA728736 2021-Feces")+grids(linetype="dashed")+border("black")
t3 <- ggpar(t2, font.title = c(10, "black")) #"bold.italic", 
ggsave("alpha-Shannon.png",t3,units = "cm",width = 8, height = 11,dpi = 1200)

#Simpson
t2 <- t1$plot_alpha(pair_compare = TRUE, measure = "Simpson", title = "PRJNA728736 2021-Feces")+grids(linetype="dashed")+border("black")
t3 <- ggpar(t2, font.title = c(10, "black")) #"bold.italic", 
ggsave("alpha-Simpson.png",t3,units = "cm",width = 8, height = 11,dpi = 1200)

#InvSimpson
t2 <- t1$plot_alpha(pair_compare = TRUE, measure = "InvSimpson", title = "PRJNA728736 2021-Feces")+grids(linetype="dashed")+border("black")
t3 <- ggpar(t2, font.title = c(10, "black")) #"bold.italic", 
ggsave("alpha-InvSimpson.png",t3,units = "cm",width = 8, height = 11,dpi = 1200)

#Fisher
t2 <- t1$plot_alpha(pair_compare = TRUE, measure = "Fisher", title = "PRJNA728736 2021-Feces")+grids(linetype="dashed")+border("black")
t3 <- ggpar(t2, font.title = c(10, "black")) #"bold.italic", 
ggsave("alpha-Fisher.png",t3,units = "cm",width = 8, height = 11,dpi = 1200)

#PD
t2 <- t1$plot_alpha(pair_compare = TRUE, measure = "PD", title = "PRJNA728736 2021-Feces")+grids(linetype="dashed")+border("black")
t3 <- ggpar(t2, font.title = c(10, "black")) #"bold.italic", 
ggsave("alpha-PD.png",t3,units = "cm",width = 8, height = 11,dpi = 1200)

####我们还使用函数cal_betadiv（）计算β分集的距离矩阵。Beta多样性分析。
#提供了四个最常用的索引：Bray-curtis，Jaccard，加权Unifrac和未加权unifrac。
# If you do not want to calculate unifrac metrics, use unifrac = FALSE
# require GUniFrac package
dataset$cal_betadiv(unifrac = TRUE)
# return dataset$beta_diversity 
class(dataset$beta_diversity)
# save dataset$beta_diversity to a directory
dir.create("beta_diversity")
dataset$save_betadiv(dirpath = "beta_diversity")

# we first create an object and select PCoA for ordination
t1 <- trans_beta$new(dataset = dataset, group = "group", measure = "bray" )#,  ordination = "PCoA" )

t1$cal_ordination(ordination = "PCoA", ncomp = 3, trans_otu = FALSE,scale_species = FALSE)#"PCoA"; "PCA", "PCoA" or "NMDS"

# t1$res_ordination is the ordination result list
class(t1$cal_ordination)
# plot the PCoA result
t1$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point" , "ellipse") )
t2 <- t1$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point" , "ellipse") )
t2
ggsave("bray-pcoa.tiff",t2,units = "cm",width = 15, height = 12,dpi = 1200)

# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")
t3 <- t1$plot_group_distance(boxplot_add = "mean")
t3
ggsave("bray-pcoa差异比较.tiff",t3,units = "cm",width = 9, height = 12,dpi = 1200)

#perMANOVA通常用于组间距离的差异测试中
# manova for all groups
t1$cal_manova(manova_all = TRUE)
t1$res_manova

t1$cal_manova(cal_manova_paired = TRUE)
t1$res_manova

#PCoA
library(vegan)
library(ggplot2)
color=c( "#3C5488B2", "#DC0000B2",
         "#00A087B2", "#8491B4B2", 
         "#F39B7FB2","#91D1C2B2", 
         "#7E6148B2","yellow", 
         "darkolivegreen1", "lightskyblue", 
         "darkgreen", "deeppink", "khaki2", 
         "firebrick", "brown1", "darkorange1", 
         "cyan1", "royalblue4", "darksalmon", 
         "darkgoldenrod1", "darkseagreen", "darkorchid")

#读取otu数据文件
otu <- read.delim('otu_table_rare.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))

#根据物种组成计算样方距离，结果以 dist 数据类型存储
bray_dis <- vegdist(otu, method = 'bray')
write.csv(as.matrix(bray_dis),'bray_dis.csv')
data=read.delim('otu_table_rare.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#样方排序坐标
pcoa <- cmdscale(bray_dis, k = (nrow(otu) - 1), eig = TRUE)
site <- data.frame(pcoa$point)[1:2]
site$name <- rownames(site)

#读取分组文件
map<-read.table('11group.txt',header=T,sep="\t",row.names=1)
#site$group <- c(rep('A', 12), rep('B', 12), rep('C', 12))
#将分组文件和数据文件以行名合并
merged=merge(site,map,by="row.names",all.x=TRUE)

#使用 wascores() 被动添加物种得分（坐标），丰度加权平均方法
species <- wascores(pcoa$points[,1:4], otu)

#计算 top10 丰度物种
abundance <- apply(otu, 2, sum)
abundance_top10 <- names(abundance[order(abundance, decreasing = TRUE)][1:10])

species_top10 <- data.frame(species[abundance_top10,1:2])
species_top10$name <- rownames(species_top10)

pcoa_exp <- pcoa$eig/sum(pcoa$eig)
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')

#ggplot2 作图
p <- ggplot(data = merged, aes(X1, X2)) +
  geom_point(aes(color = group)) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +   #添加置信椭圆，注意不是聚???
  scale_color_manual(values =color[1:length(unique(map$group))]) +
  scale_fill_manual(values =color[1:length(unique(map$group))]) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_text(data = species_top10, aes(label = name), color = "gray", size = 0) +    #??? top10 丰度物种标签
  labs(x = pcoa1, y = pcoa2, title = 'Sonja L 2020  perMANOVA=0.045*') 
p
ggsave("bray-pcoa.PDF", p, units = "cm", width = 13, height = 10,dpi = 1200)

#添加边际箱线图展示组间差异性
#加载包，对组间进行统计检验以及组合图的拼接
library(ggpubr)
library(ggsignif)
# 绘制y轴为PC2值的分组箱线图
p2 <- ggplot(merged,aes(x=group,y=X2))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(#panel.background =element_blank(), 
    panel.background = element_rect(color = 'black', fill = 'transparent'),#加边框
    axis.line=element_line(color = "black",lineend = "square"),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values=c("#3C5488B2", "#DC0000B2","#00A087B2", "#8491B4B2","#F39B7FB2","#91D1C2B2"))+
  geom_signif(comparisons = list(c("Healthy Control","NAFLD")),
              #c("HC","HBV-LC"),
              #c("HC","HBV-HCC"),
              #c("CHB","HBV-LC"),
              #c("CHB","HBV-HCC"),
              #c("HBV-LC","HBV-HCC")),
              map_signif_level = T, #TRUE显示星号，FALSE显示数字
              #test = t.test,
              test = wilcox.test,#非参数
              y_position = c(0.4),
              #c(0.55,0.6,0.67,0.41,0.35,0.48),#坐标要改
              tip_length = c(c(0,0)),
              #c(0,0),
              #c(0,0),
              #c(0,0),
              #c(0,0),
              # c(0,0)),
              size=0.5,color="black")
p2
# 绘制y轴为PC1值的分组箱线图
p3 <- ggplot(merged,aes(x=group,y=X1))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  coord_flip()+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(#panel.background =element_blank(), 
    panel.background = element_rect(color = 'black', fill = 'transparent'),#加边框
    axis.line=element_line(color = "black",lineend = "square"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values=c("#3C5488B2", "#DC0000B2","#00A087B2", "#8491B4B2","#F39B7FB2","#91D1C2B2"))+
  geom_signif(comparisons = list(c("Healthy Control","NAFLD")),
              #c("HC","HBV-LC"),
              #c("HC","HBV-HCC"),
              #c("CHB","HBV-LC"),
              #c("CHB","HBV-HCC"),
              #c("HBV-LC","HBV-HCC")),
              map_signif_level = T, #TRUE显示星号，FALSE显示数字
              #test = t.test,
              test = wilcox.test,#非参数
              y_position = c(0.5),
              #c(0.55,0.6,0.67,0.41,0.35,0.48),#坐标要改
              tip_length = c(c(0,0)),
              #c(0,0),
              #c(0,0),
              #c(0,0),
              #c(0,0),
              # c(0,0)),
              size=0.5,color="black")
p3

# ggpubr::ggarrange()函数对图进行拼接
t <- ggarrange(p3, NULL, p, p2, widths = c(5,2), heights = c(2,4), align = "hv")
t
ggsave("PCoA分析.PDF", t , units = "cm", width = 25, height = 18,dpi = 1200)
