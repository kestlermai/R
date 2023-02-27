#install.packages("pheatmap") 
library(pheatmap)

setwd("C:/Users/11146/Desktop/NAFLD2.0/17 Vincent W 2013")

# 读取数据表
otu <- read.delim("otu_table_rare.txt",header=T, row.names=1)
group <- read.delim('17group1.txt', row.names = 1, check.names = FALSE)#group排序分组
otu$sum <- rowSums(otu)# 通常我们只关注高丰度且显著差异的，按每个OTU的总丰度排序
# 按每行的和降序排列
otu_order <- otu[order(otu$sum, decreasing=TRUE), ]

# 取丰度前40的OTUs
mat <- otu_order[1:40, -ncol(otu)]

# log2 转换，通常百万比经常log2转换，数据范围由1-1000000变为1-20
scale_test <- apply(mat, 2, function(x) log2(x+1))
# scale 转换
#scale_test <- apply(mat, 2, scale)
# 读取元数据，添加样本注释
annot_data <- data.frame(row.names=colnames(mat), group=group$group)
write.csv(scale_test, 'ASV高丰度表.csv', quote = FALSE)
# 读取物种分类添加行注释 
taxonomy <- read.delim("taxonomy-ASV-genus.txt", row.names=1, check.names = FALSE)#对应高丰度表;#Unclassfied
anno_row <- data.frame(row.names = rownames(mat),Taxonomy_Genus=taxonomy$Genus)

ann_colors = list(group = c("Healthy Control" ="#6B6ECFFF", NAFLD="#ff8696")) #定义分组颜色
ann_colors = list(group = c("Healthy Control"="#6B6ECFFF", NAFLD="#FDDDAE", "NAFLD-6months later"="#ff8696")) #FDDDAE/#B6A1DE
#"#E889BD", "#B286D7", "#5189E0", "#0089CF", "#0081A1", "#007360", "#E889BD", "#B286D7", "#5189E0"
# 绘制添加样本列注释的图
pheatmap(scale_test,scale="row", treeheight_row=100,#行树高度
         annotation_col=annot_data,
         angle_col = 0,#列标签角度
         annotation_row = anno_row,
         annotation_colors= ann_colors,
         #annotation_colors=	c("#7697CB","#F47E62"),#表示行注释的颜色
         #color = rainbow(2),
         color = colorRampPalette(c("#2574AA", "white", "#ED7B79"))(100),#颜色渐变，等级随意
         gaps_col = c(22,31),
         cluster_col= F,cluster_row=T,# cluster_row/cluster_col按行列聚类
         cutree_rows = 12,#根据样品列聚类情况将热图的行方向隔开为4份
         show_rownames=T,show_colnames=F,
         cellwidth=10, cellheight=10,#更改cell的高度和宽度
         filename="pheatmap_OTU_top40_sample_color.pdf", width=20, height=10)#show_rownames/show_colnames显示行/列名

