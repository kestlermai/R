#install.packages("ggalluvial")
#install.packages("d3Network")
#install.packages("ggforce")
#install.packages("networkD3")
#读写索引失败
#setRepositories(addURLs = c(CRANxtras = "http://cran.at.r-project.org/"))
#多选几个源
#网页交互式
library(d3Network)

setwd("C:/Users/11146/Desktop")
sankey <- read.table("test3.txt",header = TRUE,sep = "\t")
nodes <- data.frame(name = unique(c(as.character(sankey$source),
                                    as.character(sankey$target))),stringsAsFactors = FALSE)
nodes$ID <- 0:(nrow(nodes)-1)
sankey <- merge(sankey,nodes,by.x = "source",by.y = "name")
sankey <- merge(sankey,nodes,by.x = "target",by.y = "name")
colnames(sankey) <- c("X","Y","value","source","target")
sankey <- subset(sankey,select = c("source","target","value"))
nodes <- subset(nodes,select = c("name"))

d3Sankey(Links = sankey,Nodes = nodes,fontsize = 20,
         Source = "source",Target = "target",Value = "value",
         NodeID = "name",file = "Sankey.html",
         width = 1200,height = 900)

#ggplot画
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(ggsci)

setwd("C:/Users/11146/Desktop")
sankey <- read.table("test2.txt", header = T, sep="\t", check.names=F)

corLodes=to_lodes_form(sankey, axes = 1:ncol(sankey), id = "Cohort")

mycol <- rep(c("#F9837B","#F5886C","#F18C5A","#EC9042","#E79419","#E09719","#DA9C19","#D49F19","#CCA319","#C4A619",
               "#BBA919","#B1AC19","#A8B019","#9CB219","#8FB519","#81B819","#70BA19","#59BC19","#30BE19","#19C043",
               "#19C25A","#19C36B","#19C47A","#19C587","#19C694","#19C6A0","#19C7AB","#19C6B6","#19C6C1","#19C5CA",
               "#19C4D4","#19C3DC","#19C0E4","#19BEEC","#19BAF2","#19B7F9","#19B3FD","#19AEFF","#58A9FF","#7AA4FF",
               "#939EFF","#A798FF","#B892FF","#FF79A3","#FF76AF","#FF73BA","#FF71C4","#FF70CE","#FC70D8","#F971E0",
               "#F474E8","#EE78F0","#E77BF7","#DE81FC","#D386FF","#C68CFF","#A953FF","#FF4374"),58)
mycol <- rep(c("#DCDCDC","#98F5FF","#FFEFD5","#FFDAB9","#E6E6FA","#76EEC6","#FFF0F5","#FFE4E1","#2F4F4F","#696969",
               "#708090","#BEBEBE","#000080","#8470FF","#B3EE3A","#FFC125","#CD5C5C","#F4A460","#B22222","#FA8072",
               "#FFA500","#CDBA96","#FF69B4","#FFC0CB","#C71585","#8A2BE2","#8B4C39"),27)

p <- ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +  
  #用aes.flow控制线调颜色，forward说明颜色和前面一致，backward说明与后面一致
  geom_flow(width = 0.2,#连线宽度
            aes.flow = "forward",#向前流动
            show.legend = TRUE) + 
  geom_stratum(alpha = 0.9,#透明度
               width = 0.3, # 方格宽度
              linetype=1,size=1.0, color = "#696969",inherit.aes=TRUE) + #画冲击图
  scale_fill_manual(values = mycol) +
  geom_text(stat = "stratum", size = 1.5,color="black") +#字体大小
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  ggtitle("") + guides(fill = FALSE)   
p
p1 <- p +	scale_color_nejm() 
p1 <- p1 + scale_fill_nejm()
p1

#输出图片至本地
ggsave('test.pdf', p, width = 10, height = 6)
#ggsave('test.png', p, width = 10, height = 5)


library(ggplot2)
library(ggforce)
setwd("C:/Users/11146/Desktop")
sankey <- read.table("test2.txt", header = T, sep="\t", check.names=F)
data$x <- factor(sankey$x,levels = c("Survived","gene","sample","age","gender","stage"),ordered = F)
p <- ggplot(data, aes(x, id = id, split = y, value = value)) +
            geom_parallel_sets(aes(fill = Survived), alpha = 0.3, axis.width = 0.1) +
            geom_parallel_sets_axes(axis.width = 0.2) +
            geom_parallel_sets_labels(colour = 'white')
p


library(networkD3)
library(Hmisc)

setwd("C:/Users/11146/Desktop/200502-3个示例学会R包networkD3和plotly的交互式桑基图")
#“env_table.txt” 记录了测量的环境变量数据
env <- read.delim('env_table.txt', row.name = 1, check.names = FALSE)
env <- t(env)
#量纲不同，作个标准化，尽管标准化前后对相关系数的计算无影响
env <- scale(env)
#环境变量间的相关性，以 spearman 秩相关系数为例
env_corr <- rcorr(t(env), type = 'spearman')

r <- env_corr$r  #相关系数 r 值
p <- env_corr$P  #相关系数的显著性 p 值
p <- p.adjust(p, method = 'BH')  #可选 p 值校正，这里使用 BH 法校正 p 值

#仅保留 r>=0.7 且 p.adj<0.01 的相关系数，且去除对角线的自相关
r[abs(r)<0.7 | p>=0.01] <- 0
r[1:6,1:6]
diag(r) <- 0

#输出相关系数矩阵
write.table(r, 'env_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#使用 corrplot 包绘制相关图可视化相关矩阵
library(corrplot)
col1 <- colorRampPalette(c('blue4', 'blue', 'white', 'orange', 'red3'))

corrplot(r, method = 'number', col = col1(21), 
         tl.cex = 0.8, number.cex = 0.8, cl.length = 5, diag = FALSE)
corrplot(r, method = 'pie',col = col1(21), 
         tl.pos = 'n', cl.pos = 'n', add = TRUE, type = 'upper', diag = FALSE)

##构建桑基图数据
#读取刚才计算的相关矩阵
r <- read.delim('env_corr.matrix.txt', check.names = FALSE)

#由于这是一个对称矩阵，仅保留下三角数值即可，把上三角的冗余值去除
r[upper.tri(r)] <- 0

#将相关系数矩阵转为“两两关系列表”的样式
env_corr_list <- reshape2::melt(r)
names(env_corr_list) <- c('source', 'target', 'spearman')
env_corr_list$source <- as.character(env_corr_list$source)
env_corr_list$target <- as.character(env_corr_list$target)

#取相关系数的绝对值赋值新列作为相关的强度，并取相关系数的符号赋值新列作为相关的方向
#以及去除 0 值
env_corr_list$weight <- abs(env_corr_list$spearman)
env_corr_list[which(env_corr_list$spearman>0),'direction'] <- '1'
env_corr_list[which(env_corr_list$spearman<0),'direction'] <- '-1'
env_corr_list <- subset(env_corr_list, !env_corr_list$direction %in% NA)
env_corr_list$label <- as.character(round(env_corr_list$spearman, 3))

#读取环境变量的属性列表，分配 id 指代环境变量，并用在后续根据分类赋值颜色等
env_group <- read.delim('env_group.txt', stringsAsFactors = FALSE, check.names = FALSE)
env_group <- subset(env_group, env %in% unique(c(env_corr_list$source, env_corr_list$target)))
env_corr_list$IDsource <- match(env_corr_list$source, env_group$env) - 1 
env_corr_list$IDtarget <- match(env_corr_list$target, env_group$env) - 1

head(env_corr_list)

#输出相关性列表
write.table(env_corr_list, 'env_corr.list.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#定义节点和连线的颜色
color <- 'd3.scaleOrdinal() .domain(["pH", "C", "N", "P", "K", "1", "-1"]) 
.range(["#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#FF00004D", "#0000FF4D"])'

#指定预定义的变量 id，按相关性方向赋值连线颜色，相关性强度赋值连线尺寸
#节点颜色按预定义的变量的分类着色
#节点、字体、边界尺寸等，详情 ?sankeyNetwork
p <- sankeyNetwork(Links = env_corr_list, Nodes = env_group,
                   Source = 'IDsource', Target = 'IDtarget',  
                   Value = 'weight', LinkGroup = 'direction', 
                   NodeID = 'env', NodeGroup = 'group', 
                   nodePadding = 50, nodeWidth = 20, fontSize = 12,
                   colourScale = color, height = 300, width = 800)

p
saveNetwork(p,file="sankey.html",selfcontained=T)#不独立保存为sankey文件夹

#install.packages('igraph',type='binary')
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("igraph")

library(networkD3)
library(ggplot2)

setwd("C:/Users/11146/Desktop")
nodes <- read.csv("A.csv",header = T, fileEncoding = "UTF-8-BOM")
links <- read.csv("B.csv",header = T, fileEncoding = "UTF-8-BOM")

links$IDsource <- match(links$source, nodes$name) -1
links$IDtarget <- match(links$target, nodes$name) -1
head(links)

my_color <- 'd3.scaleOrdinal() .domain(["c1", "c2","g1","g2"]) .range(["darkorange", "palegreen","skyblue","lightcoral"])'  

p<- sankeyNetwork(Links = links, #含有nodes(即source和target)连接信息和value值的数据框
                  Nodes= nodes, #含有nodes属性的数据框
                  Source= "IDsource", #sourceID
                  Target= "IDtarget", #targetID
                  Value= "value",
                  NodeID= "name",units = "TWh",
                  fontSize= 25, nodeWidth = 35, 
                  colourScale=my_color,
                  #sinksRight,布尔值，若为T，则将最后一个节点移动到绘图的右边界
                  NodeGroup="group", LinkGroup="group") 
p
saveNetwork(p,file="sankey.html",selfcontained=T)#不独立保存为sankey文件夹
