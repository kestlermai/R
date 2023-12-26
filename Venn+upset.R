library(VennDiagram)
library(RColorBrewer)
library(dplyr)
setwd("C:/Users/11146/Desktop/NAFLD2.0/21 Wang J 2017")

data.g <- read.table("otu_table_rare.txt" ,row.names = 1,sep = "\t")
group <- read.table("21group-rare.txt",header = T,sep = "\t")#样本的排序按照分组排列

g<-unique(group$group)
Gnum<-length(g)
number <- group %>% group_by(group) %>% count(group)
gf <- c()
k <- 1
for(i in 1:(Gnum)){
  gf[k] <- sum(number[1:i,2])
  k <- k+1
}
dg <- rownames(data.g)
k <- 1
for (i in 1:(Gnum)) {
a <- as.data.frame(apply(data.g[,c(k:gf[i])],1,sum))
dg <- cbind(dg,a)
k <- gf[i] +1
}
dg <- dg[,-1]
colnames(dg) <- g

dg <- ifelse(dg > 0, rownames(dg),NA)
dg <- as.data.frame(dg)
df <- list()
for(i in 1:length(dg)){
  a <- na.omit(dg[,i])
  df <- c(df,list(a))
}
names(df) <- colnames(dg)

p <- venn.diagram(df,filename = NULL,#height = 5400,width = 5400,
             resolution = 600,#imagetype = "tiff",
             units = "cm",
             lwd = 1,lty = 1,
             euler.d = F, scaled = F,#venn.diagram默认没有交集圆圈不重合
             fill = c("#3C5488B2","#DC0000B2"),#两变量
             #fill = c("#DC0000B2","#3C5488B2"),#两变量
             #fill = brewer.pal(length(dg),"Set2"),#三变量以上
             cex = 1,#数字大小
             fontface = 1.5,
             #cat.col = brewer.pal(length(dg),"Set2"),
             cat.col = c("#3C5488B2","#DC0000B2"),#字体颜色
             #cat.col = c("#DC0000B2","#3C5488B2"),
             #克莱因蓝#002EA6;松花黄#FFE76F;马尔斯绿#01847F;玫瑰粉#F9D2E4；
             #爱马仕橙#FF770F;深蓝色#000026;蒂芙尼蓝#80D1C8;奶酪色#F8F5D6;
             #淡黄色#FAEAD3;范戴克棕#492D22;浅卡其色#D8C7B5
             cat.dist = c(0.06,0.05),
             cat.pos = c(45,-45),
             cat.cex = 1.2,alpha = 0.8,margin = 0.05,
             cat.fontface = 1, print.mode = c("raw","percent"))
pdf("venn.pdf", width=6.5,height=6)
grid.draw(p)
dev.off()


#install.packages("ggVennDiagram")
library(ggVennDiagram)

# 查看交集详情,并导出结果
inter <- get.venn.partitions(df)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )
write.table(inter, "venn result.csv", row.names = FALSE, sep = ',', quote = FALSE)


#将venn图的数据集合可视化
library(UpSetR)
#install.packages("UpSetR")
#devtools::install_github("hms-dbmi/UpSetR")
upset(fromList(df),order.by="freq")

#个性化作图
pdf("UpSetR.pdf", width=10,height=5)
upset(fromList(df),
      sets = c("HEV", "HEV-ALF"),#指定展示的集合
      nset = 3, #集合范围的限制
      nintersects = 3, #需要绘制的交集数目
      order.by = c('degree','freq'), decreasing = c(F, T),#排序
      mb.ratio = c(0.7, 0.3),#柱状图与矩阵点图之间的比例大小
      number.angles = 0,#柱子上方数字倾斜角度
      point.size = 1.8,#矩阵中圆圈的大小
      line.size = 1, #矩阵点图中点和线的大小
      shade.color = "red",#矩阵点图阴影的颜色
      mainbar.y.label = "ASV  Intersections", #条形图的轴标签
      sets.x.label = "Set Size", #柱状图的轴标签
      main.bar.color = "#73BAD6", #柱状图柱子颜色
      sets.bar.color = "#3b7960",#条形图条形的颜色
      matrix.color = "#033250",#矩阵中交集点的颜色
      text.scale = c(1.3, 1.3, 1, 1, 1.2, 1),#文本大小设置
      keep.order = TRUE,#让集合按照 sets 参数中指定的出现的顺序排列
      queries = list(list(query = intersects,
                          params = list("HEV-ALF"),
                          active = T,color="#EF4143")))#设置查询条件 

dev.off()


#install.packages("ggvenn")
library(ggvenn)
#数据一
NAFLD <- read.delim("21-NAFLD-OTU.txt",header = T)
#数据二
HC <- read.delim("21-HC-OTU.txt",header = T)
mydata<-list(A=NAFLD$NAFLD,B=HC$HC)
#查看数列前留个数
head(mydata)

ggvenn(mydata, c("NAFLD", "Healthy Control"))
p <- ggvenn(mydata,c("NAFLD", "Healthy Control"),
         show_elements=FALSE,
         show_percentage=TRUE,
         digits=1,
         fill_color=c("#ffb2b2","#b2e7cb"),
         fill_alpha= 0.5,
         stroke_color="white",
         stroke_alpha =1,
         stroke_size = 1,
         stroke_linetype="solid",
         set_name_color=c("#ffb2b2","#b2e7cb"),
         set_name_size=7,
         text_color="black",
         text_size=3)
p
pdf(file="venn.pdf",width=8,height=8)
print(p)
dev.off()