library(vegan)
library(ggplot2)
library(ggsignif)

setwd("C:/Users/11146/Desktop/NAFLD2.0/11 Sonja L 2020")


#读取数据，一般所需是数据行名为样本名、列名为OTUxxx的数据表
otu_raw <- read.delim("otu_table_rare.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
#由于排序分析函数所需数据格式原因，需要对数据进行转置
otu <- t(otu_raw)

#计算bray_curtis距离
otu.distance <- vegdist(otu, method = 'bray')

#NMDS排序分析——vegan包中的metaMDS函数
df_nmds <- metaMDS(otu.distance, k = 2)
#结果查看——关注stress、points及species三个指标
summary(df_nmds)


#适用性检验——基于stress结果

#应力函数值（<=0.2合理）
df_nmds_stress <- df_nmds$stress
df_nmds_stress
#检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
stressplot(df_nmds)


#提取作图数据
df_points <- as.data.frame(df_nmds$points)
#添加samp1es变量
df_points$samples <- row.names(df_points)
#修改列名
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
head(df_points)

#计算每个样品的拟合度，图中圈越小表示拟合度越高。
gof <- goodness(df_nmds)
plot(df_nmds, type="t", main = "goodness of fit")
points(df_nmds, display="sites", cex=gof*100)

#绘散点图
p <- ggplot(df_points,aes(x=NMDS1, y=NMDS2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
p
ggsave("NMDS散点图.PDF", p, units = "cm", width = 13, height = 10,dpi = 1200)

#进行个性化展示

#读入分组文件
group <- read.delim("11group-rare.txt", sep='\t', header=T)
#修改列名
colnames(group) <- c("samples","group")
#将绘图数据和分组合并
df <- merge(df_points,group,by="samples")
head(df)
#使用ggplot2包绘图
color=c("#3C5488B2", "#DC0000B2","#1597A5","#FFC24B","#FEB3AE", "#8491B4B2")#颜色变量
p1<-ggplot(data=df,aes(x=NMDS1,y=NMDS2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = group), shape = 19, size=1)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  #geom_text(aes(label=group, y=NMDS2+0.03,x=NMDS1+0.03,
                #vjust=0, color = group),size=3.5, show.legend = F)+#添加数据点的标签
  stat_ellipse(data=df,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  scale_color_manual(values = color) +#点的颜色设置
  scale_fill_manual(values = color)+#椭圆颜色
  theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank())+#隐藏网格线
  ggtitle(paste('Sonja L 2020   Stress=',round(df_nmds_stress, 3)))#添加应力函数值
  #ggtitle(paste('Ayesha MK 2020   '))
p1
ggsave("NMDS分析1.PDF", p1, units = "cm", width = 13, height = 10,dpi = 1200)

#添加边际箱线图展示组间差异性
#加载包，对组间进行统计检验以及组合图的拼接
library(ggpubr)
library(ggsignif)
# 绘制y轴为PC2值的分组箱线图
p2 <- ggplot(df,aes(x=group,y=NMDS2))+
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
              y_position = c(0.8),
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
p3 <- ggplot(df,aes(x=group,y=NMDS1))+
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
              y_position = c(0.8),
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
t <- ggarrange(p3, NULL, p1, p2, widths = c(5,2), heights = c(2,4), align = "hv")
t
ggsave("NMDS分析3.PDF", t , units = "cm", width = 25, height = 15,dpi = 1200)


#anosim分析
anosim.result<-anosim(otu.distance,group$group,permutations = 999)#group要对应otu表，多余要删除
summary(anosim.result)
#结果表示R介于(-1，1)之间，R大于0，说明组间差异显著；
#R小于0，说明组内差异大于组间差异（如果R<0,那么分组也就没有太大的意义），
#统计分析的可信度用P表示，P< 0.05表示统计具有显著性。

#整理出作图数据
df1 <-data.frame(
  x=anosim.result$class.vec,
  y=anosim.result$dis.rank
)
#绘图
df2 <- ggplot(df1,aes(x=x,y=y))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
  geom_boxplot(aes(fill=x), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        legend.position="none",plot.title = element_text(size=14))+
  scale_fill_manual(values=c("darkseagreen","#3C5488B2", "#DC0000B2","#00A087B2", "#8491B4B2","#F39B7FB2","#91D1C2B2"))+ #指定颜色
  ggtitle("Bray-Curtis Anosim")+

  theme(legend.position = 'none')+
  labs(x = paste("R=",round(anosim.result$statistic,4),", ","p=", round(anosim.result$signif,4)),
       y = "Rank of Distance (Bray_Curtis)")
df2

ggsave("anosim分析.pdf", df2 , units = "cm", width = 13, height = 10 ,dpi = 1200)
